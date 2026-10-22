#!/usr/bin/env python3
"""Publish final timing results; keep reconstruction and diagnostics in Combine."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile
from datetime import datetime, timezone


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def load_json(path):
    with path.open(encoding='utf-8') as stream:
        return json.load(stream)


def relative_child(parent, name):
    value = Path(name)
    if value.is_absolute() or '..' in value.parts:
        raise ValueError(f'Unsafe relative path in execution manifest: {name}')
    candidate = parent / value
    if not candidate.resolve().is_relative_to(parent.resolve()):
        raise ValueError(f'Manifest path escapes its directory: {name}')
    return candidate


def execution_for(run, stage, execution_id=None):
    if execution_id is None:
        pointer = load_json(run / 'metadata' / f'latest_{stage}.json')
        execution_id = pointer['execution_id']
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', execution_id):
        raise ValueError(f'Invalid {stage} execution ID')
    execution = run / 'metadata' / 'executions' / execution_id
    manifest = load_json(execution / 'manifest.json')
    result = manifest.get('stages', {}).get(stage, {})
    # An all-mode run may have a successful signal stage even when its
    # independent noise stage failed. Validate the selected stage explicitly.
    if result.get('status') != 'completed' or result.get('exit_code', 0) != 0:
        raise ValueError(f'Latest {stage} stage did not complete successfully')
    config = relative_child(execution, manifest.get('config', 'config'))
    if not config.is_dir():
        raise ValueError(f'Missing {stage} execution config snapshot: {config}')
    expected_config = manifest.get('config_sha256', {})
    actual_config = {p.relative_to(config).as_posix(): sha256(p)
                     for p in config.rglob('*') if p.is_file()}
    if not expected_config or actual_config != expected_config:
        raise ValueError(f'{stage} execution config snapshot changed after analysis')
    records = manifest.get('inputs', {}).get(stage, [])
    if not records:
        raise ValueError(f'No {stage} inputs recorded by execution {execution_id}')
    expected_paths = set()
    for record in records:
        source = Path(record['path'])
        linked = Path(record.get('run_path', record['path']))
        for path in {source, linked}:
            if not path.is_file() or sha256(path) != record['sha256']:
                raise ValueError(f'{stage} input changed since analysis: {path}')
        # Calibration is also a signal input, but is not a physics acquisition.
        # Older manifests did not record roles and contained only acquisitions.
        if record.get('role', 'acquisition') == 'acquisition':
            expected_paths.add(linked.resolve())
    folder = run / 'Tracker' / ('physics' if stage == 'signal' else 'noise')
    present_paths = {p.resolve() for p in folder.iterdir()
                     if p.suffix.lower() == '.root' and p.is_file()}
    if not expected_paths or expected_paths != present_paths:
        raise ValueError(f'{stage} input file list changed since analysis: {folder}')
    outputs = {}
    for name, digest in manifest.get('outputs', {}).items():
        is_noise = name.startswith('Combine/noise/')
        if is_noise != (stage == 'noise'):
            continue
        if not name.startswith('Combine/'):
            continue
        path = relative_child(run, name)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f'{stage} output changed since analysis: {path}')
        outputs[name] = path
    if not outputs:
        raise ValueError(f'No verified {stage} output recorded in execution {execution_id}')
    if stage == 'signal':
        for name in ('Tracker-step1-raw.root', 'Tracker-step2-ana.root',
                     'Tracker-step3-dec.root', 'Tracker-step4-rec.root',
                     'Pad-step3-rec.root'):
            if f'Combine/{name}' not in outputs:
                raise ValueError(f'Signal result is incomplete: missing {name}')
    return {'id': execution_id, 'directory': execution, 'manifest': manifest,
            'config': config, 'inputs': records, 'outputs': outputs}


def copy_file(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def output_destination(name):
    """Explicit final-result allowlist; a new diagnostic is never auto-published."""
    parts = Path(name).parts
    if len(parts) < 3 or parts[:2] != ('Combine', 'time_resolution'):
        return None
    relative = Path(*parts[1:])
    tail = parts[2:]
    if tail[0] == 'snr_binned':
        if len(tail) == 3 and tail[-1] in {
                'snr_bins.tsv', 'snr_thresholds.tsv', 'snr_bands.tsv',
                'snr_timing_report.json', 'snr_threshold_scan.png',
                'snr_threshold_scan.pdf'}:
            return 'timing/' + relative.as_posix()
        if (len(tail) in (5, 6) and tail[2] == 'pairwise_details'
                and tail[-2].startswith('SNR_')
                and tail[-1] == 'detector_solution.tsv'):
            # Accept both new run-local SNR_* paths and historical run_NNN/SNR_*
            # paths; publish only the useful final table, without a second run.
            return ('timing/' + Path(*parts[1:4], 'detector_solutions',
                                     tail[-2], tail[-1]).as_posix())
        return None
    if tail[0] == 'cluster':
        tail = tail[2:]  # cluster/<execution_id>/<summary or CFD directory>
    if len(tail) == 2 and tail[0] in ('t30', 't40', 't50'):
        tail = tail[1:]
    if len(tail) == 1 and tail[0] in {
            'detector_time_resolution.tsv', 'pairwise_time_resolution.tsv',
            'time_resolution.json', 'cfd_fraction_resolution.tsv',
            'cfd_fraction_resolution.png', 'cfd_fraction_resolution.pdf'}:
        return 'timing/' + relative.as_posix()
    return None


def write_hashes(destination, root, files):
    with destination.open('w', encoding='utf-8') as stream:
        for path in sorted(files):
            stream.write(f'{sha256(path)}  {path.relative_to(root).as_posix()}\n')


def export(run, tag, run_id=None, execution_id=None):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', tag):
        raise ValueError('Invalid result tag')
    results = run / 'Results'
    target = results / tag
    if target.exists() or target.is_symlink():
        raise ValueError(f'Refusing to overwrite published result: {target}')
    signal = execution_for(run, 'signal', execution_id)
    for name in ('detector_time_resolution.tsv', 'time_resolution.json'):
        if f'Combine/time_resolution/{name}' not in signal['outputs']:
            raise ValueError(f'No verified final timing result: missing {name}')
    noise = None
    noise_status = 'unavailable: no successful noise execution'
    try:
        noise = execution_for(run, 'noise')
        noise_status = 'metadata only; calibration products remain in Combine/noise'
    except (ValueError, KeyError, OSError) as error:
        noise_status = f'unavailable: {error}'
    results.mkdir(exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f'.{tag}.export-', dir=results))
    try:
        metadata = staging / 'metadata'
        metadata.mkdir()
        published = {}
        for stage, record in (('signal', signal), ('noise', noise)):
            if record is None:
                continue
            snapshot = metadata / ('config_snapshot' if stage == 'signal' else 'noise_config_snapshot')
            shutil.copytree(record['config'], snapshot)
            copy_file(record['directory'] / 'manifest.json', metadata / f'{stage}_manifest.json')
            for name, source in record['outputs'].items():
                destination = output_destination(name)
                if destination is None:
                    continue
                copy_file(source, staging / destination)
                published[destination] = name
            # Hash the files after copying too, so a concurrently rerun source
            # cannot silently produce an export inconsistent with its manifest.
            for name, digest in record['manifest']['outputs'].items():
                destination = output_destination(name)
                if (name in record['outputs'] and destination is not None
                        and sha256(staging / destination) != digest):
                    raise ValueError(f'{stage} output changed during export: {name}')
        details = {
            'analysis_tag': tag, 'run_id': run_id or run.name, 'run_directory': str(run),
            'exported_utc': datetime.now(timezone.utc).isoformat(),
            'signal_execution_id': signal['id'],
            'noise_execution_id': noise['id'] if noise else None,
            'noise_status': noise_status,
            'configuration_source': 'successful execution snapshot',
            'result_scope': 'final timing tables, reports and summary plots',
            'intermediates_directory': str(run / 'Combine'),
            'published_files': published,
            'timing_available': True,
        }
        (metadata / 'analysis.json').write_text(json.dumps(details, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
        input_records = {'signal': signal['inputs'], 'noise': noise['inputs'] if noise else []}
        (metadata / 'inputs.json').write_text(json.dumps(input_records, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
        with (metadata / 'raw_sha256.txt').open('w', encoding='utf-8') as stream:
            for stage, records in input_records.items():
                for record in records:
                    stream.write(f"{record['sha256']}  {record.get('run_path', record['path'])}  # {stage}\n")
        config_files = list((metadata / 'config_snapshot').rglob('*'))
        config_files += list((metadata / 'noise_config_snapshot').rglob('*'))
        write_hashes(metadata / 'config_sha256.txt', staging, (p for p in config_files if p.is_file()))
        index = [f'# {run_id or run.name}: {tag}', '',
                 'Final timing results / 最终时间分辨结果。', '',
                 'Step1–4、Pad ROOT、噪声标定、残差/时间差分布和逐事件计算过程保存在',
                 '[Combine](../../Combine/)。配置快照、来源及 SHA-256 保存在 [metadata](metadata/)。', '',
                 '## 结果文件', '']
        index += [f'- [{name}]({name})' for name in sorted(published)]
        (staging / 'README.md').write_text('\n'.join(index) + '\n', encoding='utf-8')
        write_hashes(metadata / 'result_sha256.txt', staging, (p for p in staging.rglob('*')
                     if p.is_file() and p != metadata / 'result_sha256.txt'))
        if target.exists() or target.is_symlink():
            raise ValueError(f'Refusing to overwrite published result: {target}')
        staging.rename(target)
    except BaseException:
        if staging.is_dir() and staging.parent == results:
            shutil.rmtree(staging)
        raise
    print(f'Published immutable result directory: {target}')
    print(f"Signal execution: {signal['id']}; noise: {noise_status}")
    print(f'Final timing files: {len(published)}; intermediates remain in {run / "Combine"}')
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('analysis_tag')
    parser.add_argument('run_id', nargs='?')
    parser.add_argument('--execution-id', help='export this execution instead of latest_signal')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', args.analysis_tag):
        parser.error('analysis_tag must contain only letters, digits, dots, underscores and hyphens')
    try:
        run = args.run_directory.resolve(strict=True)
        if not run.is_dir():
            raise ValueError(f'Not a run directory: {run}')
        export(run, args.analysis_tag, args.run_id, args.execution_id)
    except (ValueError, KeyError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

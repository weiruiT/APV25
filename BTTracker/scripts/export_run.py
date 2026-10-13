#!/usr/bin/env python3
"""Publish an immutable, verified snapshot from a managed BTTracker run."""

import argparse
import hashlib
import json
import os
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


def execution_for(run, stage):
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
        expected_paths.add(linked.resolve())
    folder = run / 'Tracker' / ('physics' if stage == 'signal' else 'noise')
    present_paths = {p.resolve() for p in folder.iterdir()
                     if p.suffix.lower() == '.root' and p.is_file()}
    if expected_paths != present_paths:
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
    relative = name.removeprefix('Combine/')
    if relative.startswith('noise/'):
        return relative
    if relative.startswith('time_resolution/'):
        return 'timing/' + relative
    if relative.startswith('sigmoid_only/'):
        return 'timing/' + relative.removeprefix('sigmoid_only/')
    if relative == 'peak_summary.tsv' or Path(relative).name.startswith('peaks_'):
        return 'diagnostics/peaks/' + Path(relative).name
    if Path(relative).name.startswith('MMalignmentPar'):
        return 'metadata/alignment/' + Path(relative).name
    return 'analysis/' + relative


def write_hashes(destination, root, files):
    with destination.open('w', encoding='utf-8') as stream:
        for path in sorted(files):
            stream.write(f'{sha256(path)}  {path.relative_to(root).as_posix()}\n')


def export(run, tag, run_id):
    results = run / 'Results'
    target = results / tag
    if target.exists() or target.is_symlink():
        raise ValueError(f'Refusing to overwrite published result: {target}')
    signal = execution_for(run, 'signal')
    noise = None
    noise_status = 'unavailable: no successful noise execution'
    try:
        noise = execution_for(run, 'noise')
        noise_status = 'included'
    except (ValueError, KeyError, OSError) as error:
        noise_status = f'unavailable: {error}'
    results.mkdir(exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f'.{tag}.export-', dir=results))
    try:
        metadata = staging / 'metadata'
        metadata.mkdir()
        for stage, record in (('signal', signal), ('noise', noise)):
            if record is None:
                continue
            snapshot = metadata / ('config_snapshot' if stage == 'signal' else 'noise_config_snapshot')
            shutil.copytree(record['config'], snapshot)
            copy_file(record['directory'] / 'manifest.json', metadata / f'{stage}_manifest.json')
            for name, source in record['outputs'].items():
                copy_file(source, staging / output_destination(name))
            # Hash the files after copying too, so a concurrently rerun source
            # cannot silently produce an export inconsistent with its manifest.
            for name, digest in record['manifest']['outputs'].items():
                if name in record['outputs'] and sha256(staging / output_destination(name)) != digest:
                    raise ValueError(f'{stage} output changed during export: {name}')
        details = {
            'analysis_tag': tag, 'run_id': run_id or run.name, 'run_directory': str(run),
            'exported_utc': datetime.now(timezone.utc).isoformat(),
            'signal_execution_id': signal['id'],
            'noise_execution_id': noise['id'] if noise else None,
            'noise_status': noise_status,
            'configuration_source': 'successful execution snapshot',
            'bm_alignment': 'included' if (metadata / 'alignment' / 'MMalignmentParBM.dat').is_file()
                            else 'not produced (optional)',
            'peak_diagnostic_pngs': len(list((staging / 'diagnostics' / 'peaks').glob('peaks_*.png'))),
            'timing_available': ((staging / 'timing/time_resolution/bm2_tracker_time_differences.root').is_file() or
                                 (staging / 'timing/delta_histograms.root').is_file()),
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
    print(f"Peak diagnostic PNGs: {details['peak_diagnostic_pngs']}; timing: {details['timing_available']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('analysis_tag')
    parser.add_argument('run_id', nargs='?')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', args.analysis_tag):
        parser.error('analysis_tag must contain only letters, digits, dots, underscores and hyphens')
    try:
        run = args.run_directory.resolve(strict=True)
        if not run.is_dir():
            raise ValueError(f'Not a run directory: {run}')
        export(run, args.analysis_tag, args.run_id)
    except (ValueError, KeyError, OSError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

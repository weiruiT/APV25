#!/usr/bin/env python3
"""Run-local APV25 acquisition import and batch execution (Python stdlib only)."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import uuid

PROJECT = Path(__file__).resolve().parents[1]
MODES = ("signal", "noise", "all", "amplitude", "energy", "arrival")
SIGNAL_FILES = [f"Tracker-step{i}-{suffix}.root" for i, suffix in
                ((1, "raw"), (2, "ana"), (3, "dec"), (4, "rec"))] + ["Pad-step3-rec.root"]
TIMING_FILES = [
    "time_resolution/pairwise_time_resolution.tsv",
    "time_resolution/detector_time_resolution.tsv",
    "time_resolution/time_resolution.json",
    "time_resolution/bm2_tracker_gaussian_fits.tsv",
    "time_resolution/bm2_tracker_time_differences.png",
    "time_resolution/bm2_tracker_time_differences.pdf",
    "time_resolution/bm2_tracker_time_differences.root",
    "time_resolution/BM2_Tracker1_gaussian_fit.png",
    "time_resolution/BM2_Tracker2_gaussian_fit.png",
    "time_resolution/BM2_Tracker3_gaussian_fit.png",
]


def stamp():
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ") + "_" + uuid.uuid4().hex[:8]


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
    temporary.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def object_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def deep_merge(base, override):
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_merge(base[key], value)
        else:
            base[key] = value
    return base


def generated_calibration(name):
    return (name.startswith("noise-global-parameters-") or
            name.startswith("fit_parameters.") or
            name in ("Tracker-noise-channel-parameters.txt",
                     "BM-noise-channel-parameters.txt"))


def config_snapshot(run, destination):
    # Pedestal/noise products belong to Combine/noise for each run. Excluding
    # historical generated products keeps snapshots small and unambiguous.
    shutil.copytree(PROJECT / "config", destination,
                    ignore=lambda _directory, names:
                    {name for name in names
                     if generated_calibration(name) or name.endswith(".md")})
    override = run / "config"
    if override.is_dir():
        for item in sorted(override.rglob("*")):
            if not item.is_file() or item.suffix.lower() == ".md":
                continue
            if generated_calibration(item.name):
                raise ValueError(
                    f"generated noise calibration belongs in Combine/noise, not run/config: {item}")
            target = destination / item.relative_to(override)
            target.parent.mkdir(parents=True, exist_ok=True)
            if item.name == "configCR.json":
                base = json.loads(target.read_text()) if target.exists() else {}
                value = json.loads(item.read_text())
                if not isinstance(base, dict) or not isinstance(value, dict):
                    raise ValueError("configCR.json must be a JSON object")
                write_json(target, deep_merge(base, value))
            else:
                shutil.copy2(item, target)
    return {str(p.relative_to(destination)): digest(p)
            for p in sorted(destination.rglob("*")) if p.is_file()}


def inputs_for(run, mode):
    folder = run / "Tracker" / ("physics" if mode == "signal" else "noise")
    if not folder.is_dir():
        return []
    files = sorted(p for p in folder.iterdir() if p.suffix == ".root")
    if not files:
        return []
    for path in files:
        if not path.is_file():
            raise ValueError(f"missing input or broken symlink: {path}")
        if "-Tracker-" in path.name or "-Pad-" in path.name:
            raise ValueError(f"analysis output in input directory: {path}")
    records = [{"path": str(p.resolve()), "run_path": str(p), "bytes": p.stat().st_size,
                "sha256": digest(p), "role": "acquisition"} for p in files]
    if mode == "signal":
        calibration = run / "Combine/noise/noise-global-parameters-1778.txt"
        if not calibration.is_file():
            raise ValueError(
                "signal peak selection requires Combine/noise/"
                "noise-global-parameters-1778.txt; run noise mode first or use --mode all")
        records.append({"path": str(calibration.resolve()),
                        "run_path": str(calibration),
                        "bytes": calibration.stat().st_size,
                        "sha256": digest(calibration),
                        "role": "peak_noise_calibration"})
    return records


def code_hash():
    files = [PROJECT / "ReadTracker.C", PROJECT / "CMakeLists.txt"]
    for folder in ("src", "include", "scripts"):
        files.extend(p for p in (PROJECT / folder).rglob("*")
                     if p.is_file() and p.suffix in (".C", ".h", ".py"))
    return object_hash({str(p.relative_to(PROJECT)): digest(p) for p in sorted(files)})


@contextmanager
def run_lock(run):
    lockfile = run / "metadata" / "run.lock"
    lockfile.parent.mkdir(parents=True, exist_ok=True)
    with lockfile.open("a+") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"run is already processing: {run}") from exc
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def build():
    directory = Path(os.environ.get("BTTRACKER_BUILD_DIR", PROJECT / "build_apv25")).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    # Build once before batch workers; flock also protects simultaneous launchers.
    with (directory / ".build.lock").open("a+") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        subprocess.run(["cmake", "-S", str(PROJECT), "-B", str(directory)], check=True)
        subprocess.run(["cmake", "--build", str(directory), "-j",
                        os.environ.get("BTTRACKER_BUILD_JOBS", "4")], check=True)
    return directory / "ReadTracker"


def expected_outputs(run, mode):
    if mode == "signal":
        return [run / "Combine" / name for name in SIGNAL_FILES]
    if mode == "noise":
        return [p for p in sorted((run / "Combine" / "noise").glob("*")) if p.is_file()]
    return []


def successful_previous(run, mode, fingerprint):
    pointer = run / "metadata" / f"latest_{mode}.json"
    if not pointer.is_file():
        return None
    value = json.loads(pointer.read_text())
    execution = run / "metadata" / "executions" / value["execution_id"]
    previous = json.loads((execution / "manifest.json").read_text())
    stage = previous.get("stages", {}).get(mode, {})
    if stage.get("status") != "completed" or stage.get("fingerprint") != fingerprint:
        return None
    outputs = stage.get("outputs", {})
    if not outputs or any(not (run / name).is_file() or digest(run / name) != checksum
                          for name, checksum in outputs.items()):
        return None
    return previous


def run_stage(command, cwd, logfile):
    with logfile.open("wb") as stream:
        result = subprocess.run(command, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT)
    # ROOT can install a signal handler that exits zero on fatal errors.
    with logfile.open(errors="replace") as stream:
        fatal = any("segmentation violation" in line.lower() or
                    "segmentation fault" in line.lower() or
                    line.lstrip().startswith("ERROR:") for line in stream)
    return result.returncode if result.returncode else (1 if fatal else 0)


def execute_run(run, args, executable, source_hash):
    result = {"run": str(run), "status": "failed", "execution_id": None}
    try:
        with run_lock(run):
            execution_id = stamp()
            execution = run / "metadata" / "executions" / execution_id
            execution.mkdir(parents=True)
            result["execution_id"] = execution_id
            manifest = {"run": str(run), "execution_id": execution_id,
                        "mode": args.mode, "force": args.force, "status": "running",
                        "stages": {}, "inputs": {}, "outputs": {}, "code_sha256": source_hash,
                        "config": "config", "started_utc": datetime.now(timezone.utc).isoformat()}
            manifest_path = execution / "manifest.json"
            write_json(manifest_path, manifest)
            try:
                config_hashes = config_snapshot(run, execution / "config")
                manifest["config_sha256"] = config_hashes
                analysis_config = json.loads(
                    (execution / "config/configCR.json").read_text()
                )
                bm2_matched_pad_count = analysis_config.get(
                    "timing_cuts", {}
                ).get("bm2_matched_pad_count", 0)
                timing_detectors = analysis_config.get(
                    "timing_cuts", {}
                ).get(
                    "detectors",
                    ["BM1", "BM2", "Tracker1", "Tracker2", "Tracker3"],
                )
                if (
                    isinstance(bm2_matched_pad_count, bool)
                    or not isinstance(bm2_matched_pad_count, int)
                    or bm2_matched_pad_count < 0
                ):
                    raise ValueError(
                        "timing_cuts.bm2_matched_pad_count must be an integer >= 0"
                    )
                valid_timing_detectors = {
                    "BM1", "BM2", "Tracker1", "Tracker2", "Tracker3"
                }
                if (
                    not isinstance(timing_detectors, list)
                    or len(timing_detectors) < 3
                    or len(set(timing_detectors)) != len(timing_detectors)
                    or any(
                        not isinstance(name, str)
                        or name not in valid_timing_detectors
                        for name in timing_detectors
                    )
                    or "BM2" not in timing_detectors
                ):
                    raise ValueError(
                        "timing_cuts.detectors must contain BM2 and at least "
                        "two other unique known detectors"
                    )
                (run / "Combine").mkdir(exist_ok=True)
                if args.alignment_source:
                    source = Path(args.alignment_source).resolve(strict=True)
                    for name in ("MMalignmentPar3.dat", "MMalignmentPar6.dat", "MMalignmentParBM.dat"):
                        target = run / "Combine" / name
                        if (source / name).is_file() and not target.exists():
                            shutil.copy2(source / name, target)
                modes = ["noise", "signal"] if args.mode == "all" else [args.mode]
                for mode in modes:
                    print(f"[{run.name}] {mode}: starting ({execution_id})", flush=True)
                    stage = {"status": "running", "log": f"{mode}.log"}
                    manifest["stages"][mode] = stage
                    try:
                        inputs = inputs_for(run, mode) if mode in ("signal", "noise") else []
                        manifest["inputs"][mode] = inputs
                        if mode in ("signal", "noise") and not inputs:
                            raise ValueError(f"no input files for {mode}: {run}/Tracker")
                        fingerprint = object_hash({"inputs": inputs, "config": config_hashes,
                                                   "code": source_hash, "mode": mode})
                        stage["fingerprint"] = fingerprint
                        previous = None
                        if args.force == 0 and mode in ("signal", "noise"):
                            previous = successful_previous(run, mode, fingerprint)
                            present = any(p.exists() for p in expected_outputs(run, mode))
                            if present and previous is None:
                                raise ValueError("existing outputs cannot be verified against current inputs/config/code; use force=1")
                        write_json(manifest_path, manifest)
                        if previous:
                            stage["reused_execution"] = previous["execution_id"]
                            stage["exit_code"] = 0
                            (execution / f"{mode}.log").write_text(f"Verified and reused {previous['execution_id']}\n")
                        else:
                            stage["exit_code"] = run_stage(
                                [str(executable), str(run), str(args.force), mode],
                                execution, execution / f"{mode}.log")
                        if stage["exit_code"] != 0:
                            raise RuntimeError(f"{mode} failed with exit {stage['exit_code']}; see {execution / (mode + '.log')}")
                        files = expected_outputs(run, mode)
                        if mode == "signal" and any(not p.is_file() for p in files):
                            raise RuntimeError("signal stage did not produce all five required ROOT outputs")
                        if mode == "noise" and not (run / "Combine/noise/noise-global-parameters-1778.root").is_file():
                            raise RuntimeError("noise stage did not produce noise-global-parameters-1778.root")
                        if mode == "signal":
                            if not previous:
                                stage["timing_log"] = "timing.log"
                                stage["timing_exit_code"] = run_stage(
                                    [
                                        os.environ.get(
                                            "BTTRACKER_PYTHON", sys.executable
                                        ),
                                        str(PROJECT / "scripts/analyze_time_resolution.py"),
                                        str(run / "Combine/Pad-step3-rec.root"),
                                        "--output-dir",
                                        str(run / "Combine/time_resolution"),
                                        "--bm2-matched-pad-count",
                                        str(bm2_matched_pad_count),
                                        "--detectors",
                                        ",".join(timing_detectors),
                                    ],
                                    execution,
                                    execution / "timing.log",
                                )
                                if stage["timing_exit_code"] != 0:
                                    raise RuntimeError(
                                        "time-resolution analysis failed with exit "
                                        f"{stage['timing_exit_code']}; see "
                                        f"{execution / 'timing.log'}"
                                    )
                            timing_files = [
                                run / "Combine" / name for name in TIMING_FILES
                            ]
                            missing_timing = [
                                str(path) for path in timing_files if not path.is_file()
                            ]
                            if missing_timing:
                                raise RuntimeError(
                                    "time-resolution analysis is incomplete: "
                                    + ", ".join(missing_timing)
                                )
                            files += timing_files
                            files += sorted((run / "Combine").glob("MMalignmentPar*.dat"))
                        stage["outputs"] = {str(p.relative_to(run)): digest(p) for p in files}
                        stage["status"] = "completed"
                        manifest["outputs"].update(stage["outputs"])
                        print(f"[{run.name}] {mode}: completed", flush=True)
                    except Exception as exc:
                        stage.update(status="failed", error=str(exc))
                        print(f"[{run.name}] {mode}: {exc}", file=sys.stderr, flush=True)
                    write_json(manifest_path, manifest)
                if args.plots and manifest["stages"].get("signal", {}).get("status") == "completed":
                    macro = PROJECT.parent / "debug/macros/plot_peak_timing.sh"
                    rc = run_stage(["bash", str(macro), str(run / "Combine")], execution, execution / "plots.log")
                    manifest["stages"]["plots"] = {"status": "completed" if rc == 0 else "failed",
                                                    "exit_code": rc, "log": "plots.log"}
                    if rc == 0:
                        artifacts = list((run / "Combine").glob("peaks_*.png"))
                        artifacts += list((run / "Combine").glob("peak_summary.tsv"))
                        artifacts += list((run / "Combine/sigmoid_only").rglob("*"))
                        manifest["outputs"].update({str(p.relative_to(run)): digest(p) for p in artifacts if p.is_file()})
                stage_statuses = [stage["status"] for stage in manifest["stages"].values()]
                if stage_statuses and all(status == "completed" for status in stage_statuses):
                    manifest["status"] = "completed"
                elif any(status == "completed" for status in stage_statuses):
                    manifest["status"] = "partial_failed"
                else:
                    manifest["status"] = "failed"
            except Exception as exc:
                manifest.update(status="failed", error=str(exc))
            manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
            write_json(manifest_path, manifest)
            for mode in ("signal", "noise"):
                if manifest["stages"].get(mode, {}).get("status") == "completed":
                    write_json(run / "metadata" / f"latest_{mode}.json", {"execution_id": execution_id})
            result.update(status=manifest["status"], manifest=str(manifest_path),
                          stages=manifest["stages"])
            if "error" in manifest:
                result["error"] = manifest["error"]
    except Exception as exc:
        result["error"] = str(exc)
        print(f"[{run.name}] {exc}", file=sys.stderr)
    return result


def discover(paths, pattern):
    runs = set()
    for value in paths:
        path = Path(value).resolve(strict=True)
        if (path / "Tracker/physics").is_dir() or (path / "Tracker/noise").is_dir():
            runs.add(path)
        else:
            candidates = sorted(path.glob(pattern))
            for candidate in candidates:
                if not candidate.is_dir() or not (candidate / "metadata/run.json").is_file():
                    continue
                acquisition = list((candidate / "Tracker/physics").glob("*.root"))
                acquisition += list((candidate / "Tracker/noise").glob("*.root"))
                if acquisition:
                    runs.add(candidate.resolve())
    if not runs:
        raise ValueError("no runs found; pass a run directory or its batches parent")
    return sorted(runs)


def init_run(args):
    run = Path(args.directory).absolute()
    existing = run.exists()
    if existing and not (run / "metadata/run.json").is_file():
        raise ValueError(f"existing directory is not a managed run: {run}")
    if existing and ((run / "metadata/latest_signal.json").exists() or
                     (run / "metadata/latest_noise.json").exists()):
        raise ValueError(f"cannot add acquisition files after analysis has completed: {run}")
    imports = []
    for mode, sources in (("physics", args.signal), ("noise", args.noise)):
        names = set()
        for source in sources:
            source = Path(source).resolve(strict=True)
            candidates = sorted(p for p in source.iterdir() if p.suffix == ".root") if source.is_dir() else [source]
            if not candidates:
                raise ValueError(f"no decoded ROOT files in {source}")
            for path in candidates:
                if not path.is_file() or path.suffix != ".root":
                    raise ValueError(f"not a supported input: {path}")
                if path.name in names:
                    raise ValueError(f"duplicate input name in {mode}: {path.name}")
                if "-Tracker-" in path.name or "-Pad-" in path.name:
                    raise ValueError(f"analysis result is not an acquisition input: {path}")
                names.add(path.name)
                imports.append((mode, path))
    for folder in ("Tracker/physics", "Tracker/noise", "Combine/noise", "Results", "metadata", "config"):
        (run / folder).mkdir(parents=True, exist_ok=True)
    run_record_path = run / "metadata/run.json"
    run_record = json.loads(run_record_path.read_text()) if existing else {
        "run_id": run.name, "created_utc": datetime.now(timezone.utc).isoformat(),
        "notes": args.notes, "imports": []}
    planned = [(mode, source, run / "Tracker" / mode / source.name)
               for mode, source in imports]
    for _mode, _source, target in planned:
        if target.exists() or target.is_symlink():
            raise ValueError(f"refusing to overwrite run input: {target}")
    records = []
    for mode, source, target in planned:
        if args.copy:
            shutil.copy2(source, target)
        else:
            target.symlink_to(source)
        records.append({"kind": mode, "source": str(source), "destination": str(target.relative_to(run)),
                        "method": "copy" if args.copy else "symlink",
                        "imported_utc": datetime.now(timezone.utc).isoformat()})
    run_record.setdefault("imports", []).extend(records)
    if args.notes:
        run_record["notes"] = args.notes
    write_json(run_record_path, run_record)
    action = "Updated" if existing else "Created"
    print(f"{action} {run}; imported {len(records)} files ({'copy' if args.copy else 'symlink'})")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    init = sub.add_parser("init", help="create one run; no automatic pairing")
    init.add_argument("directory")
    init.add_argument("--signal", nargs="+", action="extend", default=[])
    init.add_argument("--noise", nargs="+", action="extend", default=[])
    init.add_argument("--copy", action="store_true")
    init.add_argument("--notes", default="")
    batch = sub.add_parser("batch", help="process runs, compile once, continue after failure")
    batch.add_argument("paths", nargs="+")
    batch.add_argument("--mode", choices=MODES, default="all")
    batch.add_argument("--force", type=int, choices=(0, 1), default=1)
    batch.add_argument("--jobs", type=int, default=1)
    batch.add_argument("--pattern", default="run_*")
    batch.add_argument("--dry-run", action="store_true")
    batch.add_argument("--plots", action="store_true")
    batch.add_argument("--alignment-source")
    single = sub.add_parser("run", help="old run.sh positional interface")
    single.add_argument("path")
    single.add_argument("force", type=int, choices=(0, 1), nargs="?", default=1)
    single.add_argument("mode", choices=MODES, nargs="?", default="signal")
    args = parser.parse_args(argv)
    try:
        if args.command == "init":
            return init_run(args)
        if args.command == "run":
            path = Path(args.path).resolve(strict=True)
            if path.is_file():
                if args.mode == "all":
                    raise ValueError("all mode requires a paired run directory")
                executable = build()
                with tempfile.TemporaryDirectory(prefix="bttracker_direct_") as temporary:
                    work = Path(temporary)
                    config_snapshot(work / "empty_run", work / "config")
                    return subprocess.run(
                        [str(executable), str(path), str(args.force), args.mode],
                        cwd=work).returncode
            args.paths, args.jobs, args.pattern = [str(path)], 1, "run_*"
            args.dry_run, args.plots, args.alignment_source = False, False, None
        if args.jobs < 1:
            raise ValueError("jobs must be >= 1")
        if args.plots and args.mode not in ("signal", "all"):
            raise ValueError("--plots requires signal or all mode")
        runs = discover(args.paths, args.pattern)
        if args.command == "run" and (len(runs) != 1 or runs[0] != Path(args.path).resolve()):
            raise ValueError("run.sh takes one run; use run_batch.sh for a batches directory")
        if args.dry_run:
            for run in runs:
                print(f"{run}\tmode={args.mode}\tphysics={len(list((run / 'Tracker/physics').glob('*.root')))} ROOT\tnoise={len(list((run / 'Tracker/noise').glob('*.root')))} ROOT")
            return 0
        executable = build()
        source_hash = code_hash()
        summary_base = Path(args.paths[0]).resolve()
        summary_dir = summary_base / "metadata/batch_executions" / stamp()
        summary_dir.mkdir(parents=True)
        results = []
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(execute_run, run, args, executable, source_hash) for run in runs]
            for future in as_completed(futures):
                results.append(future.result())
                write_json(summary_dir / "summary.json", sorted(results, key=lambda r: r["run"]))
        results.sort(key=lambda r: r["run"])
        (summary_dir / "summary.tsv").write_text("run\tstatus\texecution_id\n" + "".join(
            f"{r['run']}\t{r['status']}\t{r.get('execution_id') or ''}\n" for r in results))
        print(f"Batch summary: {summary_dir / 'summary.tsv'}")
        return 0 if all(r["status"] == "completed" for r in results) else 1
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())

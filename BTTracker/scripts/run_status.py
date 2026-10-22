#!/usr/bin/env python3
"""Generate a compact index of managed runs and their latest execution."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

STAGES = ("noise", "signal", "cluster_timing", "snr_timing")
COLUMNS = ("run", "latest_execution", "mode", "status", *STAGES,
           "snr_bin_edges", "results", "physics_roots", "noise_roots", "error")


def read_json(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def root_count(path):
    return sum(1 for item in path.iterdir() if item.is_file() and item.suffix == ".root") if path.is_dir() else 0


def run_row(run):
    metadata = run / "metadata"
    if not (metadata / "run.json").is_file():
        return None
    executions = metadata / "executions"
    manifests = []
    if executions.is_dir():
        for folder in executions.iterdir():
            if folder.is_dir():
                manifest = read_json(folder / "manifest.json")
                if isinstance(manifest, dict):
                    manifests.append((str(manifest.get("execution_id", folder.name)), manifest))
    execution_id, manifest = max(manifests, key=lambda pair: pair[0]) if manifests else ("", {})
    stages = manifest.get("stages", {}) if isinstance(manifest.get("stages", {}), dict) else {}
    snapshot = read_json(executions / execution_id / "config/configCR.json") if execution_id else None
    snr = snapshot.get("snr_timing", {}) if isinstance(snapshot, dict) else {}
    edges = snr.get("bin_edges", "") if isinstance(snr, dict) else ""
    pointer = read_json(metadata / "latest_results.json")
    result_path = run / pointer.get("directory", "") if isinstance(pointer, dict) else None
    results = pointer.get("execution_id", "") if isinstance(pointer, dict) and result_path and result_path.is_dir() else ""
    error = manifest.get("error", "")
    failed = [stage.get("error", "") for stage in stages.values()
              if isinstance(stage, dict) and stage.get("status") == "failed" and stage.get("error")]
    if not error and failed:
        error = "; ".join(failed)
    return {
        "run": run.name,
        "latest_execution": execution_id,
        "mode": manifest.get("mode", ""),
        "status": manifest.get("status", "not_run"),
        **{stage: stages.get(stage, {}).get("status", "not_run")
           if isinstance(stages.get(stage, {}), dict) else "not_run" for stage in STAGES},
        "snr_bin_edges": json.dumps(edges, ensure_ascii=False) if edges else "",
        "results": results,
        "physics_roots": root_count(run / "Tracker/physics"),
        "noise_roots": root_count(run / "Tracker/noise"),
        "error": str(error).replace("\t", " ").replace("\n", " "),
    }


def write_status_index(batches_dir, output=None):
    batches_dir = Path(batches_dir).resolve()
    output = Path(output).resolve() if output else batches_dir / "run_status.tsv"
    rows = [row for run in sorted(batches_dir.iterdir()) if run.is_dir()
            if (row := run_row(run)) is not None]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        stream.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            stream.write("\t".join(str(row.get(column, "")).replace("\t", " ")
                                    for column in COLUMNS) + "\n")
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batches_dir", nargs="?", default="../batches")
    parser.add_argument("--output")
    args = parser.parse_args(argv)
    output = write_status_index(args.batches_dir, args.output)
    print(f"Run status index: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

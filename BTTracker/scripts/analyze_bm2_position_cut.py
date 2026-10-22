#!/usr/bin/env python3
"""Recompute BM2 timing resolutions after a BM2 prediction-position cut.

The source is each run's Combine/Pad-step3-rec.root. The existing timing reader
provides the single-matched-pad cohort, and this script joins its event IDs to
the BM2 prediction coordinates stored in TrackArrivalEvents.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
import sys
import uuid

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import analyze_time_resolution as timing


ACTIVE_DETECTORS = ("BM2", "Tracker1", "Tracker2", "Tracker3")
X_MIN_MM, X_MAX_MM = -10.0, 40.0
Y_MIN_MM, Y_MAX_MM = -40.0, 10.0


def read_position_events(source: Path) -> tuple[dict[int, tuple[float, float]], int]:
    root_file = timing.ROOT.TFile.Open(str(source), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"cannot open {source}")
    tree = root_file.Get("TrackArrivalEvents")
    if not tree:
        raise RuntimeError(f"TrackArrivalEvents is missing in {source}")
    for name in ("event", "bm_pred_x", "bm_pred_y"):
        if not tree.GetBranch(name):
            raise RuntimeError(f"required branch {name} is missing in {source}")

    result: dict[int, tuple[float, float]] = {}
    entries = int(tree.GetEntries())
    for index in range(entries):
        tree.GetEntry(index)
        event = int(tree.event)
        if event in result:
            raise RuntimeError(f"duplicate event ID {event} in {source}")
        try:
            x = float(tree.bm_pred_x[1])
            y = float(tree.bm_pred_y[1])
        except (TypeError, IndexError) as exc:
            raise RuntimeError("BM2 prediction branches do not have two detector slots") from exc
        result[event] = (x, y)
    root_file.Close()
    return result, entries


def analyze_run(run_dir: Path, output: Path) -> dict[str, object]:
    source = run_dir / "Combine" / "Pad-step3-rec.root"
    if not source.is_file():
        raise FileNotFoundError(source)
    all_positions, tree_entries = read_position_events(source)
    records, availability, timing_tree_entries = timing.read_events(source, 1, "pad")
    if tree_entries != timing_tree_entries:
        raise RuntimeError(f"entry count changed while reading {source}")
    if len(all_positions) != tree_entries:
        raise RuntimeError(f"event-position join is not one-to-one for {source}")

    roi_events = {
        event
        for event, (x, y) in all_positions.items()
        if math.isfinite(x) and math.isfinite(y)
        and X_MIN_MM <= x <= X_MAX_MM
        and Y_MIN_MM <= y <= Y_MAX_MM
    }
    selected = [record for record in records if int(record["event"]) in roi_events]
    if not selected:
        raise RuntimeError(f"no single-pad events pass the BM2 position cut in {run_dir}")

    output.mkdir(parents=True, exist_ok=False)
    reports: dict[str, object] = {}
    for method in timing.TIMING_METHODS:
        pairs, _raw, _core, weighted = timing.analyze(
            selected, "mean", ACTIVE_DETECTORS, method
        )
        method_output = output if method == "t20" else output / method
        fits = timing.plot_bm2_tracker_gaussians(
            selected,
            method_output,
            1,
            f"BM2 ROI; {method}",
            method,
        )
        reports[method] = {
            "pairwise": pairs,
            "weighted_core": weighted,
            "gaussian_fits": fits,
        }

    report: dict[str, object] = {
        "run": run_dir.name,
        "source": str(source),
        "tree_entries": tree_entries,
        "single_matched_bm2_pad_events_before_position_cut": len(records),
        "events_after_position_cut": len(selected),
        "position_cut": {
            "detector": "BM2",
            "x_mm_inclusive": [X_MIN_MM, X_MAX_MM],
            "y_mm_inclusive": [Y_MIN_MM, Y_MAX_MM],
        },
        "timing_selection": "BM2 bmS_n[1]==1; pair-specific valid timing cohorts",
        "detector_time_valid_counts_before_position_cut": availability,
        "methods": reports,
    }
    (output / "position_cut_time_resolution.json").write_text(
        json.dumps(timing.serializable(report), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (output / "README.md").write_text(
        "# BM2 position-cut time resolution\n\n"
        f"Run: `{run_dir.name}`\n\n"
        "Selection: one valid matched BM2 pad (`bmS_n[1] == 1`) and the BM2 "
        "track-prediction position inside inclusive x=[-10, 40] mm, "
        "y=[-40, 10] mm. Pairwise distributions use pair-specific valid "
        "timing cohorts, matching the main time-resolution workflow.\n\n"
        "The top-level PNG/PDF show the t20 BM2−Tracker1/2/3 distributions "
        "with Gaussian fits. Separate `t30`, `t40`, and `t50` folders contain "
        "the corresponding plots. `position_cut_time_resolution.json` records "
        "the weighted-core detector resolutions and fit details.\n",
        encoding="utf-8",
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace", type=Path, default=SCRIPT_DIR.parents[1],
        help="APV25 workspace root (defaults to this script's repository)",
    )
    parser.add_argument(
        "--runs", nargs="+", default=["006", "008", "009", "010", "011"],
        help="run numbers, default: 006 008 009 010 011",
    )
    parser.add_argument(
        "--execution-id",
        default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + "_bm2_position_cut_" + uuid.uuid4().hex[:8],
    )
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    for number in args.runs:
        run_dir = workspace / "batches" / f"run_{number}"
        output = (
            run_dir / "Combine" / "time_resolution" / "position_cut"
            / args.execution_id
        )
        report = analyze_run(run_dir, output)
        t20 = report["methods"]["t20"]["weighted_core"]
        bm2_index = ACTIVE_DETECTORS.index("BM2")
        resolution = t20["resolutions_ns"][bm2_index]
        error = t20["resolution_uncertainties_ns"][bm2_index]
        print(
            f"{run_dir.name}: position-cut N={report['events_after_position_cut']}; "
            f"BM2 t20 weighted-core resolution={resolution:.4f} +/- {error:.4f} ns; "
            f"output={output}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

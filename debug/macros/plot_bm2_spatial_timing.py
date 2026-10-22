#!/usr/bin/env python3
"""Map BM2 event counts/resolution and plot matched-pad amplitude vs delta-t20."""
from __future__ import annotations

import argparse
from array import array
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import statistics
import sys
import uuid

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "BTTracker/scripts"))

import ROOT

from analyze_time_resolution import DETECTORS, TIMING_METHODS, analyze, read_events
from analyze_bm_snr_timing import ACTIVE, join_events, read_pad_rows

PLANE_LOW_MM = -50.
PLANE_HIGH_MM = 50.
DEFAULT_CELL_SIZE_MM = 20.
MINIMUM_EVENTS = 40


def write_tsv(path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def quantile(values, fraction):
    ordered = sorted(values)
    pos = fraction * (len(ordered) - 1)
    low = int(pos)
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (ordered[high] - ordered[low]) * (pos - low)


def json_finite(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: json_finite(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_finite(item) for item in value]
    return value


def attach_bm2_position(run, rows):
    path = run / "Combine/Pad-step3-rec.root"
    source = ROOT.TFile.Open(str(path), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"cannot open {path}")
    try:
        tree = source.Get("TrackArrivalEvents")
        if not tree:
            raise RuntimeError(f"TrackArrivalEvents missing from {path}")
        required = ("event", "bm_pred_x", "bm_pred_y")
        for name in required:
            if not tree.GetBranch(name):
                raise RuntimeError(f"TrackArrivalEvents missing branch {name}")
        tree.SetBranchStatus("*", 0)
        for name in required:
            tree.SetBranchStatus(name, 1)
        positions = {}
        for entry in tree:
            event = int(entry.event)
            if event in positions:
                raise RuntimeError(f"duplicate TrackArrivalEvents event {event}")
            positions[event] = (float(entry.bm_pred_x[1]), float(entry.bm_pred_y[1]))
    finally:
        source.Close()
    for row in rows:
        event = int(row["event"])
        if event not in positions:
            raise RuntimeError(f"event {event} is missing its BM2 prediction position")
        row["x_mm"], row["y_mm"] = positions[event]
        row["delta_t20_ns"] = float(row["mean"][1]) - float(row["mean"][3])
    return rows


def bin_events(rows, cell_size_mm):
    bins_per_axis = round((PLANE_HIGH_MM - PLANE_LOW_MM) / cell_size_mm)
    edges = tuple(PLANE_LOW_MM + i * cell_size_mm for i in range(bins_per_axis + 1))
    groups = {(ix, iy): [] for ix in range(bins_per_axis) for iy in range(bins_per_axis)}
    outside = 0
    for row in rows:
        x, y = row["x_mm"], row["y_mm"]
        if not (PLANE_LOW_MM <= x < PLANE_HIGH_MM and
                PLANE_LOW_MM <= y < PLANE_HIGH_MM):
            outside += 1
            continue
        ix = min(bins_per_axis - 1, int((x - PLANE_LOW_MM) // cell_size_mm))
        iy = min(bins_per_axis - 1, int((y - PLANE_LOW_MM) // cell_size_mm))
        groups[ix, iy].append(row)
    return groups, outside, edges


def calculate_resolution(groups, minimum, edges):
    results = []
    for (ix, iy), selected in sorted(groups.items()):
        low_x, high_x = edges[ix:ix + 2]
        low_y, high_y = edges[iy:iy + 2]
        for method in TIMING_METHODS:
            row = dict(ix=ix, iy=iy, x_low_mm=low_x, x_high_mm=high_x,
                       y_low_mm=low_y, y_high_mm=high_y,
                       x_center_mm=(low_x + high_x) / 2,
                       y_center_mm=(low_y + high_y) / 2,
                       selected_events=len(selected), timing_method=method,
                       minimum_pair_core_events=0, bm2_sigma_ns=math.nan,
                       bm2_sigma_error_ns=math.nan, status="insufficient_events",
                       reason="")
            if len(selected) >= minimum:
                try:
                    pairs, _raw, _core, weighted = analyze(selected, "mean", ACTIVE, method)
                    row["minimum_pair_core_events"] = min(p["core_events"] for p in pairs)
                    if row["minimum_pair_core_events"] < minimum:
                        row["status"] = "insufficient_core_events"
                    elif any(not math.isfinite(v) or v <= 0 for v in
                             weighted["resolutions_ns"] + weighted["resolution_uncertainties_ns"]):
                        row["status"] = "invalid_variance_solution"
                    else:
                        row.update(status="ok", bm2_sigma_ns=weighted["resolutions_ns"][0],
                                   bm2_sigma_error_ns=weighted["resolution_uncertainties_ns"][0])
                except RuntimeError as exc:
                    row.update(status="fit_failed", reason=str(exc))
            results.append(row)
    return results


def draw_count_map(groups, resolution_rows, output, root_file, bins_per_axis,
                   low_edge, high_edge, cell_size_mm):
    canvas = ROOT.TCanvas("bm2_event_count_map", "BM2 selected event count", 850, 750)
    canvas.SetRightMargin(.23)
    hist = ROOT.TH2D("bm2_event_count", "Selected common timing events on BM2;X [mm];Y [mm]",
                     bins_per_axis, low_edge, high_edge,
                     bins_per_axis, low_edge, high_edge)
    hist.SetStats(0)
    hist.GetZaxis().SetTitle(f"Events / {cell_size_mm:g} mm #times {cell_size_mm:g} mm cell")
    for (ix, iy), selected in groups.items():
        hist.SetBinContent(ix + 1, iy + 1, len(selected))
    hist.Draw("COLZ")
    labels = []
    valid = {(int(r["ix"]), int(r["iy"])): r for r in resolution_rows
             if r["timing_method"] == "t20" and r["status"] == "ok"}
    for (ix, iy), selected in groups.items():
        row = valid.get((ix, iy))
        label_text = f"{row['bm2_sigma_ns']:.2f} ns" if row else "N/A"
        x = low_edge + (ix + .5) * cell_size_mm
        y = low_edge + (iy + .5) * cell_size_mm
        label = ROOT.TLatex(x, y, label_text)
        label.SetTextAlign(22)
        label.SetTextSize(min(.030, .19 / bins_per_axis))
        label.SetTextFont(62)
        label.Draw()
        labels.append(label)
    canvas.SaveAs(str(output / "bm2_event_count_map.png"))
    canvas.SaveAs(str(output / "bm2_event_count_map.pdf"))
    root_file.cd()
    hist.Write()
    canvas.Write()
    return canvas, hist, labels


def draw_resolution_map(rows, output, root_file, minimum, bins_per_axis,
                        low_edge, high_edge, cell_size_mm):
    canvas = ROOT.TCanvas("bm2_resolution_map", "BM2 timing resolution by position", 850, 750)
    canvas.SetRightMargin(.23)
    hist = ROOT.TH2D("bm2_resolution", "BM2 weighted core resolution;X [mm];Y [mm]",
                     bins_per_axis, low_edge, high_edge,
                     bins_per_axis, low_edge, high_edge)
    hist.SetStats(0)
    hist.GetZaxis().SetTitle("BM2 time resolution [ns]")
    valid = {(int(r["ix"]), int(r["iy"])): r for r in rows
             if r["timing_method"] == "t20" and r["status"] == "ok"}
    values = [r["bm2_sigma_ns"] for r in valid.values()]
    if values:
        hist.SetMinimum(min(values) * .9)
        hist.SetMaximum(max(values) * 1.1)
    hist.Draw("COLZ")
    labels = []
    for ix in range(bins_per_axis):
        for iy in range(bins_per_axis):
            x = low_edge + (ix + .5) * cell_size_mm
            y = low_edge + (iy + .5) * cell_size_mm
            row = valid.get((ix, iy))
            if row:
                label_text = (f"#splitline{{{row['bm2_sigma_ns']:.2f}#pm"
                              f"{row['bm2_sigma_error_ns']:.2f}}}"
                              f"{{N={row['selected_events']}}}")
            else:
                count = next(r["selected_events"] for r in rows
                             if r["ix"] == ix and r["iy"] == iy and r["timing_method"] == "t20")
                label_text = f"#splitline{{N={count}}}{{<{minimum}}}"
                mask = ROOT.TBox(low_edge + cell_size_mm * ix,
                                 low_edge + cell_size_mm * iy,
                                 low_edge + cell_size_mm * (ix + 1),
                                 low_edge + cell_size_mm * (iy + 1))
                mask.SetFillColor(ROOT.kWhite)
                mask.SetFillStyle(1001)
                mask.SetLineColor(ROOT.kGray + 1)
                mask.Draw()
                labels.append(mask)
            label = ROOT.TLatex(x, y, label_text)
            label.SetTextAlign(22)
            label.SetTextSize(min(.020, .16 / bins_per_axis))
            label.Draw()
            labels.append(label)
    canvas.SaveAs(str(output / "bm2_resolution_map_t20.png"))
    canvas.SaveAs(str(output / "bm2_resolution_map_t20.pdf"))
    root_file.cd()
    hist.Write()
    canvas.Write()
    return canvas, hist, labels


def draw_amplitude_dt(rows, output, root_file, max_amplitude):
    differences = [float(row["delta_t20_ns"]) for row in rows]
    # Use a robust central range for legibility and report excluded tails.
    low = quantile(differences, .005)
    high = quantile(differences, .995)
    padding = max((high - low) * .06, .1)
    low -= padding
    high += padding
    canvas = ROOT.TCanvas("dt20_vs_maxA", "BM2 - Tracker2 t20 versus peak amplitude", 1000, 800)
    canvas.SetRightMargin(.16)
    hist = ROOT.TH2D("dt20_vs_maxA",
                     f"BM2 matched single-pad sample;Peak amplitude maxA [ADC];"
                     f"#Delta t_{{20}} = t_{{BM2}} - t_{{Tracker2}} [ns]",
                     120, 0., max_amplitude, 120, low, high)
    hist.SetStats(0)
    hist.GetZaxis().SetTitle("Events / bin")
    in_display = 0
    for row in rows:
        amplitude = float(row["maxA_adc"])
        difference = float(row["delta_t20_ns"])
        if math.isfinite(amplitude) and math.isfinite(difference):
            hist.Fill(amplitude, difference)
            in_display += int(0. <= amplitude < max_amplitude and low <= difference <= high)
    hist.Draw("COLZ")
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextSize(.035)
    label.DrawLatex(.15, .92,
                    f"N in displayed range = {in_display:,} / {len(rows):,}; BM2 single matched pad")
    objects = [hist, label]
    canvas.SaveAs(str(output / "dt20_vs_maxA.png"))
    canvas.SaveAs(str(output / "dt20_vs_maxA.pdf"))
    root_file.cd()
    hist.Write()
    canvas.Write()
    return canvas, hist, objects, (low, high), in_display


def set_heatmap_palette():
    """Use the blue-cyan-green-yellow scale from the reference, without red."""
    stops = array("d", [0.00, 0.28, 0.53, 0.77, 1.00])
    red = array("d", [0.04, 0.00, 0.00, 0.18, 1.00])
    green = array("d", [0.08, 0.18, 0.78, 0.92, 0.96])
    blue = array("d", [0.48, 0.92, 0.90, 0.28, 0.02])
    ncolors = 255
    first = ROOT.TColor.CreateGradientColorTable(
        len(stops), stops, red, green, blue, ncolors)
    ROOT.gStyle.SetPalette(ncolors, array("i", range(first, first + ncolors)))
    ROOT.gStyle.SetNumberContours(ncolors)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run")
    parser.add_argument("--min-events", type=int, default=MINIMUM_EVENTS)
    parser.add_argument("--maxA-max", type=float, default=2000.)
    parser.add_argument("--cell-size-mm", type=float, default=DEFAULT_CELL_SIZE_MM)
    args = parser.parse_args()
    bins_per_axis = ((PLANE_HIGH_MM - PLANE_LOW_MM) / args.cell_size_mm
                     if math.isfinite(args.cell_size_mm) and args.cell_size_mm > 0 else 0)
    if (args.min_events < 3 or not math.isfinite(args.maxA_max) or args.maxA_max <= 0 or
            bins_per_axis < 1 or not math.isclose(bins_per_axis, round(bins_per_axis),
                                                 rel_tol=0, abs_tol=1e-9)):
        parser.error("min-events must be >= 3, maxA-max and cell-size-mm must be positive, "
                     "and cell-size-mm must divide the 100 mm plane evenly")
    bins_per_axis = round(bins_per_axis)
    edges = tuple(PLANE_LOW_MM + i * args.cell_size_mm for i in range(bins_per_axis + 1))
    run = Path(args.run).resolve(strict=True)
    required_inputs = (run / "Combine/Tracker-step2-ana.root",
                       run / "Combine/Pad-step3-rec.root")
    missing = [str(path) for path in required_inputs if not path.is_file()]
    if missing:
        parser.error("required run input ROOT files are missing: " + ", ".join(missing))

    _all_pads, best = read_pad_rows(run)
    joined, _availability, _total = join_events(run, best)
    all_times, _availability, _total = read_events(run / "Combine/Pad-step3-rec.root", 1)
    timing_by_event = {row["event"]: row for row in all_times}
    for row in joined:
        timing = timing_by_event.get(row["event"])
        if timing is None:
            raise RuntimeError(f"event {row['event']} is missing from TrackArrivalEvents")
        for method in ("t30", "t40", "t50"):
            row["mean_" + method] = timing["mean_" + method]
            row["min_" + method] = timing["min_" + method]
    timing_keys = ("mean", "mean_t30", "mean_t40", "mean_t50")
    indices = [DETECTORS.index(name) for name in ACTIVE]
    selected = [row for row in joined if row["all_four_valid"] and all(
        math.isfinite(row[key][i]) for key in timing_keys for i in indices)]
    if len(selected) < 3:
        raise RuntimeError("fewer than three common events with valid four-detector timing")
    selected = attach_bm2_position(run, selected)
    groups, outside, _edges = bin_events(selected, args.cell_size_mm)
    resolution_rows = calculate_resolution(groups, args.min_events, edges)
    differences = [float(row["delta_t20_ns"]) for row in selected]
    dt_low, dt_high = quantile(differences, .005), quantile(differences, .995)

    execution = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    output = run / "Combine/time_resolution/spatial_binned" / execution
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "spatial_resolution.tsv", resolution_rows)
    event_rows = [dict(event=r["event"], x_bm2_mm=r["x_mm"], y_bm2_mm=r["y_mm"],
                       cell_x=int((r["x_mm"] - PLANE_LOW_MM) // args.cell_size_mm),
                       cell_y=int((r["y_mm"] - PLANE_LOW_MM) // args.cell_size_mm), maxA_adc=r["maxA"],
                       bm2_t20_ns=r["mean"][1], tracker2_t20_ns=r["mean"][3],
                       delta_t20_ns=r["delta_t20_ns"])
                  for r in selected if -50 <= r["x_mm"] < 50 and -50 <= r["y_mm"] < 50]
    write_tsv(output / "selected_events.tsv", event_rows)
    root_file = ROOT.TFile(str(output / "bm2_spatial_diagnostics.root"), "RECREATE")
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    set_heatmap_palette()
    keepalive = [draw_count_map(groups, resolution_rows, output, root_file,
                                bins_per_axis, PLANE_LOW_MM, PLANE_HIGH_MM,
                                args.cell_size_mm),
                 draw_resolution_map(resolution_rows, output, root_file, args.min_events,
                                     bins_per_axis, PLANE_LOW_MM, PLANE_HIGH_MM,
                                     args.cell_size_mm),
                 draw_amplitude_dt(event_rows, output, root_file, args.maxA_max)]
    root_file.Close()
    dt_display = keepalive[2][4]
    report = dict(run=str(run), grid=f"{bins_per_axis}x{bins_per_axis} cells; "
                  f"{args.cell_size_mm:g} mm pitch; x,y in "
                  f"[{PLANE_LOW_MM:g},{PLANE_HIGH_MM:g}) mm about BM2 centre",
                  position="TrackArrivalEvents.bm_pred_x/y[1]", timing_source="Pad-step3-rec.root/TrackArrivalEvents",
                  amplitude_source="Tracker-step2-ana.root/BMPadSignals; matched chip 13 pad in BM2 3x3 window",
                  bm2_selection="exactly one valid sigmoid-timed pad in matching window",
                  sample="common valid BM2,Tracker1,Tracker2,Tracker3 mean times for t20,t30,t40,t50",
                  count_map_events=len(selected), events_outside_spatial_area=outside,
                  dt_definition="mean_t20(BM2)-mean_t20(Tracker2), ns",
                  dt_amplitude_events=len(event_rows), maxA_axis_adc=[0, args.maxA_max],
                  dt_display_quantiles=[.005, .995], dt_display_range_ns=[dt_low, dt_high],
                  dt_events_inside_display_range=dt_display,
                  dt_events_outside_display_range=len(event_rows) - dt_display,
                  min_events_per_cell_and_pair_core=args.min_events,
                  resolution_rows=resolution_rows,
                  caveat="approximate weighted pairwise statistical errors; covariance and systematics omitted")
    (output / "analysis.json").write_text(
        json.dumps(json_finite(report), indent=2, allow_nan=False) + "\n")
    print(f"Selected common events: {len(selected)}; inside 100x100 mm: {len(event_rows)}")
    print(f"Output: {output}")


if __name__ == "__main__":
    main()

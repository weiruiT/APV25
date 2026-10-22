#!/usr/bin/env python3
"""Read existing run ROOT files and plot BM2 timing resolution versus track angle."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import statistics
import uuid

from analyze_time_resolution import ROOT, DETECTORS, TIMING_METHODS, analyze, read_events

ACTIVE = ("BM2", "Tracker1", "Tracker2", "Tracker3")


def table(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def join_angles(records, track_path, pad_path):
    """Join selected events by DAQ ID, checking timestamp and duplicate IDs."""
    wanted = {row["event"] for row in records}
    if len(wanted) != len(records):
        raise RuntimeError("duplicate event IDs in the selected timing sample")
    angles = {}
    root = ROOT.TFile.Open(str(track_path), "READ")
    if not root or root.IsZombie():
        raise RuntimeError(f"cannot open {track_path}")
    try:
        tree = root.Get("Track")
        if not tree:
            raise RuntimeError("Step4 is missing the Track tree")
        tree.SetBranchStatus("*", 0)
        for name in ("event", "timestamp", "useflag", "XfitPars", "YfitPars"):
            if not tree.GetBranch(name):
                raise RuntimeError(f"Track is missing {name}")
            tree.SetBranchStatus(name, 1)
        for entry in tree:
            event = int(entry.event)
            if event not in wanted or not bool(entry.useflag):
                continue
            if event in angles:
                raise RuntimeError(f"duplicate selected Track event {event}")
            kx, ky = float(entry.XfitPars[0]), float(entry.YfitPars[0])
            if not math.isfinite(kx) or not math.isfinite(ky):
                raise RuntimeError(f"nonfinite slope for event {event}")
            angles[event] = (math.degrees(math.atan(math.hypot(kx, ky))),
                             math.degrees(math.atan(kx)), math.degrees(math.atan(ky)),
                             int(entry.timestamp))
    finally:
        root.Close()
    missing = wanted - angles.keys()
    if missing:
        raise RuntimeError(f"{len(missing)} timing events lack a selected Step4 track; check the input pair")
    root = ROOT.TFile.Open(str(pad_path), "READ")
    try:
        tree = root.Get("TrackArrivalEvents")
        if not tree.GetBranch("timestamp"):
            raise RuntimeError("TrackArrivalEvents lacks timestamp for join verification")
        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("event", 1)
        tree.SetBranchStatus("timestamp", 1)
        for entry in tree:
            event = int(entry.event)
            if event in angles and int(entry.timestamp) != angles[event][3]:
                raise RuntimeError(f"Step4/Pad timestamp mismatch for event {event}")
    finally:
        root.Close()
    for row in records:
        row["theta_deg"], row["theta_x_deg"], row["theta_y_deg"], _ = angles[row["event"]]
    return records


def angle_edges(values, count):
    """Approximate equal populations; never split identical angles across bins."""
    values = sorted(values)
    if values[0] == values[-1]:
        raise ValueError("all selected angles are identical; cannot form angle intervals")
    edges = [values[0]]
    for index in range(1, count):
        position = len(values) * index // count
        low, high = values[position - 1], values[position]
        boundary = (low + high) / 2 if low < high else high
        if edges[-1] < boundary < values[-1]:
            edges.append(boundary)
    return edges + [values[-1]]


def calculate(records, edges, minimum):
    rows, details = [], []
    for index, (low, high) in enumerate(zip(edges, edges[1:])):
        selected = [row for row in records if low <= row["theta_deg"] and
                    (row["theta_deg"] < high or
                     (index == len(edges) - 2 and row["theta_deg"] == high))]
        mean_angle = statistics.fmean(row["theta_deg"] for row in selected) if selected else math.nan
        for method in TIMING_METHODS:
            row = dict(bin=index + 1, theta_low_deg=low, theta_high_deg=high,
                       theta_mean_deg=mean_angle,
                       theta_median_deg=statistics.median([r["theta_deg"] for r in selected])
                       if selected else (low + high) / 2,
                       selected_events=len(selected), timing_method=method,
                       minimum_core_events=0, bm2_sigma_ns=math.nan,
                       bm2_sigma_error_ns=math.nan, chi2=math.nan, ndf=2,
                       status="insufficient_events", reason="")
            if len(selected) >= minimum:
                try:
                    pairs, raw, core, weighted = analyze(selected, "mean", ACTIVE, method)
                    details.append(dict(bin=index + 1, timing_method=method, pairs=pairs,
                                        raw=raw, core=core, weighted_core=weighted))
                    row["minimum_core_events"] = min(pair["core_events"] for pair in pairs)
                    if row["minimum_core_events"] < minimum:
                        row["status"] = "insufficient_core_events"
                    elif any(not math.isfinite(value) or value <= 0 for value in
                             weighted["resolutions_ns"] + weighted["resolution_uncertainties_ns"]):
                        row["status"] = "invalid_variance_solution"
                    else:
                        row.update(status="ok", bm2_sigma_ns=weighted["resolutions_ns"][0],
                                   bm2_sigma_error_ns=weighted["resolution_uncertainties_ns"][0],
                                   chi2=weighted["chi2"], ndf=weighted["ndf"])
                except RuntimeError as exc:
                    row.update(status="fit_failed", reason=str(exc))
            rows.append(row)
    return rows, details


def draw(records, edges, rows, output, title, selected_methods):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    root = ROOT.TFile(str(output / "timing_vs_angle.root"), "RECREATE")
    canvas = ROOT.TCanvas("timing_vs_angle", title, 1100, 1000)
    canvas.Divide(1, 2)
    canvas.cd(1)
    low, high = edges[0], edges[-1]
    upper = math.nextafter(high, math.inf)
    hist = ROOT.TH1D("selected_track_theta",
                     title + ";Track angle to z [deg];Selected events / bin",
                     40, low, upper)
    for row in records:
        hist.Fill(row["theta_deg"])
    hist.SetLineColor(ROOT.kBlue + 1)
    hist.SetLineWidth(2)
    peak = max(1.0, hist.GetMaximum())
    hist.SetMaximum(peak * 1.28)
    hist.Draw("HIST")
    objects = [hist]
    for edge in edges[1:-1]:
        line = ROOT.TLine(edge, 0, edge, peak * 1.18)
        line.SetLineStyle(2)
        line.Draw()
        objects.append(line)
    # All CFD methods share the same pre-clipping sample; label each bin once.
    for row in rows:
        if row["timing_method"] != "t20":
            continue
        x = (row["theta_low_deg"] + row["theta_high_deg"]) / 2
        label = ROOT.TLatex(x, peak * 1.12, f"N = {row['selected_events']}")
        label.SetTextAlign(22)
        label.SetTextSize(.032)
        label.Draw()
        objects.append(label)
    canvas.cd(2).SetGrid()
    good = [row for row in rows if row["status"] == "ok" and
            row["timing_method"] in selected_methods]
    if good:
        lower = min(row["bm2_sigma_ns"] - row["bm2_sigma_error_ns"] for row in good)
        upper_error = max(row["bm2_sigma_ns"] + row["bm2_sigma_error_ns"] for row in good)
        span = max(upper_error - lower, 0.1)
        ymin = max(0.0, lower - .25 * span)
        ymax = upper_error + .25 * span
    else:
        ymin, ymax = 0.0, 1.0
    frame = ROOT.TH1D("resolution_axes",
                      "BM2: weighted 3-sigma core solution;Track angle to z [deg];BM2 time resolution [ns]",
                      100, low, upper)
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymax)
    frame.Draw("AXIS")
    objects.append(frame)
    legend = ROOT.TLegend(.70, .67, .90, .89)
    legend.SetBorderSize(0)
    colors = {"t20": ROOT.kBlue + 1, "t30": ROOT.kRed + 1,
              "t40": ROOT.kGreen + 2, "t50": ROOT.kMagenta + 1}
    for method in selected_methods:
        color = colors[method]
        graph = ROOT.TGraphAsymmErrors()
        graph.SetName("bm2_resolution_" + method)
        graph.SetMarkerStyle(20 + TIMING_METHODS.index(method))
        graph.SetMarkerColor(color)
        graph.SetLineColor(color)
        for i, row in enumerate(r for r in good if r["timing_method"] == method):
            x = row["theta_mean_deg"]
            graph.SetPoint(i, x, row["bm2_sigma_ns"])
            graph.SetPointError(i, x - row["theta_low_deg"], row["theta_high_deg"] - x,
                               row["bm2_sigma_error_ns"], row["bm2_sigma_error_ns"])
        if graph.GetN():
            graph.Draw("P SAME")
            legend.AddEntry(graph, method, "p")
        graph.Write()
        objects.append(graph)
    legend.Draw()
    objects.append(legend)
    for row in rows:
        if row["timing_method"] != "t20" or not math.isfinite(row["theta_mean_deg"]):
            continue
        x = row["theta_mean_deg"]
        label = ROOT.TLatex(x, ymin + .075 * (ymax - ymin),
                           "#splitline{#LT#theta#GT}{" + f"{x:.2f}" + "^{#circ}}")
        label.SetTextAlign(22)
        label.SetTextSize(.028)
        label.Draw()
        objects.append(label)
    if not good:
        label = ROOT.TLatex()
        label.SetNDC()
        label.DrawLatex(.15, .50, "No valid bins; see angle_bins.tsv")
        objects.append(label)
    canvas.SaveAs(str(output / "timing_vs_angle.png"))
    canvas.SaveAs(str(output / "timing_vs_angle.pdf"))
    hist.Write()
    canvas.Write()
    root.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run")
    parser.add_argument("--bins", type=int, default=5)
    parser.add_argument("--angle-edges", help="fixed degree edges, e.g. 0,3,6,9,12,18")
    parser.add_argument("--min-events", type=int, default=40)
    parser.add_argument("--bm2-matched-pad-count", type=int)
    parser.add_argument("--timing-method", type=int, choices=range(5), default=0,
                        help="0: draw all CFD methods; 1-4: draw t20, t30, t40, or t50 only")
    args = parser.parse_args()
    if args.bins < 1 or args.min_events < 3:
        parser.error("bins must be >= 1; min-events must be >= 3")
    run = Path(args.run).resolve(strict=True)
    pad_path = run / "Combine/Pad-step3-rec.root"
    track_path = run / "Combine/Tracker-step4-rec.root"
    report_path = run / "Combine/time_resolution/time_resolution.json"
    count = args.bm2_matched_pad_count
    selection_source = "command line"
    if count is None:
        count, selection_source = 1, "default: BM2 single matched pad"
        if report_path.is_file():
            report = json.loads(report_path.read_text())
            selection = report.get("timing_selection", {})
            if "bm2_matched_pad_count" in selection:
                count = selection["bm2_matched_pad_count"]
                selection_source = str(report_path)
    if type(count) is not int or count < 0:
        parser.error("BM2 matched pad count must be an integer >= 0")
    records, _, total = read_events(pad_path, count)
    indices = [DETECTORS.index(name) for name in ACTIVE]
    # The top histogram and all four CFD curves use exactly the same input events.
    records = [row for row in records if all(math.isfinite(row[key][i])
               for key in ("mean", "mean_t30", "mean_t40", "mean_t50") for i in indices)]
    if len(records) < 3:
        raise RuntimeError("fewer than three common BM2/Tracker1/2/3 timing events")
    records = join_angles(records, track_path, pad_path)
    if args.angle_edges:
        edges = [float(value) for value in args.angle_edges.split(",")]
        if len(edges) < 2 or any(not math.isfinite(v) or v < 0 or v > 90 for v in edges) or any(
                b <= a for a, b in zip(edges, edges[1:])):
            parser.error("angle edges must strictly increase within [0,90] degrees")
    else:
        edges = angle_edges([row["theta_deg"] for row in records], args.bins)
    outside = sum(not edges[0] <= row["theta_deg"] <= edges[-1] for row in records)
    records = [row for row in records if edges[0] <= row["theta_deg"] <= edges[-1]]
    if not records:
        raise RuntimeError("no selected events inside the angle range")
    rows, details = calculate(records, edges, args.min_events)
    execution = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    output = run / "Combine/time_resolution/angle_binned" / execution
    output.mkdir(parents=True, exist_ok=False)
    table(output / "angle_bins.tsv", rows)
    event_rows = []
    for row in records:
        event_row = {key: row[key] for key in ("event", "theta_deg", "theta_x_deg", "theta_y_deg")}
        for method in TIMING_METHODS:
            key = "mean" if method == "t20" else "mean_" + method
            event_row.update({name + "_" + method + "_ns": row[key][i] for name, i in zip(ACTIVE, indices)})
        event_rows.append(event_row)
    table(output / "selected_events.tsv", event_rows)
    provenance = dict(run=str(run), inputs=[dict(path=str(p), bytes=p.stat().st_size,
                     mtime_ns=p.stat().st_mtime_ns) for p in (pad_path, track_path)],
                     total_arrival_entries=total, common_selected_events=len(records),
                     outside_angle_range=outside, active_detectors=ACTIVE,
                     bm2_matched_pad_count=count, selection_source=selection_source,
                     angle_definition="atan(hypot(XfitPars[0],YfitPars[0])) in degrees; angle to z",
                     angle_edges_deg=edges, last_interval_includes_upper_edge=True,
                     graph_x="arithmetic mean of the common selected angles in each bin, before timing clipping",
                     minimum_core_events_per_pair=args.min_events,
                     sample="useflag; BM2 pad count cut; all four detectors valid for all four CFD times",
                     estimator="mean time, iterative 3-sigma core, weighted pairwise variance solution",
                     uncertainties="approximate propagated statistical errors; pair covariances are not included",
                     assumptions="independent detector timing errors after common trigger cancellation",
                     solutions=details)
    def finite_json(value):
        if isinstance(value, float) and not math.isfinite(value):
            return None
        if isinstance(value, dict):
            return {k: finite_json(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [finite_json(v) for v in value]
        return value
    selected_methods = TIMING_METHODS if args.timing_method == 0 else (TIMING_METHODS[args.timing_method - 1],)
    provenance["drawn_timing_methods"] = selected_methods
    (output / "analysis.json").write_text(json.dumps(finite_json(provenance), indent=2) + "\n")
    draw(records, edges, rows, output,
         run.name + ": common BM2 + Tracker1/2/3 timing sample", selected_methods)
    print(f"Selected {len(records)} events; output: {output}")
    print(f"Valid bin/method results: {sum(row['status'] == 'ok' for row in rows)}/{len(rows)}")


if __name__ == "__main__":
    main()

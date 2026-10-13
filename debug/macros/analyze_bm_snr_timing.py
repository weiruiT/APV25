#!/usr/bin/env python3
"""Read-only BM2 SNR/timing diagnostic; needs PyROOT and Python stdlib.

Join the actual single matched BM2 pad to TrackArrivalEvents, then reuse
the production pair-specific mean-t20/core variance estimator in each SNR
bin and above each threshold. Inputs and production outputs are never changed.
"""
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

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "BTTracker/scripts"))
import ROOT
from analyze_time_resolution import analyze, read_events, valid_time

ACTIVE = ("BM2", "Tracker1", "Tracker2", "Tracker3")
BINS = (5., 15., 25., 40., 60., 90., 130., math.inf)
COLORS = (ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1)


def quantile(values, fraction):
    values = sorted(v for v in values if math.isfinite(v))
    if not values:
        return math.nan
    index = fraction * (len(values) - 1)
    lo = int(index)
    hi = min(lo + 1, len(values) - 1)
    return values[lo] + (values[hi] - values[lo]) * (index - lo)


def describe(values):
    values = [v for v in values if math.isfinite(v)]
    return {"n": len(values), "mean": statistics.fmean(values) if values else math.nan,
            "median": quantile(values, .5), "q16": quantile(values, .16),
            "q84": quantile(values, .84)}


def clean(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    return value


def write_tsv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_pad_rows(run):
    source = run / "Combine/Tracker-step2-ana.root"
    f = ROOT.TFile.Open(str(source))
    if not f or f.IsZombie():
        raise RuntimeError(f"cannot open {source}")
    t = f.Get("BMPadSignals")
    if not t:
        raise RuntimeError(f"BMPadSignals missing: {source}")
    integers = ("event", "raw_chip", "position", "global_channel", "sigmoidFitStatus")
    doubles = ("maxA", "peak_noise_rms", "arrivalT_sigmoid", "sigmoid_H",
               "sigmoid_k", "sigmoidFitChi2NDF", "riseT", "TOT", "pedrms")
    buffers = {name: array("i" if name in integers else "d", [0])
               for name in integers + doubles}
    t.SetBranchStatus("*", 0)
    for name, buffer in buffers.items():
        if not t.GetBranch(name):
            raise RuntimeError(f"{source}: missing {name}")
        t.SetBranchStatus(name, 1)
        t.SetBranchAddress(name, buffer)
    all_rows, best = [], {}
    for entry in range(int(t.GetEntries())):
        t.GetEntry(entry)
        if buffers["raw_chip"][0] != 13:
            continue
        row = {name: buffer[0] for name, buffer in buffers.items()}
        row["step2_entry"] = entry
        sigma, height, k = row["peak_noise_rms"], row["sigmoid_H"], row["sigmoid_k"]
        row["snr"] = row["maxA"] / sigma if sigma > 0 else math.nan
        row["slope20_adc_ns"] = .16 * height / (25. * k) if height > 0 and k > 0 else math.nan
        row["jitter20_proxy_ns"] = sigma / row["slope20_adc_ns"] if row["slope20_adc_ns"] > 0 else math.nan
        row["sigmoid_rise10_90_ns"] = 25. * k * math.log(81.) if k > 0 else math.nan
        row["fit_residual_over_noise"] = math.sqrt(max(0., row["sigmoidFitChi2NDF"])) / sigma if sigma > 0 and row["sigmoidFitChi2NDF"] >= 0 else math.nan
        all_rows.append(row)
        key = (row["event"], row["position"])
        if key not in best or row["maxA"] >= best[key]["maxA"]:
            best[key] = row
    t.ResetBranchAddresses()
    f.Close()
    return all_rows, best


def join_events(run, best):
    records, availability, total = read_events(run / "Combine/Pad-step3-rec.root", 1)
    by_event = {r["event"]: r for r in records}
    if len(by_event) != len(records):
        raise RuntimeError(f"{run}: duplicate TrackArrivalEvents event IDs")
    f = ROOT.TFile.Open(str(run / "Combine/Pad-step3-rec.root"))
    t = f.Get("TrackArrivalEvents")
    event = array("i", [0]); ix = array("i", [0, 0]); iy = array("i", [0, 0])
    for name, buffer in (("event", event), ("bm_pred_ix", ix), ("bm_pred_iy", iy)):
        t.SetBranchAddress(name, buffer)
    joined = []
    for entry in range(int(t.GetEntries())):
        t.GetEntry(entry)
        record = by_event.get(event[0])
        if record is None:
            continue
        matches = []
        for x in range(max(0, ix[1] - 1), min(11, ix[1] + 2)):
            for y in range(max(0, iy[1] - 1), min(11, iy[1] + 2)):
                row = best.get((event[0], x + 11 * y))
                if row and valid_time(row["arrivalT_sigmoid"]):
                    matches.append(row)
        if len(matches) != 1:
            raise RuntimeError(f"{run}: event {event[0]} has {len(matches)} joined pads, expected 1")
        row = dict(matches[0])
        if abs(25. * row["arrivalT_sigmoid"] - record["mean"][1]) > 1.e-7:
            raise RuntimeError(f"{run}: event {event[0]} timing mismatch")
        row["mean"] = record["mean"]
        row["min"] = record["min"]
        row["all_four_valid"] = all(math.isfinite(v) for v in record["mean"][1:])
        row["run"] = run.name
        joined.append(row)
    t.ResetBranchAddresses(); f.Close()
    return joined, availability, total


def timing_summary(records):
    result = {"selected_events": len(records), "bm2_sigma_ns": math.nan,
              "bm2_sigma_error_ns": math.nan, "model_chi2": math.nan,
              "model_ndf": 2, "minimum_pair_events": 0, "error": ""}
    try:
        pairs, _raw, _core, weighted = analyze(records, "mean", ACTIVE)
        result.update(bm2_sigma_ns=weighted["resolutions_ns"][0],
                      bm2_sigma_error_ns=weighted["resolution_uncertainties_ns"][0],
                      model_chi2=weighted["chi2"],
                      minimum_pair_events=min(r["events"] for r in pairs))
        for row in pairs:
            key = row["first"] + "_" + row["second"]
            result[key + "_n"] = row["events"]
            result[key + "_core_sigma_ns"] = row["core_sigma_ns"]
            result[key + "_core_mu_ns"] = row["core_mean_ns"]
    except RuntimeError as exc:
        result["error"] = str(exc)
    return result


def bin_rows(run_name, events):
    rows = []
    for low, high in zip(BINS, BINS[1:]):
        selected = [r for r in events if low <= r["snr"] < high]
        if len(selected) < 40:
            continue
        row = {"run": run_name, "snr_low": low, "snr_high": high,
               "snr_median": quantile([r["snr"] for r in selected], .5),
               "amplitude_median_adc": quantile([r["maxA"] for r in selected], .5),
               "noise_median_adc": quantile([r["peak_noise_rms"] for r in selected], .5),
               "jitter20_proxy_median_ns": quantile([r["jitter20_proxy_ns"] for r in selected], .5),
               "rise10_90_median_ns": quantile([r["sigmoid_rise10_90_ns"] for r in selected], .5),
               "fit_residual_over_noise_median": quantile([r["fit_residual_over_noise"] for r in selected], .5)}
        row.update(timing_summary(selected))
        common = timing_summary([r for r in selected if r["all_four_valid"]])
        row["common_four_events"] = common["selected_events"]
        row["common_four_bm2_sigma_ns"] = common["bm2_sigma_ns"]
        row["common_four_bm2_sigma_error_ns"] = common["bm2_sigma_error_ns"]
        rows.append(row)
    return rows


def add_hist(pad, name, title, cohorts, field, bins, low, high, objects):
    pad.cd()
    legend = ROOT.TLegend(.59, .70, .88, .88)
    legend.SetBorderSize(0); legend.SetFillStyle(0)
    histograms = []
    for index, (run, records) in enumerate(cohorts.items()):
        hist = ROOT.TH1D(name + "_" + run, title, bins, low, high)
        for record in records:
            value = record[field]
            if math.isfinite(value):
                hist.Fill(value)
        if hist.GetEntries():
            hist.Scale(1. / hist.GetEntries())
        hist.SetLineColor(COLORS[index % len(COLORS)]); hist.SetLineWidth(2)
        legend.AddEntry(hist, run, "l"); histograms.append(hist)
    maximum = max(h.GetMaximum() for h in histograms)
    for index, hist in enumerate(histograms):
        hist.SetMaximum(maximum * 1.23)
        hist.Draw("HIST" if index == 0 else "HIST SAME")
    legend.Draw(); objects.extend(histograms + [legend])


def add_graph(pad, name, title, rows, xfield, yfield, errorfield, objects):
    pad.cd(); pad.SetGrid()
    multigraph = ROOT.TMultiGraph(name, title)
    legend = ROOT.TLegend(.62, .69, .89, .88)
    legend.SetBorderSize(0); legend.SetFillStyle(0)
    runs = list(dict.fromkeys(r["run"] for r in rows))
    for index, run in enumerate(runs):
        selected = [r for r in rows if r["run"] == run and math.isfinite(r.get(yfield, math.nan))]
        graph = ROOT.TGraphErrors(len(selected))
        graph.SetName(name + "_" + run)
        graph.SetLineColor(COLORS[index % len(COLORS)]); graph.SetMarkerColor(COLORS[index % len(COLORS)])
        graph.SetMarkerStyle(20 + index); graph.SetLineWidth(2)
        for i, row in enumerate(selected):
            graph.SetPoint(i, row[xfield], row[yfield])
            error = row.get(errorfield, 0.) if errorfield else 0.
            graph.SetPointError(i, 0., error if math.isfinite(error) else 0.)
        multigraph.Add(graph, "LP"); legend.AddEntry(graph, run, "lp"); objects.append(graph)
    multigraph.Draw("A"); legend.Draw(); objects.extend([multigraph, legend])


def make_plots(output, all_pads, joined, bins, thresholds):
    objects = []
    overview = ROOT.TCanvas("overview", "BM2 SNR and timing", 1500, 1120)
    overview.Divide(2, 2)
    add_hist(overview.cd(1), "amplitude", "BM2 Step2 waveforms;maxA [ADC];Fraction / bin", all_pads, "maxA", 70, 0, 1750, objects)
    add_hist(overview.cd(2), "noise", "BM2 matched single-pad timing cohort;Channel noise #sigma [ADC];Fraction / bin", joined, "peak_noise_rms", 50, 0, 25, objects)
    add_hist(overview.cd(3), "snr", "BM2 matched single-pad timing cohort;SNR = maxA / peak_noise_rms;Fraction / bin", joined, "snr", 60, 0, 180, objects)
    add_graph(overview.cd(4), "resolution_snr", "Same production variance estimator in each SNR bin;Median SNR;BM2 core resolution [ns]", bins, "snr_median", "bm2_sigma_ns", "bm2_sigma_error_ns", objects)
    overview.SaveAs(str(output / "snr_timing_comparison.png")); overview.SaveAs(str(output / "snr_timing_comparison.pdf"))
    mechanisms = ROOT.TCanvas("mechanisms", "SNR, shape, jitter and time walk", 1800, 600)
    mechanisms.Divide(3, 1)
    add_graph(mechanisms.cd(1), "rise_snr", "Fitted rise duration;Median SNR;Median sigmoid t90-t10 [ns]", bins, "snr_median", "rise10_90_median_ns", None, objects)
    add_graph(mechanisms.cd(2), "jitter_snr", "Noise / slope proxy (not total fitted-time error);Median SNR;Median #sigma_{noise} / slope_{20} [ns]", bins, "snr_median", "jitter20_proxy_median_ns", None, objects)
    add_graph(mechanisms.cd(3), "walk_snr", "Amplitude-dependent mean shift;Median SNR;Core mean t_{BM2} - t_{Tracker2} [ns]", bins, "snr_median", "BM2_Tracker2_core_mu_ns", None, objects)
    mechanisms.SaveAs(str(output / "snr_timing_mechanisms.png")); mechanisms.SaveAs(str(output / "snr_timing_mechanisms.pdf"))
    cuts = ROOT.TCanvas("cuts", "SNR threshold scan", 1500, 620)
    cuts.Divide(2, 1)
    add_graph(cuts.cd(1), "sigma_cut", "SNR threshold scan;Minimum SNR;BM2 core resolution [ns]", thresholds, "snr_min", "bm2_sigma_ns", "bm2_sigma_error_ns", objects)
    add_graph(cuts.cd(2), "fraction_cut", "Fraction of original timing cohort retained;Minimum SNR;Retained fraction", thresholds, "snr_min", "retained_fraction", None, objects)
    cuts.SaveAs(str(output / "snr_threshold_scan.png")); cuts.SaveAs(str(output / "snr_threshold_scan.pdf"))
    rf = ROOT.TFile(str(output / "snr_timing_plots.root"), "RECREATE")
    for canvas in (overview, mechanisms, cuts):
        canvas.Write()
    for obj in objects:
        if obj.InheritsFrom("TH1") or obj.InheritsFrom("TGraph"):
            obj.Write()
    rf.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, nargs="+")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output = (args.output_dir or PROJECT / "debug/reports" / ("bm2_snr_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))).resolve()
    output.mkdir(parents=True, exist_ok=False)
    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.TH1.AddDirectory(False)
    all_pads, joined, summaries, bins, thresholds, bands, event_rows = {}, {}, [], [], [], [], []
    for path in args.runs:
        run = path.resolve(); name = run.name
        pads, best = read_pad_rows(run)
        events, availability, total = join_events(run, best)
        all_pads[name] = pads; joined[name] = events
        full = timing_summary(events)
        stored = json.loads((run / "Combine/time_resolution/time_resolution.json").read_text())
        expected = stored["solutions"]["mean"]["weighted_core"]["resolutions_ns"][stored["active_detectors"].index("BM2")]
        if abs(full["bm2_sigma_ns"] - expected) > 1.e-8:
            raise RuntimeError(f"{run}: diagnostic cohort does not reproduce production resolution")
        summary = {"run": name, "step2_bm2_waveforms": len(pads), "track_events": total,
                   "availability": availability, "timing": full,
                   "common_four_detector_timing": timing_summary([r for r in events if r["all_four_valid"]]),
                   "step2_amplitude": describe([r["maxA"] for r in pads]),
                   "timing_amplitude": describe([r["maxA"] for r in events]),
                   "timing_noise": describe([r["peak_noise_rms"] for r in events]),
                   "timing_snr": describe([r["snr"] for r in events]),
                   "timing_jitter20_proxy": describe([r["jitter20_proxy_ns"] for r in events]),
                   "timing_rise10_90": describe([r["sigmoid_rise10_90_ns"] for r in events]),
                   "invalid_noise_events": sum(not math.isfinite(r["snr"]) for r in events),
                   "step2_high_amplitude_fraction_ge1400": sum(r["maxA"] >= 1400 for r in pads) / len(pads)}
        summaries.append(summary); bins.extend(bin_rows(name, events))
        for cut in (0., 15., 25., 40., 60., 90.):
            selected = [r for r in events if math.isfinite(r["snr"]) and r["snr"] >= cut]
            row = {"run": name, "snr_min": cut, "retained_fraction": len(selected) / len(events)}
            row.update(timing_summary(selected)); thresholds.append(row)
        for low, high in ((25., 90.), (40., 90.), (25., 60.)):
            selected = [r for r in events if low <= r["snr"] < high]
            row = {"run": name, "snr_min": low, "snr_max": high,
                   "retained_fraction": len(selected) / len(events)}
            row.update(timing_summary(selected)); bands.append(row)
        for row in events:
            record = {k: v for k, v in row.items() if k not in ("mean", "min")}
            for i, detector in enumerate(("BM1",) + ACTIVE):
                record[detector + "_time_ns"] = row["mean"][i]
            event_rows.append(record)
        print(name, json.dumps(clean(summary), ensure_ascii=False), flush=True)
    # Rows can have absent pair fields after an unsuccessful low-statistics fit.
    for rows in (bins, thresholds, bands):
        fields = list(dict.fromkeys(k for row in rows for k in row))
        for row in rows:
            for key in fields:
                row.setdefault(key, math.nan)
        if rows:
            rows[:] = [{key: row[key] for key in fields} for row in rows]
    write_tsv(output / "snr_bins.tsv", bins)
    write_tsv(output / "snr_thresholds.tsv", thresholds)
    write_tsv(output / "snr_bands.tsv", bands)
    write_tsv(output / "matched_bm2_events.tsv", event_rows)
    report = {"runs": summaries, "snr_definition": "maxA/peak_noise_rms",
              "cohort": "production TrackArrivalEvents, bmS_n[1]==1; exact 3x3 pad join; 25 ns/sample; pair-specific mean-t20 cohorts",
              "jitter_proxy": "25*k*noise/(0.16*H), local-crossing noise proxy, not fit covariance uncertainty",
              "uncertainty": "production weighted-core approximate pair-variance propagation; pair correlations and systematic errors not included",
              "bins": bins, "thresholds": thresholds, "bands": bands}
    (output / "snr_timing_report.json").write_text(json.dumps(clean(report), indent=2, ensure_ascii=False) + "\n")
    make_plots(output, all_pads, joined, bins, thresholds)
    print("Output:", output)


if __name__ == "__main__":
    main()

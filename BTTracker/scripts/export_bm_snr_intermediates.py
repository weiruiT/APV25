#!/usr/bin/env python3
"""Replay a saved BM2 SNR diagnostic and export its intermediate distributions.

Input is the frozen matched_bm2_events.tsv + snr_bins.tsv from
analyze_bm_snr_timing.py. No production input or result is changed.
Needs only PyROOT and the Python standard library.
"""
from __future__ import annotations

import argparse
from array import array
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import statistics

import ROOT
from analyze_bm_snr_timing import ACTIVE, clean, write_tsv
from analyze_time_resolution import analyze, clipped_stats, mean_std, solve_linear

REQUESTED_BINS = ((15., 25.), (25., 40.), (40., 60.), (60., 90.))
MINIMUM_BIN_EVENTS = 40
TIME_COLUMNS = ("BM1",) + ACTIVE


def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def load_records(path):
    groups = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for source in csv.DictReader(stream, delimiter="\t"):
            record = {"event": int(source["event"]), "snr": float(source["snr"]),
                      "maxA": float(source["maxA"]),
                      "peak_noise_rms": float(source["peak_noise_rms"]),
                      "position": int(source["position"]),
                      "global_channel": int(source["global_channel"]),
                      "step2_entry": int(source["step2_entry"]),
                      "mean": [float(source[d + "_time_ns"]) for d in TIME_COLUMNS]}
            groups.setdefault(source["run"], []).append(record)
    return groups


def clip_with_trace(values):
    """Replay production clipped_stats, preserving exact final membership."""
    center = statistics.median(values)
    mad = statistics.median(abs(value - center) for value in values)
    sigma = 1.4826 * mad
    if not math.isfinite(sigma) or sigma <= 0.:
        _, sigma = mean_std(values)
    mask = [True] * len(values)
    trace = []
    for iteration in range(12):
        if not math.isfinite(sigma) or sigma <= 0.:
            break
        low, high = center - 3. * sigma, center + 3. * sigma
        updated = [abs(value - center) <= 3. * sigma for value in values]
        selected = [value for value, keep in zip(values, updated) if keep]
        if len(selected) < 3:
            break
        new_center, new_sigma = mean_std(selected)
        stable = sum(updated) == sum(mask) and abs(new_center - center) < 1.e-10
        trace.append({"iteration": iteration + 1, "input_center_ns": center,
                      "input_sigma_ns": sigma, "selection_low_ns": low,
                      "selection_high_ns": high, "kept_events": len(selected),
                      "output_center_ns": new_center, "output_sigma_ns": new_sigma,
                      "stable": stable})
        mask = updated
        center, sigma = new_center, new_sigma
        if stable:
            break
    expected = clipped_stats(values)
    actual_center, actual_sigma = mean_std([v for v, keep in zip(values, mask) if keep])
    if (sum(mask) != expected["entries"] or
            abs(actual_center - expected["mean"]) > 1.e-10 or
            abs(actual_sigma - expected["sigma"]) > 1.e-10):
        raise RuntimeError("clipping replay does not reproduce production statistics")
    return mask, trace


def annotate(row, compact=False):
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextSize(.102 if compact else .080)
    label.DrawLatex(.06, .89, f"N={row['events']}, N_{{core}}={row['core_events']}")
    label.DrawLatex(.06, .745,
                    f"core #mu={row['core_mean_ns']:.2f}, #sigma={row['core_sigma_ns']:.2f} ns")
    label.DrawLatex(.06, .60,
                    f"Gaussian #mu={row['gaussian_mu_ns']:.2f} #pm {row['gaussian_mu_error_ns']:.2f} ns")
    label.DrawLatex(.06, .455,
                    f"Gaussian #sigma={row['gaussian_sigma_ns']:.2f} #pm {row['gaussian_sigma_error_ns']:.2f} ns")
    label.DrawLatex(.06, .31,
                    f"fit=[{row['fit_low_ns']:.1f}, {row['fit_high_ns']:.1f}] ns")
    label.DrawLatex(.06, .165,
                    f"status={row['gaussian_status']}, cov={row['covariance_status']}, #chi^{{2}}/ndf={row['gaussian_chi2']:.1f}/{row['gaussian_ndf']}")
    return label


def zoom_panel(container, suffix, row, graphic, compact=False):
    container.cd()
    plot = ROOT.TPad("plot_" + suffix, "Distribution", 0., .30, 1., 1.)
    info = ROOT.TPad("info_" + suffix, "Fit information", 0., 0., 1., .30)
    plot.SetTopMargin(.12); plot.SetBottomMargin(.13)
    plot.Draw(); info.Draw()
    plot.cd(); draw_zoom(*graphic)
    info.cd(); label = annotate(row, compact)
    return plot, info, label


def canvas_with_header(name, title, width, height, heading, columns, rows):
    canvas = ROOT.TCanvas(name, title, width, height)
    header = ROOT.TPaveText(0., .95, 1., 1., "NDC")
    header.SetBorderSize(0); header.SetFillColor(ROOT.kWhite)
    header.SetTextSize(.035); header.AddText(heading); header.Draw()
    body = ROOT.TPad(name + "_body", title, 0., 0., 1., .95)
    body.Draw(); body.cd(); body.Divide(columns, rows)
    return canvas, body, header


def draw_zoom(raw, core, gaussian, bounds):
    ROOT.gPad.SetGrid()
    raw.Draw("HIST")
    core.Draw("HIST SAME")
    gaussian.Draw("SAME")
    for boundary in bounds:
        boundary.Draw()


def export_pair(directory, root_file, records, source_row):
    row = dict(source_row)
    first = TIME_COLUMNS.index(row["first"])
    second = TIME_COLUMNS.index(row["second"])
    samples = [record for record in records
               if math.isfinite(record["mean"][first]) and math.isfinite(record["mean"][second])]
    values = [record["mean"][first] - record["mean"][second] for record in samples]
    mask, trace = clip_with_trace(values)
    key = row["first"] + "_" + row["second"]
    row["pair"] = key
    low = row["core_mean_ns"] - 6. * row["core_sigma_ns"]
    high = row["core_mean_ns"] + 6. * row["core_sigma_ns"]
    fit_low = row["core_mean_ns"] - 3. * row["core_sigma_ns"]
    fit_high = row["core_mean_ns"] + 3. * row["core_sigma_ns"]
    title = f"{row['first']} - {row['second']};#Delta t [ns];Events / bin"
    raw = ROOT.TH1D("h_" + key, title, 120, low, high)
    core = ROOT.TH1D("h_core_" + key, title, 120, low, high)
    full = ROOT.TH1D("h_full_" + key, title, 150,
                     min(low, min(values) - 1.), max(high, max(values) + 1.))
    raw.SetLineColor(ROOT.kBlue + 1); raw.SetLineWidth(2)
    core.SetLineColor(ROOT.kGreen + 2); core.SetLineStyle(2); core.SetLineWidth(2)
    full.SetLineColor(ROOT.kBlue + 1); full.SetLineWidth(2)
    data_rows = []
    for record, value, keep in zip(samples, values, mask):
        raw.Fill(value); full.Fill(value)
        if keep:
            core.Fill(value)
        data_rows.append({"event": record["event"], "snr": record["snr"],
                          "first_time_ns": record["mean"][first],
                          "second_time_ns": record["mean"][second],
                          "delta_t_ns": value, "in_core": int(keep),
                          "bm2_position": record["position"],
                          "bm2_global_channel": record["global_channel"],
                          "bm2_step2_entry": record["step2_entry"]})
    gaussian = ROOT.TF1("gaus_" + key, "gaus", fit_low, fit_high)
    gaussian.SetParameters(raw.GetMaximum(), row["core_mean_ns"], row["core_sigma_ns"])
    gaussian.SetLineColor(ROOT.kRed + 1); gaussian.SetLineWidth(3)
    fit = raw.Fit(gaussian, "QRS0")
    row.update(fit_low_ns=fit_low, fit_high_ns=fit_high,
               gaussian_status=int(fit), covariance_status=int(fit.CovMatrixStatus()),
               gaussian_mu_ns=float(gaussian.GetParameter(1)),
               gaussian_mu_error_ns=float(gaussian.GetParError(1)),
               gaussian_sigma_ns=abs(float(gaussian.GetParameter(2))),
               gaussian_sigma_error_ns=float(gaussian.GetParError(2)),
               gaussian_chi2=float(gaussian.GetChisquare()), gaussian_ndf=int(gaussian.GetNDF()),
               gaussian_probability=float(gaussian.GetProb()),
               gaussian_mu_sigma_correlation=float(fit.Correlation(1, 2)),
               zoom_low_ns=low, zoom_high_ns=high,
               zoom_underflow=int(raw.GetBinContent(0)),
               zoom_overflow=int(raw.GetBinContent(raw.GetNbinsX() + 1)))
    raw.SetMaximum(raw.GetMaximum() * 1.30)
    boundaries = []
    for boundary in (fit_low, fit_high):
        line = ROOT.TLine(boundary, 0., boundary, raw.GetMaximum())
        line.SetLineColor(ROOT.kGray + 2); line.SetLineStyle(2)
        boundaries.append(line)
    heading = f"{row['run']} | {row['snr_low']:g} <= SNR < {row['snr_high']:g} | {row['first']} - {row['second']}"
    canvas, body, header = canvas_with_header("c_" + key, key, 1600, 860, heading, 2, 1)
    body.cd(1); ROOT.gPad.SetGrid(); full.Draw("HIST")
    legend = ROOT.TLegend(.56, .74, .90, .88)
    legend.SetBorderSize(0); legend.AddEntry(full, "All pair events", "l")
    legend.Draw()
    full_label = ROOT.TLatex(); full_label.SetNDC(True); full_label.SetTextSize(.035)
    full_label.DrawLatex(.14, .91, "Full distribution (all events visible)")
    right = body.cd(2)
    panel_objects = zoom_panel(right, key, row, (raw, core, gaussian, boundaries))
    canvas.SaveAs(str(directory / (key + ".png")))
    canvas.SaveAs(str(directory / (key + ".pdf")))
    root_file.cd()
    for obj in (raw, core, full, gaussian, canvas):
        obj.Write()
    sample_tree = ROOT.TTree("samples_" + key, "Exact pair-specific events; in_core is production clipping membership")
    buffers = {name: array("i" if name in ("event", "in_core", "bm2_position", "bm2_global_channel", "bm2_step2_entry") else "d", [0])
               for name in data_rows[0]}
    for name, buffer in buffers.items():
        sample_tree.Branch(name, buffer, name + ("/I" if buffer.typecode == "i" else "/D"))
    for sample in data_rows:
        for name, buffer in buffers.items():
            buffer[0] = sample[name]
        sample_tree.Fill()
    sample_tree.Write()
    canvas.Close()
    write_tsv(directory / (key + "_events.tsv"), data_rows)
    write_tsv(directory / (key + "_clipping.tsv"), trace)
    return row, (raw, core, gaussian, boundaries), (legend, full_label, *panel_objects, sample_tree, buffers, body, header)


def matrix_details(pairs, weighted):
    design, weights, observations, equations = [], [], [], []
    for row in pairs:
        first, second = int(row["first_index"]), int(row["second_index"])
        coefficients = [int(index in (first, second)) for index in range(4)]
        observed = row["core_sigma_ns"] ** 2
        error = math.sqrt(2. * observed * observed / (row["core_events"] - 1))
        weight = 1. / (error * error)
        prediction = weighted["variances_ns2"][first] + weighted["variances_ns2"][second]
        equations.append({"first": row["first"], "second": row["second"],
                          "core_events": row["core_events"], "core_sigma_ns": row["core_sigma_ns"],
                          "observed_variance_ns2": observed, "variance_error_ns2": error,
                          "weight": weight, "predicted_variance_ns2": prediction,
                          "residual_ns2": observed - prediction,
                          "residual_over_error": (observed - prediction) / error,
                          **{detector + "_coefficient": coefficient
                             for detector, coefficient in zip(ACTIVE, coefficients)}})
        design.append(coefficients); weights.append(weight); observations.append(observed)
    normal = [[sum(weight * a[i] * a[j] for a, weight in zip(design, weights))
               for j in range(4)] for i in range(4)]
    rhs = [sum(weight * a[i] * d for a, weight, d in zip(design, weights, observations))
           for i in range(4)]
    solution = solve_linear(normal, rhs)
    if max(abs(a - b) for a, b in zip(solution, weighted["variances_ns2"])) > 1.e-9:
        raise RuntimeError("normal matrix does not reproduce production weighted solution")
    inverse_columns = [solve_linear(normal, [float(i == j) for i in range(4)]) for j in range(4)]
    covariance = [[inverse_columns[j][i] for j in range(4)] for i in range(4)]
    return equations, {"detector_order": ACTIVE, "A": design, "W_diagonal": weights,
                       "d_ns2": observations, "normal_matrix_A_transpose_W_A": normal,
                       "rhs_A_transpose_W_d": rhs, "solution_variances_ns2": solution,
                       "variance_covariance_ns4": covariance}


def export_bin(directory, run, low, high, records, reference):
    directory.mkdir(parents=True, exist_ok=False)
    pairs, raw_result, core_result, weighted = analyze(records, "mean", ACTIVE)
    expected_sigma = float(reference["bm2_sigma_ns"])
    if (len(records) != int(reference["selected_events"]) or
            abs(weighted["resolutions_ns"][0] - expected_sigma) > 1.e-9 or
            abs(weighted["resolution_uncertainties_ns"][0] - float(reference["bm2_sigma_error_ns"])) > 1.e-9):
        raise RuntimeError(f"{run} [{low},{high}): cannot reproduce saved table")
    for row in pairs:
        key = row["first"] + "_" + row["second"]
        if (int(reference[key + "_n"]) != row["events"] or
                abs(float(reference[key + "_core_sigma_ns"]) - row["core_sigma_ns"]) > 1.e-9):
            raise RuntimeError(f"{run} [{low},{high}): pair {key} differs from saved table")
    event_rows = []
    for record in records:
        event_rows.append({k: v for k, v in record.items() if k != "mean"} |
                          {detector + "_time_ns": value for detector, value in zip(TIME_COLUMNS, record["mean"])})
    write_tsv(directory / "selected_events.tsv", event_rows)
    root_file = ROOT.TFile(str(directory / "time_differences.root"), "RECREATE")
    if root_file.IsZombie():
        raise RuntimeError(f"cannot create ROOT output in {directory}")
    definition = ROOT.TNamed("AnalysisDefinition",
                            f"{run}; {low}<=SNR<{high}; BM2 single matched pad; mean sigmoid t20; 25 ns/sample; "
                            "core=production iterative 3sigma; Gaussian fit is diagnostic only, NOT used in variance solution")
    definition.Write()
    fit_rows, graphics, keep_alive = [], [], [definition]
    for pair in pairs:
        pair_metadata = dict(pair, run=run, snr_low=low, snr_high=high)
        row, graphic, keep = export_pair(directory, root_file, records, pair_metadata)
        fit_rows.append(row); graphics.append(graphic); keep_alive.extend(keep)
    heading = (f"{run} | {low:g} <= SNR < {high:g} | Selected events: {len(records)} | "
               f"BM2 weighted-core resolution: {weighted['resolutions_ns'][0]:.3f} +/- {weighted['resolution_uncertainties_ns'][0]:.3f} ns")
    canvas, body, header = canvas_with_header("six_pairs", "All six time differences", 1800, 1320, heading, 3, 2)
    keep_alive.extend((body, header))
    for index, (row, graphic) in enumerate(zip(fit_rows, graphics), 1):
        panel = body.cd(index)
        keep_alive.extend(zoom_panel(panel, "combined_" + row["pair"], row, graphic, compact=True))
    canvas.SaveAs(str(directory / "time_differences_6pairs.png"))
    canvas.SaveAs(str(directory / "time_differences_6pairs.pdf"))
    root_file.cd(); canvas.Write(); root_file.Close(); canvas.Close()
    write_tsv(directory / "pairwise_statistics.tsv", fit_rows)
    equations, matrices = matrix_details(pairs, weighted)
    write_tsv(directory / "variance_equations.tsv", equations)
    detector_rows = [{"detector": detector, "variance_ns2": weighted["variances_ns2"][i],
                      "variance_error_ns2": weighted["variance_uncertainties_ns2"][i],
                      "resolution_ns": weighted["resolutions_ns"][i],
                      "resolution_error_ns": weighted["resolution_uncertainties_ns"][i]}
                     for i, detector in enumerate(ACTIVE)]
    write_tsv(directory / "detector_solution.tsv", detector_rows)
    report = {"run": run, "snr_low": low, "snr_high": high, "selected_events": len(records),
              "reproduced_saved_table": True, "source_table_bm2_sigma_ns": expected_sigma,
              "pairs": fit_rows, "solutions": {"raw": raw_result, "core": core_result,
                                                "weighted_core": weighted}, "matrices": matrices,
              "gaussian_role": "diagnostic: fitted to raw zoom histogram over core mean +/-3 core sigma; not used in original table",
              "error_definition": "same production approximate variance-estimator propagation; pair correlations/systematics not included"}
    (directory / "solution.json").write_text(json.dumps(clean(report), indent=2) + "\n")
    print(f"{run} SNR[{low:g},{high:g}): N={len(records)}, BM2={weighted['resolutions_ns'][0]:.6f} +/- {weighted['resolution_uncertainties_ns'][0]:.6f} ns; verified", flush=True)
    return {"run": run, "snr_low": low, "snr_high": high, "selected_events": len(records),
            "bm2_sigma_ns": weighted["resolutions_ns"][0],
            "bm2_sigma_error_ns": weighted["resolution_uncertainties_ns"][0],
            "model_chi2": weighted["chi2"], "model_ndf": weighted["ndf"],
            "directory": str(directory), "verified": True,
            "gaussian_failed_pairs": sum(row["gaussian_status"] != 0 for row in fit_rows)}


def write_index(output, summary, provenance):
    completed = sum(row["status"] == "completed" for row in summary)
    lines = ["# SNR 分组时间分析中间结果", "", f"SNR 区间为左闭右开；{completed}/{len(summary)} 组达到统计门槛并完成方差解核对。", "",
             "六个探测器对：BM2−Tracker1、BM2−Tracker2、BM2−Tracker3、Tracker1−Tracker2、Tracker1−Tracker3、Tracker2−Tracker3。", "",
             "每个图的蓝线为全部有效时间差，绿虚线为生产算法保留的核心事件，红线为另行添加的高斯拟合。",
             "原表的方差求解使用核心事件的样本标准差，不使用图中的高斯拟合 sigma。", "",
             "| Run | SNR | 选中事件 | BM2 sigma (ns) | 六组时间差 | ROOT | 方差求解 |",
             "|---|---|---:|---:|---|---|---|"]
    for row in summary:
        relative = Path(row["directory"]).relative_to(output).as_posix()
        if row["status"] == "completed":
            resolution = f"{row['bm2_sigma_ns']:.3f} ± {row['bm2_sigma_error_ns']:.3f}"
            plots = f"[PNG]({relative}/time_differences_6pairs.png) / [PDF]({relative}/time_differences_6pairs.pdf)"
            root = f"[ROOT]({relative}/time_differences.root)"
            solution = f"[JSON]({relative}/solution.json)"
        else:
            resolution = row["status"]
            plots = root = solution = f"[状态]({relative}/status.json)"
        lines.append(f"| {row['run']} | [{row['snr_low']:g}, {row['snr_high']:g}) | {row['selected_events']} | "
                     f"{resolution} | {plots} | {root} | {solution} |")
    lines += ["", "## 每个子目录", "",
              "- `selected_events.tsv`：SNR cut 后的事件编号、BM2 pad/Step2 entry、SNR 和四个探测器的时间。",
              "- `<探测器对>_events.tsv`：双方时间都有效的事件、时间差、是否保留在核心中。",
              "- `<探测器对>_clipping.tsv`：每轮 3sigma 筛选范围、事件数、更新后的均值与标准差。",
              "- `<探测器对>.png/pdf`：左图展示包含尾部的完整分布，右图展示核心附近、高斯拟合及参数。",
              "- `pairwise_statistics.tsv`：各对原始数、RMS、核心数、核心均值/标准差、高斯参数/误差、范围和拟合状态。",
              "- `variance_equations.tsv`：六个方差方程的观测值、误差、权重、模型预测和残差。",
              "- `detector_solution.tsv`：四个探测器的方差、时间分辨和近似统计误差。",
              "- `solution.json`：全部数值，包括 A、W、正规矩阵、右端向量、方差协方差和原表核对状态。",
              "- `time_differences.root`：全部/核心/完整范围直方图、高斯函数、画布和逐事件 `samples_*` tree。",
              "", "## 方差求解", "",
              "令 x 为四个探测器方差，d 为六个时间差核心 sigma 的平方。每行 A 在对应的两个探测器列为 1。",
              "每对使用自己的有效事件集合。近似观测误差为 sqrt(2*d^2/(N_core-1))，权重为误差平方的倒数。",
              "解 (A^T W A)x = A^T W d；单探测器分辨为 sqrt(x)，误差由正规矩阵逆传播。",
              "这套近似不包含共享事件造成的方程间相关性和其他系统误差。", "",
              "## 来源", "", f"源目录：`{provenance['diagnostic_directory']}`", "",
              "完整输入文件 SHA256 与生成信息保存在 `provenance.json`。"]
    (output / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("diagnostic_directory", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--snr-bin-edges", help="comma-separated finite edges")
    parser.add_argument("--single-run", help="run name for a run-local output directory")
    args = parser.parse_args()
    if args.snr_bin_edges:
        try:
            edges = tuple(float(value) for value in args.snr_bin_edges.split(","))
        except ValueError:
            parser.error("--snr-bin-edges must be comma-separated numbers")
        if (len(edges) < 2 or any(not math.isfinite(value) or value < 0 for value in edges)
                or any(second <= first for first, second in zip(edges, edges[1:]))):
            parser.error("--snr-bin-edges must contain at least two increasing non-negative values")
        requested_bins = tuple(zip(edges, edges[1:]))
    else:
        requested_bins = REQUESTED_BINS
    source = args.diagnostic_directory.resolve()
    events_path, bins_path = source / "matched_bm2_events.tsv", source / "snr_bins.tsv"
    records = load_records(events_path)
    if args.single_run and records and set(records) != {args.single_run}:
        parser.error("--single-run must match the only run in matched_bm2_events.tsv")
    with bins_path.open(newline="", encoding="utf-8") as stream:
        references = {(row["run"], float(row["snr_low"]), float(row["snr_high"])): row
                      for row in csv.DictReader(stream, delimiter="\t")}
    output = (args.output_dir or source / ("pairwise_details_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))).resolve()
    output.mkdir(parents=True, exist_ok=False)
    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptFit(0)
    ROOT.TH1.AddDirectory(False)
    summary = []
    for run, cohort in records.items():
        for low, high in requested_bins:
            reference = references.get((run, low, high))
            selected = [record for record in cohort if low <= record["snr"] < high]
            parent = output if args.single_run else output / run
            directory = parent / f"SNR_{low:g}_{high:g}"
            if reference is None:
                if len(selected) >= MINIMUM_BIN_EVENTS:
                    raise RuntimeError(
                        f"reference table missing {run} [{low},{high}) despite "
                        f"{len(selected)} events"
                    )
                directory.mkdir(parents=True, exist_ok=False)
                event_rows = [{k: v for k, v in record.items() if k != "mean"} |
                              {detector + "_time_ns": value
                               for detector, value in zip(TIME_COLUMNS, record["mean"])}
                              for record in selected]
                if event_rows:
                    write_tsv(directory / "selected_events.tsv", event_rows)
                else:
                    (directory / "selected_events.tsv").write_text(
                        "event\tsnr\tstep2_entry\tBM2_time_ns\tTracker1_time_ns\t"
                        "Tracker2_time_ns\tTracker3_time_ns\n", encoding="utf-8"
                    )
                status = {"status": "insufficient_statistics",
                          "selected_events": len(selected),
                          "minimum_events_for_timing_bin": MINIMUM_BIN_EVENTS,
                          "message": "The timing diagnostic did not calculate this bin because it has fewer than the configured minimum events."}
                (directory / "status.json").write_text(json.dumps(status, indent=2) + "\n")
                summary.append({"run": run, "snr_low": low, "snr_high": high,
                                "selected_events": len(selected), "status": "insufficient_statistics",
                                "directory": str(directory), "verified": False,
                                "gaussian_failed_pairs": 0, "pair_distributions": 0})
                continue
            row = export_bin(directory, run, low, high, selected, reference)
            row["status"] = "completed"
            row["pair_distributions"] = 6
            summary.append(row)
    summary_fields = list(dict.fromkeys(key for row in summary for key in row))
    for row in summary:
        for field in summary_fields:
            row.setdefault(field, "")
    summary[:] = [{field: row[field] for field in summary_fields} for row in summary]
    write_tsv(output / "summary.tsv", summary)
    provenance = {"diagnostic_directory": str(source),
                  "created_utc": datetime.now(timezone.utc).isoformat(),
                  "source_sha256": {path.name: digest(path) for path in (events_path, bins_path)},
                  "script_sha256": digest(Path(__file__)), "ROOT_version": ROOT.gROOT.GetVersion(),
                  "groups": len(summary), "pair_distributions": sum(row["pair_distributions"] for row in summary),
                  "all_computed_results_verified": all(row["verified"] for row in summary
                                                        if row["status"] == "completed"),
                  "insufficient_statistics_groups": sum(row["status"] != "completed" for row in summary),
                  "gaussian_failed_pairs": sum(row["gaussian_failed_pairs"] for row in summary)}
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    write_index(output, summary, provenance)
    print("Output:", output, flush=True)


if __name__ == "__main__":
    main()

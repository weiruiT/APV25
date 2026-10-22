#!/usr/bin/env python3
"""Estimate per-detector timing resolution from TrackArrivalEvents.

The stored sigmoid t20/t30/t40/t50 values are in APV samples (25 ns/sample).
A Tracker
event time is the strip-count-weighted mean of its matched X/Y clusters.  The
Each detector pair uses all events available to that pair; a Tracker pools
whichever matched X/Y axes are valid.  Common trigger jitter is removed by
using pairwise time differences.  Individual variances follow

    Var(t_i - t_j) = sigma_i^2 + sigma_j^2

with an over-constrained least-squares solution over all configured detector pairs.
Both full-distribution sample RMS and an iterative 3-sigma core result are
reported; the core result is the recommended Gaussian-core estimate.
"""

from __future__ import annotations

import argparse
from array import array
import csv
import json
import math
from pathlib import Path
import statistics
import sys

try:
    import ROOT
except ImportError as exc:
    raise SystemExit(
        "PyROOT is unavailable; source ROOT's thisroot.sh before running"
    ) from exc


SAMPLE_NS = 25.0
DETECTORS = ("BM1", "BM2", "Tracker1", "Tracker2", "Tracker3")
TIMING_METHODS = ("t20", "t30", "t40", "t50")


def valid_time(value: float) -> bool:
    return math.isfinite(value) and value > -900.0


def mean_std(values: list[float]) -> tuple[float, float]:
    if len(values) < 2:
        return math.nan, math.nan
    mean = math.fsum(values) / len(values)
    variance = math.fsum((value - mean) ** 2 for value in values) / (len(values) - 1)
    return mean, math.sqrt(max(0.0, variance))


def clipped_stats(values: list[float], nsigma: float = 3.0) -> dict[str, float | int]:
    """Iterative Gaussian-core estimate, initialized with median/MAD."""
    if len(values) < 3:
        return {"entries": len(values), "mean": math.nan, "sigma": math.nan}
    center = statistics.median(values)
    mad = statistics.median(abs(value - center) for value in values)
    sigma = 1.4826 * mad
    if not math.isfinite(sigma) or sigma <= 0.0:
        _, sigma = mean_std(values)
    selected = list(values)
    for _ in range(12):
        if not math.isfinite(sigma) or sigma <= 0.0:
            break
        updated = [value for value in values if abs(value - center) <= nsigma * sigma]
        if len(updated) < 3:
            break
        new_center, new_sigma = mean_std(updated)
        stable = len(updated) == len(selected) and abs(new_center - center) < 1.0e-10
        selected = updated
        center, sigma = new_center, new_sigma
        if stable:
            break
    return {"entries": len(selected), "mean": center, "sigma": sigma}


def solve_variances(
    pair_rows: list[dict[str, object]], sigma_key: str, detector_count: int
) -> dict[str, object]:
    """Unweighted LS solution for all x_i in d_ij = x_i + x_j."""
    count = detector_count
    sums = [0.0] * count
    observations: list[tuple[int, int, float]] = []
    for row in pair_rows:
        sigma = float(row[sigma_key])
        if not math.isfinite(sigma):
            continue
        variance = sigma * sigma
        first = int(row["first_index"])
        second = int(row["second_index"])
        sums[first] += variance
        sums[second] += variance
        observations.append((first, second, variance))
    if len(observations) != count * (count - 1) // 2:
        raise RuntimeError("incomplete pairwise timing matrix")
    total_variance = math.fsum(sums) / (2.0 * (count - 1))
    variances = [(value - total_variance) / (count - 2) for value in sums]
    residuals = [
        observed - variances[first] - variances[second]
        for first, second, observed in observations
    ]
    residual_rms = math.sqrt(math.fsum(value * value for value in residuals) / len(residuals))
    return {
        "variances_ns2": variances,
        "resolutions_ns": [math.sqrt(value) if value >= 0.0 else math.nan for value in variances],
        "pair_variance_residual_rms_ns2": residual_rms,
    }


def solve_linear(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    size = len(rhs)
    work = [list(row) + [rhs[index]] for index, row in enumerate(matrix)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(work[row][column]))
        if abs(work[pivot][column]) < 1.0e-20:
            raise RuntimeError("singular timing-resolution normal matrix")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            work[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(work[row], work[column])
            ]
    return [work[index][-1] for index in range(size)]


def solve_weighted_variances(
    pair_rows: list[dict[str, object]], sigma_key: str, entries_key: str,
    detector_count: int,
) -> dict[str, object]:
    """Variance fit weighted by the approximate variance-estimator error."""
    count = detector_count
    normal = [[0.0] * count for _ in range(count)]
    rhs = [0.0] * count
    observations: list[tuple[int, int, float, float]] = []
    for row in pair_rows:
        sigma = float(row[sigma_key])
        entries = int(row[entries_key])
        if not math.isfinite(sigma) or sigma <= 0.0 or entries < 3:
            continue
        variance = sigma * sigma
        variance_error2 = 2.0 * variance * variance / (entries - 1)
        weight = 1.0 / variance_error2
        first = int(row["first_index"])
        second = int(row["second_index"])
        normal[first][first] += weight
        normal[first][second] += weight
        normal[second][first] += weight
        normal[second][second] += weight
        rhs[first] += weight * variance
        rhs[second] += weight * variance
        observations.append((first, second, variance, weight))
    if len(observations) != count * (count - 1) // 2:
        raise RuntimeError("incomplete weighted pairwise timing matrix")
    variances = solve_linear(normal, rhs)
    inverse_columns = []
    for column in range(count):
        unit = [0.0] * count
        unit[column] = 1.0
        inverse_columns.append(solve_linear(normal, unit))
    variance_errors = [
        math.sqrt(max(0.0, inverse_columns[index][index]))
        for index in range(count)
    ]
    resolutions = [math.sqrt(value) if value >= 0.0 else math.nan for value in variances]
    resolution_errors = [
        variance_errors[index] / (2.0 * resolutions[index])
        if math.isfinite(resolutions[index]) and resolutions[index] > 0.0 else math.nan
        for index in range(count)
    ]
    residuals = [
        observed - variances[first] - variances[second]
        for first, second, observed, _ in observations
    ]
    chi2 = math.fsum(
        weight * (observed - variances[first] - variances[second]) ** 2
        for first, second, observed, weight in observations
    )
    return {
        "variances_ns2": variances,
        "variance_uncertainties_ns2": variance_errors,
        "resolutions_ns": resolutions,
        "resolution_uncertainties_ns": resolution_errors,
        "pair_variance_residual_rms_ns2": math.sqrt(
            math.fsum(value * value for value in residuals) / len(residuals)),
        "chi2": chi2,
        "ndf": len(observations) - count,
    }

def read_events(
    path: Path, bm2_matched_pad_count: int = 0, bm_time_source: str = "pad"
) -> tuple[list[dict[str, object]], dict[str, int], int]:
    if bm_time_source not in ("pad", "cluster"):
        raise ValueError("bm_time_source must be pad or cluster")
    if bm_time_source == "cluster" and bm2_matched_pad_count != 0:
        raise ValueError("cluster timing must not apply a single-pad multiplicity cut")
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"cannot open {path}")
    tree = root_file.Get("TrackArrivalEvents")
    if not tree:
        raise RuntimeError(f"TrackArrivalEvents is missing in {path}")

    event = array("i", [0])
    bm_mean = {method: array("d", [0.0] * 2) for method in TIMING_METHODS}
    bm_min = {method: array("d", [0.0] * 2) for method in TIMING_METHODS}
    bm_count = {method: array("i", [0] * 2) for method in TIMING_METHODS}
    cluster_count = array("i", [0] * 2)
    tracker_mean = {method: array("d", [0.0] * 6) for method in TIMING_METHODS}
    tracker_min = {method: array("d", [0.0] * 6) for method in TIMING_METHODS}
    tracker_count = {method: array("i", [0] * 6) for method in TIMING_METHODS}
    for branch, target in (("event", event),):
        if not tree.GetBranch(branch):
            raise RuntimeError(f"required branch {branch} is missing")
        tree.SetBranchAddress(branch, target)
    bm_stem = "bmClusterS" if bm_time_source == "cluster" else "bmS"
    for method in TIMING_METHODS:
        suffix = "" if method == "t20" else method[1:]
        branches = (
            (f"{bm_stem}{suffix}_mean", bm_mean[method]),
            (f"{bm_stem}{suffix}_min", bm_min[method]),
            (f"{bm_stem}{suffix}_n", bm_count[method]),
            (f"trS{suffix}_mean", tracker_mean[method]),
            (f"trS{suffix}_min", tracker_min[method]),
            (f"trS{suffix}_n", tracker_count[method]),
        )
        for branch, target in branches:
            if not tree.GetBranch(branch):
                raise RuntimeError(
                    f"required branch {branch} is missing; rerun the signal workflow to rebuild Step2/Step3"
                )
            tree.SetBranchAddress(branch, target)
    if bm_time_source == "cluster":
        if not tree.GetBranch("bm_matched_cluster_count"):
            raise RuntimeError("cluster timing branches are missing; rerun pad reconstruction")
        tree.SetBranchAddress("bm_matched_cluster_count", cluster_count)

    availability = {name: 0 for name in DETECTORS}
    records: list[dict[str, object]] = []
    total_entries = int(tree.GetEntries())
    for entry in range(total_entries):
        tree.GetEntry(entry)
        if bm_time_source == "cluster" and (
            cluster_count[1] != 1 or bm_count["t20"][1] <= 0
        ):
            continue
        if (
            bm2_matched_pad_count > 0
            and int(bm_count["t20"][1]) != bm2_matched_pad_count
        ):
            continue
        record: dict[str, object] = {"event": int(event[0])}
        for method in TIMING_METHODS:
            mean_times = [float(bm_mean[method][0]), float(bm_mean[method][1])]
            min_times = [float(bm_min[method][0]), float(bm_min[method][1])]
            detector_valid = [
                bm_count[method][index] > 0
                and valid_time(mean_times[index])
                and valid_time(min_times[index])
                for index in range(2)
            ]
            if bm_time_source == "cluster":
                detector_valid = [
                    valid and cluster_count[index] == 1
                    for index, valid in enumerate(detector_valid)
                ]
            for tracker in range(3):
                x_index = 2 * tracker
                y_index = x_index + 1
                nx = int(tracker_count[method][x_index])
                ny = int(tracker_count[method][y_index])
                x_valid = (
                    nx > 0
                    and valid_time(float(tracker_mean[method][x_index]))
                    and valid_time(float(tracker_min[method][x_index]))
                )
                y_valid = (
                    ny > 0
                    and valid_time(float(tracker_mean[method][y_index]))
                    and valid_time(float(tracker_min[method][y_index]))
                )
                tracker_valid = x_valid or y_valid
                detector_valid.append(tracker_valid)
                if tracker_valid:
                    valid_nx = nx if x_valid else 0
                    valid_ny = ny if y_valid else 0
                    mean_times.append(
                        (
                            valid_nx * float(tracker_mean[method][x_index])
                            + valid_ny * float(tracker_mean[method][y_index])
                        )
                        / (valid_nx + valid_ny)
                    )
                    minima = []
                    if x_valid:
                        minima.append(float(tracker_min[method][x_index]))
                    if y_valid:
                        minima.append(float(tracker_min[method][y_index]))
                    min_times.append(min(minima))
                else:
                    mean_times.append(math.nan)
                    min_times.append(math.nan)
            if method == "t20":
                for name, available in zip(DETECTORS, detector_valid):
                    availability[name] += int(available)
            for index, available in enumerate(detector_valid):
                if not available:
                    mean_times[index] = math.nan
                    min_times[index] = math.nan
            record["mean" if method == "t20" else f"mean_{method}"] = [
                value * SAMPLE_NS for value in mean_times
            ]
            record["min" if method == "t20" else f"min_{method}"] = [
                value * SAMPLE_NS for value in min_times
            ]
        records.append(record)
    tree.ResetBranchAddresses()
    root_file.Close()
    return records, availability, total_entries


def analyze(
    records: list[dict[str, object]], estimator: str,
    active_detectors: tuple[str, ...], timing_method: str = "t20",
) -> tuple[
    list[dict[str, object]], dict[str, object], dict[str, object], dict[str, object]
]:
    if timing_method not in TIMING_METHODS:
        raise ValueError(f"unsupported sigmoid CFD method: {timing_method}")
    record_key = estimator if timing_method == "t20" else f"{estimator}_{timing_method}"
    rows: list[dict[str, object]] = []
    detector_indices = [DETECTORS.index(name) for name in active_detectors]
    for first_local, first in enumerate(detector_indices):
        for second_local in range(first_local + 1, len(detector_indices)):
            second = detector_indices[second_local]
            differences = [
                float(record[record_key][first]) - float(record[record_key][second])
                for record in records
                if math.isfinite(float(record[record_key][first]))
                and math.isfinite(float(record[record_key][second]))
            ]
            if len(differences) < 3:
                raise RuntimeError(f"too few events for {DETECTORS[first]}-{DETECTORS[second]}")
            raw_mean, raw_sigma = mean_std(differences)
            core = clipped_stats(differences)
            rows.append(
                {
                    "estimator": estimator,
                    "timing_method": timing_method,
                    "first_index": first_local,
                    "second_index": second_local,
                    "first": DETECTORS[first],
                    "second": DETECTORS[second],
                    "events": len(differences),
                    "mean_ns": raw_mean,
                    "rms_ns": raw_sigma,
                    "core_events": int(core["entries"]),
                    "core_mean_ns": float(core["mean"]),
                    "core_sigma_ns": float(core["sigma"]),
                }
            )
    return (
        rows,
        solve_variances(rows, "rms_ns", len(active_detectors)),
        solve_variances(rows, "core_sigma_ns", len(active_detectors)),
        solve_weighted_variances(
            rows, "core_sigma_ns", "core_events", len(active_detectors)
        ),
    )


def plot_bm2_tracker_gaussians(
    records: list[dict[str, object]], output: Path, bm2_matched_pad_count: int,
    selection_label: str | None = None, timing_method: str = "t20",
) -> list[dict[str, object]]:
    """Fit and plot BM2-versus-Tracker timing differences for one CFD level."""
    if timing_method not in TIMING_METHODS:
        raise ValueError(f"unsupported sigmoid CFD method: {timing_method}")
    output.mkdir(parents=True, exist_ok=True)
    record_key = "mean" if timing_method == "t20" else f"mean_{timing_method}"
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.TH1.AddDirectory(False)
    pairs = ((1, 2), (1, 3), (1, 4))
    rows: list[dict[str, object]] = []
    objects: list[tuple[object, object, object, object]] = []
    selection_label = selection_label or (
        f"BM2 matched pads = {bm2_matched_pad_count}"
        if bm2_matched_pad_count > 0
        else "BM2 matched pad multiplicity uncut"
    )

    for first, second in pairs:
        first_name = DETECTORS[first]
        second_name = DETECTORS[second]
        differences = [
            float(record[record_key][first]) - float(record[record_key][second])
            for record in records
            if math.isfinite(float(record[record_key][first]))
            and math.isfinite(float(record[record_key][second]))
        ]
        if len(differences) < 3:
            raise RuntimeError(
                f"too few events for Gaussian fit {first_name}-{second_name}"
            )
        raw_mean, raw_sigma = mean_std(differences)
        core = clipped_stats(differences)
        core_mean = float(core["mean"])
        core_sigma = float(core["sigma"])
        if not math.isfinite(core_sigma) or core_sigma <= 0.0:
            raise RuntimeError(f"invalid core sigma for {first_name}-{second_name}")

        key = f"{first_name}_{second_name}"
        object_suffix = "" if timing_method == "t20" else f"_{timing_method}"
        display_low = core_mean - 6.0 * core_sigma
        display_high = core_mean + 6.0 * core_sigma
        fit_low = core_mean - 3.0 * core_sigma
        fit_high = core_mean + 3.0 * core_sigma
        histogram = ROOT.TH1D(
            f"h_{key}{object_suffix}",
            (
                f"{first_name} - {second_name};"
                f"#Delta t = t_{{{first_name}}} - t_{{{second_name}}} [ns];"
                "Events / bin"
            ),
            120,
            display_low,
            display_high,
        )
        for value in differences:
            histogram.Fill(value)
        histogram.SetLineColor(ROOT.kBlue + 1)
        histogram.SetLineWidth(2)

        gaussian = ROOT.TF1(
            f"gaus_{key}{object_suffix}", "gaus", fit_low, fit_high
        )
        gaussian.SetParameters(histogram.GetMaximum(), core_mean, core_sigma)
        gaussian.SetLineColor(ROOT.kRed + 1)
        gaussian.SetLineWidth(3)
        fit_result = histogram.Fit(gaussian, "QRS0")
        fit_status = int(fit_result)
        covariance_status = int(fit_result.CovMatrixStatus())
        gaussian_mu = float(gaussian.GetParameter(1))
        gaussian_mu_error = float(gaussian.GetParError(1))
        gaussian_sigma = abs(float(gaussian.GetParameter(2)))
        gaussian_sigma_error = float(gaussian.GetParError(2))
        mu_sigma_correlation = float(fit_result.Correlation(1, 2))
        chi2 = float(gaussian.GetChisquare())
        ndf = int(gaussian.GetNDF())
        probability = float(gaussian.GetProb())

        line_low = ROOT.TLine(
            fit_low, 0.0, fit_low, histogram.GetMaximum() * 1.03
        )
        line_high = ROOT.TLine(
            fit_high, 0.0, fit_high, histogram.GetMaximum() * 1.03
        )
        for line in (line_low, line_high):
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kGray + 2)

        canvas = ROOT.TCanvas(
            f"c_{key}{object_suffix}", "time difference", 1000, 720
        )
        canvas.SetGrid()
        histogram.Draw("HIST")
        gaussian.Draw("SAME")
        line_low.Draw()
        line_high.Draw()
        label = ROOT.TLatex()
        label.SetNDC(True)
        label.SetTextSize(0.033)
        label.DrawLatex(0.14, 0.91, selection_label)
        label.DrawLatex(
            0.14,
            0.86,
            f"raw N={len(differences)}, core N={int(core['entries'])}",
        )
        label.DrawLatex(
            0.14,
            0.81,
            f"core mean={core_mean:.3f} ns, core #sigma={core_sigma:.3f} ns",
        )
        label.DrawLatex(
            0.14,
            0.76,
            (
                f"Gaussian #mu={gaussian_mu:.3f} #pm "
                f"{gaussian_mu_error:.3f} ns"
            ),
        )
        label.DrawLatex(
            0.14,
            0.71,
            (
                f"Gaussian #sigma={gaussian_sigma:.3f} #pm "
                f"{gaussian_sigma_error:.3f} ns"
            ),
        )
        label.DrawLatex(
            0.14,
            0.66,
            (
                f"fit=[{fit_low:.2f}, {fit_high:.2f}] ns, "
                f"status={fit_status}, cov={covariance_status}, "
                f"#chi^{{2}}/ndf={chi2:.1f}/{ndf}"
            ),
        )
        canvas.SaveAs(str(output / f"{key}_gaussian_fit.png"))
        canvas.Close()

        row = {
            "first": first_name,
            "second": second_name,
            "timing_method": timing_method,
            "raw_events": len(differences),
            "raw_mean_ns": raw_mean,
            "raw_rms_ns": raw_sigma,
            "core_events": int(core["entries"]),
            "core_mean_ns": core_mean,
            "core_sigma_ns": core_sigma,
            "fit_low_ns": fit_low,
            "fit_high_ns": fit_high,
            "gaussian_status": fit_status,
            "covariance_status": covariance_status,
            "gaussian_mu_ns": gaussian_mu,
            "gaussian_mu_error_ns": gaussian_mu_error,
            "gaussian_sigma_ns": gaussian_sigma,
            "gaussian_sigma_error_ns": gaussian_sigma_error,
            "mu_sigma_correlation": mu_sigma_correlation,
            "chi2": chi2,
            "ndf": ndf,
            "probability": probability,
        }
        rows.append(row)
        objects.append((histogram, gaussian, line_low, line_high))

    combined = ROOT.TCanvas(
        "bm2_tracker_combined", "BM2 minus Tracker time differences", 1800, 620
    )
    combined.Divide(3, 1)
    labels = []
    for index, (histogram, gaussian, line_low, line_high) in enumerate(objects, 1):
        combined.cd(index)
        ROOT.gPad.SetGrid()
        histogram.Draw("HIST")
        gaussian.Draw("SAME")
        line_low.Draw()
        line_high.Draw()
        row = rows[index - 1]
        label = ROOT.TLatex()
        label.SetNDC(True)
        label.SetTextSize(0.038)
        label.DrawLatex(0.14, 0.91, selection_label)
        label.DrawLatex(
            0.14,
            0.84,
            (
                f"core N={row['core_events']}, "
                f"#sigma_{{core}}={row['core_sigma_ns']:.3f} ns"
            ),
        )
        label.DrawLatex(
            0.14, 0.77, f"Gaussian #mu={row['gaussian_mu_ns']:.3f} ns"
        )
        label.DrawLatex(
            0.14,
            0.70,
            (
                f"Gaussian #sigma={row['gaussian_sigma_ns']:.3f} #pm "
                f"{row['gaussian_sigma_error_ns']:.3f} ns"
            ),
        )
        labels.append(label)
    combined.SaveAs(str(output / "bm2_tracker_time_differences.png"))
    combined.SaveAs(str(output / "bm2_tracker_time_differences.pdf"))
    combined.Close()

    root_output = ROOT.TFile(
        str(output / "bm2_tracker_time_differences.root"), "RECREATE"
    )
    if root_output.IsZombie():
        raise RuntimeError("cannot create BM2-Tracker Gaussian ROOT output")
    definition = ROOT.TNamed(
        "TimingDifferenceDefinition",
        (
            f"mean sigmoid {timing_method}; 25 ns/sample; pair-specific valid cohorts; "
            f"{selection_label}; "
            "core=iterative 3-sigma initialized by median/MAD; histogram "
            "range=core mean +/-6 sigma; Gaussian fit range=core mean +/-3 sigma"
        ),
    )
    definition.Write()
    for histogram, gaussian, _line_low, _line_high in objects:
        histogram.Write()
        gaussian.Write()
    root_output.Close()

    with (output / "bm2_tracker_gaussian_fits.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), delimiter="	"
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def serializable(value: object) -> object:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, list):
        return [serializable(item) for item in value]
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    return value


def plot_cfd_fraction_resolution(
    solutions_by_method: dict[str, dict[str, object]],
    active_detectors: tuple[str, ...],
    output: Path,
) -> None:
    """Compare weighted-core detector resolutions across sigmoid CFD levels."""
    ROOT.gROOT.SetBatch(True)
    canvas = ROOT.TCanvas("c_cfd_fraction_resolution", "CFD timing comparison", 1000, 720)
    canvas.SetGrid()
    multigraph = ROOT.TMultiGraph(
        "cfd_fraction_resolution",
        ";Sigmoid constant fraction [%];Weighted-core time resolution [ns]",
    )
    colors = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2,
              ROOT.kMagenta + 1, ROOT.kOrange + 7]
    graphs = []
    legend = ROOT.TLegend(0.72, 0.68, 0.89, 0.88)
    any_points = False
    for detector_index, detector in enumerate(active_detectors):
        graph = ROOT.TGraphErrors()
        graph.SetName(f"g_cfd_resolution_{detector}")
        graph.SetTitle(detector)
        graph.SetMarkerStyle(20 + detector_index)
        graph.SetMarkerColor(colors[detector_index % len(colors)])
        graph.SetLineColor(colors[detector_index % len(colors)])
        point = 0
        for method_index, method in enumerate(TIMING_METHODS):
            result = solutions_by_method[method]["mean"]["weighted_core"]
            sigma = float(result["resolutions_ns"][detector_index])
            error = float(result["resolution_uncertainties_ns"][detector_index])
            if not math.isfinite(sigma) or not math.isfinite(error):
                continue
            fraction = 20 + 10 * method_index
            graph.SetPoint(point, fraction, sigma)
            graph.SetPointError(point, 0.0, error)
            point += 1
        if point:
            any_points = True
            multigraph.Add(graph, "LP")
            legend.AddEntry(graph, detector, "lp")
            graphs.append(graph)
    if not any_points:
        canvas.Close()
        raise RuntimeError("no finite weighted-core resolutions for CFD comparison plot")
    multigraph.Draw("A")
    multigraph.GetXaxis().SetLimits(15.0, 55.0)
    multigraph.SetMinimum(0.0)
    legend.Draw()
    canvas.SaveAs(str(output / "cfd_fraction_resolution.png"))
    canvas.SaveAs(str(output / "cfd_fraction_resolution.pdf"))
    canvas.Close()


def compare_with_single_pad(
    source: Path, cluster_records: list[dict[str, object]],
    active_detectors: tuple[str, ...],
) -> dict[str, object]:
    """Separate changes in the estimator from changes in the selected cohort."""
    pad_records, _, _ = read_events(source, 1)
    common = ({int(row["event"]) for row in pad_records}
              & {int(row["event"]) for row in cluster_records})
    cohorts = {
        "single_pad_all": pad_records,
        "cluster_all": cluster_records,
        "single_pad_common": [row for row in pad_records if row["event"] in common],
        "cluster_common": [row for row in cluster_records if row["event"] in common],
    }
    result = {"common_events": len(common), "cohorts": {}}
    for name, records in cohorts.items():
        cohort = {"events": len(records), "estimators": {}}
        for method in TIMING_METHODS:
            cohort[method] = {}
            for estimator in ("mean", "min"):
                try:
                    pairs, raw, core, weighted = analyze(
                        records, estimator, active_detectors, method
                    )
                    cohort[method][estimator] = {
                        "pairs": pairs, "raw": raw, "core": core,
                        "weighted_core": weighted,
                    }
                except RuntimeError as exc:
                    cohort[method][estimator] = {"unavailable": str(exc)}
        result["cohorts"][name] = cohort
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pad_root", type=Path, help="Pad-step3-rec.root")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--bm-time-source", choices=("pad", "cluster"), default="pad",
        help="pad: existing window timing; cluster: entire uniquely matched cluster",
    )
    parser.add_argument(
        "--bm2-matched-pad-count",
        type=int,
        default=0,
        help=(
            "require exactly this many valid BM2 pads in the track-centered "
            "3x3 matching window; 0 disables the multiplicity cut"
        ),
    )
    parser.add_argument(
        "--detectors",
        default=",".join(DETECTORS),
        help=(
            "comma-separated detectors used in the variance solution; "
            "choices are " + ",".join(DETECTORS)
        ),
    )
    args = parser.parse_args()
    if args.bm2_matched_pad_count < 0:
        parser.error("--bm2-matched-pad-count must be >= 0")
    if args.bm_time_source == "cluster" and args.bm2_matched_pad_count != 0:
        parser.error("cluster mode requires --bm2-matched-pad-count 0")
    active_detectors = tuple(
        name.strip() for name in args.detectors.split(",") if name.strip()
    )
    if (
        len(active_detectors) < 3
        or len(set(active_detectors)) != len(active_detectors)
        or any(name not in DETECTORS for name in active_detectors)
    ):
        parser.error("--detectors requires at least three unique known detectors")
    if "BM2" not in active_detectors:
        parser.error("--detectors must include BM2 for the BM2-Tracker plots")
    source = args.pad_root.resolve()
    default_output = source.parent / "time_resolution"
    if args.bm_time_source == "cluster":
        default_output /= "cluster"
    output = (args.output_dir or default_output).resolve()
    if args.bm_time_source == "cluster" and output == (source.parent / "time_resolution").resolve():
        parser.error("cluster output must use a separate directory from single-pad results")
    output.mkdir(parents=True, exist_ok=True)

    records, availability, total_entries = read_events(
        source, args.bm2_matched_pad_count, args.bm_time_source
    )
    if not records:
        raise RuntimeError(f"BM2 {args.bm_time_source} timing selection accepted no events")
    active_indices = [DETECTORS.index(name) for name in active_detectors]
    solutions_by_method: dict[str, dict[str, object]] = {}
    gaussian_rows_by_method: dict[str, list[dict[str, object]]] = {}
    common_events_by_method: dict[str, int] = {}
    pair_rows_all: list[dict[str, object]] = []
    detector_rows_all: list[tuple[object, ...]] = []
    method_reports: dict[str, dict[str, object]] = {}
    for method in TIMING_METHODS:
        method_output = output if method == "t20" else output / method
        method_output.mkdir(parents=True, exist_ok=True)
        common_events_by_method[method] = sum(
            all(math.isfinite(float(record[
                "mean" if method == "t20" else f"mean_{method}"
            ][index])) for index in active_indices)
            for record in records
        )
        pair_rows: list[dict[str, object]] = []
        solutions: dict[str, object] = {}
        for estimator in ("mean", "min"):
            rows, raw, core, weighted_core = analyze(
                records, estimator, active_detectors, method
            )
            pair_rows.extend(rows)
            solutions[estimator] = {
                "raw": raw, "core": core, "weighted_core": weighted_core
            }
        solutions_by_method[method] = solutions
        selection_label = (
            f"BM2 cluster timing, sigmoid {method}"
            if args.bm_time_source == "cluster"
            else f"BM2 track-matched pad timing, sigmoid {method}"
        )
        gaussian_rows = plot_bm2_tracker_gaussians(
            records, method_output, args.bm2_matched_pad_count,
            selection_label, method,
        )
        gaussian_rows_by_method[method] = gaussian_rows
        pair_rows_all.extend(pair_rows)

        for estimator in ("mean", "min"):
            for distribution in ("raw", "core", "weighted_core"):
                result = solutions[estimator][distribution]
                for index, detector in enumerate(active_detectors):
                    detector_rows_all.append(
                        (
                            method,
                            estimator,
                            distribution,
                            detector,
                            result["variances_ns2"][index],
                            result["resolutions_ns"][index],
                            result.get(
                                "resolution_uncertainties_ns",
                                [math.nan] * len(active_detectors),
                            )[index],
                        )
                    )
        with (method_output / "pairwise_time_resolution.tsv").open(
            "w", newline="", encoding="utf-8"
        ) as stream:
            fields = [
                "estimator", "timing_method", "first", "second", "events",
                "mean_ns", "rms_ns", "core_events", "core_mean_ns",
                "core_sigma_ns",
            ]
            writer = csv.DictWriter(
                stream, fieldnames=fields, delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(pair_rows)

        method_report: dict[str, object] = {
            "timing_method": method,
            "fraction": int(method[1:]),
            "bm_time_branches": [
                f"{('bmClusterS' if args.bm_time_source == 'cluster' else 'bmS')}{'' if method == 't20' else method[1:]}_{field}"
                for field in ("mean", "min", "n")
            ],
            "tracker_time_branches": [
                f"trS{'' if method == 't20' else method[1:]}_{field}"
                for field in ("mean", "min", "n")
            ],
            "selected_track_events": len(records),
            "common_active_detector_events": common_events_by_method[method],
            "detector_availability": availability,
            "active_detectors": active_detectors,
            "core_definition": "iterative 3-sigma clipping initialized by median/MAD",
            "model": "Var(t_i-t_j)=sigma_i^2+sigma_j^2; unweighted LS over pair-specific cohorts",
            "bm2_tracker_gaussian_fits": gaussian_rows,
            "solutions": solutions,
        }
        method_reports[method] = method_report
        with (method_output / "time_resolution.json").open(
            "w", encoding="utf-8"
        ) as stream:
            json.dump(serializable(method_report), stream, indent=2, ensure_ascii=False)
            stream.write("\n")

    with (output / "pairwise_time_resolution.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        fields = [
            "estimator", "timing_method", "first", "second", "events",
            "mean_ns", "rms_ns", "core_events", "core_mean_ns",
            "core_sigma_ns",
        ]
        writer = csv.DictWriter(
            stream, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(pair_rows_all)

    with (output / "detector_time_resolution.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("timing_method", "estimator", "distribution", "detector",
                         "variance_ns2", "resolution_ns", "resolution_uncertainty_ns"))
        writer.writerows(detector_rows_all)

    with (output / "cfd_fraction_resolution.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("timing_method", "fraction_percent", "estimator",
                         "distribution", "detector", "variance_ns2",
                         "resolution_ns", "resolution_uncertainty_ns"))
        for method in TIMING_METHODS:
            for row in detector_rows_all:
                if row[0] == method:
                    writer.writerow((method, int(method[1:]), *row[1:]))
    plot_cfd_fraction_resolution(solutions_by_method, active_detectors, output)

    t20_solutions = solutions_by_method["t20"]
    t20_gaussian_rows = gaussian_rows_by_method["t20"]
    report: dict[str, object] = {
        "input": str(source),
        "tree": "TrackArrivalEvents",
        "timing": "sigmoid t20 (legacy primary); parallel CFD fractions t30/t40/t50",
        "timing_methods": method_reports,
        "sample_ns": SAMPLE_NS,
        "tracker_time": "strip-count-weighted mean over available matched X/Y axes",
        "total_track_events": total_entries,
        "selected_track_events": len(records),
        "timing_selection": {
            "bm2_matched_pad_count": args.bm2_matched_pad_count,
            "branch": "bmS_n[1]",
            "scope": "track-centered 3x3 pad matching window",
            "pad_time_requirement": "valid sigmoid fit; same event cohort used for t20/t30/t40/t50",
            "disabled_when": 0,
        },
        "detector_availability": availability,
        "active_detectors": active_detectors,
        "common_events": common_events_by_method["t20"],
        "core_definition": "iterative 3-sigma clipping initialized by median/MAD",
        "model": "Var(t_i-t_j)=sigma_i^2+sigma_j^2; unweighted LS over pair-specific cohorts",
        "bm2_tracker_gaussian_fits": t20_gaussian_rows,
        "solutions": t20_solutions,
        "cfd_fraction_comparison": {
            "table": "cfd_fraction_resolution.tsv",
            "plot": "cfd_fraction_resolution.png",
        },
    }
    if args.bm_time_source == "cluster":
        report["timing_selection"] = {
            "bm_time_source": "cluster",
            "matched_cluster_count": 1,
            "branch": "bm_matched_cluster_count[1]",
            "time_branches": [
                "bmClusterS_mean", "bmClusterS30_mean", "bmClusterS40_mean",
                "bmClusterS50_mean",
            ],
            "association": "any cluster member intersects the track-centered 3x3 pad window",
            "ambiguity": "count all spatial candidates, including clusters without valid time",
            "time_scope": "all unique pads of the matched cluster, including outside the window",
            "pad_time_requirement": "valid sigmoid fit; same event cohort used for t20/t30/t40/t50",
            "mean_weighting": "equal weight per valid pad; largest maxA waveform per duplicate pad",
            "unmatched_clusters": "excluded from event time",
        }
        report["single_pad_comparison"] = compare_with_single_pad(
            source, records, active_detectors
        )
    with (output / "time_resolution.json").open("w", encoding="utf-8") as stream:
        json.dump(serializable(report), stream, indent=2, ensure_ascii=False)
        stream.write("\n")

    print(f"input entries: {total_entries}")
    print(
        "selected entries: "
        f"{len(records)} (BM time source={args.bm_time_source}, "
        f"pad count={'uncut' if args.bm2_matched_pad_count == 0 else args.bm2_matched_pad_count})"
    )
    print("availability: " + ", ".join(f"{name}={availability[name]}" for name in DETECTORS))
    print(
        "active-detector common cohort "
        f"({','.join(active_detectors)}): {common_events_by_method['t20']}"
    )
    for method in TIMING_METHODS:
        print(f"sigmoid {method} weighted-core resolution [ns]:")
        result = solutions_by_method[method]["mean"]["weighted_core"]
        for index, (name, resolution) in enumerate(
            zip(active_detectors, result["resolutions_ns"])
        ):
            uncertainty = result["resolution_uncertainties_ns"][index]
            print(f"  {name}: {resolution:.6g} +/- {uncertainty:.3g}"
                  if math.isfinite(resolution) else f"  {name}: invalid")
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)

#!/usr/bin/env python3
"""Estimate per-detector timing resolution from TrackArrivalEvents.

The stored sigmoid t20 values are in APV samples (25 ns/sample).  A Tracker
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
    path: Path, bm2_matched_pad_count: int = 0
) -> tuple[list[dict[str, object]], dict[str, int], int]:
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"cannot open {path}")
    tree = root_file.Get("TrackArrivalEvents")
    if not tree:
        raise RuntimeError(f"TrackArrivalEvents is missing in {path}")

    event = array("i", [0])
    bm_mean = array("d", [0.0] * 2)
    bm_min = array("d", [0.0] * 2)
    bm_count = array("i", [0] * 2)
    tracker_mean = array("d", [0.0] * 6)
    tracker_min = array("d", [0.0] * 6)
    tracker_count = array("i", [0] * 6)
    for branch, target in (
        ("event", event),
        ("bmS_mean", bm_mean),
        ("bmS_min", bm_min),
        ("bmS_n", bm_count),
        ("trS_mean", tracker_mean),
        ("trS_min", tracker_min),
        ("trS_n", tracker_count),
    ):
        if not tree.GetBranch(branch):
            raise RuntimeError(f"required branch {branch} is missing")
        tree.SetBranchAddress(branch, target)

    availability = {name: 0 for name in DETECTORS}
    records: list[dict[str, object]] = []
    total_entries = int(tree.GetEntries())
    for entry in range(total_entries):
        tree.GetEntry(entry)
        if (
            bm2_matched_pad_count > 0
            and int(bm_count[1]) != bm2_matched_pad_count
        ):
            continue
        mean_times = [float(bm_mean[0]), float(bm_mean[1])]
        min_times = [float(bm_min[0]), float(bm_min[1])]
        detector_valid = [
            bm_count[0] > 0 and valid_time(mean_times[0]) and valid_time(min_times[0]),
            bm_count[1] > 0 and valid_time(mean_times[1]) and valid_time(min_times[1]),
        ]
        for tracker in range(3):
            x_index = 2 * tracker
            y_index = x_index + 1
            nx = int(tracker_count[x_index])
            ny = int(tracker_count[y_index])
            x_valid = (
                nx > 0
                and valid_time(float(tracker_mean[x_index]))
                and valid_time(float(tracker_min[x_index]))
            )
            y_valid = (
                ny > 0
                and valid_time(float(tracker_mean[y_index]))
                and valid_time(float(tracker_min[y_index]))
            )
            tracker_valid = x_valid or y_valid
            detector_valid.append(tracker_valid)
            if tracker_valid:
                valid_nx = nx if x_valid else 0
                valid_ny = ny if y_valid else 0
                mean_times.append(
                    (
                        valid_nx * float(tracker_mean[x_index])
                        + valid_ny * float(tracker_mean[y_index])
                    )
                    / (valid_nx + valid_ny)
                )
                minima = []
                if x_valid:
                    minima.append(float(tracker_min[x_index]))
                if y_valid:
                    minima.append(float(tracker_min[y_index]))
                min_times.append(min(minima))
            else:
                mean_times.append(math.nan)
                min_times.append(math.nan)
        for name, available in zip(DETECTORS, detector_valid):
            availability[name] += int(available)
        for index, available in enumerate(detector_valid):
            if not available:
                mean_times[index] = math.nan
                min_times[index] = math.nan
        records.append(
            {
                "event": int(event[0]),
                "mean": [value * SAMPLE_NS for value in mean_times],
                "min": [value * SAMPLE_NS for value in min_times],
            }
        )
    tree.ResetBranchAddresses()
    root_file.Close()
    return records, availability, total_entries


def analyze(
    records: list[dict[str, object]], estimator: str,
    active_detectors: tuple[str, ...],
) -> tuple[
    list[dict[str, object]], dict[str, object], dict[str, object], dict[str, object]
]:
    rows: list[dict[str, object]] = []
    detector_indices = [DETECTORS.index(name) for name in active_detectors]
    for first_local, first in enumerate(detector_indices):
        for second_local in range(first_local + 1, len(detector_indices)):
            second = detector_indices[second_local]
            differences = [
                float(record[estimator][first]) - float(record[estimator][second])
                for record in records
                if math.isfinite(float(record[estimator][first]))
                and math.isfinite(float(record[estimator][second]))
            ]
            if len(differences) < 3:
                raise RuntimeError(f"too few events for {DETECTORS[first]}-{DETECTORS[second]}")
            raw_mean, raw_sigma = mean_std(differences)
            core = clipped_stats(differences)
            rows.append(
                {
                    "estimator": estimator,
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
    records: list[dict[str, object]], output: Path, bm2_matched_pad_count: int
) -> list[dict[str, object]]:
    """Fit and plot mean-t20 differences for BM2 versus each Tracker."""
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.TH1.AddDirectory(False)
    pairs = ((1, 2), (1, 3), (1, 4))
    rows: list[dict[str, object]] = []
    objects: list[tuple[object, object, object, object]] = []
    selection_label = (
        f"BM2 matched pads = {bm2_matched_pad_count}"
        if bm2_matched_pad_count > 0
        else "BM2 matched pad multiplicity uncut"
    )

    for first, second in pairs:
        first_name = DETECTORS[first]
        second_name = DETECTORS[second]
        differences = [
            float(record["mean"][first]) - float(record["mean"][second])
            for record in records
            if math.isfinite(float(record["mean"][first]))
            and math.isfinite(float(record["mean"][second]))
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
        display_low = core_mean - 6.0 * core_sigma
        display_high = core_mean + 6.0 * core_sigma
        fit_low = core_mean - 3.0 * core_sigma
        fit_high = core_mean + 3.0 * core_sigma
        histogram = ROOT.TH1D(
            f"h_{key}",
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

        gaussian = ROOT.TF1(f"gaus_{key}", "gaus", fit_low, fit_high)
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

        canvas = ROOT.TCanvas(f"c_{key}", "time difference", 1000, 720)
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
            "mean sigmoid t20; 25 ns/sample; pair-specific valid cohorts; "
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pad_root", type=Path, help="Pad-step3-rec.root")
    parser.add_argument("--output-dir", type=Path)
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
    output = (args.output_dir or source.parent / "time_resolution").resolve()
    output.mkdir(parents=True, exist_ok=True)

    records, availability, total_entries = read_events(
        source, args.bm2_matched_pad_count
    )
    if not records:
        raise RuntimeError("BM2 matched-pad cut selected no events")
    active_indices = [DETECTORS.index(name) for name in active_detectors]
    common_events = sum(
        all(math.isfinite(float(record["mean"][index])) for index in active_indices)
        for record in records
    )
    pair_rows: list[dict[str, object]] = []
    solutions: dict[str, object] = {}
    for estimator in ("mean", "min"):
        rows, raw, core, weighted_core = analyze(
            records, estimator, active_detectors
        )
        pair_rows.extend(rows)
        solutions[estimator] = {"raw": raw, "core": core, "weighted_core": weighted_core}

    gaussian_rows = plot_bm2_tracker_gaussians(
        records, output, args.bm2_matched_pad_count
    )

    with (output / "pairwise_time_resolution.tsv").open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "estimator", "first", "second", "events", "mean_ns", "rms_ns",
            "core_events", "core_mean_ns", "core_sigma_ns",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(pair_rows)

    with (output / "detector_time_resolution.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("estimator", "distribution", "detector", "variance_ns2",
                         "resolution_ns", "resolution_uncertainty_ns"))
        for estimator in ("mean", "min"):
            for distribution in ("raw", "core", "weighted_core"):
                result = solutions[estimator][distribution]
                for index, detector in enumerate(active_detectors):
                    writer.writerow(
                        (
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

    report = {
        "input": str(source),
        "tree": "TrackArrivalEvents",
        "timing": "sigmoid t20",
        "sample_ns": SAMPLE_NS,
        "tracker_time": "strip-count-weighted mean over available matched X/Y axes",
        "total_track_events": total_entries,
        "selected_track_events": len(records),
        "timing_selection": {
            "bm2_matched_pad_count": args.bm2_matched_pad_count,
            "branch": "bmS_n[1]",
            "scope": "track-centered 3x3 pad matching window",
            "pad_time_requirement": "valid sigmoid t20",
            "disabled_when": 0,
        },
        "detector_availability": availability,
        "active_detectors": active_detectors,
        "common_events": common_events,
        "core_definition": "iterative 3-sigma clipping initialized by median/MAD",
        "model": "Var(t_i-t_j)=sigma_i^2+sigma_j^2; unweighted LS over pair-specific cohorts",
        "bm2_tracker_gaussian_fits": gaussian_rows,
        "solutions": solutions,
    }
    with (output / "time_resolution.json").open("w", encoding="utf-8") as stream:
        json.dump(serializable(report), stream, indent=2, ensure_ascii=False)
        stream.write("\n")

    print(f"input entries: {total_entries}")
    print(
        "selected entries: "
        f"{len(records)} (BM2 matched pad count "
        f"{'uncut' if args.bm2_matched_pad_count == 0 else args.bm2_matched_pad_count})"
    )
    print("availability: " + ", ".join(f"{name}={availability[name]}" for name in DETECTORS))
    print(
        "active-detector common cohort "
        f"({','.join(active_detectors)}): {common_events}"
    )
    for estimator in ("mean", "min"):
        print(f"{estimator} estimator, weighted-core per-detector resolution [ns]:")
        result = solutions[estimator]["weighted_core"]
        for index, (name, resolution) in enumerate(
            zip(active_detectors, result["resolutions_ns"])
        ):
            uncertainty = result["resolution_uncertainties_ns"][index]
            print(f"  {name}: {resolution:.6g} +/- {uncertainty:.3g}"
                  if math.isfinite(resolution) else f"  {name}: invalid")
        print(
            "  pair-variance model residual RMS [ns^2]: "
            f"{result['pair_variance_residual_rms_ns2']:.6g}"
        )
    print("BM2-Tracker Gaussian fits:")
    for row in gaussian_rows:
        print(
            f"  {row['first']}-{row['second']}: "
            f"mu={row['gaussian_mu_ns']:.6g} +/- "
            f"{row['gaussian_mu_error_ns']:.3g} ns, "
            f"sigma={row['gaussian_sigma_ns']:.6g} +/- "
            f"{row['gaussian_sigma_error_ns']:.3g} ns"
        )
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)

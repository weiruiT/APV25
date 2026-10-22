"""Exercise publication with small fixtures; no ROOT runtime or data is needed."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location(
    "export_run", Path(__file__).resolve().parents[1] / "scripts/export_run.py")
exporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)


class ExportRunTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="apv25_export_test_")
        self.addCleanup(self.temporary.cleanup)
        self.run = Path(self.temporary.name) / "run_001"
        self.execution_id = "execution_001"
        self.execution = self.run / "metadata/executions" / self.execution_id
        self.write("metadata/latest_signal.json", json.dumps({"execution_id": self.execution_id}))
        self.write("metadata/latest_noise.json", json.dumps({"execution_id": self.execution_id}))
        self.write(f"metadata/executions/{self.execution_id}/config/configCR.json", "{}")
        physics = self.write("Tracker/physics/physics.root", "COSMIC")
        noise = self.write("Tracker/noise/noise.root", "NOISE")
        calibration = self.write("Combine/noise/noise-global-parameters-1778.txt", "SIGMA")
        noise_root = self.write("Combine/noise/noise-global-parameters-1778.root", "CALIBRATION")
        core = [f"Combine/Tracker-step{i}-{suffix}.root" for i, suffix in
                ((1, "raw"), (2, "ana"), (3, "dec"), (4, "rec"))]
        core += ["Combine/Pad-step3-rec.root", "Combine/MMalignmentPar3.dat",
                 "Combine/TrackerXYResidual.png", "Combine/peaks_all.png"]
        timing = "Combine/time_resolution/"
        finals = [timing + name for name in (
            "detector_time_resolution.tsv", "pairwise_time_resolution.tsv",
            "time_resolution.json", "cfd_fraction_resolution.tsv",
            "cfd_fraction_resolution.png", "cfd_fraction_resolution.pdf",
            "t30/time_resolution.json", "t30/pairwise_time_resolution.tsv")]
        diagnostics = [timing + name for name in (
            "bm2_tracker_time_differences.root", "bm2_tracker_time_differences.png",
            "BM2_Tracker1_gaussian_fit.png", "bm2_tracker_gaussian_fits.tsv")]
        cluster = timing + "cluster/execution_001/"
        cluster_names = [cluster + "detector_time_resolution.tsv",
                         cluster + "time_resolution.json",
                         cluster + "t40/pairwise_time_resolution.tsv",
                         cluster + "BM2_Tracker1_gaussian_fit.png"]
        snr = timing + "snr_binned/execution_001/"
        snr_names = [snr + name for name in (
            "snr_bins.tsv", "snr_thresholds.tsv", "snr_bands.tsv",
            "snr_timing_report.json", "snr_threshold_scan.png", "snr_threshold_scan.pdf",
            "matched_bm2_events.tsv", "snr_timing_plots.root", "snr_timing_comparison.png",
            "pairwise_details/summary.tsv",
            "pairwise_details/SNR_15_25/detector_solution.tsv",
            "pairwise_details/SNR_15_25/solution.json",
            "pairwise_details/SNR_15_25/selected_events.tsv",
            "pairwise_details/SNR_15_25/time_differences.root")]
        outputs = {}
        for name in core + finals + diagnostics + cluster_names + snr_names:
            outputs[name] = exporter.sha256(self.write(name, name))
        noise_outputs = {str(p.relative_to(self.run)): exporter.sha256(p)
                         for p in (calibration, noise_root)}
        self.manifest = {
            "execution_id": self.execution_id, "status": "completed", "config": "config",
            "config_sha256": {"configCR.json": exporter.sha256(self.execution / "config/configCR.json")},
            "inputs": {"signal": [self.input_record(physics),
                                   self.input_record(calibration, "peak_noise_calibration")],
                       "noise": [self.input_record(noise)]},
            "stages": {
                "signal": {"status": "completed", "exit_code": 0,
                           "outputs": {n: outputs[n] for n in core + finals + diagnostics}},
                "noise": {"status": "completed", "exit_code": 0, "outputs": noise_outputs},
                "cluster_timing": {"status": "completed", "exit_code": 0,
                                   "outputs": {n: outputs[n] for n in cluster_names}},
                "snr_timing": {"status": "completed", "exit_code": 0,
                               "outputs": {n: outputs[n] for n in snr_names}},
            },
            "outputs": dict(outputs, **noise_outputs),
        }
        self.save_manifest()

    def write(self, name, content):
        path = self.run / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def input_record(self, path, role="acquisition"):
        return {"path": str(path), "run_path": str(path),
                "sha256": exporter.sha256(path), "role": role}

    def save_manifest(self):
        (self.execution / "manifest.json").write_text(json.dumps(self.manifest), encoding="utf-8")

    def test_only_final_timing_is_published_and_combine_is_preserved(self):
        original = {p: exporter.sha256(p) for p in (self.run / "Combine").rglob("*") if p.is_file()}
        # An unrecorded result must never slip into the published snapshot.
        self.write("Combine/time_resolution/t50/time_resolution.json", "UNRECORDED")
        target = exporter.export(self.run, "final_v1")
        timing = target / "timing/time_resolution"
        for name in ("detector_time_resolution.tsv", "time_resolution.json",
                     "cfd_fraction_resolution.png", "t30/time_resolution.json",
                     "cluster/execution_001/t40/pairwise_time_resolution.tsv",
                     "snr_binned/execution_001/snr_bins.tsv",
                     "snr_binned/execution_001/snr_threshold_scan.png",
                     "snr_binned/execution_001/detector_solutions/SNR_15_25/detector_solution.tsv"):
            self.assertTrue((timing / name).is_file(), name)
        self.assertFalse((timing / "t50/time_resolution.json").exists())
        self.assertFalse(list(target.rglob("*.root")))
        self.assertFalse(list(target.rglob("*gaussian*")))
        self.assertFalse(list(target.rglob("*events.tsv")))
        self.assertFalse(list(target.rglob("pairwise_details")))
        self.assertFalse((target / "analysis").exists())
        self.assertFalse((target / "noise").exists())
        self.assertTrue((target / "metadata/noise_manifest.json").is_file())
        self.assertIn("../../Combine/", (target / "README.md").read_text(encoding="utf-8"))
        self.assertEqual(original, {p: exporter.sha256(p) for p in original})
        details = json.loads((target / "metadata/analysis.json").read_text())
        self.assertTrue(details["timing_available"])
        for destination, source in details["published_files"].items():
            self.assertEqual(exporter.sha256(target / destination), self.manifest["outputs"][source])

    def test_historical_run_nested_snr_table_still_exports_flat(self):
        flat = "Combine/time_resolution/snr_binned/execution_001/pairwise_details/SNR_15_25/detector_solution.tsv"
        old = flat.replace("pairwise_details/SNR_", "pairwise_details/run_001/SNR_")
        self.manifest["outputs"].pop(flat)
        self.manifest["stages"]["snr_timing"]["outputs"].pop(flat)
        checksum = exporter.sha256(self.write(old, "HISTORICAL"))
        self.manifest["outputs"][old] = checksum
        self.manifest["stages"]["snr_timing"]["outputs"][old] = checksum
        self.save_manifest()
        target = exporter.export(self.run, "historical")
        destination = (target / "timing/time_resolution/snr_binned/execution_001/"
                       "detector_solutions/SNR_15_25/detector_solution.tsv")
        self.assertEqual(destination.read_text(), "HISTORICAL")
        self.assertFalse((target / "timing/time_resolution/snr_binned/execution_001/"
                          "detector_solutions/run_001").exists())

    def test_changed_calibration_is_still_rejected(self):
        self.write("Combine/noise/noise-global-parameters-1778.txt", "CHANGED")
        with self.assertRaisesRegex(ValueError, "input changed"):
            exporter.export(self.run, "invalid")
        self.assertFalse((self.run / "Results/invalid").exists())

    def test_extra_physics_input_is_rejected(self):
        self.write("Tracker/physics/extra.root", "EXTRA")
        with self.assertRaisesRegex(ValueError, "input file list changed"):
            exporter.export(self.run, "invalid")

    def test_changed_intermediate_is_rejected_even_though_not_copied(self):
        self.write("Combine/Tracker-step4-rec.root", "CHANGED")
        with self.assertRaisesRegex(ValueError, "output changed"):
            exporter.export(self.run, "invalid")

    def test_missing_final_table_prevents_publication(self):
        del self.manifest["outputs"]["Combine/time_resolution/detector_time_resolution.tsv"]
        self.save_manifest()
        with self.assertRaisesRegex(ValueError, "No verified final timing"):
            exporter.export(self.run, "invalid")

    def test_existing_result_is_never_overwritten(self):
        target = exporter.export(self.run, "final_v1")
        before = exporter.sha256(target / "metadata/result_sha256.txt")
        with self.assertRaisesRegex(ValueError, "Refusing to overwrite"):
            exporter.export(self.run, "final_v1")
        self.assertEqual(before, exporter.sha256(target / "metadata/result_sha256.txt"))

    def test_copy_corruption_cleans_up_unpublished_directory(self):
        original_copy = exporter.copy_file

        def corrupt_copy(source, destination):
            original_copy(source, destination)
            if source.name == "detector_time_resolution.tsv":
                destination.write_text("CORRUPTED")

        with patch.object(exporter, "copy_file", side_effect=corrupt_copy):
            with self.assertRaisesRegex(ValueError, "output changed during export"):
                exporter.export(self.run, "invalid")
        self.assertEqual(list((self.run / "Results").iterdir()), [])

    def test_legacy_acquisition_records_without_role_remain_supported(self):
        self.manifest["inputs"]["signal"] = self.manifest["inputs"]["signal"][:1]
        del self.manifest["inputs"]["signal"][0]["role"]
        self.save_manifest()
        self.assertTrue(exporter.export(self.run, "legacy").is_dir())

    def test_explicit_execution_does_not_follow_another_latest_pointer(self):
        self.write("metadata/latest_signal.json", '{"execution_id":"other"}')
        self.assertTrue(exporter.export(self.run, "exact", execution_id=self.execution_id).is_dir())


if __name__ == "__main__":
    unittest.main()

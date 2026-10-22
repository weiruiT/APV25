from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
spec = importlib.util.spec_from_file_location(
    "run_manager", Path(__file__).resolve().parents[1] / "scripts/run_manager.py")
manager = importlib.util.module_from_spec(spec)
spec.loader.exec_module(manager)


class RunManagerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="apv25_batch_test_")
        self.root = Path(self.temp.name)
        self.project = self.root / "code"
        (self.project / "config").mkdir(parents=True)
        (self.project / "config/configCR.json").write_text(
            '{"run":{"analysis_apv":true,"noise_analysis":true}}')
        (self.project / "config/map.txt").write_text("DEFAULT")
        self.project_patch = patch.object(manager, "PROJECT", self.project)
        self.project_patch.start()
        (self.project / "scripts").mkdir()
        fake_timing = self.project / "scripts/analyze_time_resolution.py"
        fake_timing.write_text(
            "import pathlib, sys\n"
            "output = pathlib.Path(sys.argv[sys.argv.index('--output-dir') + 1])\n"
            f"names = {manager.TIMING_FILES!r}\n"
            "for name in names:\n"
            "    path = output / pathlib.Path(name).name\n"
            "    path.parent.mkdir(parents=True, exist_ok=True)\n"
            "    path.write_text('TIMING')\n"
            "print('TIMING_ARGS=' + ' '.join(sys.argv[1:]))\n"
        )
        scripts = self.project / "scripts"
        (scripts / "analyze_bm_snr_timing.py").write_text(
            "import pathlib, sys\n"
            "output = pathlib.Path(sys.argv[sys.argv.index('--output-dir') + 1])\n"
            "output.mkdir(parents=True, exist_ok=False)\n"
            "for name in ('snr_bins.tsv','snr_timing_report.json',"
            "'snr_timing_comparison.png','snr_timing_plots.root'):\n"
            "    (output / name).write_text('SNR')\n"
            "print('SNR_ARGS=' + ' '.join(sys.argv[1:]))\n"
        )
        (scripts / "export_bm_snr_intermediates.py").write_text(
            "import pathlib, sys\n"
            "output = pathlib.Path(sys.argv[sys.argv.index('--output-dir') + 1])\n"
            "output.mkdir(parents=True, exist_ok=False)\n"
            "for name in ('README.md','summary.tsv','provenance.json'):\n"
            "    (output / name).write_text('DETAILS')\n"
            "print('DETAILS_ARGS=' + ' '.join(sys.argv[1:]))\n"
        )
        self.executable = self.root / "fake_tracker.py"
        self.executable.write_text('''#!/usr/bin/env python3
import pathlib, sys
run, mode = pathlib.Path(sys.argv[1]), sys.argv[3]
if (run / ("fail_" + mode)).exists(): sys.exit(7)
value = pathlib.Path("config/map.txt").read_text() + " " + mode
if mode == "signal":
    names = ["Tracker-step1-raw.root","Tracker-step2-ana.root","Tracker-step3-dec.root","Tracker-step4-rec.root","Pad-step3-rec.root"]
    output = run / "Combine"
else:
    names = ["Tracker-noise-raw.root","fit_parameters.txt","fit_parameters.root","Tracker-noise-channel-parameters.txt","BM-noise-channel-parameters.txt","noise-global-parameters-1778.txt","noise-global-parameters-1778.root"]
    output = run / "Combine/noise"
output.mkdir(parents=True, exist_ok=True)
for name in names: (output / name).write_text(value)
''')
        self.executable.chmod(0o755)

    def tearDown(self):
        self.project_patch.stop()
        self.temp.cleanup()

    def make_run(self, name, noise=True):
        run = self.root / name
        for directory in ("Tracker/physics", "Tracker/noise", "config"):
            (run / directory).mkdir(parents=True)
        (run / "Tracker/physics/input.root").write_text("COSMIC:" + name)
        if noise:
            (run / "Tracker/noise/noise.root").write_text("NOISE:" + name)
        (run / "config/map.txt").write_text(name)
        return run

    def execute(self, run, force=1, mode="all"):
        args = SimpleNamespace(force=force, mode=mode, alignment_source=None,
                               plots=False)
        return manager.execute_run(run, args, self.executable, "stable-code")

    def test_parallel_runs_isolate_noise_and_config(self):
        runs = [self.make_run("run_001"), self.make_run("run_002")]
        with ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(self.execute, runs))
        self.assertTrue(all(r["status"] == "completed" for r in results))
        for run in runs:
            self.assertEqual(
                (run / "Combine/noise/noise-global-parameters-1778.txt").read_text(),
                run.name + " noise")
            self.assertEqual((run / "Combine/Tracker-step2-ana.root").read_text(),
                             run.name + " signal")
            self.assertTrue(
                (run / "Combine/time_resolution/"
                 "bm2_tracker_time_differences.png").is_file())
        self.assertEqual((self.project / "config/map.txt").read_text(), "DEFAULT")

    def test_timing_cut_is_forwarded_from_run_config(self):
        run = self.make_run("run_single_pad")
        (run / "config/configCR.json").write_text(
            '{"timing_cuts":{"bm2_matched_pad_count":1,'
            '"detectors":["BM2","Tracker1","Tracker2","Tracker3"]}}')
        result = self.execute(run)
        timing_log = (
            run / "metadata/executions" / result["execution_id"] / "timing.log"
        ).read_text()
        self.assertIn("--bm2-matched-pad-count 1", timing_log)
        self.assertIn(
            "--detectors BM2,Tracker1,Tracker2,Tracker3", timing_log
        )

    def test_success_automatically_publishes_final_timing_only(self):
        run = self.make_run("run_results")
        first = self.execute(run)
        self.assertEqual(first["status"], "completed", first)
        published = run / "Results" / first["execution_id"]
        self.assertEqual(first["results"], str(published))
        self.assertTrue((published / "timing/time_resolution/detector_time_resolution.tsv").is_file())
        self.assertFalse(list(published.rglob("*.root")))
        self.assertFalse(list(published.rglob("*gaussian_fit*")))
        self.assertTrue((run / "Combine/Tracker-step1-raw.root").is_file())
        self.assertTrue((run / "Combine/time_resolution/BM2_Tracker1_gaussian_fit.png").is_file())
        second = self.execute(run, force=0)
        self.assertEqual(second["status"], "completed", second)
        self.assertNotEqual(first["results"], second["results"])
        self.assertTrue(published.is_dir())
        latest = json.loads((run / "metadata/latest_results.json").read_text())
        self.assertEqual(latest["execution_id"], second["execution_id"])

    def test_noise_only_does_not_publish_timing(self):
        run = self.make_run("run_noise_only")
        result = self.execute(run, mode="noise")
        self.assertEqual(result["status"], "completed", result)
        self.assertFalse((run / "metadata/latest_results.json").exists())
        self.assertFalse((run / "Results").exists())

    def test_signal_only_does_not_publish_previous_snr_results(self):
        run = self.make_run("run_snr_disabled")
        override = run / "config/configCR.json"
        override.write_text('{"snr_timing":{"enabled":true}}')
        first = self.execute(run)
        self.assertEqual(first["status"], "completed", first)
        self.assertTrue(list(Path(first["results"]).rglob("snr_bins.tsv")))
        override.write_text('{"snr_timing":{"enabled":false}}')
        second = self.execute(run, mode="signal")
        self.assertEqual(second["status"], "completed", second)
        self.assertFalse(list(Path(second["results"]).rglob("snr_bins.tsv")))
        self.assertTrue(list((run / "Combine").rglob("snr_bins.tsv")))

    def test_export_failure_is_reported_without_losing_analysis(self):
        run = self.make_run("run_export_failure")
        with patch.object(manager, "export_results", side_effect=OSError("disk full")):
            result = self.execute(run)
        self.assertEqual(result["status"], "partial_failed")
        self.assertIn("disk full", result["error"])
        self.assertEqual(result["stages"]["signal"]["status"], "completed")
        self.assertTrue((run / "metadata/latest_signal.json").is_file())
        self.assertFalse((run / "metadata/latest_results.json").exists())

    def test_cluster_timing_is_additive_and_reusable(self):
        run = self.make_run("run_cluster")
        (run / "config/configCR.json").write_text(
            '{"cluster_timing":{"enabled":true},'
            '"timing_cuts":{"bm2_matched_pad_count":1,'
            '"detectors":["BM2","Tracker1","Tracker2","Tracker3"]}}')
        result = self.execute(run)
        self.assertEqual(result["status"], "completed", result)
        stage = result["stages"]["cluster_timing"]
        self.assertEqual(stage["status"], "completed")
        execution = run / "metadata/executions" / result["execution_id"]
        log = (execution / "cluster_timing.log").read_text()
        self.assertIn("--bm-time-source cluster", log)
        self.assertNotIn("--bm2-matched-pad-count", log)
        self.assertIn("--bm2-matched-pad-count 1", (execution / "timing.log").read_text())
        self.assertTrue(all("time_resolution/cluster/" in name for name in stage["outputs"]))
        self.assertTrue((run / "Combine/time_resolution/time_resolution.json").is_file())
        self.assertTrue((run / "metadata/latest_cluster_timing.json").is_file())
        again = self.execute(run, force=0)
        self.assertEqual(again["status"], "completed", again)
        self.assertEqual(again["stages"]["cluster_timing"]["reused_execution"],
                         result["execution_id"])

    def test_cluster_timing_failure_preserves_single_pad_stage(self):
        run = self.make_run("run_cluster_failure")
        (run / "config/configCR.json").write_text('{"cluster_timing":{"enabled":true}}')
        original = manager.run_stage

        def fail_cluster(command, cwd, logfile):
            if "--bm-time-source" in command:
                logfile.write_text("Cluster analysis failure\n")
                return 9
            return original(command, cwd, logfile)

        with patch.object(manager, "run_stage", side_effect=fail_cluster):
            result = self.execute(run)
        self.assertEqual(result["status"], "partial_failed")
        self.assertEqual(result["stages"]["signal"]["status"], "completed")
        self.assertEqual(result["stages"]["cluster_timing"]["status"], "failed")
        self.assertTrue((run / "metadata/latest_signal.json").is_file())
        self.assertFalse((run / "metadata/latest_cluster_timing.json").exists())
        self.assertFalse((run / "metadata/latest_results.json").exists())
        self.assertFalse((run / "Results").exists())

    def test_snr_grouped_timing_is_a_separate_manifest_stage(self):
        run = self.make_run("run_snr")
        (run / "config/configCR.json").write_text(
            '{"snr_timing":{"enabled":true,"bin_edges":[15,25,40,60,90]}}')

        def fake_stage(command, _cwd, logfile):
            logfile.write_text("ARGS=" + " ".join(command[1:]) + "\n")
            executable = Path(command[0]) == self.executable
            script_name = Path(command[1]).name if not executable else ""
            if executable:
                mode = command[3]
                output = run / ("Combine/noise" if mode == "noise" else "Combine")
                names = (manager.expected_outputs(run, "signal") if mode == "signal" else [
                    output / "noise-global-parameters-1778.txt",
                    output / "noise-global-parameters-1778.root",
                ])
                for path in names:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(mode)
            elif script_name == "analyze_time_resolution.py":
                for relative in manager.TIMING_FILES:
                    path = run / "Combine" / relative
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text("TIMING")
            else:
                output = Path(command[command.index("--output-dir") + 1])
                output.mkdir(parents=True, exist_ok=False)
                names = (("snr_bins.tsv", "snr_timing_report.json",
                          "snr_timing_comparison.png", "snr_timing_plots.root")
                         if script_name == "analyze_bm_snr_timing.py"
                         else ("README.md", "summary.tsv", "provenance.json"))
                for name in names:
                    (output / name).write_text("SNR")
            return 0

        with patch.object(manager, "run_stage", side_effect=fake_stage):
            result = self.execute(run, mode="all")
        self.assertEqual(result["status"], "completed")
        stage = result["stages"]["snr_timing"]
        self.assertEqual(stage["status"], "completed")
        self.assertEqual(stage["bin_edges"], [15, 25, 40, 60, 90])
        self.assertTrue(any(Path(name).name == "README.md" and
                            "pairwise_details" in Path(name).parts
                            for name in stage["outputs"]))
        execution = run / "metadata/executions" / result["execution_id"]
        self.assertIn("--snr-bin-edges 15,25,40,60,90",
                      (execution / "snr_timing.log").read_text())
        self.assertIn("--snr-bin-edges 15,25,40,60,90",
                      (execution / "snr_intermediates.log").read_text())
        self.assertIn("--single-run run_snr",
                      (execution / "snr_intermediates.log").read_text())
        second = self.execute(run, force=0, mode="all")
        self.assertEqual(second["status"], "completed")
        self.assertEqual(second["stages"]["snr_timing"]["reused_execution"],
                         result["execution_id"])

    def test_missing_noise_blocks_signal_peak_selection(self):
        run = self.make_run("run_missing", noise=False)
        result = self.execute(run)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["stages"]["noise"]["status"], "failed")
        self.assertEqual(result["stages"]["signal"]["status"], "failed")
        self.assertFalse((run / "metadata/latest_signal.json").exists())
        self.assertFalse((run / "metadata/latest_noise.json").exists())

    def test_single_failed_stage_has_failed_status(self):
        run = self.make_run("run_noise_failure")
        (run / "fail_noise").touch()
        result = self.execute(run, mode="noise")
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["stages"]["noise"]["status"], "failed")

    def test_one_failure_does_not_stop_other_run(self):
        bad, good = self.make_run("run_bad"), self.make_run("run_good")
        (bad / "fail_signal").touch()
        with ThreadPoolExecutor(max_workers=2) as pool:
            results = list(pool.map(self.execute, [bad, good]))
        self.assertEqual(results[0]["status"], "partial_failed")
        self.assertEqual(results[1]["status"], "completed")

    def test_force_zero_requires_matching_fingerprint_and_hashes(self):
        run = self.make_run("run_cache")
        first = self.execute(run)
        second = self.execute(run, force=0)
        self.assertEqual(second["status"], "completed")
        self.assertEqual(second["stages"]["signal"]["reused_execution"],
                         first["execution_id"])
        (run / "Tracker/physics/input.root").write_text("CHANGED")
        third = self.execute(run, force=0)
        self.assertEqual(third["stages"]["signal"]["status"], "failed")
        self.assertEqual(third["stages"]["noise"]["status"], "completed")
        (run / "Combine/noise/noise-global-parameters-1778.root").write_text(
            "CORRUPTED")
        fourth = self.execute(run, force=0)
        self.assertEqual(fourth["stages"]["noise"]["status"], "failed")

    def test_lock_prevents_same_run_reentry(self):
        run = self.make_run("run_locked")
        with manager.run_lock(run):
            result = self.execute(run)
        self.assertEqual(result["status"], "failed")
        self.assertIn("already processing", result["error"])

    def test_json_override_deep_merges_defaults(self):
        run = self.make_run("run_override")
        (self.project / "config/noise-global-parameters-1778.root").write_text(
            "HISTORICAL")
        (run / "config/configCR.json").write_text(
            '{"run":{"noise_analysis":false}}')
        destination = self.root / "snapshot"
        manager.config_snapshot(run, destination)
        config = json.loads((destination / "configCR.json").read_text())
        self.assertEqual(config["run"],
                         {"analysis_apv": True, "noise_analysis": False})
        self.assertFalse(
            (destination / "noise-global-parameters-1778.root").exists())

    def test_parent_discovery_skips_empty_and_unmanaged_directories(self):
        empty = self.root / "run_empty"
        manager.init_run(SimpleNamespace(directory=str(empty), signal=[], noise=[],
                                         copy=False, notes=""))
        ready = self.root / "run_ready"
        source = self.root / "ready.root"
        source.write_text("READY")
        manager.init_run(SimpleNamespace(directory=str(ready), signal=[str(source)],
                                         noise=[], copy=False, notes=""))
        unmanaged = self.make_run("run_unmanaged")
        self.assertEqual(manager.discover([str(self.root)], "run_*"), [ready.resolve()])
        self.assertEqual(manager.discover([str(unmanaged)], "run_*"), [unmanaged.resolve()])

    def test_empty_init_then_safe_append(self):
        run = self.root / "run_import"
        args = SimpleNamespace(directory=str(run), signal=[], noise=[],
                               copy=False, notes="unpaired")
        manager.init_run(args)
        self.assertEqual(list((run / "Tracker/physics").iterdir()), [])
        source = self.root / "source.root"
        source.write_text("ORIGINAL")
        args.signal = [str(source)]
        manager.init_run(args)
        target = run / "Tracker/physics/source.root"
        self.assertTrue(target.is_symlink())
        self.assertEqual(source.read_text(), "ORIGINAL")
        extra = self.root / "extra.root"
        extra.write_text("EXTRA")
        args.signal = [str(extra), str(source)]
        with self.assertRaises(ValueError):
            manager.init_run(args)
        self.assertFalse((run / "Tracker/physics/extra.root").exists())

    def test_init_rejects_legacy_dat(self):
        run = self.root / "run_dat"
        source = self.root / "legacy.dat"
        source.write_text("NOT-A-SUPPORTED-INPUT")
        args = SimpleNamespace(directory=str(run), signal=[str(source)],
                               noise=[], copy=False, notes="")
        with self.assertRaisesRegex(ValueError, "not a supported input"):
            manager.init_run(args)
        self.assertFalse(run.exists())



if __name__ == "__main__":
    unittest.main()

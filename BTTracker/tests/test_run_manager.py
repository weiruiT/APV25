from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

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
            "    path = output.parent / name\n"
            "    path.parent.mkdir(parents=True, exist_ok=True)\n"
            "    path.write_text('TIMING')\n"
            "print('TIMING_ARGS=' + ' '.join(sys.argv[1:]))\n"
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

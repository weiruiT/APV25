"""ROOT integration checks for independent BM1Events/BM2Events.

Run after sourcing thisroot.sh: /usr/bin/python3 -m unittest discover
    -s BTTracker/tests -p test_pad_reconstruction.py
"""
from array import array
import importlib.util
from pathlib import Path
import tempfile
import unittest

import ROOT

timing_spec = importlib.util.spec_from_file_location(
    "pad_timing_analysis", Path(__file__).resolve().parents[1] / "scripts/analyze_time_resolution.py")
timing_analysis = importlib.util.module_from_spec(timing_spec)
timing_spec.loader.exec_module(timing_analysis)

class PadReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ROOT.gROOT.SetBatch(True)
        project = Path(__file__).resolve().parents[1]
        ROOT.gInterpreter.AddIncludePath(str(project / "include"))
        if not ROOT.gInterpreter.Declare(
                f'#include "{project / "src/PadReconstruction.C"}"'):
            raise RuntimeError("Could not load PadReconstruction.C")

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="apv25_pad_test_")
        self.directory = Path(self.temp.name)
        self.step2 = self.directory / "Tracker-step2-ana.root"
        self.step4 = self.directory / "Tracker-step4-rec.root"
        self.output = self.directory / "Pad-step3-rec.root"
        (self.directory / "MMalignmentParBM.dat").write_text("12 0 0\n13 0 0\n")

    def tearDown(self):
        self.temp.cleanup()

    def make_inputs(self, with_chip=True, extra_signals=(), track_x=0.0):
        with ROOT.TFile(str(self.step2), "RECREATE") as source:
            tree = ROOT.TTree("BMPadSignals", "Selected BM waveforms")
            integers = {name: array("i", [0]) for name in
                        ("event", "channel", "global_channel", "position", "raw_chip")}
            doubles = {name: array("d", [0.0]) for name in
                       ("ped", "pedrms", "maxA", "arrivalT", "arrivetime_a",
                        "arrivetime_e", "arrivetime_sigmoid",
                        "arrivalT_sigmoid30", "arrivalT_sigmoid40",
                        "arrivalT_sigmoid50")}
            timestamp = array("q", [0])
            wave = ROOT.std.vector("short")()
            for name, value in integers.items():
                if name != "raw_chip" or with_chip:
                    tree.Branch(name, value, f"{name}/I")
            for name, value in doubles.items():
                tree.Branch(name, value, f"{name}/D")
            tree.Branch("timestamp", timestamp, "timestamp/L")
            tree.Branch("wave", wave)

            def signal(event, chip, position, samples, time, stamp):
                for name, value in dict(event=event, raw_chip=chip, position=position,
                                        channel=0, global_channel=0).items():
                    integers[name][0] = value
                for name in ("ped", "pedrms"):
                    doubles[name][0] = 0.0
                doubles["maxA"][0] = max(abs(v) for v in samples)
                for name, value in zip(("arrivalT", "arrivetime_a", "arrivetime_e",
                                        "arrivetime_sigmoid"), time):
                    doubles[name][0] = value
                for name, offset in (("arrivalT_sigmoid30", 1.0),
                                     ("arrivalT_sigmoid40", 2.0),
                                     ("arrivalT_sigmoid50", 3.0)):
                    base = time[3]
                    doubles[name][0] = base + offset if base > -900 else base
                timestamp[0] = stamp
                wave.clear()
                for sample in samples:
                    wave.push_back(sample)
                tree.Fill()

            # Same pad in both planes, adjacency within each plane, and a
            # repeated BM1 record: no cross-plane charge or timing leakage.
            signal(1, 12, 60, [0, 10, 20], [10, 10, 11, 12], 100)
            signal(1, 12, 61, [0, 10, 30], [14, 14, 15, 16], 110)
            signal(1, 12, 60, [0, 10, 10], [18, 18, 19, 20], 105)
            signal(1, 13, 60, [0, 100, 200], [50, 50, 51, 52], 90)
            signal(1, 13, 71, [0, -30, -60], [54, 54, 55, 56], 95)
            signal(1, 15, 60, [0, 1000], [1, 1, 1, 1], 80)
            # Two diagonal pads do not cluster. Missing BM2 and invalid
            # timing must be cleared, including after the previous event.
            signal(2, 12, 0, [0, 10], [-999, -999, -999, 30], 200)
            signal(2, 12, 12, [0, 20], [-999, -999, -999, float("nan")], 210)
            signal(3, 13, 120, [0, 60], [38, 38, 39, 40], 300)
            # No selected pad: invalid mapping or <= 5 ADC must not create
            # a BM tree row, even if the event has a valid track.
            signal(4, 12, -1, [0, 100], [1, 1, 1, 1], 400)
            signal(4, 13, 60, [0, 5], [1, 1, 1, 1], 400)
            for row in extra_signals:
                signal(*row)
            tree.Write()

            # Empty tracker trees still exercise the full pad function.
            pulse_status = array("i", [-1])
            sigmoid_status = array("i", [0])
            strip_event = array("i", [0])
            strip_channel = array("i", [128])
            strip_times = {name: array("d", [-999.0]) for name in
                           ("arrivalT", "arrivalT_e", "arrivalT_sigmoid",
                            "arrivalT_sigmoid30", "arrivalT_sigmoid40",
                            "arrivalT_sigmoid50")}
            for layer in range(3):
                for direction in "XY":
                    strip = ROOT.TTree(f"Track{direction}{layer}", "Tracker signals")
                    strip.Branch("event", strip_event, "event/I")
                    strip.Branch("channel", strip_channel, "channel/I")
                    for name, value in strip_times.items():
                        strip.Branch(name, value, f"{name}/D")
                    strip.Branch("pulseFitStatus", pulse_status, "pulseFitStatus/I")
                    strip.Branch("sigmoidFitStatus", sigmoid_status, "sigmoidFitStatus/I")
                    if layer == 0 and direction == "X":
                        strip_event[0] = 1
                        for name, value in zip(strip_times, (23, 24, 25, 26, 27, 28)):
                            strip_times[name][0] = value
                        strip.Fill()
                    strip.Write()
            source.Write()

        with ROOT.TFile(str(self.step4), "RECREATE") as source:
            track = ROOT.TTree("Track", "Three-layer tracks")
            event = array("i", [0])
            timestamp = array("q", [0])
            useflag = array("b", [1])
            fit = array("d", [0.0, 0.0, 0.0])
            prediction = array("d", [128.0] * 3)
            track.Branch("event", event, "event/I")
            track.Branch("timestamp", timestamp, "timestamp/L")
            track.Branch("useflag", useflag, "useflag/O")
            xfit = array("d", [0.0, track_x, 0.0])
            track.Branch("XfitPars", xfit, "XfitPars[3]/D")
            track.Branch("YfitPars", fit, "YfitPars[3]/D")
            for name in ("predStripX", "predStripY"):
                track.Branch(name, prediction, f"{name}[3]/D")
            for number in range(1, 5):
                event[0] = number
                timestamp[0] = number * 1000
                useflag[0] = number != 2
                track.Fill()
            track.Write()

    def test_planes_have_independent_hits_clusters_charge_and_times(self):
        self.make_inputs()
        self.assertTrue(ROOT.ReconstructPadEvents(str(self.step2), str(self.step4),
                                                 str(self.output)))
        with ROOT.TFile(str(self.output), "READ") as source:
            self.assertFalse(source.Get("PadEvents"))
            bm1 = source.Get("BM1Events")
            bm2 = source.Get("BM2Events")
            self.assertEqual(bm1.GetEntries(), 2)
            self.assertEqual(bm2.GetEntries(), 2)
            self.assertEqual([int(row.event) for row in bm1], [1, 2])
            self.assertEqual([int(row.event) for row in bm2], [1, 3])
            self.assertEqual(bm1.GetBranch("event").GetEntries(), 2)
            bm1.GetEntry(0)
            bm2.GetEntry(0)
            self.assertEqual((bm1.timestamp, bm2.timestamp), (100, 90))
            self.assertEqual(list(bm1.hit_position), [60, 61])
            self.assertEqual(list(bm2.hit_position), [60, 71])
            self.assertEqual(list(bm1.hit_charge), [50, 40])
            self.assertEqual(list(bm2.hit_charge), [300, 90])
            self.assertEqual(list(bm1.hit_peak), [20, 30])
            self.assertEqual(list(bm2.hit_peak), [200, 60])
            self.assertEqual(list(bm1.cluster_size), [2])
            self.assertEqual(list(bm2.cluster_size), [2])
            self.assertEqual((bm1.primary_charge, bm2.primary_charge), (90, 390))
            self.assertAlmostEqual(bm1.primary_x, (50 * 49.5 + 40 * 58.5) / 90)
            self.assertAlmostEqual(bm2.primary_y, (300 * 49.5 + 90 * 58.5) / 390)
            self.assertEqual((bm1.arrivalT_sigmoid_mean, bm2.arrivalT_sigmoid_mean),
                             (16, 54))
            self.assertEqual((bm1.arrivalT_sigmoid_min, bm2.arrivalT_sigmoid_min),
                             (12, 52))
            self.assertEqual((bm1.n_arrivalT_sigmoid, bm2.n_arrivalT_sigmoid), (3, 2))
            self.assertEqual((bm1.arrivalT_sigmoid30_mean,
                              bm2.arrivalT_sigmoid30_mean), (17, 55))
            self.assertEqual((bm1.arrivalT_sigmoid40_mean,
                              bm2.arrivalT_sigmoid40_mean), (18, 56))
            self.assertEqual((bm1.arrivalT_sigmoid50_mean,
                              bm2.arrivalT_sigmoid50_mean), (19, 57))
            self.assertEqual(bm1.arrivetime_mean, bm1.arrivalT_a_mean)
            bm1.GetEntry(1)
            self.assertEqual(bm1.event, 2)
            self.assertEqual(bm1.nclusters, 2)
            self.assertEqual(bm1.primary_charge, 20)
            self.assertEqual(bm1.arrivalT_sigmoid_mean, 30)
            self.assertEqual(bm1.n_arrivalT_sigmoid, 1)
            self.assertEqual(bm1.arrivalT_mean, -999)
            bm2.GetEntry(1)
            self.assertEqual(bm2.event, 3)
            self.assertEqual(bm2.primary_charge, 60)

            times = source.Get("TrackArrivalEvents")
            self.assertEqual(times.GetEntries(), 3)
            times.GetEntry(0)
            self.assertEqual(list(times.bmS_n), [2, 2])
            self.assertEqual(list(times.bmS_mean), [14, 54])
            self.assertEqual(list(times.bmS30_mean), [15, 55])
            self.assertEqual(list(times.bmS40_mean), [16, 56])
            self.assertEqual(list(times.bmS50_mean), [17, 57])
            self.assertEqual(times.trS_mean[0][0], 25)
            self.assertEqual(times.trS30_mean[0][0], 26)
            self.assertEqual(times.trS40_mean[0][0], 27)
            self.assertEqual(times.trS50_mean[0][0], 28)
            self.assertEqual(list(times.bmClusterS_mean), [14, 54])
            self.assertEqual(list(times.bmClusterS30_mean), [15, 55])
            self.assertEqual(list(times.bmClusterS_n), [2, 2])
            self.assertEqual(list(times.bm_matched_cluster_count), [1, 1])
            self.assertEqual(list(times.bm_pred_ix), [5, 5])
            self.assertEqual(times.evS_n, 0)

    def test_cluster_time_includes_whole_cluster_and_excludes_other_clusters(self):
        self.make_inputs(extra_signals=(
            (1, 13, 82, [0, 30], [70] * 4, 110),
            (1, 13, 93, [0, 40], [80] * 4, 120),
            (1, 13, 0, [0, 500], [200] * 4, 130),
        ))
        self.assertTrue(ROOT.ReconstructPadEvents(str(self.step2), str(self.step4), str(self.output)))
        with ROOT.TFile(str(self.output), "READ") as source:
            bm2 = source.Get("BM2Events")
            bm2.GetEntry(0)
            self.assertEqual(list(bm2.hit_cluster_index), [0, 1, 1, 1, 1])
            self.assertEqual(list(bm2.cluster_size), [1, 4])
            self.assertEqual(list(bm2.cluster_time_sigmoid_n), [1, 4])
            self.assertEqual(list(bm2.cluster_time_sigmoid_mean), [200, 64.5])
            self.assertEqual(list(bm2.cluster_time_sigmoid30_mean), [201, 65.5])
            times = source.Get("TrackArrivalEvents")
            times.GetEntry(0)
            self.assertEqual(times.bmS_n[1], 2)
            self.assertEqual(times.bmS_mean[1], 54)  # legacy window unchanged
            self.assertEqual(times.bm_matched_cluster_count[1], 1)
            self.assertEqual(times.bm_matched_cluster_index[1], 1)
            self.assertEqual(times.bm_matched_cluster_size[1], 4)
            self.assertEqual(times.bmClusterS_n[1], 4)
            self.assertEqual(times.bmClusterS_mean[1], 64.5)
            self.assertEqual(times.bmClusterS30_mean[1], 65.5)
            self.assertEqual(times.bmClusterS_min[1], 52)
        cluster, _, total = timing_analysis.read_events(self.output, 0, "cluster")
        legacy, _, _ = timing_analysis.read_events(self.output, 1)
        self.assertEqual(total, 3)
        self.assertEqual([row["event"] for row in cluster], [1])
        self.assertEqual(cluster[0]["mean"][1], 64.5 * 25)
        self.assertNotIn(1, [row["event"] for row in legacy])
        with self.assertRaises(ValueError):
            timing_analysis.read_events(self.output, 1, "cluster")

    def test_ambiguous_spatial_clusters_are_rejected_even_with_one_valid_time(self):
        self.make_inputs(extra_signals=(
            (3, 13, 49, [0, 30], [20] * 4, 310),
            (3, 13, 61, [0, 40], [-999] * 4, 320),
        ))
        self.assertTrue(ROOT.ReconstructPadEvents(str(self.step2), str(self.step4), str(self.output)))
        with ROOT.TFile(str(self.output), "READ") as source:
            times = source.Get("TrackArrivalEvents")
            times.GetEntry(1)
            self.assertEqual(times.event, 3)
            self.assertEqual(times.bmS_n[1], 1)
            self.assertEqual(times.bm_matched_cluster_count[1], 2)
            self.assertEqual(times.bm_matched_cluster_index[1], -1)
            self.assertEqual(times.bmClusterS_n[1], 0)
            self.assertEqual(times.bmClusterS_mean[1], -999)
        cluster, _, _ = timing_analysis.read_events(self.output, 0, "cluster")
        legacy, _, _ = timing_analysis.read_events(self.output, 1)
        self.assertNotIn(3, [row["event"] for row in cluster])
        self.assertIn(3, [row["event"] for row in legacy])

    def test_prediction_outside_pad_plane_has_no_cluster_time(self):
        self.make_inputs(track_x=100.0)
        self.assertTrue(ROOT.ReconstructPadEvents(str(self.step2), str(self.step4), str(self.output)))
        with ROOT.TFile(str(self.output), "READ") as source:
            for row in source.Get("TrackArrivalEvents"):
                self.assertEqual(list(row.bm_matched_cluster_count), [0, 0])
                self.assertEqual(list(row.bmClusterS_n), [0, 0])
                self.assertEqual(list(row.bmClusterS_mean), [-999, -999])

    def test_disabled_bm_has_zero_tree_and_branch_entries(self):
        self.make_inputs()
        self.assertTrue(ROOT.ReconstructPadEvents(str(self.step2), str(self.step4),
                                                 str(self.output), False, True))
        with ROOT.TFile(str(self.output), "READ") as source:
            bm1 = source.Get("BM1Events")
            bm2 = source.Get("BM2Events")
            self.assertEqual(bm1.GetEntries(), 0)
            for branch in bm1.GetListOfBranches():
                self.assertEqual(branch.GetEntries(), 0, branch.GetName())
            self.assertEqual(bm2.GetEntries(), 2)
            for row in source.Get("TrackArrivalEvents"):
                self.assertEqual(row.bmS_n[0], 0)
                self.assertEqual(row.bmS_mean[0], -999)
                self.assertEqual(row.bm_matched_cluster_count[0], 0)
                self.assertEqual(row.bmClusterS_n[0], 0)
            self.assertIn("BM1 enabled=0; BM2 enabled=1",
                          source.Get("BMPadDetectorSelection").GetTitle())

    def test_both_disabled_still_write_empty_bm_trees(self):
        self.make_inputs()
        self.assertTrue(ROOT.ReconstructPadEvents(str(self.step2), str(self.step4),
                                                 str(self.output), False, False))
        with ROOT.TFile(str(self.output), "READ") as source:
            self.assertEqual(source.Get("BM1Events").GetEntries(), 0)
            self.assertEqual(source.Get("BM2Events").GetEntries(), 0)
            for row in source.Get("TrackArrivalEvents"):
                self.assertEqual(list(row.bmS_n), [0, 0])

    def test_missing_chip_is_rejected_without_creating_output(self):
        self.make_inputs(with_chip=False)
        self.assertFalse(ROOT.ReconstructPadEvents(str(self.step2), str(self.step4),
                                                  str(self.output)))
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()

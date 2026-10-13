#include "CRAnalysis.h"
#include "CRConfig.h"

#include "AnalysisAPV.h"
#include "DecodeAPV.h"
#include "DecodedAPVRoot2TrackRaw.h"
#include "PadReconstruction.h"
#include "ReconstructTrack.h"
#include "RunDataInput.h"

#include <functional>

namespace
{
bool RunStage(const char *label, const TString &output,
              const vector<TString> &trees, int force,
              const std::function<bool()> &run, bool nonempty = false)
{
  if (force != 1 && !gSystem->AccessPathName(output))
  {
    if (!ValidateRootTrees(output, trees, nonempty))
      return false;
    cout << " --> Reusing " << label << ": " << output << endl;
    return true;
  }
  cout << endl << " --> " << label << ", Output: " << output << endl;
  return run() && ValidateRootTrees(output, trees, nonempty);
}

bool LoadSignalSelection(const CRConfig &config, SignalSelectionConfig &cuts)
{
  cuts.amplitudeMin = config.GetDouble("amplitude_min", cuts.amplitudeMin);
  cuts.maxTMin = config.GetDouble("max_t_min", cuts.maxTMin);
  cuts.maxTMax = config.GetDouble("max_t_max", cuts.maxTMax);
  cuts.riseTimeMin = config.GetDouble("rise_time_min", cuts.riseTimeMin);
  cuts.riseTimeMax = config.GetDouble("rise_time_max", cuts.riseTimeMax);
  cuts.totMin = config.GetDouble("tot_min", cuts.totMin);
  cuts.totMax = config.GetDouble("tot_max", cuts.totMax);
  cuts.peakThresholdSigma =
      config.GetDouble("peak_threshold_sigma", cuts.peakThresholdSigma);
  cuts.requireSinglePeak = config.GetBool("require_single_peak", true);
  cuts.requireValidSigmoid = config.GetBool("require_valid_sigmoid", true);
  if (cuts.maxTMin >= cuts.maxTMax || cuts.riseTimeMin >= cuts.riseTimeMax ||
      cuts.totMin >= cuts.totMax || cuts.peakThresholdSigma <= 0.0)
  {
    cerr << "ERROR: invalid signal_cuts range in configCR.json." << endl;
    return false;
  }
  return true;
}
}

// A directory is exactly one run with Tracker/physics and Combine. The batch
// launcher invokes this once per run from that run's configuration snapshot.
bool ReadTrackerConfigured(TString fileDir, int force)
{
  CRConfig config;
  if (!config.Load())
    return false;

  TBenchmark bench;
  bench.Start("full");
  vector<TString> roots;
  TString prefix;
  const std::filesystem::path input(fileDir.Data());
  if (input.extension() == ".root")
  {
    if (!std::filesystem::is_regular_file(input))
    {
      cerr << "ERROR: input file does not exist: " << fileDir << endl;
      return false;
    }
    roots.push_back(fileDir);
    prefix = fileDir;
    prefix.Remove(prefix.Last('.'));
    prefix += "-";
  }
  else
  {
    if (!CollectRunInputs((input / "Tracker/physics").string(), roots))
    {
      cerr << "Pass one run directory; use run_batch.sh for a collection of runs."
           << endl;
      return false;
    }
    const TString combine = (input / "Combine").string();
    if (!EnsureOutputDirectory(combine))
      return false;
    prefix = combine + "/";
  }

  for (const TString &root : roots)
    if (!ValidateRootTrees(root, {"apv_EventSamples"}, true))
      return false;
  const TString raw = prefix + "Tracker-step1-raw.root";
  const TString ana = prefix + "Tracker-step2-ana.root";
  const TString dec = prefix + "Tracker-step3-dec.root";
  const TString rec = prefix + "Tracker-step4-rec.root";
  const TString pad = prefix + "Pad-step3-rec.root";
  vector<TString> anaTrees, decTrees;
  for (int layer = 0; layer < MMdetN; ++layer)
  {
    anaTrees.push_back(Form("TrackX%d", layer));
    anaTrees.push_back(Form("TrackY%d", layer));
    decTrees.push_back(Form("Track%d", layer));
  }
  anaTrees.push_back("BMPadSignals");

  cout << " --> Input: " << fileDir << endl;
  cout << " --> Input format: decoded ROOT (baseline-subtracted and 3-sigma selected)" << endl;

  const int bmSlot8LogicalAPV =
      static_cast<int>(config.GetDouble("bm_slot8_logical_chip", 12));
  const int bmSlot9LogicalAPV =
      static_cast<int>(config.GetDouble("bm_slot9_logical_chip", 13));
  cout << " --> BM slot mapping: slot8->chip" << bmSlot8LogicalAPV
       << ", slot9->chip" << bmSlot9LogicalAPV << endl;

  bool anyEnabled = false;
  if (config.GetBool("normalize_decoded_root", true))
  {
    anyEnabled = true;
    if (!RunStage("Normalising decoded ROOT", raw, {"TrackRaw"}, force, [&]() {
          DecodedAPVRoot2TrackRaw converter(
              roots, raw, 1, true, bmSlot8LogicalAPV, bmSlot9LogicalAPV);
          return converter.Init() && converter.Loop();
        }, true))
      return false;
  }
  if (config.GetBool("analysis_apv", true))
  {
    anyEnabled = true;
    if (gSystem->AccessPathName("config/map_APV.txt") ||
        !ValidateRootTrees(raw, {"TrackRaw"}, true))
    {
      cerr << "ERROR: waveform analysis requires TrackRaw and config/map_APV.txt."
           << endl;
      return false;
    }
    SignalSelectionConfig cuts;
    if (!LoadSignalSelection(config, cuts))
      return false;
    if (input.extension() == ".root")
    {
      cerr << "ERROR: peak selection now requires run-local noise calibration; "
              "pass the managed run directory instead of a standalone ROOT file."
           << endl;
      return false;
    }
    const std::filesystem::path peakNoiseParameters =
        input / "Combine/noise/noise-global-parameters-1778.txt";
    if (!std::filesystem::is_regular_file(peakNoiseParameters))
    {
      cerr << "ERROR: missing run-local peak-selection noise calibration: "
           << peakNoiseParameters
           << ". Run noise mode first, or use --mode all." << endl;
      return false;
    }
    if (!RunStage("Analyzing waveforms", ana, anaTrees, force, [&]() {
          AnalysisAPV analysis(raw, ana, 1, cuts, peakNoiseParameters.string());
          if (!analysis.Init())
            return false;
          analysis.Loop();
          return true;
        }))
      return false;
  }
  if (config.GetBool("decode_apv", true))
  {
    anyEnabled = true;
    if (!ValidateRootTrees(ana, anaTrees) ||
        !RunStage("Decoding hits", dec, decTrees, force, [&]() {
          DecodeAPV decoder(ana, dec, 1);
          if (!decoder.Init())
            return false;
          decoder.Loop();
          return true;
        }))
      return false;
  }
  if (config.GetBool("reconstruct_track", true))
  {
    anyEnabled = true;
    if (!ValidateRootTrees(dec, decTrees) ||
        !RunStage("Reconstructing three-layer tracks", rec, {"Track"}, force, [&]() {
          ReconstructTrack reconstruction(dec, rec, 1);
          if (!reconstruction.Init())
            return false;
          reconstruction.Loop();
          return true;
        }))
      return false;
  }
  if (config.GetBool("pad_reconstruction", true))
  {
    anyEnabled = true;
    if (!ValidateRootTrees(ana, {"BMPadSignals"}) ||
        !ValidateRootTrees(rec, {"Track"}) ||
        !RunStage("Reconstructing pad events and timing", pad,
                  {"BM1Events", "BM2Events", "TrackArrivalEvents"}, force, [&]() {
          return ReconstructPadEvents(ana, rec, pad,
                                      config.GetBool("bm1_enabled", true),
                                      config.GetBool("bm2_enabled", true));
        }))
      return false;
  }
  if (!anyEnabled)
  {
    cerr << "ERROR: all signal stages are disabled in configCR.json." << endl;
    return false;
  }
  cout << " --> Tracker analysis completed." << endl;
  bench.Show("full");
  return true;
}

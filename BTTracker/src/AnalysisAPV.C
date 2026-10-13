#include "AnalysisAPV.h"
#include "TFitResult.h"
#include "TNamed.h"
#include <algorithm>
#include <cmath>
#include <array>

namespace
{
const int kExpectedWaveSamples = 30;
const int kLocalPedestalMaxSamples = 6;
const double kMinLocalPedestalRms = 1.0;

enum PulsePeakFailureReason
{
  kPulsePeakOk = 0,
  kPulsePeakBadSize,
  kPulsePeakNoSampleOverThreshold,
  kPulsePeakNoLocalPeak,
  kPulsePeakMultiPeak,
  kPulsePeakEdgePeak,
  kPulsePeakLowAmplitude,
  kPulsePeakBadEdge,
  kPulsePeakUnknown,
  kPulsePeakFailureReasonCount
};

const char *PulsePeakFailureReasonName(int reason)
{
  switch (reason)
  {
  case kPulsePeakBadSize:
    return "bad_size";
  case kPulsePeakNoSampleOverThreshold:
    return "no_sample_over_peak_threshold";
  case kPulsePeakNoLocalPeak:
    return "no_local_peak";
  case kPulsePeakMultiPeak:
    return "multi_peak";
  case kPulsePeakEdgePeak:
    return "edge_peak";
  case kPulsePeakLowAmplitude:
    return "low_peak_amplitude";
  case kPulsePeakBadEdge:
    return "bad_edge";
  case kPulsePeakUnknown:
    return "unknown";
  default:
    return "ok";
  }
}

double PulseModel(double *x, double *par)
{
  const double b = par[3];
  if (b == 0)
    return par[0];

  const double term = (x[0] - par[2]) / b;
  return par[0] + par[1] * std::pow(term, 3) * std::exp(-term);
}

double SigmoidModel(double *x, double *par)
{
  const double k = par[3];
  if (k == 0.0)
    return par[0];
  return par[0] + par[1] / (1.0 + std::exp(-(x[0] - par[2]) / k));
}

double SigmoidFractionTime(double t50, double k, double fraction)
{
  return t50 + k * std::log(fraction / (1.0 - fraction));
}
}

AnalysisAPV::AnalysisAPV(TString fRawName, TString fAnaName, int f_force,
                         const SignalSelectionConfig &fSelection,
                         TString fNoiseParameterName)
{
  rawName = fRawName;
  anaName = fAnaName;
  force = f_force;
  selection = fSelection;
  noiseParameterName = fNoiseParameterName;
}

AnalysisAPV::~AnalysisAPV()
{
}

// Signal selection cut applied to every waveform (Tracker and BM pad)
// before it is written to step2:
//   maxA > amplitudeMin; timing ranges are open intervals. The decoded
//   signal ROOT has already had baseline subtraction and 3-sigma selection.
bool AnalysisAPV::Cut()
{
  return maxA > selection.amplitudeMin &&
         maxT > selection.maxTMin && maxT < selection.maxTMax &&
         riseT > selection.riseTimeMin && riseT < selection.riseTimeMax &&
         TOT > selection.totMin && TOT < selection.totMax;
}


bool AnalysisAPV::LoadBMMap()
{
  bmMap.clear();

  std::ifstream input("config/map_APV4BMM.txt");
  TString usedName = "config/map_APV4BMM.txt";

  AnalysisBMMapRow row;
  int connectorPin = -1;
  while (input >> row.board >> row.chip >> row.channelId
               >> connectorPin >> row.position)
  {
    if (row.board < 0 || row.board > 1 ||
        row.channelId < 1 || row.channelId > 128 ||
        row.position < 0 || row.position >= 121)
      continue;
    bmMap.push_back(row);
  }

  cout << " --> Loaded BM map " << usedName
       << " with " << bmMap.size() << " entries" << endl;
  return !bmMap.empty();
}



int AnalysisAPV::GlobalChannelForRaw(int rawBoard, int rawChip,
                                     int rawChannel, int &mappedChannel,
                                     int &padPosition)
{
  mappedChannel = -1;
  padPosition = -1;

  if (rawChip == 12 || rawChip == 13)
  {
    if (rawChannel < 0 || rawChannel >= 128 || bmMap.empty())
      return -1;

    const int padBoard = rawChip == 12 ? 0 : 1;
    const int mapChannelId = rawChannel + 1;
    for (const AnalysisBMMapRow &row : bmMap)
    {
      if (row.board == padBoard && row.channelId == mapChannelId)
      {
        padPosition = row.position;
        return 2 * MMdetN * 256 + padBoard * 121 + row.position + 1;
      }
    }
    return -1;
  }

  mappedChannel = IndexStrip(rawBoard, rawChip, rawChannel);
  if (mappedChannel < 0)
    return -1;

  const int trackerIndex = rawBoard * 2 + rawChip / 2;
  if (trackerIndex < 0 || trackerIndex >= 2 * MMdetN)
    return -1;
  return trackerIndex * 256 + mappedChannel + 1;
}

bool AnalysisAPV::LoadNoiseParameters()
{
  bmNoiseCalibrationLoaded = false;                                   //初始化噪声加载false
  channelCalibrationLoaded = false;
  std::vector<TString> candidates;
  if (noiseParameterName.IsNull() || noiseParameterName.Length() == 0)
  {
    cerr << "ERROR: peak selection requires the run-local 1778-channel "
            "noise parameter file." << endl;
    return false;
  }
  candidates.push_back(noiseParameterName);

  for (const TString &candidate : candidates)
  {
    std::ifstream input(candidate.Data());
    if (!input.is_open())
    {
      cerr << "ERROR: cannot open peak-selection noise calibration: "
           << candidate << endl;
      return false;
    }

    for (int i = 0; i < kGlobalNoiseChannels; ++i)
    {
      globalPedestal[i] = 0.0;
      globalRms[i] = kMinLocalPedestalRms;
    }
    for (int i = 0; i < 2 * MMdetN; ++i)
      for (int ch = 0; ch < 256; ++ch)
      {
        channelped[i][ch] = 0.0;
        channelrms[i][ch] = kMinLocalPedestalRms;
      }

    std::vector<bool> seen(kGlobalNoiseChannels, false);
    int count = 0;
    int globalChannel = -1;
    double mean = 0.0;
    double sigma = 0.0;
    std::string line;
    while (std::getline(input, line))
    {
      if (line.empty() || line[0] == '#')
        continue;

      std::istringstream iss(line);
      if (!(iss >> globalChannel >> mean >> sigma))
        continue;
      if (globalChannel < 1 || globalChannel > kGlobalNoiseChannels ||
          !std::isfinite(sigma))
        continue;

      const int index = globalChannel - 1;
      if (seen[index])
      {
        cerr << "ERROR: duplicate global channel " << globalChannel
             << " in " << candidate << endl;
        return false;
      }
      // Signal samples are already baseline-subtracted.  Import only sigma;
      // the noise pedestal mean must not be subtracted a second time.
      globalPedestal[index] = 0.0;
      globalRms[index] = std::max(std::fabs(sigma), kMinLocalPedestalRms);
      seen[index] = true;
      ++count;

      if (globalChannel <= 2 * MMdetN * 256)
      {
        const int treeIndex = index / 256;
        const int localChannel = index % 256;
        channelped[treeIndex][localChannel] = 0.0;
        channelrms[treeIndex][localChannel] = globalRms[index];
      }
    }

    if (count == kGlobalNoiseChannels)
    {
      cout << " --> Loaded peak-selection RMS from " << candidate
           << ": " << count << " global channels; pedestal means ignored" << endl;
      channelCalibrationLoaded = true;
      bmNoiseCalibrationLoaded = true;
      return true;
    }

    cout << " --> Ignore incomplete noise calibration " << candidate
         << ": " << count << " channels (need exactly "
         << kGlobalNoiseChannels << ")" << endl;
  }

  return false;
}

bool AnalysisAPV::Init()
{
  // Reading map: board chip channel connector strip.
  ifstream mapfile;
  mapfile.open("config/map_APV.txt");
  int mapin;
  while (mapfile >> mapin)
  {
    for (int kk = 0; kk < 3; kk++)
      mapfile >> mapin;
    mapfile >> mapin;
    APVmap.push_back(mapin);
  }

  frawFile = new TFile(rawName, "read");
  if (!frawFile->IsOpen())
  {
    cout << " ERROR! " << rawName << " cant open!" << endl;
    return false;
  }
  cout << " --> Analyzing waveform files from " << rawName << endl;
  frawTree = (TTree *)frawFile->Get("TrackRaw");
  if (frawTree == 0)
  {
    cout << " ERROR! TrackRaw Tree not exist in " << rawName << endl;
    return false;
  }
  frawTree->SetBranchAddress("event", &event, &b_event);
  frawTree->SetBranchAddress("board", &board, &b_board);
  frawTree->SetBranchAddress("chip", &chip, &b_chip);
  frawTree->SetBranchAddress("channel", &channel, &b_channel);
  frawTree->SetBranchAddress("wave", &wave, &b_wave);
  frawTree->SetBranchAddress("timestamp", &timestamp, &b_timestamp);

  for (int i = 0; i < 2 * MMdetN; ++i)
    for (int ch = 0; ch < 256; ++ch)
    {
      channelped[i][ch] = 0.0;
      channelrms[i][ch] = kMinLocalPedestalRms;
    }
  for (int i = 0; i < kGlobalNoiseChannels; ++i)
  {
    globalPedestal[i] = 0.0;
    globalRms[i] = kMinLocalPedestalRms;
  }

  LoadBMMap();
  if (!LoadNoiseParameters())
    return false;
  cout << " --> Preprocessed signal input: pedestal remains 0; run-local noise "
          "RMS is used only by significant-peak selection at "
       << selection.peakThresholdSigma << " sigma" << endl;

  if (force != 1)
  {
    FileStat_t fStat;
    gSystem->GetPathInfo(anaName, fStat);
    if (fStat.fSize != 0)
    {
      cout << " ERROR! " << anaName << " is exist!" << endl;
      return false;
    }
  }

  fanaFile = new TFile(anaName, "recreate");
  if (!fanaFile->IsOpen())
  {
    cout << "ERROR! " << anaName << " cant open!" << endl;
    return false;
  }

  for (int i = 0; i < 2 * MMdetN; i++)
  {
    fanaTree[i] = new TTree(
        Form("Track%s%d", XYname(i % 2).Data(), i / 2),
        Form("Track%s%d", XYname(i % 2).Data(), i / 2));
    fanaTree[i]->Branch("event", &event);
    fanaTree[i]->Branch("maxA", &maxA);
    fanaTree[i]->Branch("maxT", &maxT);
    fanaTree[i]->Branch("channel", &channel);
    fanaTree[i]->Branch("global_channel", &currentGlobalChannel);
    fanaTree[i]->Branch("ped", &ped);
    fanaTree[i]->Branch("pedrms", &pedrms);
    fanaTree[i]->Branch("peak_noise_rms", &peakNoiseRms);
    fanaTree[i]->Branch("peak_threshold_adc", &peakThresholdAdc);
    fanaTree[i]->Branch("riseT", &riseT);
    fanaTree[i]->Branch("TOT", &TOT);
    // Legacy Pause branches keep their schema, but remain invalid (-999,
    // pulseFitStatus=-1) in sigmoid-only mode. Do not relabel a as sigmoid.
    fanaTree[i]->Branch("arrivalT", &arrivalT);
    fanaTree[i]->Branch("arrivalT_a", &arrivalT_a);
    fanaTree[i]->Branch("arrivalT_e", &arrivalT_e);
    fanaTree[i]->Branch("arrivalT_sigmoid", &arrivalT_sigmoid);
    fanaTree[i]->Branch("arrivetime_sigmoid", &arrivalT_sigmoid);
    fanaTree[i]->Branch("sigmoid_t20", &arrivalT_sigmoid);
    fanaTree[i]->Branch("sigmoid_t50", &sigmoidT50);
    fanaTree[i]->Branch("sigmoid_B", &sigmoidBaseline);
    fanaTree[i]->Branch("sigmoid_H", &sigmoidHeight);
    fanaTree[i]->Branch("sigmoid_k", &sigmoidK);
    fanaTree[i]->Branch("sigmoidFitStatus", &sigmoidFitStatus);
    fanaTree[i]->Branch("sigmoidFitChi2NDF", &sigmoidFitChi2NDF);
    fanaTree[i]->Branch("sigmoid_fitStart", &sigmoidFitStart);
    fanaTree[i]->Branch("sigmoid_fitEnd", &sigmoidFitEnd);
    fanaTree[i]->Branch("arrivetime_a", &arrivalT_a);
    fanaTree[i]->Branch("arrivetime_e", &arrivalT_e);
    fanaTree[i]->Branch("pulseFitStatus", &pulseFitStatus);
    fanaTree[i]->Branch("pulseFitChi2NDF", &pulseFitChi2NDF);
    fanaTree[i]->Branch("n_significant_peaks", &pulsePeakCount);
    fanaTree[i]->Branch("wave", &waveAna);
    fanaTree[i]->Branch("timestamp", &timestamp, "timestamp/L");
  }

  bmTree = new TTree("BMPadSignals", "BM signals after configured waveform selection");
  bmTree->Branch("event", &event);
  bmTree->Branch("timestamp", &timestamp, "timestamp/L");
  bmTree->Branch("global_channel", &currentGlobalChannel);
  bmTree->Branch("raw_chip", &currentRawChip);
  bmTree->Branch("raw_channel", &currentRawChannel);
  bmTree->Branch("channel", &channel);
  bmTree->Branch("position", &currentPadPosition);
  bmTree->Branch("ped", &ped);
  bmTree->Branch("pedrms", &pedrms);
  bmTree->Branch("peak_noise_rms", &peakNoiseRms);
  bmTree->Branch("peak_threshold_adc", &peakThresholdAdc);
  bmTree->Branch("nisel", &nisel);
  bmTree->Branch("maxA", &maxA);
  bmTree->Branch("maxT", &maxT);
  bmTree->Branch("riseT", &riseT);
  bmTree->Branch("TOT", &TOT);
  // Same legacy Pause sentinel convention as the Tracker trees.
  bmTree->Branch("arrivalT", &arrivalT);
  bmTree->Branch("arrivalT_a", &arrivalT_a);
  bmTree->Branch("arrivalT_e", &arrivalT_e);
  bmTree->Branch("arrivalT_sigmoid", &arrivalT_sigmoid);
  bmTree->Branch("arrivetime_sigmoid", &arrivalT_sigmoid);
  bmTree->Branch("sigmoid_t20", &arrivalT_sigmoid);
  bmTree->Branch("sigmoid_t50", &sigmoidT50);
  bmTree->Branch("sigmoid_B", &sigmoidBaseline);
  bmTree->Branch("sigmoid_H", &sigmoidHeight);
  bmTree->Branch("sigmoid_k", &sigmoidK);
  bmTree->Branch("sigmoidFitStatus", &sigmoidFitStatus);
  bmTree->Branch("sigmoidFitChi2NDF", &sigmoidFitChi2NDF);
  bmTree->Branch("sigmoid_fitStart", &sigmoidFitStart);
  bmTree->Branch("sigmoid_fitEnd", &sigmoidFitEnd);
  bmTree->Branch("arrivetime_a", &arrivalT_a);
  bmTree->Branch("arrivetime_e", &arrivalT_e);
  bmTree->Branch("pulseFitStatus", &pulseFitStatus);
  bmTree->Branch("pulseFitChi2NDF", &pulseFitChi2NDF);
  bmTree->Branch("n_significant_peaks", &pulsePeakCount);
  bmTree->Branch("wave", &waveAna);

  return true;
}

void AnalysisAPV::Loop()
{
  // Signal pedestal remains zero.  Run-local RMS is kept in separate peak
  // selection fields so it cannot alter the existing timing or charge logic.

  ofstream op = LogFile(anaName);
  TCanvas *cTmp = new TCanvas("tmpana", "tmpana");
  TH1I *hChannel = new TH1I("Channel", "Channel", 256, 0, 256);
  hChannel->GetXaxis()->SetTitle("channel");
  hChannel->GetYaxis()->SetTitle("Counts");
  TGraph *gEvent[2 * MMdetN];
  TH1I *hHit[2 * MMdetN];
  for (int i = 0; i < 2 * MMdetN; i++)
  {
    hHit[i] = new TH1I(Form("Hit%s%d", XYname(i % 2).Data(), i / 2), Form("Hit%s%d", XYname(i % 2).Data(), i / 2), 256, 0, 256);
    hHit[i]->GetXaxis()->SetTitle("channel");
    hHit[i]->GetYaxis()->SetTitle("Counts");
  }

  /////////// begin to loop
  vector<int> eventvec[2 * MMdetN];
  vector<int> entryvec[2 * MMdetN];
  vector<int> valideventvec[2 * MMdetN];
  vector<int> allvalidevent;
  int RecordStart = -1;
  int RecordEnd = -1;
  int NumofTriggers = 0;
  int nentries = frawTree->GetEntriesFast();
  long long pulsePeakCheckedCount = 0;
  long long pulsePeakPassCount = 0;
  long long pulsePeakFailureCounts[kPulsePeakFailureReasonCount] = {0};
  // Count before the single-peak rejection, including waveforms which never
  // enter a signal tree. "All" means mapped, correctly sized raw waveforms.
  const char *peakGroups[] = {
      "All", "Tracker", "BM", "Tracker1", "Tracker2", "Tracker3",
      "BM1", "BM2", "Tracker1X", "Tracker1Y", "Tracker2X", "Tracker2Y",
      "Tracker3X", "Tracker3Y"};
  const int nPeakGroups = sizeof(peakGroups) / sizeof(peakGroups[0]);
  TH1D *peakAll[nPeakGroups], *peakAfterCut[nPeakGroups];
  for (int g = 0; g < nPeakGroups; ++g)
  {
    peakAll[g] = new TH1D(Form("PeakCount_all_%s", peakGroups[g]),
        Form("%s: mapped raw waveforms;Significant peak count;Waveforms", peakGroups[g]),
        31, -0.5, 30.5);
    peakAfterCut[g] = new TH1D(Form("PeakCount_afterCut_%s", peakGroups[g]),
        Form("%s: after basic Cut, before peak selection;Significant peak count;Waveforms", peakGroups[g]),
        31, -0.5, 30.5);
  }
  // A compact audit trail for the efficiency denominator; no rejected
  // waveform is promoted into the selected signal trees by this diagnostic.
  TTree peakDiagnostics("WaveformPeakDiagnostics",
      "Waveforms after basic Cut, BEFORE single-peak rejection");
  peakDiagnostics.SetDirectory(nullptr);
  Long64_t rawEntry = -1;
  int detectorId = -1; // 1..4 Tracker, 5=BM1, 6=BM2
  bool peakAccepted = false, writtenToSignalTree = false;
  peakDiagnostics.Branch("raw_entry", &rawEntry);
  peakDiagnostics.Branch("event", &event);
  peakDiagnostics.Branch("detector_id", &detectorId);
  peakDiagnostics.Branch("direction", &XY);
  peakDiagnostics.Branch("raw_chip", &currentRawChip);
  peakDiagnostics.Branch("raw_channel", &currentRawChannel);
  peakDiagnostics.Branch("channel", &channel);
  peakDiagnostics.Branch("global_channel", &currentGlobalChannel);
  peakDiagnostics.Branch("peak_noise_rms", &peakNoiseRms);
  peakDiagnostics.Branch("peak_threshold_adc", &peakThresholdAdc);
  peakDiagnostics.Branch("n_significant_peaks", &pulsePeakCount);
  peakDiagnostics.Branch("peak_selection_passed", &peakAccepted);
  peakDiagnostics.Branch("peak_failure_reason", &pulsePeakFailureReason);
  peakDiagnostics.Branch("written_to_signal_tree", &writtenToSignalTree);
  cout << " --> start to loop" << endl;
  for (int i = 0; i < nentries; i++)
  {
    DrawProcessbar(i, nentries);
    InitData();
    wave->clear();
    frawTree->GetEntry(i);
    //if (board == 5) continue;
    waveAna.clear();
    if (wave->size() == 30) 
    { 
      waveAna = *wave; 
    } 
    else 
    {
      cout << "Invalid waveform data at entry " << i << endl;
      continue;
    }
    ped = 0;
    pedrms = 0;
    currentRawChip = chip;
    currentRawChannel = channel;
    currentIsBM = (chip == 12 || chip == 13);
    int mappedChannel = -1;
    int padPosition = -1;
    currentGlobalChannel =
        GlobalChannelForRaw(board, chip, channel, mappedChannel, padPosition);
    currentPadPosition = padPosition;
    if (currentGlobalChannel < 1)
    {
      op << " [WARNING!] No global noise mapping: board " << board
         << " chip " << chip << " raw channel " << channel << endl;
      continue;
    }
    peakNoiseRms = globalRms[currentGlobalChannel - 1];
    peakThresholdAdc = selection.peakThresholdSigma * peakNoiseRms;

    // Tracker records use the existing detector/tree mapping. BM records
    // use the pad position and are written to BMPadSignals below.
    XY = currentIsBM ? -1 : IndexXY(board, chip);
    Det = currentIsBM ? -1 : IndexDet(board, chip);
    if (!currentIsBM &&
        (Det < 0 || Det >= MMdetN || XY < 0 || XY >= 2))
    {
      op << " [WARNING!] Ignore unmapped board: " << board
         << " chip " << chip << endl;
      continue;
    }
    if (!currentIsBM && mappedChannel < 0)
    {
      op << " [WARNING!] Invalid APV mapping: board " << board
         << " chip " << chip << " raw channel " << channel << endl;
      continue;
    }
    if (currentIsBM)
      channel = padPosition;
    else
      channel = mappedChannel;
    if (!currentIsBM)
    {
      if (i == 0)
      {
        RecordStart = event;
        RecordEnd = event;
        allvalidevent.push_back(event);
      }
      if (RecordEnd < event /*&& event - RecordEnd < 50*/)
      {
        RecordEnd = event;
        NumofTriggers++;
        allvalidevent.push_back(event);
      }
    }
    ped = 0.0;
    pedrms = 0.0;
    CalWaveInfo();
    vector<double> waveform;
    for (int i = 0; i < (int)wave->size(); i++)
    {
      waveform.push_back((double)wave->at(i) - ped);
    }
    const bool singlePeakAccepted = SelectPulseWaveform();
    peakAccepted = !selection.requireSinglePeak || singlePeakAccepted;
    writtenToSignalTree = false;
    rawEntry = i;
    detectorId = currentIsBM ? 5 + (currentRawChip - 12) : Det + 1;
    const int groups[4] = {0, currentIsBM ? 2 : 1,
        currentIsBM ? 6 + currentRawChip - 12 : 3 + Det,
        currentIsBM ? -1 : 8 + Det * 2 + XY};
    for (int g : groups)
    {
      if (g < 0) continue;
      peakAll[g]->Fill(pulsePeakCount);
    }
    riseTandTOT(waveform);
    const bool basicCutPassed = Cut();
    if (!basicCutPassed) continue;
    for (int g : groups)
      if (g >= 0) peakAfterCut[g]->Fill(pulsePeakCount);

    if (currentIsBM)
    {
      // BM pad waveforms go through the same amplitude/timing cut as the
      // Tracker waveforms before they are written to BMPadSignals.
      if (!basicCutPassed) continue;
      ++pulsePeakCheckedCount;
      if (!peakAccepted)
      {
        int reason = pulsePeakFailureReason;
        if (reason <= kPulsePeakOk || reason >= kPulsePeakFailureReasonCount)
          reason = kPulsePeakUnknown;
        ++pulsePeakFailureCounts[reason];
        peakDiagnostics.Fill();
        continue;
      }
      ++pulsePeakPassCount;
      // Keep selected BM waveforms even if a timing fit fails. The fit status
      // and the corresponding -999 timing value are saved in step2.
      // Pause a/e disabled; implementation retained for future comparison.
      // FitPulseArrivalTime();
      FitSigmoidArrivalTime();
      bmTree->Fill();
      writtenToSignalTree = true;
      peakDiagnostics.Fill();
      continue;
    }

    if (!basicCutPassed) continue;
    pulsePeakCheckedCount++;
    if (!peakAccepted)
    {
      int reason = pulsePeakFailureReason;
      if (reason <= kPulsePeakOk || reason >= kPulsePeakFailureReasonCount)
        reason = kPulsePeakUnknown;
      pulsePeakFailureCounts[reason]++;
      peakDiagnostics.Fill();
      continue;
    }
    pulsePeakPassCount++;

    // Previous Pause prerequisite retained, but disabled:
    // FitPulseArrivalTime();
    // if (pulseFitStatus != 0)
    // {
    //   peakDiagnostics.Fill();
    //   continue;
    // }
    // Use the active timing fit as the Tracker acceptance prerequisite.
    FitSigmoidArrivalTime();
    if (selection.requireValidSigmoid &&
        (sigmoidFitStatus != 0 || !std::isfinite(arrivalT_sigmoid) ||
         arrivalT_sigmoid <= -900.0))
    {
      peakDiagnostics.Fill();
      continue;
    }

    // The sigmoid fit above is independent of the disabled Pause fit.

    //if(event < 10) cout << event << "," << timestamp << endl;
    fanaTree[Det * 2 + XY]->Fill();
    writtenToSignalTree = true;
    peakDiagnostics.Fill();
    hChannel->Fill(channel);
    if (eventvec[Det * 2 + XY].size() == 0 || eventvec[Det * 2 + XY].back() != event)
    {
      eventvec[Det * 2 + XY].push_back(event);
      entryvec[Det * 2 + XY].push_back(eventvec[Det * 2 + XY].size());
    }
    hHit[Det * 2 + XY]->Fill(channel);
    if (valideventvec[Det * 2 + XY].size() == 0 || valideventvec[Det * 2 + XY].back() != event)
      valideventvec[Det * 2 + XY].push_back(event);
  }
  cout << " --> pulse peak selection summary" << endl;
  op << " --> pulse peak selection summary" << endl;
  cout << "    checked_after_cut: " << pulsePeakCheckedCount << endl;
  op << "    checked_after_cut: " << pulsePeakCheckedCount << endl;
  cout << "    pass: " << pulsePeakPassCount << endl;
  op << "    pass: " << pulsePeakPassCount << endl;
  const double pulsePeakDenom = pulsePeakCheckedCount > 0 ? (double)pulsePeakCheckedCount : 1.0;
  for (int r = kPulsePeakBadSize; r < kPulsePeakFailureReasonCount; ++r)
  {
    const double frac = 100.0 * (double)pulsePeakFailureCounts[r] / pulsePeakDenom;
    cout << "    fail_" << PulsePeakFailureReasonName(r) << ": " << pulsePeakFailureCounts[r]
         << " (" << frac << "%)" << endl;
    op << "    fail_" << PulsePeakFailureReasonName(r) << ": " << pulsePeakFailureCounts[r]
       << " (" << frac << "%)" << endl;
  }

  fanaFile->cd();
  peakDiagnostics.Write();
  TNamed("TimingFitPolicy",
         "sigmoid_only: sigmoid rising-edge t20; Pause a/e calls disabled; "
         "Pause time branches=-999 and pulseFitStatus=-1; "
         "Tracker sigmoid acceptance is controlled by config; BM failed fits "
         "retain their waveform with invalid timing.").Write();
  TNamed("SignalInputPolicy",
         "decoded ROOT only: pedestal=0; input already baseline-subtracted and "
         "3-sigma channel-selected; run-local noise RMS is used only for peak finding").Write();
  TNamed("PeakSelectionNoiseCalibration",
         Form("source=%s; threshold=%g*noise_sigma; pedestal_mean_ignored=1",
              noiseParameterName.Data(), selection.peakThresholdSigma)).Write();
  TNamed("WaveformSelection",
         Form("amplitude>%g; %g<maxT<%g; %g<riseT<%g; %g<TOT<%g; "
              "peakThreshold=%g*runNoiseRMS; singlePeak=%d; validTrackerSigmoid=%d",
              selection.amplitudeMin, selection.maxTMin, selection.maxTMax,
              selection.riseTimeMin, selection.riseTimeMax,
              selection.totMin, selection.totMax,
              selection.peakThresholdSigma,
              selection.requireSinglePeak, selection.requireValidSigmoid)).Write();
  for (int g = 0; g < nPeakGroups; ++g)
  {
    peakAll[g]->Write();
    peakAfterCut[g]->Write();
  }
  for (int i = 0; i < 2 * MMdetN; i++)
    fanaFile->WriteTObject(fanaTree[i]);
  hChannel->Write();
  if (bmTree)
    fanaFile->WriteTObject(bmTree);
  //// check lost ID
  cout << " --> check lost event ID " << endl;
  op << " --> check lost event ID " << endl;
  TCanvas *cEvent = new TCanvas("Event", "Event");
  cEvent->Divide(2 * MMdetN, 2);
  for (int i = 0; i < 2 * MMdetN; i++)
  {
    TString title = Form("Track%s%dEvent", XYname(i % 2).Data(), i / 2);
    if (eventvec[i].empty())
      gEvent[i] = new TGraph();
    else
      gEvent[i] = new TGraph(eventvec[i].size(), &entryvec[i][0], &eventvec[i][0]);
    gEvent[i]->SetTitle(title);
    gEvent[i]->SetName(title);
    gEvent[i]->GetXaxis()->SetTitle("entryID");
    gEvent[i]->GetYaxis()->SetTitle("eventID");
    cEvent->cd(i + 1);
    gEvent[i]->Draw();
    int lostcounter = eventvec[i].empty() ? 0 : CalculateLostID(eventvec[i]);
    double lostratio = eventvec[i].empty() ? 0 : (double)lostcounter / eventvec[i].size();
    cout << " --> The total numbers of " << title << ": " << eventvec[i].size() << endl
         << " --> The lost numbers of " << title << ": " << lostcounter
         << ", lostratio: " << lostratio << endl;
    op << " --> The total numbers of " << title << ": " << eventvec[i].size() << endl
       << " --> The lost numbers of " << title << ": " << lostcounter
       << ", lostratio: " << lostratio << endl;
  }
  cEvent->Write();

  //// check eff
  // cout << " --> check MM efficiency " << endl;
  // op << " --> check MM efficiency " << endl;
  // TH1I *hMMeff[3];
  // TH1I *hTrackereff[4];
  // for (int i = 0; i < 3; i++)
  // {
  //   hMMeff[i] = new TH1I(Form("MM%seff", XYname(i).Data()), Form("MM%seff", XYname(i).Data()), MMdetN, 0, MMdetN);
  //   hMMeff[i]->GetXaxis()->SetTitle("ID of MM");
  //   hMMeff[i]->GetYaxis()->SetTitle("Counts");
  //   hTrackereff[i] = new TH1I(Form("Tracker%seff", XYname(i).Data()), Form("MM%seff", XYname(i).Data()), MMdetN + 1, 0, MMdetN + 1);
  //   hTrackereff[i]->GetXaxis()->SetTitle("Number of fired MMs");
  //   hTrackereff[i]->GetYaxis()->SetTitle("Counts");
  // }
  // hTrackereff[3] = new TH1I("TrackerXorYeff", "TrackerXorYeff", 2 * MMdetN + 1, 0, 2 * MMdetN + 1);
  // hTrackereff[3]->GetXaxis()->SetTitle("Number of fired Chips");
  // hTrackereff[3]->GetYaxis()->SetTitle("Counts");
  // int matchedycounter = 0;
  // int matchedxcounter = 0;
  // int matchedandcounter = 0;
  // int matchedtotalcounter = 0;
  // vector<int> matchedyvec;
  // vector<int> matchedxvec;
  // vector<int> matchedandvec;
  // cout << " --> start to choose entries " << endl;
  // for (int s = RecordStart; s < RecordEnd; s++)
  // {
  //   DrawProcessbar(s, RecordEnd - RecordStart);
  //   bool valid = IsMatched(s, allvalidevent);
  //   if(valid)
  //   {
  //     matchedxvec.clear();
  //     matchedyvec.clear();
  //     matchedandvec.clear();
  //     bool matchedx = 0;
  //     bool matchedy = 0;
  //     for (int i = 0; i < MMdetN; i++)
  //     {
  //       matchedx = IsMatched(s, valideventvec[i * 2]);
  //       matchedy = IsMatched(s, valideventvec[i * 2 + 1]);
  //       if (matchedx)
  //       {
  //         matchedxvec.push_back(i);
  //         hMMeff[0]->Fill(i);
  //       }
  //       if (matchedy)
  //       {
  //         matchedyvec.push_back(i);
  //         hMMeff[1]->Fill(i);
  //       }
  //       if (matchedy && matchedx)
  //       {
  //       matchedandvec.push_back(i);
  //       hMMeff[2]->Fill(i);
  //       }
  //     }
  //     hTrackereff[0]->Fill(matchedxvec.size());
  //     if (matchedxvec.size())
  //     {
  //       matchedxcounter++;
  //       matchedtotalcounter++;
  //     }

  //     hTrackereff[1]->Fill(matchedyvec.size());
  //     if (matchedyvec.size())
  //     {
  //       matchedycounter++;
  //       matchedtotalcounter++;
  //     }

  //     hTrackereff[2]->Fill(matchedandvec.size());
  //     if (matchedandvec.size())
  //     {
  //       matchedandcounter++;
  //     }
  //     hTrackereff[3]->Fill(matchedxvec.size() + matchedyvec.size());
  //   }
  // }
  // int ValidTracker[4];
  // for (int k = 0; k < 4; k++)
  // {
  //   ValidTracker[k] = hTrackereff[k]->Integral(2, 100);
  // }
  // cout << " --> Trigger Start = " << RecordStart << ", Trigger End =" << RecordEnd << endl
  //      << " --> Trigger coincidence: " << NumofTriggers << endl
  //      << " --> Valid tracker x counter: " << ValidTracker[0] << ", GE: " << (double)ValidTracker[0] / NumofTriggers << endl
  //      << " --> Valid tracker y counter: " << ValidTracker[1] << ", GE: " << (double)ValidTracker[1] / NumofTriggers << endl
  //      << " --> Valid tracker (x&y) counter: " << ValidTracker[2] << ", GE: " << (double)ValidTracker[2] / NumofTriggers << endl
  //      << " --> Valid tracker counter (or): " << ValidTracker[3] << ", GE: " << (double)ValidTracker[3] / NumofTriggers << endl;
  // op << " --> Trigger Start = " << RecordStart << ", Trigger End =" << RecordEnd << endl
  //    << " --> Trigger coincidence: " << NumofTriggers << endl
  //    << " --> Valid tracker x counter: " << ValidTracker[0] << ", GE: " << (double)ValidTracker[0] / NumofTriggers << endl
  //    << " --> Valid tracker y counter: " << ValidTracker[1] << ", GE: " << (double)ValidTracker[1] / NumofTriggers << endl
  //    << " --> Valid tracker (x&y) counter: " << ValidTracker[2] << ", GE: " << (double)ValidTracker[2] / NumofTriggers << endl
  //    << " --> Valid tracker counter (or): " << ValidTracker[3] << ", GE: " << (double)ValidTracker[3] / NumofTriggers << endl;

  // TCanvas *cHit = new TCanvas("Hit", "Hit");
  // cHit->Divide(MMdetN, 2);
  // for (int i = 0; i < MMdetN; i++)
  // {
  //   cHit->cd(i + 1);
  //   hHit[2 * i]->Draw();
  //   cHit->cd(i + MMdetN + 1);
  //   hHit[2 * i + 1]->Draw();
  // }
  // cHit->Write();

  // TCanvas *cEff = new TCanvas("MMeff", "MMeff");
  // cEff->Divide(4, 2);
  // for (int s = 0; s < 3; s++)
  // {
  //   cEff->cd(s + 1);
  //   hMMeff[s]->Draw();
  //   double MMeff[MMdetN];
  //   for (int i = 0; i < MMdetN; i++)
  //   {
  //     MMeff[MMdetN] = hMMeff[s]->GetBinContent(i + 1) / NumofTriggers;
  //     TString latexbuff = Form("MMID=%d,Ratio=%.1f%%", i, MMeff[i] * 100);
  //     TLatex *la = CreatLatex(latexbuff, 0.2, 0.7 - 0.07 * i, 42);
  //     la->Draw("same");
  //   }
  // }

  // for (int s = 0; s < 4; s++)
  // {
  //   cEff->cd(s + 5);
  //   hTrackereff[s]->Draw();
  //   int total = 0;
  //   if (s < 3)
  //     total = MMdetN + 1;
  //   else
  //     total = 2 * MMdetN + 1;
  //   double Trackereff[2 * MMdetN + 1];
  //   for (int i = 0; i < total; i++)
  //   {
  //     Trackereff[i] = hTrackereff[s]->GetBinContent(i + 1) / NumofTriggers;
  //     TString latexbuff;
  //     if (s == 3)
  //       latexbuff = Form("Nchips=%d,Ratio=%.1f%%", i, Trackereff[i] * 100);
  //     else
  //       latexbuff = Form("NMM=%d,Ratio=%.1f%%", i, Trackereff[i] * 100);
  //     TLatex *la = CreatLatex(latexbuff, 0.2, 0.7 - 0.05 * i, 42);
  //     la->Draw("same");
  //   }
  // }
  // cEff->Write();

  fanaFile->Close();
  frawFile->Close();
  cout << " --> Analyzed results has been stored to " << anaName << endl;
  op << " --> Analyzed results has been stored to " << anaName << endl;
  op.close();
}

// Legacy Pause implementation retained intact; calls are commented out.
double AnalysisAPV::FitPulseArrivalTime()
{
  arrivalT = -999;
  arrivalT_a = -999;
  arrivalT_e = -999;
  arrivalT_sigmoid = -999;
  sigmoidT50 = -999;
  sigmoidBaseline = -999;
  sigmoidHeight = -999;
  sigmoidK = -999;
  sigmoidFitChi2NDF = -1;
  sigmoidFitStatus = -1;
  sigmoidFitStart = -1;
  sigmoidFitEnd = -1;
  pulseFitStatus = -1;
  pulseFitChi2NDF = -1;

  if (!wave || wave->size() < 4 || maxA <= 0 || maxT <= 0)
    return arrivalT;

  const double safePedRms = std::max(1.0, std::abs(pedrms));
  const double threshold = std::max(3.0 * safePedRms, 0.1 * maxA);
  int fitStart = 0;
  for (int i = maxT; i > 0; --i)
  {
    const double adc = (double)wave->at(i) - ped;
    const double prevAdc = (double)wave->at(i - 1) - ped;
    if (adc >= threshold && prevAdc < threshold)
    {
      fitStart = std::max(0, i - 1);
      break;
    }
  }

  int fitEnd = std::min(maxT, (int)wave->size() - 1);
  if (fitEnd - fitStart + 1 < 4)
  {
    fitStart = std::max(0, fitEnd - 3);
  }
  if (fitEnd - fitStart + 1 < 4)
    return arrivalT;

  TGraph gWave;
  int nPoints = 0;
  for (int i = fitStart; i <= fitEnd; i++)
  {
    gWave.SetPoint(nPoints++, i, (double)wave->at(i));
  }

  const double initialB = ped;
  const double initialShapeB = std::max(0.5, ((double)maxT - fitStart) / 3.0);
  const double initialArrival = (double)maxT - 3.0 * initialShapeB;
  const double initialShapeA = std::max(1.0, maxA * std::exp(3.0) / 27.0);

  TF1 fitFunc("pulseFitModel", PulseModel, fitStart, fitEnd, 4);
  fitFunc.SetParNames("B", "A", "a", "b");
  fitFunc.SetParameters(initialB, initialShapeA, initialArrival, initialShapeB);
  // Keep the Pause-fit baseline fixed to the global pedestal measured for
  // this ADC channel.  A signal waveform must not refit its own pedestal.
  fitFunc.FixParameter(0, ped);
  fitFunc.SetParLimits(1, 0.0, std::max(100.0, 10.0 * maxA));
  fitFunc.SetParLimits(2, fitStart - 10.0, (double)fitEnd);
  fitFunc.SetParLimits(3, 0.1, 20.0);

  TFitResultPtr fitResult = gWave.Fit(&fitFunc, "QSN");
  pulseFitStatus = fitResult;
  if (pulseFitStatus == 0)
  {
    const double a = fitFunc.GetParameter(2);
    const double b = fitFunc.GetParameter(3);
    const double target = 0.20 * std::pow(3.0, 3) * std::exp(-3.0);
    double lo = 0.0;
    double hi = 3.0;
    for (int iteration = 0; iteration < 100; ++iteration)
    {
      const double u = 0.5 * (lo + hi);
      const double shape = std::pow(u, 3) * std::exp(-u);
      if (shape < target)
        lo = u;
      else
        hi = u;
    }
    const double u20 = 0.5 * (lo + hi);
    arrivalT_a = a;
    arrivalT_e = a + b * u20;
    arrivalT = arrivalT_a;
    const int ndf = fitFunc.GetNDF();
    if (ndf > 0)
      pulseFitChi2NDF = fitFunc.GetChisquare() / ndf;
  }

  return arrivalT;
}

double AnalysisAPV::FitSigmoidArrivalTime()
{
  arrivalT_sigmoid = -999;
  sigmoidT50 = -999;
  sigmoidBaseline = -999;
  sigmoidHeight = -999;
  sigmoidK = -999;
  sigmoidFitChi2NDF = -1;
  sigmoidFitStatus = -1;
  sigmoidFitStart = -1;
  sigmoidFitEnd = -1;

  if (!wave || wave->size() < 5 || maxA <= 0 || maxT <= 0)
    return arrivalT_sigmoid;

  const double safePedRms = std::max(1.0, std::abs(pedrms));
  const double threshold = std::max(3.0 * safePedRms, 0.1 * maxA);
  int thresholdStart = 0;
  for (int i = maxT; i > 0; --i)
  {
    const double adc = (double)wave->at(i) - ped;
    const double prevAdc = (double)wave->at(i - 1) - ped;
    if (adc >= threshold && prevAdc < threshold)
    {
      thresholdStart = std::max(0, i - 1);
      break;
    }
  }

  sigmoidFitStart = std::max(0, thresholdStart - 3);
  sigmoidFitEnd = std::min(maxT, (int)wave->size() - 1);
  if (sigmoidFitEnd - sigmoidFitStart + 1 < 5)
    sigmoidFitStart = std::max(0, sigmoidFitEnd - 4);
  if (sigmoidFitEnd - sigmoidFitStart + 1 < 5)
    return arrivalT_sigmoid;

  double t50Initial = (double)maxT;
  for (int i = maxT; i > 0; --i)
  {
    const double adc = (double)wave->at(i) - ped;
    const double prevAdc = (double)wave->at(i - 1) - ped;
    if (adc >= 0.5 * maxA && prevAdc < 0.5 * maxA)
    {
      const double denom = adc - prevAdc;
      const double fraction = denom != 0.0 ? (0.5 * maxA - prevAdc) / denom : 0.0;
      t50Initial = (i - 1) + fraction;
      break;
    }
  }

  TGraph gWave;
  int nPoints = 0;
  for (int i = sigmoidFitStart; i <= sigmoidFitEnd; ++i)
    gWave.SetPoint(nPoints++, i, (double)wave->at(i));

  TF1 fitFunc("sigmoidRisingFit", SigmoidModel, sigmoidFitStart, sigmoidFitEnd, 4);
  fitFunc.SetParNames("B", "H", "t50", "k");
  fitFunc.SetParameters(ped, maxA, t50Initial,
                        std::max(0.2, ((double)maxT - t50Initial) / 2.0));
  // Use the global channel pedestal as a fixed baseline, consistently with
  // the Pause fit and the waveform selection cuts.
  fitFunc.FixParameter(0, ped);
  fitFunc.SetParLimits(1, 0.0, std::max(100.0, 10.0 * maxA));
  fitFunc.SetParLimits(2, sigmoidFitStart - 5.0, sigmoidFitEnd + 5.0);
  fitFunc.SetParLimits(3, 0.05, 20.0);

  TFitResultPtr fitResult = gWave.Fit(&fitFunc, "QSN");
  sigmoidFitStatus = fitResult;
  if (sigmoidFitStatus != 0)
    return arrivalT_sigmoid;

  sigmoidBaseline = fitFunc.GetParameter(0);
  sigmoidHeight = fitFunc.GetParameter(1);
  sigmoidT50 = fitFunc.GetParameter(2);
  sigmoidK = fitFunc.GetParameter(3);
  arrivalT_sigmoid = SigmoidFractionTime(sigmoidT50, sigmoidK, 0.20);
  const int ndf = fitFunc.GetNDF();
  if (ndf > 0)
    sigmoidFitChi2NDF = fitFunc.GetChisquare() / ndf;

  if (!std::isfinite(arrivalT_sigmoid))
  {
    arrivalT_sigmoid = -999;
    return arrivalT_sigmoid;
  }
  return arrivalT_sigmoid;
}
bool AnalysisAPV::SelectPulseWaveform()
{
  pulsePeakFailureReason = kPulsePeakOk;
  pulsePeakCount = -1;


  if (!wave || wave->size() != kExpectedWaveSamples)
  {
    pulsePeakFailureReason = kPulsePeakBadSize;
    return false;
  }

  bool hasSampleOverThreshold = false;
  vector<int> valid_peaks;

  for (int i = 1; i < kExpectedWaveSamples - 1; ++i)
  {
    if (((double)wave->at(i) - ped) <= peakThresholdAdc)
      continue;
    hasSampleOverThreshold = true;

    if (wave->at(i) == wave->at(i + 1))
    {
      int plateau_start = i;
      int plateau_end = i;
      for (int j = i + 1; j < kExpectedWaveSamples; ++j)
      {
        if (wave->at(j) == wave->at(i))
          plateau_end = j;
        else
          break;
      }

      if (plateau_start > 0 && plateau_end < kExpectedWaveSamples - 1)
      {
        bool left_ok = true;
        bool right_ok = true;
        for (int j = plateau_start - 1; j >= std::max(0, plateau_start - 3); --j)
        {
          if (wave->at(j) >= wave->at(i))
          {
            left_ok = false;
            break;
          }
        }
        for (int j = plateau_end + 1; j <= std::min(kExpectedWaveSamples - 1, plateau_end + 3); ++j)
        {
          if (wave->at(j) >= wave->at(i))
          {
            right_ok = false;
            break;
          }
        }
        if (left_ok && right_ok)
        {
          valid_peaks.push_back((plateau_start + plateau_end) / 2);
          i = plateau_end;
          continue;
        }
      }
    }

    bool is_local_max = true;
    const int left_bound = std::max(0, i - 3);
    const int right_bound = std::min(kExpectedWaveSamples - 1, i + 3);
    for (int j = left_bound; j <= right_bound; ++j)
    {
      if (j == i)
        continue;
      if (wave->at(j) >= wave->at(i))
      {
        is_local_max = false;
        break;
      }
    }
    if (is_local_max)
      valid_peaks.push_back(i);
  }

  pulsePeakCount = static_cast<int>(valid_peaks.size());
  if (valid_peaks.empty())
  {
    pulsePeakFailureReason = hasSampleOverThreshold ? kPulsePeakNoLocalPeak : kPulsePeakNoSampleOverThreshold;
    return false;
  }

  if (valid_peaks.size() != 1)
  {
    pulsePeakFailureReason = kPulsePeakMultiPeak;
    return false;
  }

  const int peak_pos = valid_peaks[0];
  if (peak_pos <= 0 || peak_pos >= kExpectedWaveSamples - 1)
  {
    pulsePeakFailureReason = kPulsePeakEdgePeak;
    return false;
  }
  if (((double)wave->at(peak_pos) - ped) <= peakThresholdAdc)
  {
    pulsePeakFailureReason = kPulsePeakLowAmplitude;
    return false;
  }
  if (wave->at(0) >= wave->at(peak_pos) || wave->at(kExpectedWaveSamples - 1) >= wave->at(peak_pos))
  {
    pulsePeakFailureReason = kPulsePeakBadEdge;
    return false;
  }

  return true;
}
void AnalysisAPV::InitData()
{
  maxA = 0;
  maxT = 0;
  ped = 0;
  pedrms = 0;
  peakNoiseRms = 0;
  peakThresholdAdc = 0;
  nisel = 0;
  riseT = 0;
  TOT = 0;
  arrivalT = -999;
  arrivalT_a = -999;
  arrivalT_e = -999;
  arrivalT_sigmoid = -999;
  sigmoidT50 = -999;
  sigmoidBaseline = -999;
  sigmoidHeight = -999;
  sigmoidK = -999;
  sigmoidFitChi2NDF = -1;
  sigmoidFitStatus = -1;
  sigmoidFitStart = -1;
  sigmoidFitEnd = -1;
  pulseFitStatus = -1;
  pulseFitChi2NDF = -1;
  pulsePeakFailureReason = kPulsePeakOk;
  currentGlobalChannel = -1;
  pulsePeakCount = -1;
  currentIsBM = false;
  currentRawChip = -1;
  currentRawChannel = -1;
  currentPadPosition = -1;
  //timestamp = 0;
}

int AnalysisAPV::IndexXY(int board, int chip)
{
  return chip / 2;
}

int AnalysisAPV::IndexDet(int board, int chip)
{
  return board;
}

int AnalysisAPV::IndexStrip(int board, int chip, int channel)
{
  if (board < 0 || board >= MMdetN ||
      chip < 0 || chip >= apvChipsPerTracker ||
      channel < 0 || channel >= apvChannelsPerChip)
    return -1;

  // map_APV.txt is ordered by board, chip, channel. The board offset
  // is required because every Tracker uses local chip IDs 0-3.
  const int index = (board * apvChipsPerTracker + chip) * apvChannelsPerChip + channel;
  if (index < 0 || index >= static_cast<int>(APVmap.size()))
    return -1;
  return APVmap.at(index);
}

int AnalysisAPV::CalculateLostID(vector<int> eventvec)
{
  int start = 1;
  int lostcounter = 0;
  while (start < eventvec.size())
  {
    int l = abs(eventvec.at(start) - eventvec.at(start - 1));
    if (l != 1 && l < 50) // if l>50, may be stop record
    {
      for (int i = 1; i < l; i++)
      {
        lostcounter++;
        //cout << " --> lost event " << eventvec.at(start - 1) + i << " between (" << eventvec.at(start - 1) << "," << eventvec.at(start) << ") !" << endl;
      }
    }
    start++;
  }
  return lostcounter;
}

void AnalysisAPV::riseTandTOT(vector<double> waveform)
{
  if (maxA < 0)
  {
    InitData();
    return;
  }

  double cfamp[4] = {0.2 * maxA, 0.8 * maxA, 0.5 * maxA, 0.5 * maxA};
  double stamp[4] = {double(maxT) - 1, double(maxT), double(maxT), double(maxT) + 1};

  // find start 0.2*maxA
  for (int i = maxT; i > 0; i--)
    if (waveform[i] - cfamp[0] > 0 && waveform[i - 1] - cfamp[0] < 0)
    {
      stamp[0] = i;
      break;
    }
  // find start 0.8*maxA
  for (int i = maxT; i > 0; i--)
    if (waveform[i] - cfamp[1] > 0 && waveform[i - 1] - cfamp[1] < 0)
    {
      stamp[1] = i;
      break;
    }
  // find start 0.5*maxA
  for (int i = maxT; i > 0; i--)
    if (waveform[i] - cfamp[2] > 0 && waveform[i - 1] - cfamp[2] < 0)
    {
      stamp[2] = i;
      break;
    }
  // find end 0.5*maxA
  for (int i = maxT; i < waveform.size() - 1; i++)
    if (waveform[i] - cfamp[3] > 0 && waveform[i + 1] - cfamp[3] < 0)
    {
      stamp[3] = i + 1;
      break;
    }
  double time[4];
  for (int i = 0; i < 4; i++)
  {
    double y1 = stamp[i] - 1;
    double y2 = stamp[i];
    double x1 = waveform[stamp[i] - 1];
    double x2 = waveform[stamp[i]];
    time[i] = (y2 - y1) * (cfamp[i] - x1) / (x2 - x1) + y1;
  }
  riseT = time[1] - time[0];
  TOT = time[3] - time[2];
}

void AnalysisAPV::Cal_ped()
{
  if (LoadNoiseParameters())
    return;
  cout << " --> No usable external noise calibration; estimating Tracker pedestals from raw waveforms" << endl;
  vector<short> channelData[MMdetN * 2][256];
  ifstream fitfile("config/fit_parameters.txt");
  if (fitfile.is_open()) 
  {
    string line;
    int global_ch;
    double mean, sigma;
    int count = 0;
    while (getline(fitfile, line)) 
    {
      // Skip comment lines.
      if (line.empty() || line[0] == '#') continue;

      istringstream iss(line);
      if (iss >> global_ch >> mean >> sigma) 
      {
        // //璁＄畻iTree鍜宑h (鍏ㄥ眬閫氶亾鍙蜂粠1寮€濮?
        // if (global_ch < 1025)
        // {
        //   continue;
        // }
        // if (global_ch > 1024)
        // {
        //   global_ch = global_ch - 1024;
        // }

        // if (global_ch > 1024)
        // {
        //   if (global_ch < 2049)
        //   {
        //     continue;
        //   }
        //   global_ch = global_ch - 1024;
        // }

        int index = global_ch - 1;
        int iTree = index / 256;      // 姣忎釜iTree鏈?56涓€氶亾
        int ch = index % 256;          // 閫氶亾鍦╥Tree鍐呯殑绱㈠紩

        // 纭繚绱㈠紩鍦ㄦ湁鏁堣寖鍥村唴
        
        if (iTree < MMdetN * 2 && ch < 256) 
        {
          channelped[iTree][ch] = mean;
          channelrms[iTree][ch] = sigma;
          count++;
        }
      }
    }
    fitfile.close();
    cout << count << endl;
    cout << MMdetN << endl;
    // 妫€鏌ユ槸鍚︽垚鍔熻鍙栨墍鏈夐€氶亾
    //count = 0;
    if (count == 2 * 256 * MMdetN) 
    {  
      cout << "Successfully read pedestal parameters from file." << endl;
      channelCalibrationLoaded = true;
      return;
    }
    else
    {
      cout << "read pedestal parameters fail." << endl;
    }
  }

  cout << "not parameters file." << endl;
  int nentries = frawTree->GetEntriesFast();
  for (int iEntry = 0 ; iEntry < nentries ; iEntry++)
  {
    frawTree->GetEntry(iEntry);
    const int ch = IndexStrip(board, chip, channel);
    if (ch < 0)
      continue;
    int treeIndex = board * 2 + (chip / 2);
    channelData[treeIndex][ch].push_back(wave->at(0));
  }
  for (int iTree = 0; iTree < 2 * MMdetN; iTree++)
  {
    double fit_mean;
    double fit_sigma;
    for (int ch = 0; ch < 256; ch++)
    {
      TString histName = Form("hWave0_Ch%d", ch);
      TH1F* hWave0 = new TH1F(histName, Form("Channel %d Wave0;ADC;Counts", ch),4096, -2048, 2047);
      
      for (auto val : channelData[iTree][ch]) 
      {
        hWave0->Fill(val);
      }
      int maxBin = hWave0->GetMaximumBin();
      double peakCenter = hWave0->GetBinCenter(maxBin);
      double fitLow = peakCenter - 200;
      double fitHigh = peakCenter + 200;

      fitLow = std::max(fitLow, hWave0->GetXaxis()->GetXmin());
      fitHigh = std::min(fitHigh, hWave0->GetXaxis()->GetXmax());
      double histMean = hWave0->GetMean();
      double histRMS = hWave0->GetRMS();
      TF1 *gaussFit = new TF1("gaussFit", "gaus", fitLow, fitHigh);
      gaussFit->SetParameters(hWave0->GetMaximum(), histMean, histRMS);
      hWave0->Fit(gaussFit, "QR", "", fitLow, fitHigh);
      TF1 *fitFunc = hWave0->GetFunction("gaussFit");
      if (fitFunc) 
      {
        fit_mean = fitFunc->GetParameter(1);
        fit_sigma = fitFunc->GetParameter(2);
      } 

      channelped[iTree][ch] = fit_mean;
      channelrms[iTree][ch] = fit_sigma;
      //cout << iTree << "channel" << ch << channelped[iTree][ch] ;
      delete hWave0;
      delete gaussFit;// 閲婃斁鎷熷悎鍑芥暟
    }
  }
  for (int treeIndex = 0; treeIndex < 2 * MMdetN; ++treeIndex)
    for (int localChannel = 0; localChannel < 256; ++localChannel)
    {
      const int globalIndex = treeIndex * 256 + localChannel;
      globalPedestal[globalIndex] = channelped[treeIndex][localChannel];
      globalRms[globalIndex] = std::max(channelrms[treeIndex][localChannel], kMinLocalPedestalRms);
    }
  channelCalibrationLoaded = true;
}


void AnalysisAPV::CalWavePed()
{
  ped = 0.0;
  pedrms = 1.0e30;

  if (!channelCalibrationLoaded ||
      (currentIsBM && !bmNoiseCalibrationLoaded) ||
      currentGlobalChannel < 1 ||
      currentGlobalChannel > kGlobalNoiseChannels)
    return;

  const int index = currentGlobalChannel - 1;
  ped = globalPedestal[index];
  pedrms = std::max(globalRms[index], kMinLocalPedestalRms);
}
void AnalysisAPV::CalWaveInfo()
{
  int peakstart = 0;
  int peakend = (int)wave->size();
  int threshold = -384;
  vector<int> id;
  maxA = -65536;
  for (int i = peakstart; i < peakend; i++)
  {
    if (maxA < (double)wave->at(i))
    {
      maxA = (double)wave->at(i);
      maxT = i;
    }
  }
  maxA -= ped;
  // Number of waveform samples above the 7-sigma amplitude threshold used
  // by Cut(). This is saved as the pad step2 nisel quantity.
  nisel = 0;
  const double selectionThreshold =
      rmscut * std::max(std::abs(pedrms), kMinLocalPedestalRms);
  for (int i = 0; i < peakend; ++i)
    if (static_cast<double>(wave->at(i)) - ped > selectionThreshold)
      ++nisel;
  //if(maxT < 15 && maxA > 1200)
  if(0)
  {
    bool c = 1;
    for (int s = 0 ; s < maxT + 5 ; s++)
    {
      if((double)wave->at(maxT) - (double)wave->at(s) < 20)
      {
        id.push_back(s);
      }
    }
    for (int s = 1 ; s < id.size() ; s++)
    {
      int n = id[s] - id[s-1] ;
      if (n != 1) c = 0;
    }
    if(c && id.size() > 2)
    {
      /*int plateauStart = id[0];     
      int plateauEnd = id.back();   
      ouble ct = (plateauEnd - plateauStart)/2;
      maxA = maxA + 150 * ct;*/
    
      int plateauStart = id[0];
      int plateauEnd = id.back();
      double ct = (plateauEnd - plateauStart)/2.0;
      maxA = maxA + 150 * ct;

      // --- 鏂板锛氬鎵炬嫙鍚堟暟鎹偣鑼冨洿 ---
      int startFit = 0;
      int endFit = wave->size() - 1;
      // 浠庡潽鍖哄墠瀵绘壘绗竴涓?100鐨勭偣锛屽彇鍏跺墠涓€涓偣
      for (int i = plateauStart - 1; i >= 0; --i) 
      {
        double adc = (double)wave->at(i) - ped;
        if (adc > 100) 
        {
          startFit = (i > 0) ? (i - 1) : 0;
          break;
        }
      }

      // 浠庡潽鍖哄悗瀵绘壘绗竴涓?100鐨勭偣锛屽彇鍏跺悗涓€涓偣
      for (int i = plateauEnd + 1; i < wave->size(); ++i) 
      {
        double adc = (double)wave->at(i) - ped;
        if (adc < 100) 
        {
          endFit = (i < wave->size() - 1) ? (i + 1) : wave->size() - 1;
          break;
        }
      }



      TGraph* gWave = new TGraph();
      int nPoints = 0;
      for (int i = startFit; i <= endFit; i++) 
      {
        double adc = (double)wave->at(i) - ped;
        gWave->SetPoint(nPoints, i, adc);
        nPoints++;
      }


      auto fitFunc = new TF1
      (
        "fitFunc",
        [](double* x, double* par)
        {
          double t = x[0];
          double I = par[0];
          double RC = par[1];
          double t1 = par[2];
          double t2 = par[3];
          double term = (1 - exp(-t1/RC)) * exp(t2/RC) + 1;
          double t0 = (term > 0) ? RC * log(term) : (t1 + t2)/2.0;
          double S = I * (1 - exp(-t1/RC));
          if (t < t1)
            return I * (1 - exp(-t/RC));
          else if (t <= t2)
            return S;
          else
            return I * exp(-t/RC) * (exp(t0/RC) - 1);
        },
        0, wave->size()-1, 4
      );



      fitFunc->SetParameters(
        maxA,
        100.0,
        plateauStart,
        plateauEnd
      );
      fitFunc->SetParLimits(2, plateauStart-2, plateauStart+2);
      fitFunc->SetParLimits(3, plateauEnd-2, plateauEnd+2);     // t2 鍦ㄥ潽鍖洪檮杩?




      // 4. 鎵ц鎷熷悎
      bool fitSuccess = (gWave->Fit(fitFunc, "SQ") == 0);

      if (fitSuccess) 
      {
        // 鑾峰彇鎷熷悎鍙傛暟
        double I_fit = fitFunc->GetParameter(0);  // I
        double RC_fit = fitFunc->GetParameter(1); // RC
        double t1_fit = fitFunc->GetParameter(2); // t1
        double t2_fit = fitFunc->GetParameter(3); // t2

        // 鍔ㄦ€佽绠?t0 (閫氳繃鏁板鍏紡)
        double term = (1 - exp(-t1_fit/RC_fit)) * exp(t2_fit/RC_fit) + 1;
        if (term <= 0) 
        {
          term = 1.0; // 闃叉瀵规暟璐熸暟
        }
        double t0_fit = RC_fit * log(term);
        double peakValue = I_fit * (1 - exp(-t0_fit/RC_fit));

        // 鏇存柊 maxA
        if (peakValue > maxA - 150 * ct) 
        {
          maxA = peakValue;
          maxT = t0_fit; // 鏇存柊宄板€兼椂闂达紙鍙€夛級
        }
      }  

      /*// 6. 缁樺浘淇濆瓨
      gSystem->mkdir("FitPlots", kTRUE);
      TCanvas* cFit = new TCanvas("cFit", "Fit Result", 800, 600);
    
      // 缁樺埗鍘熷鏁版嵁锛堟墍鏈夌偣锛?      TGraph* gAll = new TGraph();
      for (int i = 0; i < (int)wave->size(); i++) 
      {
        gAll->SetPoint(i, i, (double)wave->at(i) - ped); // 鎵ｉ櫎鍩虹嚎
      }
      gAll->SetTitle(Form("Tracker%d Event%d;Sample;ADC", Det, event));
      gAll->SetMarkerStyle(20);
      gAll->GetYaxis()->SetRangeUser(-200, 2000); // 鏂板浠ｇ爜
      gAll->Draw("AP");
      gAll->Draw("AP");

      // 缁樺埗鎷熷悎鏇茬嚎
      fitFunc->SetLineColor(kRed);
      fitFunc->Draw("same");
      gAll->GetYaxis()->SetRangeUser(-200, 2000); // 鑻ユ湭缁ф壙璁剧疆锛屽彲鏄惧紡璋冪敤

      // 鏍囪鍧尯鑼冨洿
      TBox* plateauBox = new TBox
      (
        plateauStart, gAll->GetYaxis()->GetXmin(),
        plateauEnd, gAll->GetYaxis()->GetXmax()
      );
      plateauBox->SetFillColorAlpha(kGray, 0.3);
      plateauBox->Draw("same");

      cFit->SaveAs(Form("FitPlots/Tracker%d_Event%d_Fit.png", Det, event));

      // 7. 娓呯悊璧勬簮
      delete gWave;
      delete gAll;
      delete fitFunc;
      delete plateauBox;
      delete cFit;*/


      delete gWave;
      delete fitFunc;
    }
  }
  // ped and pedrms remain the global per-ADC-channel calibration values.



}












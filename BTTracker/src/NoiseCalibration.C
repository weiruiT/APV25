#include "CRAnalysis.h"
#include "CRConfig.h"
#include "DecodedAPVRoot2TrackRaw.h"
#include "RunDataInput.h"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFitResult.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

namespace
{
struct NoiseAccumulator
{
  long long count = 0;
  double mean = 0.0;
  double m2 = 0.0;

  void Add(double value)
  {
    ++count;
    const double delta = value - mean;
    mean += delta / static_cast<double>(count);
    const double delta2 = value - mean;
    m2 += delta * delta2;
  }

  double Sigma() const
  {
    return count > 0 ? std::sqrt(m2 / static_cast<double>(count)) : 0.0;
  }
};

struct PadMapRow
{
  int board = -1;
  int chip = -1;
  int channelId = -1;
  int connectorPin = -1;
  int position = -1;
  int globalChannel = -1;
};

struct PadDetail
{
  int readoutBoard = -1;
  int rawChip = -1;
  int rawChannel = -1;
  int connectorPin = -1;
  int position = -1;
  int mapBoard = -1;
  int globalChannel = -1;
};

bool LoadAPVMap(vector<int> &apvMap)
{
  ifstream mapFile("config/map_APV.txt");
  if (!mapFile.is_open())
  {
    cerr << "ERROR: cannot open config/map_APV.txt" << endl;
    return false;
  }

  int board = 0;
  int chip = 0;
  int channel = 0;
  int unused = 0;
  int strip = 0;
  while (mapFile >> board >> chip >> channel >> unused >> strip)
    apvMap.push_back(strip);

  const int expected = MMdetN * apvChipsPerTracker * apvChannelsPerChip;
  if (static_cast<int>(apvMap.size()) < expected)
  {
    cerr << "ERROR: map_APV.txt has " << apvMap.size()
         << " rows, need at least " << expected << endl;
    return false;
  }
  if (static_cast<int>(apvMap.size()) > expected)
    apvMap.resize(expected);
  return true;
}

int MapRawChannel(const vector<int> &apvMap, int board, int chip, int channel)
{
  if (board < 0 || board >= MMdetN ||
      chip < 0 || chip >= apvChipsPerTracker ||
      channel < 0 || channel >= apvChannelsPerChip)
    return -1;

  const int index =
      (board * apvChipsPerTracker + chip) * apvChannelsPerChip + channel;
  if (index < 0 || index >= static_cast<int>(apvMap.size()))
    return -1;
  return apvMap[index];
}

bool LoadPadMap(vector<PadMapRow> &rows, TString &usedName)
{
  ifstream mapFile("config/map_APV4BMM.txt");
  usedName = "config/map_APV4BMM.txt";

  PadMapRow row;
  while (mapFile >> row.board >> row.chip >> row.channelId
                 >> row.connectorPin >> row.position)
  {
    if (row.position >= 0 && row.position < 121)
      rows.push_back(row);
  }
  return !rows.empty();
}

const PadMapRow *UniquePadMatch(const vector<PadMapRow> &rows,
                                int board, int chip, int channelId,
                                bool requireChip)
{
  const PadMapRow *match = nullptr;
  int position = -1;

  for (const PadMapRow &row : rows)
  {
    if (row.board != board || row.channelId != channelId)
      continue;
    if (requireChip && row.chip != chip)
      continue;

    if (!match)
    {
      match = &row;
      position = row.position;
    }
    else if (row.position != position)
    {
      // A local-chip map with more than one chip cannot be resolved without
      // an explicit raw-chip row, so reject this fallback as ambiguous.
      return nullptr;
    }
  }
  return match;
}

const PadMapRow *FindPadMapRow(const vector<PadMapRow> &rows,
                               int readoutBoard, int rawChip, int rawChannel)
{
  // In the raw stream both BM chips can be in readoutBoard 1. The two
  // map blocks therefore follow the raw chip number, not readoutBoard:
  // chip12 -> map board0, chip13 -> map board1.
  int mapBoard = -1;
  if (rawChip == 12)
    mapBoard = 0;
  else if (rawChip == 13)
    mapBoard = 1;
  else if (readoutBoard > 0)
    mapBoard = readoutBoard - 1;
  else
    mapBoard = readoutBoard;

  // map_APV4BMM channel_id is 1-based; TrackRaw channel is 0-based.
  const int mapChannelId = rawChannel + 1;
  if (mapBoard < 0 || mapChannelId <= 0)
    return nullptr;

  for (const PadMapRow &row : rows)
  {
    if (row.board != mapBoard || row.channelId != mapChannelId)
      continue;

    // Current map files use chip=0 as a local placeholder. Also accept
    // explicit raw-chip or local-chip values in future map files.
    if (row.chip == 0 || row.chip == rawChip ||
        row.chip == rawChip - 12)
      return &row;
  }

  return nullptr;
}

int PadKey(const PadDetail &detail)
{
  const int chipSlot = detail.rawChip == 13 ? 1 : 0;
  return ((detail.mapBoard * 2) + chipSlot) * apvChannelsPerChip +
         detail.rawChannel;
}

// Fitted sigma above this limit is considered unphysical for noise
// waveforms. Typical tracker noise is ~7-15 ADC; contaminated channels
// whose fit latches onto a wrong peak produce sigma of 50-360 ADC.
constexpr double kMaxNoiseSigma = 50.0;

bool HistQuantiles(TH1D *h, double &q16, double &q50, double &q84)
{
  const double total = h->Integral();
  if (total <= 0.0)
    return false;

  double acc = 0.0;
  bool got16 = false;
  bool got50 = false;
  bool got84 = false;
  q16 = q50 = q84 = 0.0;
  for (int b = 1; b <= h->GetNbinsX(); ++b)
  {
    acc += h->GetBinContent(b);
    const double frac = acc / total;
    const double x = h->GetBinCenter(b);
    if (!got16 && frac >= 0.16)
    {
      q16 = x;
      got16 = true;
    }
    if (!got50 && frac >= 0.50)
    {
      q50 = x;
      got50 = true;
    }
    if (!got84 && frac >= 0.84)
    {
      q84 = x;
      got84 = true;
    }
    if (got16 && got50 && got84)
      break;
  }
  return got16 && got50 && got84;
}

struct NoiseFitResult
{
  int status = -1;
  TString reason;
  double mean = 0.0;
  double sigma = 0.0;
  double chi2Ndf = 0.0;
  double robustMean = 0.0;
  double robustSigma = 0.0;
};

// Robust Gaussian fit of a noise waveform histogram.
//
// status codes:
//   0 = fit OK
//   1 = no data (empty channel)
//   2 = MINUIT fit failed / invalid result
//   3 = fitted sigma unphysical (<=0, NaN or above kMaxNoiseSigma)
//   4 = fitted mean far from robust center (wrong peak selected)
//   5 = contamination still present after chi2-trim refit
//   6 = degenerate distribution (robust width < 0.5 ADC)
//
// On failure the robust statistics (median and rescaled quantile width)
// are returned as the fallback mean/sigma.
NoiseFitResult FitNoiseGaussian(TH1D *h, long long sampleCount, int gc)
{
  NoiseFitResult out;
  if (sampleCount <= 0)
  {
    out.status = 1;
    out.reason = "no samples";
    return out;
  }

  double q16 = 0.0;
  double q50 = 0.0;
  double q84 = 0.0;
  if (!HistQuantiles(h, q16, q50, q84))
  {
    out.status = 1;
    out.reason = "empty histogram";
    return out;
  }

  out.robustMean = q50;
  // For a Gaussian, 0.5*(q84-q16) = 0.6745*sigma.
  const double rs = 0.5 * (q84 - q16);
  out.robustSigma = rs * 1.4826;
  if (rs < 0.5)
  {
    out.status = 6;
    out.reason = "degenerate distribution (robust width < 0.5 ADC)";
    out.mean = q50;
    out.sigma = out.robustSigma;
    return out;
  }

  // Peak search restricted around the robust center so that spikes far
  // away (stuck ADC codes, other pick-up) cannot win the amplitude seed.
  const int binLo = h->FindBin(q50 - 3.0 * rs);
  const int binHi = h->FindBin(q50 + 3.0 * rs);
  int peakBin = -1;
  double peakContent = 0.0;
  for (int b = binLo; b <= binHi; ++b)
  {
    const double c = h->GetBinContent(b);
    if (c > peakContent)
    {
      peakContent = c;
      peakBin = b;
    }
  }
  if (peakBin < 0 || peakContent <= 0.0)
  {
    out.status = 2;
    out.reason = "no peak found in robust window";
    out.mean = q50;
    out.sigma = out.robustSigma;
    return out;
  }

  const double peakX = h->GetBinCenter(peakBin);
  const double fitLow = peakX - 5.0 * rs;
  const double fitHigh = peakX + 5.0 * rs;

  // Work on a clone so that the stored QA histogram is never modified by
  // the outlier trimming below.
  TH1D *hc = static_cast<TH1D *>(h->Clone(Form("%s_trim", h->GetName())));
  TF1 *fg = new TF1(Form("robustgaus_%s", h->GetName()), "gaus",
                    fitLow, fitHigh);
  fg->SetParameters(peakContent, peakX, rs);
  fg->SetParLimits(0, 0.0, peakContent * 100.0 + 100.0);
  fg->SetParLimits(1, q50 - 8.0 * rs, q50 + 8.0 * rs);
  fg->SetParLimits(2, 0.2, 100.0);

  // Fit, then remove outlier bins (pull > 3) and refit. This removes a
  // small secondary peak sitting next to the main noise peak.
  TFitResultPtr r;
  int trimRounds = 0;
  for (int round = 0; round <= 2; ++round)
  {
    r = hc->Fit(fg, "QS", "", fitLow, fitHigh);
    if (!(r.Get() && r->IsValid()))
      break;
    if (round == 2)
      break;

    int trimmed = 0;
    for (int b = hc->FindBin(fitLow); b <= hc->FindBin(fitHigh); ++b)
    {
      const double c = hc->GetBinContent(b);
      if (c <= 0.0)
        continue;
      const double x = hc->GetBinCenter(b);
      const double e = hc->GetBinError(b);
      const double pull =
          std::fabs(c - fg->Eval(x)) / (e > 0.0 ? e : 1.0);
      if (pull > 3.0)
      {
        hc->SetBinContent(b, 0.0);
        hc->SetBinError(b, 1.0e12);
        ++trimmed;
      }
    }
    if (trimmed == 0)
      break;
    ++trimRounds;
  }

  if (!(r.Get() && r->IsValid()))
  {
    out.status = 2;
    out.reason = "MINUIT fit failed or invalid result";
    out.mean = q50;
    out.sigma = out.robustSigma;
    delete fg;
    delete hc;
    return out;
  }

  const double mean = r->Parameter(1);
  const double sigma = r->Parameter(2);
  out.chi2Ndf = r->Ndf() > 0 ? r->Chi2() / r->Ndf() : -1.0;

  TString reason;
  int status = 0;
  if (std::isnan(sigma) || sigma <= 0.0)
  {
    status = 3;
    reason = "fitted sigma unphysical (<=0 or NaN)";
  }
  else if (sigma > kMaxNoiseSigma)
  {
    status = 3;
    reason = Form("fitted sigma %.1f exceeds limit %.0f ADC",
                  sigma, kMaxNoiseSigma);
  }
  else if (std::fabs(mean - q50) > 5.0 * rs)
  {
    status = 4;
    reason = Form("fitted mean %.1f far from robust center %.1f "
                  "(wrong peak selected)", mean, q50);
  }

  if (status != 0 && trimRounds > 0)
  {
    status = 5;
    reason = TString("contamination persists after trimming: ") + reason;
  }

  if (status == 0)
  {
    out.mean = mean;
    out.sigma = sigma;

    // Attach the final fit to the stored histogram so that the fitted
    // curve is saved into the QA ROOT file together with the data. The
    // histogram owns the function. Failed fits are not attached.
    TF1 *fitDraw =
        new TF1(Form("gausfit_gc%d", gc), "gaus", fitLow, fitHigh);
    fitDraw->SetParameters(r->Parameter(0), r->Parameter(1),
                           r->Parameter(2));
    fitDraw->SetLineColor(kRed);
    fitDraw->SetLineWidth(2);
    h->GetListOfFunctions()->Add(fitDraw);
  }
  else
  {
    out.mean = q50;
    out.sigma = out.robustSigma;
  }
  out.status = status;
  out.reason = reason;

  delete fg;
  delete hc;
  return out;
}

bool BuildNoiseCalibration(TString rawName,
                           TString outputName,
                           TString detailOutputName,
                           TString padDetailOutputName,
                           TString globalOutputName)
{
  vector<int> apvMap;
  if (!LoadAPVMap(apvMap))
    return false;

  const int nTrees = 2 * MMdetN;
  const int trackerChannelCount = nTrees * 256;

  vector<PadMapRow> padMap;
  TString padMapName;
  const bool hasPadMap = LoadPadMap(padMap, padMapName);
  if (hasPadMap)
  {
    const int bmChannelsPerChip = 121;
    for (size_t i = 0; i < padMap.size(); ++i)
    {
      if (padMap[i].board < 0 || padMap[i].position < 0)
        continue;
      padMap[i].globalChannel =
          trackerChannelCount +
          padMap[i].board * bmChannelsPerChip +
          padMap[i].position + 1;
    }

    // Global BM channels follow physical position i, matching the Tracker
    // convention of writing mapped strip order rather than raw APV order.
    sort(padMap.begin(), padMap.end(),
         [](const PadMapRow &a, const PadMapRow &b) {
           return a.globalChannel < b.globalChannel;
         });

    cout << " --> Loaded BM noise map: " << padMapName
         << " (" << padMap.size() << " rows)" << endl;
  }
  else
    cerr << "WARNING: BM noise map not found; chip12/13 will be skipped."
         << endl;

  const int globalChannelCount =
      trackerChannelCount + static_cast<int>(padMap.size());
  if (globalChannelCount != ::globalNoiseChannelCount)
    cerr << "WARNING: global noise channel count is "
         << globalChannelCount << " (expected "
         << ::globalNoiseChannelCount << ")." << endl;

  TString globalQaName = globalOutputName;
  globalQaName.ReplaceAll(".txt", ".root");
  TFile globalQaFile(globalQaName, "RECREATE");
  if (globalQaFile.IsZombie())
  {
    cerr << "ERROR: cannot write " << globalQaName << endl;
    return false;
  }

  TH1D hGlobalPedestal("global_pedestal_mean",
                       "Pedestal mean;Global channel;ADC",
                       globalChannelCount, 0, globalChannelCount);
  TH1D hGlobalRMS("global_noise_rms",
                  "Noise RMS;Global channel;ADC",
                  globalChannelCount, 0, globalChannelCount);
  TH1D hGlobalPedestalRaw(
      "global_pedestal_mean_raw",
      "Raw pedestal mean (no fit or trimming);Global channel;ADC",
      globalChannelCount, 0, globalChannelCount);
  TH1D hGlobalRMSRaw(
      "global_noise_rms_raw",
      "Raw noise RMS (no fit or trimming);Global channel;ADC",
      globalChannelCount, 0, globalChannelCount);
  TH1D hGlobalEntries("global_noise_sample_count",
                      "Noise samples;Global channel;Samples",
                      globalChannelCount, 0, globalChannelCount);
  TH1D hPedestalDistribution(
      "pedestal_distribution",
      "Distribution of channel pedestal means;Pedestal mean (ADC);Channels",
      4096, -2048, 2048);
  TH1D hRMSDistribution(
      "rms_distribution",
      "Distribution of channel noise RMS;RMS (ADC);Channels",
      4096, 0, 2048);

  vector<TH1D *> waveformHistograms(globalChannelCount, nullptr);
  for (int globalChannel = 1;
       globalChannel <= globalChannelCount; ++globalChannel)
  {
    TString name = Form("pedestal_waveform_global_%04d", globalChannel);
    TString title = Form("Noise waveform samples, global channel %d;ADC;Samples",
                         globalChannel);
    waveformHistograms[globalChannel - 1] =
        new TH1D(name, title, 4096, -2048, 2048);
  }

  TTree waveformTree(
      "NoiseWaveforms",
      "Raw noise waveform vectors associated with global channels");
  int waveformGlobalChannel = -1;
  int waveformEvent = -1;
  int waveformBoard = -1;
  int waveformChip = -1;
  int waveformRawChannel = -1;
  vector<short> waveformSamples;
  const int maxWaveformEntriesPerChannel = 100;
  vector<int> waveformTreeEntries(globalChannelCount, 0);
  waveformTree.Branch("global_channel", &waveformGlobalChannel,
                      "global_channel/I");
  waveformTree.Branch("event", &waveformEvent, "event/I");
  waveformTree.Branch("board", &waveformBoard, "board/I");
  waveformTree.Branch("chip", &waveformChip, "chip/I");
  waveformTree.Branch("raw_channel", &waveformRawChannel,
                      "raw_channel/I");
  waveformTree.Branch("wave", &waveformSamples);

  TFile rawFile(rawName, "READ");
  if (rawFile.IsZombie())
  {
    cerr << "ERROR: cannot open " << rawName << endl;
    return false;
  }

  TTree *rawTree = dynamic_cast<TTree *>(rawFile.Get("TrackRaw"));
  if (!rawTree)
  {
    cerr << "ERROR: TrackRaw is missing in " << rawName << endl;
    return false;
  }

  int event = 0;
  int board = -1;
  int readoutBoard = -1;
  int chip = -1;
  int channel = -1;
  vector<short> *wave = nullptr;

  rawTree->SetBranchAddress("event", &event);
  rawTree->SetBranchAddress("board", &board);
  rawTree->SetBranchAddress("chip", &chip);
  rawTree->SetBranchAddress("channel", &channel);
  rawTree->SetBranchAddress("wave", &wave);
  if (rawTree->GetBranch("readoutBoard"))
    rawTree->SetBranchAddress("readoutBoard", &readoutBoard);

  vector<NoiseAccumulator> accum(trackerChannelCount);
  map<int, NoiseAccumulator> padAccum;
  map<int, PadDetail> padDetails;

  long long validWaveforms = 0;
  long long invalidWaveforms = 0;
  long long padRecords = 0;
  long long padMappedWaveforms = 0;
  long long padUnmappedWaveforms = 0;

  globalQaFile.cd();
  const Long64_t nEntries = rawTree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i)
  {
    rawTree->GetEntry(i);

    if (chip == 12 || chip == 13)
    {
      ++padRecords;
      if (!wave || wave->empty() || !hasPadMap || readoutBoard < 0)
      {
        ++padUnmappedWaveforms;
        continue;
      }

      const PadMapRow *mapRow =
          FindPadMapRow(padMap, readoutBoard, chip, channel);
      if (!mapRow)
      {
        ++padUnmappedWaveforms;
        continue;
      }

      PadDetail detail;
      detail.readoutBoard = readoutBoard;
      detail.rawChip = chip;
      detail.rawChannel = channel;
      detail.connectorPin = mapRow->connectorPin;
      detail.position = mapRow->position;
      detail.mapBoard = mapRow->board;
      detail.globalChannel = mapRow->globalChannel;

      padDetails[detail.globalChannel] = detail;
      NoiseAccumulator &channelAccumulator =
          padAccum[detail.globalChannel];
      for (short sample : *wave)
      {
        channelAccumulator.Add(static_cast<double>(sample));
        waveformHistograms[detail.globalChannel - 1]->Fill(sample);
      }

      if (waveformTreeEntries[detail.globalChannel - 1] <
          maxWaveformEntriesPerChannel)
      {
        waveformGlobalChannel = detail.globalChannel;
        waveformEvent = event;
        waveformBoard = readoutBoard;
        waveformChip = chip;
        waveformRawChannel = channel;
        waveformSamples = *wave;
        waveformTree.Fill();
        ++waveformTreeEntries[detail.globalChannel - 1];
      }

      ++padMappedWaveforms;
      continue;
    }

    if (!wave || wave->empty() ||
        board < 0 || board >= MMdetN ||
        chip < 0 || chip >= apvChipsPerTracker)
    {
      ++invalidWaveforms;
      continue;
    }

    const int strip = MapRawChannel(apvMap, board, chip, channel);
    if (strip < 0 || strip >= 256)
    {
      ++invalidWaveforms;
      continue;
    }

    const int treeIndex = board * 2 + chip / 2;
    const int globalChannel = treeIndex * 256 + strip + 1;
    NoiseAccumulator &channelAccumulator =
        accum[treeIndex * 256 + strip];

    for (short sample : *wave)
    {
      channelAccumulator.Add(static_cast<double>(sample));
      waveformHistograms[globalChannel - 1]->Fill(sample);
    }

    if (waveformTreeEntries[globalChannel - 1] <
        maxWaveformEntriesPerChannel)
    {
      waveformGlobalChannel = globalChannel;
      waveformEvent = event;
      waveformBoard = board;
      waveformChip = chip;
      waveformRawChannel = channel;
      waveformSamples = *wave;
      waveformTree.Fill();
      ++waveformTreeEntries[globalChannel - 1];
    }

    ++validWaveforms;
  }

  if (validWaveforms == 0 && padMappedWaveforms == 0)
  {
    cerr << "ERROR: no usable mapped waveforms in noise input." << endl;
    return false;
  }

  // Gaussian fit over ALL samples of each unified global channel. A robust
  // quantile-based seed plus a restricted fit window and chi2-trimming keep
  // the fit on the true noise peak even when contaminated by far ADC spikes
  // or small secondary peaks. On failure the robust statistics (median and
  // rescaled quantile width) are kept as fallback.
  vector<double> fitMean(globalChannelCount, 0.0);
  vector<double> fitSigma(globalChannelCount, 0.0);
  vector<int> fitStatus(globalChannelCount, -1);
  vector<double> fitChi2Ndf(globalChannelCount, 0.0);
  vector<double> robustMean(globalChannelCount, 0.0);
  vector<double> robustSigma(globalChannelCount, 0.0);
  vector<TString> fitReason(globalChannelCount, "");
  int fitOkCount = 0;
  int fitFallbackCount = 0;
  for (int gc = 1; gc <= globalChannelCount; ++gc)
  {
    const NoiseAccumulator *a = nullptr;
    if (gc <= trackerChannelCount)
      a = &accum[gc - 1];
    else
    {
      map<int, NoiseAccumulator>::const_iterator it = padAccum.find(gc);
      if (it != padAccum.end())
        a = &it->second;
    }
    if (!a || a->count == 0)
    {
      fitStatus[gc - 1] = 1;
      fitReason[gc - 1] = "no samples";
      continue;
    }

    TH1D *h = waveformHistograms[gc - 1];
    const NoiseFitResult fit = FitNoiseGaussian(h, a->count, gc);
    fitMean[gc - 1] = fit.mean;
    fitSigma[gc - 1] = fit.sigma;
    fitStatus[gc - 1] = fit.status;
    fitChi2Ndf[gc - 1] = fit.chi2Ndf;
    robustMean[gc - 1] = fit.robustMean;
    robustSigma[gc - 1] = fit.robustSigma;
    fitReason[gc - 1] = fit.reason;
    if (fit.status == 0)
      ++fitOkCount;
    else
      ++fitFallbackCount;
  }
  cout << " --> gaussian-fit channels: " << fitOkCount
       << ", fallback: " << fitFallbackCount << endl;

  // Report channels whose fit failed together with the failure reason.
  {
    int statusCounts[7] = {0, 0, 0, 0, 0, 0, 0};
    int failedPrinted = 0;
    for (int gc = 1; gc <= globalChannelCount; ++gc)
    {
      const int status = fitStatus[gc - 1];
      if (status >= 0 && status < 7)
        ++statusCounts[status];
      if (status == 0)
        continue;

      TString det;
      const NoiseAccumulator *a = nullptr;
      if (gc <= trackerChannelCount)
      {
        const int treeIndex = (gc - 1) / 256;
        det = Form("Tracker board %d chip %d-%d strip %d",
                   treeIndex / 2, (treeIndex % 2) * 2,
                   (treeIndex % 2) * 2 + 1, (gc - 1) % 256);
        a = &accum[gc - 1];
      }
      else
      {
        det = Form("BM global channel %d", gc);
        map<int, NoiseAccumulator>::const_iterator it = padAccum.find(gc);
        if (it != padAccum.end())
          a = &it->second;
        for (const PadMapRow &row : padMap)
        {
          if (row.globalChannel == gc)
          {
            det = Form("BM board %d chip %d position %d", row.board,
                       row.board == 0 ? 12 : 13, row.position);
            break;
          }
        }
      }

      cout << "  [fit FAIL] gc=" << gc << " " << det
           << " samples=" << (a ? a->count : 0)
           << " status=" << status
           << " reason: " << fitReason[gc - 1]
           << " -> fallback mean=" << fitMean[gc - 1]
           << " sigma=" << fitSigma[gc - 1] << endl;
      ++failedPrinted;
    }
    if (failedPrinted == 0)
      cout << "  --> all channels fitted successfully." << endl;
    else
    {
      cout << "  --> failed channel count: " << failedPrinted
           << " (by status:";
      for (int s = 1; s < 7; ++s)
        cout << " " << s << "=" << statusCounts[s];
      cout << ")" << endl;
    }
  }

  ofstream parameterFile(outputName.Data());
  if (!parameterFile.is_open())
  {
    cerr << "ERROR: cannot write " << outputName << endl;
    return false;
  }

  ofstream detailFile(detailOutputName.Data());
  if (!detailFile.is_open())
  {
    cerr << "ERROR: cannot write " << detailOutputName << endl;
    return false;
  }

  ofstream padDetailFile(padDetailOutputName.Data());
  if (!padDetailFile.is_open())
  {
    cerr << "ERROR: cannot write " << padDetailOutputName << endl;
    return false;
  }

  ofstream globalFile(globalOutputName.Data());
  if (!globalFile.is_open())
  {
    cerr << "ERROR: cannot write " << globalOutputName << endl;
    return false;
  }

  parameterFile << "# Channel parameters from APV noise data\n";
  parameterFile << "# Format: global_channel mean sigma\n";
  parameterFile << fixed << setprecision(9);

  detailFile << "# Per-channel Tracker noise parameters\n";
  detailFile << "# Format: board chip raw_channel mapped_strip global_channel "
                "samples mean sigma\n";
  detailFile << fixed << setprecision(9);

  padDetailFile << "# Per-channel BM pad noise parameters\n";
  padDetailFile << "# Format: readout_board raw_chip raw_channel "
                   "connector_pin position x_mm y_mm samples mean sigma\n";
  padDetailFile << fixed << setprecision(9);

  globalFile << "# Unified noise parameters for global channels 1.."
             << globalChannelCount << "\n";
  globalFile << "# Tracker: 1.." << trackerChannelCount
             << "; BM: " << trackerChannelCount + 1 << ".."
             << globalChannelCount << " in BM position order\n";
  globalFile << "# Format: global_channel mean sigma detector board chip "
                "raw_channel mapped_channel samples\n";
  globalFile << fixed << setprecision(9);

  TTree globalTree(
      "GlobalNoiseChannels",
      "Pedestal and RMS for every unified global channel");
  int globalDetectorType = -1;
  int globalBoard = -1;
  int globalChip = -1;
  int globalRawChannel = -1;
  int globalMappedChannel = -1;
  int globalChannel = -1;
  Long64_t globalSamples = 0;
  double globalMean = 0.0;
  double globalSigma = 0.0;
  double globalRawMean = 0.0;
  double globalRawSigma = 0.0;
  int globalFitStatus = -1;
  double globalFitChi2Ndf = 0.0;
  double globalRobustMean = 0.0;
  double globalRobustSigma = 0.0;

  globalTree.Branch("global_channel", &globalChannel,
                    "global_channel/I");
  globalTree.Branch("detector_type", &globalDetectorType,
                    "detector_type/I");
  globalTree.Branch("board", &globalBoard, "board/I");
  globalTree.Branch("chip", &globalChip, "chip/I");
  globalTree.Branch("raw_channel", &globalRawChannel,
                    "raw_channel/I");
  globalTree.Branch("mapped_channel", &globalMappedChannel,
                    "mapped_channel/I");
  globalTree.Branch("samples", &globalSamples, "samples/L");
  globalTree.Branch("mean", &globalMean, "mean/D");
  globalTree.Branch("sigma", &globalSigma, "sigma/D");
  // Arithmetic moments of every mapped raw ADC sample. These values are
  // deliberately independent of the robust seed, Gaussian fit and pull
  // trimming used for the production mean/sigma branches.
  globalTree.Branch("raw_mean", &globalRawMean, "raw_mean/D");
  globalTree.Branch("raw_sigma", &globalRawSigma, "raw_sigma/D");
  globalTree.Branch("fit_status", &globalFitStatus, "fit_status/I");
  globalTree.Branch("fit_chi2_ndf", &globalFitChi2Ndf, "fit_chi2_ndf/D");
  globalTree.Branch("robust_mean", &globalRobustMean, "robust_mean/D");
  globalTree.Branch("robust_sigma", &globalRobustSigma, "robust_sigma/D");

  TH1D hMean("noise_pedestal_mean", "Noise pedestal mean;Global channel;ADC",
             nTrees * 256, 0, nTrees * 256);
  TH1D hRMS("noise_pedestal_rms", "Noise pedestal sigma;Global channel;ADC",
            nTrees * 256, 0, nTrees * 256);
  TH1D hEntries("noise_sample_count", "Noise sample count;Global channel;Samples",
                nTrees * 256, 0, nTrees * 256);

  TTree channelTree("NoiseChannels",
                    "Per board/chip/raw-channel Tracker noise parameters");
  int outBoard = -1;
  int outChip = -1;
  int outRawChannel = -1;
  int outStrip = -1;
  int outGlobalChannel = -1;
  Long64_t outSamples = 0;
  double outMean = 0.0;
  double outSigma = 0.0;

  channelTree.Branch("board", &outBoard, "board/I");
  channelTree.Branch("chip", &outChip, "chip/I");
  channelTree.Branch("raw_channel", &outRawChannel, "raw_channel/I");
  channelTree.Branch("mapped_strip", &outStrip, "mapped_strip/I");
  channelTree.Branch("global_channel", &outGlobalChannel,
                     "global_channel/I");
  channelTree.Branch("samples", &outSamples, "samples/L");
  channelTree.Branch("mean", &outMean, "mean/D");
  channelTree.Branch("sigma", &outSigma, "sigma/D");

  TTree padTree("PadNoiseChannels",
                "Per BM pad channel noise parameters");
  int outReadoutBoard = -1;
  int outRawChip = -1;
  int outPadRawChannel = -1;
  int outConnectorPin = -1;
  int outPosition = -1;
  Long64_t outPadSamples = 0;
  double outX = 0.0;
  double outY = 0.0;
  double outPadMean = 0.0;
  double outPadSigma = 0.0;

  padTree.Branch("readout_board", &outReadoutBoard, "readout_board/I");
  padTree.Branch("raw_chip", &outRawChip, "raw_chip/I");
  padTree.Branch("raw_channel", &outPadRawChannel, "raw_channel/I");
  padTree.Branch("connector_pin", &outConnectorPin, "connector_pin/I");
  padTree.Branch("position", &outPosition, "position/I");
  padTree.Branch("x_mm", &outX, "x_mm/D");
  padTree.Branch("y_mm", &outY, "y_mm/D");
  padTree.Branch("samples", &outPadSamples, "samples/L");
  padTree.Branch("mean", &outPadMean, "mean/D");
  padTree.Branch("sigma", &outPadSigma, "sigma/D");

  long long calibratedChannels = 0;
  for (int treeIndex = 0; treeIndex < nTrees; ++treeIndex)
  {
    for (int strip = 0; strip < 256; ++strip)
    {
      const int globalChannel = treeIndex * 256 + strip + 1;
      const NoiseAccumulator &a = accum[treeIndex * 256 + strip];
      const double mean = a.count > 0 ? fitMean[globalChannel - 1] : 0.0;
      const double sigma = a.count > 0 ? fitSigma[globalChannel - 1] : 0.0;

      parameterFile << globalChannel << " " << mean << " "
                    << sigma << "\n";
      hMean.SetBinContent(globalChannel, mean);
      hRMS.SetBinContent(globalChannel, sigma);
      hEntries.SetBinContent(globalChannel,
                             static_cast<double>(a.count));
      if (a.count > 0)
        ++calibratedChannels;
    }
  }

  for (int boardIndex = 0; boardIndex < MMdetN; ++boardIndex)
  {
    for (int chipIndex = 0; chipIndex < apvChipsPerTracker; ++chipIndex)
    {
      const int treeIndex = boardIndex * 2 + chipIndex / 2;
      for (int rawChannel = 0;
           rawChannel < apvChannelsPerChip; ++rawChannel)
      {
        const int strip =
            MapRawChannel(apvMap, boardIndex, chipIndex, rawChannel);
        const int globalChannel =
            strip >= 0 ? treeIndex * 256 + strip + 1 : -1;

        const NoiseAccumulator *a = nullptr;
        if (strip >= 0 && strip < 256)
          a = &accum[treeIndex * 256 + strip];

        const long long samples = a ? a->count : 0;
        const double mean = a && a->count > 0 && globalChannel >= 1
                                ? fitMean[globalChannel - 1]
                                : 0.0;
        const double sigma = a && a->count > 0 && globalChannel >= 1
                                 ? fitSigma[globalChannel - 1]
                                 : 0.0;

        detailFile << boardIndex << " " << chipIndex << " "
                   << rawChannel << " " << strip << " "
                   << globalChannel << " " << samples << " "
                   << mean << " " << sigma << "\n";

        outBoard = boardIndex;
        outChip = chipIndex;
        outRawChannel = rawChannel;
        outStrip = strip;
        outGlobalChannel = globalChannel;
        outSamples = samples;
        outMean = mean;
        outSigma = sigma;
        channelTree.Fill();
      }
    }
  }

  long long calibratedPadChannels = 0;
  for (const auto &item : padDetails)
  {
    const int key = item.first;
    const PadDetail &detail = item.second;
    const NoiseAccumulator &a = padAccum[key];
    const int position = detail.position;
    const double x = (position % 11) * 9.0 + 4.5;
    const double y = (position / 11) * 9.0 + 4.5;
    const double mean = a.count > 0 ? fitMean[key - 1] : 0.0;
    const double sigma = a.count > 0 ? fitSigma[key - 1] : 0.0;

    padDetailFile << detail.readoutBoard << " "
                  << detail.rawChip << " "
                  << detail.rawChannel << " "
                  << detail.connectorPin << " "
                  << position << " " << x << " " << y << " "
                  << a.count << " " << mean << " " << sigma << "\n";

    outReadoutBoard = detail.readoutBoard;
    outRawChip = detail.rawChip;
    outPadRawChannel = detail.rawChannel;
    outConnectorPin = detail.connectorPin;
    outPosition = position;
    outX = x;
    outY = y;
    outPadSamples = a.count;
    outPadMean = mean;
    outPadSigma = sigma;
    padTree.Fill();

    if (a.count > 0)
      ++calibratedPadChannels;
  }

  for (int treeIndex = 0; treeIndex < nTrees; ++treeIndex)
  {
    for (int strip = 0; strip < 256; ++strip)
    {
      const int globalChannelNumber = treeIndex * 256 + strip + 1;
      const NoiseAccumulator &a =
          accum[treeIndex * 256 + strip];
      const double mean = a.count > 0 ? fitMean[globalChannelNumber - 1]
                                      : 0.0;
      const double sigma = a.count > 0 ? fitSigma[globalChannelNumber - 1]
                                       : 0.0;
      const double rawMean = a.count > 0 ? a.mean : 0.0;
      const double rawSigma = a.count > 0 ? a.Sigma() : 0.0;

      hGlobalPedestal.SetBinContent(globalChannelNumber, mean);
      hGlobalRMS.SetBinContent(globalChannelNumber, sigma);
      hGlobalPedestalRaw.SetBinContent(globalChannelNumber, rawMean);
      hGlobalRMSRaw.SetBinContent(globalChannelNumber, rawSigma);
      hGlobalEntries.SetBinContent(globalChannelNumber,
                                   static_cast<double>(a.count));
      if (a.count > 0)
      {
        hPedestalDistribution.Fill(mean);
        hRMSDistribution.Fill(sigma);
      }

      globalFile << globalChannelNumber << " " << mean << " "
                 << sigma << " Tracker " << (treeIndex / 2) << " "
                 << ((treeIndex % 2) * 2) << " -1 " << strip << " "
                 << a.count << "\n";

      globalChannel = globalChannelNumber;
      globalDetectorType = 0;
      globalBoard = treeIndex / 2;
      globalChip = (treeIndex % 2) * 2;
      globalRawChannel = -1;
      globalMappedChannel = strip;
      globalSamples = a.count;
      globalMean = mean;
      globalSigma = sigma;
      globalRawMean = rawMean;
      globalRawSigma = rawSigma;
      globalFitStatus = fitStatus[globalChannelNumber - 1];
      globalFitChi2Ndf = fitChi2Ndf[globalChannelNumber - 1];
      globalRobustMean = robustMean[globalChannelNumber - 1];
      globalRobustSigma = robustSigma[globalChannelNumber - 1];
      globalTree.Fill();
    }
  }

  for (const PadMapRow &row : padMap)
  {
    map<int, NoiseAccumulator>::const_iterator accumulatorIt =
        padAccum.find(row.globalChannel);
    const NoiseAccumulator *a =
        accumulatorIt != padAccum.end() ? &accumulatorIt->second : nullptr;
    const long long samples = a ? a->count : 0;
    const double mean = a && a->count > 0 ? fitMean[row.globalChannel - 1]
                                          : 0.0;
    const double sigma = a && a->count > 0 ? fitSigma[row.globalChannel - 1]
                                           : 0.0;
    const double rawMean = a && a->count > 0 ? a->mean : 0.0;
    const double rawSigma = a && a->count > 0 ? a->Sigma() : 0.0;

    PadDetail detail;
    map<int, PadDetail>::const_iterator detailIt =
        padDetails.find(row.globalChannel);
    if (detailIt != padDetails.end())
      detail = detailIt->second;
    else
    {
      detail.readoutBoard = row.board + 1;
      detail.rawChip = row.board == 0 ? 12 : 13;
      detail.rawChannel = row.channelId > 0 ? row.channelId - 1
                                            : row.channelId;
      detail.connectorPin = row.connectorPin;
      detail.position = row.position;
      detail.mapBoard = row.board;
      detail.globalChannel = row.globalChannel;
    }

    const double x = (row.position % 11) * 9.0 + 4.5;
    const double y = (row.position / 11) * 9.0 + 4.5;
    hGlobalPedestal.SetBinContent(row.globalChannel, mean);
    hGlobalRMS.SetBinContent(row.globalChannel, sigma);
    hGlobalPedestalRaw.SetBinContent(row.globalChannel, rawMean);
    hGlobalRMSRaw.SetBinContent(row.globalChannel, rawSigma);
    hGlobalEntries.SetBinContent(row.globalChannel,
                                 static_cast<double>(samples));
    if (samples > 0)
    {
      hPedestalDistribution.Fill(mean);
      hRMSDistribution.Fill(sigma);
    }

    globalFile << row.globalChannel << " " << mean << " "
               << sigma << " BM " << row.board << " "
               << detail.rawChip << " " << detail.rawChannel << " "
               << row.position << " " << samples << "\n";

    globalChannel = row.globalChannel;
    globalDetectorType = 1;
    globalBoard = row.board;
    globalChip = detail.rawChip;
    globalRawChannel = detail.rawChannel;
    globalMappedChannel = row.position;
    globalSamples = samples;
    globalMean = mean;
    globalSigma = sigma;
    globalRawMean = rawMean;
    globalRawSigma = rawSigma;
    globalFitStatus = fitStatus[row.globalChannel - 1];
    globalFitChi2Ndf = fitChi2Ndf[row.globalChannel - 1];
    globalRobustMean = robustMean[row.globalChannel - 1];
    globalRobustSigma = robustSigma[row.globalChannel - 1];
    globalTree.Fill();
  }

  globalFile.close();
  globalQaFile.cd();
  hGlobalPedestal.Write();
  hGlobalRMS.Write();
  hGlobalPedestalRaw.Write();
  hGlobalRMSRaw.Write();
  hGlobalEntries.Write();
  hPedestalDistribution.Write();
  hRMSDistribution.Write();
  for (size_t i = 0; i < waveformHistograms.size(); ++i)
    waveformHistograms[i]->Write();
  waveformTree.Write();
  globalTree.Write();
  globalQaFile.Close();

  parameterFile.close();
  detailFile.close();
  padDetailFile.close();

  TString qaName = outputName;
  qaName.ReplaceAll(".txt", ".root");
  TFile qaFile(qaName, "RECREATE");
  if (qaFile.IsZombie())
  {
    cerr << "ERROR: cannot write " << qaName << endl;
    return false;
  }
  hMean.Write();
  hRMS.Write();
  hEntries.Write();
  channelTree.Write();
  padTree.Write();
  qaFile.Close();

  cout << "Noise calibration finished." << endl;
  cout << "  raw entries: " << nEntries << endl;
  cout << "  valid Tracker waveforms: " << validWaveforms << endl;
  cout << "  invalid Tracker waveforms: " << invalidWaveforms << endl;
  cout << "  Tracker calibrated channels: " << calibratedChannels
       << " / " << nTrees * 256 << endl;
  cout << "  BM records: " << padRecords
       << ", mapped waveforms: " << padMappedWaveforms
       << ", unmapped: " << padUnmappedWaveforms << endl;
  cout << "  BM calibrated channels: " << calibratedPadChannels << endl;
  cout << "  flat parameters: " << outputName << endl;
  cout << "  Tracker detail: " << detailOutputName << endl;
  cout << "  BM detail: " << padDetailOutputName << endl;
  cout << "  QA ROOT: " << qaName << endl;
  return true;
}
}

bool CalibrateTrackerNoise(TString fileDir, int force)
{
  namespace fs = std::filesystem;
  vector<TString> rootList;
  fs::path input = fs::absolute(fileDir.Data()).lexically_normal();
  fs::path outputDirectory;
  if (input.extension() == ".root")
  {
    if (!fs::is_regular_file(input))
    {
      cerr << "ERROR: noise input does not exist: " << fileDir << endl;
      return false;
    }
    rootList.emplace_back(input.string());
    // Explicit standalone noise files have their own output directory.
    outputDirectory = input.parent_path() / (input.stem().string() + "-noise");
  }
  else
  {
    if (input.filename() == "noise" && input.parent_path().filename() == "Tracker")
      input = input.parent_path().parent_path();
    else if (input.filename() == "Tracker")
      input = input.parent_path();
    if (!CollectRunInputs((input / "Tracker/noise").string(), rootList))
      return false;
    outputDirectory = input / "Combine/noise";
  }
  if (!EnsureOutputDirectory(outputDirectory.string()))
    return false;

  const TString outputDir = outputDirectory.string();
  const TString rawName = outputDir + "/Tracker-noise-raw.root";
  const TString parameters = outputDir + "/fit_parameters.txt";
  const TString detailName = outputDir + "/Tracker-noise-channel-parameters.txt";
  const TString padDetailName = outputDir + "/BM-noise-channel-parameters.txt";
  const TString globalName = outputDir + "/noise-global-parameters-1778.txt";

  CRConfig config;
  if (!config.Load())
    return false;
  const int bmSlot8LogicalAPV =
      static_cast<int>(config.GetDouble("bm_slot8_logical_chip", 12));
  const int bmSlot9LogicalAPV =
      static_cast<int>(config.GetDouble("bm_slot9_logical_chip", 13));
  cout << " --> BM slot mapping: slot8->chip" << bmSlot8LogicalAPV
       << ", slot9->chip" << bmSlot9LogicalAPV << endl;

  if (force != 1 && !gSystem->AccessPathName(rawName))
  {
    if (!ValidateRootTrees(rawName, {"TrackRaw"}, true))
      return false;
    cout << " --> Reusing existing raw noise ROOT: " << rawName << endl;
  }
  else
  {
    cout << " --> Normalising decoded noise ROOT: " << rawName << endl;
    DecodedAPVRoot2TrackRaw converter(
        rootList, rawName, 1, false,
        bmSlot8LogicalAPV, bmSlot9LogicalAPV);
    if (!converter.Init() || !converter.Loop())
      return false;
    if (!ValidateRootTrees(rawName, {"TrackRaw"}, true))
      return false;
  }
  return BuildNoiseCalibration(rawName, parameters, detailName,
                               padDetailName, globalName);
}

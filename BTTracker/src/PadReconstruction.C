#include "CRInfo.h"
#include "PadReconstruction.h"
#include "TNamed.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace
{
struct PadMapRow
{
  int board = -1;
  int chip = -1;
  int channelId = -1;
  int connectorPin = -1;
  int position = -1;
};

struct PadHit
{
  int position = -1;
  double charge = 0.0;
  double peak = 0.0;
  double timingAmplitude = -1.0;
  std::array<double, 4> timeSigmoid = {{-999.0, -999.0, -999.0, -999.0}};
};

struct TimeAccumulator
{
  double sum = 0.0;
  double minimum = 1.0e30;
  int count = 0;

  void Add(double value)
  {
    if (!std::isfinite(value) || value <= -900.0)
      return;
    sum += value;
    minimum = std::min(minimum, value);
    ++count;
  }

  double Mean() const
  {
    return count > 0 ? sum / static_cast<double>(count) : -999.0;
  }

  double Minimum() const
  {
    return count > 0 ? minimum : -999.0;
  }
};
struct PadEvent
{
  Long64_t timestamp = -1;
  std::map<int, PadHit> hits;
  TimeAccumulator arrivalT;
  TimeAccumulator arrivalTa;
  TimeAccumulator arrivalTe;
  TimeAccumulator arrivalTSigmoid[4]; // t20, t30, t40, t50
};

// Track-matched timing geometry.
const double kBMZ[2] = {BMZpos[0], BMZpos[1]}; // chip12 (BM1), chip13 (BM2)
const int kBMChip[2] = {12, 13};
const int kPadGrid = 11;
const double kPadPitch = 9.0;
const double kPadCentre = kPadGrid / 2.0 * kPadPitch; // 49.5 mm
const double kPadDxOff = 0.0; // pad centre vs tracker offset (mm)
const double kPadDyOff = 0.0;
const double kStripPitch = 0.4;
const double kStripCentre = 256.0 / 2.0 * kStripPitch; // 51.2 mm
const double kStripWindow = 2.0; // strips around the predicted strip
const int kClusterGap = 1;       // strips with gap > 1 split clusters

struct StripWave
{
  int channel = -1;
  double timeA = -999.0;
  double timeE = -999.0;
  std::array<double, 4> timeS = {{-999.0, -999.0, -999.0, -999.0}};
};

struct PadTimeKey
{
  int event = -1;
  int chip = -1;
  int position = -1;

  bool operator<(const PadTimeKey &other) const
  {
    if (event != other.event)
      return event < other.event;
    if (chip != other.chip)
      return chip < other.chip;
    return position < other.position;
  }
};

bool LoadPadMap(std::vector<PadMapRow> &rows, TString &usedName)
{
  std::ifstream input("config/map_APV4BMM.txt");
  usedName = "config/map_APV4BMM.txt";

  PadMapRow row;
  while (input >> row.board >> row.chip >> row.channelId
               >> row.connectorPin >> row.position)
  {
    if (row.position < 0 || row.position >= 121)
      continue;
    rows.push_back(row);
  }

  if (rows.empty())
  {
    std::cerr << " ERROR! Pad map is empty: " << usedName << std::endl;
    return false;
  }

  std::cout << " --> Loaded pad map " << usedName
            << " with " << rows.size() << " entries" << std::endl;
  return true;
}

const PadMapRow *UniqueMatch(const std::vector<PadMapRow> &rows,
                             int board, int chip, int channelId,
                             bool requireChip, bool &ambiguous)
{
  const PadMapRow *match = nullptr;
  int position = -1;
  ambiguous = false;

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
      ambiguous = true;
      return nullptr;
    }
  }
  return match;
}

const PadMapRow *FindPadMapRow(const std::vector<PadMapRow> &rows,
                               int readoutBoard, int rawChip, int rawChannel)
{
  (void)readoutBoard;
  if ((rawChip != 12 && rawChip != 13) ||
      rawChannel < 0 || rawChannel >= 128)
    return nullptr;

  // The decoded-ROOT adapter normalises BM records to chip 12/13. The chip
  // selects the map board; raw APV channel is zero-based while the BMM map
  // channel_id is one-based.
  const int mapBoard = rawChip == 12 ? 0 : 1;
  const int mapChannelId = rawChannel + 1;
  const PadMapRow *match = nullptr;
  for (const PadMapRow &row : rows)
  {
    if (row.board != mapBoard || row.channelId != mapChannelId)
      continue;
    if (match && match->position != row.position)
      return nullptr;
    match = &row;
  }
  return match;
}


double PadCharge(const std::vector<short> &wave, double baseline,
                 double rms, double &peak, bool &overThreshold)
{
  peak = 0.0;
  overThreshold = false;
  if (wave.empty())
    return 0.0;

  double positivePeak = 0.0;
  double negativePeak = 0.0;
  double positiveCharge = 0.0;
  double negativeCharge = 0.0;

  for (short sample : wave)
  {
    const double delta = static_cast<double>(sample) - baseline;
    positivePeak = std::max(positivePeak, delta);
    negativePeak = std::max(negativePeak, -delta);
    positiveCharge += std::max(0.0, delta);
    negativeCharge += std::max(0.0, -delta);
  }

  peak = std::max(positivePeak, negativePeak);
  overThreshold = peak > 5.0 * std::max(rms, 1.0);
  if (!overThreshold)
    return 0.0;

  if (negativePeak > positivePeak)
    return negativeCharge;
  return positiveCharge;
}

bool Adjacent(int first, int second)
{
  const int firstX = first % 11;
  const int firstY = first / 11;
  const int secondX = second % 11;
  const int secondY = second / 11;
  return std::abs(firstX - secondX) + std::abs(firstY - secondY) == 1;
}

// One independent output buffer per physical BM plane. Pad positions are
// local to a plane, so neither hit merging nor clustering may cross planes.
struct PadPlaneOutput
{
  Long64_t timestamp = -1;
  int nhits = 0;
  int nclusters = 0;
  double timeMean[7] = {-999.0, -999.0, -999.0, -999.0,
                        -999.0, -999.0, -999.0};
  double timeMin[7] = {-999.0, -999.0, -999.0, -999.0,
                       -999.0, -999.0, -999.0};
  int timeCount[7] = {0, 0, 0, 0, 0, 0, 0};
  std::vector<int> hitPosition;
  std::vector<double> hitX, hitY, hitCharge, hitPeak;
  std::vector<int> clusterSize;
  std::vector<double> clusterX, clusterY, clusterCharge;
  std::vector<int> hitClusterIndex;
  std::array<std::vector<int>, 4> clusterTimeSigmoidN;
  std::array<std::vector<double>, 4> hitTimeSigmoid;
  std::array<std::vector<double>, 4> clusterTimeSigmoidMean;
  std::array<std::vector<double>, 4> clusterTimeSigmoidMin;
  double primaryX = -1.0;
  double primaryY = -1.0;
  double primaryCharge = 0.0;

  void Branch(TTree &tree, const TString &prefix)
  {
    const auto scalar = [&](const char *suffix, auto *value, char type) {
      const TString name = prefix + suffix;
      tree.Branch(name, value, Form("%s/%c", name.Data(), type));
    };
    scalar("timestamp", &timestamp, 'L');
    scalar("nhits", &nhits, 'I');
    scalar("nclusters", &nclusters, 'I');
    const char *timeNames[7] = {"arrivalT", "arrivalT_a", "arrivalT_e",
                                "arrivalT_sigmoid", "arrivalT_sigmoid30",
                                "arrivalT_sigmoid40", "arrivalT_sigmoid50"};
    for (int m = 0; m < 7; ++m)
    {
      scalar(Form("%s_mean", timeNames[m]), &timeMean[m], 'D');
      scalar(Form("%s_min", timeNames[m]), &timeMin[m], 'D');
      scalar(Form("n_%s", timeNames[m]), &timeCount[m], 'I');
    }
    // Preserve the spelling aliases within each plane only.
    scalar("arrivetime_mean", &timeMean[0], 'D');
    scalar("arrivetime_min", &timeMin[0], 'D');
    tree.Branch(prefix + "hit_position", &hitPosition);
    tree.Branch(prefix + "hit_x", &hitX);
    tree.Branch(prefix + "hit_y", &hitY);
    tree.Branch(prefix + "hit_charge", &hitCharge);
    tree.Branch(prefix + "hit_peak", &hitPeak);
    tree.Branch(prefix + "cluster_size", &clusterSize);
    tree.Branch(prefix + "cluster_x", &clusterX);
    tree.Branch(prefix + "cluster_y", &clusterY);
    tree.Branch(prefix + "cluster_charge", &clusterCharge);
    tree.Branch(prefix + "hit_cluster_index", &hitClusterIndex);
    const char *fractionTags[4] = {"", "30", "40", "50"};
    for (int f = 0; f < 4; ++f)
    {
      const TString tag = fractionTags[f];
      tree.Branch(prefix + "hit_time_sigmoid" + tag, &hitTimeSigmoid[f]);
      tree.Branch(prefix + "cluster_time_sigmoid" + tag + "_mean",
                  &clusterTimeSigmoidMean[f]);
      tree.Branch(prefix + "cluster_time_sigmoid" + tag + "_min",
                  &clusterTimeSigmoidMin[f]);
      tree.Branch(prefix + "cluster_time_sigmoid" + tag + "_n",
                  &clusterTimeSigmoidN[f]);
    }
    scalar("primary_x", &primaryX, 'D');
    scalar("primary_y", &primaryY, 'D');
    scalar("primary_charge", &primaryCharge, 'D');
  }

  void Set(const PadEvent &event)
  {
    // Reset even when this plane has no hits in an otherwise valid event.
    *this = PadPlaneOutput{};
    timestamp = event.timestamp;
    const TimeAccumulator *times[7] = {&event.arrivalT, &event.arrivalTa,
                                       &event.arrivalTe, &event.arrivalTSigmoid[0],
                                       &event.arrivalTSigmoid[1],
                                       &event.arrivalTSigmoid[2],
                                       &event.arrivalTSigmoid[3]};
    for (int m = 0; m < 7; ++m)
    {
      timeMean[m] = times[m]->Mean();
      timeMin[m] = times[m]->Minimum();
      timeCount[m] = times[m]->count;
    }
    for (const auto &item : event.hits)
    {
      const PadHit &hit = item.second;
      hitPosition.push_back(hit.position);
      hitX.push_back((hit.position % kPadGrid) * kPadPitch + kPadPitch / 2.0);
      hitY.push_back((hit.position / kPadGrid) * kPadPitch + kPadPitch / 2.0);
      hitCharge.push_back(hit.charge);
      hitPeak.push_back(hit.peak);
      for (int f = 0; f < 4; ++f)
        hitTimeSigmoid[f].push_back(hit.timeSigmoid[f]);
    }
    nhits = static_cast<int>(hitPosition.size());
    hitClusterIndex.assign(nhits, -1);
    std::vector<bool> visited(hitPosition.size(), false);
    for (std::size_t seed = 0; seed < hitPosition.size(); ++seed)
    {
      if (visited[seed])
        continue;
      std::vector<std::size_t> component;
      std::vector<std::size_t> pending = {seed};
      visited[seed] = true;
      while (!pending.empty())
      {
        const std::size_t current = pending.back();
        pending.pop_back();
        component.push_back(current);
        for (std::size_t candidate = 0; candidate < hitPosition.size(); ++candidate)
        {
          if (!visited[candidate] && Adjacent(hitPosition[current], hitPosition[candidate]))
          {
            visited[candidate] = true;
            pending.push_back(candidate);
          }
        }
      }
      double charge = 0.0, xCharge = 0.0, yCharge = 0.0;
      TimeAccumulator clusterTime[4];
      for (std::size_t index : component)
      {
        charge += hitCharge[index];
        xCharge += hitCharge[index] * hitX[index];
        yCharge += hitCharge[index] * hitY[index];
        hitClusterIndex[index] = static_cast<int>(clusterSize.size());
        for (int f = 0; f < 4; ++f)
          clusterTime[f].Add(hitTimeSigmoid[f][index]);
      }
      const double x = charge > 0.0 ? xCharge / charge : -1.0;
      const double y = charge > 0.0 ? yCharge / charge : -1.0;
      clusterSize.push_back(static_cast<int>(component.size()));
      clusterX.push_back(x);
      clusterY.push_back(y);
      clusterCharge.push_back(charge);
      for (int f = 0; f < 4; ++f)
      {
        clusterTimeSigmoidMean[f].push_back(clusterTime[f].Mean());
        clusterTimeSigmoidMin[f].push_back(clusterTime[f].Minimum());
        clusterTimeSigmoidN[f].push_back(clusterTime[f].count);
      }
      if (charge > primaryCharge)
      {
        primaryCharge = charge;
        primaryX = x;
        primaryY = y;
      }
    }
    nclusters = static_cast<int>(clusterSize.size());
  }
};
}

bool ReconstructPadEvents(TString step2Name, TString step4Name,
                          TString outputName, bool bm1Enabled, bool bm2Enabled)
{
  const bool bmEnabled[2] = {bm1Enabled, bm2Enabled};
  TFile input(step2Name, "READ");
  if (input.IsZombie())
  {
    std::cerr << " ERROR! Cannot open pad step2 input " << step2Name << std::endl;
    return false;
  }

  TFile trackInput(step4Name, "READ");
  if (trackInput.IsZombie())
  {
    std::cerr << " ERROR! Cannot open track input " << step4Name << std::endl;
    return false;
  }

  TTree *trackTree = dynamic_cast<TTree *>(trackInput.Get("Track"));
  if (!trackTree)
  {
    std::cerr << " ERROR! Track tree is missing in " << step4Name << std::endl;
    return false;
  }

  if (!trackTree->GetBranch("predStripX") || !trackTree->GetBranch("predStripY"))
  {
    std::cerr << " ERROR! Track-matched timing requires predStripX/Y; "
              << "cannot substitute zero predictions." << std::endl;
    return false;
  }
  TTree *step2Tree = dynamic_cast<TTree *>(input.Get("BMPadSignals"));
  if (!step2Tree)
  {
    std::cerr << " ERROR! BMPadSignals tree is missing in " << step2Name << std::endl;
    return false;
  }
  if (!step2Tree->GetBranch("raw_chip"))
  {
    std::cerr << " ERROR! Separating BM1/BM2 requires BMPadSignals.raw_chip "
                 "(normalised chip 12/13)." << std::endl;
    return false;
  }

  // Load the BM pad-plane alignment parameters (dx/dy per chip) fitted
  // by ReconstructTrack against tracker-projected tracks. They are applied
  // to the track projection below so that BM pads are matched in the
  // detector-aligned frame (chip12 -> index 0, chip13 -> index 1).
  double bmAlignDx[2] = {0.0, 0.0};
  double bmAlignDy[2] = {0.0, 0.0};
  {
    TString bmAlignName = step4Name;
    const int lastSlash = bmAlignName.Last('/');
    if (lastSlash >= 0)
      bmAlignName.Remove(lastSlash + 1);
    else
      bmAlignName = "";
    bmAlignName += "MMalignmentParBM.dat";
    std::ifstream bmAlignInput(bmAlignName.Data());
    if (bmAlignInput.is_open())
    {
      std::string line;
      while (std::getline(bmAlignInput, line))
      {
        if (line.empty() || line[0] == '#')
          continue;
        int chip = -1;
        double dx = 0.0, dy = 0.0;
        std::istringstream iss(line);
        if (!(iss >> chip >> dx >> dy))
          continue;
        if (chip == 12)
        {
          bmAlignDx[0] = dx;
          bmAlignDy[0] = dy;
        }
        else if (chip == 13)
        {
          bmAlignDx[1] = dx;
          bmAlignDy[1] = dy;
        }
      }
      std::cout << " --> Loaded BM alignment " << bmAlignName
                << ": chip12 dx=" << bmAlignDx[0] << " dy=" << bmAlignDy[0]
                << " | chip13 dx=" << bmAlignDx[1] << " dy=" << bmAlignDy[1]
                << std::endl;
    }
    else
      std::cerr << " WARNING! BM alignment file not found: " << bmAlignName
                << "; BM pads assumed unaligned." << std::endl;
  }

  int event = -1;
  int channel = -1;
  int globalChannel = -1;
  int position = -1;
  int nisel = 0;
  double ped = 0.0;
  double pedrms = 1.0;
  double arrivalT = -999.0;
  double arrivalTa = -999.0;
  double arrivalTe = -999.0;
  double arrivalTSigmoid[4] = {-999.0, -999.0, -999.0, -999.0};
  Long64_t timestamp = -1;
  std::vector<short> *wave = nullptr;
  int rawChip = -1;
  double maxA = 0.0;

  step2Tree->SetBranchAddress("event", &event);
  step2Tree->SetBranchAddress("channel", &channel);
  step2Tree->SetBranchAddress("global_channel", &globalChannel);
  step2Tree->SetBranchAddress("position", &position);
  step2Tree->SetBranchAddress("ped", &ped);
  step2Tree->SetBranchAddress("pedrms", &pedrms);
  step2Tree->SetBranchAddress("wave", &wave);
  if (step2Tree->GetBranch("timestamp"))
    step2Tree->SetBranchAddress("timestamp", &timestamp);
  if (step2Tree->GetBranch("nisel"))
    step2Tree->SetBranchAddress("nisel", &nisel);
  step2Tree->SetBranchAddress("raw_chip", &rawChip);
  if (step2Tree->GetBranch("maxA"))
    step2Tree->SetBranchAddress("maxA", &maxA);
  if (step2Tree->GetBranch("arrivalT"))
    step2Tree->SetBranchAddress("arrivalT", &arrivalT);
  if (step2Tree->GetBranch("arrivetime_a"))
    step2Tree->SetBranchAddress("arrivetime_a", &arrivalTa);
  else if (step2Tree->GetBranch("arrivalT_a"))
    step2Tree->SetBranchAddress("arrivalT_a", &arrivalTa);
  if (step2Tree->GetBranch("arrivetime_e"))
    step2Tree->SetBranchAddress("arrivetime_e", &arrivalTe);
  else if (step2Tree->GetBranch("arrivalT_e"))
    step2Tree->SetBranchAddress("arrivalT_e", &arrivalTe);
  if (step2Tree->GetBranch("arrivetime_sigmoid"))
    step2Tree->SetBranchAddress("arrivetime_sigmoid", &arrivalTSigmoid[0]);
  else if (step2Tree->GetBranch("arrivalT_sigmoid"))
    step2Tree->SetBranchAddress("arrivalT_sigmoid", &arrivalTSigmoid[0]);
  const char *sigmoidBranches[3] = {"arrivalT_sigmoid30",
                                    "arrivalT_sigmoid40",
                                    "arrivalT_sigmoid50"};
  for (int f = 1; f < 4; ++f)
    if (step2Tree->GetBranch(sigmoidBranches[f - 1]))
      step2Tree->SetBranchAddress(sigmoidBranches[f - 1], &arrivalTSigmoid[f]);

  // (event, chip, position) -> {maxA, a, e, sigmoid t20/t30/t40/t50} for BM
  // timing. Keeps the entry with the largest maxA per pad.
  std::map<PadTimeKey, std::array<double, 7>> padTimes;


  std::map<int, std::array<PadEvent, 2>> events;
  Long64_t padRecords = 0;
    arrivalT = -999.0;
    arrivalTa = -999.0;
    arrivalTe = -999.0;
    for (double &time : arrivalTSigmoid)
      time = -999.0;
  Long64_t selectedRecords = 0;
  Long64_t unmappedRecords = 0;

  const Long64_t entries = step2Tree->GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry)
  {
    step2Tree->GetEntry(entry);

    ++padRecords;
    if (!wave || wave->empty() || (rawChip != 12 && rawChip != 13))
      continue;
    // An explicitly disconnected BM is excluded from pad reconstruction
    // and matched timing. Its raw/step2 records remain available for QA.
    if (!bmEnabled[rawChip - 12])
      continue;

    if (position >= 0 && position < 121)
    {
      PadTimeKey key;
      key.event = event;
      key.chip = rawChip;
      key.position = position;
      std::array<double, 7> &slot = padTimes[key];
      if (maxA >= slot[0])
      {
        slot[0] = maxA;
        slot[1] = arrivalT > -900.0 ? arrivalT : -999.0;
        slot[2] = arrivalTe > -900.0 ? arrivalTe : -999.0;
        for (int f = 0; f < 4; ++f)
          slot[3 + f] = arrivalTSigmoid[f] > -900.0
                            ? arrivalTSigmoid[f]
                            : -999.0;
      }
    }

    if (position < 0 || position >= 121)
    {
      ++unmappedRecords;
      continue;
    }

    const double baseline = ped;
    const double rms = std::max(1.0, std::abs(pedrms));

    double peak = 0.0;
    bool overThreshold = false;
    const double charge =
        PadCharge(*wave, baseline, rms, peak, overThreshold);
    if (!overThreshold || charge <= 0.0)
      continue;

    ++selectedRecords;
    PadEvent &padEvent = events[event][rawChip - 12];
    if (padEvent.timestamp < 0 ||
        (timestamp >= 0 && timestamp < padEvent.timestamp))
      padEvent.timestamp = timestamp;

    PadHit &hit = padEvent.hits[position];
    // Cluster timing counts each physical pad once. As in the legacy
    // window timing, duplicate records use the largest-amplitude waveform.
    if (maxA >= hit.timingAmplitude)
    {
      hit.timingAmplitude = maxA;
      for (int f = 0; f < 4; ++f)
        hit.timeSigmoid[f] = arrivalTSigmoid[f];
    }
    padEvent.arrivalT.Add(arrivalT);
    padEvent.arrivalTa.Add(arrivalTa);
    padEvent.arrivalTe.Add(arrivalTe);
    for (int f = 0; f < 4; ++f)
      padEvent.arrivalTSigmoid[f].Add(arrivalTSigmoid[f]);
    hit.position = position;
    hit.charge += charge;
    hit.peak = std::max(hit.peak, peak);
  }

  // Per-tracker strip waveforms indexed by event for the track-matched
  // tracker timing: tree index = layer*2 + direction (0=X, 1=Y).
  std::map<int, std::vector<StripWave>> stripsByEvent[2 * MMdetN];
  for (int t = 0; t < 2 * MMdetN; ++t)
  {
    TString name = Form("Track%s%d", (t % 2 == 0) ? "X" : "Y", t / 2);
    TTree *stripTree = dynamic_cast<TTree *>(input.Get(name));
    if (!stripTree)
    {
      std::cerr << " WARNING! Missing step2 tree " << name
                << "; tracker track-matched times stay invalid for it"
                << std::endl;
      continue;
    }
    int stripEvent = -1;
    int stripChannel = -1;
    double stripA = -999.0;
    double stripE = -999.0;
    double stripS[4] = {-999.0, -999.0, -999.0, -999.0};
    int pulseFitStatus = -1;
    int sigmoidFitStatus = -1;
    stripTree->SetBranchAddress("event", &stripEvent);
    stripTree->SetBranchAddress("channel", &stripChannel);
    stripTree->SetBranchAddress("arrivalT", &stripA);
    if (stripTree->GetBranch("arrivalT_e"))
      stripTree->SetBranchAddress("arrivalT_e", &stripE);
    if (stripTree->GetBranch("arrivalT_sigmoid"))
      stripTree->SetBranchAddress("arrivalT_sigmoid", &stripS[0]);
    for (int f = 1; f < 4; ++f)
      if (stripTree->GetBranch(sigmoidBranches[f - 1]))
        stripTree->SetBranchAddress(sigmoidBranches[f - 1], &stripS[f]);
    stripTree->SetBranchAddress("pulseFitStatus", &pulseFitStatus);
    if (stripTree->GetBranch("sigmoidFitStatus"))
      stripTree->SetBranchAddress("sigmoidFitStatus", &sigmoidFitStatus);

    std::map<int, std::vector<StripWave>> &strips = stripsByEvent[t];
    const Long64_t stripEntries = stripTree->GetEntries();
    for (Long64_t k = 0; k < stripEntries; ++k)
    {
      stripTree->GetEntry(k);
      if (stripChannel < 0 || stripChannel >= 256)
        continue;
      StripWave w;
      w.channel = stripChannel;
      if (pulseFitStatus == 0)
      {
        if (std::isfinite(stripA) && stripA > -900.0)
          w.timeA = stripA;
        if (std::isfinite(stripE) && stripE > -900.0)
          w.timeE = stripE;
      }
      if (sigmoidFitStatus == 0 && std::isfinite(stripS[0]) &&
          stripS[0] > -900.0)
        for (int f = 0; f < 4; ++f)
          if (std::isfinite(stripS[f]) && stripS[f] > -900.0)
            w.timeS[f] = stripS[f];
      strips[stripEvent].push_back(w);
    }
  }

  TFile output(outputName, "RECREATE");
  if (output.IsZombie())
  {
    std::cerr << " ERROR! Cannot create pad output " << outputName << std::endl;
    return false;
  }

  TTree bm1Tree("BM1Events", "BM1 adjacent-pad clusters");
  TTree bm2Tree("BM2Events", "BM2 adjacent-pad clusters");
  TTree *padTrees[2] = {&bm1Tree, &bm2Tree};
  int outputEvent[2] = {-1, -1};
  PadPlaneOutput planes[2];
  for (int c = 0; c < 2; ++c)
  {
    padTrees[c]->Branch("event", &outputEvent[c], "event/I");
    planes[c].Branch(*padTrees[c], "");
  }

  Long64_t outputClusters[2] = {0, 0};
  std::map<int, std::array<PadPlaneOutput, 2>> clusterEvents;
  for (const auto &eventItem : events)
  {
    for (int c = 0; c < 2; ++c)
    {
      if (eventItem.second[c].hits.empty())
        continue;
      outputEvent[c] = eventItem.first;
      planes[c].Set(eventItem.second[c]);
      clusterEvents[eventItem.first][c] = planes[c];
      outputClusters[c] += planes[c].nclusters;
      padTrees[c]->Fill();
    }
  }

  // Track-matched event arrival times: BM pads inside the 3x3 window around
  // the track projection and tracker strips in clusters (gap <= 1) within
  // +/-2 strips of the predicted hit strip. Time definitions: a, e, sigmoid;
  // per definition the mean and min are written. The event-level value is
  // the mean/min over ALL selected strips (X and Y pooled).
  {
    TTree trackTimeTree("TrackArrivalEvents",
                        "Track-matched arrival times per event");

    int trEvent = -1;
    Long64_t trTimestamp = -1;
    bool trUseflag = false;
    double trXfit[3] = {0.0, 0.0, 0.0};
    double trYfit[3] = {0.0, 0.0, 0.0};
    double trPredStripX[MMdetN] = {0.0};
    double trPredStripY[MMdetN] = {0.0};
    trackTree->SetBranchAddress("event", &trEvent);
    trackTree->SetBranchAddress("timestamp", &trTimestamp);
    trackTree->SetBranchAddress("useflag", &trUseflag);
    trackTree->SetBranchAddress("XfitPars", trXfit);
    trackTree->SetBranchAddress("YfitPars", trYfit);
    if (trackTree->GetBranch("predStripX"))
      trackTree->SetBranchAddress("predStripX", trPredStripX);
    if (trackTree->GetBranch("predStripY"))
      trackTree->SetBranchAddress("predStripY", trPredStripY);

    double bmA_mean[2], bmA_min[2], bmE_mean[2], bmE_min[2];
    double bmS_mean[4][2], bmS_min[4][2];
    int bmA_n[2], bmE_n[2], bmS_n[4][2];
    double bm_pred_x[2], bm_pred_y[2];
    int bm_pred_ix[2], bm_pred_iy[2];
    double bmClusterS_mean[4][2], bmClusterS_min[4][2];
    int bmClusterS_n[4][2], bmMatchedClusterCount[2], bmMatchedClusterIndex[2];
    int bmMatchedClusterSize[2];
    double bmMatchedClusterDx[2], bmMatchedClusterDy[2];
    double trA_mean[MMdetN][2], trA_min[MMdetN][2];
    double trE_mean[MMdetN][2], trE_min[MMdetN][2];
    double trS_mean[4][MMdetN][2], trS_min[4][MMdetN][2];
    int trA_n[MMdetN][2], trE_n[MMdetN][2], trS_n[4][MMdetN][2];
    // Detector-level timing: pool every selected X/Y strip once. In
    // particular physical Tracker2 is index 1, NOT the all-layer ev* time.
    double trXY_mean[6][MMdetN], trXY_min[6][MMdetN];
    int trXY_n[6][MMdetN];
    double evA_mean = -999.0, evA_min = -999.0;
    double evE_mean = -999.0, evE_min = -999.0;
    double evS_mean[4] = {-999.0, -999.0, -999.0, -999.0};
    double evS_min[4] = {-999.0, -999.0, -999.0, -999.0};
    int evA_n = 0, evE_n = 0, evS_n[4] = {0, 0, 0, 0};

    trackTimeTree.Branch("event", &trEvent, "event/I");
    trackTimeTree.Branch("timestamp", &trTimestamp, "timestamp/L");
    trackTimeTree.Branch("bmA_mean", bmA_mean, "bmA_mean[2]/D");
    trackTimeTree.Branch("bmA_min", bmA_min, "bmA_min[2]/D");
    trackTimeTree.Branch("bmA_n", bmA_n, "bmA_n[2]/I");
    trackTimeTree.Branch("bmE_mean", bmE_mean, "bmE_mean[2]/D");
    trackTimeTree.Branch("bmE_min", bmE_min, "bmE_min[2]/D");
    trackTimeTree.Branch("bmE_n", bmE_n, "bmE_n[2]/I");
    const char *fractionTags[4] = {"", "30", "40", "50"};
    for (int f = 0; f < 4; ++f)
    {
      const TString bmName = Form("bmS%s", fractionTags[f]);
      trackTimeTree.Branch(bmName + "_mean", bmS_mean[f],
                           Form("%s_mean[2]/D", bmName.Data()));
      trackTimeTree.Branch(bmName + "_min", bmS_min[f],
                           Form("%s_min[2]/D", bmName.Data()));
      trackTimeTree.Branch(bmName + "_n", bmS_n[f],
                           Form("%s_n[2]/I", bmName.Data()));
      const TString clusterName = Form("bmClusterS%s", fractionTags[f]);
      trackTimeTree.Branch(clusterName + "_mean", bmClusterS_mean[f],
                           Form("%s_mean[2]/D", clusterName.Data()));
      trackTimeTree.Branch(clusterName + "_min", bmClusterS_min[f],
                           Form("%s_min[2]/D", clusterName.Data()));
      trackTimeTree.Branch(clusterName + "_n", bmClusterS_n[f],
                           Form("%s_n[2]/I", clusterName.Data()));
    }
    trackTimeTree.Branch("bm_pred_x", bm_pred_x, "bm_pred_x[2]/D");
    trackTimeTree.Branch("bm_pred_y", bm_pred_y, "bm_pred_y[2]/D");
    trackTimeTree.Branch("bm_pred_ix", bm_pred_ix, "bm_pred_ix[2]/I");
    trackTimeTree.Branch("bm_pred_iy", bm_pred_iy, "bm_pred_iy[2]/I");
    trackTimeTree.Branch("bm_matched_cluster_count", bmMatchedClusterCount,
                        "bm_matched_cluster_count[2]/I");
    trackTimeTree.Branch("bm_matched_cluster_index", bmMatchedClusterIndex,
                        "bm_matched_cluster_index[2]/I");
    trackTimeTree.Branch("bm_matched_cluster_size", bmMatchedClusterSize,
                        "bm_matched_cluster_size[2]/I");
    trackTimeTree.Branch("bm_matched_cluster_dx", bmMatchedClusterDx,
                        "bm_matched_cluster_dx[2]/D");
    trackTimeTree.Branch("bm_matched_cluster_dy", bmMatchedClusterDy,
                        "bm_matched_cluster_dy[2]/D");
    trackTimeTree.Branch("trA_mean", trA_mean,
                         Form("trA_mean[%d][2]/D", MMdetN));
    trackTimeTree.Branch("trA_min", trA_min,
                         Form("trA_min[%d][2]/D", MMdetN));
    trackTimeTree.Branch("trA_n", trA_n, Form("trA_n[%d][2]/I", MMdetN));
    trackTimeTree.Branch("trE_mean", trE_mean,
                         Form("trE_mean[%d][2]/D", MMdetN));
    trackTimeTree.Branch("trE_min", trE_min,
                         Form("trE_min[%d][2]/D", MMdetN));
    trackTimeTree.Branch("trE_n", trE_n, Form("trE_n[%d][2]/I", MMdetN));
    for (int f = 0; f < 4; ++f)
    {
      const TString trName = Form("trS%s", fractionTags[f]);
      trackTimeTree.Branch(trName + "_mean", trS_mean[f],
                           Form("%s_mean[%d][2]/D", trName.Data(), MMdetN));
      trackTimeTree.Branch(trName + "_min", trS_min[f],
                           Form("%s_min[%d][2]/D", trName.Data(), MMdetN));
      trackTimeTree.Branch(trName + "_n", trS_n[f],
                           Form("%s_n[%d][2]/I", trName.Data(), MMdetN));
    }
    const char *timeTags[6] = {"A", "E", "S", "S30", "S40", "S50"};
    for (int m = 0; m < 6; ++m)
    {
      TString name = Form("trXY%s_mean", timeTags[m]);
      trackTimeTree.Branch(name, trXY_mean[m], Form("%s[%d]/D", name.Data(), MMdetN));
      name = Form("trXY%s_min", timeTags[m]);
      trackTimeTree.Branch(name, trXY_min[m], Form("%s[%d]/D", name.Data(), MMdetN));
      name = Form("trXY%s_n", timeTags[m]);
      trackTimeTree.Branch(name, trXY_n[m], Form("%s[%d]/I", name.Data(), MMdetN));
    }
    trackTimeTree.Branch("evA_mean", &evA_mean, "evA_mean/D");
    trackTimeTree.Branch("evA_min", &evA_min, "evA_min/D");
    trackTimeTree.Branch("evA_n", &evA_n, "evA_n/I");
    trackTimeTree.Branch("evE_mean", &evE_mean, "evE_mean/D");
    trackTimeTree.Branch("evE_min", &evE_min, "evE_min/D");
    trackTimeTree.Branch("evE_n", &evE_n, "evE_n/I");
    for (int f = 0; f < 4; ++f)
    {
      const TString evName = Form("evS%s", fractionTags[f]);
      trackTimeTree.Branch(evName + "_mean", &evS_mean[f],
                           Form("%s_mean/D", evName.Data()));
      trackTimeTree.Branch(evName + "_min", &evS_min[f],
                           Form("%s_min/D", evName.Data()));
      trackTimeTree.Branch(evName + "_n", &evS_n[f],
                           Form("%s_n/I", evName.Data()));
    }

    Long64_t trackTimeEvents = 0;
    const Long64_t nTrackEntries = trackTree->GetEntries();
    for (Long64_t k = 0; k < nTrackEntries; ++k)
    {
      trackTree->GetEntry(k);
      if (!trUseflag)
        continue;

      TimeAccumulator bmAcc[2][6];
      TimeAccumulator trAcc[MMdetN][2][6];
      TimeAccumulator trXYAcc[MMdetN][6];
      TimeAccumulator evAcc[6];

      for (int c = 0; c < 2; ++c)
      {
        for (int f = 0; f < 4; ++f)
        {
          bmClusterS_mean[f][c] = bmClusterS_min[f][c] = -999.0;
          bmClusterS_n[f][c] = 0;
        }
        bmMatchedClusterCount[c] = bmMatchedClusterSize[c] = 0;
        bmMatchedClusterIndex[c] = -1;
        bmMatchedClusterDx[c] = bmMatchedClusterDy[c] = -999.0;
        const double xPred = trXfit[0] * kBMZ[c] + trXfit[1] +
                             bmAlignDx[c];
        const double yPred = trYfit[0] * kBMZ[c] + trYfit[1] +
                             bmAlignDy[c];
        bm_pred_x[c] = xPred;
        bm_pred_y[c] = yPred;
        const double xLocal = xPred + kPadCentre + kPadDxOff;
        const double yLocal = yPred + kPadCentre + kPadDyOff;
        const int ix = static_cast<int>(std::floor(
            (xLocal - kPadPitch / 2.0) / kPadPitch + 0.5));
        const int iy = static_cast<int>(std::floor(
            (yLocal - kPadPitch / 2.0) / kPadPitch + 0.5));
        bm_pred_ix[c] = ix;
        bm_pred_iy[c] = iy;
        if (ix < 0 || ix >= kPadGrid || iy < 0 || iy >= kPadGrid)
          continue;
        // Associate clusters using the same spatial window as the legacy
        // pad mode. Count spatial candidates before testing timing validity,
        // so a failed time fit cannot hide an ambiguous cluster association.
        const auto clusterEvent = clusterEvents.find(trEvent);
        if (clusterEvent != clusterEvents.end())
        {
          const PadPlaneOutput &plane = clusterEvent->second[c];
          std::set<int> matches;
          for (std::size_t h = 0; h < plane.hitPosition.size(); ++h)
          {
            const int position = plane.hitPosition[h];
            if (std::abs(position % kPadGrid - ix) <= 1 &&
                std::abs(position / kPadGrid - iy) <= 1)
              matches.insert(plane.hitClusterIndex[h]);
          }
          bmMatchedClusterCount[c] = static_cast<int>(matches.size());
          if (matches.size() == 1)
          {
            const int cluster = *matches.begin();
            bmMatchedClusterIndex[c] = cluster;
            bmMatchedClusterSize[c] = plane.clusterSize[cluster];
            bmMatchedClusterDx[c] = plane.clusterX[cluster] - xLocal;
            bmMatchedClusterDy[c] = plane.clusterY[cluster] - yLocal;
            // Use the complete matched cluster, including members outside
            // the matching window; other clusters never enter this time.
            for (int f = 0; f < 4; ++f)
            {
              bmClusterS_mean[f][c] = plane.clusterTimeSigmoidMean[f][cluster];
              bmClusterS_min[f][c] = plane.clusterTimeSigmoidMin[f][cluster];
              bmClusterS_n[f][c] = plane.clusterTimeSigmoidN[f][cluster];
            }
          }
        }
        for (int di = -1; di <= 1; ++di)
        {
          for (int dj = -1; dj <= 1; ++dj)
          {
            const int nx = ix + di;
            const int ny = iy + dj;
            if (nx < 0 || nx >= kPadGrid || ny < 0 || ny >= kPadGrid)
              continue;
            PadTimeKey key;
            key.event = trEvent;
            key.chip = kBMChip[c];
            key.position = nx + kPadGrid * ny;
            std::map<PadTimeKey, std::array<double, 7>>::const_iterator it =
                padTimes.find(key);
            if (it == padTimes.end())
              continue;
            const std::array<double, 7> &slot = it->second;
            bmAcc[c][0].Add(slot[1]);
            bmAcc[c][1].Add(slot[2]);
            for (int f = 0; f < 4; ++f)
              bmAcc[c][2 + f].Add(slot[3 + f]);
          }
        }
      }

      for (int j = 0; j < MMdetN; ++j)
      {
        for (int d = 0; d < 2; ++d)
        {
          const double pred = (d == 0) ? trPredStripX[j] : trPredStripY[j];
          if (!std::isfinite(pred) || pred < 0.0 || pred >= 256.0)
            continue;
          const int t = j * 2 + d;
          std::map<int, std::vector<StripWave>>::const_iterator it =
              stripsByEvent[t].find(trEvent);
          if (it == stripsByEvent[t].end())
            continue;

          std::vector<StripWave> sorted = it->second;
          std::sort(sorted.begin(), sorted.end(),
                    [](const StripWave &a, const StripWave &b) {
                      return a.channel < b.channel;
                    });

          std::vector<StripWave> current;
          const auto flushCluster = [&]() {
            if (current.empty())
              return;
            const int lo = current.front().channel;
            const int hi = current.back().channel;
            const bool selected = hi >= pred - kStripWindow &&
                                  lo <= pred + kStripWindow;
            if (selected)
            {
              for (const StripWave &w : current)
              {
                trAcc[j][d][0].Add(w.timeA);
                trAcc[j][d][1].Add(w.timeE);
                trXYAcc[j][0].Add(w.timeA);
                trXYAcc[j][1].Add(w.timeE);
                evAcc[0].Add(w.timeA);
                evAcc[1].Add(w.timeE);
                for (int f = 0; f < 4; ++f)
                {
                  trAcc[j][d][2 + f].Add(w.timeS[f]);
                  trXYAcc[j][2 + f].Add(w.timeS[f]);
                  evAcc[2 + f].Add(w.timeS[f]);
                }
              }
            }
            current.clear();
          };
          for (const StripWave &w : sorted)
          {
            if (!current.empty() &&
                w.channel - current.back().channel > kClusterGap)
              flushCluster();
            current.push_back(w);
          }
          flushCluster();
        }
      }

      for (int c = 0; c < 2; ++c)
      {
        bmA_mean[c] = bmAcc[c][0].Mean();
        bmA_min[c] = bmAcc[c][0].Minimum();
        bmA_n[c] = bmAcc[c][0].count;
        bmE_mean[c] = bmAcc[c][1].Mean();
        bmE_min[c] = bmAcc[c][1].Minimum();
        bmE_n[c] = bmAcc[c][1].count;
        for (int f = 0; f < 4; ++f)
        {
          bmS_mean[f][c] = bmAcc[c][2 + f].Mean();
          bmS_min[f][c] = bmAcc[c][2 + f].Minimum();
          bmS_n[f][c] = bmAcc[c][2 + f].count;
        }
      }
      for (int j = 0; j < MMdetN; ++j)
      {
        for (int d = 0; d < 2; ++d)
        {
          trA_mean[j][d] = trAcc[j][d][0].Mean();
          trA_min[j][d] = trAcc[j][d][0].Minimum();
          trA_n[j][d] = trAcc[j][d][0].count;
          trE_mean[j][d] = trAcc[j][d][1].Mean();
          trE_min[j][d] = trAcc[j][d][1].Minimum();
          trE_n[j][d] = trAcc[j][d][1].count;
          for (int f = 0; f < 4; ++f)
          {
            trS_mean[f][j][d] = trAcc[j][d][2 + f].Mean();
            trS_min[f][j][d] = trAcc[j][d][2 + f].Minimum();
            trS_n[f][j][d] = trAcc[j][d][2 + f].count;
          }
        }
      }
      evA_mean = evAcc[0].Mean();
      for (int m = 0; m < 6; ++m)
        for (int j = 0; j < MMdetN; ++j)
        {
          trXY_mean[m][j] = trXYAcc[j][m].Mean();
          trXY_min[m][j] = trXYAcc[j][m].Minimum();
          trXY_n[m][j] = trXYAcc[j][m].count;
        }
      evA_min = evAcc[0].Minimum();
      evA_n = evAcc[0].count;
      evE_mean = evAcc[1].Mean();
      evE_min = evAcc[1].Minimum();
      evE_n = evAcc[1].count;
      for (int f = 0; f < 4; ++f)
      {
        evS_mean[f] = evAcc[2 + f].Mean();
        evS_min[f] = evAcc[2 + f].Minimum();
        evS_n[f] = evAcc[2 + f].count;
      }

      trackTimeTree.Fill();
      ++trackTimeEvents;
    }

    output.WriteTObject(&bm1Tree);
    output.WriteTObject(&bm2Tree);
    output.WriteTObject(&trackTimeTree);
    TNamed detectorSelection("BMPadDetectorSelection",
                             Form("BM1 enabled=%d; BM2 enabled=%d; internal chips=12,13",
                                  bm1Enabled, bm2Enabled));
    output.WriteTObject(&detectorSelection);
    output.Close();

    std::cout << " --> Pad reconstruction: records=" << padRecords
              << ", mapped=" << selectedRecords
              << ", BM1 events=" << bm1Tree.GetEntries()
              << ", BM2 events=" << bm2Tree.GetEntries()
              << ", BM1 clusters=" << outputClusters[0]
              << ", BM2 clusters=" << outputClusters[1];
    if (unmappedRecords > 0)
      std::cout << ", unmapped=" << unmappedRecords;
    std::cout << ", output=" << outputName << std::endl;
    std::cout << " --> Track-matched arrival times: events="
              << trackTimeEvents << std::endl;
    return true;
  }
}

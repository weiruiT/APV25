#include "ReconstructTrack.h"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TMatrixD.h>
#include <TNamed.h>
#include <TFitResultPtr.h>
#include <TVectorD.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <tuple>
#include <vector>

#include "TH1D.h"

namespace
{
constexpr double kGeometryHalfWidthMm = 50.0;
constexpr double kResidualAcceptanceMm = 10.0;

struct JointAlignmentEvent
{
  double local[MMdetN][2] = {{0.0}};
};

struct JointAlignmentResult
{
  double dx[MMdetN] = {0.0};
  double dy[MMdetN] = {0.0};
  double rotationZ[MMdetN] = {0.0};
  long long trackCount = 0;
  double objectiveBefore = 0.0;
  double objectiveCoarse = 0.0;
  double objectiveAfter = 0.0;
  int status = -1;
};

struct ResidualStatistics
{
  long long entries = 0;
  double rms = 0.0;
  double quantile16 = 0.0;
  double quantile84 = 0.0;
  double central68Center = 0.0;
  double central68HalfWidth = 0.0;

  int adaptiveStatus = -1;
  double adaptiveMean = 0.0;
  double adaptiveMeanError = 0.0;
  double adaptiveSigma = 0.0;
  double adaptiveSigmaError = 0.0;
  double adaptiveRangeLow = 0.0;
  double adaptiveRangeHigh = 0.0;
  double adaptiveChi2Ndf = 0.0;
  TF1 *adaptiveFunction = nullptr;

  int fixedStatus = -1;
  double fixedMean = 0.0;
  double fixedMeanError = 0.0;
  double fixedSigma = 0.0;
  double fixedSigmaError = 0.0;
  double fixedRangeLow = -2000.0;
  double fixedRangeHigh = 2000.0;
  double fixedChi2Ndf = 0.0;
  TF1 *fixedFunction = nullptr;
};

array<double, 2> TransformTrackerPosition(double localX, double localY,
                                          int detector,
                                          const double dx[MMdetN],
                                          const double dy[MMdetN],
                                          const double rotationZ[MMdetN])
{
  const double cosine = std::cos(rotationZ[detector]);
  const double sine = std::sin(rotationZ[detector]);
  return {localX * cosine - localY * sine + dx[detector],
          localX * sine + localY * cosine + dy[detector]};
}

array<double, 2> TransformTrackerPosition(double localX, double localY,
                                          int detector,
                                          const double pars[MMdetN][6])
{
  const double cosine = std::cos(pars[detector][5]);
  const double sine = std::sin(pars[detector][5]);
  return {localX * cosine - localY * sine + pars[detector][0],
          localX * sine + localY * cosine + pars[detector][1]};
}

array<double, 2> GlobalToTrackerLocal(double globalX, double globalY,
                                      int detector,
                                      const double pars[MMdetN][6])
{
  const double cosine = std::cos(pars[detector][5]);
  const double sine = std::sin(pars[detector][5]);
  const double shiftedX = globalX - pars[detector][0];
  const double shiftedY = globalY - pars[detector][1];
  return {shiftedX * cosine + shiftedY * sine,
          -shiftedX * sine + shiftedY * cosine};
}

double ProjectionChi2(const double position[MMdetN])
{
  double meanZ = 0.0;
  double meanPosition = 0.0;
  for (int detector = 0; detector < MMdetN; ++detector)
  {
    meanZ += Zpos[detector] / MMdetN;
    meanPosition += position[detector] / MMdetN;
  }
  double szz = 0.0;
  double szp = 0.0;
  for (int detector = 0; detector < MMdetN; ++detector)
  {
    szz += (Zpos[detector] - meanZ) * (Zpos[detector] - meanZ);
    szp += (Zpos[detector] - meanZ) *
           (position[detector] - meanPosition);
  }
  if (szz <= 0.0)
    return std::numeric_limits<double>::infinity();
  const double slope = szp / szz;
  const double intercept = meanPosition - slope * meanZ;
  double chi2 = 0.0;
  for (int detector = 0; detector < MMdetN; ++detector)
  {
    const double residual = position[detector] -
                            (slope * Zpos[detector] + intercept);
    chi2 += residual * residual;
  }
  return chi2;
}

double JointAlignmentObjective(const vector<JointAlignmentEvent> &events,
                               const JointAlignmentResult &alignment)
{
  if (events.empty())
    return 1e30;
  double sum = 0.0;
  for (const JointAlignmentEvent &event : events)
  {
    double x[MMdetN] = {0.0};
    double y[MMdetN] = {0.0};
    for (int detector = 0; detector < MMdetN; ++detector)
    {
      const array<double, 2> global = TransformTrackerPosition(
          event.local[detector][0], event.local[detector][1], detector,
          alignment.dx, alignment.dy, alignment.rotationZ);
      x[detector] = global[0];
      y[detector] = global[1];
    }
    sum += ProjectionChi2(x) + ProjectionChi2(y);
  }
  return sum / events.size();
}

int FineJointAlignment(const vector<JointAlignmentEvent> &events,
                       JointAlignmentResult &alignment)
{
  if (events.empty())
    return -1;
  unique_ptr<ROOT::Math::Minimizer> minimizer(
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  if (!minimizer)
    return -2;

  auto objective = [&events](const double *parameters) {
    JointAlignmentResult trial;
    for (int detector = 1; detector < MMdetN; ++detector)
    {
      const int offset = (detector - 1) * 3;
      trial.dx[detector] = parameters[offset];
      trial.dy[detector] = parameters[offset + 1];
      trial.rotationZ[detector] = parameters[offset + 2];
    }
    return JointAlignmentObjective(events, trial);
  };

  ROOT::Math::Functor function(objective, (MMdetN - 1) * 3);
  minimizer->SetFunction(function);
  minimizer->SetPrintLevel(0);
  minimizer->SetMaxFunctionCalls(200000);
  minimizer->SetMaxIterations(20000);
  minimizer->SetTolerance(0.00005);
  for (int detector = 1; detector < MMdetN; ++detector)
  {
    const int offset = (detector - 1) * 3;
    minimizer->SetVariable(offset, Form("dx_tracker%d", detector + 1),
                           alignment.dx[detector], 0.001);
    minimizer->SetVariable(offset + 1,
                           Form("dy_tracker%d", detector + 1),
                           alignment.dy[detector], 0.001);
    minimizer->SetVariable(offset + 2,
                           Form("rotationZ_tracker%d", detector + 1),
                           alignment.rotationZ[detector], 0.0001);
  }
  minimizer->Minimize();
  for (int detector = 1; detector < MMdetN; ++detector)
  {
    const int offset = (detector - 1) * 3;
    alignment.dx[detector] = minimizer->X()[offset];
    alignment.dy[detector] = minimizer->X()[offset + 1];
    alignment.rotationZ[detector] = minimizer->X()[offset + 2];
  }
  return minimizer->Status();
}

ResidualStatistics MeasureResidual(TH1 *histogram, const string &label)
{
  ResidualStatistics result;
  if (!histogram || histogram->GetEntries() <= 0)
    return result;

  result.entries = static_cast<long long>(histogram->GetEntries());
  result.rms = histogram->GetRMS();
  double probabilities[2] = {0.16, 0.84};
  double quantiles[2] = {0.0, 0.0};
  histogram->GetQuantiles(2, quantiles, probabilities);
  result.quantile16 = quantiles[0];
  result.quantile84 = quantiles[1];
  result.central68Center = 0.5 * (quantiles[0] + quantiles[1]);
  result.central68HalfWidth = 0.5 * (quantiles[1] - quantiles[0]);

  const double peak = histogram->GetBinCenter(histogram->GetMaximumBin());
  const double initialSigma =
      std::max(50.0, std::min(2000.0, histogram->GetRMS()));
  const double initialHalfRange = std::max(300.0, 2.5 * initialSigma);
  double fitLow = std::max(histogram->GetXaxis()->GetXmin(),
                           peak - initialHalfRange);
  double fitHigh = std::min(histogram->GetXaxis()->GetXmax(),
                            peak + initialHalfRange);
  result.adaptiveFunction = new TF1((label + "_adaptive").c_str(), "gaus",
                                    fitLow, fitHigh);
  result.adaptiveFunction->SetParameters(histogram->GetMaximum(), peak,
                                         initialSigma);
  TFitResultPtr adaptiveFit =
      histogram->Fit(result.adaptiveFunction, "Q0RS");
  result.adaptiveStatus = static_cast<int>(adaptiveFit);
  double fittedSigma = std::abs(result.adaptiveFunction->GetParameter(2));
  if (std::isfinite(fittedSigma) && fittedSigma > 5.0)
  {
    const double fittedMean = result.adaptiveFunction->GetParameter(1);
    fitLow = std::max(histogram->GetXaxis()->GetXmin(),
                      fittedMean - 2.5 * fittedSigma);
    fitHigh = std::min(histogram->GetXaxis()->GetXmax(),
                       fittedMean + 2.5 * fittedSigma);
    result.adaptiveFunction->SetRange(fitLow, fitHigh);
    adaptiveFit = histogram->Fit(result.adaptiveFunction, "Q0RS+");
    result.adaptiveStatus = static_cast<int>(adaptiveFit);
  }
  result.adaptiveRangeLow = fitLow;
  result.adaptiveRangeHigh = fitHigh;
  result.adaptiveMean = result.adaptiveFunction->GetParameter(1);
  result.adaptiveMeanError = result.adaptiveFunction->GetParError(1);
  result.adaptiveSigma =
      std::abs(result.adaptiveFunction->GetParameter(2));
  result.adaptiveSigmaError = result.adaptiveFunction->GetParError(2);
  if (result.adaptiveFunction->GetNDF() > 0)
    result.adaptiveChi2Ndf = result.adaptiveFunction->GetChisquare() /
                             result.adaptiveFunction->GetNDF();
  result.adaptiveFunction->SetLineColor(kRed);
  result.adaptiveFunction->SetLineWidth(2);

  result.fixedFunction = new TF1((label + "_fixed2mm").c_str(), "gaus",
                                 result.fixedRangeLow,
                                 result.fixedRangeHigh);
  result.fixedFunction->SetParameters(histogram->GetMaximum(), peak,
                                      initialSigma);
  TFitResultPtr fixedFit = histogram->Fit(result.fixedFunction, "Q0RS+");
  result.fixedStatus = static_cast<int>(fixedFit);
  result.fixedMean = result.fixedFunction->GetParameter(1);
  result.fixedMeanError = result.fixedFunction->GetParError(1);
  result.fixedSigma = std::abs(result.fixedFunction->GetParameter(2));
  result.fixedSigmaError = result.fixedFunction->GetParError(2);
  if (result.fixedFunction->GetNDF() > 0)
    result.fixedChi2Ndf = result.fixedFunction->GetChisquare() /
                          result.fixedFunction->GetNDF();
  result.fixedFunction->SetLineColor(kBlue + 1);
  result.fixedFunction->SetLineStyle(2);
  result.fixedFunction->SetLineWidth(2);
  return result;
}

bool HasFitData(TH1 *hist, double fitLow, double fitHigh)
{
  if (!hist || hist->GetEntries() <= 0 || hist->GetNbinsX() <= 0)
    return false;

  int firstBin = hist->GetXaxis()->FindBin(fitLow);
  int lastBin = hist->GetXaxis()->FindBin(fitHigh);
  if (lastBin < firstBin)
  {
    const int tmp = firstBin;
    firstBin = lastBin;
    lastBin = tmp;
  }
  if (firstBin < 1)
    firstBin = 1;
  if (lastBin > hist->GetNbinsX())
    lastBin = hist->GetNbinsX();

  return hist->Integral(firstBin, lastBin) > 0;
}

TF1 *FitGaussianSafe(TH1 *hist, double fitLow, double fitHigh, const char *label)
{
  if (!HasFitData(hist, fitLow, fitHigh))
  {
    cout << " [WARNING!] Skip gaussian fit for " << label << ": no data in fit range." << endl;
    return nullptr;
  }

  const int fitStatus = hist->Fit("gaus", "QR", "", fitLow, fitHigh);
  TF1 *fit = hist->GetFunction("gaus");
  if (fitStatus != 0 || !fit)
  {
    cout << " [WARNING!] Gaussian fit failed for " << label << ", status=" << fitStatus << endl;
    return nullptr;
  }

  return fit;
}

double GaussianParameterSafe(TH1 *hist, int par, double fitLow, double fitHigh, double fallback, const char *label)
{
  TF1 *fit = FitGaussianSafe(hist, fitLow, fitHigh, label);
  return fit ? fit->GetParameter(par) : fallback;
}

double DivideSafe(double numerator, double denominator)
{
  return denominator != 0 ? numerator / denominator : 0;
}
}

using namespace TMath;

ReconstructTrack::ReconstructTrack(TString fdecName, TString frecName, int f_force)
{
  decName = fdecName;
  recName = frecName;
  force = f_force;
  alignmentPar3Name = MMalignmentPar3Name;
  alignmentPar6Name = MMalignmentPar6Name;
  define_const_matrix();
  define_origin_matrix();
}

ReconstructTrack::~ReconstructTrack()
{
}

bool ReconstructTrack::Init()
{
  fdecFile = new TFile(decName, "read");
  if (!fdecFile->IsOpen())
  {
    cout << " ERROR! " << decName << " cant open!" << endl;
    return false;
  }
  cout << " --> Reconstruct track from " << decName << endl;
  for (int i = 0; i < MMdetN; i++)
  {
    int q = i;
    fdecTree[i] = (TTree *)fdecFile->Get(Form("Track%d", q));
    // fdecTree[i] = (TTree *)fdecFile->Get(Form("Track%d", i));
    if (fdecTree[i] == 0)
    {
      cout << " ERROR! Track" << i << " Tree not exist in " << decName << endl;
      return false;
    }
    fdecTree[i]->SetBranchAddress("event", &event, &b_event);
    fdecTree[i]->SetBranchAddress("sig", &sig, &b_sig);
    fdecTree[i]->SetBranchAddress("sig_x", &sig_x, &b_sig_x);
    fdecTree[i]->SetBranchAddress("sig_y", &sig_y, &b_sig_y);
    fdecTree[i]->SetBranchAddress("x", &x, &b_x);
    fdecTree[i]->SetBranchAddress("y", &y, &b_y);
    fdecTree[i]->SetBranchAddress("z", &z, &b_z);
    fdecTree[i]->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
    fdecTree[i]->SetBranchAddress("hit_amp_x", &hit_amp_x, &b_hit_amp_x);
    fdecTree[i]->SetBranchAddress("hit_amp_y", &hit_amp_y, &b_hit_amp_y);
    fdecTree[i]->SetBranchAddress("hit_strip_num_x", &hit_strip_num_x, &b_hit_strip_num_x);
    fdecTree[i]->SetBranchAddress("hit_strip_num_y", &hit_strip_num_y, &b_hit_strip_num_y);
    fdecTree[i]->SetBranchAddress("x_nhits", &x_nhits, &b_x_nhits);
    fdecTree[i]->SetBranchAddress("y_nhits", &y_nhits, &b_y_nhits);
    if (!fdecTree[i]->GetBranch("cluster_count_x") ||
        !fdecTree[i]->GetBranch("cluster_count_y"))
    {
      cout << " ERROR! Track" << i
           << " is missing cluster_count_x/y required by joint alignment."
           << endl;
      return false;
    }
    fdecTree[i]->SetBranchAddress("cluster_count_x", &cluster_count_x,
                                  &b_cluster_count_x);
    fdecTree[i]->SetBranchAddress("cluster_count_y", &cluster_count_y,
                                  &b_cluster_count_y);
    fdecTree[i]->SetBranchAddress("x_other", &x_other, &b_x_other);
    fdecTree[i]->SetBranchAddress("y_other", &y_other, &b_y_other);
    fdecTree[i]->SetBranchAddress("xMaxT", &xMaxT, &b_xMaxT);
    fdecTree[i]->SetBranchAddress("yMaxT", &yMaxT, &b_yMaxT);
  }

  if (force != 1)
  {
    FileStat_t fStat;
    gSystem->GetPathInfo(recName, fStat);
    if (fStat.fSize != 0)
    {
      cout << " ERROR! " << recName << " is exist!" << endl;
      return false;
    }
  }
  frecFile = new TFile(recName, "recreate");
  if (!frecFile->IsOpen())
  {
    cout << "ERROR! " << recName << " cant open!" << endl;
    return false;
  }
  frecTree = new TTree("Track", "Track");
  frecTree->Branch("event", &event, "event/I");
  frecTree->Branch("MMx", MMx, Form("MMx[%d]/D", MMdetN));
  frecTree->Branch("MMy", MMy, Form("MMy[%d]/D", MMdetN));
  frecTree->Branch("MMz", MMz, Form("MMz[%d]/D", MMdetN));
  frecTree->Branch("MMsigx", MMsigx, Form("MMsigx[%d]/O", MMdetN));
  frecTree->Branch("MMsigy", MMsigy, Form("MMsigy[%d]/O", MMdetN));
  frecTree->Branch("XfitPars", XfitPars, "XfitPars[3]/D");
  frecTree->Branch("YfitPars", YfitPars, "YfitPars[3]/D");
  frecTree->Branch("X01fitPars", X01fitPars, "X01fitPars[3]/D");
  frecTree->Branch("Y01fitPars", Y01fitPars, "Y01fitPars[3]/D");
  frecTree->Branch("X2345fitPars", X2345fitPars, "X2345fitPars[3]/D");
  frecTree->Branch("Y2345fitPars", Y2345fitPars, "Y2345fitPars[3]/D");
  frecTree->Branch("predStripX", predStripX,
                   Form("predStripX[%d]/D", MMdetN));
  frecTree->Branch("predStripY", predStripY,
                   Form("predStripY[%d]/D", MMdetN));
  frecTree->Branch("useflag", &useflag, "useflag/O");
  frecTree->Branch("richuseflag", &richuseflag, "richuseflag/O");
  frecTree->Branch("richextraHit", richextraHit, "richextraHit[3]/D");
  frecTree->Branch("richtrackVec", richtrackVec, "richtrackVec[3]/D");
  frecTree->Branch("btofuseflag", &btofuseflag, "btofuseflag/O");
  frecTree->Branch("btofextraHit", btofextraHit, "btofextraHit[3]/D");
  frecTree->Branch("btoftrackVec", btoftrackVec, "btoftrackVec[3]/D");
  frecTree->Branch("ecaluseflag", &ecaluseflag, "ecaluseflag/O");
  frecTree->Branch("ecalextraHit", ecalextraHit, "ecalextraHit[3]/D");
  frecTree->Branch("ecaltrackVec", ecaltrackVec, "ecaltrackVec[3]/D");
  frecTree->Branch("timestamp", &timestamp, "timestamp/L");
  frecTree->Branch("twotrackangle", &twotrackangle, "twotrackangle/D");
  // get_APV2det_lists(asic2detName.Data());

  alignmentPar3Name = GetPath(recName) + "/MMalignmentPar3.dat";
  alignmentPar6Name = GetPath(recName) + "/MMalignmentPar6.dat";
  if (!JointAlignment(alignmentPar3Name.Data(), alignmentPar6Name.Data()))
  {
    cerr << " [WARNING] Joint tracker alignment failed; attempting to use "
            "existing run-local alignment files."
         << endl;
    ReadAlignmentPars(alignmentPar3Name.Data(), 3);
    ReadAlignmentPars(alignmentPar6Name.Data(), 6);
  }
  // BM pad plane alignment is separate from the tracker alignment: the
  // track parameters are taken from the (existing) tracker alignment and the
  // BM planes are aligned against tracker-projected tracks.
  {
    const TString bmAlignmentName =
        GetPath(recName) + "/MMalignmentParBM.dat";
    std::ifstream test(bmAlignmentName.Data());
    if (test.is_open())
    {
      test.close();
      ReadBMAignmentPars(bmAlignmentName);
    }
    else if (BuildBMAggregation())
    {
      if (!AlignBMPads(bmAlignmentName))
        cerr << " [WARNING] BM alignment failed; BM planes stay unaligned."
             << endl;
    }
    else
      cerr << " [WARNING] BM aggregation failed; BM planes stay unaligned."
           << endl;
  }
  return true;
}

namespace
{
const int kPadGridSize = 11;
const double kPadPitchMm = 9.0;
const double kPadCentreMm = kPadGridSize / 2.0 * kPadPitchMm; // 49.5
const double kBMChargeThresholdRms = 5.0;

// Half-wave charge integral along the dominant polarity, identical to the
// PadReconstruction PadCharge logic (single waveform, no threshold gate).
double BMPadCharge(const std::vector<short> &wave, double baseline,
                   double rms, double &peak)
{
  peak = 0.0;
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
  if (peak <= kBMChargeThresholdRms * std::max(rms, 1.0))
    return 0.0;
  return negativePeak > positivePeak ? negativeCharge : positiveCharge;
}

} // namespace

bool ReconstructTrack::BuildBMAggregation()
{
  bmClusters.clear();
  TString anaName = decName;
  anaName.ReplaceAll("step3-dec", "step2-ana");
  anaName.ReplaceAll("step3_", "step2_");
  TFile *anaFile = new TFile(anaName.Data(), "READ");
  if (!anaFile->IsOpen())
  {
    cerr << " [WARNING] Cannot open step2 file for BM aggregation: "
         << anaName << endl;
    delete anaFile;
    return false;
  }
  TTree *bmTree = dynamic_cast<TTree *>(anaFile->Get("BMPadSignals"));
  if (!bmTree)
  {
    cerr << " [WARNING] BMPadSignals missing in " << anaName << endl;
    anaFile->Close();
    delete anaFile;
    return false;
  }

  int bmEvent = -1;
  int rawChip = -1;
  int padPosition = -1;
  double maxA = 0.0;
  double ped = 0.0;
  double pedrms = 1.0;
  std::vector<short> *wave = nullptr;
  bmTree->SetBranchAddress("event", &bmEvent);
  bmTree->SetBranchAddress("raw_chip", &rawChip);
  bmTree->SetBranchAddress("position", &padPosition);
  bmTree->SetBranchAddress("maxA", &maxA);
  bmTree->SetBranchAddress("ped", &ped);
  bmTree->SetBranchAddress("pedrms", &pedrms);
  bmTree->SetBranchAddress("wave", &wave);

  // Aggregate to a single best record per (event, chip, position)
  struct BestPad
  {
    double maxA = -1.0;
    double ped = 0.0;
    double pedrms = 1.0;
    std::vector<short> wave;
  };
  std::map<std::tuple<int, int, int>, BestPad> best;
  const Long64_t padEntries = bmTree->GetEntries();
  for (Long64_t k = 0; k < padEntries; ++k)
  {
    bmTree->GetEntry(k);
    if (!wave || padPosition < 0 || padPosition >= 121 ||
        (rawChip != 12 && rawChip != 13))
      continue;
    auto &slot = best[std::make_tuple(bmEvent, rawChip, padPosition)];
    if (maxA > slot.maxA)
    {
      slot.maxA = maxA;
      slot.ped = ped;
      slot.pedrms = pedrms;
      slot.wave = *wave;
    }
  }

  // Cluster pads per (event, chip) with 8-connectivity (diagonal counts).
  // Grouped by (event, chip) first so chip12/13 (different detectors) never
  // mix; only the pad positions with the best maxA entry participate.
  std::map<std::pair<int, int>, std::vector<int>> grouped; // (ev,chip) -> positions
  for (const auto &item : best)
    grouped[{std::get<0>(item.first), std::get<1>(item.first)}]
        .push_back(std::get<2>(item.first));

  Long64_t usedClusters = 0;
  Long64_t droppedClusters = 0;
  for (const auto &group : grouped)
  {
    const int ev = group.first.first;
    const int chip = group.first.second;
    const int chipIndex = (chip == 12) ? 0 : 1;

    bool onGrid[kPadGridSize][kPadGridSize] = {{false}};
    for (int pos : group.second)
      onGrid[pos % kPadGridSize][pos / kPadGridSize] = true;

    bool visited[kPadGridSize][kPadGridSize] = {{false}};
    for (int ix = 0; ix < kPadGridSize; ++ix)
    {
      for (int iy = 0; iy < kPadGridSize; ++iy)
      {
        if (!onGrid[ix][iy] || visited[ix][iy])
          continue;
        std::vector<std::pair<int, int>> queue = {{ix, iy}};
        visited[ix][iy] = true;
        std::vector<int> clusterPositions;
        while (!queue.empty())
        {
          const int cx = queue.back().first;
          const int cy = queue.back().second;
          queue.pop_back();
          clusterPositions.push_back(cx + kPadGridSize * cy);
          for (int dx = -1; dx <= 1; ++dx)
          {
            for (int dy = -1; dy <= 1; ++dy)
            {
              if (dx == 0 && dy == 0)
                continue;
              const int nx = cx + dx;
              const int ny = cy + dy;
              if (nx < 0 || nx >= kPadGridSize || ny < 0 ||
                  ny >= kPadGridSize)
                continue;
              if (onGrid[nx][ny] && !visited[nx][ny])
              {
                visited[nx][ny] = true;
                queue.push_back({nx, ny});
              }
            }
          }
        }

        // Physical constraint: a single track footprint on a 9 mm pitch pad
        // plane covers at most a 2x2 pad share (4 pads). Larger clusters are
        // almost certainly noise or multiparticle hit - drop them.
        if (static_cast<int>(clusterPositions.size()) > 4)
        {
          ++droppedClusters;
          continue;
        }

        double chargeSum = 0.0;
        double xSum = 0.0;
        double ySum = 0.0;
        double maxW = 0.0;
        std::vector<double> weights;
        for (int pos : clusterPositions)
        {
          auto it = best.find(std::make_tuple(ev, chip, pos));
          if (it == best.end())
            continue;
          double peak = 0.0;
          const double w = BMPadCharge(it->second.wave, it->second.ped,
                                       it->second.pedrms, peak);
          weights.push_back(w);
          maxW = std::max(maxW, w);
        }
        if (maxW <= 0.0)
          continue;
        std::size_t weightIndex = 0;
        for (int pos : clusterPositions)
        {
          auto it = best.find(std::make_tuple(ev, chip, pos));
          if (it == best.end())
            continue;
          const double w = weights[weightIndex++];
          if (w < 0.25 * maxW)
            continue;
          const int px = pos % kPadGridSize;
          const int py = pos / kPadGridSize;
          const double x = px * kPadPitchMm + kPadPitchMm / 2.0 -
                           kPadCentreMm;
          const double y = py * kPadPitchMm + kPadPitchMm / 2.0 -
                           kPadCentreMm;
          chargeSum += w;
          xSum += w * x;
          ySum += w * y;
        }
        if (chargeSum <= 0.0)
          continue;
        BMCluster cluster;
        cluster.x = xSum / chargeSum;
        cluster.y = ySum / chargeSum;
        cluster.charge = chargeSum;
        cluster.npads = static_cast<int>(clusterPositions.size());
        bmClusters[ev][chipIndex].push_back(cluster);
        ++usedClusters;
      }
    }
  }

  anaFile->Close();
  delete anaFile;
  // Restore the output file as the current directory (ROOT changes
  // gDirectory when a file is opened/closed).
  if (frecFile)
    frecFile->cd();

  cout << " --> BM pad aggregation: clusters kept=" << usedClusters
       << ", dropped (size>4)=" << droppedClusters << endl;
  return !bmClusters.empty();
}

void ReconstructTrack::ReadBMAignmentPars(TString bmAlignmentName)
{
  std::ifstream input(bmAlignmentName.Data());
  if (!input.is_open())
    return;
  std::string line;
  while (std::getline(input, line))
  {
    if (line.empty() || line[0] == '#')
      continue;
    int chip = -1;
    double dx = 0.0;
    double dy = 0.0;
    std::istringstream iss(line);
    if (!(iss >> chip >> dx >> dy))
      continue;
    if (chip == 12)
    {
      bmDx[0] = dx;
      bmDy[0] = dy;
    }
    else if (chip == 13)
    {
      bmDx[1] = dx;
      bmDy[1] = dy;
    }
  }
  bmAligned = true;
  cout << " --> BM alignment: chip12 dx=" << bmDx[0] << " dy=" << bmDy[0]
       << " | chip13 dx=" << bmDx[1] << " dy=" << bmDy[1] << endl;
}

bool ReconstructTrack::AlignBMPads(TString bmAlignmentName)
{
  if (bmClusters.empty())
    return false;
  cout << " --> Begin BM alignment (tracker-projected residuals)" << endl;

  // Track resolution per plane (strip pitch 0.4 mm / sqrt(12)).
  const double sigmaTrack = 0.4 / TMath::Sqrt(12.0);
  const double zBM[2] = {BMZpos[0], BMZpos[1]};

  double dx = 0.0;
  double dy = 0.0;
  double muX = 0.0;
  double muY = 0.0;
  double sigX = 0.0;
  double sigY = 0.0;
  int useCount = 0;

  for (int round = 0; round < 3; ++round)
  {
    TH1D hRx("bmResidualX", "BM residual x;r [mm];events",
             120, -60.0, 60.0);
    TH1D hRy("bmResidualY", "BM residual y;r [mm];events",
             120, -60.0, 60.0);
    useCount = 0;
    int nAllTrackers = 0;
    int nR2Ok = 0;
    int nHasCluster = 0;

    const Long64_t nEntries = fdecTree[0]->GetEntriesFast();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
      std::vector<double> xs, ys, zs;
      int ev = -1;
      for (int j = 0; j < MMdetN; ++j)
      {
        fdecTree[j]->GetEntry(i);
        if (j == 0)
          ev = event;
        if (!(sig_x && sig_y))
          continue;
        if (!(x > 25 && x < 225 && y > 25 && y < 225))
          continue;
        xs.push_back(x * 0.4 - centre);
        ys.push_back(y * 0.4 - centre);
        zs.push_back(Zpos[j]);
      }
      if (static_cast<int>(xs.size()) < 3)
        continue;
      ++nAllTrackers;

      // Own OLS with the same formulas as fit_track (R2 cut applied below).
      const int n = static_cast<int>(xs.size());
      double Kx = 0, Bx = 0, Ky = 0, By = 0;
      double Mx = 0, My = 0, Mz = 0, Mxx = 0, Mxy = 0, Mzz = 0, Mzx = 0,
             Myy = 0, Mzy = 0;
      for (int k = 0; k < n; ++k)
      {
        Mz += zs[k] / n;
        Mx += xs[k] / n;
        My += ys[k] / n;
        Mzz += zs[k] * zs[k] / n;
        Mxx += xs[k] * xs[k] / n;
        Myy += ys[k] * ys[k] / n;
        Mzx += zs[k] * xs[k] / n;
        Mzy += zs[k] * ys[k] / n;
      }
      Kx = (Mzx - Mz * Mx) / (Mzz - Mz * Mz);
      Bx = Mx - Kx * Mz;
      Ky = (Mzy - Mz * My) / (Mzz - Mz * Mz);
      By = My - Ky * Mz;
      const double Rx = (Mzx - Mz * Mx) * (Mzx - Mz * Mx) /
                        (Mxx - Mx * Mx) / (Mzz - Mz * Mz);
      const double Ry = (Mzy - Mz * My) * (Mzy - Mz * My) /
                        (Myy - My * My) / (Mzz - Mz * Mz);
      // BM alignment only needs tracks good enough for extrapolation;
      // the strict 0.99 R2 cut of the main fit is too restrictive here.
      if (Rx < 0.90 || Ry < 0.90)
        continue;
      ++nR2Ok;

      const int chipIndex = 1; // chip13 only; chip12 has no data in this run
      auto clusterIt = bmClusters.find(ev);
      if (clusterIt == bmClusters.end() ||
          clusterIt->second[chipIndex].empty())
        continue;
      ++nHasCluster;

      const double xPred = Kx * zBM[chipIndex] + Bx + dx;
      const double yPred = Ky * zBM[chipIndex] + By + dy;

      double bestDist = 1.0e30;
      const BMCluster *bestCluster = nullptr;
      for (const BMCluster &c : clusterIt->second[chipIndex])
      {
        const double dist = std::hypot(c.x - xPred, c.y - yPred);
        if (dist < bestDist)
        {
          bestDist = dist;
          bestCluster = &c;
        }
      }
      if (!bestCluster || bestDist > 40.0)
        continue;

      hRx.Fill(xPred - bestCluster->x);
      hRy.Fill(yPred - bestCluster->y);
      ++useCount;
    }

    if (round == 0)
      cout << " --> BM alignment selection: all3=" << nAllTrackers
           << " r2ok=" << nR2Ok << " hasCluster=" << nHasCluster
           << " used=" << useCount << endl;

    if (useCount < 20)
    {
      cerr << " [WARNING] BM alignment has too few residual events ("
           << useCount << ")." << endl;
      return false;
    }

    // Core gaussian fit: +-3 sigma trimmed twice; on failure fall back to
    // the histogram mean/RMS (the residual spread is large and sparse).
    TF1 *gausX = nullptr;
    TF1 *gausY = nullptr;
    double loX = -60.0, hiX = 60.0;
    double loY = -60.0, hiY = 60.0;
    for (int it = 0; it < 2; ++it)
    {
      const int stX = hRx.Fit("gaus", "QRN", "", loX, hiX);
      gausX = hRx.GetFunction("gaus");
      const int stY = hRy.Fit("gaus", "QRN", "", loY, hiY);
      gausY = hRy.GetFunction("gaus");
      if (stX != 0 || stY != 0 || !gausX || !gausY)
      {
        gausX = nullptr;
        gausY = nullptr;
        break;
      }
      loX = gausX->GetParameter(1) - 3.0 * std::fabs(gausX->GetParameter(2));
      hiX = gausX->GetParameter(1) + 3.0 * std::fabs(gausX->GetParameter(2));
      loY = gausY->GetParameter(1) - 3.0 * std::fabs(gausY->GetParameter(2));
      hiY = gausY->GetParameter(1) + 3.0 * std::fabs(gausY->GetParameter(2));
    }
    if (!gausX || !gausY)
    {
      cout << " [BM] gaussian core fit failed; using histogram mean/RMS"
           << endl;
      muX = hRx.GetMean();
      muY = hRy.GetMean();
      sigX = hRx.GetRMS();
      sigY = hRy.GetRMS();
    }
    else
    {
      muX = gausX->GetParameter(1);
      muY = gausY->GetParameter(1);
      sigX = std::fabs(gausX->GetParameter(2));
      sigY = std::fabs(gausY->GetParameter(2));
    }

    cout << " --> BM align round " << round << ": events=" << useCount
         << " muX=" << muX << " sigX=" << sigX
         << " | muY=" << muY << " sigY=" << sigY << endl;

    dx -= muX; // constant shift: dx = -mu
    dy -= muY;
    if (std::fabs(muX) < 0.01 && std::fabs(muY) < 0.01)
      break;
  }

  const double weight = (sigmaTrack * sigmaTrack) /
                        (sigX * sigX); // relative to tracker weight 1
  cout << " --> BM alignment: dx=" << dx << " dy=" << dy
       << " | core sigma=" << sigX << " mm (" << sigY << "), weight="
       << weight << endl;
  if (std::fabs(dx) > 40.0 || std::fabs(dy) > 40.0)
    cout << " [CHECK] |dx| or |dy| > 40 mm: check BM coordinate convention."
         << endl;

  ofstream out(bmAlignmentName.Data());
  if (!out.is_open())
  {
    cerr << " [WARNING] Cannot write " << bmAlignmentName << endl;
    return false;
  }
  out << "# BM pad-plane alignment fitted against tracker-projected tracks\n";
  out << "# Format: chip dx dy  (units mm; chip12 = z" << BMZpos[0]
      << ", chip13 = z" << BMZpos[1] << ")\n";
  out << "# residual core sigma x = " << sigX
      << " mm, y = " << sigY << " mm\n";
  out << "12 " << bmDx[0] << " " << bmDy[0] << "\n";
  out << "13 " << dx << " " << dy << "\n";
  out.close();
  bmDx[1] = dx;
  bmDy[1] = dy;
  bmAligned = true;
  cout << " --> BM alignment written to " << bmAlignmentName << endl;
  return true;
}

void ReconstructTrack::Loop()
{
  ofstream op = LogFile(recName);
  TCanvas *cTmp = new TCanvas("rectmp", "rectmp");
  TH1F *prehpos[MMdetN][2];
  TH1F *hpos[MMdetN][2];
  TH1F *hposXY[MMdetN][2];
  TH1F *bhpos[MMdetN][2];
  TH1F *htheta[MMdetN][3];
  TH1F *trackangle[2];
  TH1F *dtrackangle;
  TH2F *hitmap[4];
  TH1F *fitK[4];
  TH1F *reshit[2];
  TH2F *KXY[4];
  for (int i = 0; i < MMdetN; i++)
  {
    prehpos[i][0] = new TH1F(Form("preresX%d", i), Form("preresX%d", i), 500, -4000, 4000);
    prehpos[i][0]->GetXaxis()->SetTitle("x [#mum]");
    prehpos[i][0]->GetYaxis()->SetTitle("counts");
    prehpos[i][1] = new TH1F(Form("preresY%d", i), Form("preresY%d", i), 500, -4000, 4000);
    prehpos[i][1]->GetXaxis()->SetTitle("y [#mum]");
    prehpos[i][1]->GetYaxis()->SetTitle("counts");
    hpos[i][0] = new TH1F(Form("resX%d", i), Form("resX%d", i), 2000, -10000, 10000);
    hpos[i][0]->GetXaxis()->SetTitle("x [#mum]");
    hpos[i][0]->GetYaxis()->SetTitle("counts");
    hpos[i][1] = new TH1F(Form("resY%d", i), Form("resY%d", i), 2000, -10000, 10000);
    hpos[i][1]->GetXaxis()->SetTitle("y [#mum]");
    hpos[i][1]->GetYaxis()->SetTitle("counts");
    hposXY[i][0] = new TH1F(Form("resXY_X%d", i),
                            Form("Tracker%d X residual, three-layer XY coincidence", i + 1),
                            2000, -10000, 10000);
    hposXY[i][0]->GetXaxis()->SetTitle("x [#mum]");
    hposXY[i][0]->GetYaxis()->SetTitle("counts");
    hposXY[i][1] = new TH1F(Form("resXY_Y%d", i),
                            Form("Tracker%d Y residual, three-layer XY coincidence", i + 1),
                            2000, -10000, 10000);
    hposXY[i][1]->GetXaxis()->SetTitle("y [#mum]");
    hposXY[i][1]->GetYaxis()->SetTitle("counts");
    htheta[i][0] = new TH1F(Form("thetaX%d", i), Form("thetaX%d", i), 600, 0, 100);
    htheta[i][0]->GetXaxis()->SetTitle("#theta_{x} [#circ]");
    htheta[i][0]->GetYaxis()->SetTitle("counts");
    htheta[i][1] = new TH1F(Form("thetaY%d", i), Form("thetaY%d", i), 600, 0, 100);
    htheta[i][1]->GetXaxis()->SetTitle("#theta_{y} [#circ]");
    htheta[i][1]->GetYaxis()->SetTitle("counts");
    htheta[i][2] = new TH1F(Form("thetaMM%d", i), Form("thetaMM%d", i), 600, 0, 100);
    htheta[i][2]->GetXaxis()->SetTitle("#theta [#circ]");
    htheta[i][2]->GetYaxis()->SetTitle("counts");
    if (i < 4)
    {
      bhpos[i][0] = new TH1F(Form("resbX%d", i), Form("resbX%d", i), 500, -4000, 4000);
      bhpos[i][0]->GetXaxis()->SetTitle("x [#mum]");
      bhpos[i][0]->GetYaxis()->SetTitle("counts");
      bhpos[i][1] = new TH1F(Form("resbY%d", i), Form("resbY%d", i), 500, -4000, 4000);
      bhpos[i][1]->GetXaxis()->SetTitle("y [#mum]");
      bhpos[i][1]->GetYaxis()->SetTitle("counts");
    }
  }
  for (int i = 0; i < 4; i++)
  {
    hitmap[i] = new TH2F(Form("hitmap%d", i), Form("hitmap%d", i), 800, -100, 100, 800, -100, 100);
    hitmap[i]->GetXaxis()->SetTitle("X [mm]");
    hitmap[i]->GetYaxis()->SetTitle("Y [mm]");
    fitK[i] = new TH1F(Form("k%d", i), Form("k%d", i), 4000, -0.15, 0.15);
    fitK[i]->GetXaxis()->SetTitle("k");
    fitK[i]->GetYaxis()->SetTitle("counts");
    KXY[i] = new TH2F(Form("KXY%d", i), Form("KXY%d", i), 1000, -100, 100, 4000, -0.15, 0.15);
    KXY[i]->GetXaxis()->SetTitle("xy");
    KXY[i]->GetYaxis()->SetTitle("k");
  }
  reshit[0] = new TH1F(Form("reshitX"), Form("reshitX"), 500, -6000, 6000);
  reshit[0]->GetXaxis()->SetTitle("resX [#mum]");
  reshit[0]->GetYaxis()->SetTitle("counts");
  reshit[1] = new TH1F(Form("reshitY"), Form("reshitY"), 500, -6000, 6000);
  reshit[1]->GetXaxis()->SetTitle("resY [#mum]");
  reshit[1]->GetYaxis()->SetTitle("counts");
  trackangle[0] = new TH1F(Form("track01angle"), Form("track01angle"), 1000, 0, 100);
  trackangle[0]->GetXaxis()->SetTitle("mrad");
  trackangle[0]->GetYaxis()->SetTitle("counts");
  trackangle[1] = new TH1F(Form("track2345angle"), Form("track2345angle"), 1000, 0, 100);
  trackangle[1]->GetXaxis()->SetTitle("mrad");
  trackangle[1]->GetYaxis()->SetTitle("counts");
  dtrackangle = new TH1F(Form("dtrackangle"), Form("dtrackangle"), 5000, 50, 50);
  dtrackangle->GetXaxis()->SetTitle("mrad");
  dtrackangle->GetYaxis()->SetTitle("counts");
  double res_sigma1[MMdetN][2] = {0};  // 第一个高斯分量的sigma
  double res_sigma2[MMdetN][2] = {0};  // 第二个高斯分量的sigma
  double bres_sigma1[MMdetN][2] = {0}; // 第一个高斯分量的sigma
  double bres_sigma2[MMdetN][2] = {0}; // 第二个高斯分量的sigma
  double preres[MMdetN][2] = {0};
  double res[MMdetN][2] = {0};
  double resXY[MMdetN][2] = {0};
  ResidualStatistics residualStats[2][MMdetN][2];
  double bres[MMdetN][2] = {0};
  double efficiencyX[MMdetN][3] = {0};
  double efficiencyY[MMdetN][3] = {0};
  double efficiencybX[MMdetN][3] = {0};
  double efficiencybY[MMdetN][3] = {0};
  int count[MMdetN][2] = {0};

  double area1[MMdetN][2] = {0};  // 第一个高斯函数的面积 [det][0:X,1:Y]
  double area2[MMdetN][2] = {0};  // 第二个高斯函数的面积 [det][0:X,1:Y]
  double barea1[MMdetN][2] = {0}; // 第一个高斯函数的面积 [det][0:X,1:Y]
  double barea2[MMdetN][2] = {0}; // 第二个高斯函数的面积 [det][0:X,1:Y]

  int nentries = fdecTree[0]->GetEntriesFast();
  for (int i = 0; i < nentries; i++)
  {
    DrawProcessbar(i, nentries);
    InitData();
    // frecTree->Fill();
    double localX[MMdetN] = {0};
    double localY[MMdetN] = {0};
    int eventClusterCountX[MMdetN] = {0};
    int eventClusterCountY[MMdetN] = {0};
    bool xsigy[MMdetN] = {0};
    bool ysigx[MMdetN] = {0};
    long long ftime = 0;
    for (int j = 0; j < MMdetN; j++)
    {
      fdecTree[j]->GetEntry(i);
      if (j == 0)
        ftime = timestamp;
      eventClusterCountX[j] = cluster_count_x;
      eventClusterCountY[j] = cluster_count_y;
      if (sig_x)
      {
        localX[j] = x * 0.4 - centre;
        MMsigx[j] = true;
      }
      if (sig_y)
      {
        localY[j] = y * 0.4 - centre;
        MMsigy[j] = true;
      }
      MMz[j] = Zpos[j];
    }
    vector<double> fitx, fity, fitxz, fityz;
    vector<double> fitbx, fitby, fitbxz, fitbyz;
    for (int j = 0; j < MMdetN; j++)
    {
      // A Z rotation couples X and Y. Events with both coordinates receive
      // the exact transform. A one-coordinate DUT candidate is corrected
      // later with the predicted orthogonal local coordinate.
      const array<double, 2> aligned = TransformTrackerPosition(
          MMsigx[j] ? localX[j] : 0.0,
          MMsigy[j] ? localY[j] : 0.0, j, algPar6);
      if (MMsigx[j])
        MMx[j] = aligned[0];
      if (MMsigy[j])
        MMy[j] = aligned[1];
      if (j == 1 && MMsigx[1] && MMsigy[1])
        hitmap[0]->Fill(MMx[1], MMy[1]);
    }

    vector<double> measuredX, measuredY, measuredXZ, measuredYZ;
    for (int j = 0; j < MMdetN; ++j)
    {
      if (MMsigx[j])
      {
        measuredX.push_back(MMx[j]);
        measuredXZ.push_back(MMz[j]);
      }
      if (MMsigy[j])
      {
        measuredY.push_back(MMy[j]);
        measuredYZ.push_back(MMz[j]);
      }
    }
    double measuredFitX[3] = {0.0};
    double measuredFitY[3] = {0.0};
    const bool canPredictX = fit_track2(measuredXZ, measuredX, measuredFitX);
    const bool canPredictY = fit_track2(measuredYZ, measuredY, measuredFitY);
    for (int j = 0; j < MMdetN; ++j)
    {
      if (!MMsigx[j] && canPredictX)
      {
        MMx[j] = measuredFitX[0] * MMz[j] + measuredFitX[1];
        xsigy[j] = true;
      }
      if (!MMsigy[j] && canPredictY)
      {
        MMy[j] = measuredFitY[0] * MMz[j] + measuredFitY[1];
        ysigx[j] = true;
      }
    }

    vector<double> fitx01, fity01, fitxz01, fityz01;
    vector<double> fitx23, fity23, fitxz23, fityz23;
    vector<double> fitx45, fity45, fitxz45, fityz45;
    vector<double> fitx2345, fity2345, fitxz2345, fityz2345;
    vector<double> fitallx, fitally, fitallxz, fitallyz;
    double parX01[3], parY01[3], parX2345[3], parY2345[3], parX23[3], parY23[3], parX45[3], parY45[3];

    for (int j = 0; j < MMdetN; j++)
    {
      double parx[3], pary[3];
      double parbx[3], parby[3];
      bool good_fit_x = false;
      bool good_fit_y = false;
      fitx.clear();
      fitxz.clear();
      fity.clear();
      fityz.clear();
      fitbx.clear();
      fitby.clear();
      fitbxz.clear();
      fitbyz.clear();
      if (MMsigx[j])
      {
        if (j < 2)
        {
          fitx01.push_back(MMx[j]);
          fitxz01.push_back(MMz[j]);
        }
        if (j < 6 && j > 1)
        {
          fitx2345.push_back(MMx[j]);
          fitxz2345.push_back(MMz[j]);
        }
        if (j > 1 && j < 4)
        {
          fitx23.push_back(MMx[j]);
          fitxz23.push_back(MMz[j]);
        }
        if (j > 3)
        {
          fitx45.push_back(MMx[j]);
          fitxz45.push_back(MMz[j]);
        }
        fitallx.push_back(MMx[j]);
        fitallxz.push_back(MMz[j]);
      }
      if (MMsigy[j])
      {
        if (j < 2)
        {
          fity01.push_back(MMy[j]);
          fityz01.push_back(MMz[j]);
        }
        if (j < 6 && j > 1)
        {
          fity2345.push_back(MMy[j]);
          fityz2345.push_back(MMz[j]);
        }
        if (j > 1 && j < 4)
        {
          fity23.push_back(MMy[j]);
          fityz23.push_back(MMz[j]);
        }
        if (j > 3)
        {
          fity45.push_back(MMy[j]);
          fityz45.push_back(MMz[j]);
        }
        fitally.push_back(MMy[j]);
        fitallyz.push_back(MMz[j]);
      }
      bool referencePointsAccepted = true;
      for (int k = 0; k < MMdetN; k++)
      {
        if (k == j)
          continue;
        const bool singleTwoDimensionalPoint =
            MMsigx[k] && MMsigy[k] && eventClusterCountX[k] == 1 &&
            eventClusterCountY[k] == 1;
        const bool insideReferenceGeometry =
            singleTwoDimensionalPoint &&
            Abs(MMx[k]) < kGeometryHalfWidthMm &&
            Abs(MMy[k]) < kGeometryHalfWidthMm;
        if (!insideReferenceGeometry)
        {
          referencePointsAccepted = false;
          continue;
        }
        fitx.push_back(MMx[k]);
        fitxz.push_back(MMz[k]);
        fity.push_back(MMy[k]);
        fityz.push_back(MMz[k]);
        if (k > 1 && j > 1)
        {
          fitbx.push_back(MMx[k]);
          fitbxz.push_back(MMz[k]);
          fitby.push_back(MMy[k]);
          fitbyz.push_back(MMz[k]);
        }
      }
      double fitXX = 0;
      double fitYY = 0;
      bool badfitx = false;
      bool badfity = false;
      if (referencePointsAccepted && fit_track(fitxz, fitx, parx))
      {
        good_fit_x = true;
        fitXX = parx[0] * MMz[j] + parx[1];
        for (int kk = 0; kk < 3; kk++)
          XfitPars[kk] = parx[kk];
      }
      if (referencePointsAccepted && fit_track(fityz, fity, pary))
      {
        good_fit_y = true;
        fitYY = pary[0] * MMz[j] + pary[1];
        for (int kk = 0; kk < 3; kk++)
          YfitPars[kk] = pary[kk];
      }

      const bool predictionInsideGeometry =
          good_fit_x && good_fit_y &&
          Abs(fitXX) < kGeometryHalfWidthMm &&
          Abs(fitYY) < kGeometryHalfWidthMm;
      double measuredX = 0.0;
      double measuredY = 0.0;
      bool measuredXValid = false;
      bool measuredYValid = false;
      if (predictionInsideGeometry)
      {
        efficiencyX[j][0]++;
        efficiencyY[j][0]++;
        const array<double, 2> predictedLocal =
            GlobalToTrackerLocal(fitXX, fitYY, j, algPar6);
        if (MMsigx[j])
        {
          const array<double, 2> candidate = TransformTrackerPosition(
              localX[j], MMsigy[j] ? localY[j] : predictedLocal[1], j,
              algPar6);
          measuredX = candidate[0];
          measuredXValid = true;
          if (Abs(measuredX - fitXX) < kResidualAcceptanceMm)
          {
            efficiencyX[j][1]++;
            hpos[j][0]->Fill((measuredX - fitXX) * 1000);
          }
        }
        if (MMsigy[j])
        {
          const array<double, 2> candidate = TransformTrackerPosition(
              MMsigx[j] ? localX[j] : predictedLocal[0], localY[j], j,
              algPar6);
          measuredY = candidate[1];
          measuredYValid = true;
          if (Abs(measuredY - fitYY) < kResidualAcceptanceMm)
          {
            efficiencyY[j][1]++;
            hpos[j][1]->Fill((measuredY - fitYY) * 1000);
          }
        }
        if (j == 1 && measuredXValid &&
            Abs(measuredX - fitXX) * 1000 > 350)
          badfitx = true;
        if (j == 1 && measuredYValid &&
            Abs(measuredY - fitYY) * 1000 > 350)
          badfity = true;
        if (j == 1 && badfity && badfitx)
        {
          hitmap[3]->Fill(measuredX, measuredY);
          hitmap[2]->Fill(fitXX, fitYY);
        }
      }

      // Clean common sample: every Tracker has exactly one X and one Y
      // cluster. Reference and predicted two-dimensional points must stay
      // within the requested +/-50 mm geometry window.
      bool allTrackersHaveOneXY = true;
      for (int k = 0; k < MMdetN; ++k)
        allTrackersHaveOneXY =
            allTrackersHaveOneXY && MMsigx[k] && MMsigy[k] &&
            eventClusterCountX[k] == 1 && eventClusterCountY[k] == 1;
      if (allTrackersHaveOneXY && predictionInsideGeometry &&
          measuredXValid && measuredYValid)
      {
        const double residualX = measuredX - fitXX;
        const double residualY = measuredY - fitYY;
        if (Abs(residualX) < kResidualAcceptanceMm &&
            Abs(residualY) < kResidualAcceptanceMm)
        {
          hposXY[j][0]->Fill(residualX * 1000);
          hposXY[j][1]->Fill(residualY * 1000);
        }
      }

      if (j > 1)
      {
        double bigx = 0;
        double bigy = 0;
        double bigfitx, bigfity;
        if (fit_track3(fitbxz, fitbx, parbx))
        {
          double fitX = parbx[0] * MMz[j] + parbx[1];
          efficiencybX[j - 2][0]++;
          if (MMsigx[j] && Abs(MMx[j] - fitX) < 2.5)
          {
            efficiencybX[j - 2][1]++;
            bhpos[j - 2][0]->Fill((MMx[j] - fitX) * 1000);
            bigx = (MMx[j] - fitX) * 1000;
            bigfitx = fitX;
          }
        }
        if (fit_track3(fitbyz, fitby, parby))
        {
          double fitY = parby[0] * MMz[j] + parby[1];
          efficiencybY[j - 2][0]++;
          if (MMsigy[j] && Abs(MMy[j] - fitY) < 2.5)
          {
            efficiencybY[j - 2][1]++;
            bhpos[j - 2][1]->Fill((MMy[j] - fitY) * 1000);
            bigy = (MMy[j] - fitY) * 1000;
            bigfity = fitY;
          }
        }

        if (j == 1 && bigx != 0 && bigy != 0 && Abs(bigx) > 400)
        {
          hitmap[1]->Fill(bigfitx, bigfity);

          // op << event << endl;
        }
      }
    }

    if (fit_track(fitallxz, fitallx, XfitPars) && fit_track(fitallyz, fitally, YfitPars))
    {
      useflag = 1;
    }
    else
    {
      for (int m = 0; m < 3; m++)
      {
        XfitPars[m] = 0;
        YfitPars[m] = 0;
      }
    }

    // Track-predicted hit position per layer in step2 strip units. The
    // fitted line lives in the aligned global frame; the local strip
    // position is obtained by removing the per-layer alignment offset.
    if (useflag)
    {
      for (int j = 0; j < MMdetN; ++j)
      {
        const array<double, 2> predictedLocal = GlobalToTrackerLocal(
            XfitPars[0] * MMz[j] + XfitPars[1],
            YfitPars[0] * MMz[j] + YfitPars[1], j, algPar6);
        predStripX[j] = (predictedLocal[0] + centre) / 0.4;
        predStripY[j] = (predictedLocal[1] + centre) / 0.4;
      }
    }

    if (fit_track2(fitxz01, fitx01, parX01) && fit_track2(fityz01, fity01, parY01))
    {
      if (useflag)
      {
        fitK[0]->Fill(parX01[0]);
        fitK[1]->Fill(parY01[0]);
      }
      // double ECALZ = 2000;
      // ecaluseflag = 1;
      // ecalextraHit[0] = parX01[0] * ECALZ + parX01[1];
      // ecalextraHit[1] = parY01[0] * ECALZ + parY01[1];
      // ecalextraHit[2] = ECALZ;
      // double modelecal = Sqrt(parX01[0] * parX01[0] + parY01[0] * parY01[0] + 1);
      // ecaltrackVec[0] = parX01[0] / modelecal;
      // ecaltrackVec[1] = parY01[0] / modelecal;
      // ecaltrackVec[2] = 1 / modelecal;

      // double angle01 = ATan(Sqrt(parX01[0] * parX01[0] + parY01[0] * parY01[0])) / Pi() * 180;

      // double angle01 = (TVector3(parX01[0], parY01[0], 1).Angle(TVector3(0, 0, 1))) * 1000;
      // trackangle[0]->Fill(angle01);
      // double angle23 = (TVector3((MMx[2] - MMx[3]) / (MMz[2] - MMz[3]), (MMy[2] - MMy[3]) / (MMz[2] - MMz[3]), 1).Angle(TVector3(0, 0, 1))) * 1000;
      // trackangle[1]->Fill(angle23);
      // double dangle = (TVector3(parX01[0], parY01[0], 1).Angle(TVector3(0, 0, 1))) * 1000;
      // dtrackangle->Fill(dangle);

      // fitK[0]->Fill(parX01[0]);
      // fitK[1]->Fill(parY01[0]);

      for (int m = 0; m < 3; m++)
      {
        X01fitPars[m] = parX01[m];
        Y01fitPars[m] = parY01[m];
      }

      KXY[0]->Fill(MMx[0], parX01[0]);
      KXY[1]->Fill(MMy[0], parY01[0]);
      // reshit[0]->Fill((MMx[0] - MMx[1]) * 1000);
      // reshit[1]->Fill((MMy[0] - MMy[1]) * 1000);

      // double BTOFZ = 3690;
      // btofextraHit[0] = parX01[0] * BTOFZ + parX01[1];
      // btofextraHit[1] = parY01[0] * BTOFZ + parY01[1];
      // btofextraHit[2] = BTOFZ;
      // double modelbtof = Sqrt(parX01[0] * parX01[0] + parY01[0] * parY01[0] + 1);
      // btoftrackVec[0] = parX01[0] / modelbtof;
      // btoftrackVec[1] = parY01[0] / modelbtof;
      // btoftrackVec[2] = 1 / modelbtof;

      richuseflag = true;
      // double RICHZ = 2242;
      // double RICHZ = 8192;
      double RICHZ = 9170;
      richextraHit[0] = parX01[0] * RICHZ + parX01[1];
      richextraHit[1] = parY01[0] * RICHZ + parY01[1];
      richextraHit[2] = RICHZ;
      double modelrich = Sqrt(parX01[0] * parX01[0] + parY01[0] * parY01[0] + 1);
      richtrackVec[0] = parX01[0] / modelrich;
      richtrackVec[1] = parY01[0] / modelrich;
      richtrackVec[2] = 1 / modelrich;
      // if (fit_track3(fitxz2345, fitx2345, parX2345) && fit_track3(fityz2345, fity2345, parY2345))
      // {
      //   double angle01 = (TVector3(parX01[0], parY01[0], 1).Angle(TVector3(0, 0, 1))) * 1000;
      //   trackangle[0]->Fill(angle01);
      //   double angle2345 = (TVector3(parX2345[0], parY2345[0], 1).Angle(TVector3(0, 0, 1))) * 1000;
      //   double dangle = (TVector3(parX01[0], parY01[0], 1).Angle(TVector3(parX2345[0], parY2345[0], 1))) * 1000;
      //   trackangle[1]->Fill(angle2345);
      //   // dtrackangle->Fill(dangle);
      //   fitK[2]->Fill(parX2345[0] - parX01[0]);
      //   fitK[3]->Fill(parY2345[0] - parY01[0]);
      //   if (MMx[2] != 0)
      //     KXY[2]->Fill(MMx[2], parX2345[0]);
      //   if (MMy[2] != 0)
      //     KXY[3]->Fill(MMy[2], parY2345[0]);
      //   // if (Abs(dangle) < 1)
      //   btofuseflag = 1;
      //   for (int m = 0; m < 3; m++)
      //   {
      //     X01fitPars[m] = parX01[m];
      //     Y01fitPars[m] = parY01[m];
      //     X2345fitPars[m] = parX2345[m];
      //     Y2345fitPars[m] = parY2345[m];
      //   }
      //   twotrackangle = dangle;

      //   double ECALZ = 2000;
      //   ecaluseflag = 1;
      //   ecalextraHit[0] = parX2345[0] * ECALZ + parX2345[1];
      //   ecalextraHit[1] = parY2345[0] * ECALZ + parY2345[1];
      //   ecalextraHit[2] = ECALZ;
      //   double modelecal = Sqrt(parX2345[0] * parX2345[0] + parY2345[0] * parY2345[0] + 1);
      //   ecaltrackVec[0] = parX2345[0] / modelecal;
      //   ecaltrackVec[1] = parY2345[0] / modelecal;
      //   ecaltrackVec[2] = 1 / modelecal;
      // }

      // double inbtofx = parX01[0] * 3968 + parX01[1];
      // double inbtofy = parY01[0] * 3968 + parY01[1];
      // fitx2345.push_back(inbtofx);
      // fitxz2345.push_back(3968);
      // fity2345.push_back(inbtofy);
      // fityz2345.push_back(3968);

      if (fit_track3(fitxz2345, fitx2345, parX2345) && fit_track3(fityz2345, fity2345, parY2345))
      {
        // dtrackangle->Fill(dangle);
        if (useflag)
        {
          fitK[2]->Fill(parX01[0] - parX2345[0]);
          fitK[3]->Fill(parY01[0] - parY2345[0]);
        }

        if (MMx[2] != 0)
          KXY[2]->Fill(MMx[2], parX2345[0]);
        if (MMy[2] != 0)
          KXY[3]->Fill(MMy[2], parY2345[0]);
        // if (Abs(dangle) < 1)
        btofuseflag = 1;
        for (int m = 0; m < 3; m++)
        {
          X2345fitPars[m] = parX2345[m];
          Y2345fitPars[m] = parY2345[m];
        }
        reshit[0]->Fill(((parX01[0] * 3980 + parX01[1]) - (parX2345[0] * 3980 + parX2345[1])) * 1000);
        reshit[1]->Fill(((parY01[0] * 3980 + parY01[1]) - (parY2345[0] * 3980 + parY2345[1])) * 1000);

        double ECALZ = 2000;
        ecaluseflag = 1;
        ecalextraHit[0] = parX2345[0] * ECALZ + parX2345[1];
        ecalextraHit[1] = parY2345[0] * ECALZ + parY2345[1];
        ecalextraHit[2] = ECALZ;
        double modelecal = Sqrt(parX2345[0] * parX2345[0] + parY2345[0] * parY2345[0] + 1);
        ecaltrackVec[0] = parX2345[0] / modelecal;
        ecaltrackVec[1] = parY2345[0] / modelecal;
        ecaltrackVec[2] = 1 / modelecal;
      }
    }

    // if (fit_track3(fitxz2345, fitx2345, parX2345) && fit_track3(fityz2345, fity2345, parY2345))
    // {
    //   double angle01 = (TVector3(parX01[0], parY01[0], 1).Angle(TVector3(0, 0, 1))) * 1000;
    //   trackangle[0]->Fill(angle01);
    //   double angle2345 = (TVector3(parX2345[0], parY2345[0], 1).Angle(TVector3(0, 0, 1))) * 1000;
    //   double dangle = (TVector3(parX01[0], parY01[0], 1).Angle(TVector3(parX2345[0], parY2345[0], 1))) * 1000;
    //   trackangle[1]->Fill(angle2345);
    //   dtrackangle->Fill(dangle);
    //   fitK[2]->Fill(parX01[0] - parX2345[0]);
    //   fitK[3]->Fill(parY01[0] - parY2345[0]);
    //   if (MMx[2] != 0)
    //     KXY[2]->Fill(MMx[2], parX2345[0]);
    //   if (MMy[2] != 0)
    //     KXY[3]->Fill(MMy[2], parY2345[0]);
    //   // if (Abs(dangle) < 1)
    //   btofuseflag = 1;
    //   for (int m = 0; m < 3; m++)
    //   {
    //     X01fitPars[m] = parX01[m];
    //     Y01fitPars[m] = parY01[m];
    //     X2345fitPars[m] = parX2345[m];
    //     Y2345fitPars[m] = parY2345[m];
    //   }
    //   twotrackangle = dangle;

    //   double ECALZ = 2000;
    //   ecaluseflag = 1;
    //   ecalextraHit[0] = parX2345[0] * ECALZ + parX2345[1];
    //   ecalextraHit[1] = parY2345[0] * ECALZ + parY2345[1];
    //   ecalextraHit[2] = ECALZ;
    //   double modelecal = Sqrt(parX2345[0] * parX2345[0] + parY2345[0] * parY2345[0] + 1);
    //   ecaltrackVec[0] = parX2345[0] / modelecal;
    //   ecaltrackVec[1] = parY2345[0] / modelecal;
    //   ecaltrackVec[2] = 1 / modelecal;
    // }

    // if (fit_track3(fitxz2345, fitx2345, parX2345) && fit_track3(fityz2345, fity2345, parY2345))
    // {
    //   if (fit_track2(fitxz23, fitx23, parX23) && fit_track2(fityz23, fity23, parY23))
    //   {
    //     if (fit_track2(fitxz45, fitx45, parX45) && fit_track2(fityz45, fity45, parY45))
    //     {
    //       double anglez23 = (TVector3(parX23[0], parY23[0], 1).Angle(TVector3(0, 0, 1))) * 1000;
    //       double anglez45 = (TVector3(parX45[0], parY45[0], 1).Angle(TVector3(0, 0, 1))) * 1000;
    //       double d2345 = (TVector3(parX23[0], parY23[0], 1).Angle(TVector3(parX45[0], parY45[0], 1))) * 1000;
    //       if ((anglez23 - anglez45) < 0)
    //       {
    //         d2345 = -d2345;
    //       }
    //       dtrackangle->Fill(anglez45 - anglez23);
    //     }
    //   }
    // }
    // if (fit_track3(fitxz2345, fitx2345, parX2345) && fit_track3(fityz2345, fity2345, parY2345))
    // {
    //   // dtrackangle->Fill(dangle);
    //   fitK[2]->Fill(parX2345[0] - parX01[0]);
    //   fitK[3]->Fill(parY2345[0] - parY01[0]);
    //   if (MMx[2] != 0)
    //     KXY[2]->Fill(MMx[2], parX2345[0]);
    //   if (MMy[2] != 0)
    //     KXY[3]->Fill(MMy[2], parY2345[0]);
    //   // if (Abs(dangle) < 1)
    //   btofuseflag = 1;
    //   for (int m = 0; m < 3; m++)
    //   {
    //     X2345fitPars[m] = parX2345[m];
    //     Y2345fitPars[m] = parY2345[m];
    //   }

    //   double ECALZ = 2000;
    //   ecaluseflag = 1;
    //   ecalextraHit[0] = parX2345[0] * ECALZ + parX2345[1];
    //   ecalextraHit[1] = parY2345[0] * ECALZ + parY2345[1];
    //   ecalextraHit[2] = ECALZ;
    //   double modelecal = Sqrt(parX2345[0] * parX2345[0] + parY2345[0] * parY2345[0] + 1);
    //   ecaltrackVec[0] = parX2345[0] / modelecal;
    //   ecaltrackVec[1] = parY2345[0] / modelecal;
    //   ecaltrackVec[2] = 1 / modelecal;
    // }

    // if (fit_track3(fitxz2345, fitx2345, parX2345) && fit_track3(fityz2345, fity2345, parY2345))
    // {
    //   // double angle2345 = (TVector3(parX2345[0], parY2345[0], 1).Angle(TVector3(0, 0, 1))) * 1000;
    //   // double dangle = (TVector3(parX01[0], parY01[0], 1).Angle(TVector3(parX2345[0], parY2345[0], 1))) * 1000;
    //   // trackangle[1]->Fill(angle2345);
    //   // dtrackangle->Fill(dangle);
    //   // fitK[2]->Fill(parX2345[0]);
    //   // fitK[3]->Fill(parY2345[0]);
    //   // // if (Abs(dangle) < 1)
    //   // richuseflag = true;
    //   // // double RICHZ = 2242;
    //   // // double RICHZ = 8192;
    //   // double RICHZ = 9170;
    //   // richextraHit[0] = parX2345[0] * RICHZ + parX2345[1];
    //   // richextraHit[1] = parY2345[0] * RICHZ + parY2345[1];
    //   // richextraHit[2] = RICHZ;
    //   // double modelrich = Sqrt(parX2345[0] * parX2345[0] + parY2345[0] * parY2345[0] + 1);
    //   // richtrackVec[0] = parX2345[0] / modelrich;
    //   // richtrackVec[1] = parY2345[0] / modelrich;
    //   // richtrackVec[2] = 1 / modelrich;
    //   // hitmap[2]->Fill(richextraHit[0], richextraHit[1]);
    //   double ECALZ = 9600;
    //   ecaluseflag = 1;
    //   ecalextraHit[0] = parX2345[0] * ECALZ + parX2345[1];
    //   ecalextraHit[1] = parY2345[0] * ECALZ + parY2345[1];
    //   ecalextraHit[2] = ECALZ;
    //   double modelecal = Sqrt(parX2345[0] * parX2345[0] + parY2345[0] * parY2345[0] + 1);
    //   ecaltrackVec[0] = parX2345[0] / modelecal;
    //   ecaltrackVec[1] = parY2345[0] / modelecal;
    //   ecaltrackVec[2] = 1 / modelecal;
    // }
    timestamp = ftime;
    // if (event < 10) cout << "event" << event << "timestamp" << timestamp << endl;

    frecTree->Fill();
  }
  frecTree->Write();

  for (int i = 0; i < MMdetN; i++)
  {
    residualStats[0][i][0] =
        MeasureResidual(hpos[i][0], Form("resX%d", i));
    residualStats[0][i][1] =
        MeasureResidual(hpos[i][1], Form("resY%d", i));
    residualStats[1][i][0] =
        MeasureResidual(hposXY[i][0], Form("resXY_X%d", i));
    residualStats[1][i][1] =
        MeasureResidual(hposXY[i][1], Form("resXY_Y%d", i));
    res[i][0] = residualStats[0][i][0].adaptiveSigma;
    res[i][1] = residualStats[0][i][1].adaptiveSigma;
    resXY[i][0] = residualStats[1][i][0].adaptiveSigma;
    resXY[i][1] = residualStats[1][i][1].adaptiveSigma;
    efficiencyX[i][2] = DivideSafe(efficiencyX[i][1], efficiencyX[i][0]);
    efficiencyY[i][2] = DivideSafe(efficiencyY[i][1], efficiencyY[i][0]);
  }

  frecFile->cd();
  TTree *residualSummary = new TTree(
      "TrackerResidualSummary",
      "Adaptive/fixed Gaussian and robust Tracker residual statistics");
  int summarySample = 0;
  int summaryDetector = 0;
  int summaryAxis = 0;
  long long summaryEntries = 0;
  double summaryRmsUm = 0.0;
  double summaryQuantile16Um = 0.0;
  double summaryQuantile84Um = 0.0;
  double summaryCentral68CenterUm = 0.0;
  double summaryCentral68HalfWidthUm = 0.0;
  int summaryAdaptiveStatus = -1;
  double summaryAdaptiveMeanUm = 0.0;
  double summaryAdaptiveMeanErrorUm = 0.0;
  double summaryAdaptiveSigmaUm = 0.0;
  double summaryAdaptiveSigmaErrorUm = 0.0;
  double summaryAdaptiveRangeLowUm = 0.0;
  double summaryAdaptiveRangeHighUm = 0.0;
  double summaryAdaptiveChi2Ndf = 0.0;
  int summaryFixedStatus = -1;
  double summaryFixedMeanUm = 0.0;
  double summaryFixedMeanErrorUm = 0.0;
  double summaryFixedSigmaUm = 0.0;
  double summaryFixedSigmaErrorUm = 0.0;
  double summaryFixedRangeLowUm = 0.0;
  double summaryFixedRangeHighUm = 0.0;
  double summaryFixedChi2Ndf = 0.0;
  residualSummary->Branch("sample", &summarySample, "sample/I");
  residualSummary->Branch("detector", &summaryDetector, "detector/I");
  residualSummary->Branch("axis", &summaryAxis, "axis/I");
  residualSummary->Branch("entries", &summaryEntries, "entries/L");
  residualSummary->Branch("rms_um", &summaryRmsUm, "rms_um/D");
  residualSummary->Branch("quantile16_um", &summaryQuantile16Um,
                          "quantile16_um/D");
  residualSummary->Branch("quantile84_um", &summaryQuantile84Um,
                          "quantile84_um/D");
  residualSummary->Branch("central68_center_um",
                          &summaryCentral68CenterUm,
                          "central68_center_um/D");
  residualSummary->Branch("central68_half_width_um",
                          &summaryCentral68HalfWidthUm,
                          "central68_half_width_um/D");
  residualSummary->Branch("adaptive_status", &summaryAdaptiveStatus,
                          "adaptive_status/I");
  residualSummary->Branch("adaptive_mean_um", &summaryAdaptiveMeanUm,
                          "adaptive_mean_um/D");
  residualSummary->Branch("adaptive_mean_error_um",
                          &summaryAdaptiveMeanErrorUm,
                          "adaptive_mean_error_um/D");
  residualSummary->Branch("adaptive_sigma_um", &summaryAdaptiveSigmaUm,
                          "adaptive_sigma_um/D");
  residualSummary->Branch("adaptive_sigma_error_um",
                          &summaryAdaptiveSigmaErrorUm,
                          "adaptive_sigma_error_um/D");
  residualSummary->Branch("adaptive_range_low_um",
                          &summaryAdaptiveRangeLowUm,
                          "adaptive_range_low_um/D");
  residualSummary->Branch("adaptive_range_high_um",
                          &summaryAdaptiveRangeHighUm,
                          "adaptive_range_high_um/D");
  residualSummary->Branch("adaptive_chi2_ndf", &summaryAdaptiveChi2Ndf,
                          "adaptive_chi2_ndf/D");
  residualSummary->Branch("fixed_status", &summaryFixedStatus,
                          "fixed_status/I");
  residualSummary->Branch("fixed_mean_um", &summaryFixedMeanUm,
                          "fixed_mean_um/D");
  residualSummary->Branch("fixed_mean_error_um", &summaryFixedMeanErrorUm,
                          "fixed_mean_error_um/D");
  residualSummary->Branch("fixed_sigma_um", &summaryFixedSigmaUm,
                          "fixed_sigma_um/D");
  residualSummary->Branch("fixed_sigma_error_um",
                          &summaryFixedSigmaErrorUm,
                          "fixed_sigma_error_um/D");
  residualSummary->Branch("fixed_range_low_um", &summaryFixedRangeLowUm,
                          "fixed_range_low_um/D");
  residualSummary->Branch("fixed_range_high_um", &summaryFixedRangeHighUm,
                          "fixed_range_high_um/D");
  residualSummary->Branch("fixed_chi2_ndf", &summaryFixedChi2Ndf,
                          "fixed_chi2_ndf/D");

  const TString summaryTsvName =
      GetPath(recName) + "/Tracker-step4-residual-summary.tsv";
  ofstream summaryTsv(summaryTsvName.Data());
  summaryTsv << "sample\tdetector\taxis\tentries\tadaptive_mu_um"
                "\tadaptive_sigma_um\tadaptive_sigma_error_um"
                "\tadaptive_range_low_um\tadaptive_range_high_um"
                "\tadaptive_status\tadaptive_chi2_ndf\tfixed_mu_um"
                "\tfixed_sigma_um\tfixed_sigma_error_um"
                "\tfixed_range_low_um\tfixed_range_high_um"
                "\tfixed_status\tfixed_chi2_ndf\trms_um"
                "\tcentral68_center_um\tcentral68_half_width_um"
                "\tquantile16_um\tquantile84_um\n";
  summaryTsv << std::setprecision(10);
  for (int sample = 0; sample < 2; ++sample)
  {
    for (int detector = 0; detector < MMdetN; ++detector)
    {
      for (int axis = 0; axis < 2; ++axis)
      {
        const ResidualStatistics &stats =
            residualStats[sample][detector][axis];
        summarySample = sample;
        summaryDetector = detector + 1;
        summaryAxis = axis;
        summaryEntries = stats.entries;
        summaryRmsUm = stats.rms;
        summaryQuantile16Um = stats.quantile16;
        summaryQuantile84Um = stats.quantile84;
        summaryCentral68CenterUm = stats.central68Center;
        summaryCentral68HalfWidthUm = stats.central68HalfWidth;
        summaryAdaptiveStatus = stats.adaptiveStatus;
        summaryAdaptiveMeanUm = stats.adaptiveMean;
        summaryAdaptiveMeanErrorUm = stats.adaptiveMeanError;
        summaryAdaptiveSigmaUm = stats.adaptiveSigma;
        summaryAdaptiveSigmaErrorUm = stats.adaptiveSigmaError;
        summaryAdaptiveRangeLowUm = stats.adaptiveRangeLow;
        summaryAdaptiveRangeHighUm = stats.adaptiveRangeHigh;
        summaryAdaptiveChi2Ndf = stats.adaptiveChi2Ndf;
        summaryFixedStatus = stats.fixedStatus;
        summaryFixedMeanUm = stats.fixedMean;
        summaryFixedMeanErrorUm = stats.fixedMeanError;
        summaryFixedSigmaUm = stats.fixedSigma;
        summaryFixedSigmaErrorUm = stats.fixedSigmaError;
        summaryFixedRangeLowUm = stats.fixedRangeLow;
        summaryFixedRangeHighUm = stats.fixedRangeHigh;
        summaryFixedChi2Ndf = stats.fixedChi2Ndf;
        residualSummary->Fill();

        const char *sampleName = sample == 0 ? "per_axis" : "exact_xy";
        const char *axisName = axis == 0 ? "X" : "Y";
        summaryTsv << sampleName << '\t' << detector + 1 << '\t'
                   << axisName << '\t' << stats.entries << '\t'
                   << stats.adaptiveMean << '\t' << stats.adaptiveSigma
                   << '\t' << stats.adaptiveSigmaError << '\t'
                   << stats.adaptiveRangeLow << '\t'
                   << stats.adaptiveRangeHigh << '\t'
                   << stats.adaptiveStatus << '\t'
                   << stats.adaptiveChi2Ndf << '\t' << stats.fixedMean
                   << '\t' << stats.fixedSigma << '\t'
                   << stats.fixedSigmaError << '\t' << stats.fixedRangeLow
                   << '\t' << stats.fixedRangeHigh << '\t'
                   << stats.fixedStatus << '\t' << stats.fixedChi2Ndf
                   << '\t' << stats.rms << '\t' << stats.central68Center
                   << '\t' << stats.central68HalfWidth << '\t'
                   << stats.quantile16 << '\t' << stats.quantile84 << '\n';
        op << "ResidualSummary sample=" << sampleName
           << " tracker=" << detector + 1 << " axis=" << axisName
           << " entries=" << stats.entries
           << " adaptive_mu_um=" << stats.adaptiveMean
           << " adaptive_sigma_um=" << stats.adaptiveSigma
           << " adaptive_range_um=[" << stats.adaptiveRangeLow << ","
           << stats.adaptiveRangeHigh << "]"
           << " adaptive_status=" << stats.adaptiveStatus
           << " fixed_mu_um=" << stats.fixedMean
           << " fixed_sigma_um=" << stats.fixedSigma
           << " fixed_range_um=[" << stats.fixedRangeLow << ","
           << stats.fixedRangeHigh << "]"
           << " fixed_status=" << stats.fixedStatus
           << " rms_um=" << stats.rms
           << " central68_um=[" << stats.quantile16 << ","
           << stats.quantile84 << "]" << endl;
      }
    }
  }
  summaryTsv.close();
  residualSummary->Write();

  TTree *efficiencySummary = new TTree(
      "TrackerEfficiencySummary",
      "Tracker efficiency with exact-2D references and +/-50 mm geometry");
  int efficiencyDetector = 0;
  int efficiencyAxis = 0;
  long long efficiencyDenominator = 0;
  long long efficiencyNumerator = 0;
  double efficiencyValue = 0.0;
  double referenceHalfWidthMm = kGeometryHalfWidthMm;
  double predictionHalfWidthMm = kGeometryHalfWidthMm;
  double residualAcceptanceMm = kResidualAcceptanceMm;
  efficiencySummary->Branch("detector", &efficiencyDetector,
                            "detector/I");
  efficiencySummary->Branch("axis", &efficiencyAxis, "axis/I");
  efficiencySummary->Branch("denominator", &efficiencyDenominator,
                            "denominator/L");
  efficiencySummary->Branch("numerator", &efficiencyNumerator,
                            "numerator/L");
  efficiencySummary->Branch("efficiency", &efficiencyValue,
                            "efficiency/D");
  efficiencySummary->Branch("reference_half_width_mm",
                            &referenceHalfWidthMm,
                            "reference_half_width_mm/D");
  efficiencySummary->Branch("prediction_half_width_mm",
                            &predictionHalfWidthMm,
                            "prediction_half_width_mm/D");
  efficiencySummary->Branch("residual_acceptance_mm",
                            &residualAcceptanceMm,
                            "residual_acceptance_mm/D");
  const TString efficiencyTsvName =
      GetPath(recName) + "/Tracker-step4-efficiency-summary.tsv";
  ofstream efficiencyTsv(efficiencyTsvName.Data());
  efficiencyTsv << "detector\taxis\tdenominator\tnumerator\tefficiency"
                   "\treference_half_width_mm\tprediction_half_width_mm"
                   "\tresidual_acceptance_mm\n";
  efficiencyTsv << std::setprecision(10);
  for (int detector = 0; detector < MMdetN; ++detector)
  {
    for (int axis = 0; axis < 2; ++axis)
    {
      const double (*efficiency)[3] = axis == 0 ? efficiencyX : efficiencyY;
      efficiencyDetector = detector + 1;
      efficiencyAxis = axis;
      efficiencyDenominator =
          static_cast<long long>(efficiency[detector][0]);
      efficiencyNumerator =
          static_cast<long long>(efficiency[detector][1]);
      efficiencyValue = efficiency[detector][2];
      efficiencySummary->Fill();
      const char *axisName = axis == 0 ? "X" : "Y";
      efficiencyTsv << detector + 1 << '\t' << axisName << '\t'
                    << efficiencyDenominator << '\t' << efficiencyNumerator
                    << '\t' << efficiencyValue << '\t'
                    << referenceHalfWidthMm << '\t'
                    << predictionHalfWidthMm << '\t'
                    << residualAcceptanceMm << '\n';
      op << "EfficiencySummary tracker=" << detector + 1
         << " axis=" << axisName
         << " denominator=" << efficiencyDenominator
         << " numerator=" << efficiencyNumerator
         << " efficiency=" << efficiencyValue
         << " reference_half_width_mm=" << referenceHalfWidthMm
         << " prediction_half_width_mm=" << predictionHalfWidthMm
         << " residual_acceptance_mm=" << residualAcceptanceMm << endl;
    }
  }
  efficiencyTsv.close();
  efficiencySummary->Write();

  // for (int i = 0; i < MMdetN; i++)
  // {
  //   {
  //     // 先做单高斯拟合获取初始参数
  //     hpos[i][0]->Fit("gaus", "Q0", "", -2000, 2000);
  //     TF1 *gausFit_x = hpos[i][0]->GetFunction("gaus");
  //     double amp1 = gausFit_x->GetParameter(0);
  //     double mean1 = gausFit_x->GetParameter(1);
  //     double sigma1 = gausFit_x->GetParameter(2);

  //     // 定义双高斯函数 - 使用唯一名称
  //     TF1 *dgaus_x = new TF1(Form("dgaus_x_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)",
  //                            -2000, 2000);
  //     dgaus_x->SetParameters(amp1, mean1, sigma1, amp1 / 10.0, 2.0 * sigma1);
  //     dgaus_x->SetParLimits(2, 0, 1000);
  //     dgaus_x->SetParLimits(4, 0, 2000);
  //     dgaus_x->SetParLimits(1, -100, 100);

  //     // 执行双高斯拟合
  //     hpos[i][0]->Fit(dgaus_x, "QR", "", -2000, 2000);

  //     // 获取拟合参数
  //     double A1 = dgaus_x->GetParameter(0);
  //     double mean_val = dgaus_x->GetParameter(1);
  //     double sig1 = dgaus_x->GetParameter(2);
  //     double A2 = dgaus_x->GetParameter(3);
  //     double sig2 = dgaus_x->GetParameter(4);

  //     // 保存两个sigma值
  //     res_sigma1[i][0] = sig1;
  //     res_sigma2[i][0] = sig2;

  //     // 计算有效分辨率
  //     double norm1 = A1 * sig1 * sqrt(2 * TMath::Pi());
  //     double norm2 = A2 * sig2 * sqrt(2 * TMath::Pi());
  //     double total_norm = norm1 + norm2;
  //     double moment2 = (norm1 * (mean_val * mean_val + sig1 * sig1) +
  //                       norm2 * (mean_val * mean_val + sig2 * sig2)) /
  //                      total_norm;
  //     double variance = moment2 - mean_val * mean_val;
  //     res[i][0] = sqrt(variance);
  //     area1[i][0] = norm1;
  //     area2[i][0] = norm2;
  //   }

  //   // Y方向 - 使用不同的函数名称
  //   {
  //     hpos[i][1]->Fit("gaus", "Q0", "", -2000, 2000);
  //     TF1 *gausFit_y = hpos[i][1]->GetFunction("gaus");
  //     double amp1 = gausFit_y->GetParameter(0);
  //     double mean1 = gausFit_y->GetParameter(1);
  //     double sigma1 = gausFit_y->GetParameter(2);

  //     // 定义双高斯函数 - 使用唯一名称
  //     TF1 *dgaus_y = new TF1(Form("dgaus_y_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)",
  //                            -2000, 2000);
  //     dgaus_y->SetParameters(amp1, mean1, sigma1, amp1 / 10.0, 2.0 * sigma1);
  //     dgaus_y->SetParLimits(2, 0, 1000);
  //     dgaus_y->SetParLimits(4, 0, 2000);
  //     dgaus_y->SetParLimits(1, -100, 100);

  //     hpos[i][1]->Fit(dgaus_y, "QR", "", -2000, 2000);

  //     double A1 = dgaus_y->GetParameter(0);
  //     double mean_val = dgaus_y->GetParameter(1);
  //     double sig1 = dgaus_y->GetParameter(2);
  //     double A2 = dgaus_y->GetParameter(3);
  //     double sig2 = dgaus_y->GetParameter(4);

  //     res_sigma1[i][1] = sig1;
  //     res_sigma2[i][1] = sig2;

  //     double norm1 = A1 * sig1 * sqrt(2 * TMath::Pi());
  //     double norm2 = A2 * sig2 * sqrt(2 * TMath::Pi());
  //     double total_norm = norm1 + norm2;
  //     double moment2 = (norm1 * (mean_val * mean_val + sig1 * sig1) +
  //                       norm2 * (mean_val * mean_val + sig2 * sig2)) /
  //                      total_norm;
  //     double variance = moment2 - mean_val * mean_val;
  //     res[i][1] = sqrt(variance);

  //     area1[i][1] = norm1;
  //     area2[i][1] = norm2;
  //   }

  //   efficiencyX[i][2] = efficiencyX[i][1] / efficiencyX[i][0];
  //   efficiencyY[i][2] = efficiencyY[i][1] / efficiencyY[i][0];
  // }

  for (int i = 0; i < MMdetN; i++)
  {
    bres[i][0] = GaussianParameterSafe(bhpos[i][0], 2, -2000, 2000, -1, Form("resbX%d", i));
    bres[i][1] = GaussianParameterSafe(bhpos[i][1], 2, -2000, 2000, -1, Form("resbY%d", i));
    efficiencybX[i][2] = DivideSafe(efficiencybX[i][1], efficiencybX[i][0]);
    efficiencybY[i][2] = DivideSafe(efficiencybY[i][1], efficiencybY[i][0]);
  }

  // for (int i = 0; i < 4; i++)
  // {
  //   {
  //     // 先做单高斯拟合获取初始参数
  //     bhpos[i][0]->Fit("gaus", "Q0", "", -2000, 2000);
  //     TF1 *gausFit_x = bhpos[i][0]->GetFunction("gaus");
  //     double amp1 = gausFit_x->GetParameter(0);
  //     double mean1 = gausFit_x->GetParameter(1);
  //     double sigma1 = gausFit_x->GetParameter(2);

  //     // 定义双高斯函数 - 使用唯一名称
  //     TF1 *dgaus_x = new TF1(Form("dgaus_x_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)",
  //                            -2000, 2000);
  //     dgaus_x->SetParameters(amp1, mean1, sigma1, amp1 / 10.0, 2.0 * sigma1);
  //     dgaus_x->SetParLimits(2, 0, 1000);
  //     dgaus_x->SetParLimits(4, 0, 2000);
  //     dgaus_x->SetParLimits(1, -100, 100);

  //     // 执行双高斯拟合
  //     bhpos[i][0]->Fit(dgaus_x, "QR", "", -2000, 2000);

  //     // 获取拟合参数
  //     double A1 = dgaus_x->GetParameter(0);
  //     double mean_val = dgaus_x->GetParameter(1);
  //     double sig1 = dgaus_x->GetParameter(2);
  //     double A2 = dgaus_x->GetParameter(3);
  //     double sig2 = dgaus_x->GetParameter(4);

  //     // 保存两个sigma值
  //     bres_sigma1[i][0] = sig1;
  //     bres_sigma2[i][0] = sig2;

  //     // 计算有效分辨率
  //     double norm1 = A1 * sig1 * sqrt(2 * TMath::Pi());
  //     double norm2 = A2 * sig2 * sqrt(2 * TMath::Pi());
  //     double total_norm = norm1 + norm2;
  //     double moment2 = (norm1 * (mean_val * mean_val + sig1 * sig1) +
  //                       norm2 * (mean_val * mean_val + sig2 * sig2)) /
  //                      total_norm;
  //     double variance = moment2 - mean_val * mean_val;
  //     bres[i][0] = sqrt(variance);

  //     barea1[i][0] = norm1;
  //     barea2[i][0] = norm2;
  //   }

  //   // Y方向 - 使用不同的函数名称
  //   {
  //     bhpos[i][1]->Fit("gaus", "Q0", "", -2000, 2000);
  //     TF1 *gausFit_y = bhpos[i][1]->GetFunction("gaus");
  //     double amp1 = gausFit_y->GetParameter(0);
  //     double mean1 = gausFit_y->GetParameter(1);
  //     double sigma1 = gausFit_y->GetParameter(2);

  //     // 定义双高斯函数 - 使用唯一名称
  //     TF1 *dgaus_y = new TF1(Form("dgaus_y_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)",
  //                            -2000, 2000);
  //     dgaus_y->SetParameters(amp1, mean1, sigma1, amp1 / 10.0, 2.0 * sigma1);
  //     dgaus_y->SetParLimits(2, 0, 1000);
  //     dgaus_y->SetParLimits(4, 0, 2000);
  //     dgaus_y->SetParLimits(1, -100, 100);

  //     bhpos[i][1]->Fit(dgaus_y, "QR", "", -2000, 2000);

  //     double A1 = dgaus_y->GetParameter(0);
  //     double mean_val = dgaus_y->GetParameter(1);
  //     double sig1 = dgaus_y->GetParameter(2);
  //     double A2 = dgaus_y->GetParameter(3);
  //     double sig2 = dgaus_y->GetParameter(4);

  //     bres_sigma1[i][1] = sig1;
  //     bres_sigma2[i][1] = sig2;

  //     double norm1 = A1 * sig1 * sqrt(2 * TMath::Pi());
  //     double norm2 = A2 * sig2 * sqrt(2 * TMath::Pi());
  //     double total_norm = norm1 + norm2;
  //     double moment2 = (norm1 * (mean_val * mean_val + sig1 * sig1) +
  //                       norm2 * (mean_val * mean_val + sig2 * sig2)) /
  //                      total_norm;
  //     double variance = moment2 - mean_val * mean_val;
  //     bres[i][1] = sqrt(variance);

  //     barea1[i][1] = norm1;
  //     barea2[i][1] = norm2;
  //   }

  //   efficiencybX[i][2] = efficiencybX[i][1] / efficiencybX[i][0];
  //   efficiencybY[i][2] = efficiencybY[i][1] / efficiencybY[i][0];
  // }

  // dtrackangle->Fit("gaus", "QR", "", -1, 1);
  double dtrackanglesigma = 0; // dtrackangle->GetFunction("gaus")->GetParameter(2);

  TF1 *reshitXFit = FitGaussianSafe(reshit[0], -6000, 6000, "reshitX");
  TF1 *reshitYFit = FitGaussianSafe(reshit[1], -6000, 6000, "reshitY");

  TCanvas *cMMeff = new TCanvas("Trackereff", "Trackereff");
  cMMeff->Divide(2, MMdetN);
  TLatex *la;
  for (int i = 0; i < MMdetN; i++)
  {
    cMMeff->cd(2 * i + 1);
    hpos[i][0]->GetXaxis()->SetRangeUser(-4000.0, 4000.0);
    hpos[i][0]->Draw();
    if (residualStats[0][i][0].fixedFunction)
      residualStats[0][i][0].fixedFunction->Draw("same");
    if (residualStats[0][i][0].adaptiveFunction)
      residualStats[0][i][0].adaptiveFunction->Draw("same");
    la = CreatLatex(Form("#color[2]{#sigma_{core}=%.0f#mum}",
                         residualStats[0][i][0].adaptiveSigma), 0.55, 0.82);
    la->Draw("same");
    la = CreatLatex(Form("#color[4]{#sigma_{#pm2mm}=%.0f#mum}",
                         residualStats[0][i][0].fixedSigma), 0.55, 0.74);
    la->Draw("same");
    la = CreatLatex(Form("RMS=%.0f#mum, 68%%/2=%.0f#mum",
                         residualStats[0][i][0].rms,
                         residualStats[0][i][0].central68HalfWidth),
                    0.45, 0.66);
    la->Draw("same");
    la = CreatLatex(Form("#mu_{core}=%.0f#mum, status=%d",
                         residualStats[0][i][0].adaptiveMean,
                         residualStats[0][i][0].adaptiveStatus),
                    0.45, 0.58);
    la->Draw("same");
    la = CreatLatex(Form("Eff=%.1f%%", efficiencyX[i][2] * 100), 0.16, 0.82);
    la->Draw("same");

    cMMeff->cd(2 * i + 2);
    hpos[i][1]->GetXaxis()->SetRangeUser(-4000.0, 4000.0);
    hpos[i][1]->Draw();
    if (residualStats[0][i][1].fixedFunction)
      residualStats[0][i][1].fixedFunction->Draw("same");
    if (residualStats[0][i][1].adaptiveFunction)
      residualStats[0][i][1].adaptiveFunction->Draw("same");
    la = CreatLatex(Form("#color[2]{#sigma_{core}=%.0f#mum}",
                         residualStats[0][i][1].adaptiveSigma), 0.55, 0.82);
    la->Draw("same");
    la = CreatLatex(Form("#color[4]{#sigma_{#pm2mm}=%.0f#mum}",
                         residualStats[0][i][1].fixedSigma), 0.55, 0.74);
    la->Draw("same");
    la = CreatLatex(Form("RMS=%.0f#mum, 68%%/2=%.0f#mum",
                         residualStats[0][i][1].rms,
                         residualStats[0][i][1].central68HalfWidth),
                    0.45, 0.66);
    la->Draw("same");
    la = CreatLatex(Form("#mu_{core}=%.0f#mum, status=%d",
                         residualStats[0][i][1].adaptiveMean,
                         residualStats[0][i][1].adaptiveStatus),
                    0.45, 0.58);
    la->Draw("same");
    la = CreatLatex(Form("Eff=%.1f%%", efficiencyY[i][2] * 100), 0.16, 0.82);
    la->Draw("same");
  }

  TCanvas *cXYResidual = new TCanvas("TrackerXYResidual", "TrackerXYResidual");
  cXYResidual->Divide(2, MMdetN);
  for (int i = 0; i < MMdetN; i++)
  {
    cXYResidual->cd(2 * i + 1);
    hposXY[i][0]->GetXaxis()->SetRangeUser(-4000.0, 4000.0);
    hposXY[i][0]->Draw();
    if (residualStats[1][i][0].fixedFunction)
      residualStats[1][i][0].fixedFunction->Draw("same");
    if (residualStats[1][i][0].adaptiveFunction)
      residualStats[1][i][0].adaptiveFunction->Draw("same");
    la = CreatLatex(Form("#color[2]{#sigma_{core}=%.0f#mum}",
                         residualStats[1][i][0].adaptiveSigma), 0.55, 0.82);
    la->Draw("same");
    la = CreatLatex(Form("#color[4]{#sigma_{#pm2mm}=%.0f#mum}",
                         residualStats[1][i][0].fixedSigma), 0.55, 0.74);
    la->Draw("same");
    la = CreatLatex(Form("RMS=%.0f#mum, 68%%/2=%.0f#mum",
                         residualStats[1][i][0].rms,
                         residualStats[1][i][0].central68HalfWidth),
                    0.45, 0.66);
    la->Draw("same");
    la = CreatLatex(Form("#mu_{core}=%.0f#mum, status=%d",
                         residualStats[1][i][0].adaptiveMean,
                         residualStats[1][i][0].adaptiveStatus),
                    0.45, 0.58);
    la->Draw("same");
    la = CreatLatex("3-layer exact-one X/Y", 0.16, 0.82);
    la->Draw("same");

    cXYResidual->cd(2 * i + 2);
    hposXY[i][1]->GetXaxis()->SetRangeUser(-4000.0, 4000.0);
    hposXY[i][1]->Draw();
    if (residualStats[1][i][1].fixedFunction)
      residualStats[1][i][1].fixedFunction->Draw("same");
    if (residualStats[1][i][1].adaptiveFunction)
      residualStats[1][i][1].adaptiveFunction->Draw("same");
    la = CreatLatex(Form("#color[2]{#sigma_{core}=%.0f#mum}",
                         residualStats[1][i][1].adaptiveSigma), 0.55, 0.82);
    la->Draw("same");
    la = CreatLatex(Form("#color[4]{#sigma_{#pm2mm}=%.0f#mum}",
                         residualStats[1][i][1].fixedSigma), 0.55, 0.74);
    la->Draw("same");
    la = CreatLatex(Form("RMS=%.0f#mum, 68%%/2=%.0f#mum",
                         residualStats[1][i][1].rms,
                         residualStats[1][i][1].central68HalfWidth),
                    0.45, 0.66);
    la->Draw("same");
    la = CreatLatex(Form("#mu_{core}=%.0f#mum, status=%d",
                         residualStats[1][i][1].adaptiveMean,
                         residualStats[1][i][1].adaptiveStatus),
                    0.45, 0.58);
    la->Draw("same");
    la = CreatLatex("3-layer exact-one X/Y", 0.16, 0.82);
    la->Draw("same");
  }

  // for (int i = 0; i < MMdetN; i++)
  // {

  //   cMMeff->cd(2 * i + 1);
  //   hpos[i][0]->Draw();
  //   TF1 *dgaus_x = hpos[i][0]->GetFunction(Form("dgaus_x_%d", i));
  //   if (dgaus_x)
  //   {
  //     TF1 *gaus1_x = new TF1(Form("gaus1_x_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2)",
  //                            -2000, 2000);
  //     gaus1_x->SetParameters(dgaus_x->GetParameter(0),
  //                            dgaus_x->GetParameter(1),
  //                            dgaus_x->GetParameter(2));
  //     gaus1_x->SetLineColor(kBlue);
  //     gaus1_x->SetLineStyle(1);
  //     gaus1_x->SetLineWidth(2);
  //     gaus1_x->Draw("same");

  //     TF1 *gaus2_x = new TF1(Form("gaus2_x_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2)",
  //                            -2000, 2000);
  //     gaus2_x->SetParameters(dgaus_x->GetParameter(3),
  //                            dgaus_x->GetParameter(1),
  //                            dgaus_x->GetParameter(4));
  //     gaus2_x->SetLineColor(kGreen);
  //     gaus2_x->SetLineStyle(1);
  //     gaus2_x->SetLineWidth(2);
  //     gaus2_x->Draw("same");

  //     dgaus_x->SetLineColor(kRed);
  //     dgaus_x->SetLineWidth(2);
  //     dgaus_x->Draw("same");
  //   }
  //   la = CreatLatex(Form("#color[2]{area1 = %.0f, area2 = %.0f, p = %.0f}", area1[i][0], area2[i][0], area2[i][0]/(area1[i][0] + area2[i][0])), 0.6, 0.35);
  //   la->Draw("same");
  //   cout << "x" << area2[i][0]/(area1[i][0] + area2[i][0]) << endl;
  //   TLatex *la2 = CreatLatex(Form("#color[4]{#sigma_{1}=%.0f #mum}", res_sigma1[i][0]), 0.6, 0.30);
  //   la2->Draw("same");

  //   TLatex *la3 = CreatLatex(Form("#color[6]{#sigma_{2}=%.0f #mum}", res_sigma2[i][0]), 0.6, 0.25);
  //   la3->Draw("same");
  //   la = CreatLatex(Form("Eff=%.1f%%", efficiencyX[i][2] * 100), 0.2, 0.75);
  //   la->Draw("same");

  //   cMMeff->cd(2 * i + 2);
  //   hpos[i][1]->Draw();
  //   TF1 *dgaus_y = hpos[i][1]->GetFunction(Form("dgaus_y_%d", i));
  //   if (dgaus_y)
  //   {
  //     TF1 *gaus1_y = new TF1(Form("gaus1_y_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2)",
  //                            -2000, 2000);
  //     gaus1_y->SetParameters(dgaus_y->GetParameter(0),
  //                            dgaus_y->GetParameter(1),
  //                            dgaus_y->GetParameter(2));
  //     gaus1_y->SetLineColor(kBlue);
  //     gaus1_y->SetLineStyle(1);
  //     gaus1_y->SetLineWidth(2);
  //     gaus1_y->Draw("same");

  //     TF1 *gaus2_y = new TF1(Form("gaus2_y_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2)",
  //                            -2000, 2000);
  //     gaus2_y->SetParameters(dgaus_y->GetParameter(3),
  //                            dgaus_y->GetParameter(1),
  //                            dgaus_y->GetParameter(4));
  //     gaus2_y->SetLineColor(kGreen);
  //     gaus2_y->SetLineStyle(1);
  //     gaus2_y->SetLineWidth(2);
  //     gaus2_y->Draw("same");

  //     dgaus_y->SetLineColor(kRed);
  //     dgaus_y->SetLineWidth(2);
  //     dgaus_y->Draw("same");
  //   }
  //   la = CreatLatex(Form("#color[2]{area1 = %.0f, area2 = %.0f, p = %.0f}", area1[i][1], area2[i][1], area2[i][1]/(area1[i][1] + area2[i][1])), 0.6, 0.35);
  //   cout << "y" << area2[i][1]/(area1[i][1] + area2[i][1]) << endl;
  //   la->Draw("same");
  //   la2 = CreatLatex(Form("#color[4]{#sigma_{1}=%.0f #mum}", res_sigma1[i][1]), 0.6, 0.30);
  //   la2->Draw("same");

  //   la3 = CreatLatex(Form("#color[6]{#sigma_{2}=%.0f #mum}", res_sigma2[i][1]), 0.6, 0.25);
  //   la3->Draw("same");
  //   la = CreatLatex(Form("Eff=%.1f%%", efficiencyY[i][2] * 100), 0.2, 0.75);
  //   la->Draw("same");
  // }

  TCanvas *cbMMeff = new TCanvas("bTrackereff", "bTrackereff");
  cbMMeff->Divide(2, MMdetN);
  for (int i = 0; i < MMdetN; i++)
  {
    cbMMeff->cd(2 * i + 1);
    bhpos[i][0]->Draw();
    la = CreatLatex(Form("#color[2]{#sigma=%.0f#mum}", bres[i][0]), 0.6, 0.35);
    la->Draw("same");
    la = CreatLatex(Form("Eff=%.1f%%", efficiencybX[i][2] * 100), 0.2, 0.75);
    la->Draw("same");

    cbMMeff->cd(2 * i + 2);
    bhpos[i][1]->Draw();
    la = CreatLatex(Form("#color[2]{#sigma=%.0f#mum}", bres[i][1]), 0.6, 0.35);
    la->Draw("same");
    la = CreatLatex(Form("Eff=%.1f%%", efficiencybY[i][2] * 100), 0.2, 0.75);
    la->Draw("same");
  }

  // for (int i = 0; i < 4; i++)
  // {

  //   cbMMeff->cd(2 * i + 1);
  //   bhpos[i][0]->Draw();
  //   TF1 *dgaus_x = bhpos[i][0]->GetFunction(Form("dgaus_x_%d", i));
  //   if (dgaus_x)
  //   {
  //     TF1 *gaus1_x = new TF1(Form("gaus1_x_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2)",
  //                            -2000, 2000);
  //     gaus1_x->SetParameters(dgaus_x->GetParameter(0),
  //                            dgaus_x->GetParameter(1),
  //                            dgaus_x->GetParameter(2));
  //     gaus1_x->SetLineColor(kBlue);
  //     gaus1_x->SetLineStyle(1);
  //     gaus1_x->SetLineWidth(2);
  //     gaus1_x->Draw("same");

  //     TF1 *gaus2_x = new TF1(Form("gaus2_x_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2)",
  //                            -2000, 2000);
  //     gaus2_x->SetParameters(dgaus_x->GetParameter(3),
  //                            dgaus_x->GetParameter(1),
  //                            dgaus_x->GetParameter(4));
  //     gaus2_x->SetLineColor(kGreen);
  //     gaus2_x->SetLineStyle(1);
  //     gaus2_x->SetLineWidth(2);
  //     gaus2_x->Draw("same");

  //     dgaus_x->SetLineColor(kRed);
  //     dgaus_x->SetLineWidth(2);
  //     dgaus_x->Draw("same");
  //   }
  //   la = CreatLatex(Form("#color[2]{area1 = %.0f, area2 = %.0f, p = %.0f}", barea1[i][0], barea2[i][0], barea2[i][0]/(barea1[i][0] + barea2[i][0])), 0.6, 0.35);
  //   la->Draw("same");
  //   cout << "bx" << barea2[i][0]/(barea1[i][0] + barea2[i][0]) << endl;
  //   TLatex *la2 = CreatLatex(Form("#color[4]{#sigma_{1}=%.0f #mum}", bres_sigma1[i][0]), 0.6, 0.30);
  //   la2->Draw("same");

  //   TLatex *la3 = CreatLatex(Form("#color[6]{#sigma_{2}=%.0f #mum}", bres_sigma2[i][0]), 0.6, 0.25);
  //   la3->Draw("same");
  //   la = CreatLatex(Form("Eff=%.1f%%", efficiencybX[i][2] * 100), 0.2, 0.75);
  //   la->Draw("same");

  //   cbMMeff->cd(2 * i + 2);
  //   bhpos[i][1]->Draw();
  //   TF1 *dgaus_y = bhpos[i][1]->GetFunction(Form("dgaus_y_%d", i));
  //   if (dgaus_y)
  //   {
  //     TF1 *gaus1_y = new TF1(Form("gaus1_y_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2)",
  //                            -2000, 2000);
  //     gaus1_y->SetParameters(dgaus_y->GetParameter(0),
  //                            dgaus_y->GetParameter(1),
  //                            dgaus_y->GetParameter(2));
  //     gaus1_y->SetLineColor(kBlue);
  //     gaus1_y->SetLineStyle(1);
  //     gaus1_y->SetLineWidth(2);
  //     gaus1_y->Draw("same");

  //     TF1 *gaus2_y = new TF1(Form("gaus2_y_%d", i),
  //                            "[0]*exp(-0.5*((x-[1])/[2])^2)",
  //                            -2000, 2000);
  //     gaus2_y->SetParameters(dgaus_y->GetParameter(3),
  //                            dgaus_y->GetParameter(1),
  //                            dgaus_y->GetParameter(4));
  //     gaus2_y->SetLineColor(kGreen);
  //     gaus2_y->SetLineStyle(1);
  //     gaus2_y->SetLineWidth(2);
  //     gaus2_y->Draw("same");

  //     dgaus_y->SetLineColor(kRed);
  //     dgaus_y->SetLineWidth(2);
  //     dgaus_y->Draw("same");
  //   }
  //   la = CreatLatex(Form("#color[2]{area1 = %.0f, area2 = %.0f, p = %.0f}", barea1[i][1], barea2[i][1], barea2[i][1]/(barea1[i][1] + barea2[i][1])), 0.6, 0.35);
  //   la->Draw("same");
  //   cout << "by" << barea2[i][1]/(barea1[i][1] + barea2[i][1]) << endl;
  //   la2 = CreatLatex(Form("#color[4]{#sigma_{1}=%.0f #mum}", bres_sigma1[i][1]), 0.6, 0.30);
  //   la2->Draw("same");

  //   la3 = CreatLatex(Form("#color[6]{#sigma_{2}=%.0f #mum}", bres_sigma2[i][1]), 0.6, 0.25);
  //   la3->Draw("same");
  //   la = CreatLatex(Form("Eff=%.1f%%", efficiencybY[i][2] * 100), 0.2, 0.75);
  //   la->Draw("same");
  // }

  TCanvas *tangle = new TCanvas("angle", "angle");
  tangle->Divide(2, 2);
  tangle->cd(1);
  trackangle[0]->Draw();
  tangle->cd(2);
  trackangle[1]->Draw();
  tangle->cd(3);
  FitGaussianSafe(dtrackangle, -20, 20, "dtrackangle");
  //  double scatterangel = dtrackangle->GetFunction("gaus")->GetParameter(2);
  dtrackangle->Draw();
  //  la = CreatLatex(Form("#color[2]{#sigma=%.5f}", scatterangel), 0.6, 0.35);
  //  la->Draw("same");
  // la = CreatLatex(Form("#color[2]{#sigma=%f°}", dtrackanglesigma), 0.6, 0.35);
  // la->Draw("same");

  TCanvas *hitmaps = new TCanvas("hitmap", "hitmap");
  hitmaps->Divide(2, 2);
  for (int i = 0; i < 4; i++)
  {
    hitmaps->cd(i + 1);
    hitmap[i]->Draw();
  }

  TCanvas *fitKcanvas = new TCanvas("K", "K");
  fitKcanvas->Divide(2, 2);
  for (int i = 0; i < 4; i++)
  {
    fitKcanvas->cd(i + 1);
    FitGaussianSafe(fitK[i], -0.1, 0.1, Form("fitK%d", i));
    fitK[i]->Draw();
  }

  // for (int i = 0; i < 4; i++)
  // {
  //   fitKcanvas->cd(i + 1);
  //   fitK[i]->Fit("gaus", "Q0", "", -0.1, 0.1);
  //   TF1 *gausFit_x = fitK[i]->GetFunction("gaus");
  //   double amp1 = gausFit_x->GetParameter(0);
  //   double mean1 = gausFit_x->GetParameter(1);
  //   double sigma1 = gausFit_x->GetParameter(2);

  //   // 三高斯函数定义
  //   TF1 *tgaus_x = new TF1(Form("tgaus_x_%d", i),
  //                          "[0]*exp(-0.5*((x-[1])/[2])^2) + "      // 第一个高斯
  //                          "[3]*exp(-0.5*((x-[4])/[5])^2) + "      // 第二个高斯
  //                          "[6]*exp(-0.5*((x-[7])/[8])^2)",        // 第三个高斯
  //                          -0.1, 0.1);

  //   // 设置参数初始值
  //   // 第一个高斯：使用单高斯拟合结果
  //   // 第二个高斯：较小的振幅，稍大的sigma，mean偏移一点
  //   // 第三个高斯：更小的振幅，更大的sigma，mean再偏移一点
  //   tgaus_x->SetParameters(amp1, mean1, sigma1,                    // [0],[1],[2]
  //                         amp1 * 0.3, mean1 + 0.01, sigma1 * 2, // [3],[4],[5]
  //                         amp1 * 0.1, mean1 - 0.02, sigma1 * 4); // [6],[7],[8]

  //   // 设置参数限制
  //   tgaus_x->SetParLimits(2, 0, 100000);    // 第一个sigma
  //   tgaus_x->SetParLimits(5, 0, 200000);    // 第二个sigma
  //   tgaus_x->SetParLimits(8, 0, 300000);    // 第三个sigma
  //   tgaus_x->SetParLimits(1, -100, 100);    // 第一个mean
  //   tgaus_x->SetParLimits(4, -100, 100);    // 第二个mean
  //   tgaus_x->SetParLimits(7, -100, 100);    // 第三个mean

  //   // 执行三高斯拟合
  //   fitK[i]->Fit(tgaus_x, "QR", "", -0.1, 0.1);
  //   fitK[i]->Draw();
  // }

  TCanvas *reshits = new TCanvas("reshit", "reshit");
  reshits->Divide(2, 1);
  double reshitxfitmean = reshitXFit ? reshitXFit->GetParameter(1) : 0;
  double reshityfitmean = reshitYFit ? reshitYFit->GetParameter(1) : 0;
  cout << "reshitxfitmean: " << reshitxfitmean << " reshityfitmean: " << reshityfitmean << endl;

  TCanvas *FITKXY = new TCanvas("kxy", "kxy");
  FITKXY->Divide(2, 2);
  for (int i = 0; i < 4; i++)
  {
    FITKXY->cd(i + 1);
    KXY[i]->Draw();
  }

  reshits->cd(1);
  reshit[0]->Draw();
  reshits->cd(2);
  reshit[1]->Draw();

  frecFile->cd();
  TNamed("TrackerResidualSelection",
         "No track-angle cut. For each DUT, both reference Trackers must each "
         "have exactly one X and one Y cluster and reference points must obey "
         "|x|,|y|<50 mm. The predicted DUT point must obey |x|,|y|<50 mm. "
         "Trackereff keeps a responding DUT coordinate when |residual|<10 "
         "mm. TrackerXYResidual additionally requires the DUT to have exactly "
         "one X and one Y cluster and fills X/Y as a pair. Alignment is a "
         "three-layer joint dx/dy/rotationZ Minuit2 fit with Tracker1 fixed.")
      .Write();
  for (int i = 0; i < MMdetN; ++i)
  {
    hpos[i][0]->Write();
    hpos[i][1]->Write();
    hposXY[i][0]->Write();
    hposXY[i][1]->Write();
  }
  cMMeff->SaveAs(GetPath(recName) +
                 "/Trackereff_position50_joint_alignment.png");
  cXYResidual->SaveAs(GetPath(recName) +
                      "/TrackerXYResidual_position50_joint_alignment.png");
  cMMeff->Write();
  cXYResidual->Write();
  cbMMeff->Write();
  tangle->Write();
  hitmaps->Write();
  fitKcanvas->Write();
  reshits->Write();
  FITKXY->Write();
  frecFile->Flush();
  frecFile->Close();
  op.close();
}

void ReconstructTrack::InitData()
{
  useflag = 0; // useful
  richuseflag = 0;
  btofuseflag = 0;
  ecaluseflag = 0;
  memset(MMx, 0, sizeof(MMx));
  memset(MMy, 0, sizeof(MMy));
  memset(MMz, 0, sizeof(MMz));
  memset(MMsigx, 0, sizeof(MMsigx));
  memset(MMsigy, 0, sizeof(MMsigy));
  memset(XfitPars, 0, sizeof(XfitPars));
  memset(YfitPars, 0, sizeof(YfitPars));
  memset(X01fitPars, 0, sizeof(X01fitPars));
  memset(Y01fitPars, 0, sizeof(Y01fitPars));
  memset(X2345fitPars, 0, sizeof(X2345fitPars));
  memset(Y2345fitPars, 0, sizeof(Y2345fitPars));
  for (int i = 0; i < MMdetN; ++i)
  {
    predStripX[i] = -1.0;
    predStripY[i] = -1.0;
  }
  memset(richextraHit, 0, sizeof(richextraHit));
  memset(richtrackVec, 0, sizeof(richtrackVec));
  memset(btofextraHit, 0, sizeof(btofextraHit));
  memset(btoftrackVec, 0, sizeof(btoftrackVec));
  memset(ecalextraHit, 0, sizeof(ecalextraHit));
  memset(ecaltrackVec, 0, sizeof(ecaltrackVec));
}

bool ReconstructTrack::GetBTOFposition(double kx, double ky, double bx, double by, vector<double> btofxyz)
{
  double A = 1;
  double B = 1;
  double C = 1;
  double D = 1;
  double btofx = 1;
  double btofy = 1;
  double btofz = 1;
  btofxyz.push_back(btofx);
  btofxyz.push_back(btofy);
  btofxyz.push_back(btofz);
  return true;
}

bool ReconstructTrack::fit_track(vector<double> px, vector<double> py, double par[3])
{
  int n = px.size() < py.size() ? px.size() : py.size();
  if (n < MMdetN - 1)
    return false;
  double Mx = 0, My = 0, Mxx = 0, Myy = 0, Mxy = 0; // mean
  for (int i = 0; i < n; i++)
  {
    Mx += px[i] / n;
    My += py[i] / n;
    Mxx += px[i] * px[i] / n;
    Myy += py[i] * py[i] / n;
    Mxy += px[i] * py[i] / n;
  }
  /// y = K*x+B, R2 = correlation index squre
  par[0] = (Mxy - Mx * My) / (Mxx - Mx * Mx); // K
  par[1] = My - par[0] * Mx;                  // B
  for (int i = 0; i < n; i++)
  {
    TVector3 a = TVector3(px[i], py[i], 0);
    a.RotateZ(45 * deg); // B
  }
  double Mxtemp = 0, Mytemp = 0, Mxxtemp = 0, Myytemp = 0, Mxytemp = 0;
  for (int i = 0; i < n; i++)
  {
    double tx = px[i] * Cos(45 * deg) - py[i] * Sin(45 * deg);
    double ty = px[i] * Sin(45 * deg) + py[i] * Cos(45 * deg);
    Mxtemp += tx / n;
    Mytemp += ty / n;
    Mxxtemp += tx * tx / n;
    Myytemp += ty * ty / n;
    Mxytemp += tx * ty / n;
  }
  par[2] = (Mxytemp - Mxtemp * Mytemp) * (Mxytemp - Mxtemp * Mytemp) / (Mxxtemp - Mxtemp * Mxtemp) / (Myytemp - Mytemp * Mytemp); // R2

  return par[2] > fit_R2_cut;
}

bool ReconstructTrack::fit_track3(vector<double> px, vector<double> py, double par[3])
{
  int n = px.size() < py.size() ? px.size() : py.size();
  if (n < 3)
    return false;
  double Mx = 0, My = 0, Mxx = 0, Myy = 0, Mxy = 0; // mean
  for (int i = 0; i < n; i++)
  {
    Mx += px[i] / n;
    My += py[i] / n;
    Mxx += px[i] * px[i] / n;
    Myy += py[i] * py[i] / n;
    Mxy += px[i] * py[i] / n;
  }
  /// y = K*x+B, R2 = correlation index squre
  par[0] = (Mxy - Mx * My) / (Mxx - Mx * Mx); // K
  par[1] = My - par[0] * Mx;
  for (int i = 0; i < n; i++)
  {
    TVector3 a = TVector3(px[i], py[i], 0);
    a.RotateZ(45 * deg); // B
  }
  double Mxtemp = 0, Mytemp = 0, Mxxtemp = 0, Myytemp = 0, Mxytemp = 0;
  for (int i = 0; i < n; i++)
  {
    double tx = px[i] * Cos(45 * deg) - py[i] * Sin(45 * deg);
    double ty = px[i] * Sin(45 * deg) + py[i] * Cos(45 * deg);
    Mxtemp += tx / n;
    Mytemp += ty / n;
    Mxxtemp += tx * tx / n;
    Myytemp += ty * ty / n;
    Mxytemp += tx * ty / n;
  }
  par[2] = (Mxytemp - Mxtemp * Mytemp) * (Mxytemp - Mxtemp * Mytemp) / (Mxxtemp - Mxtemp * Mxtemp) / (Myytemp - Mytemp * Mytemp); // R2

  return par[2] > fit_R2_cut;
}

bool ReconstructTrack::fit_track2(vector<double> px, vector<double> py, double par[3])
{
  int n = px.size() < py.size() ? px.size() : py.size();
  if (n < 2)
    return false;
  double Mx = 0, My = 0, Mxx = 0, Myy = 0, Mxy = 0; // mean
  for (int i = 0; i < n; i++)
  {
    Mx += px[i] / n;
    My += py[i] / n;
    Mxx += px[i] * px[i] / n;
    Myy += py[i] * py[i] / n;
    Mxy += px[i] * py[i] / n;
  }
  /// y = K*x+B, R2 = correlation index squre
  par[0] = (Mxy - Mx * My) / (Mxx - Mx * Mx);                                     // K
  par[1] = My - par[0] * Mx;                                                      // B
  par[2] = (Mxy - Mx * My) * (Mxy - Mx * My) / (Mxx - Mx * Mx) / (Myy - My * My); // R2

  return true;
}

bool ReconstructTrack::JointAlignment(const string &par3Filename,
                                      const string &par6Filename)
{
  cout << " --> Begin three-layer joint alignment (dx, dy, rotationZ)."
       << endl;
  vector<JointAlignmentEvent> tracks;
  const int nentries = fdecTree[0]->GetEntriesFast();
  tracks.reserve(nentries);
  for (int entry = 0; entry < nentries; ++entry)
  {
    bool good = true;
    JointAlignmentEvent track;
    for (int detector = 0; detector < MMdetN; ++detector)
    {
      fdecTree[detector]->GetEntry(entry);
      if (!sig_x || !sig_y || cluster_count_x != 1 ||
          cluster_count_y != 1)
      {
        good = false;
        break;
      }
      track.local[detector][0] = x * 0.4 - centre;
      track.local[detector][1] = y * 0.4 - centre;
    }
    if (good)
      tracks.push_back(track);
  }

  JointAlignmentResult alignment;
  alignment.trackCount = static_cast<long long>(tracks.size());
  alignment.objectiveBefore = JointAlignmentObjective(tracks, alignment);
  if (tracks.size() < 20)
  {
    cerr << " [ERROR] Joint tracker alignment requires at least 20 three-layer "
            "single-X/Y-cluster events, found "
         << tracks.size() << "." << endl;
    return false;
  }

  frecFile->cd();
  TDirectory *alignmentDirectory = frecFile->mkdir("alignment");
  TDirectory *coarseDirectory =
      alignmentDirectory->mkdir("coarse_alignment");
  coarseDirectory->cd();
  for (int detector = 1; detector < MMdetN; ++detector)
  {
    TH1D *differenceX = new TH1D(
        Form("coarse_dx_tracker%d", detector + 1),
        Form("Tracker%d - Tracker1 X;#DeltaX [mm];events", detector + 1),
        400, -100.0, 100.0);
    TH1D *differenceY = new TH1D(
        Form("coarse_dy_tracker%d", detector + 1),
        Form("Tracker%d - Tracker1 Y;#DeltaY [mm];events", detector + 1),
        400, -100.0, 100.0);
    for (const JointAlignmentEvent &track : tracks)
    {
      differenceX->Fill(track.local[detector][0] - track.local[0][0]);
      differenceY->Fill(track.local[detector][1] - track.local[0][1]);
    }
    TF1 *fitX = FitGaussianSafe(differenceX, -100.0, 100.0,
                                Form("coarse_dx_tracker%d", detector + 1));
    TF1 *fitY = FitGaussianSafe(differenceY, -100.0, 100.0,
                                Form("coarse_dy_tracker%d", detector + 1));
    alignment.dx[detector] =
        -(fitX ? fitX->GetParameter(1) : differenceX->GetMean());
    alignment.dy[detector] =
        -(fitY ? fitY->GetParameter(1) : differenceY->GetMean());
    differenceX->Write();
    differenceY->Write();
  }
  alignment.objectiveCoarse = JointAlignmentObjective(tracks, alignment);
  alignment.status = FineJointAlignment(tracks, alignment);
  alignment.objectiveAfter = JointAlignmentObjective(tracks, alignment);

  if (!std::isfinite(alignment.objectiveAfter))
  {
    cerr << " [ERROR] Joint tracker alignment produced a non-finite objective."
         << endl;
    frecFile->cd();
    return false;
  }

  for (int detector = 0; detector < MMdetN; ++detector)
  {
    for (int parameter = 0; parameter < 6; ++parameter)
    {
      algPar3[detector][parameter] = 0.0;
      algPar6[detector][parameter] = 0.0;
    }
    algPar3[detector][0] = alignment.dx[detector];
    algPar3[detector][1] = alignment.dy[detector];
    algPar6[detector][0] = alignment.dx[detector];
    algPar6[detector][1] = alignment.dy[detector];
    algPar6[detector][5] = alignment.rotationZ[detector];
  }

  ofstream par3(par3Filename.c_str());
  ofstream par6(par6Filename.c_str());
  if (!par3.is_open() || !par6.is_open())
  {
    cerr << " [ERROR] Cannot write run-local tracker alignment files." << endl;
    frecFile->cd();
    return false;
  }
  par3 << std::setprecision(12);
  par6 << std::setprecision(12);
  for (int detector = 0; detector < MMdetN; ++detector)
  {
    for (int parameter = 0; parameter < 6; ++parameter)
    {
      par3 << algPar3[detector][parameter]
           << (parameter == 5 ? '\n' : ' ');
      par6 << algPar6[detector][parameter]
           << (parameter == 5 ? '\n' : ' ');
    }
  }
  par3.close();
  par6.close();

  alignmentDirectory->cd();
  TTree *alignmentTree = new TTree(
      "TrackerAlignment",
      "Three-layer joint dx/dy/rotationZ alignment; Tracker1 fixed");
  int detectorNumber = 0;
  double zMm = 0.0;
  double dxMm = 0.0;
  double dyMm = 0.0;
  double rotationZRad = 0.0;
  alignmentTree->Branch("detector", &detectorNumber, "detector/I");
  alignmentTree->Branch("z_mm", &zMm, "z_mm/D");
  alignmentTree->Branch("dx_mm", &dxMm, "dx_mm/D");
  alignmentTree->Branch("dy_mm", &dyMm, "dy_mm/D");
  alignmentTree->Branch("rotation_z_rad", &rotationZRad,
                        "rotation_z_rad/D");
  alignmentTree->Branch("track_count", &alignment.trackCount,
                        "track_count/L");
  alignmentTree->Branch("objective_before", &alignment.objectiveBefore,
                        "objective_before/D");
  alignmentTree->Branch("objective_coarse", &alignment.objectiveCoarse,
                        "objective_coarse/D");
  alignmentTree->Branch("objective_after", &alignment.objectiveAfter,
                        "objective_after/D");
  alignmentTree->Branch("status", &alignment.status, "status/I");
  for (int detector = 0; detector < MMdetN; ++detector)
  {
    detectorNumber = detector + 1;
    zMm = Zpos[detector];
    dxMm = alignment.dx[detector];
    dyMm = alignment.dy[detector];
    rotationZRad = alignment.rotationZ[detector];
    alignmentTree->Fill();
  }
  alignmentTree->Write();
  TNamed("TrackerAlignmentSelection",
         "All three Trackers must each have exactly one X cluster and one Y "
         "cluster. Tracker1 is fixed; Minuit2 jointly fits Tracker2/3 dx, dy "
         "and rotationZ using the mean X/Y straight-line residual chi2.")
      .Write();

  cout << " --> Joint alignment tracks=" << alignment.trackCount
       << ", status=" << alignment.status << ", objective "
       << alignment.objectiveBefore << " -> " << alignment.objectiveCoarse
       << " -> " << alignment.objectiveAfter << endl;
  cout << " --> Tracker       dx [mm]       dy [mm]    rotationZ [mrad]"
       << endl;
  for (int detector = 0; detector < MMdetN; ++detector)
  {
    cout << " --> " << std::setw(7) << detector + 1 << std::setw(14)
         << alignment.dx[detector] << std::setw(14)
         << alignment.dy[detector] << std::setw(20)
         << alignment.rotationZ[detector] * 1000.0 << endl;
  }
  cout << " --> Alignment written to " << par6Filename << endl;
  frecFile->cd();
  return true;
}

void ReconstructTrack::ReadAlignmentPars(string filename, int Npars)
{
  for (int i = 0; i < MMdetN; i++)
    for (int j = 0; j < 6; j++)
    {
      if (Npars == 3)
        algPar3[i][j] = 0;
      if (Npars == 6)
        algPar6[i][j] = 0;
    }

  ifstream in_file(filename.c_str());
  // if(0)
  if (in_file)
  {
    cout << " --> Read alignment data from " << filename << endl;
    for (int i = 0; i < MMdetN; i++)
    {
      if (Npars == 3)
      {
        cout << " --> det =" << i;
        for (int j = 0; j < 6; j++)
        {
          in_file >> algPar3[i][j];
          cout << " " << ParName[j] << " " << algPar3[i][j];
        }
        cout << endl;
      }
      if (Npars == 6)
      {
        cout << " --> det=" << i;
        for (int j = 0; j < 6; j++)
        {
          in_file >> algPar6[i][j];
          cout << ", " << ParName[j] << "=" << algPar6[i][j];
        }
        cout << endl;
      }
    }
    cout << " --> Read Tracker alignment pars succesfully !" << endl
         << endl;
    in_file.close();
  }
  else
    Alignment(filename, Npars);
}

void ReconstructTrack::Alignment(string filename, int Npars)
{
  cout << " --> Begin to do tracker alignment !" << endl;

  double res = alg_res;                        // delta-residual resolution
  double rho = alg_rho, epsilon = alg_epsilon; // for alginment loop
  vector<double> xvec[MMdetN], yvec[MMdetN], zvec[MMdetN];
  int nentries = fdecTree[0]->GetEntriesFast();
  for (int i = 0; i < nentries; i++)
  {
    bool is_good = true;
    double posX[MMdetN], posY[MMdetN], posZ[MMdetN];
    for (int j = 0; j < MMdetN; j++)
    {
      fdecTree[j]->GetEntry(i);
      // if (j == 0)
      // {
      //   is_good = is_good && sig && x > 100 && y > 100 && x < 150 && y < 165;
      // }
      // else if (j == 1)
      // {
      //   is_good = is_good && sig && x > 100 && y > 100 && x < 150 && y < 165;
      // }
      // else if (j == 2)
      // {
      //   is_good = is_good && sig && x > 100 && y > 100 && x < 160 && y < 180;
      // }
      // else if (j == 3)
      // {
      //   is_good = is_good && sig && x > 100 && y > 100 && x < 160 && y < 180;
      // }
      // else if (j == 4)
      // {
      //   is_good = is_good && sig && x > 77 && y > 100 && x < 160 && y < 210;
      // }
      // else if (j == 5)
      // {
      //   is_good = is_good && sig && x > 70 && y > 100 && x < 160 && y < 215;
      // }
      is_good = is_good && sig && x > 25 && y > 25 && x < 225 && y < 225;
      posX[j] = x;
      posY[j] = y;
      posZ[j] = z;
    }
    if (is_good)
    {
      for (int j = 0; j < MMdetN; j++)
      {
        xvec[j].push_back(posX[j]);
        yvec[j].push_back(posY[j]);
        zvec[j].push_back(posZ[j]);
      }
    }
  }
  nentries = xvec[0].size();
  cout << " --> Select events = " << nentries << endl;

  double pxres[MMdetN][2], pyres[MMdetN][2];
  for (int i = 0; i < MMdetN; i++)
  {
    pxres[i][0] = -100 * 1000;
    pxres[i][1] = 100 * 1000;
    pyres[i][0] = -100 * 1000;
    pyres[i][1] = 100 * 1000;
  }
  double ppars[MMdetN][6] = {0};

  double pars[MMdetN][6] = {0};
  double D_pre[MMdetN][6] = {0}, Chi_pre = 0;
  double ED[MMdetN][6] = {0}, Ed[MMdetN][6] = {0};
  for (int loop = 0;; loop++)
  {
    // if (loop > 1000)
    // {
    //   alg_fit_R2_cut = 0.99;
    //   alg_theta_cut = 1;
    // }
    // if (loop > 1000)
    // {
    //   alg_fit_R2_cut = 0.999;
    //   alg_fit_res_cut = 3;
    // }
    // if (loop > 2000)
    // {
    //   alg_fit_res_cut = 3;
    //   alg_fit_R2_cut = 0.999;
    //   alg_theta_cut = 1;
    // }
    TH1F *phpos[MMdetN][2];

    for (int i = 0; i < MMdetN; i++)
    {
      phpos[i][0] = new TH1F(Form("presX%d", i), Form("presX%d", i), 500, -4000, 4000);
      phpos[i][0]->GetXaxis()->SetTitle("x [#mum]");
      phpos[i][0]->GetYaxis()->SetTitle("counts");
      phpos[i][1] = new TH1F(Form("presY%d", i), Form("presY%d", i), 500, -4000, 4000);
      phpos[i][1]->GetXaxis()->SetTitle("y [#mum]");
      phpos[i][1]->GetYaxis()->SetTitle("counts");
    }

    cout << " --> loop = " << loop << endl;
    int NChi = 0;
    double Chi2 = 0, D[MMdetN][6] = {0}, d[MMdetN][6] = {0};

    for (int entry = 0; entry < nentries; entry++)
    {
      vector<double> x0, y0, z0, xalg, yalg, zalg;
      vector<double> pxalg, pyalg, pzalg;
      for (int i = 0; i < MMdetN; i++)
      {
        x0.push_back(xvec[i][entry] * 0.4 - centre);
        y0.push_back(yvec[i][entry] * 0.4 - centre);
        z0.push_back(zvec[i][entry]);
        xalg.push_back(x0[i] + pars[i][0] - y0[i] / centre * pars[i][5]);
        yalg.push_back(y0[i] + pars[i][1] + x0[i] / centre * pars[i][5]);
        zalg.push_back(z0[i] + pars[i][2] - x0[i] / centre * pars[i][4] + y0[i] / centre * pars[i][3]);

        pxalg.push_back(x0[i] + ppars[i][0] - y0[i] / centre * ppars[i][5]);
        pyalg.push_back(y0[i] + ppars[i][1] + x0[i] / centre * ppars[i][5]);
        pzalg.push_back(z0[i] + ppars[i][2] - x0[i] / centre * ppars[i][4] + y0[i] / centre * ppars[i][3]);
      }

      bool is_good_event = true;
      // x = K*z+B, R2 = correction index squre
      double Kx[MMdetN] = {0}, Bx[MMdetN] = {0}, Rx[MMdetN] = {0};
      double Ky[MMdetN] = {0}, By[MMdetN] = {0}, Ry[MMdetN] = {0};

      double pKx[MMdetN] = {0}, pBx[MMdetN] = {0}, pRx[MMdetN] = {0};
      double pKy[MMdetN] = {0}, pBy[MMdetN] = {0}, pRy[MMdetN] = {0};

      // mean
      double Mx[MMdetN] = {0}, My[MMdetN] = {0}, Mz[MMdetN] = {0};
      double Mxx[MMdetN] = {0}, Myy[MMdetN] = {0}, Mzz[MMdetN] = {0};
      double Mzx[MMdetN] = {0}, Mzy[MMdetN] = {0};

      double pMx[MMdetN] = {0}, pMy[MMdetN] = {0}, pMz[MMdetN] = {0};
      double pMxx[MMdetN] = {0}, pMyy[MMdetN] = {0}, pMzz[MMdetN] = {0};
      double pMzx[MMdetN] = {0}, pMzy[MMdetN] = {0};

      int N = MMdetN - 1;

      for (int i = 0; i < MMdetN; i++)
      {
        for (int j = 0; j < MMdetN; j++)
        {
          if (i != j)
          {
            Mx[i] += xalg[j] / N;
            My[i] += yalg[j] / N;
            Mz[i] += zalg[j] / N;
            Mxx[i] += xalg[j] * xalg[j] / N;
            Myy[i] += yalg[j] * yalg[j] / N;
            Mzz[i] += zalg[j] * zalg[j] / N;
            Mzx[i] += zalg[j] * xalg[j] / N;
            Mzy[i] += zalg[j] * yalg[j] / N;
          }
        }
        Kx[i] = (Mzx[i] - Mz[i] * Mx[i]) / (Mzz[i] - Mz[i] * Mz[i]);
        Bx[i] = Mx[i] - Kx[i] * Mz[i];
        Rx[i] = (Mzx[i] - Mz[i] * Mx[i]) * (Mzx[i] - Mz[i] * Mx[i]) / (Mxx[i] - Mx[i] * Mx[i]) / (Mzz[i] - Mz[i] * Mz[i]);
        Ky[i] = (Mzy[i] - Mz[i] * My[i]) / (Mzz[i] - Mz[i] * Mz[i]);
        By[i] = My[i] - Ky[i] * Mz[i];
        Ry[i] = (Mzy[i] - Mz[i] * My[i]) * (Mzy[i] - Mz[i] * My[i]) / (Myy[i] - My[i] * My[i]) / (Mzz[i] - Mz[i] * Mz[i]);

        for (int j = 0; j < MMdetN; j++)
        {
          if (i != j)
          {
            pMx[i] += pxalg[j] / N;
            pMy[i] += pyalg[j] / N;
            pMz[i] += pzalg[j] / N;
            pMxx[i] += pxalg[j] * pxalg[j] / N;
            pMyy[i] += pyalg[j] * pyalg[j] / N;
            pMzz[i] += pzalg[j] * pzalg[j] / N;
            pMzx[i] += pzalg[j] * pxalg[j] / N;
            pMzy[i] += pzalg[j] * pyalg[j] / N;
          }
        }
        pKx[i] = (pMzx[i] - pMz[i] * pMx[i]) / (pMzz[i] - pMz[i] * pMz[i]);
        pBx[i] = pMx[i] - pKx[i] * pMz[i];
        pKy[i] = (pMzy[i] - pMz[i] * pMy[i]) / (pMzz[i] - pMz[i] * pMz[i]);
        pBy[i] = pMy[i] - pKy[i] * pMz[i];

        double pdx, pdy;
        pdx = pKx[i] * pzalg[i] + pBx[i] - pxalg[i];
        pdy = pKy[i] * pzalg[i] + pBy[i] - pyalg[i];
        phpos[i][0]->Fill(pdx * 1000);
        phpos[i][1]->Fill(pdy * 1000);

        is_good_event = is_good_event && Rx[i] > alg_fit_R2_cut && Ry[i] > alg_fit_R2_cut;
        is_good_event = is_good_event && Abs(ATan(Kx[i]) / Pi() * 180) < alg_theta_cut;
        is_good_event = is_good_event && Abs(ATan(Ky[i]) / Pi() * 180) < alg_theta_cut;
        is_good_event = is_good_event && Abs(Kx[i] * zalg[i] + Bx[i] - xalg[i]) < alg_fit_res_cut;
        is_good_event = is_good_event && Abs(Ky[i] * zalg[i] + By[i] - yalg[i]) < alg_fit_res_cut;
        // is_good_event = is_good_event && (pdx * 1000 > pxres[i][0] && pdx * 1000 < pxres[i][1]);
        // is_good_event = is_good_event && (pdy * 1000 > pyres[i][0] && pdy * 1000 < pyres[i][1]);
      }
      if (!is_good_event)
        continue;

      for (int i = 0; i < MMdetN; i++)
      {
        double dx, dy; // residual
        double DKxz, DKxx, DBxz, DBxx;
        double DKyz, DKyy, DByz, DByy;

        double w = 1;
        // if (i == 0)
        // {
        //   w = 1 / (0.450 * 0.450);
        // }
        // else if (i == 1)
        // {
        //   w = 0;
        // }
        // else if (i == 2)
        // {
        //   w = 1 / (0.450 * 0.450);
        // }
        // else if (i == 3)
        // {
        //   w = 1 / (0.350 * 0.350);
        // }
        // else if (i == 4)
        // {
        //   w = 1 / (0.300 * 0.300);
        // }
        // else if (i == 5)
        // {
        //   w = 1 / (0.450 * 0.450);
        // }

        NChi += 2;
        if (1)
        {
          dx = Kx[i] * zalg[i] + Bx[i] - xalg[i];
          dy = Ky[i] * zalg[i] + By[i] - yalg[i];
          Chi2 += w * Power(dx, 2) + Power(dy, 2);

          for (int j = 0; j < MMdetN; j++)
          {
            if (i != j)
            {
              DKxz = ((xalg[j] - Mx[i]) / (Mzz[i] - Mz[i] * Mz[i]) - 2 * Kx[i] * (zalg[j] - Mz[i]) / (Mzz[i] - Mz[i] * Mz[i])) / N;
              DKxx = (zalg[j] - Mz[i]) / (Mzz[i] - Mz[i] * Mz[i]) / N;
              DBxz = -Kx[i] / N - DKxz * Mz[i];
              DBxx = 1. / N - DKxx * Mz[i];
              DKyz = ((yalg[j] - My[i]) / (Mzz[i] - Mz[i] * Mz[i]) - 2 * Ky[i] * (zalg[j] - Mz[i]) / (Mzz[i] - Mz[i] * Mz[i])) / N;
              DKyy = (zalg[j] - Mz[i]) / (Mzz[i] - Mz[i] * Mz[i]) / N;
              DByz = -Ky[i] / N - DKyz * Mz[i];
              DByy = 1. / N - DKyy * Mz[i];

              D[j][0] += w * (2 * dx * (zalg[i] * DKxx + DBxx));
              D[j][1] += w * (2 * dy * (zalg[i] * DKyy + DByy));
              D[j][2] += w * (2 * dx * (zalg[i] * DKxz + DBxz) + 2 * dy * (zalg[i] * DKyz + DByz));
              D[j][3] += w * (yalg[j] / centre * (2 * dx * (zalg[i] * DKxz + DBxz) + 2 * dy * (zalg[i] * DKyz + DByz)));
              D[j][4] += w * (-xalg[j] / centre * (2 * dx * (zalg[i] * DKxz + DBxz) + 2 * dy * (zalg[i] * DKyz + DByz)));
              D[j][5] += w * (-yalg[j] / centre * 2 * dx * (zalg[i] * DKxx + DBxx) + xalg[j] / centre * 2 * dy * (zalg[i] * DKyy + DByy));
            }
            else
            {
              D[j][0] += w * (-2 * dx);
              D[j][1] += w * (-2 * dy);
              D[j][2] += w * (2 * dx * Kx[i] + 2 * dy * Ky[i]);
              D[j][3] += w * (yalg[j] / centre * (2 * dx * Kx[i] + 2 * dy * Ky[i]));
              D[j][4] += w * (-xalg[j] / centre * (2 * dx * Kx[i] + 2 * dy * Ky[i]));
              D[j][5] += w * (yalg[j] / centre * 2 * dx - xalg[j] / centre * 2 * dy);
            }
          }
        }
      }

      if (0)
      {
        const double angle_weight = 1000 * 1000 * 10;
        int nA = 2;
        double sumzA = 0, sumxA = 0, sumyA = 0, sumzzA = 0, sumzxA = 0, sumzyA = 0;
        for (int i = 0; i < 2; ++i)
        {
          sumzA += zalg[i];
          sumxA += xalg[i];
          sumyA += yalg[i];
          sumzzA += zalg[i] * zalg[i];
          sumzxA += zalg[i] * xalg[i];
          sumzyA += zalg[i] * yalg[i];
        }
        double DA = nA * sumzzA - sumzA * sumzA;
        double NxA = nA * sumzxA - sumzA * sumxA;
        double NyA = nA * sumzyA - sumzA * sumyA;
        double kxA = (DA != 0) ? NxA / DA : 0;
        double kyA = (DA != 0) ? NyA / DA : 0;

        // 计算组 B（探测器2~5）的斜率
        int nB = 4;
        double sumzB = 0, sumxB = 0, sumyB = 0, sumzzB = 0, sumzxB = 0, sumzyB = 0;
        for (int i = 2; i < MMdetN; ++i)
        {
          sumzB += zalg[i];
          sumxB += xalg[i];
          sumyB += yalg[i];
          sumzzB += zalg[i] * zalg[i];
          sumzxB += zalg[i] * xalg[i];
          sumzyB += zalg[i] * yalg[i];
        }
        double DB = nB * sumzzB - sumzB * sumzB;
        double NxB = nB * sumzxB - sumzB * sumxB;
        double NyB = nB * sumzyB - sumzB * sumyB;
        double kxB = (DB != 0) ? NxB / DB : 0;
        double kyB = (DB != 0) ? NyB / DB : 0;

        double delta_kx = kxA - kxB;
        double delta_ky = kyA - kyB;

        Chi2 += angle_weight * (delta_kx * delta_kx + delta_ky * delta_ky);

        double dkxA_dx[2], dkxA_dz[2], dkyA_dy[2], dkyA_dz[2];
        double dkxB_dx[4], dkxB_dz[4], dkyB_dy[4], dkyB_dz[4];

        for (int i = 0; i < 2; ++i)
        {
          int idx = i;
          double zi = zalg[idx];
          double xi = xalg[idx];
          double yi = yalg[idx];
          if (DA != 0)
          {
            dkxA_dx[i] = (nA * zi - sumzA) / DA;
            dkxA_dz[i] = ((nA * xi - sumxA) - 2 * kxA * (nA * zi - sumzA)) / DA;
            dkyA_dy[i] = (nA * zi - sumzA) / DA;
            dkyA_dz[i] = ((nA * yi - sumyA) - 2 * kyA * (nA * zi - sumzA)) / DA;
          }
          else
          {
            dkxA_dx[i] = dkxA_dz[i] = dkyA_dy[i] = dkyA_dz[i] = 0;
          }
        }

        for (int i = 2; i < MMdetN; ++i)
        {
          int idx = i;
          int j = i - 2;
          double zi = zalg[idx];
          double xi = xalg[idx];
          double yi = yalg[idx];
          if (DB != 0)
          {
            dkxB_dx[j] = (nB * zi - sumzB) / DB;
            dkxB_dz[j] = ((nB * xi - sumxB) - 2 * kxB * (nB * zi - sumzB)) / DB;
            dkyB_dy[j] = (nB * zi - sumzB) / DB;
            dkyB_dz[j] = ((nB * yi - sumyB) - 2 * kyB * (nB * zi - sumzB)) / DB;
          }
          else
          {
            dkxB_dx[j] = dkxB_dz[j] = dkyB_dy[j] = dkyB_dz[j] = 0;
          }
        }

        for (int j = 0; j < MMdetN; ++j)
        {
          // 原始坐标（用于求坐标对参数的导数）
          double x0j = x0[j];
          double y0j = y0[j];
          // 坐标对参数的导数（只列出非零项）
          double dx_dpar[6] = {0};
          double dy_dpar[6] = {0};
          double dz_dpar[6] = {0};
          dx_dpar[0] = 1.0;
          dx_dpar[5] = -y0j / centre;
          dy_dpar[1] = 1.0;
          dy_dpar[5] = +x0j / centre;
          dz_dpar[2] = 1.0;
          dz_dpar[3] = +y0j / centre;
          dz_dpar[4] = -x0j / centre;

          int group = (j < 2) ? 0 : 1;
          int loc = (j < 2) ? j : (j - 2);

          double dkx_dpar = 0, dky_dpar = 0;
          if (group == 0)
          { // 组A
            for (int k = 0; k < 6; ++k)
            {
              dkx_dpar += dkxA_dx[loc] * dx_dpar[k] + dkxA_dz[loc] * dz_dpar[k];
              dky_dpar += dkyA_dy[loc] * dy_dpar[k] + dkyA_dz[loc] * dz_dpar[k];
            }
          }
          else
          { // 组B
            for (int k = 0; k < 6; ++k)
            {
              dkx_dpar -= (dkxB_dx[loc] * dx_dpar[k] + dkxB_dz[loc] * dz_dpar[k]);
              dky_dpar -= (dkyB_dy[loc] * dy_dpar[k] + dkyB_dz[loc] * dz_dpar[k]);
            }
          }

          double factor = 2.0 * angle_weight;
          double dchi2_dpar = factor * (delta_kx * dkx_dpar + delta_ky * dky_dpar);

          double c_D[6];
          if (j < 2)
          {
            for (int k = 0; k < 6; k++)
            {
              c_D[k] = 2 * (delta_kx * (dkxA_dx[j] * dx_dpar[k] + dkxA_dz[j] * dz_dpar[k]) + delta_ky * (dkyA_dy[j] * dy_dpar[k] + dkyA_dz[j] * dz_dpar[k]));
            }
          }
          else
          {
            for (int k = 0; k < 6; k++)
            {
              c_D[k] = -2 * (delta_kx * (dkxB_dx[j - 2] * dx_dpar[k] + dkxB_dz[j - 2] * dz_dpar[k]) + delta_ky * (dkyB_dy[j - 2] * dy_dpar[k] + dkyB_dz[j - 2] * dz_dpar[k]));
            }
          }

          // 只累加前 Npars 个参数（3 或 6）
          for (int k = 0; k < Npars; ++k)
          {
            D[j][k] += c_D[k] * angle_weight;
          }
        }
      }

    } // entry loop

    cout << " --> Event number = " << NChi / (MMdetN) / 2 << ", Chi = " << Sqrt(Chi2 / NChi)
         << ", DeltaChi = " << Chi_pre - Sqrt(Chi2 / NChi) << endl;

    for (int i = 0; i < MMdetN; i++)
    {
      TF1 *alignXFit = FitGaussianSafe(phpos[i][0], -10000, 10000, Form("alignX%d", i));
      if (!alignXFit || Abs(alignXFit->GetParameter(1)) > 10000)
      {
        pxres[i][0] = -10000;
        pxres[i][1] = 10000;
      }
      else
      {
        pxres[i][0] = alignXFit->GetParameter(1) - 5 * alignXFit->GetParameter(2);
        pxres[i][1] = alignXFit->GetParameter(1) + 5 * alignXFit->GetParameter(2);
      }

      TF1 *alignYFit = FitGaussianSafe(phpos[i][1], -10000, 10000, Form("alignY%d", i));
      if (!alignYFit || Abs(alignYFit->GetParameter(1)) > 10000)
      {
        pyres[i][0] = -10000;
        pyres[i][1] = 10000;
      }
      else
      {
        pyres[i][0] = alignYFit->GetParameter(1) - 5 * alignYFit->GetParameter(2);
        pyres[i][1] = alignYFit->GetParameter(1) + 5 * alignYFit->GetParameter(2);
      }

      delete phpos[i][0];
      delete phpos[i][1];
      cout << pxres[i][0] << "," << pxres[i][1] << "," << pyres[i][0] << "," << pyres[i][1] << endl;
    }

    for (int i = 0; i < MMdetN; i++)
    {
      for (int j = 0; j < Npars; j++)
      {
        ppars[i][j] = pars[i][j];
      }
    }

    for (int i = 0; i < MMdetN; i++)
    {
      for (int j = 0; j < Npars; j++)
      {
        D[i][j] = D[i][j] / NChi * MMdetN * 2;
        if (D_pre[i][j] * D[i][j] < 0)
          Ed[i][j] = rho * Ed[i][j];
        ED[i][j] = rho * ED[i][j] + (1 - rho) * Power(D[i][j], 2);
        d[i][j] = -Sqrt(Ed[i][j] + epsilon) / Sqrt(ED[i][j] + epsilon) * D[i][j];
        Ed[i][j] = rho * Ed[i][j] + (1 - rho) * Power(d[i][j], 2);
        pars[i][j] += d[i][j];
        D_pre[i][j] = D[i][j];
        epsilon = epsilon * rho;
      }
      cout << " --> det=" << i;
      for (int j = 0; j < 6; j++)
        cout << ", " << ParName[j] << "=" << pars[i][j];
      cout << endl;
    }
    cout << endl;

    if (Abs(Chi_pre - Sqrt(Chi2 / NChi)) < res || loop > 5e3)
      break;
    else
      Chi_pre = Sqrt(Chi2 / NChi);
  }

  ofstream out(filename.c_str());
  for (int i = 0; i < MMdetN; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      if (j >= 3)
        pars[i][j] = pars[i][j] / centre;
      if (Npars == 3)
        algPar3[i][j] = pars[i][j];
      if (Npars == 6)
        algPar6[i][j] = pars[i][j];
      out << pars[i][j] << " ";
    }
    out << endl;
  }
  out.close();
  cout << " --> Alignment succesfully !" << endl;
}

void ReconstructTrack::ReadSecondAlignmentPars(string filename, int Npars)
{
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 6; j++)
    {
      if (Npars == 3)
        SecondalgPar3[i][j] = 0;
      if (Npars == 6)
        SecondalgPar6[i][j] = 0;
    }

  ifstream in_file(filename.c_str());
  // if(0)
  if (in_file)
  {
    cout << " --> Read alignment data from " << filename << endl;
    for (int i = 0; i < 2; i++)
    {
      if (Npars == 3)
      {
        cout << " --> det =" << i;
        for (int j = 0; j < 6; j++)
        {
          in_file >> SecondalgPar3[i][j];
          cout << " " << ParName[j] << " " << SecondalgPar3[i][j];
        }
        cout << endl;
      }
      if (Npars == 6)
      {
        cout << " --> det=" << i;
        for (int j = 0; j < 6; j++)
        {
          in_file >> SecondalgPar6[i][j];
          cout << ", " << ParName[j] << "=" << SecondalgPar6[i][j];
        }
        cout << endl;
      }
    }
    cout << " --> Read Tracker alignment pars succesfully !" << endl
         << endl;
    in_file.close();
  }
  else
    SecondAlignment(filename, Npars);
}

void ReconstructTrack::SecondAlignment(string filename, int Npars)
{
  cout << " --> Begin to do Second tracker alignment !" << endl;
  double res = alg_res;                        // delta-residual resolution
  double rho = alg_rho, epsilon = alg_epsilon; // for alginment loop
  vector<double> xvec[MMdetN], yvec[MMdetN], zvec[MMdetN];
  int nentries = fdecTree[0]->GetEntriesFast();
  for (int i = 0; i < nentries / 10; i++)
  {
    bool is_good = true;
    double posX[MMdetN], posY[MMdetN], posZ[MMdetN];
    for (int j = 0; j < MMdetN; j++)
    {
      fdecTree[j]->GetEntry(i);
      if (j == 0)
      {
        is_good = is_good && sig && x > 100 && y > 100 && x < 150 && y < 165;
      }
      else if (j == 1)
      {
        is_good = is_good && sig && x > 100 && y > 100 && x < 150 && y < 165;
      }
      else if (j == 2)
      {
        is_good = is_good && sig && x > 100 && y > 100 && x < 160 && y < 180;
      }
      else if (j == 3)
      {
        is_good = is_good && sig && x > 100 && y > 100 && x < 160 && y < 180;
      }
      else if (j == 4)
      {
        is_good = is_good && sig && x > 77 && y > 100 && x < 160 && y < 210;
      }
      else if (j == 5)
      {
        is_good = is_good && sig && x > 70 && y > 100 && x < 160 && y < 215;
      }
      posX[j] = x * 0.4 - centre;
      posY[j] = y * 0.4 - centre;
      posZ[j] = z;
    }
    if (is_good)
    {
      for (int j = 0; j < MMdetN; j++)
      {
        xvec[j].push_back(posX[j] + algPar6[j][0] - posY[j] / centre * algPar6[j][5]);
        yvec[j].push_back(posY[j] + algPar6[j][1] + posX[j] / centre * algPar6[j][5]);
        zvec[j].push_back(posZ[j] + algPar6[j][2] - posX[j] / centre * algPar6[j][4] + posY[j] / centre * algPar6[j][3]);
      }
    }
  }
  nentries = xvec[0].size();
  cout << " --> Select events = " << nentries << endl;

  double pars[2][6] = {0};
  double D_pre[2][6] = {0}, Chi_pre = 0;
  double ED[2][6] = {0}, Ed[2][6] = {0};

  for (int loop = 0;; loop++)
  {
    cout << " --> loop = " << loop << endl;
    int NChi = 0;
    double Chi2 = 0, D[2][6] = {0}, d[2][6] = {0};
    for (int entry = 0; entry < nentries; entry++)
    {
      vector<double> x0, y0, z0, xalg, yalg, zalg;
      for (int i = 0; i < 2; i++)
      {
        x0.push_back(xvec[i][entry]);
        y0.push_back(yvec[i][entry]);
        z0.push_back(zvec[i][entry]);
        xalg.push_back(x0[i] + pars[0][0] - y0[i] / centre * pars[0][5]);
        yalg.push_back(y0[i] + pars[0][1] + x0[i] / centre * pars[0][5]);
        zalg.push_back(z0[i] + pars[0][2] - x0[i] / centre * pars[0][4] + y0[i] / centre * pars[0][3]);
      }

      for (int i = 2; i < 6; i++)
      {
        x0.push_back(xvec[i][entry]);
        y0.push_back(yvec[i][entry]);
        z0.push_back(zvec[i][entry]);
        xalg.push_back(x0[i] + pars[1][0] - y0[i] / centre * pars[1][5]);
        yalg.push_back(y0[i] + pars[1][1] + x0[i] / centre * pars[1][5]);
        zalg.push_back(z0[i] + pars[1][2] - x0[i] / centre * pars[1][4] + y0[i] / centre * pars[1][3]);
      }

      bool is_good_event = true;
      // x = K*z+B, R2 = correction index squre
      double Kx[MMdetN] = {0}, Bx[MMdetN] = {0}, Rx[MMdetN] = {0};
      double Ky[MMdetN] = {0}, By[MMdetN] = {0}, Ry[MMdetN] = {0};
      // mean
      double Mx[MMdetN] = {0}, My[MMdetN] = {0}, Mz[MMdetN] = {0};
      double Mxx[MMdetN] = {0}, Myy[MMdetN] = {0}, Mzz[MMdetN] = {0};
      double Mzx[MMdetN] = {0}, Mzy[MMdetN] = {0};
      int N = MMdetN - 1;

      for (int i = 0; i < MMdetN; i++)
      {
        for (int j = 0; j < MMdetN; j++)
        {
          if (i != j)
          {
            Mx[i] += xalg[j] / N;
            My[i] += yalg[j] / N;
            Mz[i] += zalg[j] / N;
            Mxx[i] += xalg[j] * xalg[j] / N;
            Myy[i] += yalg[j] * yalg[j] / N;
            Mzz[i] += zalg[j] * zalg[j] / N;
            Mzx[i] += zalg[j] * xalg[j] / N;
            Mzy[i] += zalg[j] * yalg[j] / N;
          }
        }

        Kx[i] = (Mzx[i] - Mz[i] * Mx[i]) / (Mzz[i] - Mz[i] * Mz[i]);
        Bx[i] = Mx[i] - Kx[i] * Mz[i];
        Rx[i] = (Mzx[i] - Mz[i] * Mx[i]) * (Mzx[i] - Mz[i] * Mx[i]) / (Mxx[i] - Mx[i] * Mx[i]) / (Mzz[i] - Mz[i] * Mz[i]);
        Ky[i] = (Mzy[i] - Mz[i] * My[i]) / (Mzz[i] - Mz[i] * Mz[i]);
        By[i] = My[i] - Ky[i] * Mz[i];
        Ry[i] = (Mzy[i] - Mz[i] * My[i]) * (Mzy[i] - Mz[i] * My[i]) / (Myy[i] - My[i] * My[i]) / (Mzz[i] - Mz[i] * Mz[i]);

        is_good_event = is_good_event && Rx[i] > alg_fit_R2_cut && Ry[i] > alg_fit_R2_cut;
        is_good_event = is_good_event && Abs(ATan(Kx[i]) / Pi() * 180) < alg_theta_cut;
        is_good_event = is_good_event && Abs(ATan(Ky[i]) / Pi() * 180) < alg_theta_cut;
        is_good_event = is_good_event && Abs(Kx[i] * zalg[i] + Bx[i] - xalg[i]) < alg_fit_res_cut;
        is_good_event = is_good_event && Abs(Ky[i] * zalg[i] + By[i] - yalg[i]) < alg_fit_res_cut;
      }
      if (!is_good_event)
        continue;

      for (int i = 0; i < MMdetN; i++)
      {
        double dx, dy; // residual
        double DKxz, DKxx, DBxz, DBxx;
        double DKyz, DKyy, DByz, DByy;

        NChi += 2;
        dx = Kx[i] * zalg[i] + Bx[i] - xalg[i];
        dy = Ky[i] * zalg[i] + By[i] - yalg[i];
        Chi2 += Power(dx, 2) + Power(dy, 2);

        for (int j = 0; j < MMdetN; j++)
        {
          int j_group = (j < 2) ? 0 : 1;
          if (j_group == 0)
          {
            if (i != j)
            {
              D[0][0] += 2 * dx * (zalg[i] * DKxx + DBxx);
              D[0][1] += 2 * dy * (zalg[i] * DKyy + DByy);
              D[0][2] += 2 * dx * (zalg[i] * DKxz + DBxz) + 2 * dy * (zalg[i] * DKyz + DByz);
              D[0][3] += yalg[j] / centre * (2 * dx * (zalg[i] * DKxz + DBxz) + 2 * dy * (zalg[i] * DKyz + DByz));
              D[0][4] += -xalg[j] / centre * (2 * dx * (zalg[i] * DKxz + DBxz) + 2 * dy * (zalg[i] * DKyz + DByz));
              D[0][5] += -yalg[j] / centre * 2 * dx * (zalg[i] * DKxx + DBxx) + xalg[j] / centre * 2 * dy * (zalg[i] * DKyy + DByy);
            }
            else
            {
              D[0][0] += -2 * dx;
              D[0][1] += -2 * dy;
              D[0][2] += 2 * dx * Kx[i] + 2 * dy * Ky[i];
              D[0][3] += yalg[j] / centre * (2 * dx * Kx[i] + 2 * dy * Ky[i]);
              D[0][4] += -xalg[j] / centre * (2 * dx * Kx[i] + 2 * dy * Ky[i]);
              D[0][5] += yalg[j] / centre * 2 * dx - xalg[j] / centre * 2 * dy;
            }
          }
          else
          { // j_group == 1
            // 类似地，累加到 D[1]
            if (i != j)
            {
              D[1][0] += 2 * dx * (zalg[i] * DKxx + DBxx);
              D[1][1] += 2 * dy * (zalg[i] * DKyy + DByy);
              D[1][2] += 2 * dx * (zalg[i] * DKxz + DBxz) + 2 * dy * (zalg[i] * DKyz + DByz);
              D[1][3] += yalg[j] / centre * (2 * dx * (zalg[i] * DKxz + DBxz) + 2 * dy * (zalg[i] * DKyz + DByz));
              D[1][4] += -xalg[j] / centre * (2 * dx * (zalg[i] * DKxz + DBxz) + 2 * dy * (zalg[i] * DKyz + DByz));
              D[1][5] += -yalg[j] / centre * 2 * dx * (zalg[i] * DKxx + DBxx) + xalg[j] / centre * 2 * dy * (zalg[i] * DKyy + DByy);
            }
            else
            {
              D[1][0] += -2 * dx;
              D[1][1] += -2 * dy;
              D[1][2] += 2 * dx * Kx[i] + 2 * dy * Ky[i];
              D[1][3] += yalg[j] / centre * (2 * dx * Kx[i] + 2 * dy * Ky[i]);
              D[1][4] += -xalg[j] / centre * (2 * dx * Kx[i] + 2 * dy * Ky[i]);
              D[1][5] += yalg[j] / centre * 2 * dx - xalg[j] / centre * 2 * dy;
            }
          }
        }
      }
    }

    cout << " --> Event number = " << NChi / MMdetN / 2 << ", Chi = " << Sqrt(Chi2 / NChi)
         << ", DeltaChi = " << Chi_pre - Sqrt(Chi2 / NChi) << endl;

    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < Npars; j++)
      {
        D[i][j] = D[i][j] / NChi * MMdetN * 2;
        if (D_pre[i][j] * D[i][j] < 0)
          Ed[i][j] = rho * Ed[i][j];
        ED[i][j] = rho * ED[i][j] + (1 - rho) * Power(D[i][j], 2);
        d[i][j] = -Sqrt(Ed[i][j] + epsilon) / Sqrt(ED[i][j] + epsilon) * D[i][j];
        Ed[i][j] = rho * Ed[i][j] + (1 - rho) * Power(d[i][j], 2);
        pars[i][j] += d[i][j];
        D_pre[i][j] = D[i][j];
        epsilon = epsilon * rho;
      }
      cout << " --> det=" << i;
      for (int j = 0; j < 6; j++)
        cout << ", " << ParName[j] << "=" << pars[i][j];
      cout << endl;
    }
    cout << endl;

    if (Abs(Chi_pre - Sqrt(Chi2 / NChi)) < res || loop > 5e3)
      break;
    else
      Chi_pre = Sqrt(Chi2 / NChi);
  }

  ofstream out(filename.c_str());
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      if (j >= 3)
        pars[i][j] = pars[i][j] / centre;
      if (Npars == 3)
        SecondalgPar3[i][j] = pars[i][j];
      if (Npars == 6)
        SecondalgPar6[i][j] = pars[i][j];
      out << pars[i][j] << " ";
    }
    out << endl;
  }
  out.close();
  cout << " --> Alignment succesfully !" << endl;
}

void ReconstructTrack::Millepede_ali()
{
  define_align_matrix();
  double chi2 = 0, chi2_pre = 0;

  for (int loop = 0; loop < 100; loop++)
  {
    cout << "loop" << loop << endl;

    chi2 = get_data_chi2();
    cout << "Chi2 = " << chi2 << ", DeltaChi2 = " << chi2_pre - chi2 << endl;
    double dchi = chi2_pre - chi2;
    chi2_pre = chi2;

    print_alignment_list();
    Write_mille_file();

    system("/mnt/d/Millepede/target/pede ./pedestr.txt > log 2>1");
    read_millepede_res("./millepede.res");
  }
}

void ReconstructTrack::Write_mille_file()
{
  define_align_matrix();

  Mille *millefile = new Mille("millebin.bin");

  int nentries = fdecTree[0]->GetEntriesFast();
  for (int i = 0; i < nentries / 20; i++)
  {
    vector<double> zx_origin, x_origin, zy_origin, y_origin;
    vector<double> a_zx, a_x, a_zy, a_y;
    vector<int> layer;
    for (int j = 0; j < MMdetN; j++)
    {
      fdecTree[j]->GetEntry(i);
      if (sig)
      {
        x_origin.push_back(x * 0.4 - centre);
        zx_origin.push_back(z);
        y_origin.push_back(y * 0.4 - centre);
        zy_origin.push_back(z);

        double xx = x * 0.4 - centre;
        double yy = y * 0.4 - centre;
        double zz = z;
        local_to_global(xx, yy, zz, j);
        a_x.push_back(xx);
        a_zx.push_back(zz);
        a_y.push_back(yy);
        a_zy.push_back(zz);

        layer.push_back(j);
      }
    }

    if (layer.size() != MMdetN)
    {
      continue;
    }

    double kx, bx, ky, by;
    double chi2Dndf_x, chi2Dndf_y;

    double sigmax = 0.4 / sqrt(12);
    double sigmay = 0.4 / sqrt(12);

    // 直线拟合
    if (!least_square_fit(a_zx, a_x, sigmax, kx, bx, chi2Dndf_x) || !least_square_fit(a_zy, a_y, sigmay, ky, by, chi2Dndf_y))
    {
      continue;
    }
    define_track_matrix(kx, ky, bx, by);

    double weight = sqrt((chi2Dndf_x + chi2Dndf_y) / 2) * 1;
    if (weight > 1000)
    {
      continue;
    }

    float sigma;
    float rmeas;
    int nlc = 4;
    float derlcx[4] = {0};
    float derlcy[4] = {0};
    int ngl = 6;
    float derglx[6] = {0};
    float dergly[6] = {0};
    int label[6] = {0};

    for (size_t idx = 0; idx < layer.size(); idx++)
    {
      int k = layer[idx];
      for (int j = 0; j < 6; j++)
      {
        label[j] = (k + 1) * 100 + j;
      }

      // zx
      alignment_derlcx(derlcx, k);
      alignment_derglx(derglx, k);
      sigma = sigmax * weight;
      rmeas = (x_origin[idx] - expect_position_x(k));
      millefile->mille(nlc, derlcx, ngl, derglx, label, rmeas, sigma);

      alignment_derlcy(derlcy, k);
      alignment_dergly(dergly, k);
      sigma = sigmay * weight;
      rmeas = (y_origin[idx] - expect_position_y(k));
      millefile->mille(nlc, derlcy, ngl, dergly, label, rmeas, sigma);
    }
    millefile->end();

    // entry
  }

  delete millefile;
}

void ReconstructTrack::define_align_matrix()
{
  // Δr, ΔrT
  for (int i = 0; i < MMdetN; i++)
  {
    double Deltar_arr[3] = {alignment_list[i][0], alignment_list[i][1], alignment_list[i][2]};
    Deltar[i].ResizeTo(3, 1);
    Deltar[i] = TMatrixD(3, 1, Deltar_arr);
    DeltarT[i].ResizeTo(3, 1);
    DeltarT[i] = Deltar[i];
    DeltarT[i].T();
  }

  // ΔRx, ΔRxT
  for (int i = 0; i < MMdetN; i++)
  {
    double DeltaRx_arr[9] = {1,
                             0,
                             0,
                             0,
                             cos(alignment_list[i][3]),
                             -sin(alignment_list[i][3]),
                             0,
                             sin(alignment_list[i][3]),
                             cos(alignment_list[i][3])};
    DeltaRx[i].ResizeTo(3, 3);
    DeltaRx[i] = TMatrixD(3, 3, DeltaRx_arr);
    DeltaRxT[i].ResizeTo(3, 3);
    DeltaRxT[i] = DeltaRx[i];
    DeltaRxT[i].T();
  }

  // ΔRy, ΔRyT
  for (int i = 0; i < MMdetN; i++)
  {
    double DeltaRy_arr[9] = {cos(alignment_list[i][4]), 0, sin(alignment_list[i][4]), 0, 1, 0,
                             -sin(alignment_list[i][4]), 0, cos(alignment_list[i][4])};
    DeltaRy[i].ResizeTo(3, 3);
    DeltaRy[i] = TMatrixD(3, 3, DeltaRy_arr);
    DeltaRyT[i].ResizeTo(3, 3);
    DeltaRyT[i] = DeltaRy[i];
    DeltaRyT[i].T();
  }

  // ΔRz, ΔRzT
  for (int i = 0; i < MMdetN; i++)
  {
    double DeltaRz_arr[9] = {cos(alignment_list[i][5]),
                             -sin(alignment_list[i][5]),
                             0,
                             sin(alignment_list[i][5]),
                             cos(alignment_list[i][5]),
                             0,
                             0,
                             0,
                             1};
    DeltaRz[i].ResizeTo(3, 3);
    DeltaRz[i] = TMatrixD(3, 3, DeltaRz_arr);
    DeltaRzT[i].ResizeTo(3, 3);
    DeltaRzT[i] = DeltaRz[i];
    DeltaRzT[i].T();
  }

  // ΔR, ΔRT
  for (int i = 0; i < MMdetN; i++)
  {
    DeltaR[i].ResizeTo(3, 3);
    DeltaR[i] = DeltaRx[i] * DeltaRy[i] * DeltaRz[i];
    DeltaRT[i].ResizeTo(3, 3);
    DeltaRT[i] = DeltaRzT[i] * DeltaRyT[i] * DeltaRxT[i];
  }

  // DRx, DRxT
  for (int i = 0; i < MMdetN; i++)
  {
    double DRx_arr[9] = {0,
                         0,
                         0,
                         0,
                         -sin(alignment_list[i][3]),
                         -cos(alignment_list[i][3]),
                         0,
                         cos(alignment_list[i][3]),
                         -sin(alignment_list[i][3])};
    DRx[i].ResizeTo(3, 3);
    DRx[i] = TMatrixD(3, 3, DRx_arr);
    DRxT[i].ResizeTo(3, 3);
    DRxT[i] = DRx[i];
    DRxT[i].T();
  }

  // DRy, DRyT
  for (int i = 0; i < MMdetN; i++)
  {
    double DRy_arr[9] = {-sin(alignment_list[i][4]), 0, cos(alignment_list[i][4]), 0, 0, 0,
                         -cos(alignment_list[i][4]), 0, -sin(alignment_list[i][4])};
    DRy[i].ResizeTo(3, 3);
    DRy[i] = TMatrixD(3, 3, DRy_arr);
    DRyT[i].ResizeTo(3, 3);
    DRyT[i] = DRy[i];
    DRyT[i].T();
  }

  // DRz, DRzT
  for (int i = 0; i < MMdetN; i++)
  {
    double DRz_arr[9] = {-sin(alignment_list[i][5]),
                         -cos(alignment_list[i][5]),
                         0,
                         cos(alignment_list[i][5]),
                         -sin(alignment_list[i][5]),
                         0,
                         0,
                         0,
                         0};
    DRz[i].ResizeTo(3, 3);
    DRz[i] = TMatrixD(3, 3, DRz_arr);
    DRzT[i].ResizeTo(3, 3);
    DRzT[i] = DRz[i];
    DRzT[i].T();
  }

  // DR, DRT
  for (int i = 0; i < MMdetN; i++)
  {
    DR[i].ResizeTo(3, 3);
    DR[i] = DRx[i] * DRy[i] * DRz[i];
    DRT[i].ResizeTo(3, 3);
    DRT[i] = DRzT[i] * DRyT[i] * DRxT[i];
  }
}

void ReconstructTrack::define_const_matrix()
{
  // define X, XT, Y, YT, Z, ZT

  // X, XT
  double X_arr[3] = {1, 0, 0};
  X.ResizeTo(1, 3);
  X = TMatrixD(1, 3, X_arr);
  XT.ResizeTo(1, 3);
  XT = X;
  XT.T();

  // Y, YT
  double Y_arr[3] = {0, 1, 0};
  Y.ResizeTo(1, 3);
  Y = TMatrixD(1, 3, Y_arr);
  YT.ResizeTo(1, 3);
  YT = Y;
  YT.T();

  // Z, ZT
  double Z_arr[3] = {0, 0, 1};
  Z.ResizeTo(1, 3);
  Z = TMatrixD(1, 3, Z_arr);
  ZT.ResizeTo(1, 3);
  ZT = Z;
  ZT.T();
}

void ReconstructTrack::define_origin_matrix()
{
  // define r0, r0T, R, RT

  // r0, r0T
  for (int i = 0; i < MMdetN; i++)
  {
    double r0_arr[3] = {0, 0, Zpos[i]};
    r0[i].ResizeTo(3, 1);
    r0[i] = TMatrixD(3, 1, r0_arr);
    r0T[i].ResizeTo(3, 1);
    r0T[i] = r0[i];
    r0T[i].T();
  }

  // R, RT, z-x-z, Rz(γ)Rx(β)Rz(α)
  for (int i = 0; i < MMdetN; i++)
  {
    // TRotation的每次Rotate都相当于在前面乘一次旋转矩阵
    TRotation R_rot = TRotation()
                          .RotateZ(0)
                          .RotateX(0)
                          .RotateZ(0);
    double R_arr[9] = {R_rot.XX(), R_rot.XY(), R_rot.XZ(), R_rot.YX(), R_rot.YY(),
                       R_rot.YZ(), R_rot.ZX(), R_rot.ZY(), R_rot.ZZ()};
    R[i].ResizeTo(3, 3);
    R[i] = TMatrixD(3, 3, R_arr);
    RT[i].ResizeTo(3, 3);
    RT[i] = R[i];
    RT[i].T();
  }
}

void ReconstructTrack::set_alignment_filename(string filename)
{
  for (int i = 0; i < MMdetN; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      alignment_list[i][j] = 0;
    }
  }
  m_alignment_filename = filename;
  define_align_matrix();
}

int ReconstructTrack::local_to_global(double &xx, double &yy, double &zz, int i)
{
  // pl
  double pl_arr[3] = {xx, yy, 0};
  TMatrixD pl(3, 1, pl_arr);

  // pg
  TMatrixD pg = R[i] * DeltaR[i] * pl + r0[i] + Deltar[i];

  xx = pg[0][0];
  yy = pg[1][0];
  zz = pg[2][0];

  return 0;
}

void ReconstructTrack::define_track_matrix(double kx, double ky, double bx, double by)
{
  // define k, kT, b, bT

  // k, kT
  double k_arr[3] = {kx, ky, 1};
  k.ResizeTo(3, 1);
  k = TMatrixD(3, 1, k_arr);
  kT.ResizeTo(3, 1);
  kT = k;
  kT.T();

  // b, bT
  double b_arr[3] = {bx, by, 0};
  b.ResizeTo(3, 1);
  b = TMatrixD(3, 1, b_arr);
  bT.ResizeTo(3, 1);
  bT = b;
  bT.T();
}

int ReconstructTrack::read_millepede_res(string filename, int alg)
{
  ifstream f(filename);
  string line;
  getline(f, line);
  while (getline(f, line))
  {
    istringstream line_stream(line);
    int label;
    double val;
    line_stream >> label >> val;

    int list_row = label / 100 - 1;
    int list_col = label % 100;
    if (alg == 1)
    {
      if (list_col == 3 || list_col == 4 || list_col == 5)
      {
        alignment_list[list_row][list_col] += val * 1e0;
        // mrad to rad
        // alignment_list[list_row][list_col] = alignment_list[list_row][list_col] > 0.05 ? 0.05 : alignment_list[list_row][list_col];
        // alignment_list[list_row][list_col] = alignment_list[list_row][list_col] < -0.05 ? -0.05 : alignment_list[list_row][list_col];
      }
      else
        alignment_list[list_row][list_col] += val;
    }
    else if (alg == 2)
    {
      alignment_list[list_row][list_col] = alignment_list[list_row][list_col] * 0.8 + val * 0.2;
    }
  }

  define_align_matrix();

  return 0;
}

bool ReconstructTrack::alignment_derlcx(float *derlcx, int i)
{
  double d_derlcx[4] = {0};
  bool ret = alignment_derlcx(d_derlcx, i);
  for (int j = 0; j < 4; j++)
  {
    derlcx[j] = d_derlcx[j];
  }
  return ret;
}

bool ReconstructTrack::alignment_derlcy(float *derlcy, int i)
{
  double d_derlcy[4] = {0};
  bool ret = alignment_derlcy(d_derlcy, i);
  for (int j = 0; j < 4; j++)
  {
    derlcy[j] = d_derlcy[j];
  }
  return ret;
}

bool ReconstructTrack::alignment_derglx(float *derglx, int i)
{
  double d_derglx[6] = {0};
  bool ret = alignment_derglx(d_derglx, i);
  for (int j = 0; j < 6; j++)
  {
    derglx[j] = d_derglx[j];
  }
  return ret;
}

bool ReconstructTrack::alignment_dergly(float *dergly, int i)
{
  double d_dergly[6] = {0};
  bool ret = alignment_dergly(d_dergly, i);
  for (int j = 0; j < 6; j++)
  {
    dergly[j] = d_dergly[j];
  }
  return ret;
}

bool ReconstructTrack::alignment_derlcx(double *derlcx, int i)
{
  TMatrixD G = kT * R[i] * DeltaR[i] * ZT;
  double g = G[0][0];
  if (g != 0)
  {
    TMatrixD dxl_dkx = X * DeltaRT[i] * RT[i] * (XT * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * kT * R[i] * DeltaR[i] * ZT - k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * X * R[i] * DeltaR[i] * ZT) * (1. / (g * g));
    TMatrixD dxl_dky = X * DeltaRT[i] * RT[i] * (YT * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * kT * R[i] * DeltaR[i] * ZT - k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * Y * R[i] * DeltaR[i] * ZT) * (1. / (g * g));
    TMatrixD dxl_dbx = X * DeltaRT[i] * RT[i] * (XT - k * X * R[i] * DeltaR[i] * ZT * (1. / g));
    TMatrixD dxl_dby = X * DeltaRT[i] * RT[i] * (YT - k * Y * R[i] * DeltaR[i] * ZT * (1. / g));

    // TMatrixD dxl_dkx = X * R_inv * XT * zx;
    // TMatrixD dxl_dky = X * R_inv * YT * zx;
    // TMatrixD dxl_dbx = X * R_inv * XT;
    // TMatrixD dxl_dby = X * R_inv * YT;

    derlcx[0] = dxl_dkx[0][0];
    derlcx[1] = dxl_dky[0][0];
    derlcx[2] = dxl_dbx[0][0];
    derlcx[3] = dxl_dby[0][0];

    return true;
  }
  else
  {
    return false;
  }
}

bool ReconstructTrack::alignment_derlcy(double *derlcy, int i)
{
  TMatrixD G = kT * R[i] * DeltaR[i] * ZT;
  double g = G[0][0];
  if (g != 0)
  {
    TMatrixD dyl_dkx = Y * DeltaRT[i] * RT[i] * (XT * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * kT * R[i] * DeltaR[i] * ZT - k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * X * R[i] * DeltaR[i] * ZT) * (1. / (g * g));
    TMatrixD dyl_dky = Y * DeltaRT[i] * RT[i] * (YT * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * kT * R[i] * DeltaR[i] * ZT - k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * Y * R[i] * DeltaR[i] * ZT) * (1. / (g * g));
    TMatrixD dyl_dbx = Y * DeltaRT[i] * RT[i] * (XT - k * X * R[i] * DeltaR[i] * ZT * (1. / g));
    TMatrixD dyl_dby = Y * DeltaRT[i] * RT[i] * (YT - k * Y * R[i] * DeltaR[i] * ZT * (1. / g));

    // TMatrixD dyl_dkx = Y * R_inv * XT * zx;
    // TMatrixD dyl_dky = Y * R_inv * YT * zx;
    // TMatrixD dyl_dbx = Y * R_inv * XT;
    // TMatrixD dyl_dby = Y * R_inv * YT;

    derlcy[0] = dyl_dkx[0][0];
    derlcy[1] = dyl_dky[0][0];
    derlcy[2] = dyl_dbx[0][0];
    derlcy[3] = dyl_dby[0][0];

    return true;
  }
  else
  {
    return false;
  }
}

bool ReconstructTrack::alignment_derglx(double *derglx, int i)
{
  TMatrixD G = kT * R[i] * DeltaR[i] * ZT;
  double g = G[0][0];
  if (g != 0)
  {
    TMatrixD dxl_dx = X * DeltaRT[i] * RT[i] * (k * X * R[i] * DeltaR[i] * ZT * (1. / g) - XT);
    TMatrixD dxl_dy = X * DeltaRT[i] * RT[i] * (k * Y * R[i] * DeltaR[i] * ZT * (1. / g) - YT);
    TMatrixD dxl_dz = X * DeltaRT[i] * RT[i] * (k * Z * R[i] * DeltaR[i] * ZT * (1. / g) - ZT);
    TMatrixD dxl_dtx =
        X * DeltaRzT[i] * DeltaRyT[i] * DRxT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * (1. / g) + b - r0[i] - Deltar[i]) + X * DeltaRT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DRx[i] * DeltaRy[i] * DeltaRz[i] * ZT * kT * R[i] * DeltaR[i] * ZT - k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * kT * R[i] * DRx[i] * DeltaRy[i] * DeltaRz[i] * ZT) * (1. / (g * g));
    TMatrixD dxl_dty =
        X * DeltaRzT[i] * DRyT[i] * DeltaRxT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * (1. / g) + b - r0[i] - Deltar[i]) + X * DeltaRT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaRx[i] * DRy[i] * DeltaRz[i] * ZT * kT * R[i] * DeltaR[i] * ZT - k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * kT * R[i] * DeltaRx[i] * DRy[i] * DeltaRz[i] * ZT) * (1. / (g * g));
    TMatrixD dxl_dtz =
        X * DRzT[i] * DeltaRyT[i] * DeltaRxT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * (1. / g) + b - r0[i] - Deltar[i]) + X * DeltaRT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaRx[i] * DeltaRy[i] * DRz[i] * ZT * kT * R[i] * DeltaR[i] * ZT - k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * kT * R[i] * DeltaRx[i] * DeltaRy[i] * DRz[i] * ZT) * (1. / (g * g));

    // derivative
    // TMatrixD dxl_dx = -1. * X * R_inv * XT;
    // TMatrixD dxl_dy = -1. * X * R_inv * YT;
    // TMatrixD dxl_dz = X * R_inv * (k - ZT);
    // TMatrixD dxl_dtx = X * (R_inv * k * Z * R - I) * DRx * R_inv * (k * zx + b - r0);
    // TMatrixD dxl_dty = X * (R_inv * k * Z * R - I) * DRy * R_inv * (k * zx + b - r0);
    // TMatrixD dxl_dtz = X * (R_inv * k * Z * R - I) * DRz * R_inv * (k * zx + b - r0);

    // dergl
    derglx[0] = dxl_dx[0][0];
    derglx[1] = dxl_dy[0][0];
    derglx[2] = dxl_dz[0][0];
    derglx[3] = dxl_dtx[0][0];
    derglx[4] = dxl_dty[0][0];
    derglx[5] = dxl_dtz[0][0];

    return true;
  }
  else
  {
    return false;
  }
}

bool ReconstructTrack::alignment_dergly(double *dergly, int i)
{
  TMatrixD G = kT * R[i] * DeltaR[i] * ZT;
  double g = G[0][0];
  if (g != 0)
  {
    TMatrixD dyl_dx = Y * DeltaRT[i] * RT[i] * (k * X * R[i] * DeltaR[i] * ZT * (1. / g) - XT);
    TMatrixD dyl_dy = Y * DeltaRT[i] * RT[i] * (k * Y * R[i] * DeltaR[i] * ZT * (1. / g) - YT);
    TMatrixD dyl_dz = Y * DeltaRT[i] * RT[i] * (k * Z * R[i] * DeltaR[i] * ZT * (1. / g) - ZT);
    TMatrixD dyl_dtx =
        Y * DeltaRzT[i] * DeltaRyT[i] * DRxT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * (1. / g) + b - r0[i] - Deltar[i]) + X * DeltaRT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DRx[i] * DeltaRy[i] * DeltaRz[i] * ZT * kT * R[i] * DeltaR[i] * ZT - k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * kT * R[i] * DRx[i] * DeltaRy[i] * DeltaRz[i] * ZT) * (1. / (g * g));
    TMatrixD dyl_dty =
        Y * DeltaRzT[i] * DRyT[i] * DeltaRxT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * (1. / g) + b - r0[i] - Deltar[i]) + X * DeltaRT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaRx[i] * DRy[i] * DeltaRz[i] * ZT * kT * R[i] * DeltaR[i] * ZT - k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * kT * R[i] * DeltaRx[i] * DRy[i] * DeltaRz[i] * ZT) * (1. / (g * g));
    TMatrixD dyl_dtz =
        Y * DRzT[i] * DeltaRyT[i] * DeltaRxT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * (1. / g) + b - r0[i] - Deltar[i]) + X * DeltaRT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaRx[i] * DeltaRy[i] * DRz[i] * ZT * kT * R[i] * DeltaR[i] * ZT - k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * kT * R[i] * DeltaRx[i] * DeltaRy[i] * DRz[i] * ZT) * (1. / (g * g));

    // derivative
    // TMatrixD dyl_dx = -1. * Y * R_inv * XT;
    // TMatrixD dyl_dy = -1. * Y * R_inv * YT;
    // TMatrixD dyl_dz = Y * R_inv * (k - ZT);
    // TMatrixD dyl_dtx = Y * (R_inv * k * Z * R - I) * DRx * R_inv * (k * zx + b - r0);
    // TMatrixD dyl_dty = Y * (R_inv * k * Z * R - I) * DRy * R_inv * (k * zx + b - r0);
    // TMatrixD dyl_dtz = Y * (R_inv * k * Z * R - I) * DRz * R_inv * (k * zx + b - r0);

    // dergl
    dergly[0] = dyl_dx[0][0];
    dergly[1] = dyl_dy[0][0];
    dergly[2] = dyl_dz[0][0];
    dergly[3] = dyl_dtx[0][0];
    dergly[4] = dyl_dty[0][0];
    dergly[5] = dyl_dtz[0][0];

    return true;
  }
  else
  {
    return false;
  }
}

double ReconstructTrack::expect_position_x(int i)
{
  TMatrixD G = kT * R[i] * DeltaR[i] * ZT;
  double g = G[0][0];
  if (g != 0)
  {
    TMatrixD xl = X * DeltaRT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * (1. / g) + b - r0[i] - Deltar[i]);
    return xl[0][0];
  }
  else
  {
    return -1000;
  }
}

double ReconstructTrack::expect_position_y(int i)
{
  TMatrixD G = kT * R[i] * DeltaR[i] * ZT;
  double g = G[0][0];
  if (g != 0)
  {
    TMatrixD yl = Y * DeltaRT[i] * RT[i] * (k * (r0T[i] + DeltarT[i] - bT) * R[i] * DeltaR[i] * ZT * (1. / g) + b - r0[i] - Deltar[i]);
    return yl[0][0];
  }
  else
  {
    return -1000;
  }
}

void ReconstructTrack::print_alignment_list()
{
  for (int i = 0; i < MMdetN; i++)
  {
    cout << "alignment list " << i << ": ";
    for (int j = 0; j < 6; j++)
    {
      cout << "j" << j << " : " << alignment_list[i][j] << ", ";
    }
    cout << endl;
  }
}

bool ReconstructTrack::least_square_fit(vector<double> z, vector<double> x, double sigma, double &kx, double &bx,
                                        double &chi2Dndf)
{
  if (x.size() == z.size() && x.size() >= 3)
  {
    int N = x.size();
    double Chi2 = 0;
    int NChi2 = 0;

    // x = kx*z+bx
    double Kx = 0, Bx = 0;
    double Mz = 0, Mx = 0;
    double Mzz = 0, Mxx = 0, Mzx = 0;
    for (int i = 0; i < N; i++)
    {
      Mz += z[i];
      Mx += x[i];
      Mzz += z[i] * z[i];
      Mxx += x[i] * x[i];
      Mzx += z[i] * x[i];
    }
    Mz /= N;
    Mx /= N;
    Mzz /= N;
    Mxx /= N;
    Mzx /= N;
    if (Mzz - Mz * Mz == 0)
    {
      return false;
    }
    Kx = (Mzx - Mz * Mx) / (Mzz - Mz * Mz);
    Bx = Mx - Kx * Mz;
    for (int i = 0; i < N; i++)
    {
      Chi2 += (x[i] - Kx * z[i] - Bx) * (x[i] - Kx * z[i] - Bx) / sigma / sigma;
      NChi2 += 1;
    }

    kx = Kx;
    bx = Bx;

    NChi2 -= 2; // 由于最小二乘法带来的自由度减少
    chi2Dndf = Chi2 / NChi2;

    double fit_R2;
    fit_R2 = (Mzx - Mz * Mx) * (Mzx - Mz * Mx) / (Mzz - Mz * Mz) / (Mxx - Mx * Mx);

    if (fit_R2 < alg_fit_R2_cut)
      return false;

    return true;
  }
  else if (x.size() == z.size() && x.size() == 2)
  {
    if (z[0] == z[1])
      return false;
    kx = (x[1] - x[0]) / (z[1] - z[0]);
    bx = x[0] - kx * z[0];
    chi2Dndf = 1;
    return true;
  }
  else if (x.size() != z.size())
  {
    cout << "Error: vectors is not the same size!" << endl;
    return false;
  }
  else
  {
    // cout << "Error: vectors' size is too small!" << endl;
    return false;
  }
}

double ReconstructTrack::get_data_chi2()
{
  define_align_matrix();

  int NChi2 = 0;
  double Chi2 = 0;
  int nentries = fdecTree[0]->GetEntriesFast();

  // cout << "Init data file done. Total entries: " << entries << endl;

  for (int evt = 0; evt < nentries; evt++)
  {
    vector<double> zx_origin, x_origin, zy_origin, y_origin;
    vector<double> a_zx, a_x, a_zy, a_y;
    vector<int> layer;
    for (int i = 0; i < MMdetN; i++)
    {
      fdecTree[i]->GetEntry(evt);
      if (sig)
      {
        x_origin.push_back(x * 0.4 - centre);
        zx_origin.push_back(z);
        y_origin.push_back(y * 0.4 - centre);
        zy_origin.push_back(z);

        double xx = x * 0.4 - centre;
        double yy = y * 0.4 - centre;
        double zz = z;
        local_to_global(xx, yy, zz, i);
        a_x.push_back(xx);
        a_zx.push_back(zz);
        a_y.push_back(yy);
        a_zy.push_back(zz);

        layer.push_back(i);
      }
    }
    // 只挑选穿过所有层探测器的数据
    if (a_zx.size() != MMdetN)
    {
      continue;
    }

    // 探测器分辨率
    double sigmax = 0.4 / sqrt(12);
    double sigmay = 0.4 / sqrt(12);

    double kx, bx, ky, by;
    double chi2Dndf_x, chi2Dndf_y;

    // 直线拟合
    if (!least_square_fit(a_zx, a_x, sigmax, kx, bx, chi2Dndf_x) || !least_square_fit(a_zy, a_y, sigmay, ky, by, chi2Dndf_y))
    {
      continue;
    }
    define_track_matrix(kx, ky, bx, by);

    // 计算这条直线的对齐权重
    // double weight = 1;
    double weight = sqrt((chi2Dndf_x + chi2Dndf_y) / 2) * 1;
    if (weight > 1000)
    {
      continue;
    }

    double sigma;
    for (int idx = 0; idx < layer.size(); idx++)
    {
      int i = layer[idx];

      sigma = sigmax * weight;
      double dx = x_origin[idx] - expect_position_x(i);
      // Chi2 += dx*dx/sigma/sigma;
      Chi2 += dx * dx;

      sigma = sigmay * weight;
      double dy = y_origin[idx] - expect_position_y(i);
      // Chi2 += dy*dy/sigma/sigma;
      Chi2 += dy * dy;

      NChi2++;
    }
  }
  return Chi2 / NChi2;
}

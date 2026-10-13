#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TKey.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace
{
constexpr double kInvalidValue = -999.0;

double PausePulseModel(double *x, double *par)
{
  const double b = par[3];
  if (b == 0.0)
    return par[0];

  const double u = (x[0] - par[2]) / b;
  return par[0] + par[1] * std::pow(u, 3) * std::exp(-u);
}

double SigmoidModel(double *x, double *par)
{
  const double k = par[3];
  if (k == 0.0)
    return par[0];

  return par[0] + par[1] / (1.0 + std::exp(-(x[0] - par[2]) / k));
}

double PauseShape(double u)
{
  return std::pow(u, 3) * std::exp(-u);
}

double SolvePauseRisingFractionU(double fraction)
{
  const double target = fraction * PauseShape(3.0);
  double lo = 0.0;
  double hi = 3.0;

  for (int i = 0; i < 100; ++i)
  {
    const double mid = 0.5 * (lo + hi);
    if (PauseShape(mid) < target)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
}

double SigmoidFractionTime(double t50, double k, double fraction)
{
  return t50 + k * std::log(fraction / (1.0 - fraction));
}

struct Stats
{
  int n = 0;
  double mean = 0.0;
  double varPopulation = 0.0;
  double varSample = 0.0;
  double rms = 0.0;
};

Stats GetStats(const std::vector<double> &values)
{
  Stats s;
  s.n = (int)values.size();
  if (values.empty())
    return s;

  const double sum = std::accumulate(values.begin(), values.end(), 0.0);
  s.mean = sum / values.size();

  double ss = 0.0;
  for (double v : values)
  {
    const double d = v - s.mean;
    ss += d * d;
  }

  s.varPopulation = ss / values.size();
  s.varSample = values.size() > 1 ? ss / (values.size() - 1) : 0.0;
  s.rms = std::sqrt(s.varPopulation);
  return s;
}

TH1D *MakeHist(const char *name, const char *title, const std::vector<double> &values)
{
  if (values.empty())
    return new TH1D(name, title, 10, 0.0, 1.0);

  const auto minmax = std::minmax_element(values.begin(), values.end());
  double lo = *minmax.first;
  double hi = *minmax.second;
  if (lo == hi)
  {
    lo -= 0.5;
    hi += 0.5;
  }
  else
  {
    const double pad = 0.05 * (hi - lo);
    lo -= pad;
    hi += pad;
  }

  const int bins = std::max(20, std::min(200, (int)std::sqrt((double)values.size()) * 4));
  TH1D *hist = new TH1D(name, title, bins, lo, hi);
  for (double v : values)
    hist->Fill(v);
  return hist;
}

int FindRisingFitStart(const std::vector<short> *wave, double ped, double pedrms, double maxA, int maxT)
{
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
  return fitStart;
}

struct PauseResult
{
  double a = kInvalidValue;
  double t20 = kInvalidValue;
  double t50 = kInvalidValue;
  double b = kInvalidValue;
  double baseline = kInvalidValue;
  double shapeA = kInvalidValue;
  double peakTime = kInvalidValue;
  double peakHeight = kInvalidValue;
  double chi2ndf = -1.0;
  int status = -1;
  int refit = 0;
  int fitStart = -1;
  int fitEnd = -1;
};

bool RunPauseFit(TGraph &graph,
                 int fitStart,
                 int fitEnd,
                 double ped,
                 double pedrms,
                 double maxA,
                 int maxT,
                 bool fixBaseline,
                 PauseResult &result)
{
  const double safePedRms = std::max(1.0, std::abs(pedrms));
  const double initialBshape = std::max(0.5, ((double)maxT - fitStart) / 3.0);
  const double initialA = std::max(1.0, maxA * std::exp(3.0) / 27.0);
  const double initialArrival = (double)maxT - 3.0 * initialBshape;

  TF1 fitFunc("pausePulseFit_local", PausePulseModel, fitStart, fitEnd, 4);
  fitFunc.SetParNames("B", "A", "a", "b");
  fitFunc.SetParameters(ped, initialA, initialArrival, initialBshape);
  if (fixBaseline)
  {
    fitFunc.FixParameter(0, ped);
  }
  else
  {
    fitFunc.SetParLimits(0,
                         ped - std::max(50.0, 10.0 * safePedRms),
                         ped + std::max(50.0, 10.0 * safePedRms));
  }
  fitFunc.SetParLimits(1, 0.0, std::max(100.0, 10.0 * maxA));
  fitFunc.SetParLimits(2, fitStart - 10.0, (double)fitEnd);
  fitFunc.SetParLimits(3, 0.1, 20.0);

  TFitResultPtr fitResult = graph.Fit(&fitFunc, "QSN");
  result.status = (int)fitResult;
  if (result.status != 0)
    return false;

  result.baseline = fitFunc.GetParameter(0);
  result.shapeA = fitFunc.GetParameter(1);
  result.a = fitFunc.GetParameter(2);
  result.b = fitFunc.GetParameter(3);

  const int ndf = fitFunc.GetNDF();
  if (ndf > 0)
    result.chi2ndf = fitFunc.GetChisquare() / ndf;

  const double u20 = SolvePauseRisingFractionU(0.20);
  const double u50 = SolvePauseRisingFractionU(0.50);
  result.t20 = result.a + result.b * u20;
  result.t50 = result.a + result.b * u50;
  result.peakTime = result.a + 3.0 * result.b;
  result.peakHeight = result.shapeA * PauseShape(3.0);
  return std::isfinite(result.a) && std::isfinite(result.t20) && std::isfinite(result.t50);
}

bool FitPauseWithBaselineCheck(const std::vector<short> *wave,
                               double ped,
                               double pedrms,
                               double maxA,
                               int maxT,
                               PauseResult &result)
{
  result = PauseResult();
  if (!wave || wave->size() < 4 || maxA <= 0.0 || maxT <= 0)
    return false;

  result.fitStart = FindRisingFitStart(wave, ped, pedrms, maxA, maxT);
  result.fitEnd = std::min(maxT, (int)wave->size() - 1);
  if (result.fitEnd - result.fitStart + 1 < 4)
    result.fitStart = std::max(0, result.fitEnd - 3);
  if (result.fitEnd - result.fitStart + 1 < 4)
    return false;

  TGraph graph;
  int nPoints = 0;
  for (int i = result.fitStart; i <= result.fitEnd; ++i)
    graph.SetPoint(nPoints++, i, (double)wave->at(i));

  PauseResult first = result;
  const bool firstOk = RunPauseFit(graph, result.fitStart, result.fitEnd, ped, pedrms, maxA, maxT, false, first);
  const double sigma3 = 3.0 * std::max(1.0, std::abs(pedrms));
  if (firstOk && std::abs(first.baseline - ped) <= sigma3)
  {
    result = first;
    return true;
  }

  PauseResult second = result;
  second.refit = firstOk ? 1 : 2;
  const bool secondOk = RunPauseFit(graph, result.fitStart, result.fitEnd, ped, pedrms, maxA, maxT, true, second);
  result = secondOk ? second : first;
  return secondOk;
}

struct SigmoidResult
{
  double baseline = kInvalidValue;
  double height = kInvalidValue;
  double t20 = kInvalidValue;
  double t50 = kInvalidValue;
  double k = kInvalidValue;
  double chi2ndf = -1.0;
  int status = -1;
  int fitStart = -1;
  int fitEnd = -1;
};

bool FitSigmoidRisingEdge(const std::vector<short> *wave,
                          double ped,
                          double pedrms,
                          double maxA,
                          int maxT,
                          SigmoidResult &result)
{
  result = SigmoidResult();
  if (!wave || wave->size() < 5 || maxA <= 0.0 || maxT <= 0)
    return false;

  const double safePedRms = std::max(1.0, std::abs(pedrms));
  const int thresholdStart = FindRisingFitStart(wave, ped, pedrms, maxA, maxT);
  result.fitStart = std::max(0, thresholdStart - 3);
  result.fitEnd = std::min(maxT, (int)wave->size() - 1);
  if (result.fitEnd - result.fitStart + 1 < 5)
    result.fitStart = std::max(0, result.fitEnd - 4);
  if (result.fitEnd - result.fitStart + 1 < 5)
    return false;

  double t50Initial = (double)maxT;
  for (int i = maxT; i > 0; --i)
  {
    const double adc = (double)wave->at(i) - ped;
    const double prevAdc = (double)wave->at(i - 1) - ped;
    if (adc >= 0.5 * maxA && prevAdc < 0.5 * maxA)
    {
      const double denom = adc - prevAdc;
      const double frac = denom != 0.0 ? (0.5 * maxA - prevAdc) / denom : 0.0;
      t50Initial = (i - 1) + frac;
      break;
    }
  }

  TGraph graph;
  int nPoints = 0;
  for (int i = result.fitStart; i <= result.fitEnd; ++i)
    graph.SetPoint(nPoints++, i, (double)wave->at(i));

  TF1 fitFunc("sigmoidRisingFit_local", SigmoidModel, result.fitStart, result.fitEnd, 4);
  fitFunc.SetParNames("B", "H", "t50", "k");
  fitFunc.SetParameters(ped, maxA, t50Initial, std::max(0.2, ((double)maxT - t50Initial) / 2.0));
  fitFunc.SetParLimits(0, ped - 3.0 * safePedRms, ped + 3.0 * safePedRms);
  fitFunc.SetParLimits(1, 0.0, std::max(100.0, 10.0 * maxA));
  fitFunc.SetParLimits(2, result.fitStart - 5.0, result.fitEnd + 5.0);
  fitFunc.SetParLimits(3, 0.05, 20.0);

  TFitResultPtr fitResult = graph.Fit(&fitFunc, "QSN");
  result.status = (int)fitResult;
  if (result.status != 0)
    return false;

  result.baseline = fitFunc.GetParameter(0);
  result.height = fitFunc.GetParameter(1);
  result.t50 = fitFunc.GetParameter(2);
  result.k = fitFunc.GetParameter(3);
  result.t20 = SigmoidFractionTime(result.t50, result.k, 0.20);

  const int ndf = fitFunc.GetNDF();
  if (ndf > 0)
    result.chi2ndf = fitFunc.GetChisquare() / ndf;

  return std::isfinite(result.t20);
}

void WriteStats(std::ostream &out, const char *label, const Stats &stats)
{
  out << label << " n: " << stats.n << "\n";
  out << label << " mean: " << stats.mean << "\n";
  out << label << " variance population: " << stats.varPopulation << "\n";
  out << label << " variance sample: " << stats.varSample << "\n";
  out << label << " rms: " << stats.rms << "\n";
}
}

void analyze_refit_sigmoid_variance(const char *inputName = "Combine/Tracker-step2-ana.root",
                                    const char *outputName = "refit_sigmoid_variance.root",
                                    const char *summaryName = "refit_sigmoid_variance_summary.txt",
                                    long long maxEntriesPerTree = -1)
{
  TFile input(inputName, "READ");
  if (input.IsZombie())
  {
    std::cerr << "Cannot open input file: " << inputName << std::endl;
    return;
  }

  TFile output(outputName, "RECREATE");
  TTree resultTree("fit_results", "Pause refit and sigmoid rising-edge fit results");

  char treeName[64] = "";
  int event = 0;
  int channel = 0;
  int maxT = 0;
  int pauseStatus = -1;
  int pauseRefit = 0;
  int pauseFitStart = -1;
  int pauseFitEnd = -1;
  int sigmoidStatus = -1;
  int sigmoidFitStart = -1;
  int sigmoidFitEnd = -1;
  double maxA = 0.0;
  double ped = 0.0;
  double pedrms = 0.0;
  double storedArrivalT = kInvalidValue;
  double pauseA = kInvalidValue;
  double pauseT20 = kInvalidValue;
  double pauseT50 = kInvalidValue;
  double pauseB = kInvalidValue;
  double pauseShapeA = kInvalidValue;
  double pauseShapeB = kInvalidValue;
  double pausePeakTime = kInvalidValue;
  double pausePeakHeight = kInvalidValue;
  double pauseChi2Ndf = -1.0;
  double sigmoidB = kInvalidValue;
  double sigmoidHeight = kInvalidValue;
  double sigmoidT20 = kInvalidValue;
  double sigmoidT50 = kInvalidValue;
  double sigmoidK = kInvalidValue;
  double sigmoidChi2Ndf = -1.0;

  resultTree.Branch("tree", treeName, "tree[64]/C");
  resultTree.Branch("event", &event);
  resultTree.Branch("channel", &channel);
  resultTree.Branch("maxA", &maxA);
  resultTree.Branch("maxT", &maxT);
  resultTree.Branch("ped", &ped);
  resultTree.Branch("pedrms", &pedrms);
  resultTree.Branch("storedArrivalT", &storedArrivalT);
  resultTree.Branch("pause_a", &pauseA);
  resultTree.Branch("pause_t20", &pauseT20);
  resultTree.Branch("pause_t50", &pauseT50);
  resultTree.Branch("pause_B", &pauseB);
  resultTree.Branch("pause_A", &pauseShapeA);
  resultTree.Branch("pause_b", &pauseShapeB);
  resultTree.Branch("pause_peakTime", &pausePeakTime);
  resultTree.Branch("pause_peakHeight", &pausePeakHeight);
  resultTree.Branch("pause_status", &pauseStatus);
  resultTree.Branch("pause_refit", &pauseRefit);
  resultTree.Branch("pause_chi2ndf", &pauseChi2Ndf);
  resultTree.Branch("pause_fitStart", &pauseFitStart);
  resultTree.Branch("pause_fitEnd", &pauseFitEnd);
  resultTree.Branch("sigmoid_t20", &sigmoidT20);
  resultTree.Branch("sigmoid_t50", &sigmoidT50);
  resultTree.Branch("sigmoid_B", &sigmoidB);
  resultTree.Branch("sigmoid_H", &sigmoidHeight);
  resultTree.Branch("sigmoid_k", &sigmoidK);
  resultTree.Branch("sigmoid_status", &sigmoidStatus);
  resultTree.Branch("sigmoid_chi2ndf", &sigmoidChi2Ndf);
  resultTree.Branch("sigmoid_fitStart", &sigmoidFitStart);
  resultTree.Branch("sigmoid_fitEnd", &sigmoidFitEnd);

  std::vector<double> pauseAValues;
  std::vector<double> pauseT20Values;
  std::vector<double> pauseT50Values;
  std::vector<double> sigmoidT20Values;
  long long totalEntries = 0;
  long long goodPauseFits = 0;
  long long pauseRefitCount = 0;
  long long goodSigmoidFits = 0;

  TIter nextKey(input.GetListOfKeys());
  while (TKey *key = (TKey *)nextKey())
  {
    TObject *obj = key->ReadObj();
    if (!obj || !obj->InheritsFrom(TTree::Class()))
    {
      delete obj;
      continue;
    }

    TTree *tree = (TTree *)obj;
    const std::string name = tree->GetName();
    if (name.rfind("Track", 0) != 0)
    {
      delete obj;
      continue;
    }

    std::vector<short> *wave = nullptr;
    const bool hasStoredArrivalT = tree->GetBranch("arrivalT") != nullptr;

    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("channel", &channel);
    tree->SetBranchAddress("maxA", &maxA);
    tree->SetBranchAddress("maxT", &maxT);
    tree->SetBranchAddress("ped", &ped);
    tree->SetBranchAddress("pedrms", &pedrms);
    tree->SetBranchAddress("wave", &wave);
    if (hasStoredArrivalT)
      tree->SetBranchAddress("arrivalT", &storedArrivalT);

    const long long nentries = tree->GetEntries();
    const long long entriesToProcess = maxEntriesPerTree > 0 ? std::min(nentries, maxEntriesPerTree) : nentries;
    std::cout << "processing " << name << ": " << entriesToProcess
              << " / " << nentries << " entries" << std::endl;

    for (long long i = 0; i < entriesToProcess; ++i)
    {
      if (i > 0 && i % 5000 == 0)
      {
        std::cout << "  " << name << " entry " << i << " / " << entriesToProcess
                  << ", pause good " << goodPauseFits
                  << ", sigmoid good " << goodSigmoidFits << std::endl;
      }

      storedArrivalT = kInvalidValue;
      tree->GetEntry(i);
      ++totalEntries;
      std::snprintf(treeName, sizeof(treeName), "%s", name.c_str());

      PauseResult pause;
      const bool pauseOk = FitPauseWithBaselineCheck(wave, ped, pedrms, maxA, maxT, pause);
      pauseA = pause.a;
      pauseT20 = pause.t20;
      pauseT50 = pause.t50;
      pauseB = pause.baseline;
      pauseShapeA = pause.shapeA;
      pauseShapeB = pause.b;
      pausePeakTime = pause.peakTime;
      pausePeakHeight = pause.peakHeight;
      pauseStatus = pause.status;
      pauseRefit = pause.refit;
      pauseChi2Ndf = pause.chi2ndf;
      pauseFitStart = pause.fitStart;
      pauseFitEnd = pause.fitEnd;

      SigmoidResult sigmoid;
      const bool sigmoidOk = FitSigmoidRisingEdge(wave, ped, pedrms, maxA, maxT, sigmoid);
      sigmoidT20 = sigmoid.t20;
      sigmoidT50 = sigmoid.t50;
      sigmoidB = sigmoid.baseline;
      sigmoidHeight = sigmoid.height;
      sigmoidK = sigmoid.k;
      sigmoidStatus = sigmoid.status;
      sigmoidChi2Ndf = sigmoid.chi2ndf;
      sigmoidFitStart = sigmoid.fitStart;
      sigmoidFitEnd = sigmoid.fitEnd;

      resultTree.Fill();

      if (pauseOk)
      {
        ++goodPauseFits;
        if (pause.refit != 0)
          ++pauseRefitCount;
        pauseAValues.push_back(pause.a);
        pauseT20Values.push_back(pause.t20);
        pauseT50Values.push_back(pause.t50);
      }
      if (sigmoidOk)
      {
        ++goodSigmoidFits;
        sigmoidT20Values.push_back(sigmoid.t20);
      }
    }

    delete obj;
  }

  Stats pauseAStats = GetStats(pauseAValues);
  Stats pauseT20Stats = GetStats(pauseT20Values);
  Stats pauseT50Stats = GetStats(pauseT50Values);
  Stats sigmoidT20Stats = GetStats(sigmoidT20Values);

  TH1D *hPauseA = MakeHist("h_pause_a", "Pause fit parameter a;sample index;counts", pauseAValues);
  TH1D *hPauseT20 = MakeHist("h_pause_t20", "Pause fit 20% peak time;sample index;counts", pauseT20Values);
  TH1D *hPauseT50 = MakeHist("h_pause_t50", "Pause fit 50% peak time;sample index;counts", pauseT50Values);
  TH1D *hSigmoidT20 = MakeHist("h_sigmoid_t20", "Sigmoid fit 20% height time;sample index;counts", sigmoidT20Values);

  output.cd();
  resultTree.Write();
  hPauseA->Write();
  hPauseT20->Write();
  hPauseT50->Write();
  hSigmoidT20->Write();

  TCanvas c("c_refit_sigmoid_variance", "refit and sigmoid timing variance", 1600, 900);
  c.Divide(2, 2);
  c.cd(1);
  hPauseA->Draw();
  c.cd(2);
  hPauseT20->Draw();
  c.cd(3);
  hPauseT50->Draw();
  c.cd(4);
  hSigmoidT20->Draw();
  c.Write();
  c.SaveAs("refit_sigmoid_variance.png");
  output.Close();

  std::ofstream summary(summaryName);
  summary << "input: " << inputName << "\n";
  summary << "output: " << outputName << "\n";
  summary << "total step2 waveforms: " << totalEntries << "\n";
  summary << "good pause fits: " << goodPauseFits << "\n";
  summary << "pause refits after B-ped 3sigma/failure check: " << pauseRefitCount << "\n";
  summary << "good sigmoid fits: " << goodSigmoidFits << "\n\n";
  WriteStats(summary, "pause_a", pauseAStats);
  summary << "\n";
  WriteStats(summary, "pause_t20", pauseT20Stats);
  summary << "\n";
  WriteStats(summary, "pause_t50", pauseT50Stats);
  summary << "\n";
  WriteStats(summary, "sigmoid_t20", sigmoidT20Stats);
  summary.close();

  std::cout << "total step2 waveforms: " << totalEntries << std::endl;
  std::cout << "good pause fits: " << goodPauseFits << std::endl;
  std::cout << "pause refits: " << pauseRefitCount << std::endl;
  std::cout << "good sigmoid fits: " << goodSigmoidFits << std::endl;
  std::cout << "pause a variance population: " << pauseAStats.varPopulation << std::endl;
  std::cout << "pause 20% variance population: " << pauseT20Stats.varPopulation << std::endl;
  std::cout << "pause 50% variance population: " << pauseT50Stats.varPopulation << std::endl;
  std::cout << "sigmoid 20% variance population: " << sigmoidT20Stats.varPopulation << std::endl;
  std::cout << "wrote: " << outputName << ", " << summaryName << std::endl;
}



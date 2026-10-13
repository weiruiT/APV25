#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TKey.h"
#include "TMath.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace
{
double PulseModel(double *x, double *par)
{
  const double b = par[3];
  if (b == 0.0)
    return par[0];

  const double u = (x[0] - par[2]) / b;
  return par[0] + par[1] * std::pow(u, 3) * std::exp(-u);
}

double PulseShape(double u)
{
  return std::pow(u, 3) * std::exp(-u);
}

double SolveRisingFractionU(double fraction)
{
  const double target = fraction * PulseShape(3.0);
  double lo = 0.0;
  double hi = 3.0;

  for (int i = 0; i < 100; ++i)
  {
    const double mid = 0.5 * (lo + hi);
    if (PulseShape(mid) < target)
      lo = mid;
    else
      hi = mid;
  }
  return 0.5 * (lo + hi);
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

  int bins = std::max(20, std::min(200, (int)std::sqrt((double)values.size()) * 4));
  TH1D *hist = new TH1D(name, title, bins, lo, hi);
  for (double v : values)
    hist->Fill(v);
  return hist;
}

bool FitOneWaveform(const std::vector<short> *wave,
                    double ped,
                    double pedrms,
                    double maxA,
                    int maxT,
                    double &a,
                    double &e20,
                    double &b,
                    double &baseline,
                    double &shapeA,
                    double &peakTime,
                    double &peakHeight,
                    int &fitStatus,
                    double &chi2ndf,
                    int &fitStart,
                    int &fitEnd)
{
  a = -999.0;
  e20 = -999.0;
  b = -999.0;
  baseline = -999.0;
  shapeA = -999.0;
  peakTime = -999.0;
  peakHeight = -999.0;
  fitStatus = -1;
  chi2ndf = -1.0;
  fitStart = -1;
  fitEnd = -1;

  if (!wave || wave->size() < 4 || maxA <= 0.0 || maxT <= 0)
    return false;

  const double safePedRms = std::max(1.0, std::abs(pedrms));
  const double threshold = std::max(3.0 * safePedRms, 0.1 * maxA);

  fitStart = 0;
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

  fitEnd = std::min(maxT, (int)wave->size() - 1);
  if (fitEnd - fitStart + 1 < 4)
    fitStart = std::max(0, fitEnd - 3);
  if (fitEnd - fitStart + 1 < 4)
    return false;

  TGraph gWave;
  int nPoints = 0;
  for (int i = fitStart; i <= fitEnd; ++i)
    gWave.SetPoint(nPoints++, i, (double)wave->at(i));

  const double initialB = ped;
  const double initialBshape = std::max(0.5, ((double)maxT - fitStart) / 3.0);
  const double initialA = std::max(1.0, maxA * std::exp(3.0) / 27.0);
  const double initialArrival = (double)maxT - 3.0 * initialBshape;

  TF1 fitFunc("pulseFitModel_local", PulseModel, fitStart, fitEnd, 4);
  fitFunc.SetParNames("B", "A", "a", "b");
  fitFunc.SetParameters(initialB, initialA, initialArrival, initialBshape);
  fitFunc.SetParLimits(0,
                       initialB - std::max(50.0, 10.0 * safePedRms),
                       initialB + std::max(50.0, 10.0 * safePedRms));
  fitFunc.SetParLimits(1, 0.0, std::max(100.0, 10.0 * maxA));
  fitFunc.SetParLimits(2, fitStart - 10.0, (double)fitEnd);
  fitFunc.SetParLimits(3, 0.1, 20.0);

  TFitResultPtr fitResult = gWave.Fit(&fitFunc, "QSN");
  fitStatus = (int)fitResult;
  if (fitStatus != 0)
    return false;

  baseline = fitFunc.GetParameter(0);
  shapeA = fitFunc.GetParameter(1);
  a = fitFunc.GetParameter(2);
  b = fitFunc.GetParameter(3);

  const int ndf = fitFunc.GetNDF();
  if (ndf > 0)
    chi2ndf = fitFunc.GetChisquare() / ndf;

  const double u20 = SolveRisingFractionU(0.20);
  e20 = a + b * u20;

  peakTime = a + 3.0 * b;
  peakHeight = shapeA * PulseShape(3.0);
  return std::isfinite(a) && std::isfinite(e20);
}
}

void analyze_arrival_e_variance(const char *inputName = "Combine/Tracker-step2-ana.root",
                                const char *outputName = "arrival_e_variance.root",
                                const char *summaryName = "arrival_e_variance_summary.txt",
                                long long maxEntriesPerTree = -1)
{
  TFile input(inputName, "READ");
  if (input.IsZombie())
  {
    std::cerr << "Cannot open input file: " << inputName << std::endl;
    return;
  }

  TFile output(outputName, "RECREATE");
  TTree resultTree("fit_results", "Per-waveform pulse fit results for a and e20");

  char treeName[64] = "";
  int event = 0;
  int channel = 0;
  int maxT = 0;
  int fitStatus = -1;
  int fitStart = -1;
  int fitEnd = -1;
  double maxA = 0.0;
  double ped = 0.0;
  double pedrms = 0.0;
  double storedArrivalT = -999.0;
  double a = -999.0;
  double e20 = -999.0;
  double b = -999.0;
  double baseline = -999.0;
  double shapeA = -999.0;
  double peakTime = -999.0;
  double peakHeight = -999.0;
  double chi2ndf = -1.0;

  resultTree.Branch("tree", treeName, "tree/C");
  resultTree.Branch("event", &event);
  resultTree.Branch("channel", &channel);
  resultTree.Branch("maxA", &maxA);
  resultTree.Branch("maxT", &maxT);
  resultTree.Branch("ped", &ped);
  resultTree.Branch("pedrms", &pedrms);
  resultTree.Branch("storedArrivalT", &storedArrivalT);
  resultTree.Branch("a", &a);
  resultTree.Branch("e20", &e20);
  resultTree.Branch("b", &b);
  resultTree.Branch("baseline", &baseline);
  resultTree.Branch("shapeA", &shapeA);
  resultTree.Branch("peakTime", &peakTime);
  resultTree.Branch("peakHeight", &peakHeight);
  resultTree.Branch("fitStatus", &fitStatus);
  resultTree.Branch("chi2ndf", &chi2ndf);
  resultTree.Branch("fitStart", &fitStart);
  resultTree.Branch("fitEnd", &fitEnd);

  std::vector<double> aValues;
  std::vector<double> e20Values;
  long long totalEntries = 0;
  long long goodFits = 0;

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
    bool hasStoredArrivalT = tree->GetBranch("arrivalT") != nullptr;

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
                  << ", good fits so far " << goodFits << std::endl;
      }
      storedArrivalT = -999.0;
      tree->GetEntry(i);
      ++totalEntries;

      const bool ok = FitOneWaveform(wave, ped, pedrms, maxA, maxT,
                                     a, e20, b, baseline, shapeA,
                                     peakTime, peakHeight, fitStatus,
                                     chi2ndf, fitStart, fitEnd);
      std::snprintf(treeName, sizeof(treeName), "%s", name.c_str());
      resultTree.Fill();

      if (!ok)
        continue;

      ++goodFits;
      aValues.push_back(a);
      e20Values.push_back(e20);
    }

    delete obj;
  }

  Stats aStats = GetStats(aValues);
  Stats e20Stats = GetStats(e20Values);

  TH1D *hA = MakeHist("h_a", "Fitted arrival parameter a;sample index;counts", aValues);
  TH1D *hE20 = MakeHist("h_e20", "Fitted 20% peak time e20;sample index;counts", e20Values);

  output.cd();
  resultTree.Write();
  hA->Write();
  hE20->Write();

  TCanvas c("c_arrival_e_variance", "arrival and e20 variance", 1200, 500);
  c.Divide(2, 1);
  c.cd(1);
  hA->Draw();
  c.cd(2);
  hE20->Draw();
  c.Write();
  c.SaveAs("arrival_e_variance.png");
  output.Close();

  std::ofstream summary(summaryName);
  summary << "input: " << inputName << "\n";
  summary << "output: " << outputName << "\n";
  summary << "total step2 waveforms: " << totalEntries << "\n";
  summary << "good fits: " << goodFits << "\n";
  summary << "\n";
  summary << "a mean: " << aStats.mean << "\n";
  summary << "a variance population: " << aStats.varPopulation << "\n";
  summary << "a variance sample: " << aStats.varSample << "\n";
  summary << "a rms: " << aStats.rms << "\n";
  summary << "\n";
  summary << "e20 mean: " << e20Stats.mean << "\n";
  summary << "e20 variance population: " << e20Stats.varPopulation << "\n";
  summary << "e20 variance sample: " << e20Stats.varSample << "\n";
  summary << "e20 rms: " << e20Stats.rms << "\n";
  summary.close();

  std::cout << "total step2 waveforms: " << totalEntries << std::endl;
  std::cout << "good fits: " << goodFits << std::endl;
  std::cout << "a variance population: " << aStats.varPopulation << std::endl;
  std::cout << "e20 variance population: " << e20Stats.varPopulation << std::endl;
  std::cout << "wrote: " << outputName << ", " << summaryName << std::endl;
}




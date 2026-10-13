#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

bool pass_shape(const std::vector<short> *wave, double ped, double pedrms)
{
  const int kExpectedWaveSamples = 30;
  if (!wave || (int)wave->size() != kExpectedWaveSamples)
    return false;

  const double sigma5 = 5.0 * pedrms;
  std::vector<int> valid_peaks;

  for (int i = 1; i < kExpectedWaveSamples - 1; ++i) {
    if (((double)wave->at(i) - ped) <= sigma5)
      continue;

    if (wave->at(i) == wave->at(i + 1)) {
      int plateau_start = i;
      int plateau_end = i;
      for (int j = i + 1; j < kExpectedWaveSamples; ++j) {
        if (wave->at(j) == wave->at(i)) plateau_end = j;
        else break;
      }
      if (plateau_start > 0 && plateau_end < kExpectedWaveSamples - 1) {
        bool left_ok = true;
        bool right_ok = true;
        for (int j = plateau_start - 1; j >= std::max(0, plateau_start - 3); --j) {
          if (wave->at(j) >= wave->at(i)) { left_ok = false; break; }
        }
        for (int j = plateau_end + 1; j <= std::min(kExpectedWaveSamples - 1, plateau_end + 3); ++j) {
          if (wave->at(j) >= wave->at(i)) { right_ok = false; break; }
        }
        if (left_ok && right_ok) {
          valid_peaks.push_back((plateau_start + plateau_end) / 2);
          i = plateau_end;
          continue;
        }
      }
    }

    bool is_local_max = true;
    const int left_bound = std::max(0, i - 3);
    const int right_bound = std::min(kExpectedWaveSamples - 1, i + 3);
    for (int j = left_bound; j <= right_bound; ++j) {
      if (j == i) continue;
      if (wave->at(j) >= wave->at(i)) { is_local_max = false; break; }
    }
    if (is_local_max) valid_peaks.push_back(i);
  }

  if (valid_peaks.size() != 1) return false;
  const int peak_pos = valid_peaks[0];
  if (peak_pos <= 0 || peak_pos >= kExpectedWaveSamples - 1) return false;
  if (((double)wave->at(peak_pos) - ped) <= sigma5) return false;
  if (wave->at(0) >= wave->at(peak_pos) || wave->at(kExpectedWaveSamples - 1) >= wave->at(peak_pos)) return false;
  return true;
}

struct Counts {
  Long64_t total = 0;
  Long64_t passCut = 0;
  Long64_t passShape = 0;
  Long64_t passStatus = 0;
  Long64_t passChi2Valid = 0;
  Long64_t passChi2_30 = 0;
  Long64_t passChi2_50 = 0;
  Long64_t passChi2_100 = 0;
};

void add_counts(Counts &a, const Counts &b)
{
  a.total += b.total;
  a.passCut += b.passCut;
  a.passShape += b.passShape;
  a.passStatus += b.passStatus;
  a.passChi2Valid += b.passChi2Valid;
  a.passChi2_30 += b.passChi2_30;
  a.passChi2_50 += b.passChi2_50;
  a.passChi2_100 += b.passChi2_100;
}

double pct(Long64_t n, Long64_t d)
{
  return d > 0 ? 100.0 * (double)n / (double)d : 0;
}

void print_counts(const char *name, const Counts &c)
{
  std::cout << std::left << std::setw(10) << name
            << std::right << std::setw(12) << c.total
            << std::setw(12) << c.passCut << std::setw(9) << std::fixed << std::setprecision(2) << pct(c.passCut, c.total)
            << std::setw(12) << c.passShape << std::setw(9) << pct(c.passShape, c.total)
            << std::setw(12) << c.passStatus << std::setw(9) << pct(c.passStatus, c.total)
            << std::setw(12) << c.passChi2Valid << std::setw(9) << pct(c.passChi2Valid, c.total)
            << std::setw(12) << c.passChi2_30 << std::setw(9) << pct(c.passChi2_30, c.total)
            << std::setw(12) << c.passChi2_50 << std::setw(9) << pct(c.passChi2_50, c.total)
            << std::setw(12) << c.passChi2_100 << std::setw(9) << pct(c.passChi2_100, c.total)
            << "\n";
}

Counts scan_tree(TTree *t)
{
  Counts c;
  if (!t) return c;

  double maxA=0, ped=0, pedrms=0, riseT=0, TOT=0, chi2=0;
  int maxT=0, status=0;
  std::vector<short> *wave=nullptr;
  t->SetBranchAddress("maxA", &maxA);
  t->SetBranchAddress("maxT", &maxT);
  t->SetBranchAddress("ped", &ped);
  t->SetBranchAddress("pedrms", &pedrms);
  t->SetBranchAddress("riseT", &riseT);
  t->SetBranchAddress("TOT", &TOT);
  t->SetBranchAddress("pulseFitStatus", &status);
  t->SetBranchAddress("pulseFitChi2NDF", &chi2);
  t->SetBranchAddress("wave", &wave);

  c.total = t->GetEntries();
  for (Long64_t i = 0; i < c.total; ++i) {
    t->GetEntry(i);
    bool okCut = maxA > 7.0 * pedrms && maxT > 0 && maxT < 20 && riseT > 0 && riseT < 15 && TOT > 0 && TOT < 20;
    if (!okCut) continue;
    ++c.passCut;
    if (!pass_shape(wave, ped, pedrms)) continue;
    ++c.passShape;
    if (status != 0) continue;
    ++c.passStatus;
    if (chi2 < 0) continue;
    ++c.passChi2Valid;
    if (chi2 <= 30) ++c.passChi2_30;
    if (chi2 <= 50) ++c.passChi2_50;
    if (chi2 <= 100) ++c.passChi2_100;
  }
  return c;
}

void diagnose_strict_filter_breakdown()
{
  TFile f("/home/thoriatong/workplace/BTTracker/arrivalT_run_20260702/Combine/Tracker-step2-ana.root", "READ");

  std::cout << "Filter breakdown on old loose STEP2 waveforms\n";
  std::cout << "Columns are cumulative counts after each condition; percentages are relative to old total.\n";
  std::cout << std::left << std::setw(10) << "tree"
            << std::right << std::setw(12) << "old"
            << std::setw(12) << "Cut" << std::setw(9) << "Cut%"
            << std::setw(12) << "+shape" << std::setw(9) << "%"
            << std::setw(12) << "+status" << std::setw(9) << "%"
            << std::setw(12) << "+chi2valid" << std::setw(9) << "%"
            << std::setw(12) << "chi2<=30" << std::setw(9) << "%"
            << std::setw(12) << "chi2<=50" << std::setw(9) << "%"
            << std::setw(12) << "chi2<=100" << std::setw(9) << "%"
            << "\n";
  Counts total;
  for (int det=0; det<6; ++det) {
    for (const char *xy : {"X", "Y"}) {
      TString name = Form("Track%s%d", xy, det);
      Counts c = scan_tree((TTree*)f.Get(name));
      add_counts(total, c);
      print_counts(name, c);
    }
  }
  std::cout << std::string(121, '-') << "\n";
  print_counts("TOTAL", total);

  std::cout << "\nSequential rejection relative to previous stage (TOTAL):\n";
  std::cout << "Cut rejects        : " << std::fixed << std::setprecision(2) << (100.0 - pct(total.passCut, total.total)) << "% of old step2\n";
  std::cout << "Shape rejects      : " << (100.0 - pct(total.passShape, total.passCut)) << "% of Cut-passing waveforms\n";
  std::cout << "Fit status rejects : " << (100.0 - pct(total.passStatus, total.passShape)) << "% of shape-passing waveforms\n";
  std::cout << "NDF/chi2 invalid   : " << (100.0 - pct(total.passChi2Valid, total.passStatus)) << "% of status-good waveforms\n";
  std::cout << "chi2>30 rejects    : " << (100.0 - pct(total.passChi2_30, total.passChi2Valid)) << "% of valid-chi2 waveforms\n";
}
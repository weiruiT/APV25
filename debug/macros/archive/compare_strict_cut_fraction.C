#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

Long64_t entries(TFile &f, const char *treeName, const char *selection = "")
{
  TTree *t = (TTree*)f.Get(treeName);
  if (!t) {
    std::cout << "missing tree " << treeName << " in " << f.GetName() << std::endl;
    return 0;
  }
  if (selection && selection[0]) return t->GetEntries(selection);
  return t->GetEntries();
}

void print_row(const char *name, Long64_t oldN, Long64_t newN)
{
  double keep = oldN > 0 ? (double)newN / oldN : 0;
  double cut = oldN > 0 ? 1.0 - keep : 0;
  std::cout << std::left << std::setw(12) << name
            << std::right << std::setw(14) << oldN
            << std::setw(14) << newN
            << std::setw(12) << std::fixed << std::setprecision(2) << keep * 100.0
            << std::setw(12) << cut * 100.0 << "\n";
}

void compare_step2()
{
  const char *oldPath = "/home/thoriatong/workplace/BTTracker/arrivalT_run_20260702/Combine/Tracker-step2-ana.root";
  const char *newPath = "/home/thoriatong/workplace/BTTracker/Combine/Tracker-step2-ana.root";
  TFile oldF(oldPath, "READ");
  TFile newF(newPath, "READ");

  std::cout << "\nSTEP2 waveform entries (old loose vs new strict)\n";
  std::cout << std::left << std::setw(12) << "tree"
            << std::right << std::setw(14) << "old"
            << std::setw(14) << "new"
            << std::setw(12) << "keep%"
            << std::setw(12) << "cut%" << "\n";

  Long64_t oldTotal = 0;
  Long64_t newTotal = 0;
  for (int det = 0; det < 6; ++det) {
    for (const char *xy : {"X", "Y"}) {
      TString name = Form("Track%s%d", xy, det);
      Long64_t oldN = entries(oldF, name);
      Long64_t newN = entries(newF, name);
      oldTotal += oldN;
      newTotal += newN;
      print_row(name, oldN, newN);
    }
  }
  std::cout << "------------" << std::string(52, '-') << "\n";
  print_row("TOTAL", oldTotal, newTotal);
}

void compare_step3()
{
  const char *oldPath = "/home/thoriatong/workplace/BTTracker/arrivalT_run_20260702/Combine/Tracker-step3-dec.root";
  const char *newPath = "/home/thoriatong/workplace/BTTracker/Combine/Tracker-step3-dec.root";
  TFile oldF(oldPath, "READ");
  TFile newF(newPath, "READ");

  std::cout << "\nSTEP3 decoded hit signals, sig_x/sig_y/sig (old loose vs new strict)\n";
  std::cout << std::left << std::setw(12) << "quantity"
            << std::right << std::setw(14) << "old"
            << std::setw(14) << "new"
            << std::setw(12) << "keep%"
            << std::setw(12) << "cut%" << "\n";

  Long64_t oldSigXTot = 0, newSigXTot = 0;
  Long64_t oldSigYTot = 0, newSigYTot = 0;
  Long64_t oldSigTot = 0, newSigTot = 0;
  for (int det = 0; det < 6; ++det) {
    TString tree = Form("Track%d", det);
    Long64_t oldX = entries(oldF, tree, "sig_x");
    Long64_t newX = entries(newF, tree, "sig_x");
    Long64_t oldY = entries(oldF, tree, "sig_y");
    Long64_t newY = entries(newF, tree, "sig_y");
    Long64_t oldS = entries(oldF, tree, "sig");
    Long64_t newS = entries(newF, tree, "sig");
    oldSigXTot += oldX; newSigXTot += newX;
    oldSigYTot += oldY; newSigYTot += newY;
    oldSigTot += oldS; newSigTot += newS;
    print_row(Form("T%d sig_x", det), oldX, newX);
    print_row(Form("T%d sig_y", det), oldY, newY);
    print_row(Form("T%d sig", det), oldS, newS);
  }
  std::cout << "------------" << std::string(52, '-') << "\n";
  print_row("TOTAL sig_x", oldSigXTot, newSigXTot);
  print_row("TOTAL sig_y", oldSigYTot, newSigYTot);
  print_row("TOTAL sig", oldSigTot, newSigTot);
}

void compare_step4()
{
  const char *oldPath = "/home/thoriatong/workplace/BTTracker/arrivalT_run_20260702/Combine/Tracker-step4-rec.root";
  const char *newPath = "/home/thoriatong/workplace/BTTracker/Combine/Tracker-step4-rec.root";
  TFile oldF(oldPath, "READ");
  TFile newF(newPath, "READ");

  std::cout << "\nSTEP4 reconstructed track flags (old loose vs new strict)\n";
  std::cout << std::left << std::setw(12) << "quantity"
            << std::right << std::setw(14) << "old"
            << std::setw(14) << "new"
            << std::setw(12) << "keep%"
            << std::setw(12) << "cut%" << "\n";
  print_row("entries", entries(oldF, "Track"), entries(newF, "Track"));
  print_row("useflag", entries(oldF, "Track", "useflag"), entries(newF, "Track", "useflag"));
  print_row("richuse", entries(oldF, "Track", "richuseflag"), entries(newF, "Track", "richuseflag"));
  print_row("btofuse", entries(oldF, "Track", "btofuseflag"), entries(newF, "Track", "btofuseflag"));
  print_row("ecaluse", entries(oldF, "Track", "ecaluseflag"), entries(newF, "Track", "ecaluseflag"));
}

void compare_strict_cut_fraction()
{
  compare_step2();
  compare_step3();
  compare_step4();
}
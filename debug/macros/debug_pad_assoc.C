#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <iostream>
#include <vector>
#include <map>

void debug_pad_assoc()
{
  TFile *f4 = TFile::Open("/home/thoriatong/workplace/BTTracker/Combine/Tracker-step4-rec.root");
  TFile *f2 = TFile::Open("/home/thoriatong/workplace/BTTracker/Combine/Tracker-step2-ana.root");
  TTree *t = (TTree *)f4->Get("Track");
  TTree *p = (TTree *)f2->Get("BMPadSignals");

  int event = -1;
  double Xfit[3] = {0}, Yfit[3] = {0};
  bool useflag = false;
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("XfitPars", Xfit);
  t->SetBranchAddress("YfitPars", Yfit);
  t->SetBranchAddress("useflag", &useflag);

  int pEvent = -1, rawChip = -1, position = -1;
  p->SetBranchAddress("event", &pEvent);
  p->SetBranchAddress("raw_chip", &rawChip);
  p->SetBranchAddress("position", &position);

  std::map<int, std::vector<int>> fired;
  for (Long64_t i = 0; i < p->GetEntries(); i++)
  {
    p->GetEntry(i);
    if (rawChip != 13)
      continue;
    fired[pEvent].push_back(position);
  }
  std::cout << "chip13 events with fired pads: " << fired.size() << std::endl;

  const double zBM = 309;
  TH1D hDist("dist", "nearest fired pad distance;|dr| [mm];tracks", 200, 0, 100);
  TH1D hDx("dx2", "xpred - xfired nearest;dx [mm];tracks", 200, -50, 50);
  TH1D hDy("dy2", "ypred - yfired nearest;dy [mm];tracks", 200, -50, 50);
  int nWithPad = 0, nNear = 0;
  for (Long64_t i = 0; i < t->GetEntries(); i++)
  {
    t->GetEntry(i);
    if (!useflag)
      continue;
    const double xp = Xfit[0] * zBM + Xfit[1];
    const double yp = Yfit[0] * zBM + Yfit[1];
    auto it = fired.find(event);
    if (it == fired.end())
      continue;
    nWithPad++;
    double best = 1e9, bx = 0, by = 0;
    for (int pos : it->second)
    {
      const int nx = pos % 11, ny = pos / 11;
      const double cx = nx * 9.0 + 4.5 - 49.5;
      const double cy = ny * 9.0 + 4.5 - 49.5;
      const double d = TMath::Hypot(xp - cx, yp - cy);
      if (d < best)
      {
        best = d;
        bx = xp - cx;
        by = yp - cy;
      }
    }
    hDist.Fill(best);
    hDx.Fill(bx);
    hDy.Fill(by);
    if (best < 30)
      nNear++;
  }
  std::cout << "useflag tracks with any chip13 pad fired: " << nWithPad
            << "  nearest<30mm: " << nNear << std::endl;
  std::cout << "dx mean=" << hDx.GetMean() << " rms=" << hDx.GetRMS()
            << "  dy mean=" << hDy.GetMean() << " rms=" << hDy.GetRMS() << std::endl;
  TString macroPath = gInterpreter->GetCurrentMacroName();
  TString outDir = gSystem->DirName(macroPath);
  outDir += "/../reports";
  TCanvas c("c", "pad assoc", 1200, 500);
  c.Divide(3, 1);
  c.cd(1);
  hDist.Draw();
  c.cd(2);
  hDx.Draw();
  c.cd(3);
  hDy.Draw();
  c.SaveAs(outDir + "/pad_assoc_diag.png");
  std::cout << "saved " << outDir + "/pad_assoc_diag.png" << std::endl;
}

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TBox.h>
#include <TPad.h>

#include <iostream>
#include <vector>
#include <map>

namespace
{
// Geometry: z = 0 at the lower PMT plane, z increases upward.
const double zTracker[4] = {533, 468, 217, 154}; // Tracker1..4
const double zBM1 = 364;                          // chip12
const double zBM2 = 309;                          // chip13

const int padGrid = 11;     // 11 x 11 pads
const double padPitch = 9.0; // mm
const double padCentre = padGrid / 2.0 * padPitch; // 49.5 mm local centre

double PulseModelLocal(double *x, double *par)
{
  const double u = (x[0] - par[2]) / par[3];
  return par[0] + par[1] * u * u * u * std::exp(-u);
}

bool PadIndexFromPosition(double xGlobal, double yGlobal,
                          double dxOff, double dyOff, int &ix, int &iy)
{
  const double xLocal = xGlobal + padCentre + dxOff;
  const double yLocal = yGlobal + padCentre + dyOff;
  ix = (int)TMath::FloorNint((xLocal - padPitch / 2.0) / padPitch);
  iy = (int)TMath::FloorNint((yLocal - padPitch / 2.0) / padPitch);
  return ix >= 0 && ix < padGrid && iy >= 0 && iy < padGrid;
}
}

void plot_track_pad_waveform(TString step4File = "Combine/Tracker-step4-rec.root",
                             TString step2File = "Combine/Tracker-step2-ana.root",
                             int eventWanted = -1,
                             int bmChip = 13,
                             double dxOff = 0.0,
                             double dyOff = 0.0,
                             bool calibrate = false)
{
  TFile *f4 = TFile::Open(step4File);
  if ((!f4 || f4->IsZombie()) && !step4File.BeginsWith("/"))
  {
    f4 = TFile::Open("../" + step4File);
    if (f4 && !f4->IsZombie())
      step4File = "../" + step4File;
  }
  TFile *f2 = TFile::Open(step2File);
  if ((!f2 || f2->IsZombie()) && !step2File.BeginsWith("/"))
  {
    f2 = TFile::Open("../" + step2File);
    if (f2 && !f2->IsZombie())
      step2File = "../" + step2File;
  }
  if (!f4 || f4->IsZombie() || !f2 || f2->IsZombie())
  {
    std::cerr << "ERROR: cannot open " << step4File << " or " << step2File << std::endl;
    return;
  }

  TTree *trackTree = dynamic_cast<TTree *>(f4->Get("Track"));
  TTree *padTree = dynamic_cast<TTree *>(f2->Get("BMPadSignals"));
  if (!trackTree || !padTree)
  {
    std::cerr << "ERROR: Track or BMPadSignals tree missing" << std::endl;
    return;
  }

  int event = -1;
  double MMx[4] = {0}, MMy[4] = {0}, MMz[4] = {0};
  bool MMsigx[4] = {false}, MMsigy[4] = {false};
  double Xfit[3] = {0}, Yfit[3] = {0};
  bool useflag = false;
  trackTree->SetBranchAddress("event", &event);
  trackTree->SetBranchAddress("MMx", MMx);
  trackTree->SetBranchAddress("MMy", MMy);
  trackTree->SetBranchAddress("MMz", MMz);
  trackTree->SetBranchAddress("MMsigx", MMsigx);
  trackTree->SetBranchAddress("MMsigy", MMsigy);
  trackTree->SetBranchAddress("XfitPars", Xfit);
  trackTree->SetBranchAddress("YfitPars", Yfit);
  trackTree->SetBranchAddress("useflag", &useflag);

  int pEvent = -1;
  int rawChip = -1;
  int position = -1;
  int maxT = -1;
  int pulseFitStatus = -1;
  int sigmoidFitStatus = -1;
  int sigmoidFitStart = -1;
  int sigmoidFitEnd = -1;
  double ped = 0.0;
  double pedrms = 1.0;
  double maxA = 0.0;
  double arrivalT_a = -999.0;
  double arrivalT_e = -999.0;
  double sigmoidB = -999.0;
  double sigmoidH = -999.0;
  double sigmoidK = -999.0;
  double sigmoidT50 = -999.0;
  double sigmoidT20 = -999.0;
  std::vector<short> *wave = nullptr;
  padTree->SetBranchAddress("event", &pEvent);
  padTree->SetBranchAddress("raw_chip", &rawChip);
  padTree->SetBranchAddress("position", &position);
  padTree->SetBranchAddress("maxT", &maxT);
  padTree->SetBranchAddress("ped", &ped);
  padTree->SetBranchAddress("pedrms", &pedrms);
  padTree->SetBranchAddress("maxA", &maxA);
  padTree->SetBranchAddress("arrivalT_a", &arrivalT_a);
  padTree->SetBranchAddress("arrivalT_e", &arrivalT_e);
  padTree->SetBranchAddress("pulseFitStatus", &pulseFitStatus);
  padTree->SetBranchAddress("sigmoidFitStatus", &sigmoidFitStatus);
  padTree->SetBranchAddress("sigmoid_B", &sigmoidB);
  padTree->SetBranchAddress("sigmoid_H", &sigmoidH);
  padTree->SetBranchAddress("sigmoid_k", &sigmoidK);
  padTree->SetBranchAddress("sigmoid_t50", &sigmoidT50);
  padTree->SetBranchAddress("arrivalT_sigmoid", &sigmoidT20);
  padTree->SetBranchAddress("sigmoid_fitStart", &sigmoidFitStart);
  padTree->SetBranchAddress("sigmoid_fitEnd", &sigmoidFitEnd);
  padTree->SetBranchAddress("wave", &wave);

  // Index pad waveforms by (event, position); keep the entry with max maxA.
  std::map<std::pair<int, int>, Long64_t> padIndex;
  const Long64_t padEntries = padTree->GetEntries();
  for (Long64_t i = 0; i < padEntries; ++i)
  {
    padTree->GetEntry(i);
    if (rawChip != bmChip || position < 0 || position >= padGrid * padGrid)
      continue;
    std::pair<int, int> key(pEvent, position);
    Long64_t best = -1;
    auto it = padIndex.find(key);
    if (it != padIndex.end())
    {
      padTree->GetEntry(it->second);
      const double bestA = maxA;
      padTree->GetEntry(i);
      if (maxA <= bestA)
        continue;
    }
    padIndex[key] = i;
  }

  const double zBM = (bmChip == 12) ? zBM1 : zBM2;

  TString macroPath = gInterpreter->GetCurrentMacroName();
  TString outDir = gSystem->DirName(macroPath);
  outDir += "/../reports";
  gSystem->mkdir(outDir, kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);

  if (calibrate)
  {
    TH2D hOff("pad_track_offset",
              Form("predicted - fired pad centre (chip%d);dx [mm];dy [mm]", bmChip),
              200, -20, 20, 200, -20, 20);
    TH1D hDx("pad_track_dx", "dx residual;dx [mm];events", 200, -20, 20);
    TH1D hDy("pad_track_dy", "dy residual;dy [mm];events", 200, -20, 20);
    int nAssoc = 0;
    int nTracks = 0;
    const Long64_t nTrack = trackTree->GetEntries();
    for (Long64_t i = 0; i < nTrack; ++i)
    {
      trackTree->GetEntry(i);
      if (!useflag)
        continue;
      ++nTracks;
      const double xPred = Xfit[0] * zBM + Xfit[1];
      const double yPred = Yfit[0] * zBM + Yfit[1];
      int ix = 0, iy = 0;
      if (!PadIndexFromPosition(xPred, yPred, dxOff, dyOff, ix, iy))
        continue;
      for (int di = -1; di <= 1; ++di)
      {
        for (int dj = -1; dj <= 1; ++dj)
        {
          const int nx = ix + di;
          const int ny = iy + dj;
          if (nx < 0 || nx >= padGrid || ny < 0 || ny >= padGrid)
            continue;
          const int pos = nx + padGrid * ny;
          auto it = padIndex.find(std::make_pair(event, pos));
          if (it == padIndex.end())
            continue;
          const double cx = nx * padPitch + padPitch / 2.0 - padCentre;
          const double cy = ny * padPitch + padPitch / 2.0 - padCentre;
          const double dx = xPred - cx;
          const double dy = yPred - cy;
          hOff.Fill(dx, dy);
          hDx.Fill(dx);
          hDy.Fill(dy);
          ++nAssoc;
        }
      }
    }
    std::cout << "calibrate: tracks=" << nTracks << " pad associations=" << nAssoc << std::endl;
    if (nAssoc > 0)
    {
      const int bx = hDx.GetMaximumBin();
      const int by = hDy.GetMaximumBin();
      std::cout << "recommended dxOff (pad centre vs tracker): " << hDx.GetBinCenter(bx) * (-1)
                << "  dyOff: " << hDy.GetBinCenter(by) * (-1) << std::endl;
      std::cout << "dx RMS=" << hDx.GetRMS() << "  dy RMS=" << hDy.GetRMS() << std::endl;
      TCanvas c("cCal", "pad-track offset", 1200, 500);
      c.Divide(2, 1);
      c.cd(1);
      hOff.Draw("COLZ");
      c.cd(2);
      hDx.SetLineColor(kBlue);
      hDx.Draw();
      hDy.SetLineColor(kRed);
      hDy.Draw("same");
      TString png = Form("%s/pad_track_calibrate_chip%d.png", outDir.Data(), bmChip);
      c.SaveAs(png);
      std::cout << "saved " << png << std::endl;
    }
    f4->Close();
    f2->Close();
    return;
  }

  const Long64_t nTrack = trackTree->GetEntries();
  int plotted = 0;
  const int maxPlots = 20;

  for (Long64_t i = 0; i < nTrack; ++i)
  {
    trackTree->GetEntry(i);
    if (!useflag)
      continue;
    if (eventWanted >= 0 && event != eventWanted)
      continue;

    const double xPred = Xfit[0] * zBM + Xfit[1];
    const double yPred = Yfit[0] * zBM + Yfit[1];
    int ix = 0, iy = 0;
    const bool inRange = PadIndexFromPosition(xPred, yPred, dxOff, dyOff, ix, iy);

    TCanvas c("cPadTrack", "track and pad waveforms", 1400, 900);
    c.Divide(1, 2);
    c.cd(1);
    TPad trackPad("trackPad", "track", 0.0, 0.0, 0.62, 1.0);
    TPad hitPad("hitPad", "hitmap", 0.62, 0.0, 1.0, 1.0);
    trackPad.Draw();
    hitPad.Draw();
    trackPad.cd();
    trackPad.SetMargin(0.10, 0.05, 0.12, 0.05);

    double zmin = 100;
    double zmax = 600;
    TH2D frame("frame", Form("event %d  chip%d  zBM=%.0f mm;z [mm];x / y [mm]", event, bmChip, zBM),
               10, zmin, zmax, 10, -60, 60);
    frame.SetStats(0);
    frame.Draw();

    TGraph gX, gY;
    for (int j = 0; j < 4; ++j)
    {
      if (MMsigx[j])
        gX.SetPoint(gX.GetN(), MMz[j], MMx[j]);
      if (MMsigy[j])
        gY.SetPoint(gY.GetN(), MMz[j], MMy[j]);
    }
    gX.SetMarkerStyle(20);
    gX.SetMarkerColor(kBlue);
    gX.SetLineColor(kBlue);
    gY.SetMarkerStyle(21);
    gY.SetMarkerColor(kRed);
    gY.SetLineColor(kRed);
    gX.Draw("P");
    gY.Draw("P");

    TGraph gFitX, gFitY;
    const int nFit = 51;
    for (int k = 0; k < nFit; ++k)
    {
      const double zz = zmin + (zmax - zmin) * k / (nFit - 1);
      gFitX.SetPoint(k, zz, Xfit[0] * zz + Xfit[1]);
      gFitY.SetPoint(k, zz, Yfit[0] * zz + Yfit[1]);
    }
    gFitX.SetLineColor(kBlue);
    gFitY.SetLineColor(kRed);
    gFitX.Draw("L");
    gFitY.Draw("L");

    TLine lBM(zmin, zBM, zmax, zBM);
    lBM.SetLineStyle(2);
    lBM.Draw();
    TMarker mX(xPred, zBM, 29);
    mX.SetMarkerColor(kBlue);
    mX.Draw();
    TMarker mY(yPred, zBM, 30);
    mY.SetMarkerColor(kRed);
    mY.Draw();

    TLegend legT(0.72, 0.72, 0.95, 0.90);
    legT.SetBorderSize(0);
    legT.SetFillStyle(0);
    legT.AddEntry(&gX, "MM x hit / fit", "lp");
    legT.AddEntry(&gY, "MM y hit / fit", "lp");
    legT.AddEntry(&lBM, Form("BM%d z=%.0f", bmChip == 12 ? 1 : 2, zBM), "l");
    legT.Draw();

    hitPad.cd();
    hitPad.SetMargin(0.10, 0.13, 0.12, 0.05);
    TH2D hHit(Form("hhit_%lld", i),
              Form("BM%d (chip%d) hit map, event %d;ix;iy",
                   bmChip == 12 ? 1 : 2, bmChip, event),
              11, -0.5, 10.5, 11, -0.5, 10.5);
    hHit.SetStats(0);
    hHit.GetZaxis()->SetTitle("maxA [ADC]");
    for (const auto &kv : padIndex)
    {
      if (kv.first.first != event)
        continue;
      padTree->GetEntry(kv.second);
      hHit.Fill(kv.first.second % padGrid, kv.first.second / padGrid, maxA);
    }
    hHit.Draw("COLZ");
    if (inRange)
    {
      TBox predBox(ix - 0.5, iy - 0.5, ix + 0.5, iy + 0.5);
      predBox.SetLineColor(kMagenta);
      predBox.SetLineWidth(3);
      predBox.SetFillStyle(0);
      predBox.Draw();
    }

    c.cd(2);
    gPad->SetMargin(0.05, 0.05, 0.10, 0.05);
    TPad *grid[3][3] = {nullptr};
    const int cellW = 3;
    const int cellH = 3;
    for (int r = 0; r < 3; ++r)
    {
      for (int cc = 0; cc < 3; ++cc)
      {
        const double x1 = cc * 0.32;
        const double y1 = 1.0 - (r + 1) * 0.31;
        TString name = Form("pad_%d_%d", r, cc);
        grid[r][cc] = new TPad(name, name, x1, y1, x1 + 0.32, y1 + 0.31);
        grid[r][cc]->Draw();
      }
    }

    for (int r = 0; r < 3; ++r)
    {
      for (int cc = 0; cc < 3; ++cc)
      {
        grid[r][cc]->cd();
        const int di = cc - 1;
        const int dj = 1 - r;
        const int nx = ix + di;
        const int ny = iy + dj;
        const bool isCentre = (di == 0 && dj == 0);
        const int pos = nx + padGrid * ny;
        auto it = padIndex.find(std::make_pair(event, pos));
        Long64_t pi = -1;
        if (it != padIndex.end())
          pi = it->second;

        if (pi >= 0)
          padTree->GetEntry(pi);

        TString title = Form("pad %d (%d,%d)%s", pos, nx, ny, isCentre ? " pred" : "");
        TH2D fcell(Form("f%d_%d_%lld", r, cc, i), title, 10, 0, 30, 10, -400, 400);
        fcell.SetStats(0);
        fcell.Draw();

        if (pi >= 0 && wave && wave->size() >= 4)
        {
          const int nSamples = (int)wave->size();
          const double thr = 5.0 * std::max(1.0, std::abs(pedrms));
          TGraph gW;
          double ymin = 1e9, ymax = -1e9;
          for (int k = 0; k < nSamples; ++k)
          {
            const double yv = (double)wave->at(k) - ped;
            gW.SetPoint(k, (double)k, yv);
            ymin = std::min(ymin, yv);
            ymax = std::max(ymax, yv);
          }
          ymin = std::min(ymin, -thr);
          ymax = std::max(ymax, thr);
          gW.SetMarkerStyle(20);
          gW.SetMarkerSize(0.5);
          gW.SetMinimum(ymin * 1.2);
          gW.SetMaximum(ymax * 1.2);
          gW.SetTitle(title + ";sample;ADC-ped");
          gW.Draw("ALP");

          TLine l0(0, 0, nSamples - 1, 0);
          l0.SetLineColor(kGray + 2);
          l0.SetLineStyle(2);
          l0.Draw();
          TLine lu(0, thr, nSamples - 1, thr);
          lu.SetLineColor(kGray + 1);
          lu.SetLineStyle(7);
          lu.Draw();
          TLine ld(0, -thr, nSamples - 1, -thr);
          ld.SetLineColor(kGray + 1);
          ld.SetLineStyle(7);
          ld.Draw();

          bool sigmoidOk = false;
          (void)sigmoidOk;
          if (maxA > 0 && maxT > 0 && maxT < nSamples)
          {
            const double thresh = std::max(3.0 * std::max(1.0, std::abs(pedrms)), 0.1 * maxA);
            int start = 0;
            for (int k = maxT; k > 0; --k)
            {
              const double adc = (double)wave->at(k) - ped;
              const double prev = (double)wave->at(k - 1) - ped;
              if (adc >= thresh && prev < thresh)
              {
                start = std::max(0, k - 1);
                break;
              }
            }
            int end = std::min(maxT, nSamples - 1);
            if (end - start + 1 < 4)
              start = std::max(0, end - 3);
            if (end - start + 1 >= 4)
            {
              TF1 fGam(Form("gam_%lld_%d_%d", i, r, cc), PulseModelLocal, 0, nSamples - 1, 4);
              const double b0 = std::max(0.5, ((double)maxT - start) / 3.0);
              const double a0 = (double)maxT - 3.0 * b0;
              const double A0 = std::max(1.0, maxA * std::exp(3.0) / 27.0);
              fGam.FixParameter(0, 0.0);
              fGam.SetParameters(0.0, A0, a0, b0);
              fGam.SetParLimits(1, 0.0, std::max(100.0, 10.0 * maxA));
              fGam.SetParLimits(2, start - 10.0, (double)end);
              fGam.SetParLimits(3, 0.1, 20.0);
              TGraph gFit;
              int np = 0;
              for (int k = start; k <= end; ++k)
                gFit.SetPoint(np++, (double)k, (double)wave->at(k) - ped);
              if (gFit.Fit(&fGam, "QSN") == 0)
              {
                fGam.SetLineColor(kAzure + 2);
                fGam.SetLineWidth(2);
                fGam.Draw("same");
              }
            }
          }

          if (sigmoidFitStatus == 0 && sigmoidH > 0.0 && sigmoidFitStart <= sigmoidFitEnd)
          {
            const int sigEnd = std::max(sigmoidFitStart,
                                        std::min(nSamples - 1, maxT + 4));
            TF1 fSig(Form("sig_%lld_%d_%d", i, r, cc),
                     "[0]/(1.0+exp(-(x-[1])/[2]))", sigmoidFitStart, sigEnd);
            fSig.SetParameters(sigmoidH, sigmoidT50, sigmoidK);
            fSig.SetLineColor(kOrange + 2);
            fSig.SetLineWidth(3);
            fSig.Draw("same");
            sigmoidOk = true;
          }
          else
          {
            TLatex fail(1.0, 30.0, "sig fail");
            fail.SetTextSize(0.08);
            fail.SetTextColor(kGray + 2);
            fail.Draw();
          }
          grid[r][cc]->SetBorderMode(1);
          grid[r][cc]->SetBorderSize(3);
          grid[r][cc]->SetFillColor(kGreen - 9);
        }
        else
        {
          grid[r][cc]->SetBorderMode(1);
          grid[r][cc]->SetBorderSize(3);
          if (isCentre)
          {
            TMarker mPred(15, 0, 29);
            mPred.SetMarkerSize(2);
            mPred.SetMarkerColor(kMagenta);
            mPred.Draw();
            TLatex lt(16, 30, Form("pred (%d,%d)", nx, ny));
            lt.SetTextSize(0.08);
            lt.Draw();
          }
          if (nx < 0 || nx >= padGrid || ny < 0 || ny >= padGrid)
          {
            TLatex lb(2, 30, "outside");
            lb.SetTextSize(0.08);
            lb.Draw();
          }
        }
      }
    }

    c.cd(1);
    TString png = Form("%s/pad_track_ev%d_chip%d", outDir.Data(), event, bmChip);
    if (eventWanted < 0)
      png += Form("_%d", plotted + 1);
    png += ".png";
    c.SaveAs(png);
    std::cout << "ev " << event << " pred=(" << xPred << ", " << yPred << ")";
    if (inRange)
      std::cout << " -> pad " << ix + padGrid * iy;
    else
      std::cout << " -> outside pad plane";
    std::cout << "  saved " << png << std::endl;

    ++plotted;
    if (eventWanted >= 0)
      break;
    if (plotted >= maxPlots)
    {
      std::cout << "Reached plot limit " << maxPlots << std::endl;
      break;
    }
  }

  f4->Close();
  f2->Close();
  if (plotted == 0)
    std::cout << "No useflag track found" << std::endl;
}

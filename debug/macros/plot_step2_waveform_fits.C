#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TBox.h>

#include <iostream>
#include <vector>

double PulseModelLocal(double *x, double *par)
{
  const double u = (x[0] - par[2]) / par[3];
  return par[0] + par[1] * u * u * u * std::exp(-u);
}

void plot_step2_waveform_fits(TString file = "Combine/Tracker-step2-ana.root",
                              TString treeName = "TrackX0",
                              int eventWanted = -1,
                              int sel = -1)
{
  if (sel < 0)
  {
    std::cerr << "Usage: plot_step2_waveform_fits(file, tree, event, sel)\n"
              << "  tree : TrackX0..TrackY3 or BMPadSignals\n"
              << "  event: >=0 single event, <0 scan all waveforms\n"
              << "  sel  : strip channel for Track trees, pad position for BMPadSignals"
              << std::endl;
    return;
  }

  TFile *f = TFile::Open(file);
  if ((!f || f->IsZombie()) && !file.BeginsWith("/"))
  {
    TString up = "../" + file;
    f = TFile::Open(up);
    if (f && !f->IsZombie())
    {
      std::cout << "Opened " << up << std::endl;
      file = up;
    }
  }
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open " << file << std::endl;
    return;
  }
  TTree *t = dynamic_cast<TTree *>(f->Get(treeName));
  if (!t)
  {
    std::cerr << "ERROR: tree " << treeName << " not found in " << file << std::endl;
    f->Close();
    return;
  }

  const bool isBM = (treeName == "BMPadSignals");

  int event = -1;
  int channel = -1;
  int rawChip = -1;
  int rawChannel = -1;
  int position = -1;
  int globalChannel = -1;
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
  double pulseFitChi2NDF = -1.0;
  double sigmoidFitChi2NDF = -1.0;
  std::vector<short> *wave = nullptr;

  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("maxT", &maxT);
  t->SetBranchAddress("ped", &ped);
  t->SetBranchAddress("pedrms", &pedrms);
  t->SetBranchAddress("maxA", &maxA);
  t->SetBranchAddress("arrivalT_a", &arrivalT_a);
  t->SetBranchAddress("arrivalT_e", &arrivalT_e);
  t->SetBranchAddress("pulseFitStatus", &pulseFitStatus);
  t->SetBranchAddress("pulseFitChi2NDF", &pulseFitChi2NDF);
  t->SetBranchAddress("sigmoidFitStatus", &sigmoidFitStatus);
  t->SetBranchAddress("sigmoidFitChi2NDF", &sigmoidFitChi2NDF);
  t->SetBranchAddress("sigmoid_B", &sigmoidB);
  t->SetBranchAddress("sigmoid_H", &sigmoidH);
  t->SetBranchAddress("sigmoid_k", &sigmoidK);
  t->SetBranchAddress("sigmoid_t50", &sigmoidT50);
  t->SetBranchAddress("arrivalT_sigmoid", &sigmoidT20);
  t->SetBranchAddress("sigmoid_fitStart", &sigmoidFitStart);
  t->SetBranchAddress("sigmoid_fitEnd", &sigmoidFitEnd);
  t->SetBranchAddress("wave", &wave);
  if (isBM)
  {
    t->SetBranchAddress("position", &position);
    t->SetBranchAddress("raw_chip", &rawChip);
    t->SetBranchAddress("raw_channel", &rawChannel);
    t->SetBranchAddress("global_channel", &globalChannel);
  }
  else
  {
    t->SetBranchAddress("channel", &channel);
  }

  TString macroPath = gInterpreter->GetCurrentMacroName();
  TString outDir = gSystem->DirName(macroPath);
  outDir += "/../reports";
  gSystem->mkdir(outDir, kTRUE);

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);

  std::map<int, std::map<int, double>> bmFired;
  TTree *bmTree = dynamic_cast<TTree *>(f->Get("BMPadSignals"));
  if (bmTree)
  {
    if (bmTree != t)
    {
      bmTree->SetBranchAddress("event", &event);
      bmTree->SetBranchAddress("raw_chip", &rawChip);
      bmTree->SetBranchAddress("position", &position);
      bmTree->SetBranchAddress("maxA", &maxA);
    }
    for (Long64_t b = 0; b < bmTree->GetEntries(); ++b)
    {
      bmTree->GetEntry(b);
      if (rawChip == 13 && position >= 0 && position < 121)
        bmFired[event][position] = maxA;
    }
  }

  TCanvas c("cWaveFit", "waveform and fits", 900, 800);
  c.Divide(2, 1);
  const Long64_t nEntries = t->GetEntries();
  const int maxPlots = 20;
  int plotted = 0;

  for (Long64_t i = 0; i < nEntries; ++i)
  {
    t->GetEntry(i);

    const bool match = isBM ? (position == sel) : (channel == sel);
    if (!match)
      continue;
    if (eventWanted >= 0 && event != eventWanted)
      continue;
    if (!wave || wave->size() < 4)
      continue;

    const int nSamples = (int)wave->size();
    const double thr = 5.0 * std::max(1.0, std::abs(pedrms));

    c.cd(1);
    gPad->Clear();
    gPad->SetMargin(0.12, 0.05, 0.12, 0.05);

    TGraph gSub;
    double ymin = 1e9;
    double ymax = -1e9;
    for (int k = 0; k < nSamples; ++k)
    {
      const double y = (double)wave->at(k) - ped;
      gSub.SetPoint(k, (double)k, y);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
    }
    ymin = std::min(ymin, -thr);
    ymax = std::max(ymax, thr);

    TString label = isBM
                        ? Form("%s ev%d pos%d chip%d ch%d gc%d", treeName.Data(),
                               event, position, rawChip, rawChannel, globalChannel)
                        : Form("%s ev%d ch%d", treeName.Data(), event, channel);

    const double pad = 0.15 * std::max(1.0, ymax - ymin) + 0.5;
    gSub.SetTitle(Form("%s;sample;ADC - ped", label.Data()));
    gSub.SetMarkerStyle(20);
    gSub.SetMarkerSize(0.8);
    gSub.SetMarkerColor(kBlack);
    gSub.SetLineColor(kBlack);
    gSub.SetMinimum(ymin - pad);
    gSub.SetMaximum(ymax + pad);
    gSub.Draw("ALP");

    const double xmin = gSub.GetXaxis()->GetXmin();
    const double xmax = gSub.GetXaxis()->GetXmax();

    TLine l0(xmin, 0.0, xmax, 0.0);
    l0.SetLineColor(kGray + 2);
    l0.SetLineStyle(2);
    l0.Draw();

    TLine lUp(xmin, thr, xmax, thr);
    lUp.SetLineColor(kGray + 1);
    lUp.SetLineStyle(7);
    lUp.Draw();
    TLine lDn(xmin, -thr, xmax, -thr);
    lDn.SetLineColor(kGray + 1);
    lDn.SetLineStyle(7);
    lDn.Draw();

    TLine lT20(sigmoidT20, ymin - pad, sigmoidT20, ymax + pad);
    lT20.SetLineColor(kViolet + 1);
    lT20.SetLineStyle(3);
    lT20.Draw();

    // Legacy Pause diagnostic refit retained intact, but disabled while the
    // analysis policy is sigmoid-only.
#if 0
    double refitA = -999.0;
    int refitStatus = -1;
    bool pauseOk = false;
    TF1 fPause(Form("pause_%lld", i), PulseModelLocal, 0, nSamples - 1, 4);
    fPause.SetParNames("B", "A", "a", "b");
    if (maxA > 0 && maxT > 0 && maxT < nSamples)
    {
      const double thresh = std::max(3.0 * std::max(1.0, std::abs(pedrms)), 0.1 * maxA);
      int start = 0;
      for (int k = maxT; k > 0; --k)
      {
        const double adc = (double)wave->at(k) - ped;
        const double prevAdc = (double)wave->at(k - 1) - ped;
        if (adc >= thresh && prevAdc < thresh)
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
        const double b0 = std::max(0.5, ((double)maxT - start) / 3.0);
        const double a0 = (double)maxT - 3.0 * b0;
        const double A0 = std::max(1.0, maxA * std::exp(3.0) / 27.0);
        fPause.FixParameter(0, 0.0);
        fPause.SetParameters(0.0, A0, a0, b0);
        fPause.SetParLimits(1, 0.0, std::max(100.0, 10.0 * maxA));
        fPause.SetParLimits(2, start - 10.0, (double)end);
        fPause.SetParLimits(3, 0.1, 20.0);

        TGraph gFit;
        int np = 0;
        for (int k = start; k <= end; ++k)
          gFit.SetPoint(np++, (double)k, (double)wave->at(k) - ped);
        refitStatus = gFit.Fit(&fPause, "QSN");
        if (refitStatus == 0)
        {
          refitA = fPause.GetParameter(2);
          pauseOk = true;
        }
      }
    }
    if (pauseOk)
    {
      fPause.SetLineColor(kAzure + 2);
      fPause.SetLineWidth(2);
      fPause.Draw("same");
    }
#endif

    bool sigmoidOk = false;
    TF1 fSig(Form("sig_%lld", i), "[0]/(1.0+exp(-(x-[1])/[2]))",
             0, nSamples - 1);
    if (sigmoidFitStatus == 0 && sigmoidH > 0.0 && sigmoidFitStart <= sigmoidFitEnd)
    {
      const int sigEnd = std::max(sigmoidFitStart,
                                  std::min(nSamples - 1, maxT + 4));
      fSig.SetRange(sigmoidFitStart, sigEnd);
      fSig.SetParameters(sigmoidH, sigmoidT50, sigmoidK);
      fSig.SetLineColor(kOrange + 2);
      fSig.SetLineWidth(3);
      fSig.Draw("same");
      sigmoidOk = true;
    }

    TLegend leg(0.62, 0.56, 0.90, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.032);
    leg.AddEntry(&gSub, "wave - ped", "lp");
    leg.AddEntry(&l0, "0 (ped)", "l");
    leg.AddEntry(&lUp, Form("#pm5#sigma = #pm%.2f", thr), "l");
    if (sigmoidT20 > -900)
      leg.AddEntry(&lT20, Form("t20 = %.2f", sigmoidT20), "l");
#if 0
    if (pauseOk)
      leg.AddEntry(&fPause, "Pause fit", "l");
#endif
    if (sigmoidOk)
      leg.AddEntry(&fSig, "Sigmoid fit", "l");
    else
    {
      TLatex fail(0.02, 0.92, "sigmoid fit failed");
      fail.SetNDC();
      fail.SetTextColor(kGray + 2);
      fail.SetTextSize(0.045);
      fail.Draw();
    }

    leg.Draw();

    c.cd(2);
    gPad->Clear();
    gPad->SetMargin(0.10, 0.13, 0.12, 0.08);
    TH2D hHit(Form("hhit_%lld", i), Form("BM2 (chip13) hit map, event %d;ix;iy", event),
              11, -0.5, 10.5, 11, -0.5, 10.5);
    hHit.SetStats(0);
    hHit.GetZaxis()->SetTitle("maxA [ADC]");
    auto firedIt = bmFired.find(event);
    if (firedIt != bmFired.end())
    {
      for (const auto &kv : firedIt->second)
        hHit.Fill(kv.first % 11, kv.first / 11, kv.second);
    }
    hHit.Draw("COLZ");
    if (isBM && sel >= 0 && sel < 121)
    {
      const int sx = sel % 11;
      const int sy = sel / 11;
      TBox selBox(sx - 0.5, sy - 0.5, sx + 0.5, sy + 0.5);
      selBox.SetLineColor(kOrange + 2);
      selBox.SetLineWidth(3);
      selBox.SetFillStyle(0);
      selBox.Draw();
    }

    TString pngName = Form("%s/%s_ev%d_%s%d", outDir.Data(), treeName.Data(),
                           event, isBM ? "pos" : "ch", sel);
    if (eventWanted < 0)
      pngName += Form("_%d", plotted + 1);
    pngName += ".png";
    c.SaveAs(pngName);

    std::cout << label << " -> " << pngName << "\n"
              << "  maxA=" << maxA << " ped=" << ped << " pedrms=" << pedrms
              // Legacy Pause output retained, but disabled:
              // << "\n  stored a=" << arrivalT_a << " e=" << arrivalT_e
              // << " refit a=" << refitA << " refitStatus=" << refitStatus
              // << " pulseChi2NDF=" << pulseFitChi2NDF
              << "\n  sigmoid t50=" << sigmoidT50 << " t20=" << sigmoidT20
              << " B=" << sigmoidB << " H=" << sigmoidH << " k=" << sigmoidK
              << " status=" << sigmoidFitStatus
              << " chi2NDF=" << sigmoidFitChi2NDF << std::endl;

    ++plotted;
    if (eventWanted >= 0)
      break;
    if (plotted >= maxPlots)
    {
      std::cout << "Reached plot limit " << maxPlots << "; increase maxPlots "
                << "inside the macro to plot more." << std::endl;
      break;
    }
  }

  f->Close();

  if (plotted == 0)
    std::cout << "No matching waveform found in " << treeName
              << " (entries=" << nEntries << "). Check the sel/event values."
              << std::endl;
  else
    std::cout << "Plotted " << plotted << " waveform(s) to " << outDir << std::endl;
}

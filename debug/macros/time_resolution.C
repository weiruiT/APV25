#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TStyle.h>
#include "plot_peak_timing.C"

#include <iostream>
#include <string>
#include <vector>

struct Def
{
  std::string name;
  std::string trMeanX;
  std::string trMeanY;
  std::string trMinX;
  std::string trMinY;
  std::string padMean;
  std::string padMin;
};

TFile *OpenWithFallback(TString &name)
{
  TFile *f = TFile::Open(name);
  if ((!f || f->IsZombie()) && !name.BeginsWith("/"))
  {
    TString up = "../" + name;
    f = TFile::Open(up);
    if (f && !f->IsZombie())
      name = up;
  }
  return f;
}

void FillHist(TH1D *h, TTree *t, const char *branch)
{
  if (!t || !t->GetBranch(branch))
    return;
  for (Long64_t i = 0; i < t->GetEntries(); ++i)
  {
    t->GetEntry(i);
    const double value = t->GetLeaf(branch)->GetValue();
    if (value > -900.0)
      h->Fill(value * 25.0);
  }
}

void pad_time_resolution(TString padFile = "../Combine/Pad-step3-rec.root")
{
  TFile *fp = OpenWithFallback(padFile);
  if (!fp || fp->IsZombie())
  {
    std::cerr << "ERROR: cannot open pad file" << std::endl;
    return;
  }
  TTree *pad[2] = {dynamic_cast<TTree *>(fp->Get("BM1Events")),
                   dynamic_cast<TTree *>(fp->Get("BM2Events"))};
  if (!pad[0] || !pad[1])
  {
    std::cerr << "ERROR: BM1Events/BM2Events trees missing" << std::endl;
    fp->Close();
    return;
  }

  for (int bm = 1; bm <= 2; ++bm)
  {
    const char *defNames[3] = {"a", "e", "sigmoid"};
    const char *defTitles[3] = {"Pause a", "e (20% rising)", "sigmoid t20"};
    const char *meanBranches[3] = {"arrivalT_a_mean", "arrivalT_e_mean",
                                   "arrivalT_sigmoid_mean"};
    const char *minBranches[3] = {"arrivalT_a_min", "arrivalT_e_min",
                                  "arrivalT_sigmoid_min"};
    Color_t colors[3] = {kAzure + 2, kOrange + 2, kViolet + 2};

    TString macroPath = gInterpreter->GetCurrentMacroName();
    TString outDir = gSystem->DirName(macroPath);
    outDir += "/../reports";
    gSystem->mkdir(outDir, kTRUE);
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas(Form("cPadTR_BM%d", bm),
                             Form("BM%d pad time resolution", bm), 1400, 650);
    c->Divide(2, 1);

    std::vector<TObject *> keepAlive;
    keepAlive.push_back(c);

    TH1D *hMean[3] = {nullptr};
    TH1D *hMin[3] = {nullptr};

    // for (int k = 0; k < 3; ++k) // Pause a/e disabled.
    for (int k : PeakTimingPlots::ActiveTimingModels)
    {
      hMean[k] = new TH1D(Form("hP_BM%d_%s_mean", bm, defNames[k]), ";t [ns];events",
                          130, -25.0, 625.0);
      hMin[k] = new TH1D(Form("hP_BM%d_%s_min", bm, defNames[k]), ";t [ns];events",
                         130, -25.0, 625.0);
      keepAlive.push_back(hMean[k]);
      keepAlive.push_back(hMin[k]);

      FillHist(hMean[k], pad[bm - 1], meanBranches[k]);
      FillHist(hMin[k], pad[bm - 1], minBranches[k]);

      hMean[k]->SetLineColor(colors[k]);
      hMin[k]->SetLineColor(colors[k]);
      hMean[k]->SetLineWidth(2);
      hMin[k]->SetLineWidth(2);
    }

    c->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.05);

    // for (int k = 0; k < 3; ++k) // Pause a/e disabled.
    for (int k : PeakTimingPlots::ActiveTimingModels)
      hMean[k]->Draw(k == PeakTimingPlots::ActiveTimingModels[0] ? "" : "same");

    TLegend *legM = new TLegend(0.55, 0.60, 0.95, 0.90);
    keepAlive.push_back(legM);
    legM->SetBorderSize(0);
    legM->SetFillStyle(0);
    legM->SetTextSize(0.030);
    legM->SetMargin(0.08);
    // for (int k = 0; k < 3; ++k) // Pause a/e disabled.
    for (int k : PeakTimingPlots::ActiveTimingModels)
    {
      legM->AddEntry(hMean[k],
                     Form("%s: %.1f ns (N=%.0f)", defTitles[k],
                          hMean[k]->GetRMS(), hMean[k]->GetEntries()),
                     "l");
    }
    legM->Draw();

    TLatex *tMean = new TLatex(0.16, 0.91, Form("BM%d mean arrival time", bm));
    keepAlive.push_back(tMean);
    tMean->SetNDC();
    tMean->SetTextSize(0.040);
    tMean->SetTextFont(62);
    tMean->Draw();

    TLatex *tNoteM = new TLatex(0.16, 0.85, "fit-failure excluded (t > -900)");
    keepAlive.push_back(tNoteM);
    tNoteM->SetNDC();
    tNoteM->SetTextSize(0.028);
    tNoteM->SetTextColor(kGray + 2);
    tNoteM->Draw();

    c->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.02);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.05);

    // for (int k = 0; k < 3; ++k) // Pause a/e disabled.
    for (int k : PeakTimingPlots::ActiveTimingModels)
      hMin[k]->Draw(k == PeakTimingPlots::ActiveTimingModels[0] ? "" : "same");

    TLegend *legN = new TLegend(0.55, 0.60, 0.95, 0.90);
    keepAlive.push_back(legN);
    legN->SetBorderSize(0);
    legN->SetFillStyle(0);
    legN->SetTextSize(0.030);
    legN->SetMargin(0.08);
    // for (int k = 0; k < 3; ++k) // Pause a/e disabled.
    for (int k : PeakTimingPlots::ActiveTimingModels)
    {
      legN->AddEntry(hMin[k],
                     Form("%s: %.1f ns (N=%.0f)", defTitles[k],
                          hMin[k]->GetRMS(), hMin[k]->GetEntries()),
                     "l");
    }
    legN->Draw();

    TLatex *tMin = new TLatex(0.16, 0.91, Form("BM%d min arrival time", bm));
    keepAlive.push_back(tMin);
    tMin->SetNDC();
    tMin->SetTextSize(0.040);
    tMin->SetTextFont(62);
    tMin->Draw();

    TLatex *tNoteN = new TLatex(0.16, 0.85, "fit-failure excluded (t > -900)");
    keepAlive.push_back(tNoteN);
    tNoteN->SetNDC();
    tNoteN->SetTextSize(0.028);
    tNoteN->SetTextColor(kGray + 2);
    tNoteN->Draw();

    TString png = Form("%s/pad_time_resolution_bm%d.png", outDir.Data(), bm);
    c->SaveAs(png);
    std::cout << "saved " << png << std::endl;

    std::cout << "BM" << bm << " pad time resolution (RMS over events, ns):" << std::endl;
    // for (int k = 0; k < 3; ++k) // Pause a/e disabled.
    for (int k : PeakTimingPlots::ActiveTimingModels)
    {
      std::cout << "  " << defNames[k] << " mean: RMS = " << hMean[k]->GetRMS()
                << " ns, mean = " << hMean[k]->GetMean()
                << " ns, N = " << hMean[k]->GetEntries() << std::endl;
      std::cout << "  " << defNames[k] << " min : RMS = " << hMin[k]->GetRMS()
                << " ns, mean = " << hMin[k]->GetMean()
                << " ns, N = " << hMin[k]->GetEntries() << std::endl;
    }

    for (TObject *obj : keepAlive)
    {
      if (obj != c)
        delete obj;
    }
    delete c;
  }
  fp->Close();
}

void track_time_resolution(TString trackFile = "../Combine/Pad-step3-rec.root")
{
  TFile *ft = OpenWithFallback(trackFile);
  if (!ft || ft->IsZombie())
  {
    std::cerr << "ERROR: cannot open track arrival file" << std::endl;
    return;
  }
  TTree *ta = dynamic_cast<TTree *>(ft->Get("TrackArrivalEvents"));
  if (!ta)
  {
    std::cerr << "ERROR: TrackArrivalEvents tree missing in " << trackFile
              << std::endl;
    ft->Close();
    return;
  }

  TString macroPath = gInterpreter->GetCurrentMacroName();
  TString outDir = gSystem->DirName(macroPath);
  outDir += "/../reports";
  gSystem->mkdir(outDir, kTRUE);
  gStyle->SetOptStat(0);

  const char *defNames[3] = {"a", "e", "sigmoid"};
  const char *defTitles[3] = {"Pause a", "e (20% rising)", "sigmoid t20"};
  Color_t colors[3] = {kAzure + 2, kOrange + 2, kViolet + 2};

  double evMean[3] = {0}, evMin[3] = {0};
  double bmMean[3][2] = {{0}}, bmMin[3][2] = {{0}};
  double trMean[3][3][2] = {{{0}}}, trMin[3][3][2] = {{{0}}};
  ta->SetBranchAddress("evA_mean", &evMean[0]);
  ta->SetBranchAddress("evA_min", &evMin[0]);
  ta->SetBranchAddress("evE_mean", &evMean[1]);
  ta->SetBranchAddress("evE_min", &evMin[1]);
  ta->SetBranchAddress("evS_mean", &evMean[2]);
  ta->SetBranchAddress("evS_min", &evMin[2]);
  ta->SetBranchAddress("bmA_mean", bmMean[0]);
  ta->SetBranchAddress("bmA_min", bmMin[0]);
  ta->SetBranchAddress("bmE_mean", bmMean[1]);
  ta->SetBranchAddress("bmE_min", bmMin[1]);
  ta->SetBranchAddress("bmS_mean", bmMean[2]);
  ta->SetBranchAddress("bmS_min", bmMin[2]);
  ta->SetBranchAddress("trA_mean", trMean[0]);
  ta->SetBranchAddress("trA_min", trMin[0]);
  ta->SetBranchAddress("trE_mean", trMean[1]);
  ta->SetBranchAddress("trE_min", trMin[1]);
  ta->SetBranchAddress("trS_mean", trMean[2]);
  ta->SetBranchAddress("trS_min", trMin[2]);

  TH1D *hEvM[3] = {nullptr}, *hEvm[3] = {nullptr};
  TH1D *hBmM[3][2] = {{nullptr}}, *hBmm[3][2] = {{nullptr}};
  TH1D *hTrM[3][3][2] = {{{nullptr}}}, *hTrm[3][3][2] = {{{nullptr}}};
  // for (int k = 0; k < 3; ++k) // Pause a/e disabled.
  for (int k : PeakTimingPlots::ActiveTimingModels)
  {
    hEvM[k] = new TH1D(Form("hT_%s_ev_mean", defNames[k]),
                       Form("%s event-level mean;t [ns];events", defNames[k]),
                       130, -25.0, 625.0);
    hEvm[k] = new TH1D(Form("hT_%s_ev_min", defNames[k]),
                       Form("%s event-level min;t [ns];events", defNames[k]),
                       130, -25.0, 625.0);
    for (int c = 0; c < 2; ++c)
    {
      hBmM[k][c] = new TH1D(Form("hT_%s_bm%d_mean", defNames[k], c + 1),
                            Form("BM%d mean;t [ns];events", c + 1),
                            130, -25.0, 625.0);
      hBmm[k][c] = new TH1D(Form("hT_%s_bm%d_min", defNames[k], c + 1),
                            Form("BM%d min;t [ns];events", c + 1),
                            130, -25.0, 625.0);
    }
    for (int j = 0; j < 3; ++j)
      for (int d = 0; d < 2; ++d)
      {
        hTrM[k][j][d] = new TH1D(
            Form("hT_%s_T%d%s_mean", defNames[k], j, d == 0 ? "X" : "Y"),
            Form("T%d %s mean;t [ns];events", j, d == 0 ? "X" : "Y"),
            130, -25.0, 625.0);
        hTrm[k][j][d] = new TH1D(
            Form("hT_%s_T%d%s_min", defNames[k], j, d == 0 ? "X" : "Y"),
            Form("T%d %s min;t [ns];events", j, d == 0 ? "X" : "Y"),
            130, -25.0, 625.0);
      }
  }

  const Long64_t entries = ta->GetEntries();
  for (Long64_t i = 0; i < entries; ++i)
  {
    ta->GetEntry(i);
    // for (int k = 0; k < 3; ++k) // Pause a/e disabled.
    for (int k : PeakTimingPlots::ActiveTimingModels)
    {
      if (evMean[k] > -900.0)
      {
        hEvM[k]->Fill(evMean[k] * 25.0);
        hEvm[k]->Fill(evMin[k] * 25.0);
      }
      for (int c = 0; c < 2; ++c)
      {
        if (bmMean[k][c] > -900.0)
        {
          hBmM[k][c]->Fill(bmMean[k][c] * 25.0);
          hBmm[k][c]->Fill(bmMin[k][c] * 25.0);
        }
      }
      for (int j = 0; j < 3; ++j)
        for (int d = 0; d < 2; ++d)
        {
          if (trMean[k][j][d] > -900.0)
          {
            hTrM[k][j][d]->Fill(trMean[k][j][d] * 25.0);
            hTrm[k][j][d]->Fill(trMin[k][j][d] * 25.0);
          }
        }
    }
  }

  std::cout << "Track-matched time resolution (RMS over events, ns):"
            << std::endl;
  // for (int k = 0; k < 3; ++k) // Pause a/e disabled.
  for (int k : PeakTimingPlots::ActiveTimingModels)
  {
    TCanvas *c = new TCanvas(Form("cTrack_%s", defNames[k]),
                             Form("track-matched time resolution: %s",
                                  defNames[k]), 1200, 900);
    c->Divide(2, 2);
    c->cd(1);
    hEvM[k]->SetLineColor(colors[k]);
    hEvM[k]->Draw();
    c->cd(2);
    hEvm[k]->SetLineColor(colors[k]);
    hEvm[k]->Draw();
    c->cd(3);
    hBmM[k][0]->SetLineColor(kBlue + 1);
    hBmM[k][1]->SetLineColor(kRed + 1);
    hBmM[k][0]->Draw();
    hBmM[k][1]->Draw("same");
    c->cd(4);
    hBmm[k][0]->SetLineColor(kBlue + 1);
    hBmm[k][1]->SetLineColor(kRed + 1);
    hBmm[k][0]->Draw();
    hBmm[k][1]->Draw("same");
    TString png = Form("%s/track_time_resolution_%s.png", outDir.Data(),
                       defNames[k]);
    c->SaveAs(png);
    std::cout << "saved " << png << std::endl;
    delete c;

    std::cout << "  " << defTitles[k] << ":" << std::endl;
    std::cout << "    event mean: entries=" << hEvM[k]->GetEntries()
              << " mean=" << hEvM[k]->GetMean() << " ns RMS=" << hEvM[k]->GetRMS()
              << " ns | min: entries=" << hEvm[k]->GetEntries()
              << " mean=" << hEvm[k]->GetMean() << " ns RMS=" << hEvm[k]->GetRMS()
              << " ns" << std::endl;
    for (int c = 0; c < 2; ++c)
      std::cout << "    BM" << c + 1 << " mean: entries="
                << hBmM[k][c]->GetEntries() << " RMS=" << hBmM[k][c]->GetRMS()
                << " ns | min: entries=" << hBmm[k][c]->GetEntries()
                << " RMS=" << hBmm[k][c]->GetRMS() << " ns" << std::endl;
    for (int j = 0; j < 3; ++j)
      for (int d = 0; d < 2; ++d)
        std::cout << "    T" << j << (d == 0 ? "X" : "Y") << " mean: entries="
                  << hTrM[k][j][d]->GetEntries() << " RMS="
                  << hTrM[k][j][d]->GetRMS() << " ns | min: entries="
                  << hTrm[k][j][d]->GetEntries() << " RMS="
                  << hTrm[k][j][d]->GetRMS() << " ns" << std::endl;
  }

  // for (int k = 0; k < 3; ++k) // Pause a/e disabled.
  for (int k : PeakTimingPlots::ActiveTimingModels)
  {
    delete hEvM[k];
    delete hEvm[k];
    for (int c = 0; c < 2; ++c)
    {
      delete hBmM[k][c];
      delete hBmm[k][c];
    }
    for (int j = 0; j < 3; ++j)
      for (int d = 0; d < 2; ++d)
      {
        delete hTrM[k][j][d];
        delete hTrm[k][j][d];
      }
  }
  ft->Close();
}

void delta_time_resolution(TString trackFile = "../Combine/Pad-step3-rec.root")
{
  std::unique_ptr<TFile> ft(OpenWithFallback(trackFile));
  if (!ft || ft->IsZombie())
  {
    std::cerr << "ERROR: cannot open track arrival file" << std::endl;
    return;
  }
  TString outDir = TString(gSystem->DirName(__FILE__)) + "/../reports/bm_tracker2_timing";
  gSystem->mkdir(outDir, kTRUE);
  gStyle->SetOptStat(0);
  try
  {
    // One implementation for the normal entry point and the diagnostic run:
    // BOTH BMs, physical Tracker2 X/Y pooled; 4 sigmoid-only mean/min plots.
    PeakTimingPlots::DeltaPlots(outDir, *ft);
    std::cout << "Saved 4 sigmoid BM - Tracker2(X+Y) delta plots to "
              << outDir << "/sigmoid_only/" << std::endl;
  }
  catch (const std::exception &error)
  {
    std::cerr << "ERROR: " << error.what() << std::endl;
  }
}

void time_resolution(TString step3File = "../Combine/Tracker-step3-dec.root",
                     TString padFile = "../Combine/Pad-step3-rec.root")
{
  TString mode = step3File;
  mode.ToLower();
  if (mode == "pad" || mode == "pad_only")
  {
    pad_time_resolution(padFile);
    return;
  }
  if (mode == "track" || mode == "track_only")
  {
    track_time_resolution(padFile);
    return;
  }
  if (mode == "delta" || mode == "delta_only")
  {
    delta_time_resolution(padFile);
    return;
  }

  TFile *f3 = OpenWithFallback(step3File);
  TFile *fp = OpenWithFallback(padFile);
  if (!f3 || f3->IsZombie() || !fp || fp->IsZombie())
  {
    std::cerr << "ERROR: cannot open step3 or pad file" << std::endl;
    return;
  }

  std::vector<Def> defs = {
//      {"a", "xArrivalT", "yArrivalT", "xArrivalT_min", "yArrivalT_min",
//       "arrivalT_a_mean", "arrivalT_a_min"},
//      {"e", "xArrivalTe", "yArrivalTe", "xArrivalTe_min", "yArrivalTe_min",
//       "arrivalT_e_mean", "arrivalT_e_min"},
      {"sigmoid", "xArrivalTSigmoid", "yArrivalTSigmoid",
       "xArrivalTSigmoid_min", "yArrivalTSigmoid_min",
       "arrivalT_sigmoid_mean", "arrivalT_sigmoid_min"}};

  TString macroPath = gInterpreter->GetCurrentMacroName();
  TString outDir = gSystem->DirName(macroPath);
  outDir += "/../reports";
  gSystem->mkdir(outDir, kTRUE);
  gStyle->SetOptStat(0);

  TTree *tr[3] = {nullptr};
  for (int j = 0; j < 3; ++j)
    tr[j] = dynamic_cast<TTree *>(f3->Get(Form("Track%d", j)));
  TTree *pad[2] = {dynamic_cast<TTree *>(fp->Get("BM1Events")),
                   dynamic_cast<TTree *>(fp->Get("BM2Events"))};

  std::cout << "Time resolution summary (RMS of per-event mean/min, ns):" << std::endl;
  std::cout << "def | channel | stat | entries | mean [ns] | RMS [ns]" << std::endl;

  for (const Def &def : defs)
  {
    TCanvas *c = new TCanvas(Form("c_%s", def.name.c_str()),
                             Form("time resolution: %s", def.name.c_str()),
                             1000, 1100);
    c->Divide(2, 5);

    std::vector<TObject *> keepAlive;
    keepAlive.push_back(c);

    for (int j = 0; j < 3; ++j)
    {
      TH1D *hXM = new TH1D(Form("h_%s_T%dX_mean", def.name.c_str(), j),
                           Form("T%d X mean;t [ns]", j), 130, -25.0, 625.0);
      TH1D *hYM = new TH1D(Form("h_%s_T%dY_mean", def.name.c_str(), j),
                           Form("T%d Y mean;t [ns]", j), 130, -25.0, 625.0);
      keepAlive.push_back(hXM);
      keepAlive.push_back(hYM);
      FillHist(hXM, tr[j], def.trMeanX.c_str());
      FillHist(hYM, tr[j], def.trMeanY.c_str());

      c->cd(j + 1);
      hXM->SetLineColor(kAzure + 2);
      hYM->SetLineColor(kOrange + 2);
      hXM->Draw();
      hYM->Draw("same");
      TLegend *leg = new TLegend(0.60, 0.72, 0.88, 0.88);
      keepAlive.push_back(leg);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.06);
      leg->AddEntry(hXM, "X", "l");
      leg->AddEntry(hYM, "Y", "l");
      leg->Draw();

      TH1D *hXm = new TH1D(Form("h_%s_T%dX_min", def.name.c_str(), j),
                           Form("T%d X min;t [ns]", j), 130, -25.0, 625.0);
      TH1D *hYm = new TH1D(Form("h_%s_T%dY_min", def.name.c_str(), j),
                           Form("T%d Y min;t [ns]", j), 130, -25.0, 625.0);
      keepAlive.push_back(hXm);
      keepAlive.push_back(hYm);
      FillHist(hXm, tr[j], def.trMinX.c_str());
      FillHist(hYm, tr[j], def.trMinY.c_str());

      c->cd(5 + j + 1);
      hXm->SetLineColor(kAzure + 2);
      hYm->SetLineColor(kOrange + 2);
      hXm->Draw();
      hYm->Draw("same");
      TLegend *leg2 = new TLegend(0.60, 0.72, 0.88, 0.88);
      keepAlive.push_back(leg2);
      leg2->SetBorderSize(0);
      leg2->SetFillStyle(0);
      leg2->SetTextSize(0.06);
      leg2->AddEntry(hXm, "X", "l");
      leg2->AddEntry(hYm, "Y", "l");
      leg2->Draw();

      std::cout << def.name << " | T" << j << "X | mean | " << hXM->GetEntries()
                << " | " << hXM->GetMean() << " | " << hXM->GetRMS() << std::endl;
      std::cout << def.name << " | T" << j << "Y | mean | " << hYM->GetEntries()
                << " | " << hYM->GetMean() << " | " << hYM->GetRMS() << std::endl;
      std::cout << def.name << " | T" << j << "X | min  | " << hXm->GetEntries()
                << " | " << hXm->GetMean() << " | " << hXm->GetRMS() << std::endl;
      std::cout << def.name << " | T" << j << "Y | min  | " << hYm->GetEntries()
                << " | " << hYm->GetMean() << " | " << hYm->GetRMS() << std::endl;
    }

    if (pad[0] || pad[1])
    {
      for (int bm = 1; bm <= 2; ++bm)
      {
        if (!pad[bm - 1])
          continue;
        TH1D *hPM = new TH1D(Form("h_%s_bm%d_mean", def.name.c_str(), bm),
                             Form("BM%d mean;t [ns]", bm), 130, -25.0, 625.0);
        TH1D *hPm = new TH1D(Form("h_%s_bm%d_min", def.name.c_str(), bm),
                             Form("BM%d min;t [ns]", bm), 130, -25.0, 625.0);
        keepAlive.push_back(hPM);
        keepAlive.push_back(hPm);
        FillHist(hPM, pad[bm - 1], def.padMean.c_str());
        FillHist(hPm, pad[bm - 1], def.padMin.c_str());

        c->cd(3 + bm);
        hPM->SetLineColor(kAzure + 2);
        hPM->Draw();
        c->cd(8 + bm);
        hPm->SetLineColor(kOrange + 2);
        hPm->Draw();

        std::cout << def.name << " | BM" << bm << " | mean | " << hPM->GetEntries()
                  << " | " << hPM->GetMean() << " | " << hPM->GetRMS() << std::endl;
        std::cout << def.name << " | BM" << bm << " | min  | " << hPm->GetEntries()
                  << " | " << hPm->GetMean() << " | " << hPm->GetRMS() << std::endl;
      }
    }

    TString png = Form("%s/time_resolution_%s.png", outDir.Data(), def.name.c_str());
    c->SaveAs(png);
    std::cout << "saved " << png << std::endl;

    for (TObject *obj : keepAlive)
    {
      if (obj != c)
        delete obj;
    }
  }

  f3->Close();
  fp->Close();
}

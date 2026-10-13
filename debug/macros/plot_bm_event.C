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
#include <TMarker.h>
#include <TLatex.h>

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

namespace
{
double PulseModelLocal(double *x, double *par)
{
  const double u = (x[0] - par[2]) / par[3];
  return par[0] + par[1] * u * u * u * std::exp(-u);
}
}

TFile *OpenWithFallback2(TString &name)
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

struct PadHitEntry
{
  int position = -1;
  Long64_t entry = -1;
  double maxA = 0.0;
  double ped = 0.0;
  double pedrms = 1.0;
  double arrivalT = -999.0;
  double sigmoidT20 = -999.0;
};

void DrawPadWaveform(int ev, int r, int cc, int chip, const PadHitEntry &hit,
                     std::vector<short> *wave, int maxT, int sigmoidFitStatus,
                     int sigmoidFitStart, int sigmoidFitEnd, double sigmoidH,
                     double sigmoidT50, double sigmoidK)
{
  if (!wave || wave->size() < 4)
    return;

  const int nSamples = (int)wave->size();
  const double thr = 5.0 * std::max(1.0, std::abs(hit.pedrms));

  TGraph *gW = new TGraph();
  gW->SetMarkerStyle(20);
  gW->SetMarkerSize(0.5);
  double ymin = 1e9;
  double ymax = -1e9;
  for (int k = 0; k < nSamples; ++k)
  {
    const double yv = (double)wave->at(k) - hit.ped;
    gW->SetPoint(k, (double)k, yv);
    ymin = std::min(ymin, yv);
    ymax = std::max(ymax, yv);
  }
  ymin = std::min(ymin, -thr);
  ymax = std::max(ymax, thr);
  gW->SetMinimum(ymin * 1.25);
  gW->SetMaximum(ymax * 1.25);
  gW->SetTitle(Form("pos %d;sample;ADC-ped", hit.position));
  gW->Draw("ALP");

  TLine *l0 = new TLine(0, 0, nSamples - 1, 0);
  l0->SetLineColor(kGray + 2);
  l0->SetLineStyle(2);
  l0->Draw();
  TLine *lu = new TLine(0, thr, nSamples - 1, thr);
  lu->SetLineColor(kGray + 1);
  lu->SetLineStyle(7);
  lu->Draw();
  TLine *ld = new TLine(0, -thr, nSamples - 1, -thr);
  ld->SetLineColor(kGray + 1);
  ld->SetLineStyle(7);
  ld->Draw();

  if (sigmoidFitStatus == 0 && sigmoidH > 0.0 && sigmoidFitStart <= sigmoidFitEnd)
  {
    const int sigEnd = std::max(sigmoidFitStart, std::min(nSamples - 1, maxT + 4));
    TF1 *fSig = new TF1(Form("sig_ev%d_%d_%d_%d", ev, r, cc, chip),
                        "[0]/(1.0+exp(-(x-[1])/[2]))", sigmoidFitStart, sigEnd);
    fSig->SetParameters(sigmoidH, sigmoidT50, sigmoidK);
    fSig->SetLineColor(kOrange + 2);
    fSig->SetLineWidth(2);
    fSig->Draw("same");
  }

  if (hit.maxA > 0 && maxT > 0 && maxT < nSamples)
  {
    const double thresh = std::max(3.0 * std::max(1.0, std::abs(hit.pedrms)),
                                   0.1 * hit.maxA);
    int start = 0;
    for (int k = maxT; k > 0; --k)
    {
      const double adc = (double)wave->at(k) - hit.ped;
      const double prev = (double)wave->at(k - 1) - hit.ped;
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
      TF1 *fGam = new TF1(Form("gam_ev%d_%d_%d_%d", ev, r, cc, chip),
                          PulseModelLocal, 0, nSamples - 1, 4);
      const double b0 = std::max(0.5, ((double)maxT - start) / 3.0);
      const double a0 = (double)maxT - 3.0 * b0;
      const double A0 = std::max(1.0, hit.maxA * std::exp(3.0) / 27.0);
      fGam->FixParameter(0, 0.0);
      fGam->SetParameters(0.0, A0, a0, b0);
      fGam->SetParLimits(1, 0.0, std::max(100.0, 10.0 * hit.maxA));
      fGam->SetParLimits(2, start - 10.0, (double)end);
      fGam->SetParLimits(3, 0.1, 20.0);
      TGraph *gFit = new TGraph();
      int np = 0;
      for (int k = start; k <= end; ++k)
        gFit->SetPoint(np++, (double)k, (double)wave->at(k) - hit.ped);
      if (gFit->Fit(fGam, "QSN") == 0)
      {
        fGam->SetLineColor(kAzure + 2);
        fGam->SetLineWidth(2);
        fGam->Draw("same");
      }
      delete gFit;
    }
  }
}

void plot_bm_event(TString step2File = "../Combine/Tracker-step2-ana.root",
                   TString padFile = "../Combine/Pad-step3-rec.root",
                   int eventWanted = -1,
                   int maxEvents = 20)
{
  TFile *f2 = OpenWithFallback2(step2File);
  TFile *fp = OpenWithFallback2(padFile);
  if (!f2 || f2->IsZombie() || !fp || fp->IsZombie())
  {
    std::cerr << "ERROR: cannot open step2 or pad file" << std::endl;
    return;
  }

  TTree *bm = dynamic_cast<TTree *>(f2->Get("BMPadSignals"));
  TTree *pad[2] = {dynamic_cast<TTree *>(fp->Get("BM1Events")),
                   dynamic_cast<TTree *>(fp->Get("BM2Events"))};
  if (!bm || !pad[0] || !pad[1])
  {
    std::cerr << "ERROR: BMPadSignals or BM1Events/BM2Events tree missing" << std::endl;
    return;
  }

  int event = -1;
  int rawChip = -1;
  int position = -1;
  int maxT = -1;
  int sigmoidFitStatus = -1;
  int sigmoidFitStart = -1;
  int sigmoidFitEnd = -1;
  double ped = 0.0;
  double pedrms = 1.0;
  double maxA = 0.0;
  double arrivalT = -999.0;
  double sigmoidT20 = -999.0;
  double sigmoidH = -999.0;
  double sigmoidT50 = -999.0;
  double sigmoidK = -999.0;
  std::vector<short> *wave = nullptr;
  bm->SetBranchAddress("event", &event);
  bm->SetBranchAddress("raw_chip", &rawChip);
  bm->SetBranchAddress("position", &position);
  bm->SetBranchAddress("maxT", &maxT);
  bm->SetBranchAddress("ped", &ped);
  bm->SetBranchAddress("pedrms", &pedrms);
  bm->SetBranchAddress("maxA", &maxA);
  bm->SetBranchAddress("arrivalT", &arrivalT);
  bm->SetBranchAddress("arrivalT_sigmoid", &sigmoidT20);
  bm->SetBranchAddress("sigmoidFitStatus", &sigmoidFitStatus);
  bm->SetBranchAddress("sigmoid_fitStart", &sigmoidFitStart);
  bm->SetBranchAddress("sigmoid_fitEnd", &sigmoidFitEnd);
  bm->SetBranchAddress("sigmoid_H", &sigmoidH);
  bm->SetBranchAddress("sigmoid_t50", &sigmoidT50);
  bm->SetBranchAddress("sigmoid_k", &sigmoidK);
  bm->SetBranchAddress("wave", &wave);

  std::map<std::pair<int, int>, std::vector<Long64_t>> bmIndex;
  for (Long64_t i = 0; i < bm->GetEntries(); ++i)
  {
    bm->GetEntry(i);
    if (position >= 0 && position < 121)
      bmIndex[std::make_pair(event, rawChip)].push_back(i);
  }

  int pev[2] = {-1, -1};
  std::vector<double> *clusterX[2] = {nullptr, nullptr};
  std::vector<double> *clusterY[2] = {nullptr, nullptr};
  double primaryX[2] = {-1.0, -1.0};
  double primaryY[2] = {-1.0, -1.0};
  std::map<int, Long64_t> padIndex[2];
  for (int ci = 0; ci < 2; ++ci)
  {
    pad[ci]->SetBranchAddress("event", &pev[ci]);
    pad[ci]->SetBranchAddress("cluster_x", &clusterX[ci]);
    pad[ci]->SetBranchAddress("cluster_y", &clusterY[ci]);
    pad[ci]->SetBranchAddress("primary_x", &primaryX[ci]);
    pad[ci]->SetBranchAddress("primary_y", &primaryY[ci]);
    for (Long64_t i = 0; i < pad[ci]->GetEntries(); ++i)
    {
      pad[ci]->GetEntry(i);
      padIndex[ci][pev[ci]] = i;
    }
  }

  std::vector<int> eventList;
  if (eventWanted >= 0)
    eventList.push_back(eventWanted);
  else
  {
    for (const auto &kv : bmIndex)
      eventList.push_back(kv.first.first);
    std::sort(eventList.begin(), eventList.end());
    eventList.erase(std::unique(eventList.begin(), eventList.end()),
                    eventList.end());
  }

  TString macroPath = gInterpreter->GetCurrentMacroName();
  TString outDir = gSystem->DirName(macroPath);
  outDir += "/../reports";
  gSystem->mkdir(outDir, kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);

  const int chips[2] = {12, 13};
  const char *chipNames[2] = {"BM1 (chip12)", "BM2 (chip13)"};
  int plotted = 0;

  for (int ev : eventList)
  {
    if (plotted >= maxEvents && eventWanted < 0)
      break;

    TCanvas *c = new TCanvas(Form("cBMEv%d", ev), Form("BM event %d", ev),
                             1500, 900);
    c->Divide(1, 2);

    for (int ci = 0; ci < 2; ++ci)
    {
      const int chip = chips[ci];
      auto it = bmIndex.find(std::make_pair(ev, chip));
      const int nfired = it != bmIndex.end() ? (int)it->second.size() : 0;

      c->cd(ci + 1);
      TPad *mapPad = new TPad(Form("map_%d_%d", ci, ev), "hitmap",
                              0.02, 0.03, 0.30, 0.97);
      TPad *wavePad = new TPad(Form("wave_%d_%d", ci, ev), "waves",
                               0.34, 0.03, 0.98, 0.97);
      mapPad->Draw();
      wavePad->Draw();

      mapPad->cd();
      mapPad->SetMargin(0.10, 0.10, 0.10, 0.05);
      TH2D *hHit = new TH2D(Form("hhit_%d_%d", ci, ev),
                            Form("%s, event %d;ix;iy", chipNames[ci], ev),
                            11, -0.5, 10.5, 11, -0.5, 10.5);
      hHit->SetStats(0);
      hHit->GetZaxis()->SetTitle("maxA [ADC]");

      if (it != bmIndex.end())
      {
        for (Long64_t ei : it->second)
        {
          bm->GetEntry(ei);
          hHit->Fill(position % 11, position / 11, maxA);
        }
      }
      hHit->Draw("COLZ");

      auto pit = padIndex[ci].find(ev);
      if (pit != padIndex[ci].end())
      {
        pad[ci]->GetEntry(pit->second);
        if (clusterX[ci] && clusterY[ci])
        {
          for (size_t k = 0; k < clusterX[ci]->size(); ++k)
          {
            TMarker *mk = new TMarker(((*clusterX[ci])[k] - 4.5) / 9.0,
                                      ((*clusterY[ci])[k] - 4.5) / 9.0, 24);
            mk->SetMarkerColor(kMagenta);
            mk->SetMarkerSize(1.2);
            mk->Draw();
          }
        }
        if (primaryX[ci] >= 0.0 && primaryY[ci] >= 0.0)
        {
          TMarker *mp = new TMarker((primaryX[ci] - 4.5) / 9.0,
                                    (primaryY[ci] - 4.5) / 9.0, 29);
          mp->SetMarkerColor(kRed + 1);
          mp->SetMarkerSize(2.2);
          mp->Draw();
        }
      }

      wavePad->cd();
      if (nfired == 0)
      {
        TLatex *ln = new TLatex(0.5, 0.5, "no hits");
        ln->SetNDC();
        ln->SetTextSize(0.08);
        ln->SetTextAlign(22);
        ln->Draw();
        continue;
      }

      const int ncols = std::min(5, nfired);
      const int nrows = (nfired + ncols - 1) / ncols;
      const double cw = 0.96 / ncols;
      const double chh = 0.94 / nrows;

      std::vector<PadHitEntry> hits;
      for (int k = 0; k < nfired; ++k)
      {
        bm->GetEntry(it->second[k]);
        PadHitEntry h;
        h.position = position;
        h.entry = it->second[k];
        h.maxA = maxA;
        h.ped = ped;
        h.pedrms = pedrms;
        h.arrivalT = arrivalT;
        h.sigmoidT20 = sigmoidT20;
        hits.push_back(h);
      }
      std::sort(hits.begin(), hits.end(),
                [](const PadHitEntry &a, const PadHitEntry &b)
                { return a.position < b.position; });

      std::vector<TPad *> cells;
      for (int k = 0; k < nfired; ++k)
      {
        const int r = k / ncols;
        const int cc = k % ncols;
        const double x1 = 0.02 + cc * cw;
        const double y1 = 1.0 - (r + 1) * chh - 0.02;
        TPad *cell = new TPad(Form("cell_%d_%d_%d", ci, k, ev),
                              Form("cell_%d_%d_%d", ci, k, ev),
                              x1, y1, x1 + cw * 0.98, y1 + chh * 0.96);
        cell->Draw();
        cells.push_back(cell);
      }

      for (int k = 0; k < nfired; ++k)
      {
        cells[k]->cd();
        bm->GetEntry(hits[k].entry);
        DrawPadWaveform(ev, k / ncols, k % ncols, chip, hits[k], wave, maxT,
                        sigmoidFitStatus, sigmoidFitStart, sigmoidFitEnd,
                        sigmoidH, sigmoidT50, sigmoidK);
      }

      std::cout << "event " << ev << " " << chipNames[ci] << ": "
                << nfired << " fired pads" << std::endl;
      for (const PadHitEntry &h : hits)
      {
        std::cout << "  pos " << h.position << " maxA=" << h.maxA
                  << " pedrms=" << h.pedrms << " a=" << h.arrivalT
                  << " t20=" << h.sigmoidT20 << std::endl;
      }
    }

    TString png = Form("%s/bm_event_ev%d.png", outDir.Data(), ev);
    c->SaveAs(png);
    std::cout << "saved " << png << std::endl;

    delete c;
    ++plotted;
  }

  f2->Close();
  fp->Close();
}

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TStyle.h>
#include <TH1.h>
#include "../../BTTracker/include/CRInfo.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

namespace
{
TString noiseConfigDirectory;
TString ConfigPath(const char *name)
{
  if (!noiseConfigDirectory.IsNull())
    return noiseConfigDirectory+"/"+name;
  TString p = TString("../config/") + name;
  std::ifstream test(p.Data());
  if (test.is_open())
    return p;
  return TString("config/") + name;
}

bool MapGlobalChannel(int globalChannel, int &logicalBoard, int &chip,
                      int &rawChannel, TString &description)
{
  logicalBoard = -1;
  chip = -1;
  rawChannel = -1;
  description = "";

  if (globalChannel >= 1 && globalChannel <= trackerGlobalChannelCount)
  {
    std::ifstream in(ConfigPath("map_APV.txt"));
    if (!in.is_open())
    {
      std::cerr << "ERROR: cannot open map_APV.txt" << std::endl;
      return false;
    }
    int board, chipL, channel, connector, strip;
    while (in >> board >> chipL >> channel >> connector >> strip)
    {
      const int treeIndex = board * 2 + chipL / 2;
      const int g = treeIndex * 256 + strip + 1;
      if (g == globalChannel)
      {
        logicalBoard = board;
        chip = chipL;
        rawChannel = channel;
        description = Form("board%d chip%d ch%d (strip %d)",
                           board, chipL, channel, strip);
        return true;
      }
    }
    std::cerr << "ERROR: global channel " << globalChannel
              << " not found in map_APV.txt" << std::endl;
    return false;
  }
  if (globalChannel > trackerGlobalChannelCount && globalChannel <= globalNoiseChannelCount)
  {
    std::ifstream in(ConfigPath("map_APV4BMM.txt"));
    if (!in.is_open())
    {
      std::cerr << "ERROR: cannot open map_APV4BMM.txt" << std::endl;
      return false;
    }
    int padBoard, chipM, channelId, connectorPin, position;
    while (in >> padBoard >> chipM >> channelId >> connectorPin >> position)
    {
      const int g = trackerGlobalChannelCount + padBoard * bmChannelsPerDetector + position + 1;
      if (g == globalChannel)
      {
        logicalBoard = -1;
        chip = 12 + padBoard;
        rawChannel = channelId - 1;
        description = Form("chip%d ch%d (pad position %d)",
                           chip, rawChannel, position);
        return true;
      }
    }
    std::cerr << "ERROR: global channel " << globalChannel
              << " not found in map_APV4BMM.txt" << std::endl;
    return false;
  }
  std::cerr << "ERROR: global channel must be in 1.." << globalNoiseChannelCount << std::endl;
  return false;
}
}

void plot_noise_channel(TString noiseRawFile = "Combine/noise/Tracker-noise-raw.root",
                        int globalChannel = trackerGlobalChannelCount + 1,
                        int maxWaves = -1,
                        TString configDirectory = "", TString outputDirectory = "")
{
  noiseConfigDirectory = configDirectory;
  int logicalBoard = -1;
  int chip = -1;
  int rawChannel = -1;
  TString description;
  if (!MapGlobalChannel(globalChannel, logicalBoard, chip, rawChannel,
                        description))
    return;

  TFile *f = TFile::Open(noiseRawFile);
  if ((!f || f->IsZombie()) && !noiseRawFile.BeginsWith("/"))
  {
    TString up = "../" + noiseRawFile;
    f = TFile::Open(up);
    if (f && !f->IsZombie())
      noiseRawFile = up;
  }
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open " << noiseRawFile << std::endl;
    return;
  }
  TTree *t = dynamic_cast<TTree *>(f->Get("TrackRaw"));
  if (!t)
  {
    std::cerr << "ERROR: TrackRaw tree missing in " << noiseRawFile << std::endl;
    f->Close();
    return;
  }

  int event = -1;
  int rb = -1;
  int ch = -1;
  int chan = -1;
  std::vector<short> *wave = nullptr;
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("board", &rb);
  t->SetBranchAddress("chip", &ch);
  t->SetBranchAddress("channel", &chan);
  t->SetBranchAddress("wave", &wave);

  std::vector<std::pair<int, std::vector<short>>> waves;
  for (Long64_t i = 0; i < t->GetEntries(); ++i)
  {
    t->GetEntry(i);
    if (rb == logicalBoard && ch == chip && chan == rawChannel && wave)
      waves.push_back(std::make_pair(event, *wave));
  }

  if (waves.empty())
  {
    std::cout << "No waveforms found for global channel " << globalChannel
              << " (" << description << ")" << std::endl;
    f->Close();
    return;
  }

  double sum = 0.0;
  double sum2 = 0.0;
  long nSamples = 0;
  for (const auto &w : waves)
  {
    for (short s : w.second)
    {
      sum += s;
      sum2 += (double)s * s;
      ++nSamples;
    }
  }
  const double mean = sum / nSamples;
  const double sigma = std::sqrt(std::max(0.0, sum2 / nSamples - mean * mean));

  const int total = maxWaves > 0 ? std::min(maxWaves, (int)waves.size())
                                 : (int)waves.size();

  // Keep diagnostics with this run rather than a shared debug/reports tree.
  TString outDir = outputDirectory;
  if (outDir.IsNull())
    outDir = TString(gSystem->DirName(noiseRawFile))+"/waveforms/gc"+
             TString::Format("%d",globalChannel);
  gSystem->mkdir(outDir, kTRUE);

  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);

  TCanvas *c = new TCanvas("cNoiseWave", "noise waveform", 800, 500);
  c->SetMargin(0.12, 0.05, 0.12, 0.05);

  for (int k = 0; k < total; ++k)
  {
    c->Clear();

    const std::vector<short> &w = waves[k].second;
    const int ev = waves[k].first;
    const int n = (int)w.size();

    double ymin = 1e9;
    double ymax = -1e9;
    TGraph g;
    for (int i = 0; i < n; ++i)
    {
      g.SetPoint(i, (double)i, (double)w[i]);
      ymin = std::min(ymin, (double)w[i]);
      ymax = std::max(ymax, (double)w[i]);
    }
    ymin = std::min(ymin, mean - sigma);
    ymax = std::max(ymax, mean + sigma);

    g.SetMarkerStyle(20);
    g.SetMarkerSize(0.7);
    g.SetMarkerColor(kBlack);
    g.SetLineColor(kBlack);
    g.SetTitle(Form("gc%d (%s) wave %d/%d event %d;sample;ADC",
                    globalChannel, description.Data(), k + 1, total, ev));
    g.SetMinimum(ymin - 0.1 * (ymax - ymin + 1.0));
    g.SetMaximum(ymax + 0.1 * (ymax - ymin + 1.0));
    g.Draw("ALP");

    TLine lm(g.GetXaxis()->GetXmin(), mean, g.GetXaxis()->GetXmax(), mean);
    lm.SetLineColor(kAzure + 2);
    lm.SetLineWidth(2);
    lm.Draw();

    TLine lu(g.GetXaxis()->GetXmin(), mean + sigma,
             g.GetXaxis()->GetXmax(), mean + sigma);
    lu.SetLineColor(kGray + 1);
    lu.SetLineStyle(7);
    lu.Draw();
    TLine ld(g.GetXaxis()->GetXmin(), mean - sigma,
             g.GetXaxis()->GetXmax(), mean - sigma);
    ld.SetLineColor(kGray + 1);
    ld.SetLineStyle(7);
    ld.Draw();

    TLatex info(0.62, 0.84,
                Form("#mu = %.1f, #sigma = %.1f", mean, sigma));
    info.SetNDC();
    info.SetTextColor(kAzure + 2);
    info.SetTextSize(0.045);
    info.Draw();

    TString png = Form("%s/wave_%05d_ev%d.png", outDir.Data(), k + 1, ev);
    c->SaveAs(png);

    if ((k + 1) % 100 == 0 || k + 1 == total)
      std::cout << "saved " << k + 1 << "/" << total << std::endl;
  }

  delete c;
  f->Close();

  std::cout << "global channel " << globalChannel << " -> "
            << description << "\n"
            << "  logicalBoard=" << logicalBoard << " chip=" << chip
            << " rawChannel=" << rawChannel << "\n"
            << "  waveforms=" << waves.size()
            << (maxWaves > 0 ? Form(" (output first %d)", total) : "")
            << "\n"
            << "  pooled mean=" << mean << " sigma=" << sigma << "\n"
            << "  output: " << outDir << std::endl;
}

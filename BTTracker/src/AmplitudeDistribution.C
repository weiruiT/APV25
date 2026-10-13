#include "CRAnalysis.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>

void MakeAmplitudeDistributions(TString fileDir)
{
  TString anaName;
  TString outputDir;
  if (fileDir.EndsWith(".root"))
  {
    anaName = fileDir;
    outputDir = ".";
  }
  else
  {
    if (fileDir.EndsWith("/"))
      fileDir.Remove(fileDir.Length() - 1, 1);
    anaName = fileDir + "/Combine/Tracker-step2-ana.root";
    outputDir = fileDir + "/Combine";
  }

  TFile input(anaName, "READ");
  if (input.IsZombie())
  {
    cerr << "ERROR: cannot open analysis file " << anaName << endl;
    return;
  }

  TString outputName = outputDir + "/Tracker-amplitude.root";
  TString summaryName = outputDir + "/Tracker-amplitude-summary.tsv";
  TFile output(outputName, "RECREATE");
  if (output.IsZombie())
  {
    cerr << "ERROR: cannot write " << outputName << endl;
    return;
  }

  ofstream summary(summaryName.Data());
  if (!summary.is_open())
  {
    cerr << "ERROR: cannot write " << summaryName << endl;
    return;
  }
  summary << "tracker\taxis\tglobal_channel\tdecoded_channel\tentries\tmean_adc\trms_adc\tmax_adc\n";
  summary << fixed << setprecision(6);

  TCanvas allCanvas("amplitude_all", "Signal amplitude distributions", 1400, 1400);
  allCanvas.Divide(2, MMdetN);
  TCanvas mapCanvas("amplitude_map", "Signal amplitude by channel", 1400, 1400);
  mapCanvas.Divide(2, MMdetN);

  TH1D hTotal("amplitude_all_trackers", "All selected tracker amplitudes;Signal amplitude [ADC];Counts",
              512, 0, 4096);
  long long totalEntries = 0;

  for (int treeIndex = 0; treeIndex < 2 * MMdetN; ++treeIndex)
  {
    const int tracker = treeIndex / 2 + 1;
    const bool isX = (treeIndex % 2 == 0);
    const char *axis = isX ? "X" : "Y";
    const TString treeName = Form("Track%s%d", axis, tracker - 1);
    TTree *tree = dynamic_cast<TTree *>(input.Get(treeName));
    if (!tree)
    {
      cerr << "WARNING: missing tree " << treeName << endl;
      continue;
    }

    double maxA = 0.0;
    int channel = -1;
    tree->SetBranchAddress("maxA", &maxA);
    tree->SetBranchAddress("channel", &channel);

    TH1D hAll(Form("amplitude_T%d%s_all", tracker, axis),
              Form("Tracker%d %s signal amplitude;Signal amplitude [ADC];Counts", tracker, axis),
              512, 0, 4096);
    TH2D hMap(Form("amplitude_T%d%s_by_channel", tracker, axis),
              Form("Tracker%d %s amplitude by channel;Decoded channel;Signal amplitude [ADC]", tracker, axis),
              256, 0, 256, 512, 0, 4096);

    vector<TH1D *> perChannel(256, nullptr);
    vector<double> maxObserved(256, -numeric_limits<double>::infinity());
    for (int ch = 0; ch < 256; ++ch)
    {
      perChannel[ch] = new TH1D(
          Form("amplitude_T%d%s_ch%03d", tracker, axis, ch),
          Form("Tracker%d %s channel %d amplitude;Signal amplitude [ADC];Counts", tracker, axis, ch),
          512, 0, 4096);
    }

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
      tree->GetEntry(i);
      if (channel < 0 || channel >= 256 || !std::isfinite(maxA))
        continue;

      hAll.Fill(maxA);
      hMap.Fill(channel, maxA);
      perChannel[channel]->Fill(maxA);
      if (maxA > maxObserved[channel])
        maxObserved[channel] = maxA;
      hTotal.Fill(maxA);
      ++totalEntries;
    }

    output.cd();
    hAll.Write();
    hMap.Write();
    for (int ch = 0; ch < 256; ++ch)
    {
      perChannel[ch]->Write();
      const double observedMax = perChannel[ch]->GetEntries() > 0 ? maxObserved[ch] : -1.0;
      summary << tracker << "\t" << axis << "\t"
              << treeIndex * 256 + ch + 1 << "\t" << ch << "\t"
              << perChannel[ch]->GetEntries() << "\t"
              << perChannel[ch]->GetMean() << "\t"
              << perChannel[ch]->GetRMS() << "\t"
              << observedMax << "\n";
      delete perChannel[ch];
    }

    cout << "  Tracker" << tracker << axis << ": entries=" << hAll.GetEntries() << ", mean=" << hAll.GetMean() << ", rms=" << hAll.GetRMS() << endl;

    allCanvas.cd(treeIndex + 1);
    hAll.SetLineColor(kBlue + 1);
    hAll.Draw();
    mapCanvas.cd(treeIndex + 1);
    hMap.Draw("COLZ");
  }

  output.cd();
  hTotal.Write();
  summary.close();

  allCanvas.SaveAs(outputDir + "/Tracker-amplitude-distributions.png");
  mapCanvas.SaveAs(outputDir + "/Tracker-amplitude-by-channel.png");
  output.Close();
  input.Close();

  cout << "Amplitude distributions finished." << endl;
  cout << "  selected waveform entries: " << totalEntries << endl;
  cout << "  ROOT output: " << outputName << endl;
  cout << "  TSV summary: " << summaryName << endl;
  cout << "  plots: " << outputDir << "/Tracker-amplitude-distributions.png" << endl;
  cout << "  plots: " << outputDir << "/Tracker-amplitude-by-channel.png" << endl;
}


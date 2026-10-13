#include "CRAnalysis.h"

#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <fstream>
#include <iomanip>

namespace
{
double FitMPV(TH1D *hist, const TString &name, double fitMax)
{
  if (!hist || hist->GetEntries() < 20)
    return -1.0;

  TF1 fit(name, "landau", 0.0, fitMax);
  const int status = hist->Fit(&fit, "Q0");
  return status == 0 ? fit.GetParameter(1) : -1.0;
}
}

void MakeEnergySpectra(TString fileDir)
{
  TString decName;
  TString outputDir;
  if (fileDir.EndsWith(".root"))
  {
    decName = fileDir;
    outputDir = ".";
  }
  else
  {
    if (fileDir.EndsWith("/"))
      fileDir.Remove(fileDir.Length() - 1, 1);
    decName = fileDir + "/Combine/Tracker-step3-dec.root";
    outputDir = fileDir + "/Combine";
  }

  TFile input(decName, "READ");
  if (input.IsZombie())
  {
    cerr << "ERROR: cannot open decoded file " << decName << endl;
    return;
  }

  TString outputName = outputDir + "/Tracker-energy-spectrum-cut.root";
  TString summaryName = outputDir + "/Tracker-energy-spectrum-cut-summary.tsv";
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
  summary << "tracker\taxis\tentries\tmean_adc\trms_adc\tlandau_mpv\n";
  summary << fixed << setprecision(6);

  TH1D *hX[MMdetN];
  TH1D *hY[MMdetN];
  TCanvas canvas("energy_spectrum_cut", "Cut signal energy spectra", 1600, 800);
  canvas.Divide(MMdetN, 2);

  for (int tracker = 0; tracker < MMdetN; ++tracker)
  {
    hX[tracker] = new TH1D(Form("energy_T%dX", tracker + 1),
                            Form("Tracker%d X cluster energy;Cluster charge [ADC];Counts", tracker + 1),
                            300, 0, 15000);
    hY[tracker] = new TH1D(Form("energy_T%dY", tracker + 1),
                            Form("Tracker%d Y cluster energy;Cluster charge [ADC];Counts", tracker + 1),
                            300, 0, 15000);

    TTree *tree = dynamic_cast<TTree *>(input.Get(Form("Track%d", tracker)));
    if (!tree)
    {
      cerr << "WARNING: missing Track" << tracker << endl;
      continue;
    }

    bool sigX = false;
    bool sigY = false;
    double hitAmpX = 0.0;
    double hitAmpY = 0.0;
    tree->SetBranchAddress("sig_x", &sigX);
    tree->SetBranchAddress("sig_y", &sigY);
    tree->SetBranchAddress("hit_amp_x", &hitAmpX);
    tree->SetBranchAddress("hit_amp_y", &hitAmpY);

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry)
    {
      tree->GetEntry(entry);
      if (sigX && hitAmpX > 0.0)
        hX[tracker]->Fill(hitAmpX);
      if (sigY && hitAmpY > 0.0)
        hY[tracker]->Fill(hitAmpY);
    }

    const double mpvX = FitMPV(hX[tracker], Form("landau_T%dX", tracker + 1), 15000.0);
    const double mpvY = FitMPV(hY[tracker], Form("landau_T%dY", tracker + 1), 15000.0);

    summary << tracker + 1 << "\tX\t" << hX[tracker]->GetEntries() << "\t"
            << hX[tracker]->GetMean() << "\t" << hX[tracker]->GetRMS() << "\t" << mpvX << "\n";
    summary << tracker + 1 << "\tY\t" << hY[tracker]->GetEntries() << "\t"
            << hY[tracker]->GetMean() << "\t" << hY[tracker]->GetRMS() << "\t" << mpvY << "\n";

    cout << "Tracker" << tracker + 1 << " X: entries=" << hX[tracker]->GetEntries()
         << ", mean=" << hX[tracker]->GetMean() << ", rms=" << hX[tracker]->GetRMS()
         << ", MPV=" << mpvX << endl;
    cout << "Tracker" << tracker + 1 << " Y: entries=" << hY[tracker]->GetEntries()
         << ", mean=" << hY[tracker]->GetMean() << ", rms=" << hY[tracker]->GetRMS()
         << ", MPV=" << mpvY << endl;

    canvas.cd(tracker + 1);
    hX[tracker]->Draw();
    canvas.cd(MMdetN + tracker + 1);
    hY[tracker]->Draw();

    output.cd();
    hX[tracker]->Write();
    hY[tracker]->Write();
  }

  summary.close();
  output.cd();
  canvas.Write();
  canvas.SaveAs(outputDir + "/Tracker-energy-spectrum-cut.png");
  output.Close();
  input.Close();

  cout << "Cut X/Y energy spectra finished." << endl;
  cout << "  ROOT output: " << outputName << endl;
  cout << "  TSV summary: " << summaryName << endl;
  cout << "  PNG output: " << outputDir << "/Tracker-energy-spectrum-cut.png" << endl;
}

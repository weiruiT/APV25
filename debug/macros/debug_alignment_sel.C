#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>

void debug_alignment_sel()
{
  TFile *f = TFile::Open("/home/thoriatong/workplace/BTTracker/Combine/Tracker-step3-dec.root");
  TTree *tr[4];
  bool sig[4] = {false};
  double x[4] = {0}, y[4] = {0}, z[4] = {0};
  for (int j = 0; j < 4; j++)
  {
    tr[j] = (TTree *)f->Get(Form("Track%d", j));
    tr[j]->SetBranchAddress("sig", &sig[j]);
    tr[j]->SetBranchAddress("x", &x[j]);
    tr[j]->SetBranchAddress("y", &y[j]);
    tr[j]->SetBranchAddress("z", &z[j]);
  }

  const double centre = 256.0 / 2 * 0.4;
  const int MMdetN = 4;
  const double alg_fit_R2_cut = 0.05;
  const double alg_theta_cut = 20;
  const double alg_fit_res_cut = 5;

  int n = tr[0]->GetEntries();
  int limit = n / 5;
  int shown = 0;
  for (int entry = 0; entry < limit && shown < 12; entry++)
  {
    for (int j = 0; j < 4; j++)
      tr[j]->GetEntry(entry);
    bool is_good = true;
    for (int j = 0; j < 4; j++)
      is_good = is_good && sig[j] && x[j] > 25 && y[j] > 25 && x[j] < 225 && y[j] < 225;
    if (!is_good)
      continue;
    shown++;

    double x0[4], y0[4], z0[4];
    for (int i = 0; i < 4; i++)
    {
      x0[i] = x[i] * 0.4 - centre;
      y0[i] = y[i] * 0.4 - centre;
      z0[i] = z[i];
    }

    std::cout << "entry " << entry << " x0={" << x0[0] << "," << x0[1] << "," << x0[2] << "," << x0[3]
              << "} y0={" << y0[0] << "," << y0[1] << "," << y0[2] << "," << y0[3]
              << "} z0={" << z0[0] << "," << z0[1] << "," << z0[2] << "," << z0[3] << "}" << std::endl;

    bool good_event = true;
    int N = MMdetN - 1;
    for (int i = 0; i < MMdetN; i++)
    {
      double Mx = 0, My = 0, Mz = 0, Mxx = 0, Myy = 0, Mzz = 0, Mzx = 0, Mzy = 0;
      for (int j = 0; j < MMdetN; j++)
      {
        if (i == j)
          continue;
        Mx += x0[j] / N;
        My += y0[j] / N;
        Mz += z0[j] / N;
        Mxx += x0[j] * x0[j] / N;
        Myy += y0[j] * y0[j] / N;
        Mzz += z0[j] * z0[j] / N;
        Mzx += z0[j] * x0[j] / N;
        Mzy += z0[j] * y0[j] / N;
      }
      const double denomz = Mzz - Mz * Mz;
      const double Kx = (Mzx - Mz * Mx) / denomz;
      const double Bx = Mx - Kx * Mz;
      const double Rx = (Mzx - Mz * Mx) * (Mzx - Mz * Mx) / (Mxx - Mx * Mx) / denomz;
      const double Ky = (Mzy - Mz * My) / denomz;
      const double By = My - Ky * Mz;
      const double Ry = (Mzy - Mz * My) * (Mzy - Mz * My) / (Myy - My * My) / denomz;
      const double thx = TMath::Abs(TMath::ATan(Kx) / TMath::Pi() * 180);
      const double thy = TMath::Abs(TMath::ATan(Ky) / TMath::Pi() * 180);
      const double rx = TMath::Abs(Kx * z0[i] + Bx - x0[i]);
      const double ry = TMath::Abs(Ky * z0[i] + By - y0[i]);
      bool okx = Rx > alg_fit_R2_cut && thx < alg_theta_cut && rx < alg_fit_res_cut;
      bool oky = Ry > alg_fit_R2_cut && thy < alg_theta_cut && ry < alg_fit_res_cut;
      std::cout << "  i=" << i << " Kx=" << Kx << " Rx=" << Rx << " thx=" << thx << " rx=" << rx
                << " | Ky=" << Ky << " Ry=" << Ry << " thy=" << thy << " ry=" << ry
                << " -> " << (okx ? "Xok" : "XFAIL") << " " << (oky ? "Yok" : "YFAIL") << std::endl;
      good_event = good_event && okx && oky;
    }
    std::cout << "  event good = " << good_event << std::endl;
  }
  f->Close();
}

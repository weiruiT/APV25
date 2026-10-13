#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>

void debug_r2_dist()
{
  TFile *f = TFile::Open("/home/thoriatong/workplace/BTTracker/Combine/Tracker-step3-dec.root");
  TTree *tr[4];
  bool sx[4] = {false}, sy[4] = {false};
  double x[4] = {0}, y[4] = {0};
  for (int j = 0; j < 4; j++)
  {
    tr[j] = (TTree *)f->Get(Form("Track%d", j));
    tr[j]->SetBranchAddress("sig_x", &sx[j]);
    tr[j]->SetBranchAddress("sig_y", &sy[j]);
    tr[j]->SetBranchAddress("x", &x[j]);
    tr[j]->SetBranchAddress("y", &y[j]);
  }
  const double centre = 256.0 / 2 * 0.4;
  const double zpos[4] = {533, 468, 217, 154};
  const double dx[4] = {-1.0605, 0.456114, 1.22355, -1.05207};
  const double dy[4] = {-0.938291, 0.241082, 0.463939, -0.533685};

  TH1D hx("hx", "R2 x (all 4 layers);R2;tracks", 200, 0.9, 1.0005);
  TH1D hy("hy", "R2 y (all 4 layers);R2;tracks", 200, 0.9, 1.0005);

  int n = tr[0]->GetEntries();
  int n4x = 0, n4y = 0;
  int nJoint = 0, nJoint3 = 0;
  for (int e = 0; e < n; e++)
  {
    for (int j = 0; j < 4; j++)
      tr[j]->GetEntry(e);
    double r2v[2] = {-1, -1};
    double zs[4] = {zpos[0], zpos[1], zpos[2], zpos[3]};
    for (int axis = 0; axis < 2; axis++)
    {
      std::vector<double> px, pz;
      for (int j = 0; j < 4; j++)
      {
        bool ok = axis == 0 ? sx[j] : sy[j];
        if (!ok)
          continue;
        const double shift = axis == 0 ? dx[j] : dy[j];
        px.push_back((axis == 0 ? x[j] : y[j]) * 0.4 - centre + shift);
        pz.push_back(zs[j]);
      }
      if (px.size() < 3)
        continue;
      if (px.size() == 4)
        axis == 0 ? n4x++ : n4y++;
      double Mz = 0, Mx = 0, Mzz = 0, Mxx = 0, Mzx = 0;
      int nn = px.size();
      for (size_t i = 0; i < px.size(); i++)
      {
        Mz += pz[i] / nn;
        Mx += px[i] / nn;
        Mzz += pz[i] * pz[i] / nn;
        Mxx += px[i] * px[i] / nn;
        Mzx += pz[i] * px[i] / nn;
      }
      double dxt = 0, dyt = 0, dxxt = 0, dyyt = 0, dxyt = 0;
      for (size_t i = 0; i < px.size(); i++)
      {
        double tx = px[i] * TMath::Cos(45 * TMath::DegToRad()) - pz[i] * TMath::Sin(45 * TMath::DegToRad());
        double ty = px[i] * TMath::Sin(45 * TMath::DegToRad()) + pz[i] * TMath::Cos(45 * TMath::DegToRad());
        dxt += tx / nn;
        dyt += ty / nn;
        dxxt += tx * tx / nn;
        dyyt += ty * ty / nn;
        dxyt += tx * ty / nn;
      }
      double r2 = (dxyt - dxt * dyt) * (dxyt - dxt * dyt) / (dxxt - dxt * dxt) / (dyyt - dyt * dyt);
      (axis == 0 ? hx : hy).Fill(r2);
      r2v[axis] = r2;
    }
    if (r2v[0] > 0.9999 && r2v[1] > 0.9999)
      nJoint++;
    if (r2v[0] > 0.9999 && r2v[1] > 0.9999 && sx[0] && sx[1] && sx[2] && sx[3] && sy[0] && sy[1] && sy[2] && sy[3])
      nJoint3++;
  }
  std::cout << "joint R2>0.9999 (x and y, >=3 layers): " << nJoint
            << "  (all 4 layers): " << nJoint3 << std::endl;
  std::cout << "tracks with 4 x-layers: " << n4x << ", 4 y-layers: " << n4y << std::endl;
  for (double c : {0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99999})
    std::cout << "x R2>" << c << ": " << hx.Integral(hx.FindBin(c), hx.GetNbinsX() + 1)
              << "  y R2>" << c << ": " << hy.Integral(hy.FindBin(c), hy.GetNbinsX() + 1) << std::endl;
  f->Close();
}

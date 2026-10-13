#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

bool pass_shape(const std::vector<short> *wave, double ped, double pedrms)
{
  const int n = 30;
  if (!wave || (int)wave->size() != n) return false;
  const double sigma5 = 5.0 * pedrms;
  std::vector<int> peaks;
  for (int i=1;i<n-1;++i) {
    if ((double)wave->at(i)-ped <= sigma5) continue;
    bool is_local_max = true;
    int lb=std::max(0,i-3), rb=std::min(n-1,i+3);
    for (int j=lb;j<=rb;++j) if (j!=i && wave->at(j)>=wave->at(i)) { is_local_max=false; break; }
    if (is_local_max) peaks.push_back(i);
  }
  if (peaks.size()!=1) return false;
  int p=peaks[0];
  return p>0 && p<n-1 && (double)wave->at(p)-ped > sigma5 && wave->at(0)<wave->at(p) && wave->at(n-1)<wave->at(p);
}

void diagnose_normalized_chi2()
{
  TFile f("/home/thoriatong/workplace/BTTracker/arrivalT_run_20260702/Combine/Tracker-step2-ana.root", "READ");
  Long64_t base=0, norm1=0, norm4=0, norm9=0, norm16=0, norm25=0;
  for (int det=0; det<6; ++det) for (const char *xy : {"X","Y"}) {
    TString name=Form("Track%s%d",xy,det);
    TTree *t=(TTree*)f.Get(name);
    if(!t) continue;
    double maxA=0,ped=0,pedrms=0,riseT=0,TOT=0,chi2=0;
    int maxT=0,status=0;
    std::vector<short>* wave=nullptr;
    t->SetBranchAddress("maxA",&maxA);
    t->SetBranchAddress("maxT",&maxT);
    t->SetBranchAddress("ped",&ped);
    t->SetBranchAddress("pedrms",&pedrms);
    t->SetBranchAddress("riseT",&riseT);
    t->SetBranchAddress("TOT",&TOT);
    t->SetBranchAddress("pulseFitStatus",&status);
    t->SetBranchAddress("pulseFitChi2NDF",&chi2);
    t->SetBranchAddress("wave",&wave);
    Long64_t n=t->GetEntries();
    for(Long64_t i=0;i<n;++i){
      t->GetEntry(i);
      if (!(maxA > 7.0*pedrms && maxT>0 && maxT<20 && riseT>0 && riseT<15 && TOT>0 && TOT<20)) continue;
      if (!pass_shape(wave,ped,pedrms)) continue;
      if (status!=0 || chi2<0 || pedrms<=0) continue;
      ++base;
      double nchi2 = chi2/(pedrms*pedrms);
      if(nchi2<=1) ++norm1;
      if(nchi2<=4) ++norm4;
      if(nchi2<=9) ++norm9;
      if(nchi2<=16) ++norm16;
      if(nchi2<=25) ++norm25;
    }
  }
  std::cout << "Valid chi2 after Cut+shape+status: " << base << "\n";
  auto print=[&](const char* label, Long64_t n){ std::cout << std::setw(12) << label << std::setw(12) << n << std::setw(10) << std::fixed << std::setprecision(2) << (base?100.0*n/base:0) << "% of valid, " << (15375942?100.0*n/15375942:0) << "% of old step2\n"; };
  print("nchi2<=1", norm1);
  print("nchi2<=4", norm4);
  print("nchi2<=9", norm9);
  print("nchi2<=16", norm16);
  print("nchi2<=25", norm25);
}
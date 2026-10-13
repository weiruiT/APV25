// root -l -b -q 'debug/macros/plot_peak_timing.C("/absolute/run-directory")'
#ifndef BTTRACKER_PEAK_TIMING_PLOTS_C
#define BTTRACKER_PEAK_TIMING_PLOTS_C
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TNamed.h>
#include <TKey.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace PeakTimingPlots {
const int ActiveTimingModels[] = {
  // 0, // Pause-a disabled, retained for future comparisons.
  // 1, // Pause-e disabled, retained for future comparisons.
  2     // Sigmoid rising-edge t20 only.
};
void Require(bool ok, const std::string &what) {
  if (!ok) throw std::runtime_error(what);
}
void Text(double x,double y,const TString &s,double size=0.036) {
  TLatex l; l.SetNDC(); l.SetTextFont(42); l.SetTextSize(size);
  l.DrawLatex(x,y,s);
}
struct Moments {
  int n=0; double mean=0,m2=0;
  void Add(double x) { ++n; double d=x-mean; mean+=d/n; m2+=d*(x-mean); }
  double RMS() const { return n ? std::sqrt(std::max(0.0,m2/n)) : 0; }
};
struct EventFlags { bool multi=false,written=false; };
int PeakPlots(const TString &dir,TFile &step) {
  TTree *diag=dynamic_cast<TTree*>(step.Get("WaveformPeakDiagnostics"));
  Require(diag,"Missing WaveformPeakDiagnostics");
  // Discover diagnostics from this result; the removed fourth Tracker must
  // not introduce empty groups or a fixed expected PNG count.
  std::vector<std::string> groups;
  const std::string prefix="PeakCount_all_";
  TIter next(step.GetListOfKeys());
  while (auto *key=dynamic_cast<TKey*>(next())) {
    const std::string name=key->GetName();
    if (name.compare(0,prefix.size(),prefix)==0)
      groups.push_back(name.substr(prefix.size()));
  }
  std::sort(groups.begin(),groups.end());
  groups.erase(std::unique(groups.begin(),groups.end()),groups.end());
  Require(!groups.empty(),"Missing peak histograms");
  std::map<std::string,std::map<int,EventFlags>> events;
  std::map<std::string,Long64_t> written;
  int event,det,axis,npeaks; bool kept;
  diag->SetBranchAddress("event",&event);
  diag->SetBranchAddress("detector_id",&det);
  diag->SetBranchAddress("direction",&axis);
  diag->SetBranchAddress("n_significant_peaks",&npeaks);
  diag->SetBranchAddress("written_to_signal_tree",&kept);
  for(Long64_t i=0;i<diag->GetEntries();++i) {
    diag->GetEntry(i);
    Require(det>=1 && det<=6,"Invalid detector id");
    // The on-disk detector IDs for pads remain 5 and 6 for compatibility.
    const bool tracker=det<=4;
    Require(!tracker || axis==0 || axis==1,"Invalid Tracker direction");
    const std::string detector=tracker ? "Tracker"+std::to_string(det)
                                       : "BM"+std::to_string(det-4);
    std::vector<std::string> selected={"All",tracker?"Tracker":"BM",detector};
    if (tracker) selected.push_back(detector+(axis==0?"X":"Y"));
    for(const auto &g:selected) {
      Require(std::binary_search(groups.begin(),groups.end(),g),
              "Diagnostic detector has no histogram group: "+g);
      auto &f=events[g][event]; f.multi|=npeaks>=2; f.written|=kept;
      if(kept) ++written[g];
    }
  }
  diag->ResetBranchAddresses();
  std::ofstream summary((dir+"/peak_summary.tsv").Data());
  summary<<"group\tmapped_waveforms\tafter_basic_cut\tmulti_peak\tmulti_percent_after_cut"
         <<"\twritten_waveforms\tcandidate_events\tevents_with_multi"
         <<"\tevents_with_multi_but_no_selected_waveform\n";
  for(const auto &g:groups) {
    auto *all=dynamic_cast<TH1D*>(step.Get(Form("PeakCount_all_%s",g.c_str())));
    auto *cut=dynamic_cast<TH1D*>(step.Get(Form("PeakCount_afterCut_%s",g.c_str())));
    Require(all && cut,"Missing peak histogram");
    const double multi=cut->Integral(3,cut->GetNbinsX());
    int affected=0,empty=0;
    for(const auto &e:events[g]) {
      affected+=e.second.multi; empty+=e.second.multi && !e.second.written;
    }
    summary<<g.c_str()<<'\t'<<std::fixed<<std::setprecision(0)<<all->GetEntries()
           <<'\t'<<cut->GetEntries()<<'\t'<<multi<<'\t'<<std::setprecision(4)
           <<(cut->GetEntries()?100*multi/cut->GetEntries():0)<<'\t'
           <<written[g]<<'\t'<<events[g].size()<<'\t'<<affected<<'\t'<<empty<<'\n';
    for(int stage=0;stage<2;++stage) {
      TH1D *h=stage?cut:all;
      TCanvas c("peakCanvas","Waveform peak count",1000,720);
      c.SetLeftMargin(.12); c.SetRightMargin(.04); c.SetBottomMargin(.13); c.SetTopMargin(.20);
      c.SetLogy(h->GetEntries()>0);
      int last=std::max(4,h->FindLastBinAbove(0));
      h->GetXaxis()->SetRange(1,std::min(last+1,h->GetNbinsX()));
      h->SetTitle(""); h->SetStats(false); h->SetLineColor(kAzure+2); h->SetLineWidth(2);
      h->SetFillColorAlpha(kAzure+1,.25); h->SetMinimum(.5);
      h->SetMaximum(std::max(10.0,h->GetMaximum()*8));
      h->GetXaxis()->SetNdivisions(510); h->Draw("HIST");
      Text(.12,.95,Form("%s: significant waveform peaks",g.c_str()),.042);
      Text(.12,.90,stage?"After basic Cut; BEFORE single-peak selection":"All mapped 30-sample waveforms",.032);
      Text(.12,.845,Form("Waveforms = %.0f",h->GetEntries()),.032);
      double m=h->Integral(3,h->GetNbinsX());
      Text(.53,.845,Form("N_{peak} #geq 2: %.0f (%.2f%%)",m,
                        h->GetEntries()?100*m/h->GetEntries():0),.032);
      Text(.13,.035,"Peak: above configured run-noise threshold, maximum within #pm3 samples; not a cluster count",.024);
      c.SaveAs(Form("%s/peaks_%s_%s.png",dir.Data(),stage?"afterCut":"all",g.c_str()));
    }
  }
  return 2*static_cast<int>(groups.size());
}

void DeltaPlots(const TString &baseDir,TFile &pad) {
  // Keep corrected results separate from historical cross-estimator plots.
  // const TString dir=baseDir+"/same_estimator"; // Legacy a/e/sigmoid output.
  const TString dir=baseDir+"/sigmoid_only";
  gSystem->mkdir(dir,kTRUE);
  TTree *t=dynamic_cast<TTree*>(pad.Get("TrackArrivalEvents"));
  Require(t,"Missing TrackArrivalEvents");
  const char *tags[]={"A","E","S"},*names[]={"a","e","sigmoid"};
  const char *estimators[]={"mean","min"};
  double bmMean[3][2]={},bmMin[3][2]={};
  int bmN[3][2]={};
  const bool pooledBranches = t->GetBranch("trXYS_mean") != nullptr;
  auto *trackerLeaf=t->GetLeaf(pooledBranches?"trXYS_n":"trS_n");
  Require(trackerLeaf,"Missing Tracker timing count leaf");
  const int trackerCount=trackerLeaf->GetLenStatic()/(pooledBranches?1:2);
  Require(trackerCount>=2 && trackerCount<=4,"Unexpected Tracker timing dimensions");
  std::vector<double> trMean[3],trMin[3],axisMean[3],axisMin[3];
  std::vector<int> trN[3],axisN[3];
  for (int m:ActiveTimingModels) {
    trMean[m].resize(trackerCount); trMin[m].resize(trackerCount);
    trN[m].resize(trackerCount);
    axisMean[m].resize(2*trackerCount); axisMin[m].resize(2*trackerCount);
    axisN[m].resize(2*trackerCount);
  }
  // for(int m=0;m<3;++m) { // Legacy: a/e/sigmoid.
  for(int m: ActiveTimingModels) {
    auto bind=[&](TString name,void *p) {
      Require(t->GetBranch(name),std::string("Missing ")+name.Data());
      Require(t->SetBranchAddress(name,p)>=0,"Cannot bind timing branch");
    };
    bind(Form("bm%s_mean",tags[m]),bmMean[m]); bind(Form("bm%s_min",tags[m]),bmMin[m]);
    bind(Form("bm%s_n",tags[m]),bmN[m]);
    if (pooledBranches) {
      bind(Form("trXY%s_mean",tags[m]),trMean[m].data()); bind(Form("trXY%s_min",tags[m]),trMin[m].data());
      bind(Form("trXY%s_n",tags[m]),trN[m].data());
    } else {
      bind(Form("tr%s_mean",tags[m]),axisMean[m].data()); bind(Form("tr%s_min",tags[m]),axisMin[m].data());
      bind(Form("tr%s_n",tags[m]),axisN[m].data());
    }
  }
  struct Sample { int event; double delta; };
  std::vector<Sample> data[2][3][2];
  int event=-1; t->SetBranchAddress("event",&event);
  std::set<int> seen;
  std::ofstream values((dir+"/delta_values.tsv").Data());
  Require(values.is_open(),"Cannot create delta values table");
  values<<"event\tbm\tbm_model\ttracker_model\tbm_estimator\ttracker_estimator\tdelta_ns\n";
  const auto valid=[](double v){return std::isfinite(v)&&v>-900;};
  for(Long64_t i=0;i<t->GetEntries();++i) {
    t->GetEntry(i); Require(seen.insert(event).second,"Duplicate track event would bias histogram");
    if (!pooledBranches) {
      // Compatibility with older Pad files: pool only the matched strips of
      // each individual detector. Never fall back to all-layer ev* branches.
      // for (int m=0;m<3;++m) for (int j=0;j<trackerCount;++j) {
      for (int m: ActiveTimingModels) for (int j=0;j<trackerCount;++j) {
        int nx=axisN[m][2*j+0],ny=axisN[m][2*j+1];
        Require(nx>=0 && ny>=0,"Negative matched-strip count");
        for (int d=0;d<2;++d)
          Require(axisN[m][2*j+d]==0 || (valid(axisMean[m][2*j+d]) && valid(axisMin[m][2*j+d])),
                  "Inconsistent per-axis time/count");
        trN[m][j]=nx+ny; trMean[m][j]=trMin[m][j]=-999;
        if (nx+ny) {
          trMean[m][j]=((nx?nx*axisMean[m][2*j+0]:0)+(ny?ny*axisMean[m][2*j+1]:0))/(nx+ny);
          trMin[m][j]=std::min(nx?axisMin[m][2*j+0]:1e30,ny?axisMin[m][2*j+1]:1e30);
        }
      }
    }
    // for(int c=0;c<2;++c) for(int m=0;m<3;++m) {
    for(int c=0;c<2;++c) for(int m: ActiveTimingModels) {
      // Same timing model and estimator on BOTH sides. The mean-mean and
      // min-min pair shares one event cohort. Physical Tracker2 is index 1.
      if(bmN[m][c]<=0 || trN[m][1]<=0 || !valid(bmMean[m][c]) || !valid(bmMin[m][c])
         || !valid(trMean[m][1]) || !valid(trMin[m][1])) continue;
      for(int q=0;q<2;++q) {
        double delta=25*((q?bmMin[m][c]:bmMean[m][c])-(q?trMin[m][1]:trMean[m][1]));
        Require(std::isfinite(delta),"Nonfinite delta time");
        data[c][m][q].push_back({event,delta});
        values<<event<<'\t'<<c+1<<'\t'<<names[m]<<'\t'<<names[m]<<'\t'
              <<estimators[q]<<'\t'<<estimators[q]<<'\t'<<std::setprecision(17)<<delta<<'\n';
      }
    }
  }
  t->ResetBranchAddresses();
  TFile output(dir+"/delta_histograms.root","RECREATE");
  Require(!output.IsZombie(),"Cannot create delta output");
  TNamed definition("TimingDefinition",
      "Delta = BM track-centred 3x3 minus physical Tracker2, matched X/Y strips pooled; "
      "25 ns/sample; RMS=sqrt(mean((delta-mean(delta))^2)); no time-difference cut, "
      "no Gaussian fit, no reference subtraction. SAME model on both sides; "
      "sigmoid t20 only; mean-mean and min-min (4 histograms).");
  definition.Write();
  std::ofstream summary((dir+"/delta_summary.tsv").Data());
  Require(summary.is_open(),"Cannot create delta summary table");
  summary<<"bm\tbm_model\ttracker_model\tbm_estimator\ttracker_estimator\tevents"
         <<"\tmean_ns\trms_ns\tunderflow\toverflow\n"<<std::setprecision(12);
  // for(int c=0;c<2;++c) for(int m=0;m<3;++m) {
  for(int c=0;c<2;++c) for(int m: ActiveTimingModels) {
    double low=1e30,high=-1e30;
    for(int q=0;q<2;++q)
      for(const auto &s:data[c][m][q]) {low=std::min(low,s.delta);high=std::max(high,s.delta);}
    if(low>high) {low=-25;high=25;}
    double margin=std::max(5.0,(high-low)*.08); low-=margin;high+=margin;
    for(int q=0;q<2;++q) {
      TString key=Form("delta_BM%d_%s_%s_minus_T2_XY_%s_%s",
                      c+1,names[m],estimators[q],names[m],estimators[q]);
      TH1D h(key,";#Delta t = t_{BM} - t_{Tracker2,XY} [ns];Events",50,low,high);
      h.SetDirectory(nullptr); h.SetStats(false); h.SetLineColor(kAzure+2);
      h.SetLineWidth(2); h.SetFillColorAlpha(kAzure+1,.25);
      Moments stats;
      for(const auto &s:data[c][m][q]) {h.Fill(s.delta);stats.Add(s.delta);}
      Require(h.GetBinContent(0)==0 && h.GetBinContent(h.GetNbinsX()+1)==0,"Clipped delta events");
      Require(std::abs(h.GetRMS()-stats.RMS())<1e-7,"Histogram/unbinned RMS mismatch");
      output.cd(); h.Write();
      summary<<c+1<<'\t'<<names[m]<<'\t'<<names[m]<<'\t'<<estimators[q]<<'\t'
             <<estimators[q]<<'\t'<<stats.n<<'\t';
      if(stats.n) summary<<stats.mean<<'\t'<<stats.RMS();
      else summary<<"NA\tNA";
      summary<<"\t0\t0\n";
      TCanvas canvas("deltaCanvas","BM minus Tracker2",1100,760);
      canvas.SetLeftMargin(.12);canvas.SetRightMargin(.04);
      canvas.SetBottomMargin(.13);canvas.SetTopMargin(.17);
      h.SetMaximum(std::max(1.0,h.GetMaximum())*1.5); h.Draw("HIST");
      Text(.12,.95,Form("BM%d %s (%s) - Tracker2 X+Y %s (%s)",
                        c+1,estimators[q],names[m],estimators[q],names[m]),.036);
      Text(.12,.90,"Track-matched BM 3#times3; matched Tracker2 strips pooled",.029);
      Text(.58,.78,Form("Events = %d",stats.n));
      if(stats.n) {
        Text(.58,.73,Form("Mean = %.3f ns",stats.mean));
        Text(.58,.68,Form("RMS = %.3f ns",stats.RMS()),.041);
      } else Text(.58,.73,"RMS = N/A (no valid events)");
      if(stats.n>0 && stats.n<30) Text(.58,.62,"Small sample",.030);
      Text(.12,.035,"RMS: full valid sample; no #Delta t cut or Gaussian fit; reference contribution included",.024);
      canvas.SaveAs(dir+"/"+key+".png");
    }
  }
}
}

void plot_peak_timing(TString directory) {
  try {
    gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0);
    TFile step(directory+"/Tracker-step2-ana.root","READ");
    TFile pad(directory+"/Pad-step3-rec.root","READ");
    PeakTimingPlots::Require(!step.IsZombie()&&!pad.IsZombie(),"Missing analysis output");
    const int peakCount=PeakTimingPlots::PeakPlots(directory,step);
    PeakTimingPlots::DeltaPlots(directory,pad);
    std::cout<<"Wrote "<<peakCount<<" peak-count plots; 4 sigmoid mean-mean/min-min delta plots "
             <<"and ROOT/TSV outputs in "<<directory<<"/sigmoid_only/\n";
  } catch(const std::exception &e) {
    std::cerr<<"ERROR: "<<e.what()<<'\n'; gSystem->Exit(1);
  }
}
#endif

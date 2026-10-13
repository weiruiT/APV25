#include "time_resolution.C"

void test_timing_entry(TString workspace,TString verified) {
  try {
    gROOT->SetBatch(kTRUE);

    // Build an old-schema compatibility input from this same sigmoid-only run:
    // keep the per-axis trS* statistics, but remove the newer pooled trXYS*
    // branches. This avoids comparing different event selections/data sets.
    const TString compatibilityInput=verified+"/Pad-step3-rec.per_axis.root";
    {
      TFile source(verified+"/Pad-step3-rec.root","READ");
      auto *input=dynamic_cast<TTree*>(source.Get("TrackArrivalEvents"));
      PeakTimingPlots::Require(input,"Missing TrackArrivalEvents test input");
      input->SetBranchStatus("*",1);
      input->SetBranchStatus("trXY*",0);
      TFile compatibility(compatibilityInput,"RECREATE");
      PeakTimingPlots::Require(!compatibility.IsZombie(),"Cannot create compatibility input");
      auto *copy=input->CloneTree(-1,"fast");
      PeakTimingPlots::Require(copy,"Cannot clone per-axis compatibility tree");
      copy->Write();
      compatibility.Close();
      input->SetBranchStatus("*",1);
    }

    // Exercise the public entry and its old-file X/Y pooling fallback.
    delta_time_resolution(compatibilityInput);
    TFile expected(verified+"/sigmoid_only/delta_histograms.root","READ");
    TFile actual(workspace+"/debug/reports/bm_tracker2_timing/sigmoid_only/delta_histograms.root","READ");
    PeakTimingPlots::Require(!expected.IsZombie() && !actual.IsZombie(),"Missing compatibility outputs");

    int expectedCount=0,actualCount=0,compared=0;
    TIter expectedKeys(expected.GetListOfKeys());
    while (TObject *key=expectedKeys()) {
      auto *a=dynamic_cast<TH1D*>(expected.Get(key->GetName()));
      if (!a) continue;
      ++expectedCount;
      auto *b=dynamic_cast<TH1D*>(actual.Get(key->GetName()));
      PeakTimingPlots::Require(b && a->GetEntries()==b->GetEntries() &&
          std::abs(a->GetMean()-b->GetMean())<1e-7 &&
          std::abs(a->GetRMS()-b->GetRMS())<1e-7,
          "Pooled/per-axis sigmoid reference mismatch");
      ++compared;
    }
    TIter actualKeys(actual.GetListOfKeys());
    while (TObject *key=actualKeys())
      if (dynamic_cast<TH1D*>(actual.Get(key->GetName()))) ++actualCount;
    PeakTimingPlots::Require(expectedCount==4 && actualCount==4 && compared==4,
                             "Expected exactly 4 sigmoid delta-time histograms");

    for(int bm=1;bm<=2;++bm)
      for(const char *estimator: {"mean","min"}) {
        TString name=Form("delta_BM%d_sigmoid_%s_minus_T2_XY_sigmoid_%s",
                          bm,estimator,estimator);
        PeakTimingPlots::Require(expected.Get(name) && actual.Get(name),
                                 "Missing sigmoid same-estimator histogram");
      }
    std::cout<<"PASS public delta entry + old-file fallback: exactly 4 sigmoid "
             <<"mean-mean/min-min histograms; N/mean/RMS agree\n";
  } catch(const std::exception &error) {
    std::cerr<<"ERROR: "<<error.what()<<'\n'; gSystem->Exit(1);
  }
}

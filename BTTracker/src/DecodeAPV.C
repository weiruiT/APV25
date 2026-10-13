#include "DecodeAPV.h"
#include <algorithm>
#include <set>
#include <cmath> // 鐢ㄤ簬 std::clamp
#include "TStyle.h"  // 鐢ㄤ簬 gStyle
#include "TNamed.h"

DecodeAPV::DecodeAPV(TString fanaName, TString fdecName, int f_force)
{
  anaName = fanaName;
  decName = fdecName;
  force = f_force;
}

DecodeAPV::~DecodeAPV()
{
}


bool DecodeAPV::Init()
{
  fanaFile = new TFile(anaName, "read");
  if (!fanaFile->IsOpen())
  {
    cout << " ERROR! " << anaName << " cant open!" << endl;
    return false;
  }
  cout << " --> Decoding strips from " << anaName << endl;
  for (int i = 0; i < 2 * MMdetN; i++)
  {
    fanaTree[i] = (TTree *)fanaFile->Get(Form("Track%s%d", XYname(i % 2).Data(), i / 2));
    if (fanaTree[i] == 0)
    {
      cout << " ERROR! Track" << XYname(i % 2) << i / 2 << " Tree not exist in " << anaName << endl;
      return false;
    }
    fanaTree[i]->SetBranchAddress("event", &event, &b_event);
    fanaTree[i]->SetBranchAddress("maxA", &maxA, &b_maxA);
    fanaTree[i]->SetBranchAddress("maxT", &maxT, &b_maxT);
    fanaTree[i]->SetBranchAddress("channel", &channel, &b_channel);
    if (fanaTree[i]->GetBranch("arrivalT"))
      fanaTree[i]->SetBranchAddress("arrivalT", &arrivalT);
    if (fanaTree[i]->GetBranch("arrivalT_e"))
      fanaTree[i]->SetBranchAddress("arrivalT_e", &arrivalT_e);
    else if (fanaTree[i]->GetBranch("arrivetime_e"))
      fanaTree[i]->SetBranchAddress("arrivetime_e", &arrivalT_e);
    if (fanaTree[i]->GetBranch("arrivalT_sigmoid"))
      fanaTree[i]->SetBranchAddress("arrivalT_sigmoid", &arrivalT_sigmoid);
    else if (fanaTree[i]->GetBranch("arrivetime_sigmoid"))
      fanaTree[i]->SetBranchAddress("arrivetime_sigmoid", &arrivalT_sigmoid);
    fanaTree[i]->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
  }
  if (force != 1)
  {
    FileStat_t fStat;
    gSystem->GetPathInfo(decName, fStat);
    if (fStat.fSize != 0)
    {
      cout << " ERROR! " << decName << " is exist!" << endl;
      return false;
    }
  }
  fdecFile = new TFile(decName, "recreate");
  if (!fdecFile->IsOpen())
  {
    cout << "ERROR! " << decName << " cant open!" << endl;
    return false;
  }
  for (int i = 0; i < MMdetN; i++)
  {
    fdecTree[i] = new TTree(Form("Track%d", i), Form("Track%d", i));
    // Compact one-based event sequence shared by Track0/1/2. ROOT's
    // implicit Entry$ remains zero-based; event keeps the DAQ trigger ID.
    fdecTree[i]->Branch("entry", &entry, "entry/L");
    fdecTree[i]->Branch("event", &event);
    fdecTree[i]->Branch("sig", &sig);
    fdecTree[i]->Branch("sig_x", &sig_x);
    fdecTree[i]->Branch("sig_y", &sig_y);
    fdecTree[i]->Branch("x", &x);
    fdecTree[i]->Branch("y", &y);
    fdecTree[i]->Branch("z", &z);
    fdecTree[i]->Branch("hit_amp_x", &hit_amp_x);
    fdecTree[i]->Branch("hit_amp_y", &hit_amp_y);
    fdecTree[i]->Branch("hit_strip_num_x", &hit_strip_num_x);
    fdecTree[i]->Branch("hit_strip_num_y", &hit_strip_num_y);
    fdecTree[i]->Branch("x_nhits", &x_nhits);
    fdecTree[i]->Branch("y_nhits", &y_nhits);
    fdecTree[i]->Branch("x_other", &x_other);
    fdecTree[i]->Branch("y_other", &y_other);
    fdecTree[i]->Branch("timestamp", &timestamp, "timestamp/L");
    fdecTree[i]->Branch("xMaxT", &xMaxT);
    fdecTree[i]->Branch("yMaxT", &yMaxT);
    fdecTree[i]->Branch("xArrivalT", &xArrivalT);
    fdecTree[i]->Branch("yArrivalT", &yArrivalT);
    fdecTree[i]->Branch("xArrivalT_min", &xArrivalT_min);
    fdecTree[i]->Branch("yArrivalT_min", &yArrivalT_min);
    fdecTree[i]->Branch("xArrivalT_n", &xArrivalT_n);
    fdecTree[i]->Branch("yArrivalT_n", &yArrivalT_n);
    fdecTree[i]->Branch("xArrivalTe", &xArrivalTe);
    fdecTree[i]->Branch("yArrivalTe", &yArrivalTe);
    fdecTree[i]->Branch("xArrivalTe_min", &xArrivalTe_min);
    fdecTree[i]->Branch("yArrivalTe_min", &yArrivalTe_min);
    fdecTree[i]->Branch("xArrivalTSigmoid", &xArrivalTSigmoid);
    fdecTree[i]->Branch("yArrivalTSigmoid", &yArrivalTSigmoid);
    fdecTree[i]->Branch("xArrivalTSigmoid_min", &xArrivalTSigmoid_min);
    fdecTree[i]->Branch("yArrivalTSigmoid_min", &yArrivalTSigmoid_min);
    fdecTree[i]->Branch("cluster_count_x", &cluster_count_x);
    fdecTree[i]->Branch("cluster_count_y", &cluster_count_y);
    fdecTree[i]->Branch("cluster_peak_count_x", &cluster_peak_count_x);
    fdecTree[i]->Branch("cluster_peak_count_y", &cluster_peak_count_y);
  }
  return true;
}

void DecodeAPV::InitData()
{
  sig = false;
  sig_x = false;
  sig_y = false;
  x = 0;
  y = 0;
  z = 0;
  hit_amp_x = 0;
  hit_amp_y = 0;
  hit_strip_num_x = 0;
  hit_strip_num_y = 0;
  // timestamp = 0;
  x_nhits = 0;
  y_nhits = 0;
  xMaxT = -1;
  yMaxT = -1;
  xArrivalT = -999;
  yArrivalT = -999;
  xArrivalT_min = -999;
  yArrivalT_min = -999;
  xArrivalT_n = 0;
  yArrivalT_n = 0;
  xArrivalTe = -999;
  yArrivalTe = -999;
  xArrivalTe_min = -999;
  yArrivalTe_min = -999;
  xArrivalTSigmoid = -999;
  yArrivalTSigmoid = -999;
  xArrivalTSigmoid_min = -999;
  yArrivalTSigmoid_min = -999;
  x_other.clear();
  y_other.clear();
  cluster_count_x = 0;
  cluster_count_y = 0;
  cluster_peak_count_x = 0;
  cluster_peak_count_y = 0;
}

void DecodeAPV::Loop()
{
  ofstream op = LogFile(decName);
  TCanvas *cTmp = new TCanvas("tmpraw", "tmpraw");
  TH1F *hx[MMdetN];
  TH1F *hy[MMdetN];
  TH2F *hxy[MMdetN];
  TH1F *hampx[MMdetN];
  TH1F *hampy[MMdetN];
  TH1F *hampxy[MMdetN];
  TH1I *hstripx[MMdetN];
  TH1I *hstripy[MMdetN];
  TH2F *damp[MMdetN];
  for (int i = 0; i < MMdetN; i++)
  {
    hx[i] = new TH1F(Form("hitx%d", i), Form("hitx%d", i), 100, -100, 100);
    hx[i]->GetXaxis()->SetTitle("X [mm]");
    hx[i]->GetYaxis()->SetTitle("Counts");
    hy[i] = new TH1F(Form("hity%d", i), Form("hity%d", i), 100, -100, 100);
    hy[i]->GetXaxis()->SetTitle("Y [mm]");
    hy[i]->GetYaxis()->SetTitle("Counts");
    hxy[i] = new TH2F(Form("hitxy%d", i), Form("hitxy%d", i), 100, -100, 100, 100, -100, 100);
    hxy[i]->GetXaxis()->SetTitle("X [mm]");
    hxy[i]->GetYaxis()->SetTitle("Y [mm]");
    hampx[i] = new TH1F(Form("ampx%d", i), Form("ampx%d", i), 150, 0, 15000);
    hampx[i]->GetXaxis()->SetTitle("amplitude X [ADC]");
    hampx[i]->GetYaxis()->SetTitle("Counts");
    hampy[i] = new TH1F(Form("ampy%d", i), Form("ampy%d", i), 150, 0, 15000);
    hampy[i]->GetXaxis()->SetTitle("amplitude Y [ADC]");
    hampy[i]->GetYaxis()->SetTitle("Counts");
    hampxy[i] = new TH1F(Form("ampxy%d", i), Form("ampxy%d", i), 150, 0, 30000);
    hampxy[i]->GetXaxis()->SetTitle("amplitude X+Y [ADC]");
    hampxy[i]->GetYaxis()->SetTitle("Counts");
    hstripx[i] = new TH1I(Form("Nstripx%d", i), Form("Nstripx%d", i), 255, 0, 255);
    hstripx[i]->GetXaxis()->SetTitle("Num of fired X-strips");
    hstripx[i]->GetYaxis()->SetTitle("Counts");
    hstripy[i] = new TH1I(Form("Nstripy%d", i), Form("Nstripy%d", i), 255, 0, 255);
    hstripy[i]->GetXaxis()->SetTitle("Num of fired Y-strips");
    hstripy[i]->GetYaxis()->SetTitle("Counts");
    damp[i] = new TH2F(Form("N%d", i), Form("N%d", i), 255, 0, 255, 255, 0, 255);
    damp[i]->GetXaxis()->SetTitle("x");
    damp[i]->GetYaxis()->SetTitle("y");
  }
  vector<int> eventvec[2 * MMdetN];
  set<int> trackerEvents;
  for (int i = 0; i < 2 * MMdetN; i++)
  {
    int nentries = fanaTree[i]->GetEntriesFast();
    for (int j = 0; j < nentries; j++)
    {
      fanaTree[i]->GetEntry(j);
      eventvec[i].push_back(event);
      trackerEvents.insert(event);
    }
    cout << " --> Tracker" << XYname(i % 2) << i / 2
         << ", waveform entries: " << nentries;
    op << " --> Tracker" << XYname(i % 2) << i / 2
       << ", waveform entries: " << nentries;
    if (!eventvec[i].empty())
    {
      cout << ", trigger range: " << eventvec[i].front()
           << " to " << eventvec[i].back();
      op << ", trigger range: " << eventvec[i].front()
         << " to " << eventvec[i].back();
    }
    else
    {
      cout << ", no accepted events";
      op << ", no accepted events";
    }
    cout << endl;
    op << endl;
  }
  cout << " --> Compacting " << trackerEvents.size()
       << " unique Tracker events" << endl;
  op << " --> Compacting " << trackerEvents.size()
     << " unique Tracker events" << endl;

  entry = 0;
  Long64_t processedEvents = 0;
  for (const int sourceEvent : trackerEvents)
  {
    ++entry;
    DrawProcessbar(processedEvents++, trackerEvents.size());

    int firstPos[2 * MMdetN];
    int lastPos[2 * MMdetN];
    bool matched[2 * MMdetN];
    for (int j = 0; j < 2 * MMdetN; j++)
    {
      const auto lower = lower_bound(eventvec[j].begin(), eventvec[j].end(),
                                     sourceEvent);
      const auto upper = upper_bound(lower, eventvec[j].end(), sourceEvent);
      firstPos[j] = static_cast<int>(distance(eventvec[j].begin(), lower));
      lastPos[j] = static_cast<int>(distance(eventvec[j].begin(), upper));
      matched[j] = lower != upper;
    }
    long long ftime = 0;
    for (int j = 0; j < 2 * MMdetN; j++)
    {
      if (j % 2 == 0)
      {
        InitData();
      }
      if (matched[j])
      {
        vector<int> aget_ch;
        vector<double> aget_amp;
        vector<int> aget_t;
        vector<double> aget_arrival;
        vector<double> aget_arrival_e;
        vector<double> aget_arrival_sig;
        // vector<vector<short>> aget_waves;w
        bool firstcount = true;
        bool writeflag = false;
        bool tracker1maxflag = 0;
        for (int k = firstPos[j]; k < lastPos[j]; k++)
        {
          arrivalT = -999;
          arrivalT_e = -999;
          arrivalT_sigmoid = -999;
          fanaTree[j]->GetEntry(k);
          if (event == sourceEvent && firstcount)
          {
            ftime = timestamp;
            firstcount = false;
          }
          // if (event < 2) cout << "k" << k << "j" << j << "event" << event << "timestamp" << timestamp << endl;
          // if (j % 2 == 1) ftime = timestamp;
          // step2 is the single source of waveform selection. step3 only
          // groups the accepted strips into clusters and events.
          {
            aget_ch.push_back(channel);
            aget_amp.push_back(maxA);
            aget_t.push_back(maxT);
            aget_arrival.push_back(arrivalT);
            aget_arrival_e.push_back(arrivalT_e);
            aget_arrival_sig.push_back(arrivalT_sigmoid);
            // if (maxA > 1200 && j == 2)
            //   tracker1maxflag = 1;
          }
        }
        if (aget_ch.size() != 0)
        {
          // if (CalHit(aget_ch, aget_amp, aget_t, j))
          if (CalHitWithClusters(aget_ch, aget_amp, aget_t, aget_arrival, aget_arrival_e, aget_arrival_sig, j))
          {
            if (j % 2 == 0)
            {
              if (x > 0)
              {
                hx[j / 2]->Fill((x - 256 / 2) * 0.4); // total strips 384
                hampx[j / 2]->Fill(hit_amp_x);
              }
              hstripx[j / 2]->Fill(hit_strip_num_x);
            }
            if (j % 2 == 1)
            {
              if (y > 0)
              {
                hy[j / 2]->Fill((y - 256 / 2) * 0.4); // total strips 384
                hampy[j / 2]->Fill(hit_amp_y);
                hampxy[j / 2]->Fill(hit_amp_x + hit_amp_y);
              }
              hstripy[j / 2]->Fill(hit_strip_num_y);
            }
            sig = sig_x && sig_y;
            if (sig && x > 0 && y > 0)
            {
              hxy[j / 2]->Fill((x - 256 / 2) * 0.4, (y - 256 / 2) * 0.4);
              double x_mm = (x - 128) * 0.4;
              double y_mm = (y - 128) * 0.4;
              int region_x = static_cast<int>((x_mm + 40) / 8.0);
              int region_y = static_cast<int>((y_mm + 40) / 8.0);
            }
            z = Zpos[j / 2];
          }
          else
          {
            // cout<<" --> Decode failed Tracker"<<XYname(j%2)<<j/2<<", event ID: "<<event<< endl;
            // op<<" --> Decode failed Tracker"<<XYname(j%2)<<j/2<<", event ID: "<<event<< endl;
          }
        }
        else
        {
          event = sourceEvent;
          // cout<<" --> No record hit at Tracker"<<XYname(j%2)<<j/2<<", event ID: "<<event<< endl;
          // op<<" --> No record hit at Tracker"<<XYname(j%2)<<j/2<<", event ID: "<<event<< endl;
        }
      }
      else
      {
        event = sourceEvent;
        // cout<<" --> Lost record hit at Tracker"<<XYname(j%2)<<j/2<<", event ID: "<<event<< endl;
        op << " --> Lost record hit at Tracker" << XYname(j % 2) << j / 2 << ", event ID: " << event << endl;
      }
      if (j % 2 == 1)
      {
        event = sourceEvent;
        timestamp = ftime;
        // if (event < 10) cout << timestamp << "\tevent" << event << endl;
        fdecTree[j / 2]->Fill();
      }
    }
  }
  for (int i = 0; i < MMdetN; i++)
  {
    fdecFile->WriteTObject(fdecTree[i]);
  }
  TNamed eventIndexPolicy(
      "Step3EventIndexPolicy",
      "entry=1..N is the sorted union of event IDs present in the six Step2 "
      "Tracker direction trees; entry is shared by Track0/1/2; event retains "
      "the original extended DAQ trigger ID; ROOT Entry$ remains 0-based");
  fdecFile->WriteTObject(&eventIndexPolicy);
  TCanvas *cMMhit = new TCanvas("MMhit", "MMhit");
  cMMhit->Divide(MMdetN, 3);
  for (int i = 0; i < MMdetN; i++)
  {
    cMMhit->cd(i + 1);
    hx[i]->Draw();
    cMMhit->cd(i + MMdetN + 1);
    hy[i]->Draw();
    cMMhit->cd(i + 2 * MMdetN + 1);
    hxy[i]->Draw("COLZ");
  }
  cMMhit->Write();
  TCanvas *cMMamp = new TCanvas("MMamp", "MMamp");
  cMMamp->Divide(MMdetN, 5);
  TLatex *la;
  for (int i = 0; i < MMdetN; i++)
  {
    cMMamp->cd(i + 1);
    hampx[i]->Fit("landau", "Q");
    hampx[i]->Draw();
    double peakx = hampx[i]->GetFunction("landau")->GetParameter(1);
    la = CreatLatex(Form("MPV=%.f", peakx), 0.65, 0.5);
    la->Draw("same");
    cMMamp->cd(i + MMdetN + 1);
    hampy[i]->Fit("landau", "Q");
    hampy[i]->Draw();
    double peaky = hampy[i]->GetFunction("landau")->GetParameter(1);
    la = CreatLatex(Form("MPV=%.f", peaky), 0.65, 0.5);
    la->Draw("same");
    cMMamp->cd(i + 2 * MMdetN + 1);
    hampxy[i]->Fit("landau", "Q");
    hampxy[i]->Draw();
    double peakxy = hampxy[i]->GetFunction("landau")->GetParameter(1);
    la = CreatLatex(Form("MPV=%.f", peakxy), 0.65, 0.5);
    la->Draw("same");
    cMMamp->cd(i + 3 * MMdetN + 1);
    hstripx[i]->Draw();
    cMMamp->cd(i + 4 * MMdetN + 1);
    hstripy[i]->Draw();
  }
  cMMamp->Write();
  TCanvas *cMMamp2 = new TCanvas("MMamp2", "MMamp2");
  cMMamp2->Divide(MMdetN, 3);
  for (int i = 0; i < MMdetN; i++)
  {
    cMMamp2->cd(i + 1);
    hampx[i]->Fit("landau", "Q");
    hampx[i]->Draw();
    double peakx = hampx[i]->GetFunction("landau")->GetParameter(1);
    la = CreatLatex(Form("MPV=%.f", peakx), 0.65, 0.5);
    la->Draw("same");
    cMMamp2->cd(i + MMdetN + 1);
    hampy[i]->Fit("landau", "Q");
    hampy[i]->Draw();
    double peaky = hampy[i]->GetFunction("landau")->GetParameter(1);
    la = CreatLatex(Form("MPV=%.f", peaky), 0.65, 0.5);
    la->Draw("same");
    cMMamp2->cd(i + 2 * MMdetN + 1);
    hampxy[i]->Fit("landau", "Q");
    hampxy[i]->Draw();
    double peakxy = hampxy[i]->GetFunction("landau")->GetParameter(1);
    la = CreatLatex(Form("MPV=%.f", peakxy), 0.65, 0.5);
    la->Draw("same");
  }
  cMMamp2->Write();

  fanaFile->Close();
  fdecFile->Close();
  cout << " --> Decoding results has been stored to " << decName << endl;
  op << " --> Decoding results has been stored to " << decName << endl;
}

bool DecodeAPV::CalHit(vector<int> ch, vector<double> amp, vector<int> t, int XY)
{                             /// waves are easy to maximum
  double overthreshold = 100; // maxA should be lower than that
  if ((ch.size() == 0) || (ch.size() > 255))
    return false;
  int nn = 0;
  double temp = 0;
  int locmax = 0;
  int channelloc = 0;
  double Q1, Q2;
  double m, n;
  double Qall = 0;
  double Q0 = overthreshold;
  double A, B, C, D, E, F, G, H;
  // 鎵惧埌鏈€澶у箙搴︾殑浣嶇疆
  size_t max_idx = 0;
  double max_amp = amp[0];
  for (size_t i = 1; i < amp.size(); i++)
  {
    if (amp[i] > max_amp)
    {
      max_amp = amp[i];
      max_idx = i;
    }
  }
  // 纭畾绐楀彛杈圭晫 [max_ch - 5, max_ch + 5]
  int max_ch = ch[max_idx];
  // int window_start = max_ch - 8;
  // int window_end = max_ch + 8;
  int window_start = max_ch - 3;
  int window_end = max_ch + 3;
  if (0) // TMath::MaxElement(amp.size(), &amp.at(0)) > 1200
  {
    if (ch.size() == 1)
    {
      temp = ch.at(0);
    }
    else if (ch.size() == 2)
    {
      temp = (ch.at(0) + ch.at(1)) / 2;
    }
    else
    {
      int channelbegin = ch.at(0);
      int channelend = ch.at(ch.size() - 1);
      // locmax = TMath::LocMax(amp.size(), &amp.at(0));
      // for (int i = locmax; i < 256; i++)
      // {
      //   if (amp.at(i) < overthreshold)
      //   {
      //     channelend = i;
      //     break;
      //   }
      // }
      // for (int i = locmax; i > 0; i--)
      // {
      //   if (amp.at(i) < overthreshold)
      //   {
      //     channelbegin = i;
      //     break;
      //   }
      // }
      if (ch.size() > 1)
      {
        temp = abs(channelend - channelbegin) / 2.0;
        // A = (overthreshold - amp.at(max(0, channelbegin - 2))) / 2.0;
        // B = (overthreshold - amp.at(min(channelend + 2, 30))) / 2.0;
        // temp = (channelbegin * A - channelend * B) / (A - B);
      }
    }
  }
  else
  {
    for (int i = 0; i < ch.size(); i++)
    {
      if (ch.at(i) > window_start && ch.at(i) < window_end)
      {
        temp += amp.at(i) * ch.at(i);
        Qall += amp.at(i);
        nn++;
      }
    }
    temp /= Qall; // average
  }
  // cout<<"data : "<<endl;
  // for(int i=0;i<ch.size();i++){
  //   cout<<"channel "<<ch.at(i)<<" and amp "<<amp.at(i)<<endl;
  // }
  // cout<<"event end"<<endl;

  if (XY % 2 == 0)
  {
    sig_x = true;
    x = temp;
    temp = 0;
    // x-=128/2;
    // hit_amp_x = TMath::MaxElement(amp.size(), &amp.at(0));
    hit_amp_x = Qall;
    hit_strip_num_x = nn;
    for (int s = 1; s < ch.size() && x_nhits < 5; s++, x_nhits++)
      x_other.push_back(ch.at(s));
    Qall = 0;
    xMaxT = t[max_idx];
  }
  if (XY % 2 == 1)
  {
    sig_y = true;
    y = temp;
    temp = 0;
    // hit_amp_y = TMath::MaxElement(amp.size(), &amp.at(0));
    hit_amp_y = Qall;
    hit_strip_num_y = nn;
    for (int s = 1; s < ch.size() && y_nhits < 5; s++, y_nhits++)
      y_other.push_back(ch.at(s));
    Qall = 0;
    yMaxT = t[max_idx];
  }

  return true;
}

bool DecodeAPV::CalHitWithClusters(vector<int> ch, vector<double> amp, vector<int> t, vector<double> arrival, vector<double> arrivalE, vector<double> arrivalS, int XY)
{
  if (ch.empty() || ch.size() > 255)
    return false;
  if (amp.size() != ch.size() || t.size() != ch.size())
    return false;
  if (arrival.size() != ch.size())
    arrival.assign(ch.size(), -999.0);
  if (arrivalE.size() != ch.size())
    arrivalE.assign(ch.size(), -999.0);
  if (arrivalS.size() != ch.size())
    arrivalS.assign(ch.size(), -999.0);

  vector<vector<int>> clusters_ch;
  vector<vector<double>> clusters_amp;
  vector<vector<int>> clusters_t;
  vector<vector<double>> clusters_arrival;
  vector<vector<double>> clusters_arrivalE;
  vector<vector<double>> clusters_arrivalS;

  vector<size_t> indices(ch.size());
  for (size_t i = 0; i < indices.size(); ++i)
    indices[i] = i;
  sort(indices.begin(), indices.end(),
       [&ch](size_t i1, size_t i2) { return ch[i1] < ch[i2]; });

  vector<int> current_ch;
  vector<double> current_amp;
  vector<int> current_t;
  vector<double> current_arrival;
  vector<double> current_arrivalE;
  vector<double> current_arrivalS;

  for (size_t order = 0; order < indices.size(); ++order)
  {
    const size_t idx = indices[order];
    if (!current_ch.empty() && ch[idx] - current_ch.back() > 2)
    {
      clusters_ch.push_back(current_ch);
      clusters_amp.push_back(current_amp);
      clusters_t.push_back(current_t);
      clusters_arrival.push_back(current_arrival);
      clusters_arrivalE.push_back(current_arrivalE);
      clusters_arrivalS.push_back(current_arrivalS);
      current_ch.clear();
      current_amp.clear();
      current_t.clear();
      current_arrival.clear();
      current_arrivalE.clear();
      current_arrivalS.clear();
    }
    current_ch.push_back(ch[idx]);
    current_amp.push_back(amp[idx]);
    current_t.push_back(t[idx]);
    current_arrival.push_back(arrival[idx]);
    current_arrivalE.push_back(arrivalE[idx]);
    current_arrivalS.push_back(arrivalS[idx]);
  }

  if (!current_ch.empty())
  {
    clusters_ch.push_back(current_ch);
    clusters_amp.push_back(current_amp);
    clusters_t.push_back(current_t);
    clusters_arrival.push_back(current_arrival);
    clusters_arrivalE.push_back(current_arrivalE);
    clusters_arrivalS.push_back(current_arrivalS);
  }
  if (clusters_ch.empty())
    return false;

  size_t main_cluster_idx = 0;
  double max_total_amp = -1.0;
  for (size_t i = 0; i < clusters_amp.size(); ++i)
  {
    double total_amp = 0.0;
    for (double value : clusters_amp[i])
      total_amp += value;
    if (total_amp > max_total_amp)
    {
      max_total_amp = total_amp;
      main_cluster_idx = i;
    }
  }

  const vector<int> &main_ch = clusters_ch[main_cluster_idx];
  const vector<double> &main_amp = clusters_amp[main_cluster_idx];
  const vector<int> &main_t = clusters_t[main_cluster_idx];
  const vector<double> &main_arrival = clusters_arrival[main_cluster_idx];
  const vector<double> &main_arrivalE = clusters_arrivalE[main_cluster_idx];
  const vector<double> &main_arrivalS = clusters_arrivalS[main_cluster_idx];

  int peak_count = 0;
  for (size_t i = 0; i < main_ch.size(); ++i)
  {
    bool is_peak = true;
    for (int j = static_cast<int>(i) - 1; j >= 0; --j)
    {
      const int delta = std::abs(main_ch[j] - main_ch[i]);
      if (delta == 1 && main_amp[i] <= main_amp[j])
      {
        is_peak = false;
        break;
      }
      if (delta > 1)
        break;
    }
    if (!is_peak)
      continue;

    for (size_t j = i + 1; j < main_ch.size(); ++j)
    {
      const int delta = std::abs(main_ch[j] - main_ch[i]);
      if (delta == 1 && main_amp[i] <= main_amp[j])
      {
        is_peak = false;
        break;
      }
      if (delta > 2)
        break;
    }
    if (is_peak)
      ++peak_count;
  }

  size_t max_idx = 0;
  for (size_t i = 1; i < main_amp.size(); ++i)
  {
    if (main_amp[i] > main_amp[max_idx])
      max_idx = i;
  }

  const int max_ch = main_ch[max_idx];
  const int window_start = max_ch - 10;
  const int window_end = max_ch + 10;
  double position_numerator = 0.0;
  double charge = 0.0;
  int position_strip_count = 0;
  for (size_t i = 0; i < main_ch.size(); ++i)
  {
    if (main_ch[i] > window_start && main_ch[i] < window_end)
    {
      position_numerator += main_amp[i] * main_ch[i];
      charge += main_amp[i];
      ++position_strip_count;
    }
  }
  if (charge <= 0.0)
    return false;

  const double position = position_numerator / charge;

  // Timing statistics of the event. With useAllStripsForTime=true the mean
  // and n are computed over ALL fired strips of this event for the current
  // direction; with false the legacy main-cluster-only behaviour is kept
  // (statistics over main_arrival/main_arrivalE/main_arrivalS).
  const bool useAllStripsForTime = true;

  double sum_a = 0.0, sum_e = 0.0, sum_s = 0.0;
  double min_a = 1.0e30, min_e = 1.0e30, min_s = 1.0e30;
  int n_a = 0, n_e = 0, n_s = 0;

  if (useAllStripsForTime)
  {
    for (size_t i = 0; i < ch.size(); ++i)
    {
      const double ta = arrival[i];
      const double te = arrivalE[i];
      const double ts = arrivalS[i];
      if (std::isfinite(ta) && ta > -900.0)
      {
        sum_a += ta;
        min_a = std::min(min_a, ta);
        ++n_a;
      }
      if (std::isfinite(te) && te > -900.0)
      {
        sum_e += te;
        min_e = std::min(min_e, te);
        ++n_e;
      }
      if (std::isfinite(ts) && ts > -900.0)
      {
        sum_s += ts;
        min_s = std::min(min_s, ts);
        ++n_s;
      }
    }
  }
  else
  {
    for (size_t i = 0; i < main_ch.size(); ++i)
    {
      const double ta = main_arrival[i];
      const double te = main_arrivalE[i];
      const double ts = main_arrivalS[i];
      if (std::isfinite(ta) && ta > -900.0)
      {
        sum_a += ta;
        min_a = std::min(min_a, ta);
        ++n_a;
      }
      if (std::isfinite(te) && te > -900.0)
      {
        sum_e += te;
        min_e = std::min(min_e, te);
        ++n_e;
      }
      if (std::isfinite(ts) && ts > -900.0)
      {
        sum_s += ts;
        min_s = std::min(min_s, ts);
        ++n_s;
      }
    }
  }
  const double mean_a = n_a > 0 ? sum_a / n_a : -999.0;
  const double mean_e = n_e > 0 ? sum_e / n_e : -999.0;
  const double mean_s = n_s > 0 ? sum_s / n_s : -999.0;
  const double minOut_a = n_a > 0 ? min_a : -999.0;
  const double minOut_e = n_e > 0 ? min_e : -999.0;
  const double minOut_s = n_s > 0 ? min_s : -999.0;

  if (XY % 2 == 0)
  {
    sig_x = true;
    x = position;
    hit_amp_x = charge;
    hit_strip_num_x = position_strip_count;
    x_nhits = static_cast<int>(main_ch.size());
    for (size_t s = 1; s < main_ch.size() && x_other.size() < 5; ++s)
      x_other.push_back(main_ch[s]);
    xMaxT = main_t[max_idx];
    xArrivalT = mean_a;
    xArrivalT_min = minOut_a;
    xArrivalT_n = n_a;
    xArrivalTe = mean_e;
    xArrivalTe_min = minOut_e;
    xArrivalTSigmoid = mean_s;
    xArrivalTSigmoid_min = minOut_s;
    cluster_count_x = static_cast<int>(clusters_ch.size());
    cluster_peak_count_x = peak_count;
  }

  if (XY % 2 == 1)
  {
    sig_y = true;
    y = position;
    hit_amp_y = charge;
    hit_strip_num_y = position_strip_count;
    y_nhits = static_cast<int>(main_ch.size());
    for (size_t s = 1; s < main_ch.size() && y_other.size() < 5; ++s)
      y_other.push_back(main_ch[s]);
    yMaxT = main_t[max_idx];
    yArrivalT = mean_a;
    yArrivalT_min = minOut_a;
    yArrivalT_n = n_a;
    yArrivalTe = mean_e;
    yArrivalTe_min = minOut_e;
    yArrivalTSigmoid = mean_s;
    yArrivalTSigmoid_min = minOut_s;
    cluster_count_y = static_cast<int>(clusters_ch.size());
    cluster_peak_count_y = peak_count;
  }

  return true;
}

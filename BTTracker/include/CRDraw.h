#ifndef CRDRAW_h
#define CRDRAW_h 1

#include "CRInfo.h"

class CRDraw
{
public:
  TFile *frawFile;
  TTree * frawTree;
  TFile *fanaFile;
  TTree *fanaTree[2 * MMdetN];
  CRDraw(TString fanaName);
  ~CRDraw();
  bool Init();
  void Draw_noise_es();
  void Draw_channel_waveform(int board, int XY, int channel);
  void Draw_channel_ped();
  void Draw_all_waveforms();
  void Cal_ped();
  bool Cut();
  bool select_waveform();
  void Draw_event_waveforms(int eventID);
  void Draw_events_from_file(const char* filename);


  private:
  TString anaName;


  /// input data
  int event;
  double maxA;
  int maxT;
  double ped;
  double pedrms;
  double riseT;
  double TOT;
  int channel;
  double channelmean[2 * MMdetN][256];
  double channelrms[2 * MMdetN][256];
  vector<short> *waveAna = nullptr;
  TBranch *b_event;
  TBranch *b_maxA;
  TBranch *b_maxT;
  TBranch *b_ped;
  TBranch *b_pedrms;
  TBranch *b_channel;
  TBranch *b_riseT;
  TBranch *b_TOT;
  TBranch *b_wave;
};

#endif
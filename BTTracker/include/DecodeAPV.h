#ifndef DECODEAPV_h
#define DECODEAPV_h 1

#include "CRInfo.h"

class DecodeAPV
{
public:
  TFile *fanaFile;
  TTree *fanaTree[MMdetN * 2];
  TFile *fdecFile;
  TTree *fdecTree[MMdetN];
  DecodeAPV(TString fanaName, TString fdecName, int force);
  ~DecodeAPV();
  bool Init();
  void Loop();
  void InitData();
  bool CalHit(vector<int> ch, vector<double> amp, vector<int> t, int XY);
  double fit_track(vector<double> px, vector<double> py, double par[3]);
  bool CalHitWithClusters(vector<int> ch, vector<double> amp, vector<int> t, vector<double> arrival, vector<double> arrivalE, vector<double> arrivalS, int XY);

private:
  TString anaName;
  TString decName;
  int force;
  vector<double> x_possible[MMdetN]; // X方向所有峰值位置
  vector<double> y_possible[MMdetN]; // Y方向所有峰值位置
  vector<double> x_energies[MMdetN]; // X方向所有峰值能量
  vector<double> y_energies[MMdetN]; // Y方向所有峰值能量

  int cluster_count_x;                 // X方向的簇团数量
  int cluster_count_y;                 // Y方向的簇团数量
  vector<double> cluster_positions_x;  // 所有X簇团的位置
  vector<double> cluster_positions_y;  // 所有Y簇团的位置
  vector<double> cluster_amplitudes_x; // 所有X簇团的幅度
  vector<double> cluster_amplitudes_y; // 所有Y簇团的幅度
  int cluster_peak_count_x;  // X方向主要簇团中的极大值个数
  int cluster_peak_count_y;  // Y方向主要簇团中的极大值个数
  /// input data
  int event;
  Long64_t entry;
  double maxA;
  int maxT;
  double arrivalT;
  double arrivalT_e;
  double arrivalT_sigmoid;
  int channel;
  TBranch *b_event;
  TBranch *b_maxA;
  TBranch *b_maxT;
  TBranch *b_channel;
  TBranch *b_timestamp;
  /// ouput data
  bool sig;
  bool sig_x;
  bool sig_y;
  double x;
  double y;
  double z;
  double hit_amp_x;
  double hit_amp_y;
  int hit_strip_num_x;
  int hit_strip_num_y;
  int x_nhits;
  int y_nhits;
  int xMaxT;
  int yMaxT;
  double xArrivalT;
  double yArrivalT;
  double xArrivalT_min;
  double yArrivalT_min;
  int xArrivalT_n;
  int yArrivalT_n;
  double xArrivalTe;
  double yArrivalTe;
  double xArrivalTe_min;
  double yArrivalTe_min;
  double xArrivalTSigmoid;
  double yArrivalTSigmoid;
  double xArrivalTSigmoid_min;
  double yArrivalTSigmoid_min;
  vector<double> x_other;
  vector<double> y_other;
  long long timestamp;
};

#endif

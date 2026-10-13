#ifndef RECONSTRUCTTRACK_h
#define RECONSTRUCTTRACK_h 1

#include "CRInfo.h"
#include "Mille.h"
#include "TMatrixD.h"
#include "TRotation.h"

#include <array>
#include <map>
#include <vector>

class ReconstructTrack
{
public:
  TFile *fdecFile;
  TTree *fdecTree[MMdetN];
  TFile *frecFile;
  TTree *frecTree;
  ReconstructTrack(TString fdecName, TString frecName, int force);
  ~ReconstructTrack();
  bool Init();
  void Loop();
  void InitData();

  /// Alignment parameters of the two BM pad planes (chip12/13), fitted
  /// against tracker-projected tracks. Index: 0 = chip12, 1 = chip13.
  double bmDx[2] = {0.0, 0.0};
  double bmDy[2] = {0.0, 0.0};
  bool bmAligned = false;

  /// for efficiency cut

  //// track selection, correction, reconstruction
  bool is_single_correct = true;   // if correct single hit
  double centre = 256.0 / 2 * 0.4; // layer centre, mm
  double fit_R2_cut = 0.99;          // correlation index square
  bool fit_track(vector<double> px, vector<double> py, double par[3]);
  bool fit_track3(vector<double> px, vector<double> py, double par[3]);
  bool fit_track2(vector<double> px, vector<double> py, double par[3]);
  double hit_reconstruction(vector<int> id, vector<double> px, vector<double> &py, int xy);
  double hit_selection(int eraseN, vector<int> id, vector<double> px, vector<double> &py, int xy);
  bool track_reconstruction(vector<double> px, vector<double> py, vector<double> pz);
  bool GetBTOFposition(double kx, double ky, double bx, double by, vector<double> botfxyz);
  bool GetRICHposition(double kx, double ky, double bx, double by, vector<double> richxyz);

  /// asic channel to detector channel
  vector<int> APV2det[256];
  bool get_APV2det_lists(string lists_name);
  bool cal_APV2det_lists();
  /// detector alignment
  double alg_fit_R2_cut = 0.05;              // correlation index squre
  double alg_theta_cut = 20;                 // degree, cosmic tracks are inclined
  double alg_fit_res_cut = 20;               // mm, loose pre-alignment residual
  double alg_res = 5e-8;                     // delta-residual resolution
  double alg_rho = 0.95, alg_epsilon = 1e-4; // for alginment loop
  vector<string> ParName = {"deltaX", "deltaY", "deltaZ", "thetaX", "thetaY", "thetaZ"};
  double algPar3[MMdetN][6];
  double algPar6[MMdetN][6];
  TString alignmentPar3Name;
  TString alignmentPar6Name;
  bool JointAlignment(const string &par3Filename,
                      const string &par6Filename);
  void ReadAlignmentPars(string filename, int Npars = 6);
  void Alignment(string filename, int Npars = 6);

  void ReadSecondAlignmentPars(string filename, int Npars = 6);
  void SecondAlignment(string filename, int Npars = 6);
  double SecondalgPar3[2][6];
  double SecondalgPar6[2][6];

  void SecondAlignmentMillepede(string filename, int Npars);

  double alignment_list[MMdetN][6];
  TMatrixD X, XT, Y, YT, Z, ZT;
  TMatrixD r0[MMdetN], r0T[MMdetN], R[MMdetN], RT[MMdetN];
  TMatrixD Deltar[MMdetN], DeltarT[MMdetN], DeltaRx[MMdetN], DeltaRxT[MMdetN],
      DeltaRy[MMdetN], DeltaRyT[MMdetN], DeltaRz[MMdetN], DeltaRzT[MMdetN],
      DeltaR[MMdetN], DeltaRT[MMdetN], DRx[MMdetN], DRxT[MMdetN], DRy[MMdetN],
      DRyT[MMdetN], DRz[MMdetN], DRzT[MMdetN], DR[MMdetN], DRT[MMdetN];
  TMatrixD k, kT, b, bT;
  string m_alignment_filename = "";

  void Millepede_ali();
  void Write_mille_file();
  void define_align_matrix();
  void define_const_matrix();
  void define_origin_matrix();
  void set_alignment_filename(string filename);
  int local_to_global(double &xx, double &yy, double &zz, int i);
  void define_track_matrix(double kx, double ky, double bx, double by);
  int read_millepede_res(string filename, int alg = 2);
  bool alignment_derlcx(float *, int);
  bool alignment_derlcy(float *, int);
  bool alignment_derglx(float *, int);
  bool alignment_dergly(float *, int);
  bool alignment_derlcx(double *, int);
  bool alignment_derlcy(double *, int);
  bool alignment_derglx(double *, int);
  bool alignment_dergly(double *, int);
  double expect_position_x(int i);
  double expect_position_y(int i);
  void print_alignment_list();
  bool least_square_fit(vector<double> z, vector<double> x, double sigma, double &kx, double &bx,
                        double &chi2Dndf);

  /// BM pad plane aggregation and alignment against tracker-projected tracks.
  struct BMCluster
  {
    double x = 0.0;
    double y = 0.0;
    double charge = 0.0;
    int npads = 0;
  };
  bool BuildBMAggregation();
  bool AlignBMPads(TString bmAlignmentName);
  void ReadBMAignmentPars(TString bmAlignmentName);

  double get_data_chi2();

private:
  /// Aggregated BM clusters indexed by (event, chip index).
  std::map<int, std::array<std::vector<BMCluster>, 2>> bmClusters;

  TString decName;
  TString recName;
  int force;
  /// input data
  int event;
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
  int cluster_count_x;
  int cluster_count_y;
  int xMaxT;
  int yMaxT;
  vector<double> *x_other = 0;
  vector<double> *y_other = 0;
  TBranch *b_event;
  TBranch *b_sig;
  TBranch *b_sig_x;
  TBranch *b_sig_y;
  TBranch *b_x;
  TBranch *b_y;
  TBranch *b_z;
  TBranch *b_hit_amp_x;
  TBranch *b_hit_amp_y;
  TBranch *b_hit_strip_num_x;
  TBranch *b_hit_strip_num_y;
  TBranch *b_x_nhits;
  TBranch *b_y_nhits;
  TBranch *b_cluster_count_x;
  TBranch *b_cluster_count_y;
  TBranch *b_x_other;
  TBranch *b_y_other;
  TBranch *b_timestamp;
  TBranch *b_xMaxT;
  TBranch *b_yMaxT;
  /// output data
  double MMx[MMdetN];
  double MMy[MMdetN];
  double MMz[MMdetN];
  bool MMsigx[MMdetN];
  bool MMsigy[MMdetN];
  double XfitPars[3];
  double YfitPars[3];
  double X01fitPars[3];
  double Y01fitPars[3];
  double X2345fitPars[3];
  double Y2345fitPars[3];
  double predStripX[MMdetN];
  double predStripY[MMdetN];
  double twotrackangle;
  bool useflag;
  double richextraHit[3];
  double richtrackVec[3];
  bool richuseflag;
  double btofextraHit[3];
  double btoftrackVec[3];
  bool btofuseflag;
  double ecalextraHit[3];
  double ecaltrackVec[3];
  bool ecaluseflag;
  long long timestamp;

  double m_Ed[MMdetN][6] = {0};    // 更新量平方的移动平均
  double m_ED[MMdetN][6] = {0};    // 梯度平方的移动平均
  double m_D_pre[MMdetN][6] = {0}; // 上一次的梯度（用于符号检测）
  double m_rho = 0.95;              // 衰减因子（例如 0.9）
  double m_epsilon = 1e-4;          // 小常数，可随迭代衰减
};

#endif

#ifndef CRINFO_h
#define CRINFO_h 1

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <bitset>

#include "TSystem.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TBenchmark.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TRandom3.h"

#include "CRConst.h"

using namespace std;

const int eventrange = 0x10000; /// 65536
const TString asic2detName = "config/asic2det_64to384.csv";
const TString MMalignmentPar3Name = "config/MMalignmentPar3.dat";
const TString MMalignmentPar6Name = "config/MMalignmentPar6.dat";
// Active setup: three strip-readout MM Tracker planes. The former fourth
// Tracker is not installed in this APV25 variant.
const int MMdetN = 3;
// Decoded ROOT source apv_board layout (normalised by
// DecodedAPVRoot2TrackRaw): 0/1=T1X, 2/3=T1Y, 4/5=T2X, 6/7=T2Y,
// 8=BM1, 9=BM2, 10/11=T3X, 12/13=T3Y, 14/15=unloaded.
const int trackerReadoutBoardN = 2;
const int trackersPerReadoutBoard = 2;
const int apvChipsPerTracker = 4;
const int apvChannelsPerChip = 128;
const int bmDetectorN = 2;
const int bmChannelsPerDetector = 121;
const int trackerGlobalChannelCount = 2 * MMdetN * 256;
const int globalNoiseChannelCount =
    trackerGlobalChannelCount + bmDetectorN * bmChannelsPerDetector;

const double BETA = 3000. / TMath::Sqrt(3000. * 3000. + Mass_mu * Mass_mu); // p = 3 GeV/c
const char detName[MMdetN][10] = {"MM0", "MM1", "MM2"};
const double Xoff[MMdetN] = {0, 0, 0};
const double Yoff[MMdetN] = {0, 0, 0};
// const double Yoff[NDET + 2] = {90, 90, -35, -35, 45, 0}; // 20 degree angle
//  old position
// const double Zpos[6] = {582, 538, -120, -164, 320, 26}; // MM-0-3, T0, DTOF
// position for beam test setup==0
// const double Zpos[7] = {2207, 2067, -3947, -4147, 3357, -4922, 0}; // MM-0-3, T01,T02, DTOF
// position for beam test setup==1
//const double Zpos[7] = {0, -200, -870, -1070, 3357, -4922, 0}; // MM-0-3, T01,T02, DTOF
//const double Zpos[7] = {1224, 1424, 2852, 3052, 6205, 6405, 0};
//const double Zpos[7] = {2490, 2690, 5059, 5259, 6205, 6405, 0};
//const double Zpos[7] = {1224, 1424, 5059, 5259, 6205, 6405, 0};
//const double Zpos[7] = {3111, 3311, 4968, 5168, 6199, 6399, 0};
//const double Zpos[7] = {5059, 5259, 6205, 6405, 0, 0, 0};
//const double Zpos[7] = {3111, 3311, 6199, 6399, 0, 0, 0};
//const double Zpos[7] = {0, -200, -870, -1070, 3357, -4922, 0};
//const double Zpos[7] = {0, 251, 316, 381, 383, 448, 0};

// Preserve the original BTTracker geometry, measured from the lower
// (downstream) PMT plane with z increasing upward. The original Tracker1..4
// positions were 533, 468, 217, 154 mm; only Tracker4 at 154 mm is removed.
// BM1/BM2 retain their original positions between Tracker2 and Tracker3.
const double Zpos[MMdetN] = {533, 468, 217};
const double BMZpos[bmDetectorN] = {364, 309};
// const double Zpos[6] = {2989, 3189, 5377, 5577, 7052, 7252};
//const double Zpos[6] = {2373, 2573, 2912, 3112, 7052, 7252};
//const double Zpos[7] = {5377, 5577, 6568, 6768, 0, 0, 0};
const double rmscut = 5;                                       // for MM waveform cut
const int max_Nch_use = 10;                                    // for MM dec

const int color[21] = {1, 2, 3, 4, 6, 7, 9, 30, 50, 51, 28, 38, 41, 39, 42, 20, 40, 5, 8, 29, 33};

double Np(double WaveLength); // unit of WaveLength: nm
double Ng(double WaveLength); // unit of WaveLength: nm

TString GetPath(TString fileName);
TString GetFilename(TString fileName);
void GetFileList(TString filePath, TString filePattern, vector<TString> &fList);
ofstream LogFile(TString fileName);
void DrawProcessbar(int i, int total);
int IsMatched(int trigid, vector<int> idvec);
int IsMatched(double RecT, vector<double> RecTs);
int IsMatched(int RecT, vector<int> RecTs, int index);
TString XYname(int findex);
TLatex *CreatLatex(TString text, double x = 0.65, double y = 0.5, int textFont = 42);

Double_t CalVectorMean(vector<double> Time, int nn = 0, int bb = 0);
void SortRecT(vector<double> &TOAevent);
double SelectRecT(vector<double> &RecT, double deltaT = 1500, double w = 1);
#endif


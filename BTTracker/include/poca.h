#ifndef poca_h
#define poca_h 1

#include "CRInfo.h"
#include <iostream>
#include <cmath>
#include <TMath.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TVectorD.h>

class POCACalculator
{
public:
    POCACalculator();
    POCACalculator(TString f_frecName, TString f_fpocaName);
    POCACalculator(vector<TString> f_pocaslist, TString f_fpocaName);
    ~POCACalculator();
    bool Init();
    void InitData();
    void Loop();
    struct TrackParams 
    {
        double kx, bx, ky, by;
        TrackParams(double kx_, double bx_, double ky_, double by_) 
            : kx(kx_), bx(bx_), ky(ky_), by(by_) {}
    };
    struct POCAResult 
    {
        TVector3 poca_point;
        double distance;
        TVector3 point_on_track1;
        TVector3 point_on_track2;
        double scattering_angle; // 散射角（弧度）
        
        POCAResult() : distance(0), scattering_angle(0) {}
    };
    POCAResult CalculatePOCA(const TrackParams& track1, const TrackParams& track2);
    double CalculateScatteringAngle(const TrackParams& track1, const TrackParams& track2);
    void CalculateScatteringComponents(const TrackParams& track1, const TrackParams& track2, double& angle_x, double& angle_y);
    POCAResult HandleParallelTracks(const TrackParams& track1, const TrackParams& track2);

private:
    TString fRecName;
    TString fpocaName;
    TFile *fRecFile[4];
    TTree *fRecTree[4];
    TFile *fpocaFile;
    TTree *fpocaTree;

    vector<TString> fpocalist;

    int event;
    bool btofuseflag;
    bool useflag;
    double X01fitPars[3];
    double Y01fitPars[3];
    double X2345fitPars[3];
    double Y2345fitPars[3];
    double XfitPars[3];
    double YfitPars[3];
    double pocax;
    double pocay;
    double pocaz;
    double scattering_angle;
};



#endif
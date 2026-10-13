#ifndef ReadTrigger_h
#define ReadTrigger_h 1

#include "CRInfo.h"

class ReadTrigger
{
public:
    ReadTrigger();
    ReadTrigger(TString f_fTriggername, TString f_Triggerrootname, TString f_frecName, TString f_fmRecName, int f_force);
    ~ReadTrigger();

    bool TriggerInit();
    bool Init();
    void ReadTriggerDat();
    void Loop();
    void InitData();


private:
    TString fRecName;
    TFile* fRecFile;
    TTree* fRecTree;
    TString fTriggerName;
    TFile* fTriggerFile;
    TTree* fTriggerTree;
    TString fmRecName;
    TFile* fmRecFile;
    TTree* fmRecTree;
    TString fTriggerDat;
    TString fTriggerrootname;
    TFile* readtriggerFile;
    TTree* readtriggerTree;
    

    int force;


    int triggerID;
    long long triggerTimestamp;

    int event;
    double MMx[MMdetN];
    double MMy[MMdetN];
    double MMz[MMdetN];
    bool MMsigx[MMdetN];
    bool MMsigy[MMdetN];
    double XfitPars[3];
    double YfitPars[3];
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

    unsigned short exchangehighlow(unsigned short data);
};

#endif
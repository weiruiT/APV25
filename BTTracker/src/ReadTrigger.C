#include "ReadTrigger.h"
#include "CRInfo.h"

ReadTrigger::ReadTrigger()
{
    ;
}

ReadTrigger::ReadTrigger(TString f_fTriggername, TString f_Triggerrootname, TString f_fRecName, TString f_fmRecName, int f_force)
{
    fTriggerDat = f_fTriggername;
    fTriggerrootname = f_Triggerrootname;
    fRecName = f_fRecName;
    fmRecName = f_fmRecName;
    force = f_force;
}

ReadTrigger::~ReadTrigger()
{
}

bool::ReadTrigger::TriggerInit()
{
    fTriggerFile = new TFile(fTriggerrootname,"RECREATE");
    fTriggerTree = new TTree("trigger", "trigger");
    fTriggerTree->Branch("triggerID", &triggerID, "triggerID/I");
    fTriggerTree->Branch("triggerTimestamp", &triggerTimestamp, "triggerTimestamp/L");
    return true;
}

bool ReadTrigger::Init()
{
    // data in
    fRecFile = new TFile(fRecName, "READ");
    fRecTree = (TTree *)fRecFile->Get("Track");
    fRecTree->SetBranchAddress("event", &event);
    fRecTree->SetBranchAddress("MMx", MMx);
    fRecTree->SetBranchAddress("MMy", MMy);
    fRecTree->SetBranchAddress("MMz", MMz);
    fRecTree->SetBranchAddress("MMsigx", MMsigx);
    fRecTree->SetBranchAddress("MMsigy", MMsigy);
    fRecTree->SetBranchAddress("XfitPars", XfitPars);
    fRecTree->SetBranchAddress("YfitPars", YfitPars);
    fRecTree->SetBranchAddress("useflag", &useflag);
    fRecTree->SetBranchAddress("richuseflag", &richuseflag);
    fRecTree->SetBranchAddress("richextraHit", richextraHit);
    fRecTree->SetBranchAddress("richtrackVec", richtrackVec);
    fRecTree->SetBranchAddress("btofuseflag", &btofuseflag);
    fRecTree->SetBranchAddress("btofextraHit", btofextraHit);
    fRecTree->SetBranchAddress("btoftrackVec", btoftrackVec);
    fRecTree->SetBranchAddress("ecaluseflag", &ecaluseflag);
    fRecTree->SetBranchAddress("ecalextraHit", ecalextraHit);
    fRecTree->SetBranchAddress("ecaltrackVec", ecaltrackVec);
    fRecTree->SetBranchAddress("timestamp", &timestamp);
    // data in
    readtriggerFile = new TFile(fTriggerrootname, "READ");
    readtriggerTree = (TTree *)readtriggerFile->Get("trigger");
    readtriggerTree->SetBranchAddress("triggerID", &triggerID);
    readtriggerTree->SetBranchAddress("triggerTimestamp", &triggerTimestamp);
    // data out
    fmRecFile = new TFile(fmRecName,"RECREATE");
    fmRecTree = new TTree("Track_m", "Track_m");
    fmRecTree->Branch("event", &event, "event/I");
    fmRecTree->Branch("MMx", MMx, Form("MMx[%d]/D", MMdetN));
    fmRecTree->Branch("MMy", MMy, Form("MMy[%d]/D", MMdetN));
    fmRecTree->Branch("MMz", MMz, Form("MMz[%d]/D", MMdetN));
    fmRecTree->Branch("MMsigx", MMsigx, Form("MMsigx[%d]/O", MMdetN));
    fmRecTree->Branch("MMsigy", MMsigy, Form("MMsigy[%d]/O", MMdetN));
    fmRecTree->Branch("XfitPars", XfitPars, "XfitPars[3]/D");
    fmRecTree->Branch("YfitPars", YfitPars, "YfitPars[3]/D");
    fmRecTree->Branch("useflag", &useflag, "useflag/O");
    fmRecTree->Branch("richuseflag", &richuseflag, "richuseflag/O");
    fmRecTree->Branch("richextraHit", richextraHit, "richextraHit[3]/D");
    fmRecTree->Branch("richtrackVec", richtrackVec, "richtrackVec[3]/D");
    fmRecTree->Branch("btofuseflag", &btofuseflag, "btofuseflag/O");
    fmRecTree->Branch("btofextraHit", btofextraHit, "btofextraHit[3]/D");
    fmRecTree->Branch("btoftrackVec", btoftrackVec, "btoftrackVec[3]/D");
    fmRecTree->Branch("ecaluseflag", &ecaluseflag, "ecaluseflag/O");
    fmRecTree->Branch("ecalextraHit", ecalextraHit, "ecalextraHit[3]/D");
    fmRecTree->Branch("ecaltrackVec", ecaltrackVec, "ecaltrackVec[3]/D");
    fmRecTree->Branch("timestamp", &timestamp, "timestamp/L");
    return true;
}

void ReadTrigger::ReadTriggerDat()
{
    fstream InDat;
    InDat.open(fTriggerDat, ios::in | ios::binary);
    if (!InDat.is_open())
    {
        cout << " ERROR! " << fTriggerDat << " cant open!" << endl;
        return;
    }
    cout << " --> Reading TriggerDatnn: " << fTriggerDat << endl;
    unsigned short memblock;
    int length = sizeof(short);
    int loopnum = 0;
    int rawtriggerid = -1;
    int pretriggerid = -1;
    while (!InDat.eof())
    {
        InDat.read((char *)(&memblock), length);
        if (memblock == 0x8103)
        {
            InDat.read((char *)(&memblock), length);
            if (memblock == 0x101F)
            {
                InDat.read((char *)(&memblock), length);
                rawtriggerid = memblock;
                if (rawtriggerid - pretriggerid < -10) loopnum++;
                pretriggerid = rawtriggerid;
                triggerID = rawtriggerid + 65536 * loopnum;
                for (int i = 0; i < 9; i++)
                {
                    InDat.read((char *)(&memblock), length);
                }
                InDat.read((char *)(&memblock), length);
                Long64_t  rawtime = 0;
                Long64_t  midtime = 0;
                rawtime = memblock;
                InDat.read((char *)(&memblock), length);
                midtime = memblock;
                rawtime = midtime * 16 * 16 * 16 * 16 + rawtime;
                InDat.read((char *)(&memblock), length);
                midtime = memblock;
                rawtime = midtime * TMath::Power(16, 8) + rawtime;
                triggerTimestamp = rawtime;

                fTriggerTree->Fill();
            }
        }
    }
    fTriggerFile->WriteTObject(fTriggerTree);
    fTriggerFile->Close();
}

void ReadTrigger::Loop()
{
    int cha = 0;
    for(int i = 0; i < fRecTree->GetEntries(); i++)
    {
        fRecTree->GetEntry(i);
        fTriggerTree->GetEntry(i + cha);
        if (timestamp == 0)
        {
            InitData();
            fmRecTree->Fill();
            continue;
        }
        if (timestamp - triggerTimestamp < 6000) 
        {
            event = event + cha;
            fmRecTree->Fill();
            //cout << event << "," << triggerID << "," << timestamp - triggerTimestamp << endl;
        }
        else
        {
            event = triggerID;
            InitData();
            fmRecTree->Fill();
            cha++;
            i--;
        }
    }
    fmRecFile->WriteTObject(fmRecTree);
    fmRecFile->Close();
}

void ReadTrigger::InitData()
{
  useflag = 1; // useful
  richuseflag = 0;
  btofuseflag = 0;
  ecaluseflag = 0;
  memset(MMx, 0, sizeof(MMx));
  memset(MMy, 0, sizeof(MMy));
  memset(MMz, 0, sizeof(MMz));
  memset(MMsigx, 0, sizeof(MMsigx));
  memset(MMsigy, 0, sizeof(MMsigy));
  memset(XfitPars, 0, sizeof(XfitPars));
  memset(YfitPars, 0, sizeof(YfitPars));
  memset(richextraHit, 0, sizeof(richextraHit));
  memset(richtrackVec, 0, sizeof(richtrackVec));
  memset(btofextraHit, 0, sizeof(btofextraHit));
  memset(btoftrackVec, 0, sizeof(btoftrackVec));
  memset(ecalextraHit, 0, sizeof(ecalextraHit));
  memset(ecaltrackVec, 0, sizeof(ecaltrackVec));
}

unsigned short ReadTrigger::exchangehighlow(unsigned short data)
{

    unsigned short low = data & 0x00ff;
    unsigned short high = (data >> 8) & 0x00ff;
    unsigned short aa = low;
    aa <<= 8;
    aa |= high;
    return aa;
}
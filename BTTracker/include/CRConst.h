#ifndef CRCONST_h
#define CRCONST_h 1
// some parameters
///////////  

#include "TMath.h"

const double deg=TMath::Pi()/180;
const double mm=1;
const double ns=1;
const double MeV=1;
const double cm=10;
const double NpMax = 1.49602; // 275 nm
const double T0Np = 1.47865; // 340 nm
const double T0Ng = 1.54420; // 340 nm
const double DTOFNp = 1.47454; // 365 nm
const double DTOFNg = 1.52994; // 365 nm
const double Mass_mu = 105.6583715*MeV;


//// Z position
const double T0Z = 515*mm;
const double DTOFZ = 57.5*mm;
const double MMZ[6] = {300*mm, 220*mm,40*mm, -40*mm, -220*mm, -300*mm};
const double TRIGZ[2] = {889*mm, -475*mm}; // triger
const double LEADZ = -343*mm;

///// tracker
const double TriggerL = 22*cm;
const double TriggerW = 22*cm;
const double TriggerT = 1.5*cm;
///// lead brick
const double LeadL = 22*cm;
const double LeadW = 22*cm;
const double LeadT = 5.5*cm;

///// tracker
const double TrackerL = 20*cm;
const double TrackerW = 20*cm;
const double TrackerT = 2*cm;
const double TrackerAL = 153.6*mm;
const double TrackerAT = 0.1*mm;
const double TrackerPCBL = 200*mm;
const double TrackerPCBT = 3*mm;

////// T0
const double T0L = 22*cm;
const double T0W = 22*cm;
const double T0T = 5*cm;
const double T0BT = 3*mm; /// box thick

const double T0QL = 18*cm;  /// quartz
const double T0QW = 18*cm; ///quartz
const double T0QT = 1*cm; ///  quartz

const double T0LGR = 5*mm;  // light guide
const double T0LGT = 5*mm; ///  light guide
const double T0ST = 0.1*mm; //// sensor
const double T0SGT = 0.1*mm; //silica grease
const double T0AT = 0.1*mm; //absorber

////// DTOF
const double DTOFBT = 8*mm; /// box thick
const double DTOFBG = 2*mm; /// gap beween box and quartz

const double DTOFQH = 454*mm; // quartz high
const double DTOFQW2 = 533*mm; // quartz outer side length
const double DTOFQW1 = DTOFQW2-DTOFQH*TMath::Tan(15.*deg)*2; // quartz inner side length
const double DTOFQT = 15*mm; // quartz thick
const double DTOFQP = 240.13*mm; //quartz location, distance between system center and quatrzW2

const double DTOFSL = 27.6*mm; // sensor size
const double DTOFST = 13.1*mm; // sensor thick
const double DTOFSI = 7.6*mm; //sensor interval
const double DTOFSP = 24*mm; //sensor location, distance between sensor center and quatrzW2

const double DTOFAL = 5.28*mm; //anode size
const double DTOFAI = 0.3*mm; //anode interval
const double DTOFAT = 0.1*mm; //anode thick
const double DTOFSGT = 0.1*mm; //silica grease
const double DTOFAbT = 0.1*mm; //absorber

const int DTOFSN = 42; /// sensor number

#endif

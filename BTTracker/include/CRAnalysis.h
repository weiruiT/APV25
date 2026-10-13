#ifndef CRANALYSIS_h
#define CRANALYSIS_h 1

#include "CRInfo.h"

bool CalibrateTrackerNoise(TString fileDir, int force = 0);
bool ReadTrackerConfigured(TString fileDir, int force = 0);
void MakeAmplitudeDistributions(TString fileDir);
void MakeEnergySpectra(TString fileDir);
void CompareArrivalTimes(TString fileDir);
#endif

#ifndef DECODEDAPVROOT2TRACKRAW_H
#define DECODEDAPVROOT2TRACKRAW_H 1

#include "CRInfo.h"

// Adapter for the already-decoded APV25 ROOT files produced by
// data_reduction_baseline_v10.C. It normalises the input schema and hardware
// numbering only; it does not subtract a pedestal or apply another threshold.
class DecodedAPVRoot2TrackRaw
{
public:
  DecodedAPVRoot2TrackRaw(const vector<TString> &inputFiles,
                          TString outputName, int force = 0,
                          bool preprocessedSignal = true,
                          int bmSlot8LogicalAPV = 12,
                          int bmSlot9LogicalAPV = 13);
  bool Init();
  bool Loop();

private:
  bool ConvertFile(const TString &inputName, int fileIndex,
                   long long &lastOutputEvent);
  int LogicalAPV(int sourceAPV) const;

  vector<TString> inputs;
  TString outputName;
  int force = 0;
  bool preprocessedSignal = true;
  int bmSlot8LogicalAPV = 12;
  int bmSlot9LogicalAPV = 13;
  TFile *outputFile = nullptr;
  TTree *outputTree = nullptr;

  int event = -1;
  int board = -1;
  int readoutBoard = 1;
  int chip = -1;
  int channel = -1;
  long long timestamp = -1;
  vector<short> wave;

  int sourceFile = -1;
  int sourceEntry = -1;
  int sourceRawEventIndex = -1;
  int sourceAPVBoard = -1;
  unsigned int sourceTrigger = 0;

  long long inputFrames = 0;
  long long trackerFrames = 0;
  long long padFrames = 0;
  long long ignoredFrames = 0;
  long long invalidChannels = 0;
  long long outputWaveforms = 0;
};

#endif

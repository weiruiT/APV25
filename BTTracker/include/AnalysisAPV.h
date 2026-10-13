#ifndef AnalysisAPV_h
#define AnalysisAPV_h 1

#include "CRInfo.h"

struct AnalysisBMMapRow
{
  int board = -1;
  int chip = -1;
  int channelId = -1;
  int position = -1;
};

struct SignalSelectionConfig
{
  double amplitudeMin = 0.0;
  double maxTMin = 0.0;
  double maxTMax = 15.0;
  double riseTimeMin = 0.0;
  double riseTimeMax = 10.0;
  double totMin = 2.0;
  double totMax = 14.0;
  double peakThresholdSigma = 5.0;
  bool requireSinglePeak = true;
  bool requireValidSigmoid = true;
};



class AnalysisAPV
{
  ///// process wave to get ped, pedrms, maxA
  public:
    TFile *frawFile;
    TTree * frawTree;
    TFile *fanaFile;
    TTree * fanaTree[2 * MMdetN];
    AnalysisAPV(TString fRawName, TString fAnaName, int force=0,
                const SignalSelectionConfig &selection = SignalSelectionConfig(),
                TString fNoiseParameterName = "");
    ~AnalysisAPV();
    bool Cut();
    bool Init();
    void Loop();
    int IndexXY(int board, int chip);
    int IndexDet(int board, int chip);
    int IndexStrip(int board,int chip,int channel);
    int CalculateLostID(vector<int> eventvec);
    void InitData();
    void riseTandTOT(vector<double> waveform);
    void CalWavePed();
    void CalWaveInfo();
    void Cal_ped();
    double FitPulseArrivalTime(); // Pause a/e retained, calls disabled in sigmoid-only mode.
    double FitSigmoidArrivalTime();
    bool SelectPulseWaveform();
  private:
    static constexpr int kGlobalNoiseChannels = globalNoiseChannelCount;
    TString rawName;
    TString anaName;
    SignalSelectionConfig selection;
    // Run-local 1778-channel text calibration; only sigma enters peak finding.
    TString noiseParameterName;
    int force;
    /// input data
    int event;
    int board;
    int chip;
    int channel;
    vector<short> waveAna;
    vector<short> *wave = nullptr;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;
    TBranch *b_timestamp;
    /// output data
    double channelped[2 * MMdetN][256];
    double channelrms[2 * MMdetN][256];
    double globalPedestal[kGlobalNoiseChannels];
    double globalRms[kGlobalNoiseChannels];
    int currentGlobalChannel = -1;
    bool currentIsBM = false;
    bool bmNoiseCalibrationLoaded = false;
    int currentRawChip = -1;
    int currentRawChannel = -1;
    int currentPadPosition = -1;
    TTree *bmTree = nullptr;
    bool channelCalibrationLoaded = false;
    double maxA;
    int maxT;
    double ped;
    double pedrms;
    double peakNoiseRms;
    double peakThresholdAdc;
    int nisel;
    double riseT;
    double TOT;
    double arrivalT;
    double arrivalT_a;
    double arrivalT_e;
    double arrivalT_sigmoid;
    double sigmoidT50;
    double sigmoidBaseline;
    double sigmoidHeight;
    double sigmoidK;
    double sigmoidFitChi2NDF;
    int sigmoidFitStatus;
    int sigmoidFitStart;
    int sigmoidFitEnd;
    int pulseFitStatus;
    double pulseFitChi2NDF;
    int pulsePeakFailureReason;
    int pulsePeakCount = -1; // Significant waveform maxima, not ionisation clusters.
    int XY;
    int Det; 
    long long timestamp=0;
    vector<int> APVmap;
    vector<AnalysisBMMapRow> bmMap;

    bool LoadNoiseParameters();
    bool LoadBMMap();
    int GlobalChannelForRaw(int rawBoard, int rawChip, int rawChannel,
                            int &mappedChannel, int &padPosition);

};

#endif



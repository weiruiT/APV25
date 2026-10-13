#include "DecodedAPVRoot2TrackRaw.h"

#include "TNamed.h"
#include "TLeaf.h"

#include <algorithm>
#include <limits>

namespace
{
constexpr int kInputChannelCount = 128;
constexpr int kSamplesPerWaveform = 30;

short ToShort(int value)
{
  return static_cast<short>(std::max(
      static_cast<int>(std::numeric_limits<short>::min()),
      std::min(static_cast<int>(std::numeric_limits<short>::max()), value)));
}
}

DecodedAPVRoot2TrackRaw::DecodedAPVRoot2TrackRaw(
    const vector<TString> &inputFiles, TString name, int overwrite,
    bool signalInput, int slot8LogicalAPV, int slot9LogicalAPV)
    : inputs(inputFiles), outputName(name), force(overwrite),
      preprocessedSignal(signalInput), bmSlot8LogicalAPV(slot8LogicalAPV),
      bmSlot9LogicalAPV(slot9LogicalAPV)
{
}

int DecodedAPVRoot2TrackRaw::LogicalAPV(int sourceAPV) const
{
  // Source apv_board slots from data_reduction_baseline_v10.C:
  //   0/1=T1X, 2/3=T1Y, 4/5=T2X, 6/7=T2Y,
  //   8=BM1, 9=BM2, 10/11=T3X, 12/13=T3Y, 14/15=unloaded.
  // Move BM1/BM2 to the internal chip12/13 namespace and compact Tracker3
  // to logical APVs 8..11. Slots 14/15 remain outside every accepted range
  // and are counted as ignored frames by ConvertFile().
  const int mapping[16] = {
      0, 1, 2, 3, 4, 5, 6, 7,
      bmSlot8LogicalAPV, bmSlot9LogicalAPV,
      8, 9, 10, 11, 14, 15};
  if (sourceAPV < 0 || sourceAPV >= 16)
    return -1;
  return mapping[sourceAPV];
}

bool DecodedAPVRoot2TrackRaw::Init()
{
  if (!((bmSlot8LogicalAPV == 12 && bmSlot9LogicalAPV == 13) ||
        (bmSlot8LogicalAPV == 13 && bmSlot9LogicalAPV == 12)))
  {
    cerr << "ERROR: BM slot mapping must be either 8->12, 9->13 or "
            "8->13, 9->12."
         << endl;
    return false;
  }
  if (inputs.empty())
  {
    cerr << "ERROR: no decoded APV ROOT input files." << endl;
    return false;
  }
  if (force != 1 && !gSystem->AccessPathName(outputName))
  {
    cerr << "ERROR: output already exists: " << outputName << endl;
    return false;
  }

  outputFile = TFile::Open(outputName, "RECREATE");
  if (!outputFile || outputFile->IsZombie())
  {
    cerr << "ERROR: cannot create " << outputName << endl;
    return false;
  }

  outputTree = new TTree("TrackRaw", "Normalised preprocessed APV25 waveforms");
  outputTree->Branch("event", &event, "event/I");
  outputTree->Branch("board", &board, "board/I");
  outputTree->Branch("readoutBoard", &readoutBoard, "readoutBoard/I");
  outputTree->Branch("chip", &chip, "chip/I");
  outputTree->Branch("channel", &channel, "channel/I");
  outputTree->Branch("wave", &wave);
  outputTree->Branch("timestamp", &timestamp, "timestamp/L");
  outputTree->Branch("source_file", &sourceFile, "source_file/I");
  outputTree->Branch("source_entry", &sourceEntry, "source_entry/I");
  outputTree->Branch("source_raw_event_index", &sourceRawEventIndex,
                     "source_raw_event_index/I");
  outputTree->Branch("source_apv_board", &sourceAPVBoard,
                     "source_apv_board/I");
  outputTree->Branch("source_trigger", &sourceTrigger, "source_trigger/i");
  return true;
}

bool DecodedAPVRoot2TrackRaw::ConvertFile(const TString &inputName,
                                          int fileIndex,
                                          long long &lastOutputEvent)
{
  TFile input(inputName, "READ");
  if (input.IsZombie())
  {
    cerr << "ERROR: cannot open decoded APV ROOT input " << inputName << endl;
    return false;
  }

  TTree *samplesTree = dynamic_cast<TTree *>(input.Get("apv_EventSamples"));
  TTree *headersTree = dynamic_cast<TTree *>(input.Get("apv_EventHeaders"));
  bool legacyNoise = false;
  if (!samplesTree && !preprocessedSignal)
  {
    samplesTree = dynamic_cast<TTree *>(input.Get("apv_events"));
    legacyNoise = samplesTree != nullptr;
  }
  if (!samplesTree)
  {
    cerr << "ERROR: expected apv_EventSamples"
         << (preprocessedSignal ? "" : " or legacy noise apv_events")
         << " in " << inputName << endl;
    return false;
  }

  unsigned char apvBoard = 0;
  unsigned int trigger = 0;
  int rawEventIndex = -1;
  vector<unsigned char> *activeChannels = nullptr;
  vector<unsigned char> allChannels;
  if (legacyNoise)
    for (int channelIndex = 0; channelIndex < kInputChannelCount; ++channelIndex)
      allChannels.push_back(static_cast<unsigned char>(channelIndex));
  int samples[kInputChannelCount][kSamplesPerWaveform] = {{0}};

  if (samplesTree->SetBranchAddress("apv_board", &apvBoard) < 0 ||
      samplesTree->SetBranchAddress("trig_num", &trigger) < 0)
  {
    cerr << "ERROR: incompatible APV board/trigger schema in " << inputName << endl;
    return false;
  }
  if (samplesTree->GetBranch("raw_event_index"))
    samplesTree->SetBranchAddress("raw_event_index", &rawEventIndex);
  if (legacyNoise)
    activeChannels = &allChannels;
  else if (samplesTree->SetBranchAddress("apv_chnum", &activeChannels) < 0)
  {
    cerr << "ERROR: decoded signal schema requires apv_chnum in " << inputName << endl;
    return false;
  }
  for (int inputChannel = 0; inputChannel < kInputChannelCount; ++inputChannel)
  {
    const TString branchName = Form("samples_ch%d", inputChannel);
    TLeaf *leaf = samplesTree->GetLeaf(branchName);
    if (!leaf || leaf->GetLenStatic() != kSamplesPerWaveform ||
        samplesTree->SetBranchAddress(branchName, samples[inputChannel]) < 0)
    {
      cerr << "ERROR: expected 30 integer samples in " << branchName
           << " of " << inputName << endl;
      return false;
    }
  }

  ULong64_t inputTimestamp = 0;
  int headerRawEventIndex = -1;
  const bool timestampInSamples = legacyNoise && samplesTree->GetBranch("time_stamp");
  if (timestampInSamples)
    samplesTree->SetBranchAddress("time_stamp", &inputTimestamp);
  if (headersTree)
  {
    headersTree->SetBranchAddress("time_stamp", &inputTimestamp);
    if (headersTree->GetBranch("raw_event_index"))
      headersTree->SetBranchAddress("raw_event_index", &headerRawEventIndex);
  }

  bool firstTrigger = true;
  unsigned int previousRawTrigger = 0;
  unsigned long long triggerEpoch = 0;
  long long fileEventOffset = 0;
  const Long64_t entries = samplesTree->GetEntries();

  for (Long64_t entryIndex = 0; entryIndex < entries; ++entryIndex)
  {
    if (samplesTree->GetEntry(entryIndex) <= 0)
    {
      cerr << "ERROR: cannot read APV frame " << entryIndex << " in " << inputName << endl;
      return false;
    }
    ++inputFrames;

    if (!firstTrigger && previousRawTrigger > trigger &&
        previousRawTrigger - trigger > 32768u)
      triggerEpoch += 65536ull;
    const unsigned long long localEvent = triggerEpoch + trigger;
    if (firstTrigger)
    {
      if (lastOutputEvent >= 0 &&
          localEvent <= static_cast<unsigned long long>(lastOutputEvent))
        fileEventOffset = lastOutputEvent + 1 -
                          static_cast<long long>(localEvent);
      firstTrigger = false;
    }
    previousRawTrigger = trigger;
    const long long extendedEvent =
        static_cast<long long>(localEvent) + fileEventOffset;
    if (extendedEvent > std::numeric_limits<int>::max())
    {
      cerr << "ERROR: unwrapped trigger exceeds TrackRaw event range in "
           << inputName << endl;
      return false;
    }

    timestamp = timestampInSamples ? static_cast<long long>(inputTimestamp) : -1;
    if (headersTree && entryIndex < headersTree->GetEntries())
    {
      headersTree->GetEntry(entryIndex);
      if (headerRawEventIndex < 0 || rawEventIndex < 0 ||
          headerRawEventIndex == rawEventIndex)
        timestamp = static_cast<long long>(inputTimestamp);
      else
        cerr << "WARNING: header/sample raw_event_index mismatch at entry "
             << entryIndex << " in " << inputName << endl;
    }

    sourceFile = fileIndex;
    sourceEntry = static_cast<int>(entryIndex);
    sourceRawEventIndex = rawEventIndex;
    sourceAPVBoard = static_cast<int>(apvBoard);
    sourceTrigger = trigger;
    event = static_cast<int>(extendedEvent);
    readoutBoard = 1;

    const int logicalAPV = LogicalAPV(sourceAPVBoard);
    if (logicalAPV >= 0 && logicalAPV < MMdetN * apvChipsPerTracker)
    {
      board = logicalAPV / apvChipsPerTracker;
      chip = logicalAPV % apvChipsPerTracker;
      ++trackerFrames;
    }
    else if (logicalAPV == 12 || logicalAPV == 13)
    {
      board = -1;
      chip = logicalAPV;
      ++padFrames;
    }
    else
    {
      ++ignoredFrames;
      lastOutputEvent = std::max(lastOutputEvent, extendedEvent);
      continue;
    }

    if (!activeChannels)
    {
      ++invalidChannels;
      lastOutputEvent = std::max(lastOutputEvent, extendedEvent);
      continue;
    }

    bool seen[kInputChannelCount] = {false};
    for (unsigned char encodedChannel : *activeChannels)
    {
      const int inputChannel = static_cast<int>(encodedChannel);
      if (inputChannel < 0 || inputChannel >= kInputChannelCount ||
          seen[inputChannel])
      {
        ++invalidChannels;
        continue;
      }
      seen[inputChannel] = true;
      channel = inputChannel;
      wave.clear();
      wave.reserve(kSamplesPerWaveform);
      for (int sample = 0; sample < kSamplesPerWaveform; ++sample)
        wave.push_back(ToShort(samples[inputChannel][sample]));
      outputTree->Fill();
      ++outputWaveforms;
    }
    lastOutputEvent = std::max(lastOutputEvent, extendedEvent);
  }

  cout << " --> Adapted " << entries << " APV frames from " << inputName
       << endl;
  return true;
}

bool DecodedAPVRoot2TrackRaw::Loop()
{
  if (!outputFile || !outputTree)
    return false;

  long long lastOutputEvent = -1;
  for (size_t inputIndex = 0; inputIndex < inputs.size(); ++inputIndex)
    if (!ConvertFile(inputs[inputIndex], static_cast<int>(inputIndex),
                     lastOutputEvent))
      return false;

  outputFile->cd();
  TNamed preprocessing(
      "InputPreprocessing",
      preprocessedSignal
          ? "decoded APV25 signal ROOT; baseline already subtracted; input channels "
            "already selected at 3 sigma; no pedestal/noise threshold applied by adapter"
          : "decoded APV25 noise ROOT; original sample values retained; no baseline "
            "subtraction or threshold applied by adapter; legacy apv_events uses all 128 channels");
  preprocessing.Write();

  TTree summary("DecodedInputSummary", "Decoded ROOT adapter counters");
  summary.Branch("input_frames", &inputFrames, "input_frames/L");
  summary.Branch("tracker_frames", &trackerFrames, "tracker_frames/L");
  summary.Branch("pad_frames", &padFrames, "pad_frames/L");
  summary.Branch("ignored_frames", &ignoredFrames, "ignored_frames/L");
  summary.Branch("invalid_channels", &invalidChannels,
                 "invalid_channels/L");
  summary.Branch("output_waveforms", &outputWaveforms,
                 "output_waveforms/L");
  summary.Fill();
  summary.Write();
  outputTree->Write();
  outputFile->Close();

  cout << " --> Decoded ROOT adapter summary: frames=" << inputFrames
       << ", Tracker=" << trackerFrames << ", pad=" << padFrames
       << ", ignored=" << ignoredFrames
       << ", waveforms=" << outputWaveforms
       << ", invalid/duplicate channels=" << invalidChannels << endl;
  return true;
}

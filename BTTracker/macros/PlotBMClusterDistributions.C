#include <TCanvas.h>
#include <TColor.h>
#include <TDatime.h>
#include <TError.h>
#include <TFile.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TUUID.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

constexpr double kPredXMinMm = -10.0;
constexpr double kPredXMaxMm = 40.0;
constexpr double kPredYMinMm = -40.0;
constexpr double kPredYMaxMm = 10.0;

void DrawPeakCountStyle(TH1D &histogram, const TString &title,
                        const TString &subtitle, const TString &leftNote,
                        const TString &rightNote, const TString &xTitle,
                        const TString &yTitle, const TString &footnote,
                        const TString &outputStem)
{
  TCanvas canvas(Form("canvas_%s", histogram.GetName()), title, 1000, 750);
  canvas.SetLogy();
  canvas.SetLeftMargin(.12);
  canvas.SetRightMargin(.04);
  canvas.SetBottomMargin(.14);
  canvas.SetTopMargin(.23);

  histogram.SetTitle(Form(";%s;%s", xTitle.Data(), yTitle.Data()));
  histogram.SetStats(false);
  histogram.SetLineColor(kBlue + 1);
  histogram.SetLineWidth(2);
  histogram.SetFillColor(kAzure - 9);
  histogram.SetMinimum(.5);
  histogram.SetMaximum(std::max(1.0, histogram.GetMaximum()) * 8.0);
  histogram.Draw("HIST");

  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(.028);
  text.DrawLatex(.12, .95, title);
  text.DrawLatex(.12, .90, subtitle);
  text.DrawLatex(.12, .85, leftNote);
  text.DrawLatex(.56, .85, rightNote);
  text.SetTextSize(.021);
  text.DrawLatex(.12, .025, footnote);

  canvas.SaveAs(outputStem + ".png");
  canvas.SaveAs(outputStem + ".pdf");
}

} // namespace

// Both plots use selected tracks with X and Y signals in all three trackers,
// and with the BM prediction inside the fiducial area. Cluster number counts
// clusters in the predicted 3x3 pad window (including zero); cluster size
// counts selected events once per cluster-size bin, even if an event has
// multiple clusters of the same size.
void PlotBMClusterDistributions(const char *runDir,
                                const char *detector = "BM2")
{
  const TString detectorName(detector);
  if (detectorName != "BM1" && detectorName != "BM2")
  {
    Error("PlotBMClusterDistributions", "detector must be BM1 or BM2");
    return;
  }
  const int index = detectorName == "BM1" ? 0 : 1;
  const TString inputName = Form("%s/Combine/Pad-step3-rec.root", runDir);
  TFile input(inputName, "READ");
  if (input.IsZombie())
  {
    Error("PlotBMClusterDistributions", "Cannot open %s", inputName.Data());
    return;
  }

  TTree *trackTree = nullptr;
  TTree *clusterTree = nullptr;
  input.GetObject("TrackArrivalEvents", trackTree);
  input.GetObject(detectorName + "Events", clusterTree);
  if (!clusterTree)
  {
    // Older run outputs stored both BM planes in one PadEvents tree. Keep
    // those files usable, while making the combined-plane scope explicit.
    input.GetObject("PadEvents", clusterTree);
  }
  if (!trackTree || !clusterTree)
  {
    Error("PlotBMClusterDistributions",
          "TrackArrivalEvents and %sEvents (or legacy PadEvents) must exist in %s",
          detectorName.Data(), inputName.Data());
    return;
  }

  const TString trackerInputName =
      Form("%s/Combine/Tracker-step4-rec.root", runDir);
  TFile trackerInput(trackerInputName, "READ");
  if (trackerInput.IsZombie())
  {
    Error("PlotBMClusterDistributions", "Cannot open %s",
          trackerInputName.Data());
    return;
  }
  TTree *trackerTree = nullptr;
  trackerInput.GetObject("Track", trackerTree);
  if (!trackerTree || !trackerTree->GetBranch("event") ||
      !trackerTree->GetBranch("useflag") || !trackerTree->GetBranch("MMsigx") ||
      !trackerTree->GetBranch("MMsigy"))
  {
    Error("PlotBMClusterDistributions",
          "Track tree with event/useflag/MMsigx/MMsigy is required in %s",
          trackerInputName.Data());
    return;
  }
  int trackerEvent = -1;
  bool useTrack = false;
  bool trackerHasX[3] = {false, false, false};
  bool trackerHasY[3] = {false, false, false};
  trackerTree->SetBranchAddress("event", &trackerEvent);
  trackerTree->SetBranchAddress("useflag", &useTrack);
  trackerTree->SetBranchAddress("MMsigx", trackerHasX);
  trackerTree->SetBranchAddress("MMsigy", trackerHasY);
  std::set<int> allTrackerSignalEvents;
  for (Long64_t entry = 0; entry < trackerTree->GetEntries(); ++entry)
  {
    trackerTree->GetEntry(entry);
    if (!useTrack)
      continue;
    bool allHaveXY = true;
    for (int layer = 0; layer < 3; ++layer)
      allHaveXY = allHaveXY && trackerHasX[layer] && trackerHasY[layer];
    if (allHaveXY)
      allTrackerSignalEvents.insert(trackerEvent);
  }
  trackerTree->ResetBranchAddresses();
  if (allTrackerSignalEvents.empty())
  {
    Error("PlotBMClusterDistributions",
          "No tracks have X and Y signals in all three tracker layers");
    return;
  }

  int matchedClusterCount[2] = {0, 0};
  double predictedX[2] = {-999.0, -999.0};
  double predictedY[2] = {-999.0, -999.0};
  for (const char *branch : {"bm_pred_x", "bm_pred_y"})
  {
    if (!trackTree->GetBranch(branch))
    {
      Error("PlotBMClusterDistributions",
            "TrackArrivalEvents is missing branch %s", branch);
      return;
    }
  }
  const bool hasMatchedCounts =
      trackTree->GetBranch("bm_matched_cluster_count") != nullptr;
  int trackEvent = -1;
  if (!trackTree->GetBranch("event"))
  {
    Error("PlotBMClusterDistributions", "TrackArrivalEvents is missing event");
    return;
  }
  trackTree->SetBranchAddress("event", &trackEvent);
  if (hasMatchedCounts)
    trackTree->SetBranchAddress("bm_matched_cluster_count", matchedClusterCount);
  trackTree->SetBranchAddress("bm_pred_x", predictedX);
  trackTree->SetBranchAddress("bm_pred_y", predictedY);

  std::map<int, int> legacyEventClusterCounts;
  if (!hasMatchedCounts)
  {
    int padEvent = -1;
    int eventClusterCount = 0;
    if (!clusterTree->GetBranch("event") || !clusterTree->GetBranch("nclusters"))
    {
      Error("PlotBMClusterDistributions",
            "Legacy PadEvents is missing event or nclusters");
      return;
    }
    clusterTree->SetBranchAddress("event", &padEvent);
    clusterTree->SetBranchAddress("nclusters", &eventClusterCount);
    for (Long64_t entry = 0; entry < clusterTree->GetEntries(); ++entry)
    {
      clusterTree->GetEntry(entry);
      legacyEventClusterCounts[padEvent] = eventClusterCount;
    }
    clusterTree->ResetBranchAddresses();
  }

  gROOT->cd();
  TH1D numberHistogram("bm_cluster_number", "", 10, -.5, 9.5);
  numberHistogram.SetDirectory(nullptr);
  for (int bin = 1; bin <= 9; ++bin)
    numberHistogram.GetXaxis()->SetBinLabel(bin, Form("%d", bin - 1));
  numberHistogram.GetXaxis()->SetBinLabel(10, "#geq9");
  Long64_t acceptedTracks = 0;
  Long64_t tracksWithZeroClusters = 0;
  Long64_t tracksWithTwoOrMore = 0;
  std::set<int> acceptedEventIds;
  for (Long64_t entry = 0; entry < trackTree->GetEntries(); ++entry)
  {
    trackTree->GetEntry(entry);
    if (allTrackerSignalEvents.find(trackEvent) == allTrackerSignalEvents.end())
      continue;
    if (!std::isfinite(predictedX[index]) || !std::isfinite(predictedY[index]) ||
        predictedX[index] < kPredXMinMm || predictedX[index] > kPredXMaxMm ||
        predictedY[index] < kPredYMinMm || predictedY[index] > kPredYMaxMm)
      continue;
    int count = 0;
    if (hasMatchedCounts)
      count = matchedClusterCount[index];
    else
    {
      const auto found = legacyEventClusterCounts.find(trackEvent);
      if (found != legacyEventClusterCounts.end())
        count = found->second;
    }
    if (count < 0)
      continue;
    numberHistogram.Fill(std::min(count, 9));
    ++acceptedTracks;
    acceptedEventIds.insert(trackEvent);
    if (count == 0)
      ++tracksWithZeroClusters;
    if (count >= 2)
      ++tracksWithTwoOrMore;
  }
  trackTree->ResetBranchAddresses();
  if (acceptedTracks == 0)
  {
    Error("PlotBMClusterDistributions",
          "No selected tracks inside the %s fiducial area", detectorName.Data());
    return;
  }

  std::vector<int> *clusterSizes = nullptr;
  if (!clusterTree->GetBranch("cluster_size"))
  {
    Error("PlotBMClusterDistributions", "%s is missing cluster_size",
          clusterTree->GetName());
    return;
  }
  clusterTree->SetBranchAddress("cluster_size", &clusterSizes);
  int clusterEvent = -1;
  if (!clusterTree->GetBranch("event"))
  {
    Error("PlotBMClusterDistributions", "%s is missing event",
          clusterTree->GetName());
    return;
  }
  clusterTree->SetBranchAddress("event", &clusterEvent);
  TH1D sizeHistogram("bm_cluster_size", "", 20, .5, 20.5);
  sizeHistogram.SetDirectory(nullptr);
  for (int bin = 1; bin <= 19; ++bin)
    sizeHistogram.GetXaxis()->SetBinLabel(bin, Form("%d", bin));
  sizeHistogram.GetXaxis()->SetBinLabel(20, "#geq20");
  Long64_t totalClusters = 0;
  std::map<int, std::set<int>> sizesByEvent;
  for (Long64_t entry = 0; entry < clusterTree->GetEntries(); ++entry)
  {
    clusterTree->GetEntry(entry);
    if (acceptedEventIds.find(clusterEvent) == acceptedEventIds.end())
      continue;
    if (!clusterSizes)
      continue;
    for (const int size : *clusterSizes)
    {
      if (size < 1)
        continue;
      ++totalClusters;
      sizesByEvent[clusterEvent].insert(std::min(size, 20));
    }
  }
  clusterTree->ResetBranchAddresses();
  if (totalClusters == 0)
  {
    Error("PlotBMClusterDistributions", "%s contains no reconstructed clusters",
          detectorName.Data());
    return;
  }

  Long64_t eventsWithClusters = 0;
  Long64_t eventsWithSizeTwoOrMore = 0;
  for (const auto &eventEntry : sizesByEvent)
  {
    const auto &sizes = eventEntry.second;
    if (sizes.empty())
      continue;
    ++eventsWithClusters;
    bool hasLargeCluster = false;
    for (const int size : sizes)
    {
      sizeHistogram.Fill(size);
      hasLargeCluster = hasLargeCluster || size >= 2;
    }
    if (hasLargeCluster)
      ++eventsWithSizeTwoOrMore;
  }

  const double trackFraction = 100.0 * tracksWithTwoOrMore / acceptedTracks;
  const double zeroTrackFraction = 100.0 * tracksWithZeroClusters / acceptedTracks;
  const double sizeEventFraction =
      100.0 * eventsWithSizeTwoOrMore / acceptedTracks;
  TDatime now;
  TUUID uuid;
  const TString uuidString(uuid.AsString());
  const TString shortUuid = uuidString(0, 8);
  const TString outputDir = Form("%s/Combine/cluster_diagnostics/%04d%02d%02dT%02d%02d%02d_%s",
                                 runDir, now.GetYear(), now.GetMonth(), now.GetDay(),
                                 now.GetHour(), now.GetMinute(), now.GetSecond(),
                                 shortUuid.Data());
  if (gSystem->mkdir(outputDir, true) != 0)
  {
    Error("PlotBMClusterDistributions", "Cannot create %s", outputDir.Data());
    return;
  }

  const TString numberTitle = hasMatchedCounts
      ? Form("%s: matched clusters per track", detectorName.Data())
      : "BM1 + BM2: clusters per selected track";
  const TString clusterSubtitle = hasMatchedCounts
      ? "Selected events: all Tracker X/Y; BM2 x[-10,40], y[-40,10] mm"
      : "Legacy BM1 + BM2 clusters; same Tracker and prediction selection";
  DrawPeakCountStyle(
      numberHistogram, numberTitle,
      hasMatchedCounts
          ? Form("All three Trackers have X/Y; %s prediction x[-10,40], y[-40,10] mm",
                 detectorName.Data())
          : Form("All three Trackers have X/Y; %s prediction x[-10,40], y[-40,10] mm",
                 detectorName.Data()),
      Form("Tracks = %lld", acceptedTracks),
      Form("N_{cluster}=0: %lld (%.2f%%);  N_{cluster} #geq 2: %lld (%.2f%%)",
           tracksWithZeroClusters, zeroTrackFraction, tracksWithTwoOrMore,
           trackFraction),
      "BM clusters per selected track", "Tracks",
      hasMatchedCounts
          ? "Clusters in predicted 3x3 pad window; zero-count tracks included"
          : "Legacy BM1 + BM2 count; zero-count tracks included",
      outputDir + Form("/%s_cluster_number", detectorName.Data()));
  DrawPeakCountStyle(
      sizeHistogram, hasMatchedCounts
          ? Form("%s: reconstructed cluster size", detectorName.Data())
          : "BM1 + BM2: reconstructed cluster size",
      clusterSubtitle,
      Form("Events with clusters: %lld / %lld", eventsWithClusters, acceptedTracks),
      Form("Size #geq2 pads: %lld events (%.2f%%)",
           eventsWithSizeTwoOrMore, sizeEventFraction),
      "Cluster size [pads]", "Events",
      "Each event counted once per size bin; final bin includes sizes #geq20",
      outputDir + Form("/%s_cluster_size", detectorName.Data()));

  const TString rootName = outputDir + "/bm_cluster_distributions.root";
  TFile output(rootName, "RECREATE");
  numberHistogram.Write();
  sizeHistogram.Write();
  output.Close();
  std::cout << "Tracks in fiducial area: " << acceptedTracks
            << "; tracks with >=2 matched clusters: " << tracksWithTwoOrMore
            << " (" << trackFraction << "%)" << std::endl;
  std::cout << "Reconstructed clusters: " << totalClusters
            << "; events with a cluster >=2 pads: " << eventsWithSizeTwoOrMore
            << " (" << sizeEventFraction << "% of selected tracks)" << std::endl;
  std::cout << "Output: " << outputDir << std::endl;
}

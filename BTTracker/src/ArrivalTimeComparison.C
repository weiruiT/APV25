#include "CRAnalysis.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TTree.h"

namespace
{
struct Hit {
  int strip = -1;
  double arrival = -999.0;
  double charge = 0.0;
};

struct Cluster {
  std::vector<Hit> hits;
  double earliest = std::numeric_limits<double>::infinity();
  double mean = 0.0;
  double charge = 0.0;
};

struct Moments {
  long long n = 0;
  double mean = 0.0;
  double m2 = 0.0;

  void Add(double value) {
    ++n;
    const double delta = value - mean;
    mean += delta / static_cast<double>(n);
    m2 += delta * (value - mean);
  }

  double Variance() const {
    return n > 1 ? m2 / static_cast<double>(n - 1) : 0.0;
  }

  double StdDev() const { return std::sqrt(Variance()); }
};

std::vector<Cluster> BuildClusters(std::vector<Hit> hits, int maxStripGap) {
  std::sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
    return a.strip < b.strip;
  });

  std::vector<Cluster> clusters;
  for (const Hit &hit : hits) {
    if (clusters.empty() ||
        hit.strip - clusters.back().hits.back().strip > maxStripGap) {
      clusters.emplace_back();
    }
    clusters.back().hits.push_back(hit);
  }

  for (Cluster &cluster : clusters) {
    double sumTime = 0.0;
    for (const Hit &hit : cluster.hits) {
      cluster.earliest = std::min(cluster.earliest, hit.arrival);
      sumTime += hit.arrival;
      cluster.charge += hit.charge;
    }
    cluster.mean = sumTime / static_cast<double>(cluster.hits.size());
  }
  return clusters;
}

const Cluster *FindFirstArrivalCluster(const std::vector<Cluster> &clusters) {
  if (clusters.empty()) {
    return nullptr;
  }
  return &*std::min_element(
      clusters.begin(), clusters.end(),
      [](const Cluster &a, const Cluster &b) {
        return a.earliest < b.earliest;
      });
}
} // namespace

void CompareArrivalTimes(TString fileDir) {
  const TString inputName = fileDir + "/Combine/Tracker-step2-ana.root";
  const TString outputName =
      fileDir + "/Combine/Tracker-arrival-time-comparison.root";
  const TString tableName =
      fileDir + "/Combine/Tracker-arrival-time-comparison.tsv";
  const TString overlayName =
      fileDir + "/Combine/Tracker-arrival-time-comparison.png";
  const TString differenceName =
      fileDir + "/Combine/Tracker-arrival-time-difference.png";

  TFile input(inputName, "READ");
  if (input.IsZombie()) {
    std::cerr << "Cannot open analysis file: " << inputName << std::endl;
    return;
  }
  TFile output(outputName, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Cannot create output file: " << outputName << std::endl;
    return;
  }

  constexpr int maxStripGap = 2;
  TH1D *hEarliest[2 * MMdetN] = {};
  TH1D *hMean[2 * MMdetN] = {};
  TH1D *hDifference[2 * MMdetN] = {};

  TTree comparison("arrival_time_comparison",
                   "Event arrival time comparison");
  int tracker = -1;
  int axis = -1;
  int event = -1;
  int nClusters = 0;
  int firstClusterSize = 0;
  double firstClusterCharge = 0.0;
  double earliestArrival = -999.0;
  double clusterMeanArrival = -999.0;
  double difference = -999.0;

  comparison.Branch("tracker", &tracker, "tracker/I");
  comparison.Branch("axis", &axis, "axis/I");
  comparison.Branch("event", &event, "event/I");
  comparison.Branch("nClusters", &nClusters, "nClusters/I");
  comparison.Branch("firstClusterSize", &firstClusterSize, "firstClusterSize/I");
  comparison.Branch("firstClusterCharge", &firstClusterCharge,
                    "firstClusterCharge/D");
  comparison.Branch("earliestArrival", &earliestArrival,
                    "earliestArrival/D");
  comparison.Branch("clusterMeanArrival", &clusterMeanArrival,
                    "clusterMeanArrival/D");
  comparison.Branch("difference", &difference, "difference/D");

  std::ofstream table(tableName.Data());
  table << "tracker\taxis\tevents\tmean_earliest\tvariance_earliest"
        << "\tstd_earliest\tmean_cluster_mean\tvariance_cluster_mean"
        << "\tstd_cluster_mean\tmean_difference\tvariance_difference"
        << "\tstd_difference\n";

  Moments allEarliest;
  Moments allMean;
  Moments allDifference;

  std::cout << "Arrival-time comparison (sample variance, arrivalT in APV samples)\n";
  std::cout << "Cluster rule: strip gap <= " << maxStripGap
            << "; first cluster = cluster with the earliest waveform arrival\n";

  for (int det = 0; det < MMdetN; ++det) {
    for (int xy = 0; xy < 2; ++xy) {
      const int index = det * 2 + xy;
      const TString axisName = xy == 0 ? "X" : "Y";
      const TString treeName = Form("Track%s%d", axisName.Data(), det);

      TTree *tree = nullptr;
      input.GetObject(treeName, tree);
      if (!tree) {
        std::cerr << "Missing tree: " << treeName << std::endl;
        continue;
      }

      int treeEvent = -1;
      int channel = -1;
      double maxA = 0.0;
      double arrivalT = -999.0;
      int sigmoidFitStatus = -1;
      tree->SetBranchAddress("event", &treeEvent);
      tree->SetBranchAddress("channel", &channel);
      tree->SetBranchAddress("maxA", &maxA);
      tree->SetBranchAddress("arrivalT_sigmoid", &arrivalT);
      tree->SetBranchAddress("sigmoidFitStatus", &sigmoidFitStatus);

      std::map<int, std::vector<Hit>> eventHits;
      const Long64_t entries = tree->GetEntries();
      for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        if (sigmoidFitStatus != 0 || arrivalT <= -900.0 ||
            !std::isfinite(arrivalT)) {
          continue;
        }
        if (channel < 0 || channel >= 256 || treeEvent < 0) {
          continue;
        }
        eventHits[treeEvent].push_back({channel, arrivalT, maxA});
      }

      hEarliest[index] =
          new TH1D(Form("hArrivalEarliest_T%d%s", det + 1, axisName.Data()),
                   Form("Tracker%d %s;earliest waveform arrivalT;Events",
                        det + 1, axisName.Data()),
                   120, -5.0, 35.0);
      hMean[index] =
          new TH1D(Form("hArrivalMean_T%d%s", det + 1, axisName.Data()),
                   Form("Tracker%d %s;first-cluster mean arrivalT;Events",
                        det + 1, axisName.Data()),
                   120, -5.0, 35.0);
      hDifference[index] = new TH1D(
          Form("hArrivalDifference_T%d%s", det + 1, axisName.Data()),
          Form("Tracker%d %s;earliest - first-cluster mean;Events", det + 1,
               axisName.Data()),
          120, -15.0, 15.0);

      Moments earliestMoments;
      Moments meanMoments;
      Moments differenceMoments;

      for (const auto &eventEntry : eventHits) {
        const std::vector<Cluster> clusters =
            BuildClusters(eventEntry.second, maxStripGap);
        const Cluster *first = FindFirstArrivalCluster(clusters);
        if (!first) {
          continue;
        }

        tracker = det + 1;
        axis = xy;
        event = eventEntry.first;
        nClusters = static_cast<int>(clusters.size());
        firstClusterSize = static_cast<int>(first->hits.size());
        firstClusterCharge = first->charge;
        earliestArrival = first->earliest;
        clusterMeanArrival = first->mean;
        difference = earliestArrival - clusterMeanArrival;
        comparison.Fill();

        earliestMoments.Add(earliestArrival);
        meanMoments.Add(clusterMeanArrival);
        differenceMoments.Add(difference);
        allEarliest.Add(earliestArrival);
        allMean.Add(clusterMeanArrival);
        allDifference.Add(difference);
        hEarliest[index]->Fill(earliestArrival);
        hMean[index]->Fill(clusterMeanArrival);
        hDifference[index]->Fill(difference);
      }

      hEarliest[index]->Write();
      hMean[index]->Write();
      hDifference[index]->Write();

      table << det + 1 << '\t' << axisName << '\t' << earliestMoments.n << '\t'
            << std::setprecision(10) << earliestMoments.mean << '\t'
            << earliestMoments.Variance() << '\t'
            << earliestMoments.StdDev() << '\t' << meanMoments.mean << '\t'
            << meanMoments.Variance() << '\t' << meanMoments.StdDev() << '\t'
            << differenceMoments.mean << '\t' << differenceMoments.Variance()
            << '\t' << differenceMoments.StdDev() << '\n';

      std::cout << "Tracker" << det + 1 << ' ' << axisName
                << ": N=" << earliestMoments.n
                << ", earliest variance=" << earliestMoments.Variance()
                << " (std=" << earliestMoments.StdDev()
                << "), cluster-mean variance=" << meanMoments.Variance()
                << " (std=" << meanMoments.StdDev()
                << "), mean(earliest-mean)=" << differenceMoments.mean << '\n';
    }
  }

  comparison.Write();
  table << "ALL\tXY\t" << allEarliest.n << '\t'
        << std::setprecision(10) << allEarliest.mean << '\t'
        << allEarliest.Variance() << '\t' << allEarliest.StdDev() << '\t'
        << allMean.mean << '\t' << allMean.Variance() << '\t'
        << allMean.StdDev() << '\t' << allDifference.mean << '\t'
        << allDifference.Variance() << '\t' << allDifference.StdDev() << '\n';
  table.close();

  TCanvas overlay("arrival_time_comparison", "Arrival time comparison", 1600, 900);
  overlay.Divide(MMdetN, 2);
  for (int index = 0; index < 2 * MMdetN; ++index) {
    overlay.cd(index + 1);
    if (!hEarliest[index]) {
      continue;
    }
    hEarliest[index]->SetLineColor(kBlue + 1);
    hEarliest[index]->SetLineWidth(2);
    hMean[index]->SetLineColor(kRed + 1);
    hMean[index]->SetLineWidth(2);
    hEarliest[index]->Draw();
    hMean[index]->Draw("SAME");
    TLegend legend(0.52, 0.75, 0.95, 0.92);
    legend.AddEntry(hEarliest[index], "earliest waveform", "l");
    legend.AddEntry(hMean[index], "first-cluster mean", "l");
    legend.Draw();
  }
  overlay.SaveAs(overlayName);

  TCanvas differenceCanvas("arrival_time_difference", "Arrival time difference",
                           1600, 900);
  differenceCanvas.Divide(MMdetN, 2);
  for (int index = 0; index < 2 * MMdetN; ++index) {
    differenceCanvas.cd(index + 1);
    if (hDifference[index]) {
      hDifference[index]->Draw();
    }
  }
  differenceCanvas.SaveAs(differenceName);

  output.Write();
  output.Close();
  input.Close();

  std::cout << "All planes: earliest variance=" << allEarliest.Variance()
            << ", cluster-mean variance=" << allMean.Variance()
            << ", mean difference=" << allDifference.mean << '\n';
  std::cout << "Wrote: " << outputName << '\n'
            << "       " << tableName << '\n'
            << "       " << overlayName << '\n'
            << "       " << differenceName << std::endl;
}

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TString.h"

#include <iostream>
#include <vector>

void show_bm_wave(
    const char *fileName = "../../Combine/Tracker-step2-ana.root",
    Long64_t entry = 0)
{
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Cannot open file: " << fileName << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("BMPadSignals");
    if (!tree)
    {
        std::cerr << "BMPadSignals tree not found" << std::endl;
        return;
    }

    int event = -1;
    int globalChannel = -1;
    int rawChip = -1;
    int rawChannel = -1;
    int channel = -1;
    int position = -1;
    int maxT = -1;
    int nisel = 0;
    int pulseFitStatus = -1;

    double ped = 0.0;
    double pedrms = 0.0;
    double maxA = 0.0;
    double arrivalT = -999.0;

    std::vector<short> *wave = nullptr;

    tree->SetBranchAddress("event", &event);
    tree->SetBranchAddress("global_channel", &globalChannel);
    tree->SetBranchAddress("raw_chip", &rawChip);
    tree->SetBranchAddress("raw_channel", &rawChannel);
    tree->SetBranchAddress("channel", &channel);
    tree->SetBranchAddress("position", &position);
    tree->SetBranchAddress("ped", &ped);
    tree->SetBranchAddress("pedrms", &pedrms);
    tree->SetBranchAddress("maxA", &maxA);
    tree->SetBranchAddress("maxT", &maxT);
    tree->SetBranchAddress("nisel", &nisel);
    tree->SetBranchAddress("arrivalT", &arrivalT);
    tree->SetBranchAddress("pulseFitStatus", &pulseFitStatus);
    tree->SetBranchAddress("wave", &wave);

    if (entry < 0 || entry >= tree->GetEntries())
    {
        std::cerr << "Invalid entry: " << entry << std::endl;
        return;
    }

    tree->GetEntry(entry);

    std::cout << "entry = " << entry << std::endl;
    std::cout << "event = " << event << std::endl;
    std::cout << "global_channel = " << globalChannel << std::endl;
    std::cout << "raw_chip = " << rawChip
              << ", raw_channel = " << rawChannel << std::endl;
    std::cout << "position = " << position << std::endl;
    std::cout << "ped = " << ped
              << ", pedrms = " << pedrms << std::endl;
    std::cout << "maxA = " << maxA
              << ", maxT = " << maxT << std::endl;
    std::cout << "arrivalT = " << arrivalT
              << ", pulseFitStatus = " << pulseFitStatus << std::endl;

    if (!wave || wave->empty())
    {
        std::cerr << "Empty waveform" << std::endl;
        return;
    }

    TGraph rawGraph;
    TGraph correctedGraph;

    for (int i = 0; i < (int)wave->size(); ++i)
    {
        rawGraph.SetPoint(i, i, wave->at(i));
        correctedGraph.SetPoint(i, i, wave->at(i) - ped);
    }

    TCanvas canvas("canvasBM", "BM waveform", 900, 700);
    canvas.Divide(1, 2);

    canvas.cd(1);
    rawGraph.SetTitle("BM raw waveform;Sample;ADC");
    rawGraph.SetMarkerStyle(20);
    rawGraph.Draw("ALP");

    canvas.cd(2);
    correctedGraph.SetTitle("BM pedestal-subtracted waveform;Sample;ADC-ped");
    correctedGraph.SetMarkerStyle(20);
    correctedGraph.Draw("ALP");

    TString outputName =
        Form("bm_wave_entry%lld_global%d.png", entry, globalChannel);

    canvas.SaveAs(outputName);
    std::cout << "Saved: " << outputName << std::endl;
}
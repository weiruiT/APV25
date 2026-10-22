#include <TError.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <string>

namespace {
std::string BM2SpatialQuote(const char *value)
{
    std::string result = "'";
    for (const char c : std::string(value))
        result += c == '\'' ? "'\\''" : std::string(1, c);
    return result + "'";
}
}

void PlotBM2SpatialTiming(const char *run, int minEvents = 40,
                          double maxAMax = 2000., double cellSizeMm = 20.)
{
    TString script = gSystem->DirName(__FILE__);
    script += "/plot_bm2_spatial_timing.py";
    gSystem->ExpandPathName(script);
    const char *python = gSystem->Getenv("BTTRACKER_PYTHON");
    if (!python || !*python)
        python = gSystem->AccessPathName("/usr/bin/python3") ? "python3" : "/usr/bin/python3";
    const std::string command = BM2SpatialQuote(python) + " " +
        BM2SpatialQuote(script.Data()) + " " + BM2SpatialQuote(run) +
        " --min-events " + std::to_string(minEvents) +
        " --maxA-max " + std::to_string(maxAMax) +
        " --cell-size-mm " + std::to_string(cellSizeMm);
    const int status = gSystem->Exec(command.c_str());
    if (status != 0) {
        Error("PlotBM2SpatialTiming", "Plotting failed (status %d); inspect the Python error above", status);
        if (gROOT->IsBatch()) gSystem->Exit(1);
    }
}

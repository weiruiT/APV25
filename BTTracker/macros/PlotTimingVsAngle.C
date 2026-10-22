// ROOT entry point. Reuses the Python production timing estimator via PyROOT.
#include <TError.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <string>

namespace {
std::string AngleShellQuote(const char *text)
{
    std::string result = "'";
    for (const char c : std::string(text))
        result += c == '\'' ? "'\\''" : std::string(1, c);
    return result + "'";
}
}

void PlotTimingVsAngle(const char *run, int bins = 5,
                       const char *angleEdges = "", int minEvents = 40,
                       int timingMethod = 0)
{
    TString script = gSystem->DirName(__FILE__);
    script += "/../scripts/plot_timing_vs_angle.py";
    gSystem->ExpandPathName(script);
    const char *python = gSystem->Getenv("BTTRACKER_PYTHON");
    if (!python || !*python)
        python = gSystem->AccessPathName("/usr/bin/python3") ? "python3" : "/usr/bin/python3";
    std::string command = AngleShellQuote(python) + " " + AngleShellQuote(script.Data()) +
        " " + AngleShellQuote(run) + " --bins " + std::to_string(bins) +
        " --min-events " + std::to_string(minEvents) +
        " --timing-method " + std::to_string(timingMethod);
    if (angleEdges && *angleEdges)
        command += " --angle-edges " + AngleShellQuote(angleEdges);
    const int status = gSystem->Exec(command.c_str());
    if (status != 0) {
        Error("PlotTimingVsAngle", "Plotting failed (status %d); inspect the Python error above", status);
        if (gROOT->IsBatch()) gSystem->Exit(1);
    }
}

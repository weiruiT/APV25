#include "CRAnalysis.h"
#include "CRConfig.h"

#include <cstdlib>
#include <iostream>
#include <string>


namespace
{
bool Enabled(const CRConfig &config, const std::string &key)
{
  return config.GetBool(key, true);
}

bool NoiseAnalysisEnabled(const CRConfig &config)
{
  return config.GetBool("noise_analysis", true);
}
}

int RunCommand(int argc, char **argv)
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0]
              << " <data-directory-or-root-file> [force] [mode]" << std::endl;
    return 1;
  }

  const TString fileDir = argv[1];
  if (argc > 2 && std::string(argv[2]) != "0" && std::string(argv[2]) != "1")
  {
    std::cerr << "ERROR: force must be 0 or 1." << std::endl;
    return 2;
  }
  const int force = argc > 2 ? std::atoi(argv[2]) : 0;
  const std::string mode = argc > 3 ? argv[3] : "signal";

  CRConfig config;
  if (!config.Load())
    return 2;

  if (mode == "noise")
  {
    if (NoiseAnalysisEnabled(config))
      return CalibrateTrackerNoise(fileDir, force) ? 0 : 1;
    else
    {
      std::cerr << "ERROR: requested noise mode is disabled in configCR.json." << std::endl;
      return 2;
    }
  }
  else if (mode == "amplitude")
  {
    if (Enabled(config, "amplitude_distribution"))
      MakeAmplitudeDistributions(fileDir);
    else
    {
      std::cerr << "ERROR: requested amplitude distribution is disabled in configCR.json." << std::endl;
      return 2;
    }
  }
  else if (mode == "energy")
  {
    if (Enabled(config, "energy_spectrum"))
      MakeEnergySpectra(fileDir);
    else
    {
      std::cerr << "ERROR: requested energy spectrum is disabled in configCR.json." << std::endl;
      return 2;
    }
  }
  else if (mode == "arrival")
  {
    if (Enabled(config, "arrival_time_comparison"))
      CompareArrivalTimes(fileDir);
    else
    {
      std::cerr << "ERROR: requested arrival-time comparison is disabled in configCR.json." << std::endl;
      return 2;
    }
  }
  else if (mode == "signal")
  {
    return ReadTrackerConfigured(fileDir, force) ? 0 : 1;
  }
  else
  {
    std::cerr << "ERROR: unknown mode '" << mode
              << "'. Expected signal, noise, amplitude, energy or arrival."
              << std::endl;
    return 1;
  }

  return 0;
}

int main(int argc, char **argv)
{
  try
  {
    return RunCommand(argc, argv);
  }
  catch (const std::exception &error)
  {
    std::cerr << "ERROR: " << error.what() << std::endl;
    return 1;
  }
}

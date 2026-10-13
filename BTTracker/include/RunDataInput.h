#ifndef RUNDATAINPUT_H
#define RUNDATAINPUT_H 1

#include "CRInfo.h"

#include <algorithm>
#include <filesystem>

// Only enumerate the requested acquisition directory. In particular, never
// recursively combine different runs or mix Tracker/physics and Tracker/noise.
inline bool CollectRunInputs(const TString &directory,
                             vector<TString> &roots)
{
  namespace fs = std::filesystem;
  std::error_code error;
  if (!fs::is_directory(directory.Data(), error))
  {
    cerr << "ERROR: input directory does not exist: " << directory << endl;
    return false;
  }
  for (fs::directory_iterator it(directory.Data(), error), end;
       it != end && !error; it.increment(error))
  {
    if (!it->is_regular_file(error))
      continue;
    const fs::path path = it->path();
    if (path.extension() == ".root")
      roots.emplace_back(path.string());
  }
  if (error)
  {
    cerr << "ERROR: cannot enumerate " << directory << ": "
         << error.message() << endl;
    return false;
  }
  std::sort(roots.begin(), roots.end());
  if (roots.empty())
  {
    cerr << "ERROR: no decoded ROOT input in " << directory << endl;
    return false;
  }
  return true;
}

inline bool ValidateRootTrees(const TString &path,
                              const vector<TString> &treeNames,
                              bool requireEntries = false)
{
  TFile input(path, "READ");
  if (input.IsZombie() || input.TestBit(TFile::kRecovered))
  {
    cerr << "ERROR: missing, invalid or incomplete ROOT file: " << path << endl;
    return false;
  }
  for (const TString &name : treeNames)
  {
    TTree *tree = dynamic_cast<TTree *>(input.Get(name));
    if (!tree || (requireEntries && tree->GetEntries() == 0))
    {
      cerr << "ERROR: missing or empty required tree " << name
           << " in " << path << endl;
      return false;
    }
  }
  return true;
}

inline bool EnsureOutputDirectory(const TString &directory)
{
  std::error_code error;
  std::filesystem::create_directories(directory.Data(), error);
  if (error)
    cerr << "ERROR: cannot create " << directory << ": "
         << error.message() << endl;
  return !error;
}

#endif

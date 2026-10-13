#include "CRInfo.h"

using namespace std;

TString GetPath(TString fileName)
{
  // cout << " --> The fileName is:" << fileName << endl;
  TString prename;
  prename.Replace(0, 1000, fileName, fileName.Last('/'));
  // cout << " --> Get Path:" << prename << endl;
  return prename;
}

TString GetFilename(TString fileName)
{
  // cout << " --> The fileName is:" << fileName << endl;
  TString prename;
  prename = fileName.Copy().Remove(0, fileName.Last('/') + 1);
  prename = prename.Remove(prename.Last('.'), prename.Length());
  // cout << " --> Get Filename:" << prename << endl;
  return prename;
}

void GetFileList(TString filePath, TString filePattern, vector<TString> &fList)
{
  if (filePath.EndsWith("/"))
    filePath.Remove(filePath.Length() - 1, 1);

  FILE *fp = gSystem->OpenPipe("ls " + filePath + "/*" + filePattern, "r");
  if (!fp)
  {
    cout << " ----> NO " << filePattern << " data files exists in "
         << filePath << "!" << endl;
    return;
  }

  char line[1000];
  fList.clear();
  while (fgets(line, sizeof(line), fp))
  {
    TString s(line);
    if (s.Index(filePattern) == -1)
      continue;
    fList.push_back(s.ReplaceAll("\n", ""));
  }
}

ofstream LogFile(TString fileName)
{
  TString fname = GetFilename(fileName);
  TString fpath = GetPath(fileName);
  ofstream op(Form("%s/%s-log.dat", fpath.Data(), fname.Data()));
  return op;
}

void DrawProcessbar(int i, int total)
{
  static int j = 0;
  /// return ;
  if (i == 0)
    j = 0;
  if (total && (j + 1 == i * 100 / total || !j || i == total - 1))
  {
    if (i == total - 1)
      j++;
    cout << "\r --> ==================== " << j << "% =====================";
    if (i != total - 1)
      j++;
    fflush(stdout);
    if (j == 100)
    {
      j = 0;
      cout << endl;
    }
  }
}

int IsMatched(int trigid, vector<int> idvec)
{
  if (idvec.size() < 1)
    return false;
  int low, high, mid;
  low = 0;
  high = idvec.size() - 1;
  while (low <= high)
  {
    mid = (low + high) / 2;
    if (idvec[mid] > trigid)
      high = mid - 1;
    if (idvec[mid] < trigid)
      low = mid + 1;
    if (idvec[mid] == trigid)
      return mid + 1;
  }
  return 0;
}

int IsMatched(double RecT, vector<double> RecTs)
{
  if (RecTs.size() < 1)
    return false;
  int low, high, mid;
  low = 0;
  high = RecTs.size() - 1;
  while (low <= high)
  {
    mid = (low + high) / 2;
    if (RecTs[mid] > RecT)
      high = mid - 1;
    if (RecTs[mid] < RecT)
      low = mid + 1;
    if (RecTs[mid] == RecT)
      return mid + 1;
  }
  return 0;
}

int IsMatched(int RecT, vector<int> RecTs, int index)
{
  if (RecTs.size() < 1)
    return false;

  if (RecTs[index] == RecT)
  {
    return index+1;
  }
  int low, high, mid;
  low = 0;
  high = RecTs.size() - 1;
  while (low <= high)
  {
    mid = (low + high) / 2;
    if (RecTs[mid] > RecT)
      high = mid - 1;
    if (RecTs[mid] < RecT)
      low = mid + 1;
    if (RecTs[mid] == RecT)
      return mid + 1;
  }
  return 0;
}


TString XYname(int findex)
{
  TString xy[] = {"X", "Y", "X&Y", "XorY"};
  return xy[findex];
}

TLatex *CreatLatex(TString text, double x, double y, int textFont)
{
  TLatex *latex = new TLatex(x, y, text.Data());
  latex->SetNDC();
  latex->SetTextFont(textFont);
  latex->SetTextSize(0.07);
  latex->SetTextColor(1);
  return latex;
}

double Np(double WaveLength) // unit of WaveLength: nm
{
  WaveLength = WaveLength / 1000.; // nm to um
  double A1 = 0.6961663, B1 = 0.0684043;
  double A2 = 0.4079426, B2 = 0.1162414;
  double A3 = 0.8974794, B3 = 9.8961610;
  return TMath::Sqrt(1. + A1 * WaveLength * WaveLength / (WaveLength * WaveLength - B1 * B1) + A2 * WaveLength * WaveLength / (WaveLength * WaveLength - B2 * B2) + A3 * WaveLength * WaveLength / (WaveLength * WaveLength - B3 * B3));
}

double Ng(double WaveLength) // unit of WaveLength: nm
{
  double np = Np(WaveLength);
  double Dnp = (Np(WaveLength + 0.001) - Np(WaveLength - 0.001)) / 0.002;
  return np / (1 + WaveLength / np * Dnp);
}

Double_t CalVectorMean(vector<double> Time, int nn , int bb )
{
  int N = Time.size();
  double TMean = 0;
  if (N == 0)
    return 0;
  for (int iN = 0 + bb; iN < N - nn; iN++)
  {
    TMean += Time.at(iN);
  }
  TMean /= (N - nn - bb);
  return TMean;
}
void SortRecT(vector<double> &TOAevent)
{

  for (int i = 0; i < TOAevent.size() - 1; i++)
  {
    for (int j = 0; j < TOAevent.size() - i - 1; j++)
    {
      if (TOAevent.at(j) > TOAevent.at(j + 1))
      {
        swap(TOAevent.at(j), TOAevent.at(j + 1));
        // swap(channelevent.at(j), channelevent.at(j + 1));
        // swap(TOTevent.at(j), TOTevent.at(j + 1));
      }
    }
  }
}
double SelectRecT(vector<double> &RecT, double deltaT , double w) // deltaT=1000
{
  int N = RecT.size();
  if (N < 2)
    return 0;
  double mean = 0;
  for (int i = 0; i < N; i++)
    mean += RecT[i] / N;
  if (RecT[N - 1] - RecT[0] < deltaT)
    return mean;
  if ((RecT[N - 1] - mean) * w > mean - RecT[0])
  {
    RecT.erase(RecT.begin() + N - 1);
    // channel.erase(channel.begin() + N - 1);
    // TOT.erase(TOT.begin() + N - 1);
  }

  else
  {
    RecT.erase(RecT.begin());
    // channel.erase(channel.begin());
    // TOT.erase(TOT.begin());
  }

  return SelectRecT(RecT, deltaT, w);
}
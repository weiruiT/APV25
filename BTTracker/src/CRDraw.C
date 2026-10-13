#include "CRDraw.h"


CRDraw::CRDraw(TString fanaName)
{
  anaName = fanaName;
}

CRDraw::~CRDraw()
{
  if (waveAna) delete waveAna;
}


bool CRDraw::Init()
{
  fanaFile = new TFile(anaName, "read");
  if (!fanaFile->IsOpen())
  {
    cout << " ERROR! " << anaName << " cant open!" << endl;
    return false;
  }
  cout << " --> Data from " << anaName << endl;
  for (int i = 0; i < 2 * MMdetN; i++)
  {
    fanaTree[i] = (TTree *)fanaFile->Get(Form("Track%s%d", XYname(i % 2).Data(), i / 2));
    if (fanaTree[i] == 0)
    {
      cout << " ERROR! Track" << XYname(i % 2) << i / 2 << " Tree not exist in " << anaName << endl;
      return false;
    }
    fanaTree[i]->SetBranchAddress("event", &event, &b_event);
    fanaTree[i]->SetBranchAddress("maxT", &maxT, &b_maxT);
    fanaTree[i]->SetBranchAddress("channel", &channel, &b_channel);
    fanaTree[i]->SetBranchAddress("ped", &ped, &b_ped);
    fanaTree[i]->SetBranchAddress("pedrms", &pedrms, &b_pedrms);
    fanaTree[i]->SetBranchAddress("riseT", &riseT, &b_riseT);
    fanaTree[i]->SetBranchAddress("TOT", &TOT, &b_TOT);
    fanaTree[i]->SetBranchAddress("wave", &waveAna, &b_wave);
    fanaTree[i]->SetBranchAddress("maxA", &maxA, &b_maxA);
  }




  return true;
}

bool CRDraw::select_waveform()
{
  const double t = 3 * pedrms;
  double b = 5 * pedrms;
  const auto& w = *waveAna;
  int maxposition = 0;

  vector<int> peak_positions;
  // 检测所有极大值（包含平台型）
  for (size_t i = 1; i < w.size() - 1; ++i) 
  {
    // 普通极大值
    if (w[i] > w[i-1] && w[i] > w[i+1] && w[i] - ped > b) 
    {
      peak_positions.push_back(i);
    }
      // 平台型极大值
    else if (i < w.size()-2 && w[i] == w[i+1] && w[i] > w[i-1] && w[i] > w[i+2] && w[i] - ped > t) 
    {
        peak_positions.push_back(i);
        ++i; // 跳过重复点
    }
    else if (i < w.size()-3 && w[i] == w[i+1] && w[i] == w[i+2] && w[i] > w[i-1] && w[i] > w[i+3] && w[i] - ped > t) 
    {
        peak_positions.push_back(i);
        i = i + 2; // 跳过重复点
    }
  }



  if (!peak_positions.empty()) 
  {
    // 寻找全局最大峰值
    int max_peak_idx = 0;
    int max_peak_value = w[peak_positions[0]];
    
    // 遍历所有峰值找到最大值
    for (size_t p = 1; p < peak_positions.size(); ++p) 
    {
      if (w[peak_positions[p]] > max_peak_value) 
      {
        max_peak_value = w[peak_positions[p]];
        max_peak_idx = p;
        maxposition = p;
      }
    }
    if (peak_positions[max_peak_idx] > 17 || peak_positions[max_peak_idx] < 7 || w[0] >= w[peak_positions[max_peak_idx]] || w[peak_positions[max_peak_idx]] - ped < t) 
    {
      return false;
    }
  }
  else
  {
    return false;
  }


  // 检查相邻峰值间的最小值差
  if (peak_positions.size() == 1) return true;
  if (peak_positions.size() > 1)
  {
    for (size_t i = 1; i < peak_positions.size(); ++i) 
    {
      int start = peak_positions[i-1];
      int end = peak_positions[i];
      
      int min_val = INT_MAX;
      for (int j = start; j < end + 1; ++j) 
      {
        if (w[j] < min_val) 
        {
          min_val = w[j];
        }
      }
      /*if (i == 1 && maxposition != peak_positions[0])
      {
        if (w[peak_positions[i - 1]] - min_val >= b) 
        {
          return false;
        }
      }*/
      if (i == maxposition) continue;
      if (w[peak_positions[i]] - min_val >= b) 
      {
        return false;
      }
    }
  }

  return true;
}

bool CRDraw::Cut()
{
  return maxA > pedrms * rmscut && maxT > 0 && maxT < 20 && riseT > 0 && riseT < 15 && TOT > 2 && TOT < 14;
}


void CRDraw::Cal_ped()
{
  // 创建输出文件
  TFile* outFile = new TFile("Pedestal_Analysis.root", "RECREATE");
  
  // 数据结构：存储每个通道的wave0值 [8 trees][256 channels]
  vector<short> channelped[2 * MMdetN][256];

  // 第一阶段：收集所有通道的wave0数据
  for (int iTree = 0; iTree < 2 * MMdetN; iTree++) 
  {
    if (!fanaTree[iTree]) {
      cerr << "Error: Tree " << iTree << " not initialized!" << endl;
      continue;
    }

    int nentries = fanaTree[iTree]->GetEntriesFast();
    for (int iEntry = 0; iEntry < nentries; iEntry++) 
    {
      fanaTree[iTree]->GetEntry(iEntry);
      
      // 数据有效性检查
      if (!waveAna || waveAna->empty()) {
        cerr << "Invalid waveform data at tree " << iTree 
             << " entry " << iEntry << endl;
        continue;
      }

      // 通道号范围检查
      if (channel < 0 || channel >= 256) {
        cerr << "Invalid channel " << channel << " at tree " << iTree 
             << " entry " << iEntry << endl;
        continue;
      }

      // 存储第一个采样点值
      channelped[iTree][channel].push_back(waveAna->at(0));
    }
  }

  // 第二阶段：计算统计量并保存结果
  for (int iTree = 0; iTree < 2 * MMdetN; iTree++) 
  {
    // 获取探测器信息
    int trackerID = iTree / 2;       // 0-3对应tracker0-3
    int direction = iTree % 2;       // 0:X, 1:Y
    TString dirName = Form("Tracker%d%s", trackerID, (direction == 0) ? "X" : "Y");
    
    // 创建目录结构
    outFile->mkdir(dirName);
    outFile->cd(dirName);

    // 创建直方图容器
    TH1F* hPedMean = new TH1F("hPedMean", "Pedestal Mean;Channel;Mean [ADC]",
                              256, 0, 256);
    TH1F* hPedRMS = new TH1F("hPedRMS", "Pedestal RMS;Channel;RMS [ADC]",
                             256, 0, 256);

    // 通道级处理
    for (int ch = 0; ch < 256; ch++) 
    {
      // 跳过空通道
      if (channelped[iTree][ch].empty()) continue;

      // 计算统计量
      double sum = 0.0, sumSq = 0.0;
      for (auto val : channelped[iTree][ch]) {
        sum += val;
        sumSq += val * val;
      }
      double mean = sum / channelped[iTree][ch].size();
      double rms = sqrt(sumSq/channelped[iTree][ch].size() - mean*mean);
      // 填充统计直方图
      hPedMean->SetBinContent(ch+1, mean);
      hPedRMS->SetBinContent(ch+1, rms);

      // 创建并保存波形直方图
      TString histName = Form("hWave0_Ch%d", ch);
      TH1F* hWave0 = new TH1F(histName, Form("Channel %d Wave0;ADC;Counts", ch),
                              819, -2048, 2047);
      
      for (auto val : channelped[iTree][ch]) 
      {
        hWave0->Fill(val);
      }
      int maxBin = hWave0->GetMaximumBin();           // 最高计数bin的索引
      double peakCenter = hWave0->GetBinCenter(maxBin); // 对应的中心值
      
      // 2. 计算动态范围（峰值±200）
      double fitLow = peakCenter - 200;
      double fitHigh = peakCenter + 200;

      // 3. 边界保护（不超过直方图范围）
      fitLow = std::max(fitLow, hWave0->GetXaxis()->GetXmin());
      fitHigh = std::min(fitHigh, hWave0->GetXaxis()->GetXmax());
      double histMean = hWave0->GetMean();
      double histRMS = hWave0->GetRMS();
      // 4. 执行拟合
      TF1 *gaussFit = new TF1("gaussFit", "gaus", fitLow, fitHigh);
      gaussFit->SetParameters(
          hWave0->GetMaximum(), // 初始幅度（峰值高度）
          histMean,             // 初始均值（直方图平均值）
          histRMS               // 初始标准差（直方图RMS）
      );
      hWave0->Fit(gaussFit, "QR", "", fitLow, fitHigh); // 指定范围和静默模式
      TF1 *fitFunc = hWave0->GetFunction("gaussFit");
      if (fitFunc) 
      {
        double fit_mean = fitFunc->GetParameter(1);
        double fit_sigma = fitFunc->GetParameter(2);
        channelmean[iTree][ch] = fit_mean;
        channelrms[iTree][ch] = fit_sigma;
        hPedMean->SetBinContent(ch+1, fit_mean);
        hPedRMS->SetBinContent(ch+1, fit_sigma);
        fitFunc->SetLineColor(kRed); // 设置拟合曲线颜色
        fitFunc->SetLineWidth(2);        // 线宽设为2像素
        fitFunc->SetNpx(1000);          // 提高曲线分辨率
      } 
      else 
      {
        hPedMean->SetBinContent(ch+1, hWave0->GetMean());
        hPedRMS->SetBinContent(ch+1, hWave0->GetRMS());
        cerr << "Fit failed for channel " << ch << endl;
      }
      hWave0->Write();
      delete hWave0;
    }

    // 保存统计直方图
    outFile->cd(dirName);
    hPedMean->Write();
    hPedRMS->Write();
    
    // 清理内存
    delete hPedMean;
    delete hPedRMS;
  }

  // 关闭文件
  outFile->Close();
  delete outFile;
  cout << "Pedestal analysis saved to Pedestal_Analysis.root" << endl;
}

void CRDraw::Draw_noise_es()
{
  // 创建ROOT文件（所有Tree的直方图保存到同一个文件）
  TFile *outFile = new TFile("0wave_Distributions.root", "RECREATE");

  ofstream paramFile("fit_parameters.txt");
  if (!paramFile.is_open()) {
    cerr << "Error: Could not open fit_parameters.txt for writing!" << endl;
    return;
  }
  paramFile << "# Channel parameters (mean, sigma)\n";
  paramFile << "# Format: global_channel mean sigma\n";

  TH1F *hAllMean = new TH1F("hAllMean", "All Channels - Gaussian Mean;Global Channel;Mean", 12*256, 0, 12*256); // 8个Tree * 256通道
  TH1F *hAllSigma = new TH1F("hAllSigma", "All Channels - Gaussian Sigma;Global Channel;Sigma", 12*256, 0, 12*256);

  hAllMean->Sumw2(kFALSE);             // 关键：关闭误差计算
  hAllSigma->Sumw2(kFALSE);

  TH1F *hMeanDist = new TH1F("hMeanDist", "Distribution of Gaussian Mean;Mean;Counts", 100, -1000, 1000); // 假设均值范围在-1000到1000
  TH1F *hSigmaDist = new TH1F("hSigmaDist", "Distribution of Gaussian Sigma;Sigma;Counts", 100, 0, 500);     // 假设标准差范围在0到500

  for (int iTree = 0; iTree < 2 * MMdetN; iTree++)
  {

    int trackerID = iTree / 2;      // Tracker0, Tracker1, ...
    int direction = iTree % 2;      // 0=X, 1=Y
    TString dirName = Form("Tracker%d%s", trackerID, (direction == 0) ? "X" : "Y");

    // 在 ROOT 文件中创建目录
    outFile->mkdir(dirName);
    outFile->cd(dirName);

    int nentries = fanaTree[iTree]->GetEntriesFast();
    cout << "Processing Tree " << iTree << " (" << nentries << " entries)" << endl;

    // 创建当前Tree的直方图数组
    TH1F *hwave[256];
    for (int j = 0; j < 256; j++) 
    {
      // 唯一命名：Tree索引+通道号
      TString histName = Form("hwave_Ch%d", j);
      TString histTitle = Form("Channel%d: wave Distribution;wave;Counts", j);
      hwave[j] = new TH1F(histName, histTitle, 4096, -2048, 2048); // 调整bin范围和数量
    }

    // 填充数据
    for (int iEntry = 0; iEntry < nentries; iEntry++) 
    {
      DrawProcessbar(iEntry, nentries);
      fanaTree[iTree]->GetEntry(iEntry);
      int ch = channel;
      if (!waveAna || waveAna->size() != 30) 
      {
        cerr << "Error: Invalid waveform data in entry " << iEntry << endl;
        continue;
      }
    // 填充数据
    for (int b = 0; b < 30; b++) 
    {
        hwave[ch]->Fill(waveAna->at(b));
    }
    }

    for (int j = 0; j < 256; j++) 
    {
      if (hwave[j]->GetEntries() < 300) 
      {
        // 对于无效通道，写入占位值
        int globalChannel = iTree * 256 + j + 1;
        paramFile << globalChannel << " " << 0.0 << " " << 0.0 << "\n";
        continue;
      }

      // 定义高斯函数并拟合
      outFile->cd(dirName); // 显式切换回目标目录
            // 动态确定拟合范围：最大值附近±300
      int maxBin = hwave[j]->GetMaximumBin();
      double maxCenter = hwave[j]->GetBinCenter(maxBin);
      double fitMin = maxCenter - 300;
      double fitMax = maxCenter + 300;
      fitMin = std::max(fitMin, hwave[j]->GetXaxis()->GetXmin());
      fitMax = std::min(fitMax, hwave[j]->GetXaxis()->GetXmax());
      
      TF1 *gaussFit = new TF1("gaussFit", "gaus", fitMin, fitMax);
      //TF1 *gaussFit = new TF1("gaussFit", "gaus", hwave[j]->GetXaxis()->GetXmin(), hwave[j]->GetXaxis()->GetXmax());
      gaussFit->SetParameters(hwave[j]->GetMaximum(), hwave[j]->GetMean(), hwave[j]->GetRMS());
      hwave[j]->Fit(gaussFit, "RQ"); // 关键：拟合结果自动绑定到直方图

      // 设置拟合曲线样式
      gaussFit->SetLineColor(kRed);
      gaussFit->SetLineWidth(2);
      int globalChannel = iTree * 256 + j + 1; // 关键公式
      double mean = gaussFit->GetParameter(1);
      double sigma = gaussFit->GetParameter(2);
      paramFile << globalChannel << " " << mean << " " << sigma << "\n";
      hAllMean->SetBinContent(globalChannel, gaussFit->GetParameter(1));
      hAllMean->SetBinError(globalChannel, 0);
      hAllSigma->SetBinContent(globalChannel, gaussFit->GetParameter(2));
      hAllSigma->SetBinError(globalChannel, 0); // 关键：误差设为0
      hMeanDist->Fill(mean);
      hSigmaDist->Fill(sigma);
      // 写入直方图时自动包含拟合函数
      hwave[j]->Write();
      delete hwave[j]; 
      delete gaussFit; 
    }
    
    // 将当前Tree的直方图写入文件
    // 可选：保留PDF输出（若需要）
    // TCanvas* canvas = new TCanvas(...);
    // ...绘制代码...
    // canvas->SaveAs(Form("wave%d.pdf", iTree));
    // delete canvas;
  }

  outFile->cd();
  hAllMean->Write();
  hAllSigma->Write();
  hMeanDist->Write();
  hSigmaDist->Write();

  // 关闭ROOT文件并释放资源
  outFile->Close();
  delete outFile;
  paramFile.close();
  cout << "All histograms saved to wave_Distributions.root" << endl;
}

void CRDraw::Draw_channel_waveform(int board, int XY, int Dchannel)
{
  TString fileName = Form("Tracker%d%s_channel%d.root", board, (XY == 0) ? "X" : "Y", Dchannel);
  TFile *outFile = new TFile(fileName, "RECREATE");
  int nentries = fanaTree[board*2 + XY]->GetEntriesFast();
  cout << "nentries:" << nentries << endl;
  vector<int> effevent;
  int eventnum = 0;
  for(int i = 0; i < nentries; i++)
  {
    fanaTree[board*2 + XY]->GetEntry(i);
    if (channel == Dchannel && waveAna->size() == 30)
    {
      eventnum++;
      effevent.push_back(event);
      cout << event << endl;
    }
    else
      continue;
  }
  cout << "effevent:" << eventnum <<endl;
  TString dirName = Form("Tracker%d%s_channel%d", board, (XY == 0) ? "X" : "Y", Dchannel);
  outFile->mkdir(dirName);
  outFile->cd(dirName);
  TGraph** gwave = new TGraph*[eventnum];  // TGraph指针数组
  
  // 第二遍循环：填充波形数据
  int e = 0;
  for(int i = 0; i < nentries; i++) {
    fanaTree[board*2 + XY]->GetEntry(i);
    if (channel == Dchannel && waveAna->size() == 30) {
      // 创建x,y数组（横坐标1-30）
      double x[30], y[30];
      for(int b = 0; b < 30; b++) {
        x[b] = b + 1;         // 横坐标从1开始
        y[b] = waveAna->at(b); // 纵坐标为ADC值
      }
      
      // 创建TGraph并设置属性
      TString graphName = Form("Tracker%d%s_channel%d_event%d", 
                              board, (XY == 0) ? "X" : "Y", 
                              Dchannel, effevent[e]);
      TString graphTitle = Form("Tracker%d%s_channel%d_event%d: Waveform;Sample Number;ADC Value",
                               board, (XY == 0) ? "X" : "Y", 
                               Dchannel, effevent[e]);
                               
      gwave[e] = new TGraph(30, x, y);
      gwave[e]->SetName(graphName);
      gwave[e]->SetTitle(graphTitle);
      gwave[e]->SetLineColor(kBlue);  // 设置连线颜色
      e++;
    }
  }

  // 写入文件并清理内存
  for(int j = 0; j < eventnum; j++) {
    gwave[j]->Write();
    delete gwave[j];
  }
  delete[] gwave;
  outFile->cd();
  outFile->Close();
  delete outFile;
}

void CRDraw::Draw_all_waveforms()
{
  Cal_ped();
  TFile* outFile = new TFile("AllTrackers_Waveforms_slc.root", "RECREATE");
  const int totalTrackers = MMdetN;
  
  // 总进度：跟踪处理完成的Tracker数量
  cout << "==== Total Progress ====" << endl;
  for (int board = 0; board < totalTrackers; board++) {
    // 创建Tracker目录
    TString trackerTopDir = Form("Tracker%d", board);
    outFile->mkdir(trackerTopDir);
    
    // 处理X/Y方向
    for (int XY = 0; XY < 2; XY++) {
      int treeIndex = board * 2 + XY;
      if (!fanaTree[treeIndex]) continue;

      // 预分类数据结构
      std::map<int, std::vector<std::pair<int, int>>> channelData;
      
      // 第一阶段：收集有效事件
      int nentries = fanaTree[treeIndex]->GetEntries();
      cout << "\n==== Tracker" << board << (XY == 0 ? "X" : "Y") 
           << " 开始处理, 总事件数: " << nentries << " ====" << endl;
      
      for (int ientry = 0; ientry < nentries; ientry++) {
        fanaTree[treeIndex]->GetEntry(ientry);
        if (waveAna->size() != 30) continue;
        //if (!Cut()) continue;
        if (!select_waveform()) continue;
        int currentChannel = channel;
        if (currentChannel < 0 || currentChannel >= 256) continue;
        
        channelData[currentChannel].emplace_back(ientry, event);
      }
      cout << endl; // 结束当前进度条

      // 第二阶段：生成波形数据
      TString dirXY = (XY == 0) ? "X" : "Y";
      TString dirPath = Form("%s/%s", trackerTopDir.Data(), dirXY.Data());
      outFile->mkdir(dirPath);

      // 通道进度统计
      int totalChannels = channelData.size();
      int processedChannels = 0;
      cout << "==== 开始生成波形数据 ====" << endl;
      
      for (auto& [channel, entries] : channelData) {
        // 更新通道级进度
        cout << "\n  Channel " << channel << " (" << entries.size() << " events) ";
        TString chDir = Form("%s/Channel%d", dirPath.Data(), channel);
        outFile->mkdir(chDir);
        outFile->cd(chDir);

        // 事件处理进度
        int totalEvents = entries.size();
        for (size_t ievent = 0; ievent < entries.size(); ++ievent) {
          auto& [ientry, event] = entries[ievent];
          fanaTree[treeIndex]->GetEntry(ientry);
          
          double x[30], y[30];
          for (int b = 0; b < 30; b++) {
            x[b] = b + 1;
            y[b] = waveAna->at(b) - channelmean[treeIndex][channel];
          }

          TString graphName = Form("Event%d", event);
          TGraph* g = new TGraph(30, x, y);
          g->SetName(graphName);
          g->SetTitle(Form("Tracker%d%s-Ch%d Event%d;Sample;ADC", 
                         board, dirXY.Data(), channel, event));
          g->Write();
          delete g;

        }
        cout << endl; // 结束通道进度条
        processedChannels++;
        
        // 更新通道处理进度
        cout << "  ==> Processed channels: ";
        DrawProcessbar(processedChannels, totalChannels);
      }
    }
    
    // 更新总进
  }
  
  outFile->Close();
  delete outFile;
  cout << "\n==== 全部处理完成 ====" << endl;
}

void CRDraw::Draw_channel_ped()
{
  TFile* outFile = new TFile("channel_ped.root", "RECREATE");
  const int numTrackers = MMdetN;

  for (int board = 0; board < numTrackers; board++) {     // 循环Tracker
    for (int XY = 0; XY < 2; XY++) {                      // 循环X/Y方向
      int treeIndex = board * 2 + XY;
      if (!fanaTree[treeIndex]) continue;

      // 创建Tracker目录 (例如: Tracker0/X)
      TString trackerDir = Form("Tracker%d/%s", board, (XY == 0) ? "X" : "Y");
      outFile->mkdir(trackerDir);
      if (!outFile->cd(trackerDir)) {
        cerr << "Error creating directory: " << trackerDir << endl;
        continue;
      }

      // 预创建所有通道的直方图
      TH1F* hPed[256];
      for (int ch = 0; ch < 256; ch++) {
        TString histName = Form("hPed_Tracker%d%s_Ch%d", board, (XY==0?"X":"Y"), ch);
        TString histTitle = Form("Tracker%d%s Ch%d;Ped;Counts", board, (XY==0?"X":"Y"), ch);
        hPed[ch] = new TH1F(histName, histTitle, 200, -2000, 0); // 假设ped范围±500
      }

      // 填充数据
      int nEntries = fanaTree[treeIndex]->GetEntries();
      for (int iEntry = 0; iEntry < nEntries; iEntry++) {
        fanaTree[treeIndex]->GetEntry(iEntry);
        if (channel < 0 || channel >= 256) continue;
        hPed[channel]->Fill(ped); // 填充当前通道的ped值
      }

      // 写入并清理内存
      outFile->cd(trackerDir);
      for (int ch = 0; ch < 256; ch++) {
        if (hPed[ch]->GetEntries() > 0) {
          hPed[ch]->Write();
        }
        delete hPed[ch];
      }
    }
  }
  outFile->Close();
  delete outFile;
  cout << "Ped histograms saved to channel_ped.root" << endl;
}

void CRDraw::Draw_event_waveforms(int eventID)
{
  // 创建输出文件
  TString fileName = Form("event%d.root", eventID);
  TFile* outFile = new TFile(fileName, "RECREATE");
  
  cout << "==== 开始处理事件 " << eventID << " ====" << endl;
  
  // 遍历所有board
  for (int board = 0; board < MMdetN; board++) {
    // 创建board目录
    TString boardDir = Form("board%d", board);
    outFile->mkdir(boardDir);
    
    cout << "处理 Board " << board << "..." << endl;
    
    // 遍历X和Y方向
    for (int XY = 0; XY < 2; XY++) {
      int treeIndex = board * 2 + XY;
      if (!fanaTree[treeIndex]) continue;
      
      // 创建通道波形存储容器
      std::map<int, TGraph*> channelGraphs;  // 通道号 -> TGraph
      
      // 搜索指定事件的所有通道
      int nentries = fanaTree[treeIndex]->GetEntries();
      for (int ientry = 0; ientry < nentries; ientry++) {
        fanaTree[treeIndex]->GetEntry(ientry);
        
        // 检查是否为指定事件
        if (event != eventID) continue;
        
        // 检查波形数据是否有效
        if (!waveAna || waveAna->size() != 30) continue;
        
        int currentChannel = channel;
        if (currentChannel < 0 || currentChannel >= 256) continue;
        
        // 创建波形图
        double x[30], y[30];
        for (int b = 0; b < 30; b++) {
          x[b] = b + 1;
          y[b] = waveAna->at(b);
        }
        
        TString graphName = Form("ch%d", currentChannel);
        TGraph* graph = new TGraph(30, x, y);
        graph->SetName(graphName);
        graph->SetTitle(Form("Board%d-Channel%d;Sample;ADC", board, currentChannel));
        
        channelGraphs[currentChannel] = graph;
      }
      
      // 将找到的通道波形保存到文件
      if (!channelGraphs.empty()) {
        TString xyDir = Form("%s/%s", boardDir.Data(), (XY == 0) ? "X" : "Y");
        outFile->mkdir(xyDir);
        outFile->cd(xyDir);
        
        for (auto& [ch, graph] : channelGraphs) {
          // 创建通道目录
          TString chDir = Form("%s/channel%d", xyDir.Data(), ch);
          outFile->mkdir(chDir);
          outFile->cd(chDir);
          
          graph->Write();
          delete graph;
          
          cout << "  -> Board" << board << (XY == 0 ? "X" : "Y") 
               << " Channel" << ch << " 波形已保存" << endl;
        }
      }
    }
  }
  
  outFile->Close();
  delete outFile;
  
  cout << "==== 事件 " << eventID << " 的所有波形已保存到 " << fileName << " ====" << endl;
}

void CRDraw::Draw_events_from_file(const char* filename)
{
  // 读取事件ID文件
  ifstream eventFile(filename);
  if (!eventFile.is_open()) {
    cerr << "Error: Cannot open event file " << filename << endl;
    return;
  }

  vector<int> eventList;
  int eventID;
  while (eventFile >> eventID) {
    eventList.push_back(eventID);
  }
  eventFile.close();

  cout << "Found " << eventList.size() << " events in " << filename << endl;

  // 创建统一的输出文件
  TString outFileName = "events_waveforms_select.root";
  TFile* outFile = new TFile(outFileName, "RECREATE");
  
  TH1F *maxamp;
  maxamp = new TH1F(Form("maxA"), Form("maxA"), 500, 0, 4000);
  maxamp->GetXaxis()->SetTitle("maxA");
  maxamp->GetYaxis()->SetTitle("counts");
  // 进度显示
  int totalEvents = eventList.size();
  int processedEvents = 0;

  cout << "==== 开始处理所有事件 ====" << endl;

  // 处理每个事件
  int prenentries[2 * MMdetN] = {0};
  for (auto targetEvent : eventList) {
    cout << "处理事件 " << targetEvent << " (" << ++processedEvents << "/" << totalEvents << ")..." << endl;

    // 创建事件目录
    TString eventDirName = Form("event%d", targetEvent);
    outFile->mkdir(eventDirName);
    
    // 遍历所有board
    
    for (int board = 0; board < MMdetN; board++) {
      // 创建board目录
      TString boardDirName = Form("%s/board%d", eventDirName.Data(), board);
      outFile->mkdir(boardDirName);
      
      // 遍历X和Y方向
      for (int XY = 0; XY < 2; XY++) {
        int treeIndex = board * 2 + XY;
        if (!fanaTree[treeIndex]) continue;

        // 搜索指定事件的所有通道
        int nentries = fanaTree[treeIndex]->GetEntries();
        std::map<int, TGraph*> channelGraphs;  // 通道号 -> TGraph
        bool searchevent = 0;
        for (int ientry = prenentries[treeIndex]; ientry < nentries; ientry++) {
          fanaTree[treeIndex]->GetEntry(ientry);
          prenentries[treeIndex] = ientry;
          // 检查是否为指定事件
          if ((event) != targetEvent && searchevent) break;
          if ((event) != targetEvent) continue;
          
          // 检查波形数据是否有效
          if (!waveAna || waveAna->size() != 30) continue;
          
          int currentChannel = channel;
          if (currentChannel < 0 || currentChannel >= 256) continue;
          
          // 创建波形图
          searchevent = 1;
          if (!select_waveform()) continue;
          if (treeIndex == 2) maxamp->Fill(maxA);
          double x[30], y[30];
          for (int b = 0; b < 30; b++) {
            x[b] = b + 1;
            y[b] = waveAna->at(b) - ped;
          }
          
          TString graphName = Form("ch%d", currentChannel);
          TGraph* graph = new TGraph(30, x, y);
          graph->SetName(graphName);
          graph->SetTitle(Form("Event%d-Board%d-%s-Channel%d;Sample;ADC", 
                             targetEvent, board, (XY == 0) ? "X" : "Y", currentChannel));
          
          channelGraphs[currentChannel] = graph;
        }
        
        // 将找到的通道波形保存到文件
        if (!channelGraphs.empty()) {
          TString xyDir = Form("%s/%s", boardDirName.Data(), (XY == 0) ? "X" : "Y");
          outFile->mkdir(xyDir);
          outFile->cd(xyDir);
          
          for (auto& [ch, graph] : channelGraphs) {
            // 创建通道目录
            TString chDir = Form("%s/channel%d", xyDir.Data(), ch);
            outFile->mkdir(chDir);
            outFile->cd(chDir);
            
            graph->Write();
            delete graph;
          }
        }
      }
    }
    
    // 更新进度
    cout << "  -> 事件 " << targetEvent << " 完成 (" << processedEvents << "/" << totalEvents << ")" << endl;
  }
  maxamp->Write();
  outFile->Close();
  delete outFile;
  
  cout << "==== 所有事件处理完成 ====" << endl;
  cout << "结果保存在: " << outFileName << endl;
}
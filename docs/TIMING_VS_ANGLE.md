# BM2 时间分辨率随径迹角度变化

## 直接运行

在 WSL/Linux 的终端中，从 APV25 根目录执行：

```bash
cd /home/thoriatong/workplace/APV25
source /home/thoriatong/root/bin/thisroot.sh
root -l -b -q 'BTTracker/macros/PlotTimingVsAngle.C("batches/run_010")'
```

宏调用 `BTTracker/scripts/plot_timing_vs_angle.py`，复用正式时间分析中的事件读取、迭代 3σ 核心估计及加权两两方差求解。使用已有的 `Combine/Tracker-step4-rec.root` 和 `Combine/Pad-step3-rec.root`，无需修改 ROOT 分支或重跑重建。默认 Python 为 `/usr/bin/python3`；可用 `BTTRACKER_PYTHON` 指定能导入 PyROOT 的解释器。

默认分成近似等事件数的 5 个角度区间，每个探测器对在核心筛选后至少保留 40 个事件。相同角度的事件不会被拆开，因此重复角度较多时区间数可能少于请求值。改变分组数量：

```bash
root -l -b -q 'BTTracker/macros/PlotTimingVsAngle.C("batches/run_010",8)'
```

指定固定角度边界（单位：度），用于多个 run 之间比较：

```bash
root -l -b -q 'BTTracker/macros/PlotTimingVsAngle.C("batches/run_010",5,"0,3,6,9,12,18,25",40)'
```

第五个参数选择下图绘制的 CFD 定时方法：`0` 表示画全部四条曲线（默认行为），`1`、`2`、`3`、`4` 分别只画 `t20`、`t30`、`t40`、`t50`。例如只画 `t20`：

```bash
root -l -b -q 'BTTracker/macros/PlotTimingVsAngle.C("batches/run_010",5,"",40,1)'
```

这个编号只控制图中显示的曲线；`angle_bins.tsv` 和 `analysis.json` 仍保存四种 CFD 的计算结果。

若同时使用固定角度边界，仍按原顺序传参，再把定时方法编号放在第五个参数：

```bash
root -l -b -q 'BTTracker/macros/PlotTimingVsAngle.C("batches/run_010",5,"0,3,6,9,12,18,25",40,1)'
```

第三个参数非空时使用固定边界，忽略第二个参数的区间数。区间左闭右开，最后一个区间包含右端点。范围外事件数量记录在 `analysis.json`。统计不足、矩阵不可解或出现非正方差的分组保留状态，不在图上画成零。

也可以直接使用 Python，显式指定 BM2 有效匹配 pad 数目：

```bash
/usr/bin/python3 BTTracker/scripts/plot_timing_vs_angle.py batches/run_010 \
  --bins 5 --min-events 40 --bm2-matched-pad-count 1 --timing-method 1
```

## 样本与图像含义

- 角度定义为相对 z 轴的夹角：`theta = atan(hypot(XfitPars[0], YfitPars[0]))`，转换成度。另存有符号的 x、y 投影角。
- Step4 只取 `useflag=true` 的径迹，按 `event` 与时间树关联，并检查 `timestamp`；缺失、重复关联或时间戳不匹配时停止。
- BM2 pad 数目 cut 默认读取当前主时间分析的 `Combine/time_resolution/time_resolution.json`。没有可用记录时默认 `bmS_n[1]==1`。本工具分析单 pad/window 时间，不分析 cluster 时间。
- 为使上方角度分布和四条 CFD 曲线可直接比较，使用 BM2、Tracker1–3 在 t20/t30/t40/t50 均有有效平均时间的共同事件。不会额外施加 SNR cut。这个共同样本可能小于主分析中各探测器对分别使用的有效样本，因此整体分辨数值不必与主分析完全一致。
- 每个角度区间使用六组时间差，复用正式分析的 median/MAD 初始化、迭代 3σ 核心选择和加权方差求解，得到 BM2 分辨率。不是把单个事件的时间差当成分辨率。
- 上图为共同样本的角度分布，虚线是角度区间边界。下图横坐标为区间内共同有效事件的角度算术平均值（时间差核心筛选前），每个区间底部标注平均值，单位为度、保留两位小数；`angle_bins.tsv` 的 `theta_mean_deg` 保存完整精度。横向误差棒从平均值延伸到区间边界；纵向误差棒是现有加权方差模型传播的近似统计误差，未包含探测器对之间的协方差及系统误差。
- 上图每个区间标注的 `N` 是四种 CFD 共用、尚未进行两两时间差核心筛选的事件数。下图纵轴按所有有效点的误差棒上下端自动取范围，并留出空白；不再固定从 0 开始。

解释角度依赖时，应同时考虑 SNR 和 pad 多重性随角度变化，以及探测器时间误差独立的方差模型假设。若需区分这些影响，可后续在相同 SNR 范围内作比较。

这里的“Tracker1–3 时间共同有效”不等于三层 X/Y 六个方向全部有空间响应。
`useflag` 来自两个空间投影的拟合，每个投影允许只使用两层测量点；每层 Tracker 时间又允许
X/Y 至少一个方向有效。脚本没有额外要求三层 X/Y 全有信号、三层点均在 ±50 mm 内，
也没有自动施加 17° 上限或 BM2 矩形位置 cut。因此不能直接把假定三层完整接受区得到的
几何角度上限当作该图已有的 cut。各入口的具体条件见 [分析样本与 cut 对照](ANALYSIS_SELECTIONS.md)。

## 输出

每次创建一个新的目录：

```text
run_NNN/Combine/time_resolution/angle_binned/<执行编号>/
├── timing_vs_angle.png       # 上下两幅汇总图
├── timing_vs_angle.pdf
├── timing_vs_angle.root      # TCanvas、角度直方图和所选分辨率曲线
├── angle_bins.tsv            # 每组统计量、分辨率、误差和状态
├── selected_events.tsv       # event、角度与各探测器的四种时间
└── analysis.json             # 输入记录、样本定义、两两统计量及方差解
```

这是独立绘图入口，尚未接入 run manager 的自动执行和 Results 白名单；生成的文件保存在 Combine。既有 Results 归档不受此操作影响。

在 ROOT 窗口中查看保存的画布：

```cpp
auto f = TFile::Open("batches/run_010/Combine/time_resolution/angle_binned/<执行编号>/timing_vs_angle.root");
auto c = (TCanvas*)f->Get("timing_vs_angle");
c->Draw();
```

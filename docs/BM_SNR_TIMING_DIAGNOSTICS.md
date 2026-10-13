# BM2 信噪比与时间分辨诊断

从 APV25 项目根目录执行（只需要 PyROOT 和 Python 标准库）：

```bash
source /home/thoriatong/root/bin/thisroot.sh
python3 debug/macros/analyze_bm_snr_timing.py \
    batches/run_006 batches/run_008 batches/run_009
```

默认结果写入新的 `debug/reports/bm2_snr_<UTC时间>/`，也可以用
`--output-dir /absolute/new/directory` 指定一个尚不存在的目录。
诊断读取现有分析结果；不会修改输入、生产 ROOT、配置或主分析的 cut。

## 指标

- `SNR_peak = maxA / peak_noise_rms`，线性幅度比，无单位。
- `peak_noise_rms` 是 Step2 当时读取的本 run、对应 `global_channel` 的 noise sigma；
  宇宙线已扣除基线，`ped=0`、`pedrms=0`，不能以 `pedrms` 为分母。
- 当前 sigmoid 模型为 `B + H/(1+exp(-(s-t50)/k))`，s 和 k 以 APV 采样间隔为单位，
  每个间隔 25 ns。20% 恒比点的斜率为 `0.16*H/(25*k)` ADC/ns。
- `jitter20_proxy = 25*k*peak_noise_rms/(0.16*H)` ns，是局部过阈的噪声/斜率代理量，
  不是 sigmoid 多点拟合的协方差误差，也不是整个 BM2 时间分辨。
- `sigmoid_rise10_90_ns = 25*k*ln(81)`，用于检查波形形状是否随幅度变化。
- 当前 sigmoid 使用未指定逐点噪声的 `TGraph` 拟合，`sigmoidFitChi2NDF` 不能直接当作
  按 noise sigma 归一化的约化卡方。脚本额外保存 `sqrt(sigmoidFitChi2NDF)/peak_noise_rms`
  作为拟合残差/噪声的诊断量，但不据此额外筛选。

在 ROOT 中快速画出 Step2 的 BM2 波形 SNR：

```cpp
auto t = (TTree*)_file0->Get("BMPadSignals");
t->Draw("maxA/peak_noise_rms", "raw_chip==13 && peak_noise_rms>0");
```

`BMPadSignals` 的条目是波形，不是事件；两个 BM 都可能在此树中，必须按 `raw_chip`
区分。上述绘图不包含 Step4 径迹及 pad 匹配，不能直接当作时间分辨事件样本。

## 时间样本和方差求解

脚本复用 `BTTracker/scripts/analyze_time_resolution.py` 的读取和估计函数：

1. 读取 `TrackArrivalEvents`，使用生产条件 `bmS_n[1]==1`。
2. 根据 `bm_pred_ix[1]`、`bm_pred_iy[1]` 还原生产分析的 3×3 pad 窗口。
3. Step2 中相同 `(event, position)` 保留 maxA 最大的 chip13 波形，规则与生产代码一致。
4. 查找唯一具有有效 t20 的匹配 pad，并核对其 t20×25 与 `bmS_mean[1]` 精确相等。
5. 按该波形的 SNR 分组或增加 SNR 区间 cut，对所有六个探测器对都使用这一事件集合。
6. Tracker 时间继续按有效 X/Y strip 数量加权合并，各对使用其双方有时间的事件。
7. 各对使用 median/MAD 初始化的迭代 3σ 核心宽度，复用生产加权方差解
   `sigma_ij² = sigma_i² + sigma_j²`，从 BM2 和三个 Tracker 的六个方程估计四个方差。
8. 核对无附加 SNR cut 时的 BM2 分辨与已有 `time_resolution.json` 一致。

`snr_bins.tsv` 另外保存要求四个探测器同时有效的结果，便于检查各对样本不同带来的影响。
脚本输出误差与生产分析相同，是基于各对方差估计误差的近似传播；不包含共享事件造成的
探测器对相关性，也不包含噪声标定、拟合模型、对齐等系统误差。

## 结果文件

- `snr_timing_comparison.png/pdf`：归一化幅度、noise sigma、实际时间样本 SNR、各 SNR 组的 BM2 分辨。
- `snr_timing_mechanisms.png/pdf`：上升沿宽度、noise/20%斜率代理量、BM2−Tracker2 时间差均值随 SNR 的变化。
- `snr_threshold_scan.png/pdf`：SNR 下限与 BM2 分辨、保留事件比例。
- `snr_bins.tsv`、`snr_thresholds.tsv`、`snr_bands.tsv`：分组、下限 cut 和区间 cut 的详细结果。
- `matched_bm2_events.tsv`：逐事件关联到的 Step2 波形条目、pad、global_channel、幅度、噪声、SNR 和时间。
- `snr_timing_report.json`：定义、样本计数、分布统计和数值结果。
- `snr_timing_plots.root`：可继续查看的画布、直方图和曲线。

比较时区分：低 SNR 的随机定时抖动；幅度/形状相关的时间均值偏移；高幅区可能的峰顶
非线性或截断；以及气体、场、增益、漂移过程造成的真实时间差异。
SNR 分组相关性不能单独证明单面结构改变是因果来源，也不能保证 SNR cut 越高越好。
尤其高幅区出现均值移动或宽度回升时，需要查看波形后再判断是否饱和，不能只凭幅度谱的凸起定性。

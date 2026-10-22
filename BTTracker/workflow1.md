# 当前处理流程

## 输入与映射

宇宙线输入树为 `apv_EventSamples`/`apv_EventHeaders`，已经基线扣除和 3σ 零压缩。noise 输入树为 `apv_events`，保留原始 ADC 波形。

| 解包 ROOT `apv_board` | 物理连接 | 归一化 `TrackRaw` |
|---:|---|---|
| 0、1 | Tracker1 X | board 0，chip 0、1 |
| 2、3 | Tracker1 Y | board 0，chip 2、3 |
| 4、5 | Tracker2 X | board 1，chip 0、1 |
| 6、7 | Tracker2 Y | board 1，chip 2、3 |
| 8 | DUT / BM1 | board -1，chip 12 |
| 9 | DUT / BM2 | board -1，chip 13 |
| 10、11 | Tracker3 X | board 2，chip 0、1 |
| 12、13 | Tracker3 Y | board 2，chip 2、3 |
| 14、15 | 空载 | 忽略 |

因此源 `apv_board=12/13` 是 Tracker3 Y；只有适配后的内部 `chip=12/13` 才表示 BM1/BM2。表中是公共默认映射；接线交换必须通过该 run 的 `detector_layout.bm_slot8_logical_chip` / `bm_slot9_logical_chip` 覆盖，并会同时用于 noise 与 signal。run_005 当前覆盖为 slot8→chip13（BM2）、slot9→chip12（BM1，signal 空载）。

## 几何关系

沿 z 从上游到下游为 Tracker1=533、Tracker2=468、BM1=364、BM2=309、Tracker3=217 mm。原 BTTracker 的 Tracker4 位于 154 mm；当前配置只删除这一层，其余探测器位置和相对关系不变。所有重建和 DUT 径迹外推统一读取 `CRInfo.h` 中的 `Zpos`/`BMZpos`。

## signal

```text
DecodedAPVRoot2TrackRaw
  -> TrackRaw：slot 映射、trigger 展宽、多文件连续事件号、时间戳
AnalysisAPV
  -> TrackX0/Y0、TrackX1/Y1、TrackX2/Y2、BMPadSignals
  -> pedestal=0；保留波形形状、单峰条件，并从一次 sigmoid 拟合计算 t20/t30/t40/t50
DecodeAPV
  -> Track0、Track1、Track2：按六个 Tracker 方向的 event 并集压缩，重建 strip cluster 和二维 hit
  -> 显式 entry=1..N 在三棵树中同步；event 保留原始 DAQ trigger 编号
ReconstructTrack
  -> 三层唯一二维 cluster 联合对齐：固定 Tracker1，拟合 Tracker2/3 的 dx、dy、rotationZ
  -> Track：应用联合对齐后的三层直线拟合、残差、效率和预测 strip
PadReconstruction
  -> BM1Events、BM2Events、TrackArrivalEvents：独立 BM 事件表和径迹匹配时间
analyze_time_resolution.py
  -> 相同事件 cut 下并行分析 sigmoid t20/t30/t40/t50；默认选择 bmS_n[1]==1
  -> BM2 的径迹中心 3×3 窗口内恰好一个有效 sigmoid pad
analyze_time_resolution.py --bm-time-source cluster（独立的 cluster_timing 阶段）
  -> 唯一匹配的空间 cluster，允许多个 pad，使用完整 cluster 的有效时间
  -> Combine/time_resolution/cluster/<execution_id>/；原单 pad 输出继续保留
  -> 每种恒比时刻各自的 BM2−Tracker 时间差分布、核心选择、高斯拟合和时间分辨
  -> `timing_cuts.detectors` 中实际接入探测器的两两方差联合解与逐探测器时间分辨
`scripts/analyze_bm_snr_timing.py`（默认启用，可由 `snr_timing.enabled` 关闭）
  -> 用 Step2 的 maxA / peak_noise_rms 为实际匹配到的 BM2 pad 计算 SNR
  -> 按 `snr_timing.bin_edges` 对 BM2+三层 Tracker 的时间样本分组并解方差
`scripts/export_bm_snr_intermediates.py`
  -> 导出每组六个探测器对的 Δt 分布、核心筛选过程、Gaussian 诊断拟合和方差方程
```

输出在 `run/Combine/`：

```text
Tracker-step1-raw.root
Tracker-step2-ana.root
Tracker-step3-dec.root
Tracker-step4-rec.root
Pad-step3-rec.root
time_resolution/
├── BM2_Tracker1_gaussian_fit.png
├── BM2_Tracker2_gaussian_fit.png
├── BM2_Tracker3_gaussian_fit.png
├── bm2_tracker_time_differences.png/.pdf/.root
├── bm2_tracker_gaussian_fits.tsv
├── pairwise_time_resolution.tsv
├── detector_time_resolution.tsv
├── cfd_fraction_resolution.tsv/.png/.pdf
├── time_resolution.json
├── t30/、t40/、t50/（各比例的两两分布、拟合和时间分辨）
└── snr_binned/<execution_id>/
    ├── snr_bins.tsv
    ├── snr_timing_comparison.png/.pdf
    └── pairwise_details/（各 SNR 区间的图、逐事件表和方差解）
```

### Step4 对齐、几何接受区和残差统计

每次运行 signal 链都会用当前 run 的 `Tracker-step3-dec.root` 重新做三层联合对齐，不复用公共配置中的历史参数。对齐样本要求 Tracker1–3 的 X/Y 均恰好一个 cluster；Tracker1 固定为坐标基准，Minuit2 同时优化 Tracker2/3 的 `dx`、`dy` 和绕 Z 轴转角。结果写入：

- `Combine/MMalignmentPar3.dat`：仅平移的兼容格式；
- `Combine/MMalignmentPar6.dat`：实际使用的 `dx/dy/rotationZ`；
- `Tracker-step4-rec.root/alignment/TrackerAlignment`：参数、样本数、目标函数和 Minuit2 状态。

逐层留一效率要求两个参考 Tracker 都是唯一二维点，且参考点和预测 DUT 点分别满足 `|x|<50 mm`、`|y|<50 mm`；不施加径迹角度 cut。DUT 对应方向有信号且 `|residual|<10 mm` 时计入分子。`TrackerXYResidual` 进一步要求 DUT 也为唯一二维点。

每个 Tracker 的 X/Y、以及普通逐方向和严格 `exact_xy` 两种样本都会同时输出：

- 两次自适应核心高斯的 `mu`、`sigma`、误差、最终拟合范围、状态和 `chi2/ndf`；
- 固定 `[-2,+2] mm` 高斯的同类结果；
- 全分布 RMS、16%/84% 分位点和中心 68% 半宽。

结构化结果位于 `Tracker-step4-rec.root/TrackerResidualSummary`、`TrackerEfficiencySummary`、`Tracker-step4-residual-summary.tsv` 和 `Tracker-step4-efficiency-summary.tsv`。图像为 `Trackereff_position50_joint_alignment.png` 与 `TrackerXYResidual_position50_joint_alignment.png`。

## noise

`apv_events` 的每个 APV frame 取 128 通道、每通道 30 个采样点，适配为 `TrackRaw` 后完成稳健 pedestal/noise 拟合。输出全部位于 `run/Combine/noise/`：

```text
Tracker-noise-raw.root
fit_parameters.txt/.root
Tracker-noise-channel-parameters.txt
BM-noise-channel-parameters.txt
noise-global-parameters-1778.txt/.root
```

1778 = 3 × 2 × 256 个 Tracker strip + 2 × 121 个有效 pad。noise 标定用于噪声研究，其逐通道 sigma 还用于 Step2 显著峰最低高度；不参与当前预处理宇宙线信号的二次基线扣除或通道级 3σ 筛选。旧 `.dat` 解包入口已删除。

## 配置和可追溯性

每次执行在 `metadata/executions/<ID>/` 保存：

- `config/`：公共静态配置与 `run/config/` 覆盖合成后的快照；
- `signal.log`、`noise.log`、可选 `plots.log`；
- `manifest.json`：输入、输出、配置和代码哈希及各 stage 状态。

生成的 pedestal/noise 文件只存在各 run 的 `Combine/noise/`，不会复制进配置快照。`latest_signal.json`、`latest_noise.json` 只指向对应 stage 最近一次成功执行。

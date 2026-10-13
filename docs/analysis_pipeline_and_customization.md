# APV25/BTTracker 分析流程与二次开发指南

本文描述当前 APV25 变体的实际执行路径。这里的“0 压缩/zero-suppressed 模式”是宇宙线 signal 模式：输入已经由读出程序完成解码、基线扣除和 3σ 通道筛选。当前 C++ 框架不再接收 `.dat`，也不会在 signal 模式引用 noise pedestal mean 做第二次基线扣除；Step2 只引用本 run 的 noise sigma 设置显著峰最低高度。

## 1. 一眼看懂数据层级

| 层级 | 文件/树 | 一条记录代表什么 | 主要用途 |
|---|---|---|---|
| decoded ROOT | signal: `apv_EventSamples`；noise: `apv_events` | 一块 APV 的一个 frame | 当前唯一采集输入 |
| step1 | `Tracker-step1-raw.root/TrackRaw` | 一个物理通道的一条 30 点波形 | 统一硬件编号和事件号 |
| step2 | `Tracker-step2-ana.root/TrackX0..TrackY2`、`BMPadSignals` | 一个通过波形 cut 的 strip/pad 信号 | 通道级幅度和时间信息 |
| Tracker step3 | `Tracker-step3-dec.root/Track0..2` | 某一事例在某层 Tracker 的 X/Y 主簇和二维 hit | cluster/event 级 Tracker 信息 |
| Tracker step4 | `Tracker-step4-rec.root/Track` | 一个触发事例的三层径迹 | 直线拟合、残差、预测位置 |
| BM event | `Pad-step3-rec.root/BM1Events`、`BM2Events` | 各 BM 独立的事件表 | BM 空间、电荷和事件时间 |
| 径迹匹配时间 | `Pad-step3-rec.root/TrackArrivalEvents` | 一条有效三层径迹及其附近 BM/Tracker 时间 | BM 时间分辨研究 |

selection 只在 step2 做一次。step3 不再重复做 TOT、幅度或单峰 cut，只负责把 step2 已接受的通道组织成 cluster 和 event。

## 2. 电子学映射和几何

源 ROOT 的 `apv_board` 到物理探测器的关系为：

| 源 `apv_board` | 探测器 | step1 内部编号 |
|---:|---|---|
| 0、1 | Tracker1 X | board 0，chip 0、1 |
| 2、3 | Tracker1 Y | board 0，chip 2、3 |
| 4、5 | Tracker2 X | board 1，chip 0、1 |
| 6、7 | Tracker2 Y | board 1，chip 2、3 |
| 8 | BM1 | board -1，chip 12 |
| 9 | BM2 | board -1，chip 13 |
| 10、11 | Tracker3 X | board 2，chip 0、1 |
| 12、13 | Tracker3 Y | board 2，chip 2、3 |
| 14、15 | 空载 | 丢弃 |

注意：源 `apv_board=12/13` 是 Tracker3 Y；只有经过 step1 适配后的内部 `chip=12/13` 才是 BM1/BM2。

z 位置保持原 BTTracker 几何，只删除 Tracker4：Tracker1=533、Tracker2=468、BM1=364、BM2=309、Tracker3=217 mm。共享常量在 `include/CRInfo.h` 的 `Zpos` 和 `BMZpos`。

## 3. signal（zero-suppressed）流程

### 3.1 decoded ROOT → step1 通道波形

入口：`DecodedAPVRoot2TrackRaw`。

1. 验证输入存在 `apv_EventSamples`。
2. 读取每个 APV frame 的有效通道列表；只展开上游已经通过 3σ 的通道。
3. 每个有效通道写一条 `TrackRaw`：`event, board, readoutBoard, chip, channel, wave[30], timestamp`，同时保存 `source_file/source_entry/source_raw_event_index/source_apv_board/source_trigger` 以便追溯。
4. 16 位 trigger 回绕时扩展事件号；多个输入 ROOT 的事件号保持单调。
5. slot 14/15 计入 ignored 计数但不写出。

`DecodedInputSummary` 保存输入 frame、Tracker frame、BM frame、空载 frame、非法通道和输出波形数。`InputPreprocessing` 明确记录 signal 已完成基线扣除和 3σ 筛选。

### 3.2 step1 → step2 通道级波形参数

入口：`AnalysisAPV`。

Tracker 通道通过 `config/map_APV.txt` 映射到 0–255 strip；BM 通道通过 `config/map_APV4BMM.txt` 映射到 11×11 pad 位置。由于 signal 已预处理，框架仍固定使用 `ped=0, pedrms=0`，不重复做基线扣除；同时按 `global_channel` 从本 run 的 `noise-global-parameters-1778.txt` 读取 sigma。该 sigma 保存在独立的 `peak_noise_rms` 中，只用于计算 `peak_threshold_adc = peak_threshold_sigma × peak_noise_rms`，不会改变 `nisel`、sigmoid 拟合起点或后续 BM 电荷阈值。

每条 30 点波形计算：

- `maxA`、`maxT`：峰值和峰值采样位置；
- `riseT`：20% 到 80% 上升时间；
- `TOT`：50% 上升交点到 50% 下降交点，即离散波形的 FWHM 近似；
- `arrivalT_sigmoid`/`sigmoid_t20`：上升沿 sigmoid 拟合得到的 20% 时间；
- `sigmoid_t50/B/H/k/status/chi2/start/end`：拟合参数和质量；
- `peak_noise_rms`、`peak_threshold_adc`：本通道噪声 sigma 和实际显著峰 ADC 阈值；
- `n_significant_peaks`、`wave`、`timestamp`：峰数诊断、原波形和时间戳。

默认 step2 cut 是开区间：

```text
maxA > 0
0 < maxT < 15
0 < riseT < 10
2 < TOT < 14
peak_threshold_sigma = 5
require_single_peak = true
require_valid_sigmoid = true   # 只作为 Tracker 接受条件
```

BM 也应用幅度/时间和单峰条件；BM 的 sigmoid 拟合失败不会删除波形，而是保留 `status` 和无效时间值，便于单独研究效率。Tracker 默认要求有效 sigmoid 拟合。

输出树为 `TrackX0/Y0`、`TrackX1/Y1`、`TrackX2/Y2` 和 `BMPadSignals`。`WaveformPeakDiagnostics` 保留基本 cut 后、单峰拒绝前的诊断分母。ROOT 中的 `WaveformSelection` 元数据保存本次实际 cut。

### 3.3 step2 → Tracker step3 cluster/event

入口：`DecodeAPV`。

1. 按 `event + Tracker层 + X/Y` 收集 step2 strip。
2. 按 strip 编号排序；相邻记录的 strip 差大于 2 时分成新 cluster。
3. 以 cluster 总幅度最大的簇作为主簇。
4. 主簇位置用幅度加权；保存主簇总幅度、strip 数、`maxT`、时间均值/最早值/有效数、簇内局部峰数。
5. 其他 cluster 的位置/能量写入 `x_other/y_other`；总簇数写入 `cluster_count_x/y`。
6. 对六棵 Step2 Tracker 方向树的原始 `event` 取全局并集并排序；完全没有 Tracker Step2 信号的 DAQ event 不写入 Step3。
7. 并集中的每个 event 在 `Track0/1/2` 各写一条，缺测层/方向通过 `sig_x/sig_y/sig` 表示。三棵树的显式 `entry` branch 同步取值 `1..N`，原始 `event` 保持不变，继续作为 Step2、BM 和下游结果的关联键。ROOT 隐式 `Entry$` 仍为 `0..N-1`。

step3 不再调用波形 cut。这保证所有 cut 只有 step2/config 一个事实来源。

### 3.4 Tracker step3 → step4 径迹

入口：`ReconstructTrack`。

`Track0/1/2` 按显式 `entry` 和原始 `event` 逐行对齐；`entry` 是从 1 开始的紧凑重建序号，`event` 是不重编号的 DAQ trigger。每次 signal 分析都从当前 run 的 step3 收集“三层 X/Y 均恰好一个 cluster”的共同样本，固定 Tracker1 为坐标基准，用 Minuit2 联合优化 Tracker2/3 的 `dx`、`dy` 和 `rotationZ`。实际重建使用完整二维旋转变换，而不是只读取旧的三参数平移。参数写入 run 自己的 `Combine/MMalignmentPar3.dat`、`MMalignmentPar6.dat`，并在 `Tracker-step4-rec.root/alignment/TrackerAlignment` 保存样本数、目标函数变化和最小化状态。

对齐后对三层 X(z)、Y(z) 分别做直线拟合。`Track` 保存 `MMx/MMy/MMz`、每层有效标志、`XfitPars/YfitPars`、`useflag`、时间戳，以及供 BM 匹配使用的 `predStripX/predStripY`；预测 strip 会先对完整 `dx/dy/rotationZ` 变换做逆变换。

step4 的 `Trackereff` 留一效率不施加径迹角度 cut。对于每个 DUT，另外两层必须各自形成唯一二维 cluster，两个参考点都要满足 `|x|<50 mm`、`|y|<50 mm`，外推到 DUT 的预测点也要满足同样的 ±50 mm 几何接受区。DUT 对应方向存在主 cluster 且 `|residual|<10 mm` 时计入分子。`TrackerXYResidual` 进一步要求 DUT 也恰好有一个 X cluster 和一个 Y cluster，并要求两方向残差同时位于 ±10 mm 后成对填充，因此同一层 X/Y 图的 Entries 相同。独立直方图为 `resXY_X0/Y0`、`resXY_X1/Y1`、`resXY_X2/Y2`。

每个逐方向和严格 `exact_xy` 残差分布都同时计算：自适应两次核心高斯、固定 `[-2,+2] mm` 高斯、全分布 RMS、16%/84% 分位点和中心 68% 半宽。两种高斯均保存 `mu`、`sigma`、误差、拟合范围、状态和 `chi2/ndf`。结构化结果在 `TrackerResidualSummary`、`TrackerEfficiencySummary`、`Tracker-step4-residual-summary.tsv`、`Tracker-step4-efficiency-summary.tsv`；对应 PNG 为 `Trackereff_position50_joint_alignment.png` 和 `TrackerXYResidual_position50_joint_alignment.png`。

`useflag=true` 表示可用于后续径迹匹配的有效三层径迹。

### 3.5 BM 通道 → pad cluster 和 event 时间

入口：`ReconstructPadEvents`，同时读取 step2 `BMPadSignals` 和 step4 `Track`。

`BM1Events`、`BM2Events` 的处理是：

1. 按 `(event, raw_chip)` 分开处理 BM1（内部 chip12）和 BM2（内部 chip13），对各自 pad 的波形计算正/负极性中占优一侧的积分电荷；当前 pad 内部安全阈值是峰值大于 5 ADC（signal 中 `pedrms=0`，代码用最小 RMS=1）。`raw_chip` 必须存在，不能用输入 ROOT 的 slot 编号替代。
2. 同一 BM 内同一 pad 的重复记录累加电荷、保留最大峰值；两块 BM 的相同 `position` 独立统计。
3. 11×11 网格中共享边的 pad（曼哈顿距离 1）归为同一 cluster。
4. cluster 位置为电荷加权中心；主 cluster 是总电荷最大的簇。
5. 分别写 `BM1Events`、`BM2Events`，每棵树只保存本块 BM 有入选命中的事件。两棵树 Entries 独立，不为缺失 BM 写零值占位行。

两棵树使用相同的无 BM 前缀 branch 名称：

| branch（每棵 BM 树都有） | 含义 |
|---|---|
| `event` | 原始展宽后的 DAQ trigger 编号，用于与 Tracker/另一块 BM 按事件编号关联 |
| `timestamp` | 本块 BM 入选记录的最小有效时间戳，无有效时间戳为 `-1` |
| `nhits`, `nclusters` | 本块 BM 的不同命中 pad 数、cluster 数 |
| `hit_position`, `hit_x/y`, `hit_charge`, `hit_peak` | pad 编号、中心坐标（mm）、波形积分（ADC·sample）、峰值（ADC） |
| `cluster_size`, `cluster_x/y`, `cluster_charge` | 簇大小、电荷加权重心（mm）、总积分 |
| `primary_x/y`, `primary_charge` | 本块 BM 最大电荷簇的重心和总积分 |
| `arrivalT[_a/_e/_sigmoid]_mean/min`, `n_arrivalT[_a/_e/_sigmoid]` | 本块 BM 各方法的有效时间平均/最早值及记录数，时间单位为 sample |
| `arrivetime_mean/min` | 本块 BM 的 `arrivalT_mean/min` 别名 |

未接入 BM 时，在该 run 的 `config/configCR.json` 中显式设置：

```json
{
  "detector_layout": {
    "bm1_enabled": false,
    "bm2_enabled": true
  }
}
```

公共默认值为两块 BM 均启用。禁用 BM1 后，即使输入 slot8 有空载异常波形，`BM1Events` 和其中 `event` 等所有 branch 的 Entries 均为 `0`，不产生 `event=0` 的占位记录；BM1 的径迹匹配时间计数也为 `0`。此开关只作用于 pad 重建和 BM 匹配时间，保留 step1/step2 波形用于诊断，噪声标定仍覆盖全部 1778 个通道。`BMPadDetectorSelection` 保存本次 BM 启用状态。`timing_cuts.detectors` 仍单独决定时间分辨联合解的参与探测器。

原共享 `PadEvents` 不再写入。旧 ROOT 需重跑 signal 才能更新；正式 `Results/<tag>` 导出不修改。两棵 BM 树不能按 entry 下标直接关联，应使用 `event` 编号。

分别查看两块 BM 的最大簇电荷谱：

```cpp
BM1Events->Draw("primary_charge", "primary_charge>0");
BM2Events->Draw("primary_charge", "primary_charge>0");
```

这是未做能量标定的积分信号谱。两棵 BM 事件树使用各自全部入选 pad，不施加径迹匹配；径迹匹配时间保存在下面的 `TrackArrivalEvents` 中。

`TrackArrivalEvents` 只写 `useflag=true` 的径迹：

- BM：取径迹投影 pad 周围 3×3 窗口；`bmS_mean/min/n[2]` 是两块 BM 的 sigmoid t20；
- Tracker：在预测 strip ±2 范围内按 gap≤1 聚类；`trS_*[3][2]` 是三层 X/Y 时间；
- `evS_mean/min/n` 是所有选中 Tracker strip 合并后的 event 参考时间；
- A/E 是已停用 Pause 拟合的兼容分支，当前应使用后缀 `S` 的 sigmoid 量。

### 3.6 自动时间分辨与 BM2−Tracker 高斯图

每次成功的 signal 链在 PadReconstruction 后自动运行
`scripts/analyze_time_resolution.py`。BM2 事件时间是径迹中心 3×3 pad 窗口内
有效 sigmoid t20 的算术平均；每层 Tracker 时间是预测 strip ±2 范围内匹配 cluster
的有效 X/Y strip t20 按 strip 数合并后的平均。每一探测器对使用自己的共同有效事例，
时间差定义为 `Δt=t_first-t_second`，并乘 25 ns/sample。

默认主流程在构造所有时间差样本前先要求 `bmS_n[1] == timing_cuts.bm2_matched_pad_count`；当前默认值为 1，即 BM2 在径迹中心 3×3 窗口内恰好一个具有有效 sigmoid t20 的唯一 pad。这是匹配窗口的有效时间 pad 多重性，不是整块 BM2 的原始过阈通道多重性。设为 0 关闭此 cut。cut 在时间分析阶段统一筛选整个事例，因此 BM2–Tracker 与 Tracker–Tracker 分布来自同一单 pad 事例集合，不会删除 step2/step3 中间 ROOT 里的多 pad 信息。

两两分布先用 median/MAD 初始化并迭代做 3σ 核心选择；直方图显示核心中心 ±6σ，
高斯在核心中心 ±3σ 内拟合。正式输出位于 `Combine/time_resolution/`，包括
BM2−Tracker1/2/3 的单图和合并图、`bm2_tracker_gaussian_fits.tsv`、
`bm2_tracker_time_differences.root`，以及两两方差联合求解的
`pairwise_time_resolution.tsv`、`detector_time_resolution.tsv` 和
`time_resolution.json`。这些文件属于 signal 必需输出，写入 execution manifest，
导出时进入 `Results/<tag>/timing/time_resolution/`。

TSV 中 `gaussian_mu_error_ns` 和 `gaussian_sigma_error_ns` 是 ROOT 最小化器在最优点
由拟合参数协方差矩阵对角元开方得到的条件统计误差；应同时检查
`gaussian_status==0` 和 `covariance_status==3`。它们不包含 binning、拟合范围、
time-walk、通道刻度和事件相关性等系统误差。

## 4. noise 分析流程

noise 模式与 signal 模式共享 step1 适配器，但策略不同：

1. 输入必须是 decoded ROOT，并提供 `apv_events`；每个 APV frame 展开全部 128 通道，而不是只展开 3σ 通道。
2. 原始 ADC 不扣基线、不做阈值，写入 `Combine/noise/Tracker-noise-raw.root/TrackRaw`。
3. `BuildNoiseCalibration` 把每个通道的所有波形采样累计成分布。
4. 先用稳健统计给出拟合初值，再做受限高斯拟合；拟合失败时保留 fallback 统计和失败状态，而不是静默丢通道。
5. Tracker 1536 通道加 BM 242 个有效 pad，共 1778 个全局通道。

输出位于 `run/Combine/noise/`：

```text
Tracker-noise-raw.root
fit_parameters.txt/.root
Tracker-noise-channel-parameters.txt
BM-noise-channel-parameters.txt
noise-global-parameters-1778.txt/.root
```

`noise-global-parameters-1778.root` 的关键树为 `GlobalNoiseChannels`（mean/sigma、raw/robust 统计、fit status/χ²）和 `NoiseWaveforms`。这些结果不进入 signal 基线扣除；`noise-global-parameters-1778.txt` 的 sigma 按 global channel 进入 Step2 显著峰判定，默认阈值为 `5×sigma`。

## 5. 用 config 控制流程和 cut

公共默认值是 `BTTracker/config/configCR.json`。单个 run 的覆盖写在 `run/config/configCR.json`；批处理会递归合并两者，并把冻结快照保存到 `metadata/executions/<ID>/config/`。

例如只对一个 run 改成 `2<TOT<14`（当前默认本来就是这一范围）：

```json
{
  "signal_cuts": {
    "tot_min": 2.0,
    "tot_max": 14.0
  }
}
```

可配置的 step2 项为：`amplitude_min`、`max_t_min/max`、`rise_time_min/max`、`tot_min/max`、`require_single_peak`、`require_valid_sigmoid`。

阶段开关为：`normalize_decoded_root`、`analysis_apv`、`decode_apv`、`reconstruct_track`、`pad_reconstruction`、`noise_analysis`，以及独立统计/绘图 mode 的开关。

如果 step1 已有效，只想用新 cut 重跑 step2 以后，可以在 run 覆盖中设置：

```json
{
  "run": {
    "normalize_decoded_root": false,
    "analysis_apv": true,
    "decode_apv": true,
    "reconstruct_track": true,
    "pad_reconstruction": true
  },
  "signal_cuts": {
    "tot_min": 2.0,
    "tot_max": 14.0
  }
}
```

然后用 `force=1`。关闭某阶段并不会伪造输入；后续阶段会验证所需中间 ROOT 是否存在且包含正确树。

## 6. 增加某个 step 的 cut 或统计

先判断改动属于 selection 还是 measurement：

- selection 会改变哪些通道进入后续重建，例如 TOT cut。应放在 step2 `AnalysisAPV::Cut()`，参数放 `SignalSelectionConfig`，由 `CRPipeline.C` 从 config 读取；不要再在 step3 复制一次。
- measurement 只统计已有信息，例如 BM1−BM2 时间差。应新建只读模块读取固定 step 输出，写独立 ROOT/TSV/PNG，不要修改上游树或把统计 cut 混入重建。

增加一个新的 step2 cut 时，完整改动点是：

1. 在 `include/AnalysisAPV.h::SignalSelectionConfig` 增加字段和默认值。
2. 在 `src/CRPipeline.C::LoadSignalSelection` 从 `CRConfig` 读取并验证范围。
3. 在 `src/AnalysisAPV.C::Cut` 或相应选择位置应用。
4. 把实际值加入 `WaveformSelection` 元数据。
5. 更新 `config/configCR.json`、`config/configCR.md` 和本文件。
6. 用一个边界样例验证 cut 是 `<` 还是 `<=`，再用 `force=1` 重跑 step2 以后。

增加一个统计 mode 时，建议：

1. 新建 `include/BMTimeResolution.h`、`src/BMTimeResolution.C`，只读 `Pad-step3-rec.root/TrackArrivalEvents`。
2. 事件选择至少要求 `bmS_n[0]>0 && bmS_n[1]>0`；可再对径迹质量、pad cluster 大小或参考时间计数加 cut。
3. 统计 `Δt=bmS_mean[0]-bmS_mean[1]`（或用 `min` 定义），把选择条件和时间定义写入输出元数据。
4. 若假设两块 BM 独立且分辨相同，单块分辨可估为 `σ(Δt)/sqrt(2)`；不满足该假设时需要 Tracker/外部参考并做方差去卷积，不能直接把 Δt 宽度叫单块时间分辨。
5. 在 `CRAnalysis.h` 声明函数，在 `ReadTracker.C` 增加 mode 和 config 开关，在 `scripts/run_manager.py::MODES` 加同名 mode。
6. 输出建议放 `run/Combine/BM-time-resolution.root/.tsv`，并由 run manager 纳入 manifest，避免覆盖 step2–step4。

现有 `arrival` mode 已改为读取 Tracker 的 `arrivalT_sigmoid/sigmoidFitStatus`，适合检查 Tracker 首簇最早时间与簇平均时间；BM 事件时间应直接使用 `TrackArrivalEvents`。

## 7. 日常操作命令

创建一个 run（默认符号链接，`--copy` 才复制）：

```bash
cd /home/thoriatong/workplace/APV25/BTTracker
./init_run.sh ../batches/run_005 \
  --signal /path/to/cosmic-decoded.root \
  --noise /path/to/noise-decoded.root
```

分别运行或连续运行：

```bash
./run.sh ../batches/run_005 1 noise
./run.sh ../batches/run_005 1 signal
./run.sh ../batches/run_005 1 all
```

批量运行：

```bash
./run_batch.sh ../batches --mode all --force 1 --jobs 2
```

`force=0` 只有在输入 SHA-256、配置快照、代码指纹和输出哈希全部与上次成功执行一致时才复用；改过 cut 或代码后应使用 `force=1`。

常用检查：

```bash
rootls -t ../batches/run_005/Combine/Tracker-step2-ana.root
rootls -t ../batches/run_005/Combine/Tracker-step3-dec.root
rootls -t ../batches/run_005/Combine/Pad-step3-rec.root
```

执行日志和可追溯信息在 `run/metadata/executions/<ID>/`，包括 `manifest.json`、冻结配置和各 mode 日志。正式结果用 `export_results.sh` 导出，不要手工把生成的 noise 参数放回公共 `config/`。

# BM2 信噪比与时间分辨诊断

## 主流程接入

### 分组配置的来源

managed run 的 signal 阶段成功后，默认自动执行本诊断及逐区间中间量导出。配置位于
`BTTracker/config/configCR.json` 的 `snr_timing`；单 run 可在
`run/config/configCR.json` 覆盖：

```json
{
  "snr_timing": {
    "enabled": true,
    "bin_edges": [10, 30, 50, 70, 90, 120]
  }
}
```

当前公共配置的边界生成以下左闭右开区间：`[10,30)`、`[30,50)`、`[50,70)`、`[70,90)`、`[90,120)`。
单 run 的 `run/config/configCR.json` 若设置了 `snr_timing.bin_edges`，则以该 run 的边界为准。
run manager 按对象递归合并配置，数组整体替换；它向两个 SNR 脚本显式传入合成后的边界。
当 run 覆盖中仍有旧 `bin_edges` 时，只修改公共配置不会改变该 run 的分组。
本次实际配置以 `metadata/executions/<execution_id>/config/configCR.json` 为准。

独立运行这两个 Python 脚本时，它们不会自动读取公共或 run 配置：

| 入口 | 未传 `--snr-bin-edges` 时的行为 |
|---|---|
| `analyze_bm_snr_timing.py` | 使用脚本内边界 `5,15,25,40,60,90,130,+∞`。 |
| `export_bm_snr_intermediates.py` | 请求旧的四个区间 `[15,25)`、`[25,40)`、`[40,60)`、`[60,90)`。 |

因此独立运行时应像下方示例一样显式传边界；导出中间量时，边界还必须与源诊断的
`snr_bins.tsv` 对应，不能靠修改导出参数给旧结果重新分组。

### 阶段与输出

关闭时设 `enabled` 为 `false`。诊断固定使用 BM2 与 Tracker1–3 的四探测器方差解，
并要求 BM2 预测 3×3 pad 窗口中恰有一个有效时间 pad；这不会修改正式时间分析或 Step2 cut。
每组不足 40 个事件时仅保存事件表和 `status.json`，不计算时间分辨；其余组正常完成。

每次 execution 的结果单独写入
`run/Combine/time_resolution/snr_binned/<execution_id>/`，完整记录在 manifest，
所有请求阶段成功后，最终分组表、报告、`snr_threshold_scan.png/pdf` 自动归档到
`Results/<execution_id>/timing/time_resolution/snr_binned/`；各组逐探测器结果表
`detector_solution.tsv` 归到其中的 `detector_solutions/`。
逐事件表、ROOT、分布/机制诊断图和 `pairwise_details/` 的计算过程保留在 Combine。
也可通过 `export_results.sh` 使用自定义结果 tag。
下文的手动命令仍适用于单独诊断已有 run，不会改写主流程结果。

从 APV25 项目根目录执行（只需要 PyROOT 和 Python 标准库）：

```bash
cd /home/thoriatong/workplace/APV25
source /home/thoriatong/root/bin/thisroot.sh
python3 BTTracker/scripts/analyze_bm_snr_timing.py \
    batches/run_006 batches/run_008 batches/run_009 \
    --snr-bin-edges 10,30,50,70,90,120
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

- `snr_timing_comparison.png/pdf`：归一化幅度、各 SNR 区间的时间样本事例数、实际时间样本 SNR、各 SNR 组的 BM2 分辨。第二幅图的事例数按 `snr_bins.tsv` 中 `selected_events` 统计，表示进入 BM2 匹配单 pad 时间样本的事例数；每个区间为左闭右开 `[low, high)`。
  最后一幅图标题为 `time resolution in each SNR bin`，横坐标为区间内 SNR 中位数；数值来自 `weighted_core` 的
  `resolutions_ns` 和 `resolution_uncertainties_ns`；纵轴为分辨率 σ（ns），不是方差 σ²（ns²）。
- `snr_timing_mechanisms.png/pdf`：各 SNR 区间内上升沿宽度和 noise/20%斜率代理量的算术平均值，以及 BM2−Tracker2 的核心均值，随区间内平均 SNR 的变化。
- `snr_threshold_scan.png/pdf`：按配置的 SNR 区间绘制 BM2 分辨率及各区间样本占比；为兼容旧文件名保留 `threshold_scan` 名称。真正的累计下限扫描数值仍保存在 `snr_thresholds.tsv`。
- `snr_bins.tsv`、`snr_thresholds.tsv`、`snr_bands.tsv`：分组、下限 cut 和区间 cut 的详细结果。
- `matched_bm2_events.tsv`：逐事件关联到的 Step2 波形条目、pad、global_channel、幅度、噪声、SNR 和时间。
- `snr_timing_report.json`：定义、样本计数、分布统计和数值结果。
- `snr_timing_plots.root`：可继续查看的画布、直方图和曲线。

比较时区分：低 SNR 的随机定时抖动；幅度/形状相关的时间均值偏移；高幅区可能的峰顶
非线性或截断；以及气体、场、增益、漂移过程造成的真实时间差异。
SNR 分组相关性不能单独证明单面结构改变是因果来源，也不能保证 SNR cut 越高越好。
尤其高幅区出现均值移动或宽度回升时，需要查看波形后再判断是否饱和，不能只凭幅度谱的凸起定性。

## 导出 SNR 分组计算的中间数据

针对已经按当前边界生成的诊断目录，可以直接回放保存的逐事件数据，导出五个区间
`[10,30)`、`[30,50)`、`[50,70)`、`[70,90)`、`[90,120)` 的计算过程：

```bash
source /home/thoriatong/root/bin/thisroot.sh
cd /home/thoriatong/workplace/APV25
python3 BTTracker/scripts/export_bm_snr_intermediates.py \
    "batches/run_006/Combine/time_resolution/snr_binned/<执行编号>" \
    --snr-bin-edges 10,30,50,70,90,120
```

运行前将 `<执行编号>` 替换为真实目录名。上例要求源目录已经使用这五个区间；旧四区间
诊断应使用其原边界导出，或先重新运行诊断生成新边界的结果。导出器会核对原表，
有足够事件却没有对应区间参考行时会报错。

输入为该目录的 `matched_bm2_events.tsv` 和 `snr_bins.tsv`，不会重读后来可能改变的
生产数据。结果默认写入该诊断目录下新的 `pairwise_details_<UTC时间>/`；也可用
`--output-dir` 指定尚不存在的目录。顶层 `README.md` 提供所有 run/区间的图像和 ROOT 索引。
managed 单 run 分析会传入 `--single-run run_NNN`，使每个区间直接位于
`Combine/time_resolution/snr_binned/<execution_id>/pairwise_details/SNR_下限_上限/`；
独立处理多个 run 的诊断仍按 run 分目录，避免不同 run 的同名 SNR 区间相互覆盖。

每个 `SNR_下限_上限/` 子目录（多 run 诊断时为 `run_NNN/SNR_下限_上限/`）包含：

- 六个探测器对的独立 PNG/PDF：左图为包含尾部的完整分布，右图为核心附近与高斯拟合。
- `time_differences_6pairs.png/pdf`：六个时间差分布的汇总图。
- `time_differences.root`：全部样本、核心样本、完整范围直方图、拟合函数、画布与 `samples_*` tree。
- `selected_events.tsv` 和 `<探测器对>_events.tsv`：区间事件及每对实际使用的时间差；`in_core` 标明是否进入核心。
- `<探测器对>_clipping.tsv`：每轮迭代 3sigma 的范围、事件数、均值和标准差。
- `pairwise_statistics.tsv`：原始 RMS、核心标准差、高斯参数及误差、范围、拟合状态。
- `variance_equations.tsv`：六个方差方程的观测值、近似误差、权重、模型预测和残差。
- `detector_solution.tsv` 与 `solution.json`：四个探测器的结果、A/W 矩阵、正规矩阵、右端向量与方差协方差。

新增的高斯拟合只用于诊断，拟合原始时间差直方图的核心均值 ±3 倍核心标准差范围。
原表使用的是迭代 3sigma 后事件的样本标准差来求解方差，不使用这些新增高斯 sigma。
程序逐组核对原表的选中事件数、每对事件数与核心标准差，以及最终 BM2 sigma/误差。
`provenance.json` 保存输入文件 SHA256、ROOT 版本和核对状态。

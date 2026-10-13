# configCR.json

默认配置控制 C++ 分析阶段：

```json
{
  "run": {
    "normalize_decoded_root": true,
    "analysis_apv": true,
    "decode_apv": true,
    "reconstruct_track": true,
    "pad_reconstruction": true,
    "noise_analysis": true,
    "amplitude_distribution": false,
    "energy_spectrum": false,
    "arrival_time_comparison": false
  },
  "detector_layout": {
    "bm1_enabled": true,
    "bm2_enabled": true
  },
  "signal_cuts": {
    "amplitude_min": 0.0,
    "max_t_min": 0.0,
    "max_t_max": 15.0,
    "rise_time_min": 0.0,
    "rise_time_max": 10.0,
    "tot_min": 2.0,
    "tot_max": 14.0,
    "peak_threshold_sigma": 5.0,
    "require_single_peak": true,
    "require_valid_sigmoid": true
  },
  "timing_cuts": {
    "bm2_matched_pad_count": 1,
    "detectors": ["BM1", "BM2", "Tracker1", "Tracker2", "Tracker3"]
  }
}
```

字符串 `"True"`/`"False"` 和 JSON 布尔值都可识别。

`signal_cuts` 在 step2 统一生效，数值范围使用开区间。step3 不会重复做这些波形 cut。`peak_threshold_sigma` 控制显著峰候选的最低高度：`峰值 > peak_threshold_sigma × 本 run 对应 global_channel 的 noise sigma`。默认值为 5；noise pedestal mean 不参与信号基线扣除。实际生效值同时写入 `Tracker-step2-ana.root` 的 `WaveformSelection` 和 `PeakSelectionNoiseCalibration` 元数据。

每次运行会把公共默认配置复制到 `metadata/executions/<ID>/config/`，再应用 `run/config/`：

```json
{
  "run": {
    "pad_reconstruction": false
  }
}
```

上述 run 覆盖只改变一个键，其余保持公共默认值。执行期间 C++ 从冻结快照读取 `config/configCR.json` 和 mapping。

如果某个 run 的 BM 接线发生交换，可只在该 run 的 `config/configCR.json` 中覆盖输入 slot 到内部 BM chip 的映射。例如 run_005 使用：

```json
{
  "detector_layout": {
    "bm_slot8_logical_chip": 13,
    "bm_slot9_logical_chip": 12
  },
  "timing_cuts": {
    "detectors": ["BM2", "Tracker1", "Tracker2", "Tracker3"]
  }
}
```

默认值仍是 slot8→chip12（BM1）、slot9→chip13（BM2）。映射覆盖同时作用于 noise 和 signal 归一化；交换后 slot8 使用 BM2 的 pad map 与 z=309 mm，slot9 使用 BM1 的 pad map 与 z=364 mm。`timing_cuts.detectors` 控制方差联合解包含哪些实际接入的探测器；未接入的 BM1 应从列表中移除，但 BM2 必须保留。

`detector_layout.bm1_enabled` / `bm2_enabled` 默认均为 `true`，控制 Pad Step3 的物理 BM 启用状态。未接入 BM1 的 run 应设置 `bm1_enabled=false`，重跑 signal 后 `BM1Events` 及其 `event` 等所有 branch 的 Entries 均为 0，且不写占位 `event=0`。BM2 独立保存于 `BM2Events`。该开关也清空禁用 BM 的径迹匹配时间；step1/step2 波形和全 1778 通道噪声标定保留用于诊断。此开关与 `timing_cuts.detectors` 的方差联合解选择相互独立。

显式请求某 mode 但对应开关为 false 时返回错误，避免批处理把“未执行”记作成功。关闭中间阶段时，后续开启阶段要求该 run 已有可验证的中间 ROOT。

signal 的输入已预处理，仍保持 `pedestal=0`，不再次做基线扣除或通道 3σ 筛选；但 Step2 显著峰判定会读取本 run 的 `Combine/noise/noise-global-parameters-1778.txt`，仅使用其中的 sigma。noise mode 从 `Tracker/noise/` 读取原始波形，并把 1778 通道标定写到 `Combine/noise/`。因此首次处理应使用 `--mode all`，或先运行 noise 再运行 signal。

`timing_cuts.bm2_matched_pad_count` 控制 signal 结束后自动执行的时间分辨分析。默认值 `1` 要求 `TrackArrivalEvents.bmS_n[1] == 1`；这里统计的是径迹预测点周围 3×3 匹配窗口内、具有有效 sigmoid t20 的 BM2 唯一 pad 数，不是整块 BM2 平面上的原始过阈通道数。设为 `0` 可关闭该多重性 cut。该条件作用于整套时间差样本，使 BM2–Tracker 和用于方差联合解的 Tracker–Tracker 对使用同一个单 pad 事例集合。`timing_cuts.detectors` 默认包含 BM1、BM2 和三层 Tracker；缺少某个探测器数据时，应在 run 覆盖中只列出实际参与联合解的探测器。

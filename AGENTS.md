# 项目背景与约束

本项目分析宇宙线实验中的 Micromegas 数据，是 BTTracker 的三层 Tracker 变体。

## 探测器和输入

- 3 块 10 cm × 10 cm 的二维 strip Tracker；Tracker4 已移除。
- 每层 X、Y 各 256 条 strip，每个方向由两块 128 通道 APV25 读出。
- 另有两块 pad-readout Micromegas Detector，重点分析其时间及与 Tracker 径迹的关联。
- P1、P2 Cherenkov/PMT 符合作为宇宙线总触发。
- 几何沿用原 BTTracker：沿 z 依次为 Tracker1=533、Tracker2=468、BM1=364、BM2=309、Tracker3=217 mm；只移除原来位于 154 mm 的 Tracker4。

宇宙线输入是含 `apv_EventSamples`、`apv_EventHeaders` 的解码 ROOT，已完成基线扣除和 3σ 通道筛选。普通信号分析不得再次读取 noise pedestal mean 做基线扣除，也不得用 noise RMS 重做通道级 3σ 筛选；唯一例外是 Step2 显著峰判定读取本 run 的逐通道 noise sigma，以 `peak_threshold_sigma × sigma` 作为局部峰最低高度。

专门的噪声 ROOT 使用原始波形树 `apv_events`，所有 128 个 APV 通道参与标定。`baseline_tree` 一类统计汇总 ROOT 不是噪声原始输入。

## APV slot 映射

| 输入 `apv_board` slot | 物理用途 | `TrackRaw` 内部编号 |
|---:|---|---|
| 0、1 | Tracker1 X | board 0，chip 0、1 |
| 2、3 | Tracker1 Y | board 0，chip 2、3 |
| 4、5 | Tracker2 X | board 1，chip 0、1 |
| 6、7 | Tracker2 Y | board 1，chip 2、3 |
| 8 | DUT / BM1 | board -1，chip 12 |
| 9 | DUT / BM2 | board -1，chip 13 |
| 10、11 | Tracker3 X | board 2，chip 0、1 |
| 12、13 | Tracker3 Y | board 2，chip 2、3 |
| 14、15 | 空载 | 不写入 `TrackRaw` |

这里的输入 slot 是解包 ROOT 中的 `apv_board`；它与归一化后的 `TrackRaw.chip` 不是同一编号空间。channel 到 strip/pad 的映射以 `BTTracker/config/` 为准。

## run 是数据管理的最小单位

一组宇宙线和与其对应的 noise 数据共同属于一个 run：

```text
batches/run_NNN/
├── Tracker/physics/   # 本 run 宇宙线 ROOT
├── Tracker/noise/     # 本 run 原始噪声 ROOT
├── config/            # 可选 run 覆盖配置
├── Combine/           # 宇宙线中间结果
│   └── noise/         # 本 run 噪声标定
├── metadata/          # 导入记录、日志、配置快照、哈希与状态
└── Results/           # 不可覆盖的正式导出
```

不得递归搜集 `apv25/preprocess/`，不得自动猜测 cosmic/noise 配对，也不得把不同 run 的文件放到同一输入目录。使用 `BTTracker/init_run.sh` 显式导入；默认创建符号链接，不移动原始数据。

每次执行从公共静态配置和 `run/config/` 生成独立配置快照。噪声标定只能写入该 run 的 `Combine/noise/`，不得写回公共 `BTTracker/config/`。旧 2290 通道标定仅作历史参考，不进入新快照。

包含 signal 的 managed run 在本次所有请求阶段成功后，自动生成不可覆盖的
`Results/<execution_id>/`，仅归档最终时间分辨数值、误差、汇总图及配置/来源记录。
step1–4、Pad ROOT、残差/时间差分布、噪声标定和逐事件计算过程保存在 `Combine/`；
时间分析工作副本也保留在 Combine 供复用校验。不得把整个 Combine 复制进 Results。

## 当前分析链

```text
Tracker/physics/*.root
  -> DecodedAPVRoot2TrackRaw
  -> AnalysisAPV（pedestal=0；run noise sigma 只用于显著峰最低高度）
  -> DecodeAPV
  -> ReconstructTrack（三层）
  -> PadReconstruction

Tracker/noise/*.root
  -> DecodedAPVRoot2TrackRaw（保留原始采样）
  -> BuildNoiseCalibration
  -> Combine/noise/*1778*
```

`.dat` 路径仅用于历史兼容。新增代码和文档应围绕 run 目录、解码 ROOT、三层 Tracker 和 1778 个噪声通道设计。

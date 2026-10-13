# BTTracker 当前源码导览

## 入口与调度

- `init_run.sh` → `scripts/run_manager.py init`：创建/更新受管理 run，显式导入数据。
- `run.sh` → `scripts/run_manager.py run`：单 run 或直接 decoded ROOT 入口。
- `run_batch.sh` → `scripts/run_manager.py batch`：编译一次，并行调度多个 run，保存状态和哈希。
- `ReadTracker.C`：C++ mode 分发，失败返回非零。
- `src/CRPipeline.C`：signal 各阶段和输出树验证。
- `src/NoiseCalibration.C`：run-local noise 适配、拟合与 QA。

批量管理器从每次执行的配置快照目录启动 `ReadTracker`，因此 C++ 中的 `config/...` 总是指向这一执行的冻结配置，而不是可变的公共目录。

## ROOT 输入适配

`DecodedAPVRoot2TrackRaw` 支持两种明确策略：

1. signal：读取 `apv_EventSamples` 的 `apv_chnum`，只保留上游 3σ 筛选后的通道，写入元数据说明已扣基线；
2. noise：读取 `apv_events`，遍历 APV 的 128 个物理通道，保留原始 30 点采样，写入元数据说明未扣基线、未做阈值。

两者统一生成 `TrackRaw(event, board, chip, channel, wave, timestamp, source_*)`。16 位 trigger 会展宽，多文件事件号保持单调。

源 `apv_board` 的硬件顺序为：0/1=Tracker1 X、2/3=Tracker1 Y、4/5=Tracker2 X、6/7=Tracker2 Y、8=BM1、9=BM2、10/11=Tracker3 X、12/13=Tracker3 Y、14/15=空载。适配器将 Tracker3 的源 slot 10–13 归一化为 board 2 的 chip 0–3，将源 slot 8/9 改为内部 chip 12/13，并丢弃空载 slot 14/15。不要把源 slot 12/13 误认为 BM chip12/13。

几何同样继承原 BTTracker：`Zpos={533,468,217}` 是原 `{533,468,217,154}` 删除 Tracker4 后的前三项，`BMZpos={364,309}` 保持不变。Tracker 与 BM 的径迹投影均使用这组共享常量。

## 三层 signal 流程

`MMdetN=3`，每层四块 APV：chip0/1→X，chip2/3→Y。

- `AnalysisAPV`：signal 输入设置 pedestal/RMS 为 0；执行波形基本 cut、显著单峰选择和 sigmoid t20 拟合，写 6 个 Tracker 树与 `BMPadSignals`。
- `DecodeAPV`：先取六个 Tracker 方向的原始 event 并集并压缩为三棵树共享的显式 `entry=1..N`，再按 strip 相邻关系聚类、选择主 cluster，写 `Track0..2`；原始 `event` 不重编号，供 BM 和下游关联。
- `ReconstructTrack`：三层直线拟合、alignment、残差、效率和每层预测 strip，写 `Track`。
- `PadReconstruction`：按 BM 分别聚类并写 `BM1Events`、`BM2Events`；把 BM 3×3 窗口和 Tracker 预测 strip 的时间与径迹关联，写 `TrackArrivalEvents`。未接入 BM 可在 run 配置中关闭，其事件树保持零条目。

## noise 流程

`CalibrateTrackerNoise` 只读取 `run/Tracker/noise/`，输出到 `run/Combine/noise/`。`BuildNoiseCalibration` 对 1536 个 Tracker 通道和 242 个 pad 通道累计原始 ADC 样本，使用稳健初值和受限高斯拟合，失败时保存诊断状态和 fallback 统计。

`noise-global-parameters-1778.root` 包含 `GlobalNoiseChannels`、逐通道波形直方图、拟合曲线及 raw/fit 统计。它不写回公共配置，也不参与已预处理 signal 的二次扣除；对应 TXT 中的 sigma 会按 global channel 传给 Step2，只用于显著峰最低高度。

## 防止跨 run 污染

`RunDataInput.h` 只枚举指定输入目录第一层，并且只接受 decoded ROOT。`run_manager.py` 另外提供：

- run 级非阻塞文件锁；
- 一次构建和可控并行；
- 输入/配置/代码 fingerprint；
- 输出哈希验证后才允许 `force=0` 复用；
- 某 stage 失败后记录状态并继续独立 stage/其他 run。

## 绘图和导出

`debug/macros/plot_peak_timing.C` 从 ROOT 中动态发现 `PeakCount_all_*`，三层配置当前产生 14 组、28 张峰数图，并固定产生 4 张 sigmoid Δt 图。数量不再写死为旧四层的 34 张。

`scripts/export_run.py` 只导出 manifest 中记录且哈希仍匹配的文件，配置来自成功执行快照。旧的四层 Tracker 重跑/验证宏已经删除，当前诊断绘图入口为 `debug/macros/plot_peak_timing.sh`。

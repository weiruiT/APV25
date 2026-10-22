# BTTracker/APV25 分析工程

当前工程以 run 为单位管理实验数据：一组宇宙线数据和对应 noise 数据放在同一个 run，所有标定、日志和结果都留在该 run 内。现有 `apv25/preprocess/zerozip0917` 与 `baselinefolder` 没有被自动配对或移动。

## 建立一个 run

先建空目录，稍后再显式导入：

```bash
cd /home/thoriatong/workplace/APV25/BTTracker
./init_run.sh ../batches/run_001 --notes "实验条件说明"

./init_run.sh ../batches/run_001 \
  --signal /path/to/cosmic.root \
  --noise /path/to/raw-noise.root
```

也可以一次完成创建和导入。默认使用符号链接；需要复制原始数据时加 `--copy`。导入器拒绝覆盖同名输入，也拒绝在已有成功分析后追加采集文件。

标准目录：

```text
run_001/
├── Tracker/physics/       # 已扣基线、已做 3σ 的宇宙线 ROOT
├── Tracker/noise/         # 含 apv_events 的原始噪声 ROOT
├── config/                # 可选：本 run 的映射/开关覆盖
├── Combine/               # step1–step4、Pad
│   └── noise/             # 1778 通道噪声标定和 QA
├── metadata/              # 执行清单、日志、配置快照、哈希
└── Results/               # 自动归档的最终时间分辨结果（不可覆盖）
```

## 运行

单个 run 的宇宙线链：

```bash
./run.sh ../batches/run_001 1 signal
```

`run.sh` 会检查当前 Python 是否能导入 PyROOT；如果 Conda 的 `python3` 不兼容，它会在加载 ROOT 环境后自动改用 `/usr/bin/python3`。也可用 `BTTRACKER_PYTHON=/path/to/python` 显式指定。

批量处理 `batches/run_*`，每个 run 先做 noise、再做 signal：

```bash
./run_batch.sh ../batches --mode all --force 1 --jobs 2
```

父目录自动发现只选择由导入工具管理且已有采集文件的 run；空模板和旧的未管理目录会跳过。

先检查将处理哪些 run：

```bash
./run_batch.sh ../batches --mode all --dry-run
```

每次成功的 signal 分析都会自动读取 `Pad-step3-rec.root/TrackArrivalEvents`，在
`Combine/time_resolution/` 并行生成 sigmoid 恒比 20%、30%、40%、50%（t20/t30/t40/t50）
的 BM2−Tracker 时间差分布、拟合参数和逐探测器时间分辨。旧的 t20 文件仍保留在该目录，
t30、t40、t50 的完整结果分别位于同名子目录；`cfd_fraction_resolution.png/.tsv` 汇总各比例的
逐探测器时间分辨。全部正式结果作为 signal 输出记录进 manifest。
默认时间样本要求 `bmS_n[1] == 1`，即 BM2 在径迹预测点周围 3×3 匹配窗口内恰好一个
具有有效 sigmoid 拟合的 pad；四个恒比时刻来自同一波形拟合，样本选择保持一致。可用
`timing_cuts.bm2_matched_pad_count=0` 关闭。若某块探测器未接入，用 run 配置中的
`timing_cuts.detectors` 将其从逐探测器方差联合解中排除。

另外，`cluster_timing.enabled=true` 会在独立的 `cluster_timing` 阶段分析与径迹匹配的
cluster 时间，结果写入 `Combine/time_resolution/cluster/<execution_id>/`，不覆盖上面的
单 pad 结果。空间聚类仍使用同一 BM 内四邻接 pad；只要 cluster 的任一成员位于原来的
3×3 匹配窗口内，它就是候选。要求 BM2 恰好一个空间候选 cluster，且其中至少一个 pad
具有有效 sigmoid 拟合；多个候选（包括没有有效时间的候选）会被排除。计时使用整个 cluster 的
有效 pad，各 pad 等权平均，并另存最早时间；同一 pad 的重复记录使用 maxA 最大者。
允许多 pad cluster，不施加 `bmS_n[1]==1`。这仍是窗口关联，不额外施加 cluster 质心残差 cut。

`BM1Events/BM2Events` 新增 `hit_cluster_index`、`hit_time_sigmoid[30,40,50]` 及
`cluster_time_sigmoid[30,40,50]_{mean,min,n}`。`TrackArrivalEvents` 新增
`bmClusterS[30,40,50]_{mean,min,n}`、
`bm_matched_cluster_{count,index,size,dx,dy}`；数组下标 0/1 对应 BM1/BM2，时间单位仍是
25 ns 采样点，dx/dy 为 cluster 质心减去对齐后的预测位置（mm）。旧 event 时间和 `bmS_*`
含义不变。新增分支的数值单位仍是 25 ns 采样点。新模式的 `time_resolution.json.single_pad_comparison` 同时记录各自全样本与
共同事件样本的结果，便于区分时间定义和样本选择的影响。

对包含新增分支的 Pad ROOT，也可从 APV25 项目根目录单独执行：

```bash
cd /home/thoriatong/workplace/APV25
source /home/thoriatong/root/bin/thisroot.sh
/usr/bin/python3 BTTracker/scripts/analyze_time_resolution.py batches/run_009/Combine/Pad-step3-rec.root \
  --bm-time-source cluster --detectors BM2,Tracker1,Tracker2,Tracker3 \
  --output-dir batches/run_009/Combine/time_resolution/cluster/manual
```

旧 Pad ROOT 需要先重新执行 PadReconstruction；单独运行脚本不能从缺少成员时间的旧
cluster 分支恢复 cluster 时间。正式运行请使用 `run.sh`，以记录配置快照、哈希和阶段状态。
`--plots` 只用于额外的峰数诊断和历史兼容 Δt 图。某个 run 或某个 stage 失败不会阻止其他 run；
批次最终返回非零，并把汇总写到 `metadata/batch_executions/<ID>/summary.tsv`。

signal 完成后默认还会运行独立的 `snr_timing` 诊断阶段：按
合成配置中的 `snr_timing.bin_edges` 分组（公共配置、run 覆盖及独立脚本的区别见
[SNR 配置说明](../docs/BM_SNR_TIMING_DIAGNOSTICS.md#分组配置的来源)），
计算 BM2 与三层 Tracker 的逐探测器时间分辨，并导出每组六个探测器对的时间差分布和逐事件中间量。
结果位于 `Combine/time_resolution/snr_binned/<execution_id>/`，同样写入 execution manifest，
每组计算过程直接位于其 `pairwise_details/SNR_<下限>_<上限>/` 子目录。
成功分析后，仅最终数值和分辨率汇总图进入 `Results/<execution_id>/timing/time_resolution/snr_binned/`；
逐事件数据、时间差分布和拟合过程保留在 Combine。该诊断固定选择恰有一个
有效匹配 BM2 pad 的事件，不改变 Step2 或主时间分辨的 selection；低于 40 个事件的区间会标记为统计量不足。
可在 run-local `config/configCR.json` 设置 `snr_timing.enabled=false` 关闭。
`--force 0` 会按诊断配置、输入/正式分析指纹、诊断代码和产物哈希复用该 stage；
`--force 1` 则为当前执行重新生成。

`--force 0` 只复用输入文件、配置快照、代码指纹和输出哈希全部匹配的成功执行；存在无法验证的旧输出时会要求 `--force 1`。

## 输入格式和分析策略

宇宙线 ROOT 必须有 `apv_EventSamples`，并已在上游完成基线扣除和 3σ 通道筛选。信号分析不会引用 noise pedestal mean 再次扣基线，但 Step2 会读取本 run 的 1778 通道 noise sigma，并以默认 `5×sigma` 作为显著峰最低高度。噪声 ROOT 必须提供逐 APV 原始波形；当前 `apv_events` 格式会保留 128 个通道的 30 点采样。

三层 Tracker 的 slot 映射摘要见 [workflow1.md](workflow1.md)；各 step 的树、分支、cut、聚类和扩展方法见 [完整分析与二次开发指南](../docs/analysis_pipeline_and_customization.md)；批次和元数据规则见 [批量数据管理](../docs/batch_data_management.md)。

Pad Step3 中 `BM1Events`、`BM2Events` 分别只保存本块 BM 有入选命中的事件，两棵树的 Entries 独立。未接入 BM1 时，在该 run 的 `config/configCR.json` 中设置 `detector_layout.bm1_enabled=false`（BM2 对应 `bm2_enabled`），重跑 signal 后 `BM1Events` 及其所有 branch 的 Entries 为 0。该开关作用于 pad 重建和 BM 径迹匹配时间；step1/step2 的原始记录及 1778 通道噪声标定仍保留。各树内使用 `event`、`nhits`、`primary_charge` 等无 BM 前缀的字段。

Tracker Step3 的 `Track0/1/2` 只保存六个 Tracker 方向中至少一个有 Step2 信号的原始 event。三棵树共享从 1 开始的显式 `entry` branch；`event` 仍是原始展宽后的 DAQ trigger 编号，供 Step2、BM 和下游结果关联。ROOT 自带的 `Entry$` 始终从 0 开始，不等同于显式 `entry`。

## 结果导出

managed run 的 `signal` / `all` 在本次所有请求阶段成功后，自动生成
`Results/<execution_id>/`；`metadata/latest_results.json` 指向最近一次成功归档。
只运行 noise 或有阶段失败时不自动发布。导出失败会报告 `partial_failed`，
已成功的分析及 Combine 文件保留，可以单独重试导出，无需重跑物理分析。

```text
Results/<execution_id>/
├── README.md                    # 最终结果文件索引
├── timing/time_resolution/      # 时间分辨、误差、CFD/SNR 汇总
└── metadata/                    # 配置快照、来源 manifest、输入和结果哈希
```

最终结果包括 `detector_time_resolution.tsv`、`pairwise_time_resolution.tsv`、
`time_resolution.json` 和 `cfd_fraction_resolution.tsv/.png/.pdf`；开启 cluster 或 SNR
分析时，归档本次执行对应的最终数值与汇总图。SNR 的逐探测器结果表放在
`snr_binned/<source_execution_id>/detector_solutions/`。
step1–4、Pad ROOT、噪声标定、残差/时间差分布、高斯诊断图和逐事件计算过程只保存在
`Combine/`。时间分析的工作副本仍保留在 `Combine/time_resolution/`，用于缓存校验和诊断。

也可对最近成功的 signal 手动命名导出：

```bash
cd /home/thoriatong/workplace/APV25/BTTracker
./export_results.sh ../batches/run_001 sigmoid_t20_v1
```

导出器验证执行的输入列表与 SHA-256、配置快照及记录的中间输出哈希，再按最终结果白名单复制文件。
同名标签不会被覆盖，已有 Results 不会自动迁移。可加 `--execution-id <id>` 指定历史执行，
其输入和输出必须仍匹配 manifest。noise 参数作为依赖单独校验，不计入宇宙线 ROOT 文件列表。

## 构建

脚本默认使用 `BTTracker/build_apv25/`，一个批次只编译一次：

```bash
source /home/thoriatong/root/bin/thisroot.sh
cd /home/thoriatong/workplace/APV25/BTTracker
cmake -S . -B build_apv25
cmake --build build_apv25 -j4
```

可用 `BTTRACKER_BUILD_DIR` 和 `BTTRACKER_BUILD_JOBS` 覆盖构建位置与并行数。

## 独立分析入口

主时间分辨、角度分组、空间分格和 cluster 分布的选择条件不同，见
[分析样本与 cut 对照](../docs/ANALYSIS_SELECTIONS.md)。
[BM2 位置 cut 时间分析](../docs/BM2_POSITION_CUT.md) 和
[BM 簇数/簇大小分布](../docs/BM_CLUSTER_DISTRIBUTIONS.md)
提供完整命令及输出说明；这些独立工具的产物保存在对应 run 的 Combine 中。

运行状态索引位于 ../batches/run_status.tsv，每次导入或批处理后自动刷新；也可从 BTTracker/ 执行 python3 scripts/run_status.py ../batches 手动重建。项目根目录的 [README](../README.md) 提供统一导航。

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
└── Results/               # 正式不可覆盖结果
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
`Combine/time_resolution/` 生成 BM2−Tracker1/2/3 的 sigmoid t20 时间差高斯拟合图、
拟合参数表、逐探测器时间分辨和 ROOT 对象，并作为 signal 输出记录进 manifest。
默认时间样本要求 `bmS_n[1] == 1`，即 BM2 在径迹预测点周围 3×3 匹配窗口内恰好一个
具有有效 sigmoid t20 的 pad；可用 `timing_cuts.bm2_matched_pad_count=0` 关闭。若某块探测器未接入，用 run 配置中的 `timing_cuts.detectors` 将其从逐探测器方差联合解中排除。
`--plots` 只用于额外的峰数诊断和历史兼容 Δt 图。某个 run 或某个 stage 失败不会阻止其他 run；
批次最终返回非零，并把汇总写到 `metadata/batch_executions/<ID>/summary.tsv`。

`--force 0` 只复用输入文件、配置快照、代码指纹和输出哈希全部匹配的成功执行；存在无法验证的旧输出时会要求 `--force 1`。

## 输入格式和分析策略

宇宙线 ROOT 必须有 `apv_EventSamples`，并已在上游完成基线扣除和 3σ 通道筛选。信号分析不会引用 noise pedestal mean 再次扣基线，但 Step2 会读取本 run 的 1778 通道 noise sigma，并以默认 `5×sigma` 作为显著峰最低高度。噪声 ROOT 必须提供逐 APV 原始波形；当前 `apv_events` 格式会保留 128 个通道的 30 点采样。

三层 Tracker 的 slot 映射摘要见 [workflow1.md](workflow1.md)；各 step 的树、分支、cut、聚类和扩展方法见 [完整分析与二次开发指南](../docs/analysis_pipeline_and_customization.md)；批次和元数据规则见 [批量数据管理](../docs/batch_data_management.md)。

Pad Step3 中 `BM1Events`、`BM2Events` 分别只保存本块 BM 有入选命中的事件，两棵树的 Entries 独立。未接入 BM1 时，在该 run 的 `config/configCR.json` 中设置 `detector_layout.bm1_enabled=false`（BM2 对应 `bm2_enabled`），重跑 signal 后 `BM1Events` 及其所有 branch 的 Entries 为 0。该开关作用于 pad 重建和 BM 径迹匹配时间；step1/step2 的原始记录及 1778 通道噪声标定仍保留。各树内使用 `event`、`nhits`、`primary_charge` 等无 BM 前缀的字段。

Tracker Step3 的 `Track0/1/2` 只保存六个 Tracker 方向中至少一个有 Step2 信号的原始 event。三棵树共享从 1 开始的显式 `entry` branch；`event` 仍是原始展宽后的 DAQ trigger 编号，供 Step2、BM 和下游结果关联。ROOT 自带的 `Entry$` 始终从 0 开始，不等同于显式 `entry`。

## 结果导出

```bash
./export_results.sh ../batches/run_001 sigmoid_t20_v1
```

导出器验证最近成功执行的输入列表与 SHA-256、配置快照、所有 ROOT 输出哈希，再创建 `Results/sigmoid_t20_v1/`。同名标签不会被覆盖。新 signal 执行把 run 内 noise 参数作为输入依赖并记录哈希，因此必须先完成 noise 标定。

## 构建

脚本默认使用 `BTTracker/build_apv25/`，一个批次只编译一次：

```bash
source /home/thoriatong/root/bin/thisroot.sh
cmake -S . -B build_apv25
cmake --build build_apv25 -j4
```

可用 `BTTRACKER_BUILD_DIR` 和 `BTTRACKER_BUILD_JOBS` 覆盖构建位置与并行数。

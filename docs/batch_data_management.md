# run 与批量数据管理

## 目录约定

一组实验条件下相互对应的宇宙线和 noise 数据称为一个 run。代码只维护一份，数据按 run 隔离：

```text
/home/thoriatong/workplace/APV25/
├── BTTracker/                 # 唯一代码工程
├── batches/
│   ├── README.md
│   ├── run_001/
│   │   ├── Tracker/physics/
│   │   ├── Tracker/noise/
│   │   ├── config/
│   │   ├── Combine/noise/
│   │   ├── metadata/
│   │   └── Results/
│   └── run_002/
└── apv25/preprocess/          # 现有数据来源，保持原位置
```

程序只枚举指定 run 输入目录的第一层，不再递归扫描 `apv25/preprocess`。采集输入只接受 decoded `.root`；旧 `.dat` 解包不再属于本框架。

## 创建与导入

创建空模板不会建立任何科学配对：

```bash
cd /home/thoriatong/workplace/APV25/BTTracker
./init_run.sh ../batches/run_001 --notes "gas=..., HV=..., trigger=..."
```

确认数据属于同一次 run 后显式导入：

```bash
./init_run.sh ../batches/run_001 \
  --signal /data/cosmic_part1.root /data/cosmic_part2.root \
  --noise /data/noise.root
```

可多次使用 `--signal`/`--noise` 追加不同文件，或传入一个目录导入其中第一层的 decoded ROOT。默认建立绝对符号链接，原始文件不移动；`--copy` 才复制。导入记录写入 `metadata/run.json`。

已有成功分析、同名目标或非本工具创建的现存目录都会被拒绝，避免静默改变 run 的数据定义。

`apv25/preprocess/baselinefolder/.../run4.root` 的 `apv_events` 可作为原始噪声候选；小型 `baseline.root` 只有 `baseline_tree` 汇总，不能作为噪声输入。现有 zerozip 与 baseline 数据尚未确认对应关系，因此项目没有自动生成 `run_001`。

## 单 run 与批量运行

```bash
# 单 run，只跑宇宙线
./run.sh ../batches/run_001 1 signal

# 单 run，只跑噪声
./run.sh ../batches/run_001 1 noise

# 所有 run：noise + signal
./run_batch.sh ../batches --mode all --force 1 --jobs 2

# 指定若干 run
./run_batch.sh ../batches/run_001 ../batches/run_003 \
  --mode signal --force 1 --jobs 2 --plots
```

可用 `--pattern 'run_2026*'` 控制父目录发现规则。父目录发现只选择带 `metadata/run.json` 且已有采集文件的 run，空模板和未管理目录会跳过；显式传入某个 run 时仍会严格检查。`--dry-run` 仅列出 run 和输入数量，不编译、不分析。

`all` 对每个 run 先执行 noise，再执行 signal。signal 的显著峰阈值依赖本 run 的 `Combine/noise/noise-global-parameters-1778.txt`；noise 缺失或失败时 signal 会因缺少该输入而失败，批处理整体返回非零。其他 run 不受影响。

同一 run 使用文件锁，不能被两个进程同时处理。构建也有独立锁，多个批次启动器不会同时改写构建目录。

## force 和污染防护

`--force 1` 重新生成当前 stage 输出。`--force 0` 只有在以下项目全部一致时复用：

- 输入文件列表、大小、路径和 SHA-256；
- 合成后的配置快照；
- 当前核心源码指纹；
- 上一次成功输出的 SHA-256。

存在输出但找不到匹配成功记录时拒绝复用。输入目录中的分析产物文件名也会被拒绝。

## 配置

公共 `BTTracker/config/` 提供静态 mapping 和默认 `configCR.json`。可在 `run/config/` 放覆盖文件：

- `configCR.json` 按 JSON 对象递归合并，可只写需要变化的键；
- mapping 等其他文件按相对路径替换。

生成的 pedestal/noise 文件不能放入 `run/config/`，它们属于 `Combine/noise/`。每次执行快照保存于 `metadata/executions/<ID>/config/`。

机械位置不变时可用 `--alignment-source /path/to/Combine` 为尚无 alignment 的 run 复制已有参数；目标已存在时不覆盖。

## 日志、汇总和发布

每个 run 的执行日志和 `manifest.json` 位于：

```text
metadata/executions/<execution_id>/
```

批次汇总位于第一个父目录参数下：

```text
metadata/batch_executions/<batch_id>/summary.tsv
metadata/batch_executions/<batch_id>/summary.json
```

正式发布：

```bash
./export_results.sh ../batches/run_001 sigmoid_t20_v1
```

导出器重新验证输入、配置和输出哈希，复制最近成功的 signal，若有独立成功的 noise 也一并包含。绘图只有在执行时加过 `--plots` 且被该 manifest 记录时才导出。同名 `Results/<tag>` 拒绝覆盖。

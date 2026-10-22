# APV25 三层 Tracker 分析

本项目处理 Micromegas 宇宙线数据：三层 strip Tracker 与 BM1/BM2 pad 探测器。数据按 run 隔离；代码、配置、运行记录和正式结果各有固定入口。

## 从哪里开始

- [BTTracker 使用指南](BTTracker/README.md)：创建 run、导入数据、运行分析和导出结果。
- [文档索引](docs/README.md)：分析流程、数据管理、SNR 诊断和宏说明。
- [batches 说明](batches/README.md)：run 目录结构及运行状态索引。
- [debug 说明](debug/README.md)：临时诊断脚本和历史参考。

## 常用操作

在 `BTTracker/` 目录运行：

```bash
./init_run.sh ../batches/run_001 --signal /data/cosmic.root --noise /data/noise.root
./run_batch.sh ../batches --mode all --force 1 --jobs 2
python3 scripts/run_status.py ../batches
```

采集数据和分析产物留在各自 run 内；备份代码时不纳入实验数据、`Combine/` 或 `Results/`。

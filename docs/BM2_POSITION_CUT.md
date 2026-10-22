# BM2 预测位置 cut 后的时间分辨

入口为 `BTTracker/scripts/analyze_bm2_position_cut.py`。它读取已有的 `Combine/Pad-step3-rec.root/TrackArrivalEvents`，按 event 关联 BM2 预测位置和时间，不需要重跑 Step1–4。

## 运行

从 APV25 项目根目录执行，使用能导入 PyROOT 的 Python：

```bash
cd /home/thoriatong/workplace/APV25
source /home/thoriatong/root/bin/thisroot.sh
/usr/bin/python3 BTTracker/scripts/analyze_bm2_position_cut.py \
  --runs 006 008 009 010 011
```

只分析 run_006：

```bash
/usr/bin/python3 BTTracker/scripts/analyze_bm2_position_cut.py --runs 006
```

`--runs` 接受保留前导零的编号，如 `006`，不是 `run_006` 或目录路径；省略时默认处理上述五个 run。`--workspace` 可指定另一份 APV25 根目录；`--execution-id` 可指定新的输出编号，已有同名目录会拒绝覆盖。缺少文件、重复 event 或样本不足导致异常时，当前调用会停止；已经完成的 run 输出保留。

## 样本和计算

- 固定使用单 pad 时间：`bmS_n[1]==1`。
- 位置来自 `bm_pred_x[1]`、`bm_pred_y[1]`，两个坐标都必须有限，并满足闭区间 `x∈[-10,40] mm、y∈[-40,10] mm`。
- Tracker 时间允许 X/Y 至少一个方向有效；每个时间差仅要求其两端有效，不强制三层 Tracker 同时有效或三层 X/Y 全部有空间响应。
- 固定使用 BM2、Tracker1、Tracker2、Tracker3，逐一计算 t20/t30/t40/t50 的平均时间、六对核心宽度和 `weighted_core` 方差解。
- BM2−Tracker1/2/3 的直方图范围、核心选择和高斯拟合复用主时间分析函数。BM2 分辨来自核心方差解，高斯拟合用于展示时间差分布。

矩形范围目前写在脚本顶部的 `X_MIN_MM/X_MAX_MM/Y_MIN_MM/Y_MAX_MM` 常量中，没有命令行范围参数。该工具不读取 `timing_cuts` 来改变其单 pad 条件或探测器列表，也不修改主时间分辨样本。其他入口的区别见 [cut 对照](ANALYSIS_SELECTIONS.md)。

## 输出

每个 run 分别生成：

```text
Combine/time_resolution/position_cut/<execution_id>/
├── README.md
├── position_cut_time_resolution.json
├── BM2_Tracker1_gaussian_fit.png
├── BM2_Tracker2_gaussian_fit.png
├── BM2_Tracker3_gaussian_fit.png
├── bm2_tracker_time_differences.png/.pdf/.root
├── bm2_tracker_gaussian_fits.tsv
└── t30/、t40/、t50/       # 各方法的高斯图、ROOT 和拟合表
```

顶层图片对应 t20。JSON 的 `methods.t20.weighted_core.resolutions_ns[0]` 是 BM2 分辨，`resolution_uncertainties_ns[0]` 是其近似统计误差；其他 CFD 存在对应方法键下。`events_after_position_cut` 是通过单 pad 与位置 cut 的总事件数，具体每对有效事件数和核心事件数见 `methods.<method>.pairwise`。

这是独立诊断入口，未接入 run manager 的自动阶段、manifest 和 Results 白名单；结果保存在 Combine。

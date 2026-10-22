# 当前宏与诊断工具

主分析、noise 标定和批量调度统一使用 `BTTracker/run.sh`、`run_batch.sh`，不需要手工调用 ROOT 宏。BM2−Tracker1/2/3 的正式时间分辨和按 SNR 分组的时间分辨/中间分布都在 signal 主流程自动生成，分别位于 `Combine/time_resolution/` 和 `Combine/time_resolution/snr_binned/<execution_id>/`。

## 额外诊断绘图入口

时间分辨率随径迹角度的分组图可用 `BTTracker/macros/PlotTimingVsAngle.C` 读取现有 Step4 和 Pad ROOT 生成，详见 [角度分组绘图命令与样本定义](TIMING_VS_ANGLE.md)。

BM 簇数和簇大小分布可用 `BTTracker/macros/PlotBMClusterDistributions.C` 绘制，命令、矩形位置 cut 及按事件去重的 size 计数规则见 [BM cluster 分布](BM_CLUSTER_DISTRIBUTIONS.md)；输出写入对应 run 的 `Combine/cluster_diagnostics/<execution_id>/`。

加入 BM2 预测位置矩形 cut 后重新计算时间分辨和 BM2−Tracker1/2/3 高斯图，使用
`BTTracker/scripts/analyze_bm2_position_cut.py`，见 [位置 cut 分析](BM2_POSITION_CUT.md)。
所有入口的样本条件集中列在 [分析样本与 cut 对照](ANALYSIS_SELECTIONS.md)。

批量运行时直接加：

```bash
cd /home/thoriatong/workplace/APV25/BTTracker
./run_batch.sh ../batches --mode all --force 1 --jobs 2 --plots
```

对已有成功输出单独绘图：

```bash
source /home/thoriatong/root/bin/thisroot.sh
cd /home/thoriatong/workplace/APV25/BTTracker
bash ../debug/macros/plot_peak_timing.sh \
  /home/thoriatong/workplace/APV25/batches/run_001/Combine
```

`plot_peak_timing.C` 会动态发现当前 ROOT 中的峰数分组。三层 Tracker 当前为 14 组、28 张 PNG；同时在 `Combine/sigmoid_only/` 生成 4 张历史兼容的 BM−Tracker2 sigmoid Δt 图、`delta_histograms.root`、`delta_summary.tsv` 和 `delta_values.tsv`。

`--plots` 只控制这些额外诊断。主流程的 `Combine/time_resolution/` 不需要该参数，
会随 signal 自动生成并记录在 manifest；其中最终数值与分辨率汇总图进入正式导出，
时间差分布及高斯诊断图保留在 Combine。

SNR 分组由合成配置的 `snr_timing.bin_edges` 决定。公共边界、run 覆盖示例及独立脚本的
显式传参命令统一见 [SNR 配置说明](BM_SNR_TIMING_DIAGNOSTICS.md#分组配置的来源)。

区间是相邻边界形成的左闭右开范围。该诊断固定要求 BM2 预测 pad 窗口中恰有一个有效 pad，
不改变正式时间分辨的 cut。完成后先查看 `snr_bins.tsv`，再从
`pairwise_details/README.md` 打开每区间的六组时间差图；每个 execution 使用独立目录，
避免覆盖先前诊断结果。低于 40 个匹配事件的区间只写状态和事件表。

## 噪声 QA

当前噪声 QA 文件位于：

```text
run_001/Combine/noise/noise-global-parameters-1778.root
```

通道范围为 1–1778：Tracker 1–1536，两个 pad Detector 共 242 通道。旧四层的 2290 通道结果已从公共 `config/` 清理。

## 参考宏

旧四层 Tracker 的重跑/验证宏和冻结输出已删除。`time_resolution.C`、`debug_alignment_sel.C` 等仍有算法参考价值，但不属于正式批处理入口；复用时必须改为 `MMdetN=3`、run-local 路径和当前 ROOT 输入。调试输出应写入 run 的临时诊断子目录，不得修改 `Tracker/` 输入和 `Results/` 正式结果。

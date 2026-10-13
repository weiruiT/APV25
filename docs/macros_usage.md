# 当前宏与诊断工具

主分析、noise 标定和批量调度统一使用 `BTTracker/run.sh`、`run_batch.sh`，不需要手工调用 ROOT 宏。BM2−Tracker1/2/3 的时间差高斯拟合和逐探测器时间分辨已经是 signal 主流程的自动步骤，输出在 `Combine/time_resolution/`。

## 额外诊断绘图入口

批量运行时直接加：

```bash
cd /home/thoriatong/workplace/APV25/BTTracker
./run_batch.sh ../batches --mode all --force 1 --jobs 2 --plots
```

对已有成功输出单独绘图：

```bash
source /home/thoriatong/root/bin/thisroot.sh
bash ../debug/macros/plot_peak_timing.sh \
  /home/thoriatong/workplace/APV25/batches/run_001/Combine
```

`plot_peak_timing.C` 会动态发现当前 ROOT 中的峰数分组。三层 Tracker 当前为 14 组、28 张 PNG；同时在 `Combine/sigmoid_only/` 生成 4 张历史兼容的 BM−Tracker2 sigmoid Δt 图、`delta_histograms.root`、`delta_summary.tsv` 和 `delta_values.tsv`。

`--plots` 只控制这些额外诊断。主流程的 `Combine/time_resolution/` 不需要该参数，
会随 signal 自动生成、记录在 manifest，并进入正式导出。

## 噪声 QA

当前噪声 QA 文件位于：

```text
run_001/Combine/noise/noise-global-parameters-1778.root
```

通道范围为 1–1778：Tracker 1–1536，两个 pad Detector 共 242 通道。旧四层的 2290 通道结果已从公共 `config/` 清理。

## 参考宏

旧四层 Tracker 的重跑/验证宏和冻结输出已删除。`time_resolution.C`、`debug_alignment_sel.C` 等仍有算法参考价值，但不属于正式批处理入口；复用时必须改为 `MMdetN=3`、run-local 路径和当前 ROOT 输入。调试输出应写入 run 的临时诊断子目录，不得修改 `Tracker/` 输入和 `Results/` 正式结果。

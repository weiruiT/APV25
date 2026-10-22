# debug 目录

这里保存临时诊断宏、一次性检查和历史报告，不是主分析入口。正式流程调用的分析脚本统一放在 `BTTracker/scripts/`。

- `macros/plot_peak_timing.sh` 与 `plot_peak_timing.C`：当前支持的 run-local 峰数和 sigmoid Δt 绘图；
- `BTTracker/scripts/analyze_bm_snr_timing.py`：正式 SNR 时间分析入口，由 run manager 调用；
- `BTTracker/scripts/export_bm_snr_intermediates.py`：正式 SNR 中间量导出入口，由 run manager 调用；
- `macros/PlotBM2SpatialTiming.C` 与 `plot_bm2_spatial_timing.py`：按 BM2 平面位置统计事例数和 t20 分辨率，并绘制匹配 pad 的 `maxA` 与 Δt20；用法见 `docs/BM2_SPATIAL_TIMING.md`；
- 旧的四层 Tracker 重跑/验证宏和冻结报告已删除；仍有参考价值的算法宏只作为源码资料保留；
- `logs/`：临时日志。

生产处理使用 `BTTracker/run_batch.sh`。调试输出不得写入 `Tracker/physics`、`Tracker/noise` 或覆盖 `Results/<tag>`。

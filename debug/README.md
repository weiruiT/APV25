# debug 目录

这里保存诊断宏、一次性检查和历史报告，不是主分析入口。

- `macros/plot_peak_timing.sh` 与 `plot_peak_timing.C`：当前支持的 run-local 峰数和 sigmoid Δt 绘图；
- `macros/analyze_bm_snr_timing.py`：读取指定 run 的 Step2、Pad Step3，精确关联 BM2 单匹配 pad，绘制 SNR 与时间分辨的比较、分组和 cut 扫描；见 `docs/BM_SNR_TIMING_DIAGNOSTICS.md`；
- 旧的四层 Tracker 重跑/验证宏和冻结报告已删除；仍有参考价值的算法宏只作为源码资料保留；
- `logs/`：临时日志。

生产处理使用 `BTTracker/run_batch.sh`。调试输出不得写入 `Tracker/physics`、`Tracker/noise` 或覆盖 `Results/<tag>`。

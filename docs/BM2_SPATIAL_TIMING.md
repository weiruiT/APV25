# BM2 空间分区时间分析

本工具按整个 100 mm × 100 mm 区域分格；固定矩形 ROI 后的整体时间分辨见
[BM2 位置 cut 分析](BM2_POSITION_CUT.md)，各入口选择条件见 [cut 对照](ANALYSIS_SELECTIONS.md)。

从 APV25 项目根目录运行：

```bash
source /home/thoriatong/root/bin/thisroot.sh
root -l -b -q 'debug/macros/PlotBM2SpatialTiming.C("batches/run_010")'
```

默认格子边长为 20 mm、每格至少要求 40 个事例，`maxA` 轴画到 2000 ADC。可调整统计门槛、幅度上限和格子边长：

```bash
root -l -b -q 'debug/macros/PlotBM2SpatialTiming.C("batches/run_010",20,2000,10)'
```

宏参数依次为 run 路径、每格最低事例数、`maxA` 轴上限（ADC）和格子边长（mm）。格子边长必须是正数，并能整除 BM2 的 100 mm 边长，例如 5、10、20 或 25 mm；默认 20 mm 对应 5×5，10 mm 对应 10×10。分格越细，每格统计越少；最低事例数只是结果有效性门槛，不是对整个样本施加的物理筛选。

分析把 BM2 平面视作以中心为原点的 `[-50,50) mm × [-50,50) mm` 区域，按参数给出的边长等宽分格。每个入选事件根据 `TrackArrivalEvents.bm_pred_x[1]`、`bm_pred_y[1]` 落入一个格子。采用轨迹外推到 BM2 的位置，而非击中 pad 中心。

事例先要求 BM2 匹配窗口恰有一个有效 sigmoid pad，再将 Step2 中该事件、该 BM2 pad 的 `maxA` 与 Pad Step3 时间关联；脚本核对该 pad 的 t20 与写入的 BM2 事件时间一致。地图及振幅图共同使用 BM2、Tracker1、Tracker2、Tracker3 在四种 CFD 下均有有效时间的事件。单格的时间分辨率继续用六组探测器时间差联合求解 BM2 方差，核心样本不足或解无效的格子标明原因，不报告零分辨率。

会生成三张 PNG/PDF：

- `bm2_event_count_map`：空间事例数热度图，颜色表示事例数，每格中央标出 BM2 的 t20 时间分辨率（ns）；统计不足的格子显示 `N/A`。所有热度图使用低值深蓝、再经青色和绿色过渡到高值黄色的配色，不使用红色。
- `bm2_resolution_map_t20`：每格 BM2 t20 分辨率热度图，标注分辨率误差和 N；统计不足的格子显示为空白并标出 N。
- `dt20_vs_maxA`：二维事例数密度图，横轴为实际匹配 pad 的 Step2 `maxA`（ADC），纵轴为 `mean_t20(BM2)-mean_t20(Tracker2)`（ns）。为突出主体，绘图区使用 Δt 分布的 0.5%–99.5% 分位范围；图中标注绘图区内和总事例数，被裁出的尾部事件仍保留在 `selected_events.tsv` 并计入 `analysis.json` 的范围统计。

输出目录为 `run_NNN/Combine/time_resolution/spatial_binned/<执行编号>/`，包含图像、ROOT 文件、`spatial_resolution.tsv`、逐事件 `selected_events.tsv` 及样本/范围说明 `analysis.json`。默认门槛不是物理 cut；应检查每格事例数与分布，必要时调整分组、门槛或坐标显示范围。

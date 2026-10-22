# 分析样本与 cut 对照

本文描述当前入口的实际选择条件。相同 run 中不同图的事例数可能不同，比较结果时需要同时确认输入 ROOT 的生成版本和分析入口。

## 空间径迹与有效时间

Step4 的 `useflag=true` 表示 X(z)、Y(z) 两个投影都通过 `fit_track`。三层配置下每个投影至少需要两层测量点，并通过代码中的拟合质量要求；这不等于三层 X/Y 六个方向全部有信号，也不等于三个 Tracker 都有有效时间。

`TrackArrivalEvents` 保存 `useflag=true` 的径迹及匹配时间。时间分析读取器把每个 Tracker 的有效 X/Y 时间按有效 strip 数合并；X 或 Y 至少一个方向有效即可形成该 Tracker 时间。空间信息用于重建径迹，时间有效性是另外一层选择。

## 各入口的条件

表中“矩形 ROI”指预测位置的闭区间 `x∈[-10,40] mm、y∈[-40,10] mm`，不是实测 pad 中心的位置。

| 入口 | BM 条件 | Tracker 时间/空间条件 | 额外位置或分组条件 |
|---|---|---|---|
| 主时间分辨 `analyze_time_resolution.py`，pad 模式 | managed run 默认 `bmS_n[1]==1`，可由 `timing_cuts.bm2_matched_pad_count` 改变 | 每个探测器对分别要求双方时间有效；不强制全部 Tracker 同时有效，也不额外要求三层 X/Y 全有空间信号 | 无额外矩形 ROI |
| cluster 时间分辨，`--bm-time-source cluster` | BM2 唯一匹配空间簇，其中至少一个 pad 时间有效；允许多 pad | 沿用逐对有效时间；参与方差解的探测器由 `timing_cuts.detectors` 指定 | 无额外矩形 ROI |
| [SNR 时间诊断](BM_SNR_TIMING_DIAGNOSTICS.md) | 固定 BM2 单个有效匹配 pad，并关联到有可用 noise sigma 的 Step2 波形 | 六对分别使用双方有效的时间；另外保存四探测器共同有效的对照结果 | 按实际匹配 pad 的 SNR 分组，无额外矩形 ROI |
| [角度分组图](TIMING_VS_ANGLE.md) | 匹配 pad 数默认读取已有主分析报告，没有可用设置时取 1 | BM2、Tracker1–3 在四种 CFD 下全部有有效时间；Step4 `useflag=true`，按 event 关联并核对 timestamp | 按 `atan(hypot(kx,ky))` 分组；不额外要求三层 X/Y 全有空间信号，不自动施加 17° 上限或矩形 ROI |
| [空间分格图](BM2_SPATIAL_TIMING.md) | BM2 单个有效匹配 pad，并关联 Step2 的 maxA | BM2、Tracker1–3 在四种 CFD 下全部有有效时间 | 预测位置在 `[-50,50)×[-50,50) mm²` 内分格 |
| [BM cluster 分布](BM_CLUSTER_DISTRIBUTIONS.md) | 按空间簇统计，允许零个或多个匹配簇，不要求 BM2 单 pad 时间 | Step4 `useflag=true` 且三层 Tracker 的 X/Y 六个方向全部有空间信号；不额外要求全部有效时间 | 所选 BM 的预测位置满足矩形 ROI |
| [BM2 位置 cut 时间分析](BM2_POSITION_CUT.md) | 固定 `bmS_n[1]==1` | 四探测器的六个对分别要求双方有效时间；不额外要求三层 X/Y 全有空间信号 | BM2 预测位置满足矩形 ROI |

主时间分析脚本单独运行时，`--bm2-matched-pad-count` 的默认值是 0；managed run 会显式传入合成配置中的值。因此手动复现默认生产样本时，应明确传 `--bm2-matched-pad-count 1` 和参与的探测器列表。

## 计数与时间分辨

“通过单 pad 和位置 cut 的事件数”“某探测器对的有效时间事件数”和“该对迭代 3σ 后的核心事件数”属于不同阶段，不能直接互换。一个事件进入总体样本，不代表它会进入每一个探测器对或每个核心分布。

`weighted_core` 使用各对核心样本标准差的平方，按近似方差估计误差加权求解探测器方差。BM2 与三个 Tracker 对应六个方程、四个未知方差。最后开平方得到各探测器的时间分辨 σ（ns）。高斯图里的拟合 σ 是时间差分布宽度，不能直接当作 BM2 自身的分辨率。

现有误差传播不包含共享事件造成的探测器对相关性及对齐、标定等系统误差。独立的角度、空间、cluster 分布和位置 cut 工具写入 Combine，未接入 managed run 的自动 Results 导出白名单。

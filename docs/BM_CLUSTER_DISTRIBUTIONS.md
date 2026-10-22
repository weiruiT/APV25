# BM 簇数与簇大小分布

## 运行与输入

从 APV25 项目根目录运行：

```bash
cd /home/thoriatong/workplace/APV25
source /home/thoriatong/root/bin/thisroot.sh
root -l -b -q 'BTTracker/macros/PlotBMClusterDistributions.C("batches/run_006","BM2")'
```

第二个参数默认为 `BM2`，也可指定 `BM1`。宏读取该 run 的 `Combine/Tracker-step4-rec.root/Track`、`Combine/Pad-step3-rec.root/TrackArrivalEvents` 和所选 BM 的事件树。

## 共同事件选择

要求 `useflag=true`，并要求三个 Tracker 的 `MMsigx`、`MMsigy` 共六个空间信号标志全部为真。随后按 event 关联 `TrackArrivalEvents`，所选 BM 的预测位置必须满足闭区间 `x∈[-10,40] mm、y∈[-40,10] mm`。范围写在宏顶部常量中，目前没有范围参数；选择 BM1 时也采用这组范围。

不额外要求 BM 或所有 Tracker 的有效定时时间，不施加 `bmS_n[1]==1`。它与主时间分辨及位置 cut 时间分析的样本不同，见 [cut 对照](ANALYSIS_SELECTIONS.md)。

## 两幅图的计数含义

| 图 | 计数方式 |
|---|---|
| `BM2_cluster_number` | 每条入选径迹填一次，使用 `bm_matched_cluster_count[1]`：统计与预测 3×3 pad 窗口相交的空间簇。包含零匹配簇，最后一格合并 ≥9。 |
| `BM2_cluster_size` | 对入选 event 读取 `BM2Events.cluster_size` 中的全部 BM2 簇，不再限制其是否匹配预测窗口。同一事件在同一个 size bin 最多填一次；最后一格合并 ≥20，并按合并后的 bin 去重。 |

例如一个事件的簇大小为 `[1,1,2]`，size 图在 1 和 2 各增加 1；若为 `[20,21]`，只在 ≥20 增加 1。因而 size 图的纵轴是 Events，所有 bin 的计数之和既不等于簇总数，也不一定等于入选事件总数。它与直接 `Draw("cluster_size")` 得到的逐簇分布不同。

number 图标注总 Tracks、`N_cluster=0` 和 `N_cluster≥2`；size 图标注的“大簇事件”比例以全部入选径迹数为分母。零匹配簇只表示预测窗口没有关联到空间簇，不等于整块 BM 没有响应；这样的事件若在其他位置有簇，仍可能贡献 size 图。

旧 ROOT 缺少独立 BM 事件树时会尝试共享 `PadEvents`；缺少匹配簇计数时回退到事件树的 `nclusters`。旧格式的物理范围和匹配含义可能不同，不能直接与新格式结果比较；使用当前重建生成的独立 BM 树可保持明确的探测器定义。

## 输出

输出位于 `run_NNN/Combine/cluster_diagnostics/<execution_id>/`，包含所选探测器的 `BM2_cluster_number.png/.pdf`、`BM2_cluster_size.png/.pdf` 和 `bm_cluster_distributions.root`；选择 BM1 时图片使用 BM1 前缀。图像采用浅蓝填充、蓝色边框和对数纵轴。

这是独立诊断宏，未接入 managed run 的自动阶段和 Results 白名单。

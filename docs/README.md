# 文档入口

日常操作从 [项目总览](../README.md) 和 [BTTracker 使用指南](../BTTracker/README.md) 开始。本文档集保留专题细节：

- [run 与批量数据管理](batch_data_management.md)：run 生命周期、导入、状态和发布规则。
- [分析流程与二次开发](analysis_pipeline_and_customization.md)：步骤、树结构、选择条件与扩展方式。
- [分析样本与 cut 对照](ANALYSIS_SELECTIONS.md)：空间径迹、有效时间、共同样本和逐对样本的区别。
- [BM2 预测位置 cut 时间分辨](BM2_POSITION_CUT.md)：矩形 ROI、五组 run 的运行命令及高斯诊断图。
- [BM 簇数与簇大小分布](BM_CLUSTER_DISTRIBUTIONS.md)：三层 X/Y 响应要求、位置 cut 和按事件计数规则。
- [BM2 SNR 与时间分辨诊断](BM_SNR_TIMING_DIAGNOSTICS.md)：配置、独立运行和中间量导出。
- [时间分辨率与径迹角度](TIMING_VS_ANGLE.md)：ROOT 宏命令、角度分组、事件选择及输出说明。
- [BM2 空间分区时间分析](BM2_SPATIAL_TIMING.md)：可调格子边长的事例数/分辨率地图与 Δt20 对 maxA 二维分布。
- [宏使用说明](macros_usage.md)：临时宏及其输入输出约定。
- [BTTracker 源码导览](BTTracker_src_code_walkthrough.md)：代码结构参考。

run 目录约定和自动生成的状态表见 [batches/README.md](../batches/README.md)。

# 配置目录

`BTTracker/config/` 现在是公共静态配置来源：

- `configCR.json`：阶段开关；
- `map_APV.txt`：APV channel 到三层 Tracker strip；
- `map_APV4BMM.txt`：两块 pad Detector 的 242 个有效位置；
- `channelmap.txt`、`asic2det_64to384.csv`（若存在）：兼容映射；
- `configCR.md`：覆盖规则和开关说明。

公共配置目录不再存放任何生成的 noise 结果或旧四层 2290 通道产物。

每个 run 的新噪声标定写入：

```text
run/Combine/noise/noise-global-parameters-1778.txt/.root
```

当前宇宙线 ROOT 已完成基线扣除和 3σ 通道筛选。signal 不读取 noise pedestal mean，也不做二次基线扣除或通道筛选；Step2 只读取本 run 标定中的 sigma，用作显著峰最低高度阈值。缺少完整 1778 通道标定时 signal 会明确失败。

run 特有配置放在 `run/config/`。`configCR.json` 与公共默认值递归合并，其他同相对路径文件替换公共文件。生成的 pedestal/noise 文件不得放入 `run/config/`。

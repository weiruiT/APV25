# batches

每个 `run_*` 代表一组明确配对的宇宙线与 noise 数据。项目不会根据日期或文件夹名称自动配对现有数据。

```bash
cd ../BTTracker
./init_run.sh ../batches/run_001 --notes "实验条件"
./init_run.sh ../batches/run_001 --signal /path/cosmic.root --noise /path/noise.root
./run_batch.sh ../batches --mode all --force 1 --jobs 2
python3 scripts/run_status.py ../batches
```

每次包含 signal 的分析在所有请求阶段成功后，自动生成 `run_NNN/Results/<execution_id>/`，
只归档最终时间分辨数值、误差、汇总图和来源记录。step1–4、Pad ROOT、残差分布、
噪声标定与逐事件过程保存在 `run_NNN/Combine/`。
`run_NNN/metadata/latest_results.json` 指向最近成功归档。

完整规则见 `../docs/batch_data_management.md`。

`run_status.tsv` 汇总每个受管理 run 的最近一次执行、阶段状态、配置中的 SNR 区间、结果导出编号和输入 ROOT 数量。初始化 run 和批处理结束时会自动刷新；也可在 `BTTracker/` 手动运行 `python3 scripts/run_status.py ../batches`。这是本地生成的导航索引，不纳入代码备份。

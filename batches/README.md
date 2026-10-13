# batches

每个 `run_*` 代表一组明确配对的宇宙线与 noise 数据。项目不会根据日期或文件夹名称自动配对现有数据。

```bash
cd ../BTTracker
./init_run.sh ../batches/run_001 --notes "实验条件"
./init_run.sh ../batches/run_001 --signal /path/cosmic.root --noise /path/noise.root
./run_batch.sh ../batches --mode all --force 1 --jobs 2
```

完整规则见 `../docs/batch_data_management.md`。

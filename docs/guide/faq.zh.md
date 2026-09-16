# 常见问题

## 常见问题

- 样本表解析错误：确保第一列为 `sample_id`，使用绝对路径
- 缺少 `${sample}.mcds`：检查 `ALLCOOLS_BAM_TO_ALLC` 是否生成了单细胞 `*_allc.gz`，以及 `chrom_size_path` 是否正确
- 卡在 Bismark 步骤：验证参考基因组索引以及 `params.bismark_ref` 在容器中是否可见
- 恢复运行：使用 `-resume` 并指定相同的 `-w` 工作目录

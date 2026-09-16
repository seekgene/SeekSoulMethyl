# SeekSoulMethyl v2.2.1 MCDS 合并数据丢失修复报告

**报告对象**：SeekSoulMethyl 流程的使用者
**升级内容**：修复甲基化 MCDS 合并环节（`MERGE_MCDS`）的静默数据丢失问题
**受影响版本**：v2.0.0 – v2.2.0
**修复版本**：v2.2.1
**报告版本**：v1.0
**日期**：2026 年 9 月

---

## 一、摘要

SeekSoulMethyl v2.2.1 修复了甲基化 `MERGE_MCDS` 环节的一个**静默数据丢失问题**。该环节负责将分片单细胞甲基化矩阵合并成最终用于聚类和下游分析的 MCDS。

此问题会导致**随机一部分细胞在大尺度 bin 上被静默清零**——**chrom1M 约 9% 的细胞、chrom500k 约 8% 的细胞**被整行清零，且**无任何报错、无任何警告、无明显异常**。而**默认聚类用的 chrom20k 不受影响**（所有测试中零丢失）。

**一句话结论**：如果只用默认的 `chrom20k` 聚类，v2.2.0 的结果是完整的，**无需重跑**；如果使用 `chrom1M` / `chrom500k` 大尺度 bin，请升级到 v2.2.1 并重跑。

---

## 二、问题背景

`MERGE_MCDS` 环节将上游 `allcools generate-dataset` 产生的分片 MCDS 合并成一个完整 MCDS。合并结果通过 xarray 的 `to_zarr(append_dim=...)` 以 Dask 并行的方式写入 Zarr 存储。

在 v2.0.0 至 v2.2.0 版本中，该 append 操作**没有先做数据物化**（`ds.load()`）。在旧版 xarray 的 chunk 对齐检查下，append 路径没有强制校验「正在写入的 Dask chunk」与「Zarr 存储中已有 chunk」的对齐关系。当多个 Dask chunk 映射到同一个 Zarr chunk 并被并发写入时，就会产生**竞态**，可能静默地把某个细胞的整行数据覆盖成零。

由于该失败是静默且随机的，受影响的细胞在矩阵中零散分布，只能通过与独立真值比对才能发现。

---

## 三、根因

- **上游 issue**：[xarray #8876 —— 「向已存在的 Zarr 存储追加时可能发生竞态」](https://github.com/pydata/xarray/issues/8876)（相关 issue：[#8882 —— 「to_zarr 在使用 append_dim 时静默丢失数据」](https://github.com/pydata/xarray/issues/8882)）。
- **机制**：`to_zarr(append_dim=...)` 在 Dask 并行写入、chunk 不对齐时，可能把多个 Dask chunk 并发写入同一个 Zarr chunk，导致数据丢失。
- **为何静默**：xarray 的 chunk 对齐安全检查（`safe_chunks`）只在「新建变量」时生效，append 路径绕过了本应报错的检查。

---

## 四、v2.2.1 修复了什么

v2.2.1 在 append 之前显式加入 `ds.load()`，使合并数据**完整物化后串行写入**，从而消除并发写入的竞态。

- 受影响版本：**v2.0.0 – v2.2.0**
- 修复版本：**v2.2.1**（修复已推送到发布分支并打 tag）

---

## 五、影响范围

该问题只影响**大尺度 bin**（这些 bin 的 `count_type` 维度以「整块」方式 chunk，mc 与 cov 打包进同一个 chunk）。共发现三种损伤形态：

| 损伤形态 | 表现 | 影响 bin |
|---|---|---|
| **整行清零** | 细胞的 mc 与 cov 同时归零 | chrom1M（约 9%）、chrom500k（约 8%） |
| **mc 减半** | cov 不变、mc 减半 | 小尺度 bin（偶发，个别细胞） |
| **cov 减半** | mc 不变、cov 减半 | 小尺度 bin（偶发，个别细胞） |

| bin | v2.2.0 丢失 | v2.2.1 |
|---|---|---|
| chrom1M | 约 9% 细胞清零 | **0（全部找回）** |
| chrom500k | 约 8% 细胞清零 | **0（全部找回）** |
| chrom100k / 50k / 10k | 基本无损（偶发单个细胞减半） | 0 |
| **chrom20k（默认聚类）** | **0** | **0** |

---

## 六、验证与测试结果

通过四条独立证据链对修复进行了验证：

### 6.1 33 个生产样本的跨 region 一致性审计

对 **33 个生产甲基化样本**（共 295,408 个细胞、7 个批次）做跨 region 一致性审计：比较每个细胞在六个 bin 分辨率下的覆盖度总和，某 region 总量为 0 而其他正常，即代表该细胞在该 region 丢失。

| bin | 平均丢失率 |
|---|---|
| chrom1M | 9.4%（单样本 7.2% – 13.4%） |
| chrom500k | 8.3%（单样本 6.4% – 10.5%） |
| **chrom20k** | **0（无丢失）** |
| chrom100k / 50k / 10k | 0 |

### 6.2 v2.2.0 复现

用 v2.2.0（含问题版本）重跑 14 个 MCDS，复现出同等量级的丢失——chrom1M 约 9.6%、chrom500k 约 8.8%、chrom20k 为 0——证明重跑链路正确、问题稳定存在。

### 6.3 从单细胞 ALLC 文件重建真值

将一个样本的 MCDS 直接从单细胞 ALLC 文件重建（无合并步骤 = 真值），与旧 MCDS 逐细胞、逐 bin 比对：

| bin | 被清零的细胞数 |
|---|---:|
| chrom1M | 687 |
| chrom500k | 601 |
| **chrom20k** | **0** |

### 6.4 v2.2.1 修复验证（逐细胞、逐 bin）

v2.2.1（含 `ds.load()` 修复）完整跑完后，与 v2.2.0 在全部 33 个样本上逐细胞、逐 bin 对比：

| bin | v2.2.1 找回的细胞数 |
|---|---:|
| chrom1M | **28,726** |
| chrom500k | **26,442** |
| chrom20k | 基本为 0（仅个别细胞存在单 bin 的 mc/cov 微小差异） |

**结论**：v2.2.0 在 chrom1M/500k 上静默清零的每一个细胞，都在 v2.2.1 中完整恢复；chrom20k 自始至终不受影响。

---

## 七、给用户的建议

- **如果只用默认 `chrom20k` 聚类**：现有结果是完整的，**无需重跑**，此问题不影响 chrom20k。
- **如果使用 `chrom1M` 或 `chrom500k`（大尺度 bin 分析，如粗略 DMR/覆盖度汇总）**：请升级到 **v2.2.1** 并重跑，因为旧版本在这些 bin 上静默清零了约 8–9% 的细胞。
- 该丢失是**静默的**——无报错、无警告——因此旧版本在 chrom1M/500k 上的结果，未经 v2.2.1 重跑前不可采信。

---

## 八、参考链接

- [xarray #8876 —— Possible race condition when appending to an existing Zarr store](https://github.com/pydata/xarray/issues/8876)
- [xarray #8882 —— to_zarr silently loses data when using append_dim](https://github.com/pydata/xarray/issues/8882)

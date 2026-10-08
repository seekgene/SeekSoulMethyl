# 算法与流程原理

## 流程详情

### 转录组处理工作流
转录组数据使用 [SeekSoulTools](https://github.com/seekgene/SeekSoulMethyl/tree/nf_rna_methy/dependence/seeksoultools) 进行分析。详细步骤请参阅官方[算法概述](https://seeksoul.online/cloudplatform-doc/en/document/Software/3_SeekSoul_Tools/1_v1.3.0a/3_Use_guide/1_rna/4_Algorithm.html)。下游甲基化文库使用的细胞由转录组文库的细胞 barcode 确定。

### 甲基化处理工作流
#### 步骤 1：预处理与 Barcode 解析
**Barcode 提取与纠错**

根据设计的文库结构，我们在 read 中定位 barcode 并提取相应序列。如果提取的 barcode 存在于白名单中，则计为有效 barcode；否则，SeekSoulMethyl 将尝试进行 barcode 纠错，前提是该 barcode 与白名单中的某个条目仅有一个碱基的差异（汉明距离 = 1）：
* 如果恰好有一个白名单候选匹配：将无效 barcode 纠正为该白名单 barcode。
* 如果有多个白名单候选匹配：纠正为 read 计数最高的候选。

如果纠错失败，该 read 将被丢弃，视为最终无效 barcode read。

**UMI 提取**

根据设计结构中的 UMI 位置直接提取 UMI 序列，不进行纠错。

**Forward 和 Reverse reads 判定**

在 17L 和 ME 对应位置中，有 7 个碱基可能为 C 或转换后的 T（如下高亮所示）。我们使用第一个和最后两个 C/T 位置来判定 forward 与 reverse reads：如果三个位置均为 C，则判定为 reverse read；否则判定为 forward read。

- Forward：**T**gt**TT**gt**T**gttg**T**t**T**gtAGATGTGTATAAGAGA**T**

- Reverse：**C**gt**CC**gt**C**gttg**C**t**C**gtAGATGTGTATAAGAGA**C**

Reverse reads 对应甲基化术语中的 CTOT/CTOB（原始链的反向互补链）；forward reads 对应 OT/OB（原始链）。
"forward"或"reverse"的判定结果标注在 read name 中。

<figure style="text-align: center;">
<img src="../images/fr_determinate.png" alt="Forward 和 Reverse reads 判定" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">图 3. Forward 和 Reverse reads 判定</figcaption>
</figure>

**C–T 转换率**

我们利用 17L 和 ME 序列中原始 C 位置来计算 C-to-T 转换率。由于这些是容易出现测序错误的固定序列，我们仅对结构验证通过的 reads 进行计算：

 - 在 DD-MET3 中，7F 序列必须为 `TTGCTGT` 或 `TTGTTGT`，17L 和 ME 跨区域的序列必须为 `GTAGATGTGTATAAGAGA`，且第一个和最后两个原始 C 位置的碱基必须为 T。
 
 - 在 DD-MET5 中，17L 和 ME 跨区域的序列必须为 `GTAGATGTGTATAAGAGA`，且第一个和最后两个原始 C 位置的碱基必须为 T。

对于保留的 reads，我们提取相应位置的碱基来计算 C-to-T 转换率：

<figure style="text-align: center;">
<img src="../images/CT_conversion.png" alt="CT 转换率" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">图 4. CT 转换率</figcaption>
</figure>

> [!NOTE]
> 以上过滤步骤仅用于计算 C-to-T 转换率；不满足这些条件的 reads 不会从最终输出的 FASTQ 文件中过滤掉。

**人工序列去除**

根据文库设计中的预定义位置，从 Read1 中去除 TSO/7F 连接序列、17L 和 ME 序列。

**接头修剪**

使用 cutadapt 去除因 R1/R2 read-through 事件（双端测序 reads 重叠）引入的 ME 接头序列。

**修剪 Tn5 转座酶引入的 9 bp 间隔**

在去除接头和其他人工序列后，我们额外修剪 Tn5 转座引入的插入片段两侧的 9 bp 间隔区域。这些 9 bp 区域可能携带人工甲基化信号，并在插入边界附近虚假地提高 CH 甲基化水平，因此在下游分析前予以去除。

**过滤含有过多非 CpG 甲基化 C 碱基的 reads（可选）**

根据 read pair 中非 CpG 甲基化 C 碱基的数量进行过滤。默认情况下，read1/read2 中检测到 > 2 个非 CpG 甲基化 C 的 read pair 将被去除。如果不需要启用此过滤，请将 filter_ch 设置为 0。

> [!NOTE]
> 此过滤策略基于 Lu 等人 [1] 的研究发现，该研究表明合成接头中的切口（nick）可触发 Bst 聚合酶的切口平移活性。该活性会掺入 5-甲基-dCTP，导致完全未转换的 reads，表现为人工甲基化信号。

**过滤过短 reads**
经过前述过滤和接头修剪步骤后，如果 read pair 中 R1 的长度小于 20 bp 或 R2 的长度小于 60 bp，则该 read pair 将被过滤掉。

#### 步骤 2：Bismark 比对与按名称排序
**比对与标签添加**

我们使用 Bismark 进行甲基化比对。我们修改的 [Bismark](https://github.com/seekgene/Bismark.git) 添加了 `--add_barcode` 和 `--add_umi` 参数，通过 read name 在 BAM 文件中标记 CB（纠错后 barcode）和 UR（原始 UMI）。对于 forward reads，我们使用 `-X 1000` 允许最大插入片段长度为 1000 bp；对于 reverse reads，我们使用 `--pbat` 和 `-X 1000`。
比对完成后，使用 `samtools sort -n` 按 read name 排序 BAM 文件；按名称排序的 BAM 文件作为下游分析的输入。

#### 步骤 3：ALLCools 分析

**按细胞 barcode 拆分**

将按名称排序的 BAM 文件按 RNA 来源的细胞 barcode 拆分为单细胞 BAM 文件，每个文件包含一个细胞的唯一比对 reads。

**生成 ALLC 文件**

将每个单细胞 BAM 文件按位置排序，并使用 ALLCools `bam-to-allc` 转换为 ALLC 格式。我们修改的 ALLCools 基于 UR 标签对每个 C 位点进行 UMI 纠错和去重。

<figure style="text-align: center;">
<img src="../images/umi_correction_detailed_diagram_en.png" alt="UMI 纠错详细流程图" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">图 5. UMI 纠错详细流程图</figcaption>
</figure>

详细信息请参阅 [SeekGene ALLCools 仓库](https://github.com/seekgene/ALLCools)。

**生成 MCDS**

运行 `allcools generate-dataset` 对基因组进行分箱（chrom10k/20k/50k/100k/500k/1M/geneslop2k），并计算单细胞甲基化矩阵。geneslop2k 分箱定义为每个基因两侧各延伸 2k bp 的区域。

#### 步骤 4：降维与聚类
默认使用 ALLCools 对 chrom20k 分箱进行 LSI 降维，随后进行 UMAP 可视化。

### 系统要求
如果使用 `sc_methy_workflow.sh`，系统要求如下：
- **CPU**：64 核（推荐）
- **内存**：128GB RAM（推荐）
- **操作系统**：Linux（推荐 Ubuntu 18.04+ 或 CentOS 7+）

## Nextflow 分步详情


本节描述每个 Nextflow 进程的输入、核心逻辑、关键参数和输出，便于故障排除和结果解读。工作流定义在 `nf/main.nf` 中，进程实现在 `nf/modules/*.nf` 中。

### 步骤 1：预处理与 Barcode 解析（step1.nf）
- 计算全基因组 CpG 位点数：`COMPUTE_CPG_SITES`
  - 输入：`params.genomefa`、`params.chrom_size_path`
  - 操作：运行 `count_cg_sites.py` 统计 CpG 位点
  - 输出：`genome_cg_info.json`

- 转录组 FASTQ 质控（多组）：`FASTP_EXPRESSION_MULTI`
  - 输入：每个样本的双端 FASTQ（G1/G2/... 组）
  - 操作：`fastp` 修剪和质控
  - 输出：清洗后的 `*_expression_clean_R1/2.fastq.gz`、`*.html`、`*.json`

- 甲基化 FASTQ 质控（多组）：`FASTP_METHYLATION_MULTI`
  - 输入：每个样本的双端 FASTQ（G1/G2/... 组）
  - 操作：`fastp` 质控（禁用接头检测，按流程设置进行修剪）
  - 输出：清洗后的 `*_methylation_clean_R1/2.fastq.gz`、`*.html`、`*.json`

- RNA 比对与定量：`SEEKSOULTOOLS_RNA`
  - 输入：清洗后的转录组 R1/R2 列表
  - 操作：`seeksoultools rna run`（STAR 比对、计数、过滤、聚类、差异表达分析）
  - 输出：`Analysis/step3/filtered_feature_bc_matrix/` 等

- 甲基化 barcode 解析：`METHYLATION_BARCODE_EXTRACTION`
  - 输入：清洗后的甲基化 R1/R2 列表，白名单 `params.methy_barcode_wl`
  - 操作：运行 `barcode_cs_multi.py`，根据 `params.chemistry` 解析/纠错细胞 barcode 和 UMI；可选 `--split_fastq n` 根据 barcode 前 n 个碱基分片
  - 输出：`step1/${sample}_forward_*{1,2}.fq.gz`、`step1/${sample}_reverse_*{1,2}.fq.gz`、`${sample}_methy_summary.json`、可选 `${sample}_barcode_stats.txt`

- 构建 forward/reverse 配对列表：`PARSE_FASTQ_FILES`
  - 输入：forward/reverse FASTQ 文件集
  - 操作：扫描 `step1/` 并按标识符配对文件
  - 输出：`forward_pairs.txt`、`reverse_pairs.txt`

- Barcode 提取后质控：`FASTP_METHYLATION_BARCODE_EXTRACT`
  - 输入：配对的子 FASTQ
  - 操作：`fastp` 质控
  - 输出：每对 `*.html`、`*.json`

### 步骤 2：Bismark 比对与 BAM 排序（step2.nf）
- Forward 链比对：`BISMARK_ALIGNMENT_FORWARD`
  - 关键参数：`--add_barcode`、`--add_umi`；最大插入片段长度 `-X 1000`
  - 输出：`*_bismark_bt2_pe.bam`、`*_bismark_bt2_PE_report.txt`

- Reverse（PBAT）比对：`BISMARK_ALIGNMENT_REVERSE`
  - 关键参数：`--pbat`、`--add_barcode`、`--add_umi`
  - 输出：同上

- 按 read name 排序：`SORT_BAM_BY_NAME`
  - 操作：`samtools sort -n`
  - 输出：`*_sortbyname.bam`

### 步骤 3：单细胞拆分、ALLC 生成/合并、数据集构建（step3.nf）
- 按细胞 barcode 拆分 BAM：`SPLIT_BAM_FILES`
  - 输入：按名称排序的 BAM 和 GEX barcodes `barcodes.tsv.gz`
  - 操作：运行 `step3_split_bams.py`，按细胞 barcode 拆分并仅保留共有细胞
  - 输出：单细胞 BAM 目录、`*_filtered_barcode`、`*_filtered_barcode_reads_counts.csv`

- 合并 forward/reverse 单细胞 BAM 和计数：`MERGE_BISMARK_BAM`
  - 操作：对匹配的细胞使用 `samtools merge -n`；合并/去重 barcodes 和 read 计数
  - 输出：`*_merged_fr_bam/`、`*_merge_filtered_barcode`、`*_merge_filtered_barcode_reads_counts.csv`

- 生成单细胞 ALLC：`ALLCOOLS_BAM_TO_ALLC`
  - 操作：运行 `step3_bam_to_allc.py`（自定义 ALLCools），基于 UR 标签进行去重和甲基化定量
  - 输出：单细胞 `*_allc.gz` 和索引

- 合并细胞指标：`MERGE_FILTERED_BARCODE_READS_COUNTS`
  - 操作：合并 barcodes 和 read 计数，创建 `${sample}_cells.csv` 和 `.json`
  - 输出：`filtered_barcode`、`filtered_barcode_reads_counts.csv`、`${sample}_cells.{csv,json}`

- 构建多尺度数据集：`ALLCOOLS_GENERATE_DATASETS`
  - 操作：`allcools generate-dataset`，对 chrom10k/20k/50k/100k/500k/1M 进行分箱，计算 `count` 和 `hypo-score` 等指标
  - 输出：`${sample}.mcds`

- 合并 ALLC（分片时）：`ALLCOOLS_SUBMERGE`、`ALLCOOLS_MERGE`
  - 操作：合并每个分片/每个样本的 ALLC
  - 输出：`${sample}_merge_allc.gz` 和索引

- 提取 CG 上下文 ALLC：`ALLCOOLS_EXTRACT`
  - 操作：`allcools extract-allc --mc_contexts CGN --strandness merge`
  - 输出：`*.CGN-Merge*`

（仅甲基化工作流 `methy_only.nf` 额外包含 `COUNTS_MAPPED_READS` 和 `ESTIMATED_CELLS`，用于基于 read 计数的细胞估算和 barcode 过滤）

### 步骤 4：汇总、降维与联合报告（step4.nf）
- 甲基化汇总：`METHYLATION_SUMMARY`
  - 操作：`step4_wgs_summary.py` 聚合 Bismark 报告、细胞指标和 CpG 统计数据，生成 `${sample}_methy_summary.json` 和 `${sample}_wgs_summary.csv`

- LSI/PCA 聚类与可视化：`METHYLATION_LSI_PCA_CLUSTERING`
  - 操作：`step4_allcools_PCA_cluster.py --var_dim chrom20k --reduc lsi`
  - 输出：`*.h5ad`、`*.pdf`、`*.png`

- 转录组 + 甲基化联合报告：`MULTI_REPORT`
  - 操作：`step4_report_rna_met.py` 整合转录组和甲基化输出
  - 输出：`${sample}_rna_methyl_report.html`、`${sample}_rna_met.json`

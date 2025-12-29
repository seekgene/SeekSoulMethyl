# 如何对单细胞 BAM 文件进行去重

## 简介

在进行拷贝数变异（CNV）分析或其他对 Read 计数敏感的下游分析时，必须确保输入数据的准确性。单细胞甲基化文库构建过程中引入了随机打断和 PCR 扩增，为了消除 PCR 偏好性对分析结果的影响，我们需要对 BAM 文件进行去重。

本此去重流程包含两个关键步骤：
1. **Flag 矫正 (Flag Correction)**：修正 Bismark 比对产生的 Flag 异常。
2. **UMI 去重 (UMI Deduplication)**：基于 UMI 和比对位置进行去重。

## 准备工作

### 输入数据
已经拆分好的单细胞 BAM 文件。

### 依赖工具
请确保您的环境中安装了以下工具：
- **Python 3**
- **[umi_tools](https://umi-tools.readthedocs.io/en/latest/QUICK_START.html)**
- **samtools** (依赖于 pysam)

安装 `umi_tools`:
```shell
pip install umi_tools
```

### 必需脚本
您需要下载以下两个脚本：
- [`step1_swap_flags_on_reverse.py`](../nf/bin/utils/step1_swap_flags_on_reverse.py)
- [`step2_umi_tools_dedup.py`](../nf/bin/utils/step2_umi_tools_dedup.py)

## 详细步骤

### 步骤 1：调整 Reads 的 Flag 值 (Flag Correction)

**背景**：
Bismark 比对工具在处理 OT/OB（原始链）和 CTOT/CTOB（互补链）时，会[自动调换 CTOT/CTOB reads 的 Flag 值](https://github.com/FelixKrueger/Bismark/issues/41)。具体表现为 R1 被标记为 "second in pair"，而 R2 被标记为 "first in pair"。为了确保后续去重工具（如 umi_tools）能正确识别 R1 和 R2，我们需要将这些 Flag 恢复正常。

**操作命令**：

```shell
# 配置参数
outdir="./"
input_dir="/path/to/step3/split_bams/merged/"  # 单细胞 BAM 所在目录
workers=16                                      # 并行数

# 1. 生成文件列表
find -L "$input_dir" -name "*.bam" | sort > "$outdir/bam_list.txt"

# 2. 创建输出目录
mkdir -p $outdir/flag_reverse/

# 3. 执行 Flag 矫正
python step1_swap_flags_on_reverse.py \
  --bam_list $outdir/bam_list.txt \
  --out_dir $outdir/flag_reverse \
  --max_workers $workers
```

**输出**：
调整 Flag 后的 BAM 文件将生成在 `$outdir/flag_reverse/` 目录下，文件名后缀为 `.swapflags.bam`。

### 步骤 2：UMI 去重 (UMI Deduplication)

**背景**：
由于同一个文库片段末端是随机打断的，我们采用 **R1 比对位置 + UMI 序列** 作为去重的唯一标识（Unique Identifier）。

**操作命令**：
使用 `umi_tools` 进行去重，关键参数为 `--extract-umi-method tag --umi-tag UR --paired --ignore-tlen`。

```shell
# 1. 生成矫正后的 BAM 列表
find "$outdir/flag_reverse" -name "*.bam" | sort > "$outdir/swapflags_bam_list.txt"

# 2. 创建输出目录
mkdir -p $outdir/umi_tools_swapflags_dedup_bam/

# 3. 执行去重
python step2_umi_tools_dedup.py \
  --bam_list $outdir/swapflags_bam_list.txt \
  --out_dir $outdir/umi_tools_swapflags_dedup_bam/ \
  --max_workers $workers \
  --skip_existing
```

**输出**：
去重后的 BAM 文件将存储在 `$outdir/umi_tools_swapflags_dedup_bam/` 目录下，文件名后缀为 `.swapflags.dedup.bam`。这些文件即为最终去重后的单细胞 BAM 文件。

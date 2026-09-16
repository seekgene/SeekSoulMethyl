# 使用方法

## 使用方法


### 激活环境
```bash
conda activate seeksoulmethyl
```

### 运行双组学分析（Shell 脚本）
```bash
# sc_methy_workflow.sh 位于您克隆的 SeekSoulMethyl 目录下
bash sc_methy_workflow.sh \
/path/to/expression_R1.fastq.gz \
/path/to/expression_R2.fastq.gz \
/path/to/methy_R1.fastq.gz \
/path/to/methy_R2.fastq.gz \
--sample WTJW880 \
--outdir /path/to/results \
--database_dir /path/to/human-reference-GRCh38 \
--chemistry DD-MET5 \
--core 64 \
--filter_ch 2
```
如果一个样本有多组数据，请用逗号分隔文件路径。请确保各组数据的 FASTQ 文件按正确顺序排列。
```shell
bash sc_methy_workflow.sh \
/path/to/WTJW969_E_L003_R1.fq.gz,/path/to/WTJW969_E_L004_R1.fq.gz \
/path/to/WTJW969_E_L003_R2.fq.gz,/path/to/WTJW969_E_L004_R2.fq.gz \
/path/to/WTJW969_Met_L000_R1.fq.gz,/path/to/WTJW969_Met_L001_R1.fq.gz,/path/to/WTJW969_Met_L002_R1.fq.gz,/path/to/WTJW969_Met_L003_R1.fq.gz,/path/to/WTJW969_Met_L004_R1.fq.gz \
/path/to/WTJW969_Met_L000_R2.fq.gz,/path/to/WTJW969_Met_L001_R2.fq.gz,/path/to/WTJW969_Met_L002_R2.fq.gz,/path/to/WTJW969_Met_L003_R2.fq.gz,/path/to/WTJW969_Met_L004_R2.fq.gz \
--sample WTJW969 \
--outdir /path/to/results \
--database_dir /path/to/human-reference-GRCh38 \
--chemistry DD-MET5 \
--core 64 \
--filter_ch 2
```

### 输入参数：

- **$1**：单细胞转录组 Read1 fastq 文件路径。
- **$2**：单细胞转录组 Read2 fastq 文件路径。
- **$3**：单细胞甲基化 Read1 fastq 文件路径。
- **$4**：单细胞甲基化 Read2 fastq 文件路径。
- **sample**：样本名称。
- **outdir**：输出目录路径。
- **database_dir**：参考基因组数据库路径。
- **chemistry**：甲基化建库方案（DD-MET3 或 DD-MET5；默认值：DD-MET5）。DD-MET3 为双标签数据集，表示同一细胞的 RNA 和 DNA 甲基化数据 barcode 不同；DD-MET5 为单标签，表示同一细胞的 RNA 和 DNA 甲基化数据 barcode 相同。
- **core**：CPU 核心数。
- **filter_ch**：过滤包含 > n 个 CH 甲基化位点的 reads。如果不需要启用此过滤，请将 filter_ch 设置为 0。

## 使用 Nextflow 运行（推荐）


如果需要批量处理样本并获取流程级报告，请使用 Nextflow：

<figure style="text-align: center;">
<img src="./images/nf_SeekSoulMethyl_workflow.png" alt="SeekSoulMethyl 流程" width="600" style="max-width: 100%; height: auto;" />
<figcaption style="font-size: 0.95em; color: #666; margin-top: 4px;">图 6. SeekSoulMethyl Nextflow 流程工作流</figcaption>
</figure>

1. 安装 Nextflow：
```bash
conda install -n seeksoulmethyl -c bioconda nextflow
```
2. 准备输入样本表
```
cat > samplelist.csv << EOF
sample_id,expression_r1,expression_r2,methylation_r1,methylation_r2
XYRD-WTJW880,/path/to/XYRD-WTJW880-E_S1_L005_R1_001.fastq.gz,/path/to/XYRD-WTJW880-E_S1_L005_R2_001.fastq.gz,/path/to/XYRD-WTJW880-MET_S01_L001_R1_001.fastq.gz,/path/to/XYRD-WTJW880-MET_S01_L001_R2_001.fastq.gz
EOF

# expression_r1：单细胞转录组 Read1 fastq 文件
# expression_r2：单细胞转录组 Read2 fastq 文件
# methylation_r1：单细胞甲基化 Read1 fastq 文件
# methylation_r2：单细胞甲基化 Read2 fastq 文件
```

如果单个样本有多组数据，且转录组和甲基化的 FASTQ 数量不匹配（例如 WTJW969），请在样本表中添加多行，每行代表一组数据。
```text
sample_id,expression_r1,expression_r2,methylation_r1,methylation_r2
WTJW969,/path/to/WTJW969_E_L003_R1.fq.gz,/path/to/WTJW969_E_L003_R2.fq.gz,/path/to/WTJW969_Met_L000_R1.fq.gz,/path/to/WTJW969_Met_L000_R2.fq.gz
WTJW969,/path/to/WTJW969_E_L004_R1.fq.gz,/path/to/WTJW969_E_L004_R2.fq.gz,/path/to/WTJW969_Met_L001_R1.fq.gz,/path/to/WTJW969_Met_L001_R2.fq.gz
WTJW969,,,/path/to/WTJW969_Met_L002_R1.fq.gz,/path/to/WTJW969_Met_L002_R2.fq.gz
WTJW969,,,/path/to/WTJW969_Met_L003_R1.fq.gz,/path/to/WTJW969_Met_L003_R2.fq.gz
WTJW969,,,/path/to/WTJW969_Met_L004_R1.fq.gz,/path/to/WTJW969_Met_L004_R2.fq.gz
```

3. 运行流程：
```bash
nextflow run -bg SeekSoulMethyl/nf/main.nf \
--outdir /path/to/tiny_demo/results/ \
--samplesheet samplelist.csv \
-w /path/to/tiny_demo/results/work \
-c SeekSoulMethyl/nf/nextflow.config \
-profile aliyun_k8s \
--database_dir /path/to/human-reference-GRCh38/ \
--split_fastq 1 \
--filter_ch 2 \
--chemistry DD-MET5 > methy.log

# --outdir：最终结果目录
# --samplesheet：输入样本表文件
# -w：Nextflow 工作目录
# -c：Nextflow 配置文件。必须根据您的服务器配置进行修改。
# -profile：阿里云 K8s 的 Nextflow profile
# --database_dir：参考基因组数据库目录
# --split_fastq：为加速分析过程，根据 barcode 前 n 个碱基拆分 fastq。默认值为 4。
# --filter_ch：过滤包含 > 2 个 CH 甲基化位点的 reads。
# --expected_cell_num：methy_only 细胞估计使用的预期细胞数。默认值为 3000。
# --chemistry：甲基化建库方案（DD-MET3 或 DD-MET5；默认值：DD-MET5）
```

### nextflow.config 说明

`nextflow.config` 声明了执行环境、资源配额和运行策略。您必须根据自己的服务器或集群进行自定义配置。

- 位置：`SeekSoulMethyl/nf/nextflow.config`（运行时通过 `-c SeekSoulMethyl/nf/nextflow.config` 指定）。
- 需要根据您的基础设施调整的关键配置项：
  - 执行器：`process.executor`（例如 `local`、`slurm`、`pbs`、`k8s`、`awsbatch`）。
  - 资源：`process.cpus`、`process.memory`、`process.time`，或通过 `withLabel`/`withName` 进行细粒度配置。
  - 工作目录：`workDir`（也可通过 `-w` 设置）；确保该目录可写且有足够空间。
  - 环境：根据您的服务器支持情况，启用 `conda.enabled`、`docker.enabled` 或 `singularity.enabled` 之一。

示例配置（请将路径和参数替换为您服务器上的有效值）：

```groovy
profiles {
  // 本地机器
  local {
    process.executor = 'local'
    workDir          = '/path/to/work'
    process.cpus     = 8
    process.memory   = '32 GB'
    conda.enabled    = true
    // 或使用容器：singularity.enabled = true / docker.enabled = true
  }

  // Slurm 集群（HPC）
  slurm {
    process.executor   = 'slurm'
    workDir            = '/lustre/work/your_user'
    process.cpus       = 8
    process.memory     = '32 GB'
    process.queue      = 'normal'
    process.clusterOptions = '-A your_account --qos=normal'

    withLabel: 'high_mem' {
      cpus   = 16
      memory = '64 GB'
    }
  }

  // Kubernetes（例如阿里云 ACK）
  aliyun_k8s {
    process.executor   = 'k8s'
    workDir            = '/mnt/nf-work'    // 持久化存储卷路径
    k8s {
      namespace        = 'your-namespace'
      storageClaimName = 'your-pvc'
      cpu              = 4
      memory           = '16 GB'
    }
    // 如果全局使用容器镜像：docker.enabled = true
  }
}
```

提示：
- 选择与您环境匹配的 `-profile`（例如 `slurm`、`local`、`aliyun_k8s`），然后调整相应参数。
- 将 `workDir` 设置在具有足够容量的高速存储上；最终结果目录由 `--outdir` 控制。
- 如果使用 README 中的 conda 环境，建议启用 `conda.enabled`；如果您的集群使用容器，请使用 Docker/Singularity。

参考文档：
- Nextflow 配置与 profiles：https://www.nextflow.io/docs/latest/config.html
- 执行器（local、Slurm、K8s 等）：https://www.nextflow.io/docs/latest/executor.html
- Kubernetes 指南：https://www.nextflow.io/docs/latest/kubernetes.html

## 仅甲基化工作流（测试版本，目前不推荐使用）

当您仅有甲基化数据时，可使用简化工作流：
```bash
nextflow run SeekSoulMethyl/nf/main.nf \
  --workflow methy_only \
  --outdir /path/to/results \
  --samplesheet samplelist.csv \
  -w /path/to/work \
  -c SeekSoulMethyl/nf/nextflow.config \
  -profile aliyun_k8s \
  --database_dir /path/to/reference \
  --split_fastq 4 \
  --filter_ch 2 \
  --chemistry DD-MET5
```

## 关键参数与提示

- `--database_dir`：参考数据库目录，包含 `fasta/genome.fa`、`genes/genes.gtf`、`star/`、`bed/chr_nochrM.bed`
- `--chemistry`：`DD-MET3` 或 `DD-MET5`（默认值：`DD-MET5`）；同时设置转录组建库方案和 barcode 白名单
- `--split_fastq`：根据 barcode 前 n 个碱基进行分片（默认值 4；增加并行度但会增加调度/合并开销）
- `--filter_ch`：过滤包含 > n 个 CH 甲基化位点的 reads（默认值 2）。如果不需要启用此过滤，请将 `filter_ch` 设置为 0。
- 样本表表头必须包含 `sample_id`

## 执行环境与资源

- `-profile docker` 使用公开容器镜像 `ghcr.io/seekgene/seeksoulmethyl_docker:v1.0.0`。
- `-profile aliyun_k8s` 和 `-profile aliyun_k8s_argo` 使用内部阿里云镜像 `seekgene-registry-vpc.cn-beijing.cr.aliyuncs.com/seekgene/seeksoulmethyl:1.1.2`，用于 SeekGene Kubernetes 环境。
- 请根据实际基础设施选择 profile；如需镜像或替换容器，请修改 `nf/nextflow.config`。
- 资源密集型步骤：Bismark/ALLCools 需要大量 CPU/内存；默认值在 `withName` 块中设置，如有需要请相应调大

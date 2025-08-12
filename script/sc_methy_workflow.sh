#!/bin/bash
# 寻因单细胞转录组+甲基化双组学分析流程
# 适用于分析北京寻因生物有限公司单细胞转录组+甲基化双组学试剂盒产品
# 输入参数说明：
# $1: 样本名称
# $2: 单细胞转录组Read1的fastq文件
# $3: 单细胞转录组Read2的fastq文件  
# $4: 单细胞甲基化Read1的fastq文件
# $5: 单细胞甲基化Read2的fastq文件
# $6: 结果输出路径
# $7: 参考基因组数据库路径
# $8: CPU数量

sample=${1}
exp_fq1=${2}
exp_fq2=${3}
methy_fq1=${4}
methy_fq2=${5}
outdir=${6}
database_dir=${7}
core=${8}

# 设置bismark使用的CPU数量，建议为总CPU的1/8
# 注意：如果bismark --parallel参数设置为8，确保机器的CPU数量至少为64
# 500G数据量，--parallel设置为8，64个CPU机器大概运行30个小时
bismark_core=$((core/8))

# 获取脚本路径和数据库路径
shell_path=`dirname "$0"`;
methy_pipe_path=${shell_path}/..
script_path=${methy_pipe_path}/script

# 设置参考基因组相关文件路径
genomeDir=${database_dir}/star          # STAR索引目录
genomefa=$database_dir/fasta/genome.fa  # 参考基因组fasta文件
gtf=$database_dir/genes/genes.gtf       # 基因注释文件
genomebed=$database_dir/bed/chr_len.bed # 染色体长度文件
chrom_size_path=$database_dir/bed/chr_nochrM.bed  # 不含线粒体的染色体长度文件
annofile=$database_dir/bed/chr_100kbins_anno.bed  # 100kb bins注释文件
bismark_genome=$database_dir/fasta      # bismark基因组索引目录

# 设置barcode白名单文件
U3CB_methylation=${script_path}/barcodes/U3CB_methylation.txt  # 甲基化文库barcode白名单
cbcsv=${script_path}/barcodes/bUCB3_whitelist.csv             # 甲基化和转录组文库barcode对应关系

# 设置输出路径
exp_outdir=${outdir}/${sample}_exp    # 转录组分析结果目录
methy_dir=${outdir}/${sample}_methy   # 甲基化分析结果目录
# 转录组分析得到的细胞barcode文件，用于甲基化数据的细胞判定
gexcb=$exp_outdir/${sample}/Analysis/step3/filtered_feature_bc_matrix/barcodes.tsv.gz

########################## 数据获取 ##########################
# 判断输入文件是否为OSS云存储路径，如果是则下载到本地
case "${methy_fq1}" in
  oss://*)
    mkdir -p ${outdir}/Rawdata &&\
    ${script_path}/ossutil64 cp ${methy_fq1} ${outdir}/Rawdata/
    ${script_path}/ossutil64 cp ${methy_fq2} ${outdir}/Rawdata/
    ${script_path}/ossutil64 cp ${exp_fq1} ${outdir}/Rawdata/
    ${script_path}/ossutil64 cp ${exp_fq2} ${outdir}/Rawdata/
    exp_fq1_local=${outdir}/Rawdata/$(basename ${exp_fq1})
    exp_fq2_local=${outdir}/Rawdata/$(basename ${exp_fq2})
    methy_fq1_local=${outdir}/Rawdata/$(basename ${methy_fq1})
    methy_fq2_local=${outdir}/Rawdata/$(basename ${methy_fq2})  
    ;;
  *)
    exp_fq1_local=${exp_fq1}
    exp_fq2_local=${exp_fq2}
    methy_fq1_local=${methy_fq1}
    methy_fq2_local=${methy_fq2}
    ;;
esac

########################## 数据质控 (fastp) ##########################
# 使用fastp软件对原始数据进行质控，过滤掉质量低的reads，得到clean data
# fastp参数说明：
# --reads_to_process 0: 处理所有reads
# --cut_tail: 从reads尾部切除低质量碱基
# --cut_tail_window_size 1: 滑动窗口大小为1
# --cut_tail_mean_quality 3: 平均质量阈值
# --disable_adapter_trimming: 禁用adapter修剪
# --length_required 40: 最小reads长度要求
mkdir -p ${outdir}/fastp
exp_fq1_clean=${outdir}/fastp/$(basename ${exp_fq1})
exp_fq2_clean=${outdir}/fastp/$(basename ${exp_fq2})
fastp \
    -i ${exp_fq1_local} \
    -I ${exp_fq2_local} \
    -o ${exp_fq1_clean} \
    -O ${exp_fq2_clean} \
    --reads_to_process 0  --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 --disable_adapter_trimming \
    -j ${outdir}/fastp/$(basename ${exp_fq1})_fastp.json \
    -h ${outdir}/fastp/$(basename ${exp_fq1})_fastp.html \
    --length_required 40 --thread ${core}

# 对甲基化数据进行质控
methy_fq1_clean=${outdir}/fastp/$(basename ${methy_fq1})
methy_fq2_clean=${outdir}/fastp/$(basename ${methy_fq2})
fastp \
    -i ${methy_fq1_local} \
    -I ${methy_fq2_local} \
    -o ${methy_fq1_clean} \
    -O ${methy_fq2_clean} \
    --reads_to_process 0 --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 --disable_adapter_trimming \
    -j ${outdir}/fastp/$(basename ${methy_fq1})_fastp.json \
    -h ${outdir}/fastp/$(basename ${methy_fq1})_fastp.html \
    --length_required 40 --thread ${core}

########################## 转录组文库分析 (SeekSoulTools) ##########################
# 使用SeekSoulTools对转录组文库进行分析
# 后续甲基化文库的细胞基于转录组文库判定的细胞barcode
mkdir -p ${exp_outdir}
seeksoultools rna run \
	--fq1 ${exp_fq1_clean} \
	--fq2 ${exp_fq2_clean} \
	--samplename ${sample} \
	--genomeDir ${genomeDir} \
	--gtf $gtf \
	--chemistry DDV2 \
	--core ${core} \
	--include-introns \
	--outdir ${exp_outdir}

########################## 甲基化文库分析 ##########################
########################## step1: barcode识别和adapter去除 ##########################
# 根据Read1的结构设计和参数对barcode/UMI进行提取和处理
# 主要功能包含三部分：
# 1. 识别和输出包含白名单中barcode的reads
# 2. 识别adapter并去除
# 3. 去掉TN5酶相关的9bp碱基，包括5'端和3'端
# 考虑到reads长度对比对的影响，只输出read1长度大于等于20bp，同时read2长度大于等于60bp的reads
mkdir -p ${methy_dir}
mkdir -p ${methy_dir}/step1
python ${script_path}/barcode_cs_multi.py \
	--fq1 ${methy_fq1_clean} \
	--fq2 ${methy_fq2_clean} \
	--samplename ${sample} \
	--outdir ${methy_dir} \
	--barcode ${U3CB_methylation} \
	--chemistry DD-M \
	--core ${core}

# 对barcode识别后的fastq文件进行质控
fastp \
    -i ${methy_dir}/step1/${sample}_1.fq.gz \
    -I ${methy_dir}/step1/${sample}_2.fq.gz \
    -j ${methy_dir}/step1/${sample}_barcode_fastp.json \
    -h ${methy_dir}/step1/${sample}_barcode_fastp.html \
    --thread ${core} --reads_to_process 0 --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 --disable_adapter_trimming

########################## step2: bismark比对和去重 ##########################
mkdir -p ${methy_dir}/step2/bismark

# 向read name中添加UMI信息，方便后续使用deduplicate_bismark结合UMI信息进行去重
python ${script_path}/add_umi_to_fastq_end_of_read_name.py \
    --fq1 ${methy_dir}/step1/${sample}_1.fq.gz \
    --fq2 ${methy_dir}/step1/${sample}_2.fq.gz \
    --samplename ${sample} \
    --outdir ${methy_dir}/step2/bismark

# bismark比对，实际使用的CPU应该是设置的core的4-5倍
# 参数说明：
# --parallel ${bismark_core}: 并行线程数
# --genome ${bismark_genome}: 参考基因组路径
# -X 1000: 最大插入片段长度
# --temp_dir: 临时文件目录
# --non_directional: 非链特异性文库
bismark \
	--parallel ${bismark_core} \
	--genome ${bismark_genome} \
	-1 ${methy_dir}/step2/bismark/${sample}_1_rename.fq.gz \
	-2 ${methy_dir}/step2/bismark/${sample}_2_rename.fq.gz \
	-o ${methy_dir}/step2/bismark/ \
	-X 1000 \
	--temp_dir ${methy_dir}/step2/bismark/ \
	--non_directional > ${methy_dir}/step2/bismark/bismark.log 2>&1

# deduplicate_bismark要求bam文件按照read名称排序，使用samtools进行排序
samtools sort \
	-@ ${core} \
	-n \
	-o ${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.bam \
	${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe.bam;

# deduplicate_bismark去重，使用barcode和UMI信息进行PCR重复去除
deduplicate_bismark \
	-p \
	--output_dir ${methy_dir}/step2/bismark/ \
	--barcode \
	--bam ${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.bam

# allcools要求bam文件按照染色体位置排序，使用samtools进行排序
samtools sort \
    -@ ${core} \
    -o ${methy_dir}/step2/bismark/${sample}_1_deduplicated_sort.bam \
    ${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.deduplicated.bam

# allcools统计CNN位点的methylation reads和coverage reads
# 输出文件格式：
# 第1列：chromosome
# 第2列：position(1-based)
# 第3列：strand
# 第4列：sequence context
# 第5列：count of reads supporting methylation
# 第6列：read coverage
# 第7列：indicator of significant methylation (1 if no test is performed)
allcools bam-to-allc \
	--bam_path ${methy_dir}/step2/bismark/${sample}_1_deduplicated_sort.bam \
    --convert_bam_strandness \
    --reference_fasta ${genomefa} \
    --output_path ${methy_dir}/step2/bismark/${sample}_allc

########################## step2: 基因注释和覆盖度统计 ##########################
mkdir -p ${methy_dir}/step2/wgs

# 使用featureCounts对bam进行基因注释
# 参数说明：
# -T ${core}: 线程数
# -t gene: 特征类型为gene
# -s 0: 链特异性为0（非链特异性）
# -M -O: 多映射和重叠reads处理
# -g gene_id: 使用gene_id作为基因标识符
# --fracOverlap 0.5: 重叠比例阈值
# -p: 配对末端reads
# --countReadPairs: 计数read pairs
featureCounts \
	-T ${core} -t gene -s 0 -M -O -g gene_id --fracOverlap 0.5 \
	-a ${gtf} \
	-p \
	--countReadPairs \
	-o ${methy_dir}/step2/wgs/wgs_gene_counts.txt \
	-R BAM \
	${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.bam

# 对bam文件按位置排序
samtools sort \
	-@ ${core} \
	-o ${methy_dir}/step2/wgs/${sample}_sort.bam \
	${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.bam

# 使用bedtools统计基因组覆盖度
bedtools genomecov \
	-ibam ${methy_dir}/step2/wgs/${sample}_sort.bam > ${methy_dir}/step2/wgs/bedtools_genomecoverage.txt

# 统计样本在基因组的覆盖度，并将比对信息和覆盖度信息写入json文件中
python ${script_path}/wgs_coverage_depth.py \
	--outdir ${methy_dir}/step2/wgs/ \
	--samplename ${sample} \
	--align_summary ${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_PE_report.txt \
	--summary_json ${methy_dir}/${sample}_summary.json \
	--coveragefile ${methy_dir}/step2/wgs/bedtools_genomecoverage.txt

########################## step3: 细胞判定和单细胞分析 ##########################
# 对featureCounts注释的bam文件按read名称排序
samtools sort \
        -n -O BAM -@ ${core} \
        -o ${methy_dir}/step2/wgs/${sample}_featureCounts_SortByName.bam \
        ${methy_dir}/step2/wgs/${sample}_1_rename_bismark_bt2_pe_sortbyname.bam.featureCounts.bam

# 统计每个barcode含有的基因数、umi数、reads数
# 输出文件：
# counts.xls: 第1列是barcode，第2列是feature，第3列是umi数量，第4列是reads数量
python ${script_path}/step3dnam3.py \
	--bam ${methy_dir}/step2/wgs/${sample}_featureCounts_SortByName.bam \
	--outdir ${methy_dir}/step3

# 甲基化数据根据转录组判定的细胞barcode来过滤得到细胞的umicounts
# 输出文件：
# raw_umicounts.xls: 所有barcode按umi数排序的结果
# filter_umicounts.xls: 判定为细胞的barcode按umi数排序的结果
python ${script_path}/wgs_umi_count_cs.py  \
	--infile ${methy_dir}/step3/counts.xls \
	--rawcsv ${methy_dir}/step3/raw_umicounts.xls \
	--filtercsv ${methy_dir}/step3/filter_umicounts.xls \
	--gexcb ${gexcb} \
	--cbcsv ${cbcsv} \
	--outdir ${methy_dir} \
	--samplename ${sample}

# 将判断为细胞的barcode输出到filtered_barcode文件
cat ${methy_dir}/step3/filter_umicounts.xls |cut -f 1 |grep -v 'Barcode' > ${methy_dir}/step3/filtered_barcode

# 将去重后的bam按照filtered_barcode进行拆分为单个细胞的bam文件
python ${script_path}/split_bams.py \
    --bam ${methy_dir}/step2/bismark/${sample}_1_rename_bismark_bt2_pe_sortbyname.deduplicated.bam \
    --outdir ${methy_dir}/step3 \
	--samplename ${sample} \
	--filtered_barcode ${methy_dir}/step3/filtered_barcode

# 对每个细胞的bam文件，使用samtools进行排序，使用ALLCools统计每个细胞每个CNN位点的甲基化水平
# 输出目录：
# allcools: 每个细胞的甲基化水平文件
# allcools_generate_datasets: 用于下游分析的mcds文件
python ${script_path}/step3_run_allcools_and_generate_datasets.py \
	--indir ${methy_dir}/step3/split_bams \
	--samplename ${sample} \
	--outdir ${methy_dir}/step3/ \
	--genomefa ${genomefa} \
	--chrom_size_path ${chrom_size_path} \
	--filtered_barcode ${methy_dir}/step3/filtered_barcode

# 统计每个细胞的覆盖度
# 输出文件：
# bed目录: 每个细胞的bed文件
# cb_cov目录: 每个细胞的覆盖度文件
# ${sample}_cb_genomecov.xls: 第1列为barcode，第2列为1X深度下的覆盖度
python ${script_path}/wgs_sort_ml_anno_dev.py \
	--filtercb ${methy_dir}/step3/filtered_barcode \
	--core ${bismark_core} \
	--outdir ${methy_dir}/step3/ \
	--samplename ${sample} \
	--genomefa ${genomefa} \
	--allcpath ${methy_dir}/step3/allcools

# 统计每个细胞的CpG数量
# 输出文件：filter_gcov_umi_sort.xls
python ${script_path}/wgs_mk_martix_bins.py \
	--cbfile ${methy_dir}/step3/filtered_barcode \
	--outdir ${methy_dir}/step3 \
	--samplename ${sample} \
	--outcsv ${methy_dir}/step3/${sample}_CPGnum_cb.xls \
	--allcpgfile ${methy_dir}/step2/bismark/${sample}_allc.gz \
	--indir ${methy_dir}/step3/bed \
	--summary_json ${methy_dir}/${sample}_summary.json

# 生成最终的质控表
python ${script_path}/wgs_summary_yf.py \
	--outdir ${methy_dir} \
	--samplename ${sample} \
	--summary_json ${methy_dir}/${sample}_summary.json

########################## step4: 单细胞甲基化数据降维聚类 ##########################
# 利用ALLCools软件进行单细胞甲基化数据降维聚类
# 输出目录step4包含：
# tsne、umap降维聚类图
# h5ad文件
mkdir -p ${methy_dir}/step4
python ${script_path}/step4_allcools_PCA_cluster.py \
	--mcds_path ${methy_dir}/step3/allcools_generate_datasets/${sample}.mcds \
	--samplename ${sample} \
	--var_dim chrom1M \
	--filtered_barcode_file ${methy_dir}/step3/filtered_barcode \
	--outdir ${methy_dir}/step4

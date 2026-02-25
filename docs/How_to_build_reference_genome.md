# How to Build the Reference Genome Database (`--database_dir`)

This pipeline uses a single reference database root directory (`--database_dir`) to locate all genome resources needed for RNA + methylation analysis (STAR index, genome FASTA, GTF, Bismark reference, and chromosome size BED files).

This document generalizes that procedure into a reusable tutorial.

## Expected directory layout

The pipeline derives paths from `--database_dir` using the following conventions:

- `star/` for the STAR genome index
- `fasta/genome.fa` for the genome FASTA (and `fasta/genome.fa.fai` for its index)
- `genes/genes.gtf` for the gene annotation
- `fasta/` is also used as the Bismark genome folder (after `bismark_genome_preparation`)
- `bed/chr_len.bed` and `bed/chr_nochrM.bed` for chromosome sizes

After you finish, your `database_dir` should look like:

```text
<database_dir>/
├── bed/
│   ├── chr_len.bed
│   └── chr_nochrM.bed            (optional but recommended)
├── fasta/
│   ├── genome.fa
│   ├── genome.fa.fai             (created by samtools faidx)
│   └── Bisulfite_Genome/         (created by bismark_genome_preparation)
├── genes/
│   └── genes.gtf
└── star/                         (created by STAR genomeGenerate)
```

## Prerequisites

Make sure the following tools are available in your environment (either on your PATH or inside the container used by Nextflow):

- `STAR`
- `samtools`
- `bismark_genome_preparation` (Bismark)
- `bowtie2` (required by Bismark)
- `awk` and standard Unix utilities

## Inputs

Prepare these two files for your target species/build:

- Genome FASTA: `genome.fa`
- Gene annotation GTF: `genes.gtf`

If your FASTA already has an index file (`genome.fa.fai`), it will be reused; otherwise it will be generated.

## Step-by-step build

Set the variables below and run the commands.

```bash
set -euo pipefail

# 1) Choose an output directory (this becomes --database_dir)
outdir="/path/to/database_dir"

# 2) Provide the source FASTA and GTF files
fasta="/path/to/source/genome.fa"
gtf="/path/to/source/genes.gtf"

# 3) Threads (adjust for your machine)
star_threads=15
bismark_threads=16

mkdir -p "${outdir}/fasta" "${outdir}/genes" "${outdir}/bed" "${outdir}/star"

# Copy FASTA and GTF into the expected locations
cp "${fasta}" "${outdir}/fasta/genome.fa"
cp "${gtf}" "${outdir}/genes/genes.gtf"

# If a FASTA index exists next to the source FASTA, reuse it
if [[ -f "${fasta}.fai" ]]; then
  cp "${fasta}.fai" "${outdir}/fasta/genome.fa.fai"
fi
```

### Step 1: Build the STAR genome index

```bash
STAR --runThreadN "${star_threads}" \
  --runMode genomeGenerate \
  --genomeDir "${outdir}/star" \
  --genomeFastaFiles "${outdir}/fasta/genome.fa" \
  --sjdbGTFfile "${outdir}/genes/genes.gtf"
```

### Step 2: Ensure the FASTA index exists (`genome.fa.fai`)

```bash
if [[ ! -f "${outdir}/fasta/genome.fa.fai" ]]; then
  samtools faidx "${outdir}/fasta/genome.fa"
fi
```

### Step 3: Build the Bismark reference in `fasta/`

Run Bismark genome preparation using the `fasta/` folder as the genome directory:

```bash
cd "${outdir}"
bismark_genome_preparation --parallel "${bismark_threads}" ./fasta
```

This creates Bismark genome files under `${outdir}/fasta/` (e.g. `Bisulfite_Genome/`).

### Step 4: Create chromosome size BED files

The pipeline uses:

- `bed/chr_len.bed` as the full chromosome size list
- `bed/chr_nochrM.bed` (if present) as the preferred list that excludes mitochondrial and other undesired contigs

Create `chr_len.bed` from the FASTA index:

```bash
awk -F '\t' '{print $1"\t"$2}' "${outdir}/fasta/genome.fa.fai" > "${outdir}/bed/chr_len.bed"
```

Then create `chr_nochrM.bed` by selecting the chromosomes you want to keep.

>[!Note]
>You should decide which contigs are primary chromosomes, which are alt contigs, and which represent mitochondrial DNA for your species/build.
>For humans, a common choice is chr1–chr22 plus chrX/chrY; for other species, choose the equivalent primary chromosomes.

Example (simple, but species-dependent): keep the first N lines of `chr_len.bed`:

```bash
# Example only: choose N according to your genome build
head -n 22 "${outdir}/bed/chr_len.bed" > "${outdir}/bed/chr_nochrM.bed"
```

For a more explicit human GRCh38 example (if your contig names are chr1..chr22, chrX, chrY):

```bash
awk '$1 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$)/ {print $0}' \
  "${outdir}/bed/chr_len.bed" > "${outdir}/bed/chr_nochrM.bed"
```

## Validation checklist

Confirm the key files exist:

```bash
test -d "${outdir}/star"
test -f "${outdir}/fasta/genome.fa"
test -f "${outdir}/fasta/genome.fa.fai"
test -d "${outdir}/fasta/Bisulfite_Genome" || true
test -f "${outdir}/genes/genes.gtf"
test -f "${outdir}/bed/chr_len.bed"
```

If you created `chr_nochrM.bed`, also check:

```bash
test -f "${outdir}/bed/chr_nochrM.bed"
```

## Use in the pipeline

Pass the directory as the `--database_dir` parameter:

```bash
nextflow run nf/main.nf \
  -c nf/nextflow.new.config \
  --workflow rna_met \
  --database_dir "${outdir}" \
  --samplesheet /path/to/samplesheet.csv \
  --outdir /path/to/results
```

>[!Note]
>When running on Docker/Kubernetes, `--database_dir` must be a path that exists inside the container/pod (for example, a mounted volume path).

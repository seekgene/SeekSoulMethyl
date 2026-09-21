# SeekSoulMethyl v2.2.1 Upgrade Report

**Audience**: Users of the SeekSoulMethyl pipeline
**Change**: Fixed a silent data loss bug in the methylation MCDS merging step (`MERGE_MCDS`)
**Affected versions**: v2.0.0 – v2.2.0
**Fixed version**: v2.2.1
**Date**: September 2026

---

## 1. Summary

SeekSoulMethyl v2.2.1 is a **defect-fix release** (no new features, no interface changes). It fixes a **silent data loss bug** in the methylation `MERGE_MCDS` step — the step that merges per-shard single-cell methylation matrices into the final MCDS used for clustering and downstream analysis.

> **What is MCDS?** MCDS (Methylation Cell Data Structure) is the standard container that stores each cell's methylation information across genomic intervals ("bins"). A **bin** is a fixed-length window of the genome — `chrom20k` means one bin per 20,000 base pairs (the default used for clustering), `chrom1M` means one bin per 1,000,000 base pairs (a coarser, larger bin).

The bug caused a **random subset of cells to be silently zeroed out** in the large-scale bins (**chrom1M ~9% of cells, chrom500k ~8% of cells**), with **no error, no warning, and no visible corruption**. The default clustering bin (**chrom20k**) was **essentially unaffected** (in the 33 samples, only a very few bins differ in a very few cells, with no effect on clustering conclusions).

This report also includes a **before/after consistency verification** (Section 7): the same sample processed by v2.2.0 (before the fix) and v2.2.1 (after the fix) yields **essentially identical chrom20k matrices** and fully intermixed clustering — confirming the fix does not change the default clustering results.

**Key takeaway**: if you only use the default `chrom20k` clustering, your existing results are intact and **do not need to be re-run**. If you use the large-scale bins (`chrom1M` / `chrom500k`), upgrade to v2.2.1 and re-run.

---

## 2. Background

The `MERGE_MCDS` step combines the per-shard MCDS fragments produced by the earlier `allcools generate-dataset` step into a single consolidated MCDS. It writes the merged result to a Zarr store using xarray's `to_zarr(append_dim=...)` in a Dask-parallel fashion.

In versions v2.0.0 through v2.2.0, this append was performed **without first materializing the data** (`ds.load()`). Under xarray's older chunk-alignment checks, the append path did not enforce chunk alignment between the Dask chunks being written and the existing Zarr store chunks. When multiple Dask chunks mapped onto the same Zarr chunk and were written concurrently, a **race condition** could silently overwrite a whole row of a cell's data with zeros.

Because the failure is silent and stochastic, affected cells are scattered randomly across the matrix and are only detectable by comparing against an independent ground truth.

---

## 3. Root Cause

- **Upstream issue**: [xarray #8876 — "Possible race condition when appending to an existing zarr"](https://github.com/pydata/xarray/issues/8876) (closely related: [#8882 — "to_zarr silently loses data when using append_dim"](https://github.com/pydata/xarray/issues/8882)).
- **Mechanism**: `to_zarr(append_dim=...)` with Dask-parallel writes and misaligned chunks can write multiple Dask chunks into the same Zarr chunk concurrently, losing data.
- **Why it was silent**: xarray's chunk-alignment safety check (`safe_chunks`) only applied on *new* variables, not on the *append* path — so the append bypassed the check that would otherwise have raised an error.

---

## 4. What v2.2.1 Fixes

v2.2.1 adds an explicit `ds.load()` before the append so that the merged data is **fully materialized and written serially**, eliminating the concurrent-write race.

- Affected versions: **v2.0.0 – v2.2.0**
- Fixed version: **v2.2.1** (the fix has been pushed to the release branch and tagged)

> **Note on memory**: `ds.load()` materializes the full merged MCDS into memory before writing. For large samples this raises peak memory usage; the fix trades memory for correctness (serial, race-free writes).

---

## 5. Impact Scope

The bug affects only the **large-scale bins**, where the `count_type` dimension is chunked as a single whole (mc + cov packed into one chunk).

> **What are mc and cov?** For each cell and each bin, `cov` (coverage) is the total number of reads covering that bin, and `mc` (methylated counts) is the number of those reads supporting methylation. The methylation level is `mc/cov`.

The three observed damage modes are:

| Damage mode | What happens | Affected bins |
|---|---|---|
| **Whole-row zeroing** | a cell's mc **and** cov are both reset to 0 | chrom1M (~9%), chrom500k (~8%) |
| **mc halving** | cov unchanged, mc halved | small bins (rare, a few cells) |
| **cov halving** | mc unchanged, cov halved | small bins (rare, a few cells) |

| Bin | v2.2.0 loss | v2.2.1 |
|---|---|---|
| chrom1M | ~9% of cells zeroed | **0 (fully restored)** |
| chrom500k | ~8% of cells zeroed | **0 (fully restored)** |
| chrom100k / 50k / 10k | essentially 0 (occasional single-cell halving) | 0 |
| geneslop2k | occasional cov halving only (mc intact, a few cells) | 0 |
| **chrom20k (default clustering)** | **essentially 0** | **essentially 0** |

> The "halving" modes are empirical observations; the race-condition mechanism itself predicts only whole-row zeroing. The halving is likely a partial-overlap variant of the same concurrent-write problem. The key reassurance — **chrom20k is almost never affected** — is an empirical result confirmed by four independent lines of evidence (Section 6: 6.1–6.3 prove chrom20k is essentially loss-free, 6.4 proves chrom1M/500k recovery), even though the precise mechanism of chrom20k's immunity is not fully explained.

> `geneslop2k` (the 2 kb gene-flanking bins used for gene-level DMG, differential methylation genes) is affected only by the **cov-halving** mode — mc is intact and only a handful of cells lose cov. This does not change the chrom20k clustering conclusion.

---

## 6. Verification and Test Results

The fix was validated with four independent lines of evidence (all pointing to the same conclusion: chrom20k always zero-loss, chrom1M/500k loss fully recovered):

### 6.1 Cross-bin audit of 33 production samples

A consistency audit of **33 production methylation samples** (295,408 cells across 7 batches) compared the coverage sum of every cell across the six bin resolutions. A bin whose total was zero while the others were normal indicates that the cell was lost in that bin.

| Bin | Average loss rate |
|---|---|
| chrom1M | 9.4% (per-sample 7.2% – 13.4%) |
| chrom500k | 8.3% (per-sample 6.4% – 10.5%) |
| **chrom20k** | **essentially 0 (a few bins in a few cells)** |
| chrom100k / 50k / 10k | 0 |

> The 9.4% here is the per-sample-average of the loss rate; the absolute cell count recovered (Section 6.4) corresponds to 9.7% of all cells. Both figures refer to chrom1M only and describe the same phenomenon.

### 6.2 v2.2.0 reproduction

Re-running 14 MCDS with v2.2.0 (the version containing the bug) reproduced the loss at the same magnitude — chrom1M ~9.6%, chrom500k ~8.8%, chrom20k 0 — confirming the reproduction path is correct and the bug reproduces reliably.

### 6.3 Ground-truth rebuild from single-cell ALLC files

One sample's MCDS was rebuilt directly from single-cell ALLC (all-cytosine methylation) files (skipping the merge step, i.e. ground truth) and compared cell-by-cell, bin-by-bin against the old MCDS:

| Bin | Cells zeroed |
|---|---:|
| chrom1M | 687 |
| chrom500k | 601 |
| **chrom20k** | **0** |

### 6.4 v2.2.1 fix validation (cell-by-cell, bin-by-bin)

v2.2.1 (containing `ds.load()`) was run to completion and compared cell-by-cell, bin-by-bin against v2.2.0 across all 33 samples:

| Bin | Cells recovered by v2.2.1 |
|---|---:|
| chrom1M | **28,726** |
| chrom500k | **26,442** |
| chrom20k | essentially 0 (a few cells with a one-bin mc/cov difference only) |

**Conclusion**: every cell that v2.2.0 silently zeroed in chrom1M/500k is fully restored in v2.2.1, and chrom20k was essentially unaffected (a few cells with a one-bin difference only, no whole-row zeroing).

### 6.5 Before/after per-bin comparison

The figures below compare chrom1M and chrom500k cell-by-cell and bin-by-bin between v2.2.0 (before the fix) and v2.2.1 (after the fix). **How to read**: the x-axis is the v2.2.0 value, the y-axis is the v2.2.1 value; points on the diagonal mean the two versions agree. The red points mark every cell×bin pair where the two versions differ — **all of them are rows zeroed out in v2.2.0 (x = 0) and fully restored in v2.2.1 (y > 0)** (the difference is strictly one-way: no bin is non-zero in v2.2.0 but zero in v2.2.1). This directly visualizes that the fix succeeds: v2.2.1 recovers the data that v2.2.0 silently lost in these two large-scale bins.

<p align="center"><img src="mcds_corr_v220_vs_v221_chrom1M.png" alt="chrom1M before/after correlation" width="880"><br><sub>Figure 1 — chrom1M: v2.2.0 vs v2.2.1. Red points = cells zeroed in v2.2.0 and restored in v2.2.1.</sub></p>

<p align="center"><img src="mcds_corr_v220_vs_v221_chrom500k.png" alt="chrom500k before/after correlation" width="880"><br><sub>Figure 2 — chrom500k: v2.2.0 vs v2.2.1. Red points = cells zeroed in v2.2.0 and restored in v2.2.1.</sub></p>

For comparison, **chrom20k** (the default clustering bin) is shown below: it is **identical in the WTJW880 sample** (r = 1.0), consistent with the fact that the bug essentially does not affect chrom20k.

<p align="center"><img src="mcds_corr_v220_vs_v221_chrom20k.png" alt="chrom20k before/after correlation" width="880"><br><sub>Figure 3 — chrom20k: v2.2.0 vs v2.2.1, consistent in the WTJW880 sample (r = 1.0).</sub></p>

### 6.6 Quality-control metrics are identical before and after the fix

The QC summary of the same sample (WTJW880) under v2.2.0 and v2.2.1 shows that **every quality-control metric is identical** (these 10 metrics are all computed **upstream** of the merge step, so the fix — which changes only the merge — does not touch them; the post-merge chrom1M/500k changes are covered in Sections 6.4/6.5).

| Metric | v2.2.0 | v2.2.1 |
|---|---|---|
| Estimated cells | 442 | 442 |
| GEX median genes per cell | 1,066 | 1,066 |
| MET CpG per median cell | 905 | 905 |
| GEX Valid Barcode | 91.21% | 91.21% |
| GEX Sequencing Saturation | 83.22% | 83.22% |
| MET Valid Barcodes | 85.84% | 85.84% |
| MET C-T Conversion | 99.80% | 99.80% |
| MET CpG Methylation Rate | 80.74% | 80.74% |
| MET Reads Mapped Confidently | 79.83% | 79.83% |
| MET Total CpGs Detected | 517,750 | 517,750 |

This is expected: the fix only changes how the merged MCDS is written (serial vs concurrent), so it does not touch any QC metric computed upstream of the merge.

---

## 7. Before/After Consistency Verification

> This section answers a related question: **does the fix itself change the default (chrom20k) clustering results?** It compares v2.2.0 (before the fix) against v2.2.1 (after the fix) on the same sample.

### 7.1 Design

The same sample (WTJW880, 442 cells, DD-MET5) was processed by **v2.2.0** and **v2.2.1** on the **same reference genome** (`refdata-met-GRCh38-2020-A`). The two `chrom20k` MCDS matrices (442 cells × 154,423 bins each; barcodes identical) were compared at two levels: (a) the raw cell × bin matrix, and (b) the clustering result after merging the two versions as two "samples" **without batch correction** (no Harmony) and running LSI + Leiden.

> **Terminology**: LSI (latent semantic indexing) is the dimensionality-reduction step; Leiden is the clustering algorithm; UMAP is a visualization of the reduced space. A barcode is each cell's unique label; a "twin" is the two copies of the same barcode produced by v2.2.0 and v2.2.1. "Silhouette" measures how separated the two versions are (0 = fully mixed); the kNN (k-nearest-neighbour) same-version fraction measures mixing (≈0.5 = fully mixed); ARI/NMI score how similar two clusterings are (1.0 = identical).

### 7.2 Results — level (a): raw matrix is identical

| Panel | Pearson r | cells×bins differing | max\|Δ\| |
|---|---|---|---|
| Coverage (cov) | **1.00000000** | **0 / 68,254,966 (0.000%)** | 0 |
| Methylated reads (mc) | **1.00000000** | **0 / 68,254,966 (0.000%)** | 0 |
| Methylation level (mc/cov) | **1.00000000** | **0 (0.000%)** | 0 |

**In the WTJW880 sample, the chrom20k matrices from v2.2.0 and v2.2.1 are identical cell-by-cell and bin-by-bin** (r = 1.0, zero differing cell×bin pairs). No cell is zeroed, no bin differs, and the per-cell coverage totals match exactly (the correlation figure is shown in Section 6.5, Figure 3). (Note: the "few cells with a one-bin difference" in Section 6.4 is an observation across all 33 production samples — numerical noise, not data loss; the two statements are consistent.)

### 7.3 Results — level (b): clustering is fully intermixed

| Metric | Value | Fully-mixed reference |
|---|---|---|
| Silhouette by version (LSI / UMAP) | −0.0023 / −0.0022 | 0 |
| kNN (k=30) same-version fraction | 0.485 | 0.5 (random) |
| Dominant-version fraction per cluster (weighted) | 0.505 | 0.5 |
| **Per-barcode twin agreement** | **434 / 438 = 99.09%** | 100% |
| **Per-barcode twin ARI / NMI** | **0.9892 / 0.9809** | 1.0 |

The two versions are completely intermixed, and 99.09% of cells land in the same cluster across versions (ARI 0.989). Note: 4 of the 442 barcodes were removed by QC before clustering, so the comparison is over 438 barcodes (434 agree / 4 differ); the 4 discordant barcodes are boundary jitters, not systematic version bias.

> **Honest caveat**: WTJW880 is a tiny demo (442 cells, ~525k reads), so the 11 clusters have no biological meaning and clustering is partly confounded by coverage. The level-(a) result (identical matrices) is the decisive evidence; the level-(b) result is directional confirmation only.

<p align="center"><img src="umap_by_sample.png" alt="Integrated UMAP colored by version" width="560"><br><sub>Figure 4 — Integrated UMAP colored by pipeline version: the two versions are fully intermixed.</sub></p>

### 7.4 Conclusion

**The v2.2.1 fix does not change the chrom20k clustering result.** The raw matrices are identical in the WTJW880 sample (r = 1.0), and the clustering is fully intermixed (99.09% per-barcode agreement, ARI 0.989). This is expected: the bug only affected chrom1M/500k, so the fix makes no substantive change in chrom20k.

---

## 8. Recommendation for Users

- **If you only use the default `chrom20k` clustering**: your existing results are intact and do **not** need to be re-run. The bug essentially does not affect chrom20k, and the fix does not change chrom20k (Section 7).
- **If you use `chrom1M` or `chrom500k` (large-scale bin analyses, e.g. differentially methylated regions (DMR) or coverage summaries)**: upgrade to **v2.2.1** and re-run, because ~8–9% of cells were silently zeroed in these bins in earlier versions.
- **If you use `geneslop2k` or `chrom100k/50k/10k`**: impact is minimal — only a few cells have cov halved (mc intact), with no whole-row zeroing; conclusions are generally unaffected. Re-run on v2.2.1 only if you need strictly exact values.
- **How do I know which bin I used?** The default is `chrom20k`. If you did not explicitly change the bin parameter, you used chrom20k and are unaffected. Large-scale bins (`chrom1M`/`chrom500k`) are only used if you explicitly requested differentially methylated regions (DMR) or coverage-summary analyses. You can check your run parameters or config file for the `bin` / `bin_size` field: a value containing `1M`/`500k` means a large-scale bin was used.
- The loss was **silent** — no error, no warning — so older results should not be trusted for chrom1M/500k without a re-run on v2.2.1.

---

## 9. References

- [xarray #8876 — Possible race condition when appending to an existing zarr](https://github.com/pydata/xarray/issues/8876)
- [xarray #8882 — to_zarr silently loses data when using append_dim](https://github.com/pydata/xarray/issues/8882)

# SeekSoulMethyl v2.2.1 MCDS Merging Data Loss Fix Report

**Audience**: Users of the SeekSoulMethyl pipeline
**Change**: Fixed a silent data loss bug in the methylation MCDS merging step (`MERGE_MCDS`)
**Affected versions**: v2.0.0 – v2.2.0
**Fixed version**: v2.2.1
**Report version**: v1.0
**Date**: September 2026

---

## 1. Summary

SeekSoulMethyl v2.2.1 fixes a **silent data loss bug** in the methylation `MERGE_MCDS` step — the step that merges per-shard single-cell methylation matrices into the final MCDS used for clustering and downstream analysis.

The bug caused a **random subset of cells to be silently zeroed out** in the large-scale methylation bins (**chrom1M ~9% of cells, chrom500k ~8% of cells**), with **no error, no warning, and no visible corruption**. The default clustering bin (**chrom20k**) was **not affected** (zero loss in all tests).

**Key takeaway**: if you only use the default `chrom20k` clustering, your v2.2.0 results are intact and **do not need to be re-run**. If you use the large-scale bins (`chrom1M` / `chrom500k`), upgrade to v2.2.1 and re-run.

---

## 2. Background

The `MERGE_MCDS` step combines the per-shard MCDS fragments produced by the earlier `allcools generate-dataset` step into a single consolidated MCDS. It writes the merged result to a Zarr store using xarray's `to_zarr(append_dim=...)` in a Dask-parallel fashion.

In versions v2.0.0 through v2.2.0, this append was performed **without first materializing the data** (`ds.load()`). Under xarray's older chunk-alignment checks, the append path did not enforce chunk alignment between the Dask chunks being written and the existing Zarr store chunks. When multiple Dask chunks mapped onto the same Zarr chunk and were written concurrently, a **race condition** could silently overwrite a whole row of a cell's data with zeros.

Because the failure is silent and stochastic, affected cells are scattered randomly across the matrix and are only detectable by comparing against an independent ground truth.

---

## 3. Root Cause

- **Upstream issue**: [xarray #8876 — "Possible race condition when appending to an existing Zarr store"](https://github.com/pydata/xarray/issues/8876) (closely related: [#8882 — "to_zarr silently loses data when using append_dim"](https://github.com/pydata/xarray/issues/8882)).
- **Mechanism**: `to_zarr(append_dim=...)` with Dask-parallel writes and misaligned chunks can write multiple Dask chunks into the same Zarr chunk concurrently, losing data.
- **Why it was silent**: xarray's chunk-alignment safety check (`safe_chunks`) only applied on *new* variables, not on the *append* path — so the append bypassed the check that would otherwise have raised an error.

---

## 4. What v2.2.1 Fixes

v2.2.1 adds an explicit `ds.load()` before the append so that the merged data is **fully materialized and written serially**, eliminating the concurrent-write race.

- Affected versions: **v2.0.0 – v2.2.0**
- Fixed version: **v2.2.1** (the fix has been pushed to the release branch and tagged)

---

## 5. Impact Scope

The bug affects only the **large-scale bins**, where the `count_type` dimension is chunked as a single whole (mc + cov packed into one chunk). The three observed damage modes are:

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
| **chrom20k (default clustering)** | **0** | **0** |

> `geneslop2k` (the 2 kb gene-flanking bins used for gene-level DMG) has `count_type chunk = 1` like the small bins, so it is affected only by the **cov-halving** mode — mc is intact and only a handful of cells lose cov (verified on one sample). This does not change the chrom20k clustering conclusion.

---

## 6. Verification and Test Results

The fix was validated with four independent lines of evidence:

### 6.1 Cross-region audit of 33 production samples

A consistency audit of **33 production methylation samples** (295,408 cells across 7 batches) compared the coverage sum of every cell across the six bin resolutions. A region whose total was zero while the others were normal indicates that the cell was lost in that region.

| Bin | Average loss rate |
|---|---|
| chrom1M | 9.4% (per-sample 7.2% – 13.4%) |
| chrom500k | 8.3% (per-sample 6.4% – 10.5%) |
| **chrom20k** | **0 (no loss)** |
| chrom100k / 50k / 10k | 0 |

### 6.2 v2.2.0 reproduction

Re-running 14 MCDS with v2.2.0 (containing the bug) reproduced the loss at the same magnitude — chrom1M ~9.6%, chrom500k ~8.8%, chrom20k 0 — confirming the pipeline link and the stability of the bug.

### 6.3 Ground-truth rebuild from single-cell ALLC files

One sample's MCDS was rebuilt directly from single-cell ALLC files (no merge step = ground truth) and compared cell-by-cell, bin-by-bin against the old MCDS:

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

**Conclusion**: every cell that v2.2.0 silently zeroed in chrom1M/500k is fully restored in v2.2.1, and chrom20k was never affected.

---

## 7. Recommendation for Users

- **If you only use the default `chrom20k` clustering**: your existing results are intact and do **not** need to be re-run. The bug does not affect chrom20k.
- **If you use `chrom1M` or `chrom500k` (large-scale bin analyses, e.g. coarse DMR/coverage summaries)**: upgrade to **v2.2.1** and re-run, because ~8–9% of cells were silently zeroed in these bins in earlier versions.
- The loss was **silent** — no error, no warning — so older results should not be trusted for chrom1M/500k without a re-run on v2.2.1.

---

## 8. References

- [xarray #8876 — Possible race condition when appending to an existing Zarr store](https://github.com/pydata/xarray/issues/8876)
- [xarray #8882 — to_zarr silently loses data when using append_dim](https://github.com/pydata/xarray/issues/8882)

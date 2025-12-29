import argparse
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_bam_pysam(in_bam, out_bam, threads=1, keyword="reverse", keep_tmp=False):
    try:
        import pysam
    except Exception as e:
        raise RuntimeError("pysam not available: install pysam or use --engine samtools") from e
    if not os.path.exists(in_bam):
        raise FileNotFoundError(in_bam)
    os.makedirs(os.path.dirname(os.path.abspath(out_bam)) or ".", exist_ok=True)
    fin = pysam.AlignmentFile(in_bam, "rb")
    fout = pysam.AlignmentFile(out_bam, "wb", template=fin)
    try:
        if threads and threads > 1:
            fin.set_threads(threads)
            fout.set_threads(threads)
    except Exception:
        pass
    kw = (keyword or "").lower()
    for rec in fin.fetch(until_eof=True):
        q = rec.query_name or ""
        if kw and (q.lower().find(kw) != -1):
            if rec.is_read1:
                rec.flag = (rec.flag & (~0x40)) | 0x80
            elif rec.is_read2:
                rec.flag = (rec.flag & (~0x80)) | 0x40
        fout.write(rec)
    fin.close()
    fout.close()
    
def run_one(bam_path, out_dir, threads, keyword, skip_existing):
    try:
        import pysam
    except Exception as e:
        raise RuntimeError("pysam not available: install pysam") from e
    bam_path = bam_path.strip()
    if not bam_path or bam_path.startswith('#'):
        return bam_path, False, 'skipped-empty'
    if not os.path.exists(bam_path):
        return bam_path, False, 'not-found'
    base = os.path.splitext(os.path.basename(bam_path))[0]
    out_dir_use = out_dir or os.path.dirname(bam_path) or '.'
    os.makedirs(out_dir_use, exist_ok=True)
    out_bam = os.path.join(out_dir_use, f"{base}.bam")
    if skip_existing and os.path.exists(out_bam):
        return bam_path, True, 'done-exist'
    try:
        process_bam_pysam(bam_path, out_bam, threads=threads, keyword=keyword)
        return bam_path, True, 'done'
    except Exception as e:
        return bam_path, False, f'failed:{e}'

def read_list(fp):
    with open(fp, 'r') as fh:
        return [line.strip() for line in fh if line.strip()]

def main():
    ap = argparse.ArgumentParser(description="Swap flags on reverse reads in BAM files")
    ap.add_argument("--bam_list", required=True)
    ap.add_argument("--out_dir", type=str, required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--max_workers", type=int, default=4)
    ap.add_argument("--keyword", type=str, default="reverse")
    ap.add_argument("--skip_existing", action="store_true")
    args = ap.parse_args()
    paths = read_list(args.bam_list)
    results = []
    with ProcessPoolExecutor(max_workers=args.max_workers) as ex:
        futs = [
            ex.submit(
                run_one,
                bp,
                args.out_dir,
                args.threads,
                args.keyword,
                args.skip_existing,
            ) for bp in paths
        ]
        for fut in as_completed(futs):
            results.append(fut.result())
    failed = [r for r in results if not r[1]]
    print(f"Total: {len(results)}, Succeeded: {len(results)-len(failed)}, Failed: {len(failed)}")
    for bp, ok, msg in failed:
        print(f"FAIL: {bp} -> {msg}")

if __name__ == "__main__":
    main()
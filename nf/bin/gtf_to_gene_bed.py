#!/usr/bin/env python3

import sys
import os
import re
import gzip
import click


def open_maybe_gz(path: str):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


_ATTR_RE = re.compile(r"([A-Za-z0-9_]+)\s+\"([^\"]*)\"\s*;")


def parse_attributes(attr_text: str):
    attrs = {}
    for m in _ATTR_RE.finditer(attr_text):
        key, val = m.group(1), m.group(2)
        attrs[key] = val
    return attrs


@click.command()
@click.option('--gtf', required=True, type=click.Path(exists=True), help='gtf path')
@click.option('--out', required=True, type=click.Path(), help='gene.bed path')
def main(gtf: str, out: str):
    """convert gtf to bed4 format
    Output:
        gene.bed: BED4 format file, 4 columns: chrom, start(0-based), end, name
    """

    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    kept = 0
    total = 0
    with open_maybe_gz(gtf) as fh, open(out, 'w') as out_fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            total += 1
            chrom, source, feature, start, end, score, strand, frame, attrs_raw = parts
            if feature != 'gene':
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except Exception:
                continue

            # Convert to BED 0-based start, end stays 1-based inclusive -> BED end is exclusive
            bed_start = max(0, start_i - 1)
            bed_end = end_i  # BED expects end as exclusive; GTF end is inclusive

            attrs = parse_attributes(attrs_raw)
            gene_id = attrs.get('gene_id')
            gene_name = attrs.get('gene_name')
            if gene_id and gene_name:
                name = f"{gene_id}_{gene_name}"
            elif gene_id:
                name = gene_id
            elif gene_name:
                name = gene_name
            else:
                name = '.'

            out_fh.write(f"{chrom}\t{bed_start}\t{bed_end}\t{name}\n")
            kept += 1

    click.echo(f'total lines: {total}, gene lines: {kept}, output: {out}')


if __name__ == '__main__':
    main()
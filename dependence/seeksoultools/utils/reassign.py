from collections import defaultdict
import csv
import multiprocessing as mp
from tqdm import tqdm
import pandas as pd
import click
def update_cdr3_file(input_cdr3, output_cdr3, result_dict, buffer_size=16*1024**2):
    with open(input_cdr3, 'r', buffering=buffer_size) as fin, \
         open(output_cdr3, 'w', buffering=buffer_size) as fout:
        output_buffer = []
        buffer_limit = 10000  
        
        for line in fin:
            fields = line.strip().split('\t')
            if len(fields) >= 11:
                contig = fields[0]
                if contig in result_dict:
                    fields[10] = str(result_dict[contig])
                    output_buffer.append('\t'.join(fields) + '\n')
                    if len(output_buffer) >= buffer_limit:
                        fout.writelines(output_buffer)
                        output_buffer = []
        if output_buffer:
            fout.writelines(output_buffer)

def process_chunk(reads_chunk):
    barcode_umi_dict = defaultdict(lambda: defaultdict(int))
    
    for read in reads_chunk:
        read_name, contig = read
        parts = read_name.split('_')
        barcode = parts[0]
        umi = parts[1]
        barcode_umi = f"{barcode}_{umi}"
        barcode_umi_dict[barcode_umi][contig] += 1
    return dict(barcode_umi_dict)

def process_assign_file(input_file, num_processes=16, chunk_size=1000000):
    global_barcode_umi_dict = defaultdict(lambda: defaultdict(int))
    pool = mp.Pool(processes=num_processes)
    def chunk_reader():
        chunk = []
        with open(input_file, 'r') as f:
            for line in f:
                read_name, contig = line.strip().split('\t')
                chunk.append((read_name, contig))
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
            if chunk:
                yield chunk

    total_chunks = sum(1 for _ in chunk_reader())
    
    for chunk_result in pool.imap_unordered(process_chunk, chunk_reader()):
        for barcode_umi, contig_counts in chunk_result.items():
            for contig, count in contig_counts.items():
                global_barcode_umi_dict[barcode_umi][contig] += count

    pool.close()
    pool.join()

    final_contig_barcodes = defaultdict(set)
    
    for barcode_umi, contig_counts in global_barcode_umi_dict.items():
        sorted_counts = sorted(contig_counts.items(), key=lambda x: x[1], reverse=True)
        if len(sorted_counts) == 1:
            final_contig_barcodes[sorted_counts[0][0]].add(barcode_umi)
        else:
            ratio = sorted_counts[1][1] / sorted_counts[0][1]
            if ratio < 0.15:
                final_contig_barcodes[sorted_counts[0][0]].add(barcode_umi)
            else:
                for contig, _ in sorted_counts:
                    final_contig_barcodes[contig].add(barcode_umi)


    df = pd.DataFrame({
        'barcode_umis': [list(v) for v in final_contig_barcodes.values()]
        }, index=final_contig_barcodes.keys())



    final_counts = {contig: len(barcodes) for contig, barcodes in final_contig_barcodes.items()}
    return final_counts

##@click.command(context_settings=dict(help_option_names=['-h', '--help']))
##@click.option('--assign', help="sample_assign.out")
##@click.option('--cdr3', help="sample_cdr3.out")
##@click.option('--new', help="new.sample_cdr3.out")
##@click.option('--threads', default=None, type=int)
def main(assign, cdr3, new, threads):
    if not threads:
        num_processes = mp.cpu_count()
    else:
        num_processes = threads
    result = process_assign_file(assign, num_processes)
    update_cdr3_file(cdr3, new, result, buffer_size=16*1024**2)
##if __name__ == '__main__':
##    main()

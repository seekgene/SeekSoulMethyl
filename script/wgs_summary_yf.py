import json
import os
import click


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option("--outdir", help="outdir.")
@click.option("--samplename", help="samplename.")
@click.option("--summary_json", help="summary json.")
def outcsv(outdir, samplename, summary_json):
    os.makedirs(outdir, exist_ok = True)
    with open(summary_json,'r') as fh:
        summary = json.load(fh)
    
    rawreads = summary["stat"]["total"]
    vaildreads = summary["stat"]["valid"]
    vaildratio = vaildreads / rawreads
    rate_7fratio = float(summary["stat"]["rate_7f"])
    rate_17lmeratio = float(summary["stat"]["rate_17lme"])
    rate_7f17lmeratio = float(summary["stat"]["rate_7f17lme"])
    conversion = float(summary["stat"]["ct_mean"])
    cc_ratio = float(summary["stat"]["cc_mean"])

    raw = f'{summary["stat"]["total"]}'
    vaild = f'{vaildratio:.2%}'
    rate_7f = f'{rate_7fratio:.2%}'
    rate_17lme = f'{rate_17lmeratio:.2%}'
    rate_7f17lme = f'{rate_7f17lmeratio:.2%}'

    mapgenome = f'{summary["mapping"]["Reads Mapped to Genome"]:.2%}'
    confidently = f'{summary["mapping"]["Reads Mapped Confidently to Genome"]:.2%}'
    coverage = f'{summary["coverage"]["Genome Coverage rate"]:.2%}'
    total_cpgs = f'{summary["cells"]["Total CPGs Detected"]}'
    cov_max_cell = f'{summary["cells"]["Genome Coverage rate of max cell"]:.2%}'
    cpgs_max_cell = f'{summary["cells"]["CPGs of max cell"]}'
    umi_max_cell = f'{summary["cells"]["UMIs of max cell"]}'
    reads_max_cell = f'{summary["cells"]["Reads of max cell"]}'
    saturation_max_cell = f'{summary["cells"]["Saturation of max cell"]:.2%}'
    cov_median_cell = f'{summary["cells"]["Genome Coverage rate of median cell"]:.2%}'
    cpgs_median_cell = f'{summary["cells"]["CPGs of median cell"]}'
    umi_median_cell = f'{summary["cells"]["UMIs of median cell"]}'
    reads_median_cell = f'{summary["cells"]["Reads of median cell"]}'
    saturation_median_cell = f'{summary["cells"]["Saturation of median cell"]:.2%}'
    cellnum = f'{summary["cells"]["Estimated Number of Cells"]}'
    fraction = f'{summary["cells"]["Fraction Reads in Cells"]:.2%}'
    


    header=('Samplename,Estimated_Number_of_Cells,Number_of_Reads,valid_barcode_ratio,vaild_7F_reads_rate,vaild_17LME_reads_rate,vaild_7F17LME_reads_rate,C-T_conversion,C-C_ratio,'
            'Reads_Mapped_to_Genome,Reads_Mapped_Confidently_to_Genome,Total_Genome_Coverage_rate,'
            'Total_CPGs_Detected,Genome_Coverage_rate_of_max_cell,CPGs_of_max_cell,UMIs_of_max_cell,Reads_of_max_cell,Saturation_of_max_cell,'
            'Genome_Coverage_rate_of_median_cell,CPGs_of_median_cell,UMIs_of_median_cell,Reads_of_median_cell,Saturation_of_median_cell,Fraction_Reads_in_Cells')

    summary_data = [
             samplename,
             cellnum,
             raw,
             vaild,
             rate_7f,
             rate_17lme,
             rate_7f17lme,
             conversion,
             cc_ratio,
             mapgenome,
             confidently,
             coverage,
             total_cpgs,
             cov_max_cell,
             cpgs_max_cell,
             umi_max_cell,
             reads_max_cell,
             saturation_max_cell,
             cov_median_cell,
             cpgs_median_cell,
             umi_median_cell,
             reads_median_cell,
             saturation_median_cell,
             fraction
           ]

    with open(os.path.join(outdir, f'{samplename}_wgs_summary.csv'), 'w') as fh:
        fh.write(header + '\n')
        fh.write(','.join(str(_).replace(',', '') for _ in summary_data)+ '\n')

if __name__ == "__main__":
    outcsv()
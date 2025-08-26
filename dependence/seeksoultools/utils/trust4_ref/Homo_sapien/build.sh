# https://github.com/liulab-dfci/TRUST4/tree/master?tab=readme-ov-file#build-custom-vjc-gene-database-files-for--f-and---ref
# https://www.imgt.org//download/V-QUEST/IMGT_V-QUEST_reference_directory/

BuildImgtAnnot.pl Homo_sapien > IMGT+C.fa
grep ">" IMGT+C.fa | cut -f2 -d'>' | cut -f1 -d'*' | sort | uniq > bcr_tcr_gene_name.txt
BuildDatabaseFa.pl refdata-gex-GRCh38-2020-A/fasta/genome.fa refdata-gex-GRCh38-2020-A/genes/genes.gtf bcr_tcr_gene_name.txt > bcrtcr.fa

# https://www.imgt.org/vquest/refseqh.html#refdir


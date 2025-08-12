import os

__srcdir = os.path.dirname(os.path.abspath(__file__))

ADAPTERS = {
    "MM": [["AAAAAAAAAAAA", "3"],  ],
    "DD": [["AAGCAGTGGTATCAACGCAGAGTACATGGG", "5"], ],
}

R1_MINLEN = 20
R2_MINLEN = 60

CHEMISTRY = {
    '__SO01V3':{
        'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8L8B8L10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'linker': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'Linker1.txt'),
                   os.path.join(__srcdir, 'barcode', 'SO01V3', 'Linker2.txt'),),
    },
    '__nolinker':{
        'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
    },
    "__P3CBGB":{
        'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
    },
    "DDV1":{
       'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'match_type': (1,),
        'adapter1': [["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"], ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ["TTTTTTTTTTTT", "5"]],
        'adapter2': [["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ["AAAAAAAAAAAA", "3"], ["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"]],
    },
    "DDV2":{
        'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        'match_type': (1,),
        'adapter1': [["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"], ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ["TTTTTTTTTTTT", "5"] ], 
        'adapter2': [["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ["AAAAAAAAAAAA", "3"], ["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"]],
    },
    "DDVS":{
        'shift': False,
        'structure': 'B17U12X24B10',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),
                    os.path.join(__srcdir, 'barcode', 'P3CBGB', 'sample.barcode.txt.gz'),),
        'match_type': (1,3,),
        'adapter1': [["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"], ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ["TTTTTTTTTTTT", "5"] ],
        'adapter2': [["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"], ["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ["AAAAAAAAAAAA", "3"]],
    },
    "DD5V1":{
       'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        'match_type': (1,),
        'sc5p': True,
        'adapter1': [["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"] ], 
        'adapter2': [["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ["TTTTTTTTTTTT", "5"], ["CCCATATAAGAAA", "3"], ["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"] ],
    },
    "MM":{
        'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'match_type': (1,),
        'adapter2': [["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ["AAAAAAAAAAAA", "3"], ["AAAAAAAAAAAA", "5"], ["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"]], 
        'adapter1': [["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"],["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"], ["TTTTTTTTTTTT", "5"]], 
    },
    "MM-D":{
        'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'match_type': (1,),
    },
    "DD-Q":{
        'shift': False,
        'structure': 'B17U12X17',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        # TSO:  AAGCAGTGGTATCAACGCAGAGTACATGG
	#ACACTCTTTCCCTACACGACGCTCTTCCGATCT
        'adapter1': [["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"], ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ], ## SP2 SP1
        'adapter2': [["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"], ["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ], ## TSO SP2 reverse complement
        'sc5p': None,
        'match_type': (1,),
    }
}

REF = {
    "human": {
        "ref": os.path.join(__srcdir, 'IMGT', 'human', 'vdj_IMGT_human', 'fasta', 'regions.fa'),
        "internal_data": os.path.join(__srcdir, 'IMGT', 'human', 'vdj_IMGT_human_igblast'),
        "auxiliary_data": os.path.join(__srcdir, 'igblast', 'optional_file', 'human_gl.aux'),
        "changeo_igblast_database": os.path.join(__srcdir, 'changeo', 'igblast', 'database'),
        "igblast_germline_path": os.path.join(__srcdir, 'changeo', 'germlines', 'imgt', 'human', 'vdj', 'imgt_human_*.fasta'),

        # "igblast_imgt_db_path": os.path.join(__srcdir, 'IMGT', 'human', 'vdj_IMGT_human_igblast')
    },
    "mouse": {
        "ref": os.path.join(__srcdir, 'IMGT', 'mouse', 'vdj_IMGT_mouse', 'fasta', 'regions.fa'),
        "internal_data": os.path.join(__srcdir, 'IMGT', 'mouse', 'vdj_IMGT_mouse_igblast'),
        "auxiliary_data": os.path.join(__srcdir, 'igblast', 'optional_file', 'mouse_gl.aux'),
        "changeo_igblast_database": os.path.join(__srcdir, 'changeo', 'igblast', 'database'),
        "igblast_germline_path": os.path.join(__srcdir, 'changeo', 'germlines', 'imgt', 'mouse', 'vdj', 'imgt_mouse_*.fasta'),

    },
    "monkey": {
        "ref": os.path.join(__srcdir, 'IMGT', 'monkey', 'vdj_IMGT_monkey', 'fasta', 'regions.fa'),
        "internal_data": os.path.join(__srcdir, 'IMGT', 'monkey', 'vdj_IMGT_monkey_igblast'),
        "auxiliary_data": os.path.join(__srcdir, 'igblast', 'optional_file', 'rhesus_monkey_gl.aux'),
        "changeo_igblast_database": os.path.join(__srcdir, 'changeo', 'igblast', 'database'),
        "igblast_germline_path": os.path.join(__srcdir, 'changeo', 'germlines', 'imgt', 'monkey', 'vdj', 'imgt_monkey_*.fasta'),

    },
    "rabbit": {
        "ref": os.path.join(__srcdir, 'IMGT', 'rabbit', 'vdj_IMGT_rabbit', 'fasta', 'regions.fa'),
        "internal_data": os.path.join(__srcdir, 'IMGT', 'rabbit', 'vdj_IMGT_rabbit_igblast'),
        "auxiliary_data": os.path.join(__srcdir, 'igblast', 'optional_file', 'rabbit_gl.aux'),
        "changeo_igblast_database": os.path.join(__srcdir, 'changeo', 'igblast', 'database'),
        "igblast_germline_path": os.path.join(__srcdir, 'changeo', 'germlines', 'imgt', 'rabbit', 'vdj', 'imgt_rabbit_*.fasta'),

    },
    "rat": {
        "ref": os.path.join(__srcdir, 'IMGT', 'rat', 'vdj_IMGT_rat', 'fasta', 'regions.fa'),
        "internal_data": os.path.join(__srcdir, 'IMGT', 'rat', 'vdj_IMGT_rat_igblast'),
        "auxiliary_data": os.path.join(__srcdir, 'igblast', 'optional_file', 'rat_gl.aux'),
        "changeo_igblast_database": os.path.join(__srcdir, 'changeo', 'igblast', 'database'),
        "igblast_germline_path": os.path.join(__srcdir, 'changeo', 'germlines', 'imgt', 'rat', 'vdj', 'imgt_rat_*.fasta'),
    },
}

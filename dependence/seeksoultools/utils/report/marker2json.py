import json
from collections import defaultdict


d = defaultdict(set)

with open("marker.tsv") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        else:
            # cell_type, genes_str = line.strip().split("\t")
            cell_type, genes_str = line.strip().split(":")
            for g in genes_str.split(","):
                d[g.replace(" ", "")].add(cell_type)

d = {k: ",".join(list(v)) for k, v in d.items()}

with open("marker.json", "w") as fh:
    json.dump(d, fp=fh, indent=2)


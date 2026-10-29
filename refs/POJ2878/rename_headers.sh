#!/bin/bash
# Rename NCBI GenBank headers (>CM190499.1 ...) to the names used in Sh_POJ2878_v1.gff3.gz (>Chr01g1),
# using column 5 (GenBank accession) -> column 1 (sequence name) of the NCBI assembly report.
# Fails if any header has no match.
set -euo pipefail
cd "$(dirname "$0")"

REPORT=GCA_060231685.1_POJ2878_1.0_assembly_report.txt
IN=GCA_060231685.1_POJ2878_1.0_genomic.fna.gz
OUT=POJ2878_1.0.genome.fa.gz

grep -v '^#' "$REPORT" | awk -F'\t' '{print $5"\t"$1}' > accession_to_name.tsv

zcat "$IN" | awk -F'\t' '
    NR == FNR { name[$1] = $2; next }
    /^>/ {
        acc = substr($1, 2); sub(/ .*/, "", acc)
        if (!(acc in name)) { print "no name for " acc > "/dev/stderr"; exit 1 }
        print ">" name[acc]; next
    }
    { print }
' accession_to_name.tsv - | gzip -c > "$OUT"

echo "wrote $OUT"

#!/bin/zsh

OUTDIR=/Users/mm/Documents/project/metagenome/3_MAG_gene_cluster/mag_SRR24442557
DEPTH_FILE="$OUTDIR/depth.txt"

if [ -f "$DEPTH_FILE" ]; then
    echo "Depth calculation completed successfully!"
    head -n 10 "$DEPTH_FILE"
else
    echo "Error: Depth calculation failed."
    exit 1
fi


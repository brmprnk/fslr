#!/bin/bash

# Directory containing BAM files
BAM_DIR="data/cristiano"

# Loop through all .bam files in the specified directory
for bamfile in "$BAM_DIR"/*.bam; do
    # Check if the corresponding .bam.bai file does not exist
    if [ ! -f "${bamfile}.bai" ]; then
        echo "Indexing $bamfile"
        # Use samtools to create the index
        samtools index "$bamfile"
    else
        echo "Index already exists for $bamfile"
    fi
done

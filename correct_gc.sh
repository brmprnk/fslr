#!/bin/bash

# Directory containing BAM files
input_directory="data/emc/first/bams/"

# Directory to save corrected BAM files
output_directory="data/emc/first/bams_corrected/"

# Create output directory if it doesn't exist
mkdir -p "$output_directory"

# Loop through all BAM files in the input directory
for bam_file in "$input_directory"/*.bam; do
    # Extract filename without extension
    filename=$(basename "$bam_file" .bam)

    echo "Processing $bam_file"

    # Compute GC bias
    computeGCBias -b "$bam_file" -o "data/emc/first/bams/${filename}_gc_bias.txt" --effectiveGenomeSize 2805636331 --genome data/reference_genomes/hg38.2bit  -bl beds/hg38-blacklist.v2.bed

    echo "Correcting $filename..."

    # Correct GC bias
    correctGCBias -b "$bam_file" -o "${output_directory}/${filename}_corrected.bam" --gcBiasCorrectedFile "data/emc/first/bams/${filename}_gc_bias.txt"
done

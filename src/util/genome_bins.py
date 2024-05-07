
def bins_to_bed(bin_size, chrom_sizes_path, output_path, exclude_chroms=[]):
    # Read the file
    chrom_sizes = {}
    with open(chrom_sizes_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if parts[0] in ["chrX", "chrY"] or parts[0][3:].isdigit():
                chrom_sizes[parts[0]] = int(parts[1])


    # Sort chromosomes as per the requirement: 1-22, X, Y
    chrom_names = list(chrom_sizes.keys())
    chrom_names.sort(key=lambda x: int(x[3:]) if x[3:].isdigit() else (23 if x == "chrX" else 24))
    sorted_chrom_names = chrom_names

    # Remove excluded chromosomes
    for chrom in exclude_chroms:
        sorted_chrom_names.remove(chrom)

    # Now, create the BED file again with the full chromosome sizes and in the specified order
    with open(output_path, "w") as bed_file:
        for chr_name in sorted_chrom_names:
            chr_size = chrom_sizes[chr_name]
            for start in range(0, chr_size, bin_size):
                end = start + bin_size
                if end > chr_size:
                    end = chr_size
                bed_file.write(f"{chr_name}\t{start}\t{end}\n")


if __name__ == '__main__':
    
    # Let's read the uploaded chromosome sizes file to understand its format
    chrom_sizes_path = "data/reference_genomes/hg38.chrom.sizes"
    bin_size = 100000
    output_path = "beds/hg38_100kb.bed"

    bins_to_bed(bin_size, chrom_sizes_path, output_path)


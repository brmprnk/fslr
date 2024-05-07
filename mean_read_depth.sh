directory='./data/cristiano_stein/crc'

# Check if the provided path is a directory
if [ ! -d "$directory" ]; then
    echo "Error: $directory is not a directory."
    exit 1
fi

# Use a loop to iterate over files in the directory and print their names
for file in "$directory"/*; do

    if [ -f "$file" ]; then

        # Check if the file is a bam file
        if [[ $file == *.bam ]]; then

            # Check if the file is not a .bai file
            if [[ $file == *.bai ]]; then
                continue
            fi

            echo "Running big on $file"
            bedtools coverage -a beds/chrom_arms.tsv -b "$file" -mean > "${file%$suffix}_chrom_arms.bedgraph"
       fi

    fi
done
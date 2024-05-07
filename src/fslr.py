"""
"""
import sys
import os
import argparse
import yaml
from datetime import datetime
from time import time
import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
from util import logger
from util.Dataset import Dataset

# This is the Project Root Directory, assuming fslr.py is in the src folder
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PARSER = argparse.ArgumentParser(prog='fslr.py', description="Calculate Delfi log2(short/long) ratios for a given dataset.")
PARSER.add_argument('--config', '-c',
                    dest='config_file',
                    metavar='FILE',
                    help="path to the config file",
                    default='configs/hk.yaml')
PARSER.add_argument('--experiment', '-e',
                    help="Name of experiment",
                    default="experiment")
PARSER.add_argument('--results-path', '-r',
                    dest='results_path',
                    help="where to store the results",
                    default="results")
                    

def main():
    # Get the arguments from the command line
    args = PARSER.parse_args()

    # Get the config file
    try:
        with open(args.config_file, 'r') as file:
            config = yaml.load(file, Loader=yaml.FullLoader)
    except FileNotFoundError as e:
        print(f'Could not find config file: {e}', flush=True)
        sys.exit(1)

    # Create the output directory
    dt_string = datetime.now().strftime("%m-%d_%H-%M")
    experiment_name = config['GLOBAL_PARAMS']['experiment']
    if args.experiment != 'experiment':
        experiment_name = args.experiment
    output_dir = os.path.join(ROOT_DIR, args.results_path, dt_string + '_' + experiment_name)
    args.output_dir = output_dir
    os.makedirs(output_dir, exist_ok=True)

    # Set the output file
    logger.OUTPUT_FILE = output_dir
    logger.success(f"Running FSLR tool for {experiment_name} and saving results to {output_dir}")

    # Get the start time
    start_time = time()

    # Merge command line arguments with config file, command line arguments have priority
    args = vars(args)
    args = {**args, **config['GLOBAL_PARAMS']}
    args['output_dir'] = output_dir
    logger.info(str(args))

    data = Dataset(args)
    data.files = data.files[:2]

    calculate_fslr_parallel(args, data)

    # Print the total time
    logger.success(f'FSLR tool completed. Total time: {time() - start_time} seconds')

def calculate_fslr_parallel(args: dict, data: Dataset) -> None:
    """

    Returns
    -------
    None
    """
    nr_of_processes = min(mp.cpu_count(), len(data.files), 12)
    logger.info(f'Calculating fslr for {len(data.files)} bam files using {nr_of_processes} processes. Bed file used: {data.bed_files[0].name}')

    # Use a Manager object to share the result array between processes
    with mp.Manager() as manager:
        flsr_manager = manager.list()

        with mp.Pool(processes=nr_of_processes) as pool:
            pool.starmap(calculate_fslr_worker, [(flsr_manager, data.labels[i], data.files[i], data.bed_files[0], args) for i in range(len(data.files))])

        # Sort the manager by bam file index (alphabetically)
        flsr_manager = sorted(flsr_manager, key=lambda x: x[0])

    all_fslr = np.array([np.array(c[0][1]) for c in flsr_manager])
    all_short_longs = np.array([np.array(c[1][1]) for c in flsr_manager])
    logger.info("We've calculated ratios now just saving")
        
    file_labels = np.array([np.array(c[0][0]) for c in flsr_manager])

    np.save(os.path.join(args['output_dir'], f"labels.npy"), file_labels)
    np.save(os.path.join(args['output_dir'], f"{data.bed_files[0].name}_fslr.npy"), all_fslr.astype(np.float32))
    np.save(os.path.join(args['output_dir'], f"{data.bed_files[0].name}_short_longs.npy"), all_short_longs.astype(np.float32))

    return all_fslr

def calculate_fslr_worker(flsr_manager: mp.Manager, bam_file_index: int, bam_file_path: str, bed: pd.DataFrame, args: dict) -> None:
    """

    """
    logger.info(f"BAM {bam_file_index}: {bam_file_path} -- Calculating flsr for {bed.shape[0]} bins.")

    # Delfi definitions
    short_frag = (100, 150)
    long_frag = (151, 220)

    # Open the BAM file
    try:
        bam = pysam.AlignmentFile(bam_file_path, 'rb')
    except Exception as e:
        print(f"Something went wrong while opening the BAM file: {bam_file_path}")
        print(e)
        return
    
    # Store ratios of short to long fragments
    fslr = np.zeros((bed.shape[0]), dtype=np.float32)
    total_short = np.zeros((bed.shape[0]), dtype=np.uint32)
    total_long = np.zeros((bed.shape[0]), dtype=np.uint32)

    # Store fragment lengths
    fragment_lengths = np.zeros(args['max_frag_size'] + 1, dtype=np.int64)
    
    # From GCparagon:
    exclude_flags = np.uint32(3852)  # = 256 + 2048 + 512 + 1024 + 4 + 8
    exclude_flags_binary = bin(exclude_flags)

    # Save the bin size for easier access at blacklisted regions removal
    bin_size = 0

    for genome_bin in bed.itertuples():
        index, chrom, start, end = genome_bin

        if bin_size == 0:  # Now the code in the blacklist removals knows the bin size
            bin_size = end - start

        # Make sure the end is not larger than the chromosome length
        chrom_length = bam.get_reference_length(chrom)
        if end > chrom_length:
            print(f"End of bin {end} is larger than the chromosome length {chrom_length}, setting end to chrom_length - 1.")
            end = chrom_length - 1

        # Get all the reads in the bin.
        filtered_alignments = filter(lambda a:
                                            # bin(np.uint32(a.flag) & paired_flag) == paired_flag_binary and
                                            a.is_paired and
                                            bin(~np.uint32(a.flag) & exclude_flags) == exclude_flags_binary and
                                            a.is_forward != a.mate_is_forward and
                                            (args['min_frag_size'] <= a.template_length <= args['max_frag_size']) and not a.is_unmapped,
                                            bam.fetch(contig=chrom, start=start, stop=end))
        
        short = 0
        long = 0
        for read in filtered_alignments:
            if read.template_length <= short_frag[1] and read.template_length >= short_frag[0]:
                short += 1
            elif read.template_length <= long_frag[1] and read.template_length >= long_frag[0]:
                long += 1

            fragment_lengths[read.template_length] += 1

        fslr[index] = np.log2(max(short, 1) / max(long, 1))  # Avoid division by zero
        total_short[index] = short
        total_long[index] = long


    np.save(os.path.join(args['output_dir'], f"{bam_file_index}_fslr_nobl.npy"), fslr)

    print(f"{bam_file_index}: {bam_file_path} -- Now removing blacklisted reads.")
    # Remove reads in blacklisted regions
    if bin_size != 0:
        blacklist = pd.read_csv('beds/hg38-blacklist.v2.bed', sep='\t', header=None, names=['chr', 'start', 'end', 'reason'])

        for blacklisted_region in blacklist.itertuples():
            _, chrom, region_start, region_end, _ = blacklisted_region

            # print(f"Blacklisted region: {chrom}:{region_start}-{region_end}")
            
            # We need to calculate at what index the blacklisted region is in the bed file
            # Round bottom down and top up
            genome_bin_region_start = int(np.floor(region_start / bin_size)) * bin_size
            genome_bin_region_end = int(np.ceil(region_end / bin_size)) * bin_size
            # print("Genome bins: ", genome_bin_region_start, genome_bin_region_end)
            # print("Now finding the index in the bed file.")

            region_start_bed_row = bed[(bed['chrom'] == chrom) & (bed['start'] == genome_bin_region_start)]
            if region_start_bed_row.empty:
                # print("Start row empty, this should never happen")
                continue
            index_start = region_start_bed_row.index[0]
            region_end_bed_row = bed[(bed['chrom'] == chrom) & (bed['end'] == genome_bin_region_end)]
            if region_end_bed_row.empty:
                # print("End row empty, this could happen if the have reached the end of a chromosome")
                chr_end_bin = bed[bed['chrom'] == chrom].iloc[-1]
                index_end = int(chr_end_bin.name) if not chr_end_bin.empty else index_start
            else:
                index_end = region_end_bed_row.index[0]

            # print("Start index: ", index_start)
            # print("End index: ", index_end)
            start, end = region_start, region_end

            for bin_loop_idx, bed_idx in enumerate(range(index_start, index_end + 1)):
                # print(f"Loop : {bin_loop_idx}, Blacklisted region: {chrom}:{region_start}-{region_end} is in bin {bed_idx}, because {genome_bin_region_start}-{genome_bin_region_end} == bed.iloc[{bed_idx}]")

                if index_start != index_end:  # Logic for when the blacklisted region spans multiple bins
                    if bin_loop_idx == 0:  # First bin is from start to end of first bin
                        # print("First bin, so only remove from start")
                        start = region_start
                        end = start + bin_size
                    else:
                        if bin_loop_idx == index_end - index_start:  # Last bin is from start of last bin to end
                            # print("Last bin, so only remove from end")
                            start = int(np.floor(end / bin_size)) * bin_size
                            end = region_end
                        else:  # Middle bins are from start to end of bin
                            # print("One of the middle bins, so update start to be inside middle bin.")
                            start = region_start + bin_loop_idx * bin_size
                            end = start + bin_size
                            start = int(np.floor(start / bin_size)) * bin_size
                            end = int(np.floor(end / bin_size)) * bin_size
                # print("Start and end: ", start, end)
                # Get all the reads in the bin.
                filtered_alignments = filter(lambda a:
                                                # bin(np.uint32(a.flag) & paired_flag) == paired_flag_binary and
                                                a.is_paired and
                                                bin(~np.uint32(a.flag) & exclude_flags) == exclude_flags_binary and
                                                a.is_forward != a.mate_is_forward and
                                                (args['min_frag_size'] <= a.template_length <= args['max_frag_size']) and not a.is_unmapped,
                                                bam.fetch(contig=chrom, start=start, stop=end))
                
                short = 0
                long = 0
                read_ctr = 0
                for read in filtered_alignments:
                    read_ctr += 1
                    if read.template_length <= short_frag[1] and read.template_length >= short_frag[0]:
                        short += 1
                    elif read.template_length <= long_frag[1] and read.template_length >= long_frag[0]:
                        long += 1

                    # Remove the read from the fragment lengths
                    fragment_lengths[read.template_length] -= 1


                # print(f"In the region {chrom}:{start}-{end}, removing s: {short}, l: {long}, reads: {read_ctr}")
                total_short[index] -= short
                total_long[index] -= long
                fslr[index] = np.log2(max(total_short[index], 1) / max(total_long[index], 1))  # Avoid division by zero

    np.save(os.path.join(args['output_dir'], f"{bam_file_index}_fslr.npy"), fslr)
    np.save(os.path.join(args['output_dir'], f"{bam_file_index}_lengths.npy"), fragment_lengths)

    # Close the BAM file
    bam.close()

    # Zip total_short and total_long such that they form a tuple
    short_longs = list(zip(total_short, total_long))

    flsr_manager.append([(bam_file_index, fslr), (bam_file_index, short_longs)])

if __name__ == '__main__':
    main()
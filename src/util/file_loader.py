"""
Module for loading files, which for now includes:
- BED files
"""
import pickle
import pandas as pd
from typing import List

def read_bed_file(bed_path: str, header=False) -> pd.DataFrame:
    """
    Read a bed file and return a pandas dataframe with the columns:
    chrom, start, end, meta

    Parameters
    ----------
    bed_path : str
        Path to the bed file
    header : bool
        Whether or not the bed file has a header

    Returns
    -------
    pd.DataFrame
        A pandas dataframe with the columns chrom, start, end, meta
    """
    with open(bed_path, 'rt', encoding='utf-8') as file:
        if header:
            file.readline()  # skip header

        # Get each bed line region and put it into a pandas dataframe
        return_data = []
        for bed_line in file.readlines():
            if bed_line not in ('', '\n', None):
                chrom, region_start, *meta_info = bed_line.strip().split()
                return_data.append((str(chrom), int(region_start), meta_info))

        return pd.DataFrame(return_data, columns=['chrom', 'start', 'meta'])

def filter_bam_files(bam_files, args) -> List:
    """
    Filter the bam files based on cancer type.
    """
    with open(args['labels'], 'rb') as file:
        labels = pickle.load(file)

    filtered_bam_files = []

    for bam_file in bam_files:
        ctype = labels[bam_file]
        if args['ov'] and ctype == 'ovarian cancer':
            filtered_bam_files.append(bam_file)
        elif args['pr'] and ctype == 'prostate cancer':
            filtered_bam_files.append(bam_file)
        elif args['he'] and ctype == 'healthy':
            filtered_bam_files.append(bam_file)

    return filtered_bam_files

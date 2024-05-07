"""
Create a class that abstracts the dataset such that dataset specific settings and functions
are hidden from the main script.
"""
import os
import glob
import pandas as pd
import numpy as np
from util.file_loader import read_bed_file
from typing import List

class Dataset:

    def __init__(self, args: dict):
        """
        Initialize the dataset class.

        Parameters
        ----------
        args : dict
            Dictionary containing the arguments for the coverage tool, either from the config file or from the command line.
        """
        self.args = args
        self.files = self.get_bam_files()

        print(f'Number of files: {len(self.files)}')
        self.labels = self.get_labels()

        self.bed_files = self.get_bed_files()

    def get_bigwig_files(self) -> List:

        bigwig_files = []
        dir = self.args['dir']

        if 'cristiano' in dir:
            if 'hea' not in dir and 'crc' not in dir:
                bigwig_files = glob.glob(os.path.join(dir + '/hea', '*.bw'))
                bigwig_files += glob.glob(os.path.join(dir + '/crc', '*.bw'))

                # For now only get 100 binsize
                # bigwig_files = [bigwig_file for bigwig_file in bigwig_files if '{}.bw'.format(self.args['bin_size']) in bigwig_file]

                bigwig_files = sorted(bigwig_files)

        # if os.path.isdir(dir):
        else:
            bigwig_files = glob.glob(os.path.join(dir, '*.bw'))
            bigwig_files = sorted(bigwig_files)

        return bigwig_files

    def get_bam_files(self) -> List:
        """
        Get the bam files from the config file.

        Returns
        -------
        List
            List of bam files.
        """
        bam_files = []
        bam_dir = self.args['dir']

        if os.path.isdir(bam_dir):  # Get all the bam files in the directory
            bam_files = glob.glob(os.path.join(bam_dir + '/**/', '*.bam'), recursive=True)
            bam_files = sorted(bam_files)

            # Reverse the order
            # bam_files = bam_files[::-1]

        else:  # Single bam file or a list of bam files
            bam_files = [bam_dir]

        # Check if the bam file is not yet present in results/04-23_11-57_cris_hea_30k
        # known_files = glob.glob('results/04-23_23-13_cris_CRC_30k/*coverage.npy')
        # known_files = [file.split('/')[-1].split('_')[0] for file in known_files]

        # print(len(bam_files), 'before filtering')
        # bam_files = [bam_file for bam_file in bam_files if bam_file.split('/')[-1].split('.')[0] not in known_files]
        # print(len(bam_files), 'after filtering')


        return bam_files
    
    def get_labels(self) -> List:
        """
        Get the labels for the bam files.

        Returns
        -------
        List
            List of labels.
        """
        labels = []

        if 'labels' not in self.args or self.args['labels'] == '':
            for i in range(len(self.files)):
                labels.append(self.files[i].split('/')[-1].split('.')[0])

        else:
            labels = pd.read_csv(self.args['labels'], index_col=0)

            # @TODO: If the labels are not in the same order as the bam files, then sort them
    
        return labels
    
    def get_bed_files(self) -> List:
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
        # Check if the bed file is a directory
        if os.path.isdir(self.args['bed']):
            bed_files = glob.glob(os.path.join(self.args['bed'], '*.bed'))
            bed_files = sorted(bed_files)
        else:
            bed_files = [self.args['bed']]
        
        bed_dfs = []
        for bed_file in bed_files:
            bed = pd.read_csv(bed_file, sep='\t', header=None)
            # Keep only the first 3 columns and rename them
            bed = bed.iloc[:, :3]
            bed.columns = ['chrom', 'start', 'end']

            bed.name = bed_file.split('/')[-1].split('.')[0]
            bed_dfs.append(bed)

        return bed_dfs

        # with open(self.args['bed'], 'rt', encoding='utf-8') as file:
        #     # Get each bed line region and put it into a pandas dataframe
        #     return_data = []
        #     for bed_line in file.readlines():
        #         if bed_line not in ('', '\n', None):
        #             chrom, region_start, *meta_info = bed_line.strip().split()
        #             return_data.append((str(chrom), int(region_start), meta_info))

        #     return pd.DataFrame(return_data, columns=['chrom', 'start', 'meta'])
    
    def __str__(self) -> str:
        """
        Return the string representation of the dataset.

        Returns
        -------
        str
            String representation of the dataset.
        """
        return f"Dataset object with {len(self.files)} bam files from {self.args['dir']}. \n {self.args}"

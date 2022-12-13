#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 2022
@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import numpy as np

def main(args):
    
    bed_fname = args.i
    df_bed = pd.read_csv(bed_fname, sep='\t', header=None, names = ['chrom', 'start', 'end', 'cn', 'cell_id'])
    
    df_bed = df_bed.sort_values('chrom')

    cell_list = sorted(list(set(df_bed['cell_id'])))

    breakpoint_data = []
    for chrom in df_bed['chrom'].unique():
        df_bed_chrom = df_bed[df_bed['chrom'] == chrom]

        breakpoint_list = sorted(list(set(list(df_bed_chrom['start']) + list(df_bed_chrom['end']))))

        for cell in cell_list:
            df_bed_cell = df_bed_chrom[df_bed_chrom['cell_id'] == cell]

            for idx, breakpoint in enumerate(breakpoint_list):
                df_selected = df_bed_cell[(df_bed_cell['start'] < breakpoint) & (breakpoint <= df_bed_cell['end'])]
                # print(idx, breakpoint, sep='\t')
                if len(df_selected) == 0:
                    breakpoint_data.append([cell, chrom, idx, idx, 2])
                elif len(df_selected) == 1:
                    breakpoint_data.append([cell, chrom, idx, idx, df_selected.iloc[0]['cn']])
                else:
                    raise Exception(f"multiple entries for chrom {chrom}, cell {cell}, breakpoint idx {idx}, breakpoint {breakpoint} in input")
    
    df_breakpoint = pd.DataFrame(breakpoint_data, columns = ['node', 'chrom', 'start', 'end', 'cn_a'])
    
    df_breakpoint.to_csv(f"{args.o}_cn_profiles.csv", index=False)

    df_breakpoint = pd.DataFrame(breakpoint_data, columns = ['sample_id', 'chrom', 'start', 'end', 'cn_a'])
    df_breakpoint['sample_id'] = 'sample_' + df_breakpoint['sample_id'].astype(str)
    df_breakpoint.to_csv(f"{args.o}_cn_profiles_medicc2.tsv", index=False, sep="\t")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input bed file name', required=True)
    parser.add_argument('-o', type=str, help='output prefix', required=True)
    
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)

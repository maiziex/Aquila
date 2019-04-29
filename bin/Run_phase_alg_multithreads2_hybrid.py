#!/usr/bin/env python
import pdb
#pdb.set_trace()
from collections import defaultdict
from multiprocessing import Pool,cpu_count,active_children,Manager
import time
import pickle
import os
import numpy as np
from Assign_phase_block_v3 import * 
import sys
from scipy.special import comb
import glob
from Molecule_phase_alg2_withProbModel_v3_MT2_hybrid import *
from argparse import ArgumentParser
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 

parser = ArgumentParser(description="phase molecules:")
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=22)
parser.add_argument('--overlap_threshold','-t1',type=int,help="overlap threshold", default=3)
parser.add_argument('--support_threshold','-t2',type=int,help="support threshold", default=5)
parser.add_argument('--sample_name','-s',help="sample name")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs")
parser.add_argument('--h5_dir','-h_dir', help="Directory to store h5 file")

args = parser.parse_args()


def Run_phasing_all(chr_start,chr_end,overlap_threshold,support_threshold,output_dir,h5_dir,sample_name):
    pool = Pool(chr_end - chr_start + 1)
    for chr_num in range(chr_start,chr_end+1):
        pool.apply_async(Phase_start,(output_dir,h5_dir,sample_name,chr_num,chr_num,overlap_threshold,support_threshold,"xin"))
        #Phase_start(output_dir,h5_dir,sample_name,chr_num,chr_num,overlap_threshold,support_threshold,"xin")
    pool.close()
    while len(active_children()) > 1:
        time.sleep(0.5)
    pool.join()
    print("All Done~")




if __name__ == "__main__":
    output_dir = args.out_dir
    h5_dir = args.h5_dir
    if os.path.exists(output_dir):
        print("using existing output folder: " + output_dir)
    else:
        os.makedirs(output_dir)
    overlap_threshold = args.overlap_threshold
    support_threshold = args.support_threshold
    chr_start = args.chr_start
    chr_end = args.chr_end
    sample_name = args.sample_name
    Run_phasing_all(chr_start,chr_end,overlap_threshold,support_threshold,output_dir,h5_dir,sample_name)

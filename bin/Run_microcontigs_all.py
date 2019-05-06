#!/usr/bin/env python
import pdb
#pdb.set_trace()
import os
import numpy as np
from argparse import ArgumentParser
from Concatenate_contigs_all_v4_extend_for_HCbk import *
from multiprocessing import Pool,cpu_count,active_children,Manager
import time
script_path = os.path.dirname(os.path.abspath( __file__ ))

code_path = script_path + "/" 

parser = ArgumentParser(description="Run microcontigs all:")
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--in_dir','-i_dir', help="Directory for inputs")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs", default='results_/')
parser.add_argument('--output_prefix','-output_prefix', help="output file")

args = parser.parse_args()


def Run_all(chr_start,chr_end,in_dir,out_dir,output_file):
    pool = Pool(chr_end - chr_start + 1)
    for chr_num in range(chr_start,chr_end+1):
        output_file = output_prefix + "_chr" + str(chr_num) + ".fasta"
        in_dir_2 = in_dir + "chr" + str(chr_num) + "_files_cutPBHC/"
        pool.apply_async(Concatenate_start,(in_dir_2,out_dir,output_file,"xin"))
    pool.close()
    while len(active_children()) > 1:
        time.sleep(0.5)
    pool.join()

    print("all done~")

   


if __name__ == "__main__":
    out_dir = args.out_dir
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    chr_start = args.chr_start
    chr_end = args.chr_end
    in_dir = args.in_dir
    out_dir = args.out_dir
    output_prefix= args.output_prefix
    Run_all(chr_start,chr_end,in_dir,out_dir,output_prefix)
 

#!/usr/bin/env python
import pdb
#pdb.set_trace()
import os
import sys
from argparse import ArgumentParser
from multiprocessing import Pool,cpu_count,active_children,Manager
from subprocess import Popen
import time
from Make_supercontig_based_on_HCbk_from_sam_v5 import *
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 

parser = ArgumentParser(description="Run all:")
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=22)
parser.add_argument('--ref_file','-ref',help="ref file")
parser.add_argument('--ref_dir','-dir',help="ref dir")
parser.add_argument('--cut_threshold','-ct',type=int,help="cut threshold: 100000")
parser.add_argument('--out_dir','-o_dir', help="Directory for outputs")
parser.add_argument('--phase_cut_folder','-pc_dir', help="Directory to store phase info")
args = parser.parse_args()
print(code_path)


def align_cmd(ref_file,microcontig_fasta,microcontig_sam,xin):
    use_cmd = "minimap2 -a " + ref_file + " " + microcontig_fasta + " > " + microcontig_sam 
    Popen(use_cmd,shell=True).wait() 


def align_microcontigs(chr_start,chr_end,ref_file,out_dir):
    pool = Pool(chr_end - chr_start + 1)
    for chr_num in range(chr_start,chr_end + 1):
        microcontig_fasta = out_dir + "XinMagic_cutPBHC_minicontig_chr" + str(chr_num) + ".fasta"
        microcontig_sam = out_dir + "XinMagic_cutPBHC_minicontig_chr" + str(chr_num) + ".sam"
        pool.apply_async(align_cmd,(ref_file,microcontig_fasta,microcontig_sam,"xin"))
    pool.close()
    while len(active_children()) > 1:
        time.sleep(0.5)
    pool.join()
    print("alignment all done~")


def Concatenate_microcontigs_all(chr_start,chr_end,cut_threshold,out_dir,phase_cut_folder,ref_dir):
    pool = Pool(chr_end - chr_start + 1)
    for chr_num in range(chr_start,chr_end + 1):
        pool.apply_async(Contig_start,(chr_num,cut_threshold,ref_dir,out_dir,phase_cut_folder,"xin"))
    pool.close()
    while len(active_children()) > 1:
        time.sleep(0.5)
    pool.join()


def Concatenate_wgs(chr_start,chr_end,out_dir):
    for chr_num in range(chr_start,chr_end + 1):
        one_contig = out_dir + "XinMagic_Contig_chr" + str(chr_num) + ".fasta"
        cat_cmd = "cat " + one_contig  + " >> " + out_dir + "XinMagic_contig.fasta"
        Popen(cat_cmd,shell=True).wait()
    print("all done~")






if __name__ == "__main__":
    chr_start = args.chr_start
    chr_end = args.chr_end
    ref_file = args.ref_file
    ref_dir = args.ref_dir
    out_dir = args.out_dir
    cut_threshold = args.cut_threshold
    phase_cut_folder = args.phase_cut_folder
    align_microcontigs(chr_start,chr_end,ref_file,out_dir)
    Concatenate_microcontigs_all(chr_start,chr_end,cut_threshold,out_dir,phase_cut_folder,ref_dir)
    Concatenate_wgs(chr_start,chr_end,out_dir)

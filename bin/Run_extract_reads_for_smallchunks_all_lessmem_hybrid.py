#!/usr/bin/env python
import pdb
#pdb.set_trace()
import glob
import os
import numpy as np
from argparse import ArgumentParser
from collections import defaultdict
from multiprocessing import Pool,cpu_count,active_children
from Extract_qname_from_phased_molecule_cut_phase_blocks_v3_hybrid import *
import time
from subprocess import Popen

parser = ArgumentParser(description="Extract reads for all small chunks:")
parser.add_argument('--phase_cut_folder','-cut_folder',help="folder to save cutting profile")
parser.add_argument('--cut_threshold','-c',type=int,help="cutting threshold", default=100000)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--num_threads','-n',type=int,help="number of threads", default=8)
parser.add_argument('--h5_folder','-h5_folder',help="folder to save h5")
parser.add_argument('--sample_name','-s',help="sample name")
parser.add_argument('--mole_len_dict_prev_file','-file1',help="file1")
parser.add_argument('--mole_len_dict_file','-file2',help="file2")
parser.add_argument('--fastq_folder','-f',help="raw fastqs folder")
parser.add_argument('--out_dir','-o_dir',help="out_dir folder")
parser.add_argument("-mole_len_dict_prev_file", "-dict_prev")
parser.add_argument("-mole_len_dict_file", "-dict")

args = parser.parse_args()
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 

def check_mole_len(cur_mole_file):
    cat_cmd = "cat " + cur_mole_file + " | wc -l "
    mole_len_curr = Popen(cat_cmd,shell=True).wait()
    return mole_len_curr


def Run_extract_reads_for_smallchunks_all(out_dir,phase_cut_folder,cut_threshold,h5_folder,sample_name,fastq_folder,chr_start,chr_end,num_threads,mole_len_dict_prev,mole_len_dict):
    count = 1
    round_times = int(np.ceil((chr_end - chr_start + 1)/num_threads))
    old_list = defaultdict(list)
    for ii in range(0,num_threads):
        for jj in range(chr_start + round_times*ii,chr_start + round_times*ii + round_times):
            if jj <= chr_end:
                old_list[ii+1].append(jj)

    new_list = defaultdict(list)
    times = 1
    for key,val in old_list.items():
        if times%2 == 1:
            new_list[key] = val
        elif times%2 == 0:
            new_list[key] = val[::-1] 
        times += 1
    
    run_list = defaultdict(list)
    for key,val in new_list.items():
        for ii in range(round_times):
            if ii < len(val):
                run_list[ii+1].append(val[ii])

    for key, val in run_list.items():
        pool = Pool(processes=num_threads)
        mole_num_add = 0
        for chr_num in val:
            phased_h5_file = phase_cut_folder + "chr" + str(chr_num) + ".phased_final_cut_by_" +str(cut_threshold)
            PS_flag_dict_cut_file = phase_cut_folder + "chr" + str(chr_num) + ".phased_final_cut_by_" + str(cut_threshold) + "_phase_blocks.p"
            output_dir = out_dir + "chr" + str(chr_num) + "_files_cutPBHC/" 
            mole_qname_dict_file = h5_folder + sample_name + "_chr" + str(chr_num) +  "_qname.p"
            qname_pos_dict_file = h5_folder + sample_name + "_chr" + str(chr_num) + "_qname_pos.p"
            chr_fastq = fastq_folder + "fastq_by_Chr_" + str(chr_num)
            mole_num_start = 1
            mole_num_end = mole_len_dict[chr_num]
            mole_num_add = mole_len_dict_prev[chr_num]
            print("mole_num_start,mole_num_end,mole_num_add:")
            print(mole_num_start,mole_num_end,mole_num_add)
            if os.path.exists(output_dir):
                print("output folder: " + output_dir + ", existing..., delete old folder first")
                rm_cmd = "rm -rf " + output_dir
                Popen(rm_cmd,shell=True).wait()
                os.makedirs(output_dir)
            else:
                os.makedirs(output_dir)
            pool.apply_async(Extract_start,(output_dir,chr_num,phased_h5_file,PS_flag_dict_cut_file,mole_qname_dict_file,qname_pos_dict_file,chr_fastq,mole_num_start,mole_num_end,mole_num_add,"xin"))
            #Extract_start(output_dir,chr_num,phased_h5_file,PS_flag_dict_cut_file,mole_qname_dict_file,qname_pos_dict_file,chr_fastq,mole_num_start,mole_num_end,mole_num_add,"xin")
            mole_num_add += mole_num_end
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()
    print("All Done~")




if __name__ == "__main__":
    out_dir = args.out_dir
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    phase_cut_folder = args.phase_cut_folder
    cut_threshold = args.cut_threshold
    h5_folder= args.h5_folder
    fastq_folder = args.fastq_folder
    chr_start = args.chr_start
    chr_end = args.chr_end
    num_threads = args.num_threads
    sample_name = args.sample_name
    mole_len_dict_prev_file = args.mole_len_dict_prev_file
    mole_len_dict_prev = pickle.load(open(mole_len_dict_prev_file,"rb"))
    mole_len_dict_file = args.mole_len_dict_file
    mole_len_dict = pickle.load(open(mole_len_dict_file,"rb"))
    Run_extract_reads_for_smallchunks_all(out_dir,phase_cut_folder,cut_threshold,h5_folder,sample_name,fastq_folder,chr_start,chr_end,num_threads,mole_len_dict_prev,mole_len_dict)

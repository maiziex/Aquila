#!/usr/bin/env python
import pdb
#pdb.set_trace()
from subprocess import Popen,PIPE,STDOUT
from argparse import ArgumentParser
import os
import sys
from multiprocessing import Pool,cpu_count,active_children,Manager
import time
from collections import defaultdict
import pickle
import numpy as np
from Merge_reads_for_PB_of_merged_libs_v2 import *


script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
__author__ = "Xin Zhou@Stanford"
parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--bam_file_list','-bam',help="Required parameter, BAM file list, each BAM file is seperately by comma \",\". For example: 1.bam,2.bam" ,required=True)
parser.add_argument('--vcf_file_list','-v',help="Required parameter, VCF file list, each VCF file is seperately by comma \",\" . For example: 1.vcf,2.vcf",required=True)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from, default = 1", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by, default = 23", default=23)
parser.add_argument('--sample_name_list','-name',help="Required parameter, sample name list, each sample name is seperately by comma \",\". For example: S24385_lysis_2,S24385_lysis_2H",required=True)
parser.add_argument('--out_dir','-o', help="Directory to store Aquila assembly results, default = ./Assembly_results",default="./Asssembly_results")
parser.add_argument('--uniq_map_dir','-uniq_dir', help="Required parameter, Directory for 100-mer uniqness, you can run ./install.sh to download \"Uniqness_map\" for GRCh38",required=True)
parser.add_argument('--num_threads','-t_chr',type=int,help="number of threads, default = 8 (recommended)", default=8)
parser.add_argument('--num_threads_for_samtools_sort','-t',type=int,help="number of threads for samtools sort, default = 20", default=20)
parser.add_argument('--block_threshold','-bt',type=int,help="phase block threshold, default = 200000(200kb)",default=200000)
parser.add_argument('--block_len_use','-bl',type=int,help="phase block len threshold, default = 100000(100kb)",default=100000)

args = parser.parse_args()
bam_list = [item for item in args.bam_file_list.split(',')]
vcf_list = [item for item in args.vcf_file_list.split(',')]
sample_list = [item for item in args.sample_name_list.split(',')]

def Get_fragment_files(bam_file,vcf_file,chr_start,chr_end,h5_dir,sample_name,num_threads):
    use_cmd = "python3 " + code_path + "Run_h5_all_multithreads.py" + " --bam_file " + bam_file + " --vcf_file " + vcf_file + " --sample_name " + sample_name + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --mbq 13 --mmq 20 --boundary 50000 " + " --num_threads " + str(num_threads) + " --out_dir " + h5_dir
    Popen(use_cmd,shell=True).wait()


def Get_highconf_profile(bam_file,chr_start,chr_end,HighConf_file_dir,uniq_map_dir):
    use_cmd = "python3 " + code_path + "Generate_highconf_cut_profile_v2.py" + " --bam_file " +  bam_file + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --mbq 13 --mmq 20 " + " --out_dir " + HighConf_file_dir + " --uniq_map_dir " + uniq_map_dir
    Popen(use_cmd,shell=True).wait()
  

def Haplotying_fragments(chr_start,chr_end,phased_file_dir,h5_dir,sample_name):
    use_cmd = "python3 " + code_path + "Run_phase_alg_multithreads2_hybrid.py" + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --overlap_threshold 3 --support_threshold 5 " + " --out_dir " + phased_file_dir  + " --h5_dir " + h5_dir + " --sample_name " + sample_name
    Popen(use_cmd,shell=True).wait()


def Cut_phase_blocks(chr_start,chr_end,block_threshold,block_len_use,phase_blocks_cut_highconf_dir,phased_file_dir,HighConf_file_dir):
    use_cmd = "python3 " + code_path + "Cut_phaseblock_for_phased_h5_v4.0_highconf_v2.py" + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --block_threshold " + str(block_threshold) + " --block_len_use " + str(block_len_use)  + " --out_dir " +  phase_blocks_cut_highconf_dir + " --in_dir " + phased_file_dir + " --highconf_profile_dir " + HighConf_file_dir
    Popen(use_cmd,shell=True).wait()


def Get_fastq_files_total(bam_file,chr_start,chr_end,num_threads,Raw_fastqs_dir,Sorted_bam_dir,xin):
    use_cmd = "python3 " + code_path + "Read_fastqs_from_sortedbam_v2.py " + " --num_threads " + str(num_threads) + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Raw_fastqs_dir + " --bam_file " + bam_file + " --bam_dir " + Sorted_bam_dir
    Popen(use_cmd,shell=True).wait()


def Extract_reads_for_small_chunks(chr_start,chr_end,h5_dir,phase_blocks_cut_highconf_dir,Local_Assembly_dir,Raw_fastqs_dir,block_len_use,sample_name,num_threads,mole_len_dict_prev_file,mole_len_dict_file):
    use_cmd = "python3 " + code_path + "Run_extract_reads_for_smallchunks_all_lessmem_hybrid.py" + " --phase_cut_folder " + phase_blocks_cut_highconf_dir + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Local_Assembly_dir + " --h5_folder " + h5_dir +  " --fastq_folder " + Raw_fastqs_dir + " --cut_threshold " + str(block_len_use) + " --sample_name " + sample_name + " --num_threads " +   str(num_threads) + " --mole_len_dict_prev_file " + mole_len_dict_prev_file + " --mole_len_dict_file " + mole_len_dict_file 
    Popen(use_cmd,shell=True).wait()


def merge_h5_file(h5_file,fw,count):
    f = open(h5_file,"r")
    for line in f:
        data = line.rsplit()
        data[6] = str(count)
        count += 1
        fw.writelines("\t".join(data) + "\n")
    f.close() 
    return count


def extract_bam_by_chr(bam_file,chr_num,sample_bam_chr,xin):
    try:
        extract_bam_chr_cmd = "samtools view -b " + bam_file + " chr" + str(chr_num) + " > " + sample_bam_chr 
    except:
        extract_bam_chr_cmd = code_path + "samtools/" + "samtools view -b " + bam_file + " chr" + str(chr_num) + " > " + sample_bam_chr 

    Popen(extract_bam_chr_cmd,shell=True).wait()


def check_mole_len(cur_mole_file_prefix,chr_start,chr_end):
    mole_len_dict = defaultdict(int)
    for chr_num in range(chr_start,chr_end + 1): 
        cur_mole_file = cur_mole_file_prefix + "_chr" + str(chr_num)
        cat_cmd = "cat " + cur_mole_file + " | wc -l "
        p = Popen(cat_cmd,shell=True,stdout=PIPE, stderr=STDOUT)
        stdout, stderr = p.communicate()
        mole_len_dict[chr_num] = int(stdout.decode().rsplit()[0])
    return mole_len_dict


def add_two_dicts(mole_len_dict_prev,mole_len_dict):
    mole_len_dict_old = mole_len_dict_prev.copy()
    for key, val in mole_len_dict.items():
        mole_len_dict_prev[key] = val + mole_len_dict_old[key]
    return mole_len_dict_prev




def main():
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_step1.py -h",shell=True).wait()
    else:
        chr_start = args.chr_start
        chr_end = args.chr_end
        block_len_use = args.block_len_use
        block_threshold = args.block_threshold
        uniq_map_dir = args.uniq_map_dir + "/"
        num_threads = args.num_threads
        num_threads_for_samtools_sort = args.num_threads_for_samtools_sort
        HighConf_file_dir = args.out_dir + "/HighConf_file/"
        phased_file_dir = args.out_dir + "/results_phased_probmodel/"
        phase_blocks_cut_highconf_dir = args.out_dir + "/phase_blocks_cut_highconf/"
        Raw_fastqs_dir = args.out_dir + "/Raw_fastqs/"
        Local_Assembly_dir = args.out_dir + "/Local_Assembly_by_chunks/"
        merged_h5_dir = args.out_dir + "/H5_for_molecules/"
        if os.path.exists(merged_h5_dir):
            print("using existing output folder: " + merged_h5_dir)
        else:
            os.makedirs(merged_h5_dir)
        if os.path.exists(HighConf_file_dir):
            print("using existing output folder: " + HighConf_file_dir)
        else:
            os.makedirs(HighConf_file_dir)
        if os.path.exists(Local_Assembly_dir):
            print("using existing output folder: " + Local_Assembly_dir)
        else:
            os.makedirs(Local_Assembly_dir)


        #### Get fragment file and Merge for multi-libs ####
        _num = 0
        h5_file_prefix_list = []
        for _num in range(len(bam_list)):
            bam_file = bam_list[_num]
            vcf_file = vcf_list[_num]
            sample_name = sample_list[_num]
            _num += 1
            h5_dir = args.out_dir + "/H5_for_molecules_" + sample_name + "/"
            Get_fragment_files(bam_file,vcf_file,chr_start,chr_end,h5_dir,sample_name,num_threads)
            h5_file_prefix_list.append(h5_dir + sample_name + "_chr")

        for chr_num in range(chr_start,chr_end + 1):
            merged_h5_file = merged_h5_dir + "merged_chr" + str(chr_num)
            merged_h5_file_sorted = merged_h5_dir + "merged_chr" + str(chr_num) + "_sorted.h5"
            fw = open(merged_h5_file,"w")
            count =  1
            for one_file_prefix in h5_file_prefix_list:
                one_file = one_file_prefix + str(chr_num)
                count = merge_h5_file(one_file,fw,count)
            fw.close()
            sorted_cmd = "cat " + merged_h5_file + " | sort -k2n > " + merged_h5_file_sorted 
            print(sorted_cmd)
            Popen(sorted_cmd,shell=True).wait()
        
        
        #### HighConf file ####
        for chr_num in range(chr_start,chr_end + 1):
            merge_bam_chr = HighConf_file_dir + "merged_chr" + str(chr_num) + ".bam"  
            try:
                merge_bam_cmd = "samtools merge " + merge_bam_chr + " "
            except:
                merge_bam_cmd = code_path + "samtools/" + "samtools merge " + merge_bam_chr + " "

            _num = 0
            pool = Pool(len(bam_list))
            for _num in range(len(bam_list)):
                bam_file = bam_list[_num]
                sample_name = sample_list[_num]
                sample_bam_chr = HighConf_file_dir + sample_name + "_chr" + str(chr_num) + ".bam" 
                pool.apply_async(extract_bam_by_chr,(bam_file,chr_num,sample_bam_chr,"xin"))
                merge_bam_cmd += sample_bam_chr + " "
                _num += 1
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()

            Popen(merge_bam_cmd,shell=True).wait()
            try:
                indx_bam_chr_cmd = "samtools index " + merge_bam_chr
            except:
                indx_bam_chr_cmd = code_path + "samtools/" + "samtools index " + merge_bam_chr

            Popen(indx_bam_chr_cmd,shell=True).wait()
        
        for chr_num in range(chr_start,chr_end + 1):
            merge_bam_chr = HighConf_file_dir + "merged_chr" + str(chr_num) + ".bam"  
            Get_highconf_profile(merge_bam_chr,chr_num,chr_num,HighConf_file_dir,uniq_map_dir)

        #### Haplotying,Cuting PBs ####
        _num = 0
        all_depth_median = []
        for _num in range(len(bam_list)):
            sample_name = sample_list[_num]
            _num += 1
            h5_dir_sample_depth = args.out_dir + "/H5_for_molecules_" + sample_name + "/median_depth_for_var.txt"
            depth_file = open(h5_dir_sample_depth, 'r')
            depth_median = depth_file.readline()
            depth_file.close()
            all_depth_median.append(float(depth_median.rsplit()[0]))
            _num += 1
        depth_avg_libs = np.mean(all_depth_median)
        fw = open(merged_h5_dir + "median_depth_for_var.txt","w")
        fw.write(str(depth_avg_libs) + "\n")
        fw.close()

        sample_name = "merged"
        Haplotying_fragments(chr_start,chr_end,phased_file_dir,merged_h5_dir,sample_name)
        Cut_phase_blocks(chr_start,chr_end,block_threshold,block_len_use,phase_blocks_cut_highconf_dir,phased_file_dir,HighConf_file_dir)

        #### Get fastqs ####
        _num = 0
        pool = Pool(len(bam_list))
        for _num in range(len(bam_list)):
            bam_file = bam_list[_num]
            sample_name = sample_list[_num]
            Raw_fastqs_dir_sample = args.out_dir + "/Raw_fastqs_" + sample_name + "/"
            Bam_dir_sample = args.out_dir + "/sorted_bam_" + sample_name + "/"
            pool.apply_async(Get_fastq_files_total,(bam_file,chr_start,chr_end,num_threads_for_samtools_sort,Raw_fastqs_dir_sample,Bam_dir_sample,"xin"))
            _num += 1
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()

        #### Extract for small chunks ####
        _num = 0
        mole_len_dict_prev = defaultdict(int)
        for chr_num in range(chr_start,chr_end + 1):
            mole_len_dict_prev[chr_num] = 0
        mole_len_dict_prev_file = args.out_dir + "/mole_len_dict_prev.p"
        pickle.dump(mole_len_dict_prev,open(mole_len_dict_prev_file,"wb"))
        for _num in range(len(bam_list)):
            bam_file = bam_list[_num]
            sample_name = sample_list[_num]
            Raw_fastqs_dir_sample = args.out_dir + "/Raw_fastqs_" + sample_name + "/"
            h5_dir_sample = args.out_dir + "/H5_for_molecules_" + sample_name + "/"
            Local_Assembly_dir_sample = args.out_dir + "/Local_Assembly_by_chunks_" + sample_name  + "/"
            mole_len_dict = check_mole_len(h5_dir_sample + sample_name,chr_start,chr_end) 
            mole_len_dict_file = args.out_dir +  "/mole_len_dict.p"
            pickle.dump(mole_len_dict,open(mole_len_dict_file,"wb"))
            Extract_reads_for_small_chunks(chr_start,chr_end,h5_dir_sample,phase_blocks_cut_highconf_dir,Local_Assembly_dir_sample,Raw_fastqs_dir_sample,block_len_use,sample_name,12,mole_len_dict_prev_file,mole_len_dict_file)
            mole_len_dict_prev = add_two_dicts(mole_len_dict_prev,mole_len_dict)
            pickle.dump(mole_len_dict_prev,open(mole_len_dict_prev_file,"wb"))
            _num += 1


        #### Merge ####
        pool = Pool(chr_end - chr_start + 1)
        for chr_num in range(chr_start,chr_end + 1):
            _num = 0
            in_dir_list = ""
            local_assembly_dir_chr = Local_Assembly_dir + "chr" + str(chr_num) +  "_files_cutPBHC/" 
            for _num in range(len(bam_list)):
                sample_name = sample_list[_num]
                in_dir = args.out_dir + "/Local_Assembly_by_chunks_" + sample_name  + "/"  + "chr" + str(chr_num) +  "_files_cutPBHC/" + ","
                in_dir_list += in_dir
                _num += 1
            pool.apply_async(Merge_reads_for_PB_of_merged_libs(in_dir_list,local_assembly_dir_chr,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()



if __name__ == "__main__":
    main()



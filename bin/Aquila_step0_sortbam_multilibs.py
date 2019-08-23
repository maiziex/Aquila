#!/usr/bin/env python
from subprocess import Popen
from argparse import ArgumentParser
import os
import sys
from multiprocessing import Pool,cpu_count,active_children,Manager
import time
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
__author__ = "Xin Zhou@Stanford"

parser = ArgumentParser(description="sort bam by qname:")
parser.add_argument('--bam_file_list','-bam',help="Required Parameter, BAM file list, each BAM file is seperately by comma \",\". For example: 1.bam,2.bam",required=True)
parser.add_argument('--out_dir','-o', help="Directory to store Aquila assembly results, default =./Assembly_results", default="./Asssembly_results")
parser.add_argument('--sample_name_list','-sl', help='Required Parameter, The sample names list, each sample name is seperately by comma ",". For example: S24385_lysis_2,S24385_lysis_2H', type=str,required=True)
parser.add_argument('--num_threads_for_samtools_sort','-t',type=int, help=" The number of threads you can define for samtools sort, default = 20",default=20)
args = parser.parse_args()

def sort_start(sort_bam_cmd,idx_bam_cmd,xin):
    Popen(sort_bam_cmd,shell=True).wait()
    Popen(idx_bam_cmd,shell=True).wait()


def main():
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_sortbam_hybrid.py -h",shell=True).wait()
    else:
        sample_list = [item for item in args.sample_name_list.split(',')]
        bam_list = [item for item in args.bam_file_list.split(',')]
        out_dir = args.out_dir
        num_threads = int(args.num_threads_for_samtools_sort)
        if os.path.exists(out_dir):
            print("using existing output folder: " + out_dir)
        else:
            os.makedirs(out_dir)

        count = 1
        total_num = len(sample_list)
        pool = Pool(processes=total_num)
        for ii in range(len(sample_list)):
            sample_name = sample_list[ii]
            bam_file = bam_list[ii]
            bam_file_dir = out_dir + "/sorted_bam_" + sample_name + "/"
            if os.path.exists(bam_file_dir):
                print("using existing output folder: " + bam_file_dir)
            else:
                os.makedirs(bam_file_dir)

            bam_sorted_file = bam_file_dir + "sorted_bam.bam"
            try:
                sort_bam_cmd = "samtools sort -@ " + str(num_threads) + " -n " + bam_file + " -o " + bam_sorted_file
            except:
                sort_bam_cmd = code_path + "samtools/" + "samtools sort -@ " + str(num_threads) + " -n " + bam_file + " -o " + bam_sorted_file
                
            idx_bam_cmd = "touch "  + bam_file_dir + "finish_bam.txt"
            pool.apply_async(sort_start,(sort_bam_cmd,idx_bam_cmd,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()
        print("Sorting bam finished...")



if __name__ == "__main__":
    main()


#!/usr/bin/env python
from subprocess import Popen
from argparse import ArgumentParser
import os
import sys
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
parser = ArgumentParser(description="sort bam by qname:")
parser.add_argument('--bam_file','-b', help="bam file",required=True)
parser.add_argument('--out_dir','-o', help="output folder")
parser.add_argument('--num_threads_for_bwa_mem','-t',type=int, help=" The number of threads you can define for samtoos sort",default=20)
args = parser.parse_args()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_sortbam.py -h",shell=True).wait()
    else:
        bam_file = args.bam_file
        out_dir = args.out_dir
        num_threads = int(args.num_threads_for_bwa_mem)
        if os.path.exists(out_dir):
            print("using existing output folder: " + out_dir)
        else:
            os.makedirs(out_dir)

        if os.path.exists(out_dir + "/sorted_bam/"):
            print("using existing output folder: " + out_dir + "/sorted_bam/")
        else:
            os.makedirs(out_dir + "/sorted_bam/")

        bam_sorted_file = out_dir + "/sorted_bam/sorted_bam.bam"
        try:
            sort_bam_cmd = "samtools sort -@ " + str(num_threads) + " -n " + bam_file + " -o " + bam_sorted_file
        except:
            sort_bam_cmd = code_path + "samtools/" + "samtools sort -@ " + str(num_threads) + " -n " + bam_file + " -o " + bam_sorted_file

        Popen(sort_bam_cmd,shell=True).wait()
        idx_bam_cmd = "touch "  + out_dir + "finish_bam.txt"
        Popen(idx_bam_cmd,shell=True).wait()

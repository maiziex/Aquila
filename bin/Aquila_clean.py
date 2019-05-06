import pdb
#pdb.set_trace()
from argparse import ArgumentParser
import os
import sys
import time
from multiprocessing import Pool,cpu_count,active_children,Manager
from subprocess import Popen
import glob
__author__ = "Xin Zhou@Stanford"
parser = ArgumentParser(description="Clean Assembly Data by Aquila:")
parser.add_argument('--assembly_dir','-i', help="assembly folder",required=True)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--num_of_threads','-nt', help="number of threads",default= 10)
args = parser.parse_args()

def del_one_file(one_file,xin):
    rm_cmd = "rm -rf " + one_file
    Popen(rm_cmd,shell=True).wait()


def Aquila_clean(in_dir,chr_start,chr_end,num_of_threads):
    local_assembly_dir = in_dir + "Local_Assembly_by_chunks/"
    for chr_num in range(chr_start,chr_end + 1):
        one_dir = local_assembly_dir + "chr" + str(chr_num)  + "_files_cutPBHC/"
        all_files = sorted(glob.glob(one_dir +  "fastq_by_*_assembly"),key=os.path.getsize,reverse=True)
        total_num = len(all_files)
        pool = Pool(num_of_threads)
        count = 1
        for one_file in all_files:
            count += 1
            pool.apply_async(del_one_file,(one_file,"xin"))
            if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
                pool.close()
                while len(active_children()) > 1:
                    time.sleep(0.5)
                pool.join()
                if (count - 1) == total_num:
                    print("finished chr" + str(chr_num))
                else:
                    pool = Pool(num_of_threads)
        rm_cmd = "rm -rf " + one_dir  + "*.fastq"
        Popen(rm_cmd,shell=True).wait()

    print("Cleaning local assembly by chunks files done!")

    Raw_fastqs_dir = in_dir + "Raw_fastqs/"
    rm_cmd = "rm -rf " + Raw_fastqs_dir
    Popen(rm_cmd,shell=True).wait()
    print("Cleaning Raw_fastqs done!")

    bam_dir = in_dir + "sorted_bam/"
    rm_cmd = "rm -rf " + bam_dir
    Popen(rm_cmd,shell=True).wait()
    print("Cleaning sorted_bam done!")

    hc_dir = in_dir + "HighConf_file/"
    rm_cmd_1 = "rm -rf " + hc_dir + "*.bam" 
    rm_cmd_2 = "rm -rf " + hc_dir + "*.txt" 
    rm_cmd_3 = "rm -rf " + hc_dir + "*pos_cf.p" 
    Popen(rm_cmd_1,shell=True).wait()
    Popen(rm_cmd_2,shell=True).wait()
    Popen(rm_cmd_3,shell=True).wait()
    print("Cleaning HighConf file done!")
 


if __name__ == "__main__":
    in_dir = args.assembly_dir + "/"
    num_of_threads = int(args.num_of_threads)
    chr_start = args.chr_start
    chr_end = args.chr_end
    local_assembly_dir = in_dir + "Local_Assembly_by_chunks/"
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_clean.py -h",shell=True).wait()
    else:
        Aquila_clean(in_dir,chr_start,chr_end,num_of_threads)



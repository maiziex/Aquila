#!/usr/bin/env python
import pdb
#pdb.set_trace()
import glob
import os
from argparse import ArgumentParser
from multiprocessing import Pool,cpu_count,active_children,Manager
from subprocess import Popen
import time
from Concatenate_contigs_all_v4_extend_for_HCbk import *
import multiprocessing
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
parser = ArgumentParser(description="run local assembly by spades:")
parser.add_argument('--num_threads','-nt', help="number of threads")
parser.add_argument('--out_dir','-o', help="out dir")
parser.add_argument('--minicontig_dir','-minicontig_o', help="minicontig dir")
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=22)
args = parser.parse_args()


def use_spades(one_file_fastq,out_dir,xin):
    try:
        use_cmd = code_path + "SPAdes-3.13.0-Linux/" + "spades.py -t 5 --only-assembler --12 " + one_file_fastq + " -o " + out_dir 
    except:
        use_cmd = "spades.py -t 5 --only-assembler --12 " + one_file_fastq + " -o " + out_dir 

    Popen(use_cmd,shell=True).wait()


def del_one_dir_1(out_dir,xin):
    temp_dir = out_dir + "_temp/"
    use_cmd_1 = "mv " + out_dir + " " + temp_dir
    Popen(use_cmd_1,shell=True).wait()
    use_cmd_2 = "mkdir " + out_dir
    Popen(use_cmd_2,shell=True).wait()
    contig_file = temp_dir + "contigs.fasta"
    if os.path.exists(contig_file):
        use_cmd_3 = "mv " + contig_file + " " + out_dir
        Popen(use_cmd_3,shell=True).wait()
    rm_cmd = "rm -rf " + temp_dir
    Popen(rm_cmd,shell=True).wait()


def del_one_dir_2(out_dir,xin):
    rm_cmd = "rm -rf " + out_dir
    Popen(rm_cmd,shell=True).wait()


def run_spades_all(chr_start,chr_end,output_dir,num_of_threads,minicontig_dir):
    num_of_threads = multiprocessing.cpu_count() - 10
    for chr_num in range(chr_start, chr_end+1):
        in_dir = output_dir + "chr" + str(chr_num) + "_files_cutPBHC/" 
        count = 1
        fastq_files_all = sorted(glob.glob(in_dir +  "fastq_by_*.fastq"),key=os.path.getsize,reverse=True)
        total_num = len(fastq_files_all)
        pool = Pool(num_of_threads)
        out_dir_list = []
        for one_file_fastq in fastq_files_all:
            one_file = one_file_fastq[:-6]
            out_dir = one_file + "_spades_assembly"
            spades_contig_file = out_dir + "/" + "contigs.fasta"
            if os.path.exists(spades_contig_file):
                count += 1
                #print("using existing " + spades_contig_file)
                out_dir_list.append(out_dir)
            else:
                count += 1
                pool.apply_async(use_spades,(one_file_fastq,out_dir,"xin"))
                out_dir_list.append(out_dir)
            if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
                pool.close()
                while len(active_children()) > 1:
                    time.sleep(0.5)
                pool.join()
                
                # delete temp files
                """
                if out_dir_list != []:
                    pool = Pool(num_of_threads)
                    for one_dir in out_dir_list:
                        pool.apply_async(del_one_dir_1,(one_dir,"xin"))
                    pool.close()
                    while len(active_children()) > 1:
                        time.sleep(0.5)
                    pool.join()
                """

                if (count - 1) == total_num:
                    print("finished chr" + str(chr_num))
                else:
                    pool = Pool(num_of_threads)
                    out_dir_list = []
        output_file = "Aquila_cutPBHC_minicontig" + "_chr" + str(chr_num) + ".fasta"
        Concatenate_start(in_dir,minicontig_dir,output_file,"xin")   
        
        # delete assembly files
        """
        time.sleep(5)
        pool = Pool(num_of_threads)
        count = 1
        for one_file_fastq in fastq_files_all:
            one_file = one_file_fastq[:-6]
            out_dir = one_file + "_spades_assembly"
            count += 1
            pool.apply_async(del_one_dir_2,(out_dir,"xin"))
            if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
                pool.close()
                while len(active_children()) > 1:
                    time.sleep(0.5)
                pool.join()
                if (count - 1) == total_num:
                    print("finished chr" + str(chr_num))
                else:
                    pool = Pool(num_of_threads)
        time.sleep(5)
        """
    print("All Done~")




if __name__ == "__main__":
    output_dir = args.out_dir
    minicontig_dir = args.minicontig_dir
    chr_start = args.chr_start
    chr_end = args.chr_end
    num_of_threads = int(args.num_threads)
    run_spades_all(chr_start,chr_end,output_dir,num_of_threads,minicontig_dir)


#!/usr/bin/env python
import pdb
#pdb.set_trace()
from subprocess import Popen
from argparse import ArgumentParser
import os
import sys
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
__author__ = "Xin Zhou@Stanford"
parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--bam_file','-bam',help="Required parameter; BAM file, called by longranger align",required=True)
parser.add_argument('--vcf_file','-v',help="Required parameter; VCF file, called by FreeBayes",required=True)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from, default = 1", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by,default = 23", default=23)
parser.add_argument('--sample_name','-name',help="Required parameter; Sample Name you can define by yourself, for example: S12878",required=True)
parser.add_argument('--out_dir','-o', help="Directory to store assembly results, default = ./Assembly_results",default="./Asssembly_results")
parser.add_argument('--uniq_map_dir','-uniq_dir', help="Required Parameter; Directory for 100-mer uniqness, run ./install to download \"Uniquess_map\" for hg38",required=True)
parser.add_argument('--num_threads','-t_chr',type=int,help="number of threads, default = 8 (recommended)", default=8)
parser.add_argument('--num_threads_for_bwa_mem','-t',type=int,help="number of threads for bwa-mem, default = 20", default=20)
parser.add_argument('--num_threads_for_extract_reads','-t_extract',type=int,help="number of threads for extracting raw reads, default = 8 (recommended)", default=8)
parser.add_argument('--block_threshold','-bt',type=int,help="phase block threshold, default = 200000",default=200000)
parser.add_argument('--block_len_use','-bl',type=int,help="phase block len threshold, default = 100000",default=100000)
parser.add_argument('--mbq_threshold','-mbq',type=int,help="phred-scaled quality score for the assertion made in ALT, default = 13",default=13)
parser.add_argument('--boundary','-bound',type=int,help="boundary for long fragments with the same barcode, default = 50000 (recommended)",default=50000)

args = parser.parse_args()


def Get_fragment_files(bam_file,vcf_file,chr_start,chr_end,h5_dir,num_threads,mbq,boundary):
    use_cmd = "python3 " + code_path + "Run_h5_all_multithreads.py" + " --bam_file " + bam_file + " --vcf_file " + vcf_file + " --sample_name " + sample_name + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --mbq " + str(mbq) + "  --mmq 0  --boundary " + str(boundary) + " --num_threads " + str(num_threads) + " --out_dir " + h5_dir
    Popen(use_cmd,shell=True).wait()


def Get_highconf_profile(bam_file,chr_start,chr_end,HighConf_file_dir,uniq_map_dir):
    use_cmd = "python3 " + code_path + "Generate_highconf_cut_profile_v2.py" + " --bam_file " +  bam_file + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --mbq 13 --mmq 20 " + " --out_dir " + HighConf_file_dir + " --uniq_map_dir " + uniq_map_dir
    Popen(use_cmd,shell=True).wait()
  

def Haplotying_fragments(chr_start,chr_end,phased_file_dir,h5_dir,sample_name):
    use_cmd = "python3 " + code_path + "Run_phase_alg_multithreads2.py" + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --overlap_threshold 3 --support_threshold 5 " + " --out_dir " + phased_file_dir  + " --h5_dir " + h5_dir + " --sample_name " + sample_name
    Popen(use_cmd,shell=True).wait()


def Cut_phase_blocks(chr_start,chr_end,block_threshold,block_len_use,phase_blocks_cut_highconf_dir,phased_file_dir,HighConf_file_dir):
    use_cmd = "python3 " + code_path + "Cut_phaseblock_for_phased_h5_v4.0_highconf_v2.py" + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --block_threshold " + str(block_threshold) + " --block_len_use " + str(block_len_use)  + " --out_dir " +  phase_blocks_cut_highconf_dir + " --in_dir " + phased_file_dir + " --highconf_profile_dir " + HighConf_file_dir
    Popen(use_cmd,shell=True).wait()


def Get_fastq_files_total(bam_file,chr_start,chr_end,num_threads,Raw_fastqs_dir,Sorted_bam_dir):
    use_cmd = "python3 " + code_path + "Read_fastqs_from_sortedbam_v2.py " + " --num_threads " + str(num_threads) + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Raw_fastqs_dir + " --bam_file " + bam_file + " --bam_dir " + Sorted_bam_dir
    Popen(use_cmd,shell=True).wait()


def Extract_reads_for_small_chunks(chr_start,chr_end,h5_dir,phase_blocks_cut_highconf_dir,Local_Assembly_dir,Raw_fastqs_dir,block_len_use,sample_name,num_threads):
    use_cmd = "python3 " + code_path + "Run_extract_reads_for_smallchunks_all_lessmem.py" + " --phase_cut_folder " + phase_blocks_cut_highconf_dir + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Local_Assembly_dir + " --h5_folder " + h5_dir +  " --fastq_folder " + Raw_fastqs_dir + " --cut_threshold " + str(block_len_use) + " --sample_name " + sample_name + " --num_threads " +   str(num_threads)
    Popen(use_cmd,shell=True).wait()




if __name__ == "__main__":
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_step1.py -h",shell=True).wait()
    else:
        bam_file = args.bam_file
        vcf_file = args.vcf_file
        chr_start = args.chr_start
        chr_end = args.chr_end
        block_len_use = args.block_len_use
        block_threshold = args.block_threshold
        mbq = args.mbq_threshold
        boundary = args.boundary
        uniq_map_dir = args.uniq_map_dir + "/"
        num_threads = int(args.num_threads)
        num_threads_for_bwa_mem = int(args.num_threads_for_bwa_mem)
        num_threads_for_extract_reads = int(args.num_threads_for_extract_reads)
        sample_name = args.sample_name
        h5_dir = args.out_dir + "/H5_for_molecules/"
        HighConf_file_dir = args.out_dir + "/HighConf_file/"
        phased_file_dir = args.out_dir + "/results_phased_probmodel/"
        phase_blocks_cut_highconf_dir = args.out_dir + "/phase_blocks_cut_highconf/"
        Raw_fastqs_dir = args.out_dir + "/Raw_fastqs/"
        Sorted_bam_dir = args.out_dir + "/sorted_bam/"
        Local_Assembly_dir = args.out_dir + "/Local_Assembly_by_chunks/"
        
        Get_fragment_files(bam_file,vcf_file,chr_start,chr_end,h5_dir,num_threads,mbq,boundary)
        Get_highconf_profile(bam_file,chr_start,chr_end,HighConf_file_dir,uniq_map_dir)
        Haplotying_fragments(chr_start,chr_end,phased_file_dir,h5_dir,sample_name)
        Cut_phase_blocks(chr_start,chr_end,block_threshold,block_len_use,phase_blocks_cut_highconf_dir,phased_file_dir,HighConf_file_dir)
        
        Get_fastq_files_total(bam_file,chr_start,chr_end,num_threads_for_bwa_mem,Raw_fastqs_dir,Sorted_bam_dir)
        Extract_reads_for_small_chunks(chr_start,chr_end,h5_dir,phase_blocks_cut_highconf_dir,Local_Assembly_dir,Raw_fastqs_dir,block_len_use,sample_name,num_threads_for_extract_reads)
    
    

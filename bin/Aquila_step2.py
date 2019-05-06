#!/usr/bin/env python
import pdb
#pdb.set_trace()
from subprocess import Popen
from argparse import ArgumentParser
import os
import sys
import pickle
import os.path
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
__author__ = "Xin Zhou@Stanford"
parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--out_dir','-o', help="Directory to store output",default="./Results")
parser.add_argument('--reference','-ref', help="reference fasta file",required=True)
parser.add_argument('--num_threads','-t',type=int,help="number of threads", default=30)
parser.add_argument('--num_threads_spades','-t',type=int,help="number of threads for spades", default=5)
parser.add_argument('--block_len_use','-bl',type=int,help="phase block len threshold",default=100000)

args = parser.parse_args()

def read_ref(fasta_file,chr_num,out_dir):
    f = open(fasta_file,"r")
    count = 0
    ref_seq = ""
    for line in f:
        if count > 0:
            data = line.rsplit()
            ref_seq += data[0]
        count += 1
    print("total_len for chr" + str(chr_num))
    print(len(ref_seq))
    pickle.dump(ref_seq, open(out_dir + "ref_seq_chr" + str(chr_num) +  ".p","wb"))


def extract_ref_chr(ref_file,chr_num,out_dir):
    fw = open(out_dir + "genome_ref_chr" + str(chr_num) + ".fasta","w")
    f = open(ref_file,"r")
    flag = 0
    total_len = 0
    for line in f:
        data = line.rsplit()
        if data[0] == ">chr" + str(chr_num):
            fw.writelines(">" + str(chr_num) + "\n")
            flag = 1
        elif data[0][0] == ">" and flag == 1:
            break
        else:
            if flag == 1:
                total_len += len(data[0])
                fw.writelines(data[0] + "\n")
    print("chr" + str(chr_num) + ":")
    print(total_len)


def local_assembly_for_small_chunks(chr_start,chr_end,num_threads,num_threads_spades,Local_Assembly_dir,Assembly_Contigs_dir):
    use_cmd = "python " + code_path + "Run_spades_final_MT_2_all_noec_deltemp.py" + " --num_threads " + str(num_threads) + " --num_threads_spades " + str(num_threads_spades) + " --chr_start " + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Local_Assembly_dir + " --minicontig_dir " + Assembly_Contigs_dir
    Popen(use_cmd,shell=True).wait()


def Complete_contiguity(chr_start,chr_end,Assembly_Contigs_dir,phase_blocks_cut_highconf_dir,cut_threshold,ref_dir,ref_idx_file):
    use_cmd = "python " + code_path + "Concatenate_contigs_from_microcontigs.py" + " --chr_start "  + str(chr_start) + " --chr_end " + str(chr_end) + " --out_dir " + Assembly_Contigs_dir + " --phase_cut_folder " + phase_blocks_cut_highconf_dir + " --cut_threshold " + str(cut_threshold) + " --ref_file " + ref_idx_file + " --ref_dir " + ref_dir
    Popen(use_cmd,shell=True).wait()
  


if __name__ == "__main__":
    if len(sys.argv) == 1:
        Popen("python " + "Aquila_step2.py -h",shell=True).wait()
    else:
        chr_start = args.chr_start
        chr_end = args.chr_end
        ref_file = args.reference
        cut_threshold = args.block_len_use
        num_threads = int(args.num_threads)
        num_threads_spades = int(args.num_threads_spades)
        ref_dir = args.out_dir + "/ref_dir/"
        if os.path.exists(ref_dir):
            print("using existing output folder: " + ref_dir)
        else:
            os.makedirs(ref_dir)
        if ~os.path.exists(ref_dir + "ref.mmi"):
            try:
                mk_ref_idx = "minimap2 -d " + ref_dir + "ref.mmi "  + ref_file
            except:
                mk_ref_idx = code_path + "minimap2/" + "minimap2 -d " + ref_dir + "ref.mmi "  + ref_file

            Popen(mk_ref_idx,shell=True).wait()
        for chr_num in range(chr_start,chr_end + 1):
            extract_ref_chr(ref_file,chr_num,ref_dir)
        extract_ref_chr(ref_file,"X",ref_dir)
        for chr_num in range(chr_start,chr_end + 1):
            read_ref(ref_dir + "genome_ref_chr" + str(chr_num) + ".fasta",chr_num,ref_dir)
        read_ref(ref_dir + "genome_ref_chrX.fasta",23,ref_dir)

        ref_idx_file = ref_dir + "ref.mmi"
        Assembly_Contigs_dir = args.out_dir + "/Assembly_Contigs_files/"
        phase_blocks_cut_highconf_dir = args.out_dir + "/phase_blocks_cut_highconf/"
        Local_Assembly_dir = args.out_dir + "/Local_Assembly_by_chunks/"

        local_assembly_for_small_chunks(chr_start,chr_end,num_threads,num_threads_spades,Local_Assembly_dir,Assembly_Contigs_dir)
        Complete_contiguity(chr_start,chr_end,Assembly_Contigs_dir,phase_blocks_cut_highconf_dir,cut_threshold,ref_dir,ref_idx_file)
    

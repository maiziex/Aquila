import pdb
from subprocess import Popen
from argparse import ArgumentParser
import os
import sys

parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--assembly_dir','-i', help="Directory to store assembly results",required=True)

args = parser.parse_args()

def Cat_Contigs_for_wgs(assembly_dir,chr_start,chr_end):
    in_dir = assembly_dir + "Assembly_Contigs_files/"
    for chr_num in range(chr_start,chr_end + 1):
        one_contig = in_dir + "Aquila_Contig_chr" + str(chr_num) + ".fasta"
        cat_cmd = "cat " + one_contig  + " >> " + in_dir + "Aquila_Contig_WGS.fasta"
        Popen(cat_cmd,shell=True).wait()
    print("all done~")



if __name__ == "__main__":
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_Contigs_WGS.py -h",shell=True).wait()
    else:
        chr_start = args.chr_start
        chr_end = args.chr_end
        assembly_dir = args.assembly_dir + "/"
        Cat_Contigs_for_wgs(assembly_dir,chr_start,chr_end)










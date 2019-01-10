#!/usr/bin/env python
import pdb
#pdb.set_trace()
import pysam
from collections import defaultdict
from subprocess import Popen
from argparse import ArgumentParser
import os

parser = ArgumentParser(description="extract fastq files from bam:")
parser.add_argument('--out_dir','-o', help="output folder")
parser.add_argument('--bam_file','-bam', help="bam file")
parser.add_argument('--chr_start','-start',type=int,help="chr start", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chr end", default=23)
parser.add_argument('--num_threads','-t',type=int,help="number of threads", default=20)
args = parser.parse_args()
out_dir = args.out_dir
chr_start = int(args.chr_start)
chr_end = int(args.chr_end)

if os.path.exists(out_dir):
    print("using existing output folder: " + out_dir)
else:
    os.makedirs(out_dir)

chr_dict = {"chr1":1,"chr2":2,"chr3":3,"chr4":4,"chr5":5,"chr6":6,"chr7":7,"chr8":8,"chr9":9,"chr10":10,"chr11":11,"chr12":12,"chr13":13,"chr14":14,"chr15":15,"chr16":16,"chr17":17,"chr18":18,"chr19":19,"chr20":20,"chr21":21,"chr22":22,"chrX":23}

fw_curr = defaultdict()
for chr_num in range(chr_start, chr_end+1):
    fw_curr[chr_num] = open(out_dir + "fastq_by_Chr_" + str(chr_num),"w")

tab = str.maketrans("ACTGN", "TGACN")

def reverse_complement(seq):
    return seq.translate(tab)[::-1]


def write_pair_reads(prev_read,read,use_chr_num):
    prev_read_seq = prev_read.seq
    prev_read_qual = prev_read.qual
    if prev_read.is_reverse:
        prev_read_seq = reverse_complement(prev_read_seq)
        prev_read_qual = prev_read_qual[::-1]
    read_seq = read.seq
    read_qual = read.qual
    if read.is_reverse:
        read_seq = reverse_complement(read_seq)
        read_qual = read_qual[::-1]
    fw_use = fw_curr[use_chr_num]
    if prev_read.is_read1:
        cur_line = "@" + prev_read.qname + "\n" + prev_read_seq + "\n" + "+\n" + prev_read_qual + "\n" + "@" + read.qname + "\n" + read_seq + "\n" + "+\n" + read_qual + "\n"
        fw_use.write(cur_line)
    elif prev_read.is_read2:
        cur_line = "@" + read.qname + "\n" + read_seq + "\n" + "+\n" + read_qual + "\n" + "@" + prev_read.qname + "\n" + prev_read_seq + "\n" + "+\n" + prev_read_qual + "\n"
        fw_use.write(cur_line)


def read_fastqs_from_sorted_bam(sorted_bam_file,chr_start,chr_end):
    sam_file = pysam.AlignmentFile(sorted_bam_file, "rb")
    count = 0
    count_unmapped = 0
    count_line = 0
    use_chr_dict = defaultdict(int)
    for ii in range(chr_start,chr_end + 1):
        use_chr_dict[ii] == 1
    for read in sam_file.fetch(until_eof=True):
        if not read.is_secondary:
            count += 1
            if count == 2:
                if prev_read.is_unmapped and read.is_unmapped:
                    count_unmapped += 1
                else:
                    chr_name_1 = prev_read.reference_name
                    chr_name_2 = read.reference_name
                    if chr_name_1 == chr_name_2:
                        if chr_name_1 in chr_dict:
                            use_chr_num = chr_dict[chr_name_1]
                            if use_chr_num in use_chr_dict:
                                write_pair_reads(prev_read,read,use_chr_num)
                    else:
                        if chr_name_1 in chr_dict:
                            use_chr_num = chr_dict[chr_name_1]
                            if use_chr_num in use_chr_dict:
                                write_pair_reads(prev_read,read,use_chr_num)
                        elif chr_name_2 in chr_dict:
                            use_chr_num = chr_dict[chr_name_2]
                            if use_chr_num in use_chr_dict:
                                write_pair_reads(prev_read,read,use_chr_num)
                count = 0
            if count == 1:
                prev_read = read


def close_file(chr_start,chr_end):
    for chr_num in range(chr_start, chr_end+1):
        fw_curr[chr_num].close()




if __name__ == "__main__":
    bam_file = args.bam_file
    num_threads = int(args.num_threads)
    bam_sorted_file = out_dir + "sorted_bam.bam"
    sort_bam_cmd = "samtools sort -@ " + str(num_threads) + " -n " + bam_file + " -o " + bam_sorted_file
    #Popen(sort_bam_cmd,shell=True).wait()
    read_fastqs_from_sorted_bam(bam_sorted_file,chr_start,chr_end)
    close_file(chr_start,chr_end)


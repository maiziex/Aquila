import pdb
#pdb.set_trace()
import glob
import os
import numpy as np
from argparse import ArgumentParser
import pickle


def concatenate_contigs(contig_file,fw,contig_num,contigs_all,hp_flag,PS_flag,PS_flag_2):
    f = open(contig_file,"r")
    contig_ = ""
    flag = 0
    for line in f:
        data = line.rsplit()
        if data[0][0] == ">":
            if flag == 1:
                fw.writelines(">" + str(contig_num) + "_PS" +  str(PS_flag) + ":" + str(PS_flag_2) + "_hp" + str(hp_flag) + "\n")
                fw.writelines(contig_ + "\n")
                contigs_all.append(len(contig_))
                contig_num += 1
                contig_ = ""
                flag = 0
            flag = 1
        else:
            string_ = data[0]
            contig_ += string_

    fw.writelines(">" + str(contig_num) + "_PS" +  str(PS_flag) + ":" + str(PS_flag_2) + "_hp" + str(hp_flag) + "\n")
    fw.writelines(contig_ + "\n")
    contig_num += 1
    return (fw,contig_num,contigs_all)


def concatenate_all(input_dir,output_filename):
    count = 0
    count_2 = 0
    count_3 = 0
    contig_num = 1
    contigs_all = []

    fw_contigs_all = open(output_filename,"w")
    for one_file in glob.glob(input_dir + "fastq_by_*_spades_assembly"):
        #print(one_file)
        count += 1
        spades_contig_file = one_file + "/" + "contigs.fasta"
        if os.path.exists(spades_contig_file):
            if "hp1" in one_file:
                hp_flag = 1
            else:
                hp_flag = 2
            count_2 += 1
            PS_flag = int(one_file[len(input_dir):].split("fastq_by_")[1].split("_")[0])
            PS_flag_2 = int(one_file[len(input_dir):].split("fastq_by_")[1].split("_")[1])
            fw_contigs_all, contig_num,contigs_all = concatenate_contigs(spades_contig_file,fw_contigs_all,contig_num,contigs_all,hp_flag,PS_flag,PS_flag_2)
        else:
            count_3 += 1

    print(count,count_2,count_3)
    print(np.mean(contigs_all))
    print(np.median(contigs_all))
    print(np.max(contigs_all))
    print(np.min(contigs_all))

    contigs_all_sorted = sorted(contigs_all,reverse=True)
    total_contigs_len = sum(contigs_all)
    total_contigs_len_half = float(total_contigs_len)/2
    cumu_len = 0
    for one_len in contigs_all_sorted:
        cumu_len += one_len
        if cumu_len >= total_contigs_len_half:
            n50_len = one_len
            print("------results:-------")
            print("n50 minicontigs: " + str(n50_len))
            break

    fw_contigs_all.close()    


def Concatenate_start(input_dir,output_dir,output_file,xin): 
    if os.path.exists(output_dir):
        print("using existing output folder: " + output_dir)
    else:
        os.makedirs(output_dir)

    output_filename = output_dir + output_file
    concatenate_all(input_dir,output_filename)


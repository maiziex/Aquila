import pdb
#pdb.set_trace()
import glob
from collections import defaultdict
import os
from subprocess import Popen


def Merge_reads_for_PB_of_merged_libs(in_dir_list,out_dir,xin):
    if os.path.exists(out_dir):
        Popen("rm -rf " + out_dir,shell=True).wait()
        os.makedirs(out_dir)
    else:
        os.makedirs(out_dir)
    in_dir_list_raw = in_dir_list.split(",")
    in_dir_1 = in_dir_list_raw[0]
    fastq_files_all_1 = glob.glob(in_dir_1 +  "fastq_by_*.fastq")
    file_all_1 = []
    for one_file_raw in fastq_files_all_1:
        one_file = one_file_raw.split("/")[-1]
        file_all_1.append(one_file)
    print(len(file_all_1))
    merge_flag = defaultdict(int)
    count_merge = 0
    count_uniq_1 = 0 
    count_uniq_2 = 0
    for in_dir in in_dir_list_raw[1:]:
        fastq_files_all = glob.glob(in_dir +  "fastq_by_*.fastq")
        file_all = []
        for one_file_raw in fastq_files_all:
            one_file = one_file_raw.split("/")[-1]
            file_all.append(one_file)
            if one_file in file_all_1:
                rn_cmd = "mv " + in_dir_1 + one_file + " " + in_dir_1 + "temp.fastq"
                Popen(rn_cmd,shell=True).wait()
                cat_cmd = "cat " + in_dir_1 + "temp.fastq" + " " + in_dir + one_file + " > " + in_dir_1 + one_file
                Popen(cat_cmd,shell=True).wait()
            else:
                mv_cmd = "mv " + in_dir + one_file + " " + in_dir_1
                Popen(mv_cmd,shell=True).wait()

        ### update "in_dir_1" ###
        fastq_files_all_1 = glob.glob(in_dir_1 +  "fastq_by_*.fastq")
        file_all_1 = []
        for one_file_raw in fastq_files_all_1:
            one_file = one_file_raw.split("/")[-1]
            file_all_1.append(one_file)
        print(len(file_all_1))

    ### move final merged files ###
    mv_cmd = "mv " + in_dir_1 + "*" + " " + out_dir 
    Popen(mv_cmd,shell=True).wait()
    rm_cmd = "rm " + out_dir + "temp.fastq"
    Popen(rm_cmd,shell=True).wait()

    print("Merging Done~")






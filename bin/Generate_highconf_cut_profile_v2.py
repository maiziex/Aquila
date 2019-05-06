#!/usr/bin/env python
import pdb
#pdb.set_trace()
import os
import numpy as np
import subprocess
from argparse import ArgumentParser
from multiprocessing import Pool,cpu_count,active_children
import time
import pickle
from collections import defaultdict
import sys
from subprocess import Popen
script_path = os.path.dirname(os.path.abspath( __file__ ))

code_path = script_path + "/" 

parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--bam_file','-bam',help="input sorted bam file")
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--mbq','-bq',type=int,help="minimum base quality to consider a base for haplotype fragment, default 13", default=13)
parser.add_argument('--mmq','-mq',type=int,help="minimum read mapping quality to consider a read for phasing, default 20", default=20)
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs")
parser.add_argument('--uniq_map_dir','-uniq_dir', help="Directory to 100-mer uniqness")
args = parser.parse_args()


def get_read_depth(bam_file,mbq,mmq,chr_num,out_dir,xin):
    if chr_num == 23:
        chr_num = "X"
    else:
        chr_num = int(chr_num)
    noclips_bam = out_dir + "chr" + str(chr_num) + "_noclips.bam"
    depth_file_noHclip = out_dir + "chr" + str(chr_num) + "_depth_noHclip.txt"
    depth_file = out_dir  + "chr" + str(chr_num) + "_depth.txt"
    try:
        rm_hardclip_cmd = "samtools view -h " + bam_file + " chr" + str(chr_num) + " | awk '$6 !~ /H|S/{print}' | samtools view -bS - > " + noclips_bam
        depth_cmd = "samtools depth " + noclips_bam + " -q " + str(mbq_threshold) +  " -Q " + str(mmq_threshold) + " > " + depth_file_noHclip 
        depth_cmd_2 = "samtools depth " + bam_file + " -q " + str(mbq_threshold)  + " -Q " +str(mmq_threshold) + " -r chr" + str(chr_num) + " > " + depth_file 
    except:
        rm_hardclip_cmd = code_path + "samtools/" + "samtools view -h " + bam_file + " chr" + str(chr_num) + " | awk '$6 !~ /H|S/{print}' | " + code_path + "samtools/" + "samtools view -bS - > " + noclips_bam
        depth_cmd =  code_path + "samtools/" + "samtools depth " + noclips_bam + " -q " + str(mbq_threshold) +  " -Q " + str(mmq_threshold) + " > " + depth_file_noHclip 
        depth_cmd_2 = code_path + "samtools/" + "samtools depth " + bam_file + " -q " + str(mbq_threshold)  + " -Q " +str(mmq_threshold) + " -r chr" + str(chr_num) + " > " + depth_file 

    print(rm_hardclip_cmd)
    print(depth_cmd)
    print(depth_cmd_2)
    subprocess.Popen(rm_hardclip_cmd,shell=True).wait()
    subprocess.Popen(depth_cmd,shell=True).wait()
    subprocess.Popen(depth_cmd_2,shell=True).wait()



def get_coverage_per_pos(depth_file,output_file,xin):
    f = open(depth_file,"r")
    coverage_pos = defaultdict(int)
    for line in f:
        data = line.rsplit()
        pos_ = int(data[1]) + 1 # use "1" coordinates
        depth = int(data[2])
        coverage_pos[pos_] = depth
    print(len(coverage_pos))
    pickle.dump(coverage_pos,open(output_file,"wb"))


def get_global_track_for_breakpoints(cov_dict,uniq_map_dict,output_file,xin):
    global_track = defaultdict(list)
    cov_list = sorted(cov_dict.values())
    mean_cov = np.mean(cov_list)
    cov_threshold = mean_cov*0.8
    print(mean_cov, cov_threshold)
    for _pos, _cov in cov_dict.items():
        if _cov > cov_threshold:
            if _pos in uniq_map_dict:
                global_track[_pos] = [_cov]
    pickle.dump(global_track,open(output_file,"wb"))


def Run_all(bam_file,chr_start,chr_end,mbq_threshold,mmq_threshold,out_dir,uniq_map_dir):
    num_threads = 24
    round_times = int(np.ceil((chr_end - chr_start + 1)/num_threads))
    old_list = defaultdict(list)
    for ii in range(0,num_threads):
        for jj in range(chr_start + round_times*ii,chr_start + round_times*ii + round_times):
            if jj <= chr_end:
                old_list[ii+1].append(jj)

    new_list = defaultdict(list)
    times = 1
    for key,val in old_list.items():
        if times%2 == 1:
            new_list[key] = val
        elif times%2 == 0:
            new_list[key] = val[::-1] 
        times += 1
    
    run_list = defaultdict(list)
    for key,val in new_list.items():
        for ii in range(round_times):
            if ii < len(val):
                run_list[ii+1].append(val[ii])

    for key, val in run_list.items():
        pool = Pool(processes=num_threads)
        for chr_num in val:
            print("processing chr" + str(chr_num))
            pool.apply_async(get_read_depth,(bam_file,mbq_threshold,mmq_threshold,chr_num,out_dir,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()
        print("done")

    ###########
    num_threads = 12
    round_times = int(np.ceil((chr_end - chr_start + 1)/num_threads))
    old_list = defaultdict(list)
    for ii in range(0,num_threads):
        for jj in range(chr_start + round_times*ii,chr_start + round_times*ii + round_times):
            if jj <= chr_end:
                old_list[ii+1].append(jj)

    new_list = defaultdict(list)
    times = 1
    for key,val in old_list.items():
        if times%2 == 1:
            new_list[key] = val
        elif times%2 == 0:
            new_list[key] = val[::-1] 
        times += 1
    
    run_list = defaultdict(list)
    for key,val in new_list.items():
        for ii in range(round_times):
            if ii < len(val):
                run_list[ii+1].append(val[ii])

    for key, val in run_list.items():
        pool = Pool(processes=num_threads)
        for chr_num in val:
            if chr_num == 23:
                depth_file = out_dir + "chrX_depth_noHclip.txt"
            else:
                depth_file = out_dir + "chr" + str(chr_num) + "_depth_noHclip.txt"
            output_file = out_dir + "chr" + str(chr_num) + "_pos_cf.p"
            pool.apply_async(get_coverage_per_pos,(depth_file,output_file,"xin"))

        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()
        print("done here")
    print("Extract all depth file, Done~")

   ###################
    for chr_num in range(chr_start,chr_end + 1):
        cov_file = out_dir + "chr" + str(chr_num) + "_pos_cf.p" 
        cov_dict = pickle.load(open(cov_file,"rb"))
        uniq_map_file = uniq_map_dir + "uniq_map_chr" + str(chr_num) + ".p"
        uniq_map_dict = pickle.load(open(uniq_map_file,"rb"))
        output_file = out_dir + "chr" + str(chr_num) + "_global_track.p"
        get_global_track_for_breakpoints(cov_dict,uniq_map_dict,output_file,"xin")
    print("Extract all depth file, Done~")


    print("Finishing generating highconf profile, all done~")

   


if __name__ == "__main__":
    out_dir = args.out_dir
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    bam_file = args.bam_file
    chr_start = args.chr_start
    chr_end = args.chr_end
    mbq_threshold = args.mbq
    mmq_threshold = args.mmq
    out_dir = args.out_dir
    uniq_map_dir = args.uniq_map_dir
    Run_all(bam_file,chr_start,chr_end,mbq_threshold,mmq_threshold,out_dir,uniq_map_dir)
    #Popen("rm " + out_dir + "*.bam",shell=True).wait()
    #Popen("rm " + out_dir + "*.txt",shell=True).wait()
    #Popen("rm " + out_dir + "*_pos_cf.p",shell=True).wait()

#!/usr/bin/env python
import pdb
#pdb.set_trace()
from collections import defaultdict
import numpy as np
import pickle
from argparse import ArgumentParser
import os
from multiprocessing import Pool,cpu_count,active_children
import time
from subprocess import Popen
parser = ArgumentParser(description="Cut phase block for local assembly:")
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end at", default=23)
parser.add_argument('--block_threshold','-t',type=int,help="phase block threshold")
parser.add_argument('--block_len_use','-l',type=int,help="phase block len threshold")
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs")
parser.add_argument('--in_dir','-i_dir', help="Directory for inputs")
parser.add_argument('--highconf_profile_dir','-hc_dir', help="Directory for highconf profile")
args = parser.parse_args()


def Finalize_phase_block(output_file_raw, output_file,use_fake_PS_flag):
    PS_flag_dict = defaultdict(int)
    count_mole = defaultdict(int)
    f = open(output_file_raw, "r")
    fw = open(output_file, "w")
    for line in f:
        data = line.rsplit()
        var_dict = []
        if data[-1] == "hp1" or data[-1] == "hp2":
            PS_tag = data[-2]
            PS_flag = int(PS_tag.split(":")[1])
            HP_tag = data[-1]
            if use_fake_PS_flag[PS_flag] == 1:
                var_list = data[7:-2]
                count_ = 0
                for var in var_list:
                    var_pos = int(var.split(":")[0])
                    if count_ == 0:
                        min_pos = var_pos
                    else:
                        if min_pos > var_pos:
                            min_pos = var_pos
                    count_ += 1

                if count_mole[PS_flag] == 0:
                    PS_flag_dict[PS_flag] = min_pos
                else:
                    if min_pos < PS_flag_dict[PS_flag]:
                        PS_flag_dict[PS_flag] = min_pos

                count_mole[PS_flag] += 1

    f.close()
    f = open(output_file_raw, "r")
    for line in f:
        data = line.rsplit()
        if data[-1] == "hp1" or data[-1] == "hp2":
            PS_tag = data[-2]
            PS_flag = int(PS_tag.split(":")[1])
            if use_fake_PS_flag[PS_flag] == 1:
                HP_tag = data[-1]
                data[-2] = "PS:" + str(PS_flag_dict[PS_flag])
                fw.writelines("\t".join(data) + "\n")
            else:
                fw.writelines(line)
        else:
            fw.writelines(line)

    f.close()
    fw.close()


def Check_mole_pos_in_stop_end_list(mole_start,mole_end,stop_end_list):
    count = 0
    use_flags = []
    for stop_end in stop_end_list:
        step = stop_end
        if count > 0:
            if mole_start <= prev_step:
                if mole_end >= prev_step:
                    use_flags.append(count - 1)
            else:
                if mole_start <= step:
                    use_flags.append(count - 1)

        prev_step = step
        count += 1
    return use_flags


def Cut_phaseblock_for_phased_h5(h5_phased_file,chr_num,out_file,block_len_use,block_threshold,output_dir,bed_file,phase_block_file,global_track,HC_breakpoint_file,xin):
    HC_breakpoint = defaultdict(list)
    fw_HC_bk = open(output_dir + "chr" + str(chr_num) + "_HC_breakpoint.bed","w")
    fw_bed_file = open(bed_file,"w")
    count_mole_max = defaultdict(int)
    PS_flag_dict_max = defaultdict(int)
    f = open(h5_phased_file,"r")
    curr = 0
    for line in f:
        curr += 1
        data = line.rsplit()
        mole_len = int(data[3])
        var_dict = []
        if data[-1] == "hp1" or data[-1] == "hp2":
            PS_tag = data[-2]
            PS_flag = int(PS_tag.split(":")[1])
            HP_tag = data[-1]
            var_list = data[7:-2]
            count_ = 0
            for var in var_list:
                var_pos = int(var.split(":")[0])
                if count_ == 0:
                    max_pos = var_pos
                else:
                    if max_pos < var_pos:
                        max_pos = var_pos
                count_ += 1

            if count_mole_max[PS_flag] == 0:
                PS_flag_dict_max[PS_flag] = max_pos
            else:
                if max_pos > PS_flag_dict_max[PS_flag]:
                    PS_flag_dict_max[PS_flag] = max_pos

            count_mole_max[PS_flag] += 1

    f.close()
    all_len = []
    cut_phase_block = defaultdict(int)
    for key, val in PS_flag_dict_max.items():
        len_block = val - key
        all_len.append(len_block)
        fw_bed_file.writelines("chr" + str(chr_num) + "\t" + str(key)+ "\t" + str(val) + "\n")
        if len_block >= block_threshold:
            cut_phase_block[key] = val
      
    #pickle.dump(cut_phase_block,open("cut_phase_block.p","wb"))
    #pickle.dump(PS_flag_dict_max,open("PS_flag_dict_max.p","wb"))
    #cut_phase_block = pickle.load(open("cut_phase_block.p","rb"))
    #PS_flag_dict_max = pickle.load(open("PS_flag_dict_max.p","rb"))

    cut_phase_block_2 = defaultdict(list)
    cut_phase_block_3 = defaultdict(list)
    """  use 100kb for each phase block  """
    global_track_dict =  pickle.load(open(global_track,"rb"))
    count_break_highconf = 0
    count_break_highconf_2 = 0
    count_break_lowconf = 0
    all_cut_block_len = []
    for key, val in cut_phase_block.items():
        cut_phase_block_2[key].append(key)
        phase_block_len = val - key
        num_of_blocks_cut = int(phase_block_len/block_len_use)
        for ii in range(1,num_of_blocks_cut):
            cut_phase_block_2[key].append(key + block_len_use*ii)
        cut_phase_block_2[key].append(val)

        break_point_list = cut_phase_block_2[key]
        cut_phase_block_3[key].append(break_point_list[0])
        for ii in range(1,len(break_point_list)-1):
            break_point = break_point_list[ii]
            break_flag = 0
            if break_point in global_track_dict:
                HC_breakpoint[(key,val)].append(break_point)
                fw_HC_bk.writelines("chr" + str(chr_num) + "\t" + str(break_point) + "\t" + str(break_point+1) + "\n")
                cut_phase_block_3[key].append(break_point)
                count_break_highconf_2 += 1
            else:
                break_point_before = break_point - 5000
                break_point_after = break_point + 5000
                for jj in range(break_point_before, break_point_after+1):
                    if jj in global_track_dict:
                        HC_breakpoint[(key,val)].append(jj)
                        fw_HC_bk.writelines("chr" + str(chr_num) + "\t" + str(jj) + "\t" + str(jj+1) + "\n")
                        cut_phase_block_3[key].append(jj)
                        break_flag = 1
                        count_break_highconf += 1
                        all_cut_block_len.append(jj - cut_phase_block_3[key][-2])
                        break
                if break_flag == 0:
                    cut_phase_block_3[key].append(break_point)
                    count_break_lowconf += 1
                    all_cut_block_len.append(break_point - cut_phase_block_3[key][-2])

        cut_phase_block_3[key].append(break_point_list[-1])                   

    print(np.max(all_cut_block_len))
    print(np.min(all_cut_block_len))
    print(np.mean(all_cut_block_len))
    print(count_break_highconf_2,count_break_highconf,count_break_lowconf)
    count_test_1 = 0
    count_test_2 = 0
    for key, val_list in cut_phase_block_3.items():
        for jj in range(len(val_list)-1):
            val_1 = val_list[jj]
            val_2 = val_list[jj+1]
            count_test_1 += 1
            fw_bed_file.writelines("chr" + str(chr_num) + "\t" + str(val_1)+ "\t" + str(val_2) + "\n")

    for key, val_list in cut_phase_block_2.items():
        for jj in range(len(val_list)-1):
            val_1 = val_list[jj]
            val_2 = val_list[jj+1]
            count_test_2 += 1
            
    fw_bed_file.close()
    temp_file = output_dir + "temp_file_chr" + str(chr_num)
    fw = open(temp_file,"w")
    f = open(h5_phased_file,"r")
    define_PS_flag = defaultdict(lambda: defaultdict(int))
    final_PS_block_cut = defaultdict(int)
    for line in f:
        data = line.rsplit()
        if data[-1] == "hp1" or data[-1] == "hp2":
            PS_tag = data[-2]
            PS_flag = int(PS_tag.split(":")[1])
            HP_tag = data[-1]
            if PS_flag in cut_phase_block_3.keys(): # if the phase block length >= threshold_length
                stop_end_list = cut_phase_block_3[PS_flag]
                mole_start = int(data[1])
                mole_end = int(data[2])
                use_flags = Check_mole_pos_in_stop_end_list(mole_start,mole_end,stop_end_list)
                if use_flags == []:
                    use_flags = [0]
                for _pos in use_flags:
                    if _pos == 0:  # keep the first block PS flag
                        fw.writelines(line)
                        final_PS_block_cut[PS_flag] = stop_end_list[_pos+1]
                   
                    else:  # reassign the fake phase block number for the rest
                        cur_PS_flag = stop_end_list[_pos]
                        data[-2] = "PS:" + str(cur_PS_flag)
                        final_PS_block_cut[cur_PS_flag] = stop_end_list[_pos+1]
                        fw.writelines("\t".join(data) + "\n")
            else:
                fw.writelines(line)
                final_PS_block_cut[PS_flag] = PS_flag_dict_max[PS_flag]
        else:
            fw.writelines(line)

    fw.close()
    f.close()
    fw_HC_bk.close()
    if os.path.exists(out_file):
        Popen("rm " + out_file,shell=True).wait()
    Popen("mv " + temp_file + " " + out_file,shell=True).wait()
    pickle.dump(final_PS_block_cut,open(phase_block_file,"wb"))
    pickle.dump(HC_breakpoint,open(HC_breakpoint_file,"wb"))
    print("done~")    





if __name__ == "__main__":
    input_dir = args.in_dir
    hc_dir = args.highconf_profile_dir
    output_dir = args.out_dir
    if os.path.exists(output_dir):
        print("using existing output folder: " + output_dir)
    else:
        os.makedirs(output_dir)
    chr_start = args.chr_start
    chr_end = args.chr_end
    block_len_use = args.block_len_use
    block_threshold = args.block_threshold
    #pool = Pool(chr_end - chr_start + 1)
    for chr_num in range(chr_start,chr_end+1):
        file_name = input_dir + "chr" + str(chr_num) + ".phased_final"
        out_file = output_dir + "chr" + str(chr_num) + ".phased_final_cut_by_" + str(block_len_use) 
        bed_file = output_dir + "chr" + str(chr_num) + ".phased_final_cut_by_" + str(block_len_use) + "_phase_blocks.bed"
        phase_block_file = output_dir + "chr" + str(chr_num) + ".phased_final_cut_by_" + str(block_len_use) + "_phase_blocks.p"
        HC_breakpoint_file = output_dir + "chr" + str(chr_num) + ".phased_final_cut_by_" + str(block_len_use) + "_HC_breakpoint_2.p"
        global_track = hc_dir + "chr" + str(chr_num) + "_global_track.p"
        Cut_phaseblock_for_phased_h5(file_name,chr_num,out_file,block_len_use,block_threshold,output_dir,bed_file,phase_block_file,global_track,HC_breakpoint_file,"xin")
        #pool.apply_async(Cut_phaseblock_for_phased_h5,(file_name,chr_num,out_file,block_len_use,block_threshold,output_dir,bed_file,phase_block_file,global_track,HC_breakpoint_file,"xin"))
   
    """
    pool.close()
    while len(active_children()) > 1:
        time.sleep(0.5)
    pool.join()
    """
    print("all done~")



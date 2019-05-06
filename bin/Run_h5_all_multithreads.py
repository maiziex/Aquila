#!/usr/bin/env python
import pdb
#pdb.set_trace()
import glob
import os
import numpy as np
from argparse import ArgumentParser
from collections import defaultdict
import time
import multiprocessing as mp
import pickle
import sys
import pysam
from subprocess import Popen
from multiprocessing import Pool,cpu_count,active_children
parser = ArgumentParser(description="Extract reads for all small chunks:")
parser.add_argument('--bam_file','-bam',help="bam file")
parser.add_argument('--vcf_file','-v',help="vcf file")
parser.add_argument('--mbq_threshold','-bq',type=int,help="mbq threshold", default=13)
parser.add_argument('--mmq_threshold','-mq',type=int,help="mmq threshold", default=20)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--boundary','-b',type=int,help="cut boundary", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--num_threads','-n',type=int,help="number of threads", default=8)
parser.add_argument('--sample_name','-s',help="sample name")
parser.add_argument('--out_dir','-o',help="output folder")

args = parser.parse_args()
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 


def flatten(listoflist):
    flat_list = []
    for alist in listoflist:
        for val in alist:
            flat_list.append(val)
    return flat_list


def Cal_snp_ratio_vs_depth(vcf_file,mbq_threshold):
    f = open(vcf_file,"r")
    curr = 0
    variant_dp_dict =  defaultdict(int)
    curr = 0
    count_total = 0
    total_depth = []
    for line in f:
        data = line.rsplit()
        curr += 1
        if data[0][:2] != "##" and data[0][:2] != "#C":
            try:
                qual = float(data[5])
                chr_num = int(data[0][3:])
                pos = int(data[1]) - 1
                ref = data[3]
                alt = data[4]
                GT = data[9].split(":")[0]
                if (GT == "0/1"  or GT == "1/0" or GT == "0|1" or GT == "1|0") and len(ref) == 1 and len(alt) == 1 and qual >= mbq_threshold:
                    break_flag = 0
                    _format = data[8].split(":")
                    AO_idx = _format.index("AO")
                    RO_idx = _format.index("RO")
                    DP_idx = _format.index("DP")
                    _format_info = data[9].split(":")
                    ao_depth = float(_format_info[AO_idx])
                    ro_depth = float(_format_info[RO_idx])
                    dp_depth = float(_format_info[DP_idx])
                    if ao_depth >= ro_depth:
                        ratio = ro_depth/ao_depth
                    else:
                        ratio = ao_depth/ro_depth
                    for ii in range(0,11):
                        if dp_depth <= 10*(ii + 1) and dp_depth > 10*ii: 
                            for jj in range(0,20):
                                if ratio <= 0.05*(jj+1) and ratio > 0.05*jj:
                                    variant_dp_dict[(ii,jj)] += 1
                                    count_total += 1
                                    break_flag = 1
                                    total_depth.append(dp_depth)
                                    break
                    if break_flag == 0:
                        if dp_depth > 110:
                            for jj in range(0,20):
                                if ratio <= 0.05*(jj+1) and ratio > 0.05*jj:
                                    variant_dp_dict[(11,jj)] += 1
                                    count_total += 1
                                    total_depth.append(dp_depth)
            except:
                if data[0] == "chrX":
                    chr_num = "X"
                qual = float(data[5])
                pos = int(data[1]) - 1
                ref = data[3]
                alt = data[4]
                GT = data[9].split(":")[0]
                if (GT == "0/1"  or GT == "1/0" or GT == "0|1" or GT == "1|0") and len(ref) == 1 and len(alt) == 1 and qual >= mbq_threshold:
                    break_flag = 0
                    _format = data[8].split(":")
                    AO_idx = _format.index("AO")
                    RO_idx = _format.index("RO")
                    DP_idx = _format.index("DP")
                    _format_info = data[9].split(":")
                    ao_depth = float(_format_info[AO_idx])
                    ro_depth = float(_format_info[RO_idx])
                    dp_depth = float(_format_info[DP_idx])
                    if ao_depth >= ro_depth:
                        ratio = ro_depth/ao_depth
                    else:
                        ratio = ao_depth/ro_depth
                    for ii in range(0,11):
                        if dp_depth <= 10*(ii + 1) and dp_depth > 10*ii: 
                            for jj in range(0,20):
                                if ratio <= 0.05*(jj+1) and ratio > 0.05*jj:
                                    variant_dp_dict[(ii,jj)] += 1
                                    count_total += 1
                                    break_flag = 1
                                    total_depth.append(dp_depth)
                                    break
                    if break_flag == 0:
                        if dp_depth > 110:
                            for jj in range(0,20):
                                if ratio <= 0.05*(jj+1) and ratio > 0.05*jj:
                                    variant_dp_dict[(11,jj)] += 1
                                    count_total += 1
                                    total_depth.append(dp_depth)


    #pickle.dump(variant_dp_dict,open("variant_dp_dict.p","wb"))
    avg_depth = np.mean(total_depth)
    median_depth = np.median(total_depth)
    print(count_total)
    print(np.mean(total_depth))
    print(np.median(total_depth))
    print(np.max(total_depth))
    print(np.min(total_depth))
    return (avg_depth,median_depth)


def save_variant_dict(vcf_file,mbq_threshold,avg_depth,out_dir):
    f = open(vcf_file,"r")
    variant_dict =  defaultdict(list)
    for line in f:
        data = line.rsplit()
        if data[0][:2] != "##" and data[0][:2] != "#C":
            try:
                qual = float(data[5])
                chr_num = int(data[0][3:])
                pos = int(data[1]) - 1
                ref = data[3]
                alt = data[4]
                GT = data[9].split(":")[0]
                if (GT == "0/1"  or GT == "1/0") and len(ref) == 1 and len(alt) == 1 and qual >= mbq_threshold:
                    _format = data[8].split(":")
                    AO_idx = _format.index("AO")
                    RO_idx = _format.index("RO")
                    DP_idx = _format.index("DP")
                    _format_info = data[9].split(":")
                    ao_depth = float(_format_info[AO_idx])
                    ro_depth = float(_format_info[RO_idx])
                    dp_depth = float(_format_info[DP_idx])
                    if ao_depth >= ro_depth:
                        ratio = ro_depth/ao_depth
                    else:
                        ratio = ao_depth/ro_depth
                    if ratio >= 0.25 and dp_depth >= 0.1*avg_depth and dp_depth <= 0.9* avg_depth:
                        variant_dict[(chr_num,pos)] = [ref,alt,GT]
            except:
                if data[0] == "chrX":
                    chr_num = "X"
                qual = float(data[5])
                pos = int(data[1]) - 1
                ref = data[3]
                alt = data[4]
                GT = data[9].split(":")[0]
                if (GT == "0/1"  or GT == "1/0") and len(ref) == 1 and len(alt) == 1 and qual >= mbq_threshold:
                    _format = data[8].split(":")
                    AO_idx = _format.index("AO")
                    RO_idx = _format.index("RO")
                    DP_idx = _format.index("DP")
                    _format_info = data[9].split(":")
                    ao_depth = float(_format_info[AO_idx])
                    ro_depth = float(_format_info[RO_idx])
                    dp_depth = float(_format_info[DP_idx])
                    if ao_depth >= ro_depth:
                        ratio = ro_depth/ao_depth
                    else:
                        ratio = ao_depth/ro_depth
                    if ratio >= 0.25 and dp_depth >= 0.1*avg_depth and dp_depth <= 0.9* avg_depth:
                        variant_dict[(chr_num,pos)] = [ref,alt,GT]
       
    pickle.dump(variant_dict, open(out_dir + "variant_dict_heterozygous.p", "wb"))
    return variant_dict     


def get_match_num(cigar):
    used_num = ""
    cumu_num = 0
    for i in range(len(cigar)):
        cur_num = cigar[i]
        try:
            first_num = int(cur_num)
            used_num += cur_num
        except ValueError:
            if cur_num != "M":
                cumu_num += int(used_num)
                used_num = ""
            else:
                match_num = int(used_num)
                break
    return (cumu_num,match_num)


def get_match_num_revised(cigar):
    cigar_len = len(cigar)
    cigar_list = []
    num_string = ""
    for i in range(cigar_len):
        letter = cigar[i]
        if letter.isalpha():
            cigar_list.append(num_string)
            num_string = ""
            cigar_list.append(letter)
        else:
            num_string += letter

    indices = [i for i, x in enumerate(cigar_list) if x == "M"]
    cumu_num = 0
    cumu_start = 0
    cumu_num_list = []
    cumu_start_list = []
    match_num_list = []
    for idx in indices:
        match_num = int(cigar_list[idx - 1])
        match_num_list.append(match_num)
        for i in range(0,idx,1):
            parameter = cigar_list[i+1]
            if parameter in ["M", "S", "H", "I"]:
                cumu_num += int(cigar_list[i])
            if parameter in ["M", "D"]:
                cumu_start += int(cigar_list[i])

        cumu_num_list.append(cumu_num-match_num)
        cumu_start_list.append(cumu_start-match_num)
        cumu_num = 0
        cumu_start = 0

    return (cumu_num_list,cumu_start_list,match_num_list)


def check_read_hp(hp,barcode):
    if len(set(hp)) == 1:
        hp_use = True
        hp_flag = hp[0]
    else:
        count0 = hp.count("0")
        count1 = hp.count("1")
        if count0 >= 4*count1:
            hp_use = True
            hp_flag = "0"
        elif count1 >= 4*count0:
            hp_use = True
            hp_flag = "1"
        else:
            hp_use = False
            hp_flag = "none"
    return (hp_use,hp_flag)


def get_mole_variant(chr_num,mole_dict_3,mole_dict_4,mole_dict_5,barcode):
    mole_variant = defaultdict(list)
    for qname,cigar_all in mole_dict_5[barcode].items():
        for ii in range(len(mole_dict_5[barcode][qname])):
            cigar = cigar_all[ii]
            if cigar != None:
                curr_string =  mole_dict_4[barcode][qname][ii]
                cumu_num, match_num = get_match_num(cigar)
                start_locus = int(mole_dict_3[barcode][qname][ii])
                match_string = curr_string[cumu_num:cumu_num+match_num]
                end_locus = len(match_string) - 1
                for jj in range(0,end_locus-1):
                    idx_j = jj + start_locus - 1
                    val = variant_dict[(chr_num,idx_j)]
                    if val != []:
                        ref = val[0]
                        alt = val[1]
                        GT = val[2]
                        match_str = match_string[jj:jj+len(ref)]
                        if match_str == ref:
                            mole_variant[idx_j].append("0")
                        elif match_str == alt:
                            mole_variant[idx_j].append("1")

    return mole_variant


def get_mole_variant_revised(chr_num,mole_dict_3,mole_dict_4,mole_dict_5,barcode):
    mole_variant = defaultdict(list)
    for qname,cigar_all in mole_dict_5[barcode].items():
        for ii in range(len(mole_dict_5[barcode][qname])):
            cigar = cigar_all[ii]
            if cigar != None:
                curr_string =  mole_dict_4[barcode][qname][ii]
                cumu_num_list,cumu_start_list, match_num_list = get_match_num_revised(cigar)
                for kk in range(len(cumu_num_list)):
                    cumu_num = cumu_num_list[kk]
                    cumu_start = cumu_start_list[kk]
                    match_num = match_num_list[kk]
                    start_locus = int(mole_dict_3[barcode][qname][ii]) + cumu_start
                    match_string = curr_string[cumu_num:cumu_num+match_num]
                    end_locus = len(match_string) - 1
                    for jj in range(0,end_locus-1):
                        idx_j = jj + start_locus - 1
                        val = variant_dict[(chr_num,idx_j)]
                        if val != []:
                            ref = val[0]
                            alt = val[1]
                            GT = val[2]
                            match_str = match_string[jj:jj+len(ref)]
                            if match_str == ref:
                                mole_variant[idx_j].append("0")
                                #if kk > 0:
                                    #print("use here", ref, alt, match_str,match_num, match_num_list[0])
                            elif match_str == alt:
                                mole_variant[idx_j].append("1")
                                #if kk > 0:
                                    #print("use here", ref, alt, match_str,match_num, match_num_list[0])
                            #else:
                                #if kk > 0:
                                    #print("wrong here",ref,alt,match_str)

    return mole_variant


def save_pickle_file(dict1,fname):
    for value in dict1:
        dict1[value] = dict(dict1[value])
    my_dict = dict(dict1)
    with open(fname, "wb") as f:
        pickle.dump(my_dict,f) 


def process_sorted_bam(bam_file,output_file,qname_file,qname_pos_file,variant_dict,threshold_start,threshold_end,boundary,mmq_threshold,xin):
    #qname_pos = defaultdict(int)
    qname_pos = defaultdict(list)
    mole_qname_dict = defaultdict(lambda: defaultdict(list))
    sam_file = pysam.AlignmentFile(bam_file, "rb")
    threshold_pairs = 2
    fw= open(output_file,"w")
    mole_dict = defaultdict(list)
    mole_dict_2 = defaultdict(lambda: defaultdict(list))
    mole_dict_3 = defaultdict(lambda: defaultdict(list))
    mole_dict_4 = defaultdict(lambda: defaultdict(list))   # read string
    mole_dict_5 = defaultdict(lambda: defaultdict(list))   # cigar
    curr = 0
    count_mole = 1
    chr_begin = threshold_start
    if threshold_start == 23:
        use_chr_num = "chrX"
    else:
        use_chr_num = "chr" + str(threshold_start)
    for read in sam_file.fetch(use_chr_num):
        curr += 1
        raw_chr_num = read.reference_name
        if use_chr_num == "chrX":
            chr_num = "X"
        else:
            chr_num = int(read.reference_name[3:])
        tags = read.get_tags()
        barcode_field = [s for s in tags if "BX" in s]
        mq = read.mapping_quality
        if barcode_field != [] and mq >= mmq_threshold:
            barcode =  barcode_field[0][1].split("-")[0]
            start_pos = read.pos + 1   # use "1" coordinate
            qname = read.qname
            if read.is_read1:
                read_pair = 1
            elif read.is_read2:
                read_pair = 2
            #qname_pos[(qname,read_pair)] = start_pos
            qname_pos[qname].append(start_pos)
            if len(mole_dict[barcode]) == 0:
                mole_dict[barcode].append(start_pos) 
                mole_dict_2[barcode][qname].append(read_pair)
                mole_dict_3[barcode][qname].append(start_pos) 
                mole_dict_4[barcode][qname].append(read.seq) 
                mole_dict_5[barcode][qname].append(read.cigarstring) 
            elif len(mole_dict[barcode]) > 0:
                dist = start_pos - mole_dict[barcode][-1]
                if dist < boundary:   
                    mole_dict[barcode].append(start_pos) 
                    mole_dict_2[barcode][qname].append(read_pair)
                    mole_dict_3[barcode][qname].append(start_pos)
                    mole_dict_4[barcode][qname].append(read.seq)
                    mole_dict_5[barcode][qname].append(read.cigarstring)
                else:
                    count_del = 0
                    for key,value in mole_dict_2[barcode].items():
                        if len(value) == 1:
                            count_del += 1

                    all_pos = mole_dict_3[barcode].values()
                    num_of_pairs = len(all_pos) - count_del

                    if num_of_pairs >= threshold_pairs :  # 2 pairs
                        poss = flatten(all_pos)
                        start_ = min(poss)
                        end_ = max(poss)
                        mole_len = end_ - start_  + 150
                        mole_qname_dict[count_mole] = mole_dict_2[barcode].copy()
                        mole_variant = get_mole_variant_revised(chr_num,mole_dict_3,mole_dict_4,mole_dict_5,barcode)
                        if mole_variant == {}:
                            fw.write(str(chr_num) + "\t" + str(start_) + "\t" + str(end_) + "\t" + str(mole_len) + "\t" + str(len(poss)) + "\t" + str(barcode) + "\t" + str(count_mole) + "\n")
                        else:
                            fw.write(str(chr_num) + "\t" + str(start_) + "\t" + str(end_) + "\t" + str(mole_len) + "\t" + str(len(poss)) + "\t" + str(barcode) + "\t" + str(count_mole) + "\t")
                            count_var = 0
                            for locus_start, hp in mole_variant.items():
                                hp_use, hp_flag = check_read_hp(hp,barcode)
                                if count_var == len(mole_variant) - 1:
                                    if hp_use:
                                        fw.write(str(locus_start) + ":" + hp_flag + "\n")
                                        #count_1 += 1
                                    else:
                                        fw.write("\n")
                                        #count_2 += 1
                                else:
                                    if hp_use:
                                        fw.write(str(locus_start) + ":" + hp_flag + "\t")
                                        #count_1 += 1
                                    #else:
                                        #count_2 += 1
                                count_var += 1

                        mole_dict[barcode] = [] 
                        mole_dict_2[barcode] = defaultdict(list)
                        mole_dict_3[barcode] = defaultdict(list) 
                        mole_dict_4[barcode] = defaultdict(list) 
                        mole_dict_5[barcode] = defaultdict(list) 
                        mole_dict[barcode].append(start_pos)
                        mole_dict_2[barcode][qname].append(read_pair)
                        mole_dict_3[barcode][qname].append(start_pos)
                        mole_dict_4[barcode][qname].append(read.seq)
                        mole_dict_5[barcode][qname].append(read.cigarstring)
                        count_mole += 1
                    else:
                        mole_dict[barcode] = []
                        mole_dict_2[barcode] = defaultdict(list)
                        mole_dict_3[barcode] = defaultdict(list)
                        mole_dict_4[barcode] = defaultdict(list)
                        mole_dict_5[barcode] = defaultdict(list)
                        mole_dict[barcode].append(start_pos)
                        mole_dict_2[barcode][qname].append(read_pair)
                        mole_dict_3[barcode][qname].append(start_pos)
                        mole_dict_4[barcode][qname].append(read.seq)
                        mole_dict_5[barcode][qname].append(read.cigarstring)


    if threshold_end == 23:
        chr_num = "X"
    else:
        chr_num = threshold_end

    for barcode,value in mole_dict.items():
        if len(value) > 0:
            count_del = 0
            for qname,num_reads in mole_dict_2[barcode].items():
                if len(num_reads) == 1:
                    count_del += 1

            all_pos = mole_dict_3[barcode].values()
            num_of_pairs = len(all_pos) - count_del
            if num_of_pairs >= threshold_pairs : # 2 pairs
                poss = flatten(all_pos)
                start_ = min(poss)
                end_ = max(poss)
                mole_len = end_ - start_  + 150
                mole_qname_dict[count_mole] = mole_dict_2[barcode].copy()
                mole_variant = get_mole_variant_revised(chr_num,mole_dict_3,mole_dict_4,mole_dict_5,barcode)
                if mole_variant == {}:
                    fw.write(str(chr_num) + "\t" + str(start_) + "\t" + str(end_) + "\t" + str(mole_len) + "\t" + str(len(poss)) + "\t" + str(barcode) + "\t" + str(count_mole) + "\n")
                else:
                    fw.write(str(chr_num) + "\t" + str(start_) + "\t" + str(end_) + "\t" + str(mole_len) + "\t" + str(len(poss)) + "\t" + str(barcode) + "\t" + str(count_mole) + "\t")
                    count_var = 0
                    for locus_start, hp in mole_variant.items():
                        hp_use, hp_flag = check_read_hp(hp,barcode)
                        if count_var == len(mole_variant) - 1:
                            if hp_use:
                                fw.write(str(locus_start) + ":" + hp[0] + "\n")
                                #count_1 += 1
                            else:
                                fw.write("\n")
                                #count_2 += 1
                        else:
                            if hp_use:
                                fw.write(str(locus_start) + ":" + hp[0] + "\t")
                                #count_1 += 1
                            #else:
                                #count_2 += 1
                        count_var += 1

                count_mole += 1

    
    fw.close()
    sam_file.close()
    save_pickle_file(mole_qname_dict,qname_file)
    pickle.dump(qname_pos,open(qname_pos_file,"wb"))
    mole_qname_dict.clear()
    qname_pos.clear()



def Run_h5_all(bam_file,vcf_file,chr_start,chr_end,sample_name,mbq_threshold,mmq_threshold,boundary,num_threads,out_dir,variant_dict):
    count = 1
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
        pool = mp.Pool(processes=num_threads)
        for chr_num in val:
            print("processing chr" + str(chr_num))
            h5_file = out_dir + sample_name + "_chr" + str(chr_num)
            h5_qname_file = out_dir + sample_name + "_chr" + str(chr_num) + "_qname.p"
            h5_qname_pos_file = out_dir + sample_name + "_chr" + str(chr_num) + "_qname_pos.p"
            pool.apply_async(process_sorted_bam,(bam_file,h5_file,h5_qname_file,h5_qname_pos_file,variant_dict,chr_num,chr_num,boundary,mmq_threshold,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()
        print("done here")

    # sort h5 file by coordinates
    for chr_num in range(chr_start,chr_end + 1):
        one_h5 = out_dir + sample_name + "_chr" + str(chr_num)
        one_h5_sorted = out_dir + sample_name + "_chr" + str(chr_num) + "_sorted.h5"
        sort_cmd = "cat " + one_h5 + " | sort -k2n > " + one_h5_sorted
        Popen(sort_cmd,shell=True).wait()
    print("All Done~")




if __name__ == "__main__":
    out_dir = args.out_dir
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    bam_file = args.bam_file
    vcf_file = args.vcf_file
    mbq_threshold = args.mbq_threshold
    mmq_threshold = args.mmq_threshold
    boundary = args.boundary
    sample_name = args.sample_name
    avg_depth,median_depth = Cal_snp_ratio_vs_depth(vcf_file,mbq_threshold)
    variant_dict = save_variant_dict(vcf_file,mbq_threshold,avg_depth,out_dir)
    var_depth_file = out_dir + "median_depth_for_var.txt"
    if ~os.path.isfile(var_depth_file):
        Popen("touch " + var_depth_file,shell=True).wait()
        f = open(var_depth_file,"w")
        f.writelines(str(median_depth) + "\n")
        f.close()
    #variant_dict = pickle.load(open("variant_dict_heterozygous.p","wb"))
    chr_start = args.chr_start
    chr_end = args.chr_end
    num_threads = args.num_threads
    Run_h5_all(bam_file,vcf_file,chr_start,chr_end,sample_name,mbq_threshold,mmq_threshold,boundary,num_threads,out_dir,variant_dict)

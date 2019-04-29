import pdb
#pdb.set_trace()
import pickle
from collections import defaultdict
from argparse import ArgumentParser
import os
import sys


def Check_contig_in_HCbk_dict(HCbk_dict,PS_flag_1,PS_flag_2):
    break_flag = 0
    for key,val in HCbk_dict.items():
        PS_flag_1_curr = key[0]
        PS_flag_2_curr = key[1]
        if PS_flag_1 >= PS_flag_1_curr and PS_flag_2 <= PS_flag_2_curr:
            use_PS_flag_1 = PS_flag_1_curr
            use_PS_flag_2 = PS_flag_2_curr
            break_flag = 1
            break
    if break_flag == 1:
        return (use_PS_flag_1,use_PS_flag_2)
    else:
        return (0,0)


def Make_supercontig_based_on_HCbk(HCbk_file,contig_sam_file,ref_fasta_dict,use_chr_num):
    merge_dict = defaultdict(lambda: defaultdict(list))
    merge_len = defaultdict(list)
    merge_contig_num = defaultdict(int)
    f = open(contig_sam_file,"rb")
    HCbk_dict = pickle.load(open(HCbk_file,"rb"))
    curr = 0
    for line in f:
        curr += 1
        data = line.decode().rsplit()
        if data[0][0] != "@":
            mapq = int(data[4])
            chr_num = data[2]
            if chr_num == "chrX":
                chr_num = "chr23"
            cigar = data[5]
            if "H" in cigar:
                H_clip_flag = 1
            else:
                H_clip_flag = 0
            if mapq >= 20 and chr_num == use_chr_num and H_clip_flag == 0:
                data_split = data[0].split("_")
                contig_num = int(data_split[0])
                PS_flag_1 = int(data_split[1].split("PS")[1].split(":")[0])
                PS_flag_2 = int(data_split[1].split("PS")[1].split(":")[1])
                HP_flag = data_split[2]
                locus = int(data[3])
                contig_seq = data[9]
                use_PS_flag_1,use_PS_flag_2 = Check_contig_in_HCbk_dict(HCbk_dict,PS_flag_1,PS_flag_2)
                if use_PS_flag_1 != 0 and use_PS_flag_2 != 0:
                    merge_contig_num[contig_num] = 1
                    merge_dict[(use_PS_flag_1,use_PS_flag_2,HP_flag)][contig_num] = [contig_num,chr_num,contig_seq, str(PS_flag_1) + ":" + str(PS_flag_2) + "_" + HP_flag]

    count_1 = 0
    count_2 = 0
    count_3 = 0
    count_supercontig = 1
    count_no_supercontig = 0
    supercontig_dict = defaultdict(list)
    remove_contig_dict = defaultdict(int)
    use_for_supercontig_dict = defaultdict(int)
    ### For every big phase block, starting merging
    for key, contig_info_dict in merge_dict.items():
        max_contig_num = sorted(contig_info_dict.keys())[-1]
        new_contig_num = max_contig_num
        HP_flag = key[2]
        HCbk_list = HCbk_dict[(key[0],key[1])]
        for HC_bk in HCbk_list[1:-1]:
            HC_bk_string = ref_fasta_dict[HC_bk:(HC_bk+20)]
            idx_all = []
            contig_merge = []
            for contig_num,one_contig in contig_info_dict.items():
                break_flag = 0
                contig_seq = one_contig[2]
                try:
                    idx = contig_seq.index(HC_bk_string)
                    idx_all.append(idx)
                    contig_merge.append(one_contig)
                except:
                    pass

            idx_dict = defaultdict(int)
            contig_num_dict = defaultdict(int)
            contig_seq_dict = defaultdict(int)
            contig_PS_dict = defaultdict(str)
            before_bk_len = defaultdict(int)
            after_bk_len = defaultdict(int)
            if len(contig_merge) == 2:
                count = 0
                print("----")
                for one_contig in contig_merge:
                    idx_dict[count] = idx_all[count]
                    contig_num_dict[count] = one_contig[0]
                    contig_seq_dict[count] = one_contig[2]
                    contig_PS_dict[count] = one_contig[-1]
                    before_bk_len[count] = len(contig_seq_dict[count][:idx_dict[count]])
                    after_bk_len[count] = len(contig_seq_dict[count][idx_dict[count]:])
                    #print(contig_num,len(contig_seq),idx,before_bk_len,after_bk_len)
                    #print(contig_seq[idx:idx+50])
                    count += 1 
                count_2 += 1

                ### Make Supercontig 
                if before_bk_len[0] < before_bk_len[1] and after_bk_len[0] > after_bk_len[1]:
                    before_bk_string = contig_seq_dict[1][:idx_dict[1]]
                    after_bk_string = contig_seq_dict[0][idx_dict[0]:]
                    supercontig = before_bk_string + after_bk_string
                    if contig_num_dict[0] <= max_contig_num:
                        use_for_supercontig_dict[contig_num_dict[0]] = 1
                    if contig_num_dict[1] <= max_contig_num:
                        use_for_supercontig_dict[contig_num_dict[1]] = 1
                    del contig_info_dict[contig_num_dict[0]]
                    del contig_info_dict[contig_num_dict[1]]
                    new_contig_num += 1
                    PS_str = contig_PS_dict[0] + "-" + contig_PS_dict[1]
                    contig_info_dict[new_contig_num] = [new_contig_num,chr_num,supercontig,PS_str]
                    print(len(contig_seq_dict[0]),len(contig_seq_dict[1]),len(supercontig))
                elif before_bk_len[1] < before_bk_len[0] and after_bk_len[1] > after_bk_len[0]:
                    before_bk_string = contig_seq_dict[0][:idx_dict[0]]
                    after_bk_string = contig_seq_dict[1][idx_dict[1]:]
                    supercontig = before_bk_string + after_bk_string
                    if contig_num_dict[0] <= max_contig_num:
                        use_for_supercontig_dict[contig_num_dict[0]] = 1
                    if contig_num_dict[1] <= max_contig_num:
                        use_for_supercontig_dict[contig_num_dict[1]] = 1
                    del contig_info_dict[contig_num_dict[0]]
                    del contig_info_dict[contig_num_dict[1]]
                    new_contig_num += 1
                    PS_str = contig_PS_dict[0] + "-" + contig_PS_dict[1]
                    contig_info_dict[new_contig_num] = [new_contig_num,chr_num,supercontig,PS_str]
                    #print(len(contig_seq_dict[0]),len(contig_seq_dict[1]),len(supercontig))
                elif before_bk_len[0] >= before_bk_len[1] and after_bk_len[0] >= after_bk_len[1]:
                    if contig_num_dict[1] <= max_contig_num:
                        remove_contig_dict[contig_num_dict[1]] = 1
                elif before_bk_len[1] >= before_bk_len[0] and after_bk_len[1] >= after_bk_len[0]:
                    if contig_num_dict[0] <= max_contig_num:
                        remove_contig_dict[contig_num_dict[0]] = 1
                else:
                    print("wrong")
                    count_no_supercontig += 1

            elif len(contig_merge) == 1:
                count_1 += 1
            elif len(contig_merge) > 2:
                count_3 += 1

        print("finishe one big phase block")
        for _contig_num,_contig_one in contig_info_dict.items():
            if _contig_num > max_contig_num:
                supercontig_dict[count_supercontig,key] = _contig_one
                count_supercontig += 1

    print("finished~")
    return (supercontig_dict,use_for_supercontig_dict,remove_contig_dict)


def Finalize_contig_fasta_file(contig_fasta_file,supercontig_dict,use_for_supercontig_dict,remove_contig_dict,output_file):
    fw = open(output_file,"w")
    f = open(contig_fasta_file,"r")
    count = 1
    remove_flag = 0
    count =  1
    count_contig =  1
    curr = 0
    for line in f:
        curr += 1
        data = line.rsplit()
        if count%2 == 1:
            contig_num = int(data[0].split("_")[0][1:])
            PS_flag_1 = int(data[0].split("_")[1].split("PS")[1].split(":")[0])
            PS_flag_2 = int(data[0].split("_")[1].split("PS")[1].split(":")[1])
            HP_flag = data[0].split("_")[2]
            if contig_num in use_for_supercontig_dict or contig_num in remove_contig_dict:
                remove_flag = 1
        elif count%2 == 0:
            if remove_flag == 1:
                remove_flag = 0
            else:
                fw.writelines(">" + str(count_contig) + "_PS" + str(PS_flag_1) + ":" + str(PS_flag_2) + "_" + HP_flag +  "\n")
                fw.writelines(line)
                count_contig += 1
                remove_flag = 0
        count += 1
    
    for key,val in supercontig_dict.items():
        contig_seq = val[2]
        PS_flag_1 = key[1][0]
        PS_flag_2 = key[1][1]
        HP_flag = key[1][2]
        stitch_info = val[-1]
        fw.writelines(">" + str(count_contig) + "_PS" + str(PS_flag_1) + ":" + str(PS_flag_2) + "_" + HP_flag + "_merge" + stitch_info + "\n")
        fw.writelines(contig_seq + "\n")
        count_contig += 1
 
    f.close()
    fw.close()

    print("done")


def Contig_start(chr_num,cut_threshold,ref_dir,out_dir,phase_cut_folder,xin):
    HCbk_file = phase_cut_folder + "chr" + str(chr_num) +  ".phased_final_cut_by_" + str(cut_threshold) + "_HC_breakpoint_2.p"
    contig_fasta_file = out_dir + "Aquila_cutPBHC_minicontig_chr" + str(chr_num) + ".fasta"
    contig_sam_file = out_dir + "Aquila_cutPBHC_minicontig_chr" + str(chr_num) + ".sam"
    ref_fasta_file = ref_dir + "ref_seq_chr" + str(chr_num) + ".p"
    ref_fasta_dict = pickle.load(open(ref_fasta_file,"rb"))
    output_file = out_dir + "Aquila_Contig_chr" + str(chr_num) + ".fasta"

    supercontig_dict,use_for_supercontig_dict,remove_contig_dict = Make_supercontig_based_on_HCbk(HCbk_file,contig_sam_file,ref_fasta_dict,"chr" + str(chr_num))
    Finalize_contig_fasta_file(contig_fasta_file,supercontig_dict,use_for_supercontig_dict,remove_contig_dict,output_file)


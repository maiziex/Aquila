import pdb
#pdb.set_trace()
from collections import defaultdict
import pickle
import math
import os
import numpy as np
from Assign_phase_block_v3 import * 
import sys
import operator
from scipy.special import comb
import glob
from subprocess import Popen



def split_h5(h5_file,chr_num,output_dir):
    f = open(h5_file,"r")
    file_num = 1
    filename = output_dir + "chr" + str(chr_num) + "_" +  str(file_num)
    fw = open(filename,"w")
    mole_end_max = 0
    count = 0
    curr = 0
    for line in f:
        #print(curr)
        curr += 1
        data = line.rsplit()
        if len(data) >= 9:
            mole_start = int(data[1])
            mole_end = int(data[2])
            if mole_start > mole_end_max and count > 0:
                fw.close()
                if count > 0 and count <= 10000 and file_num > 1:
                    prev_filename = output_dir + "chr" + str(chr_num) + "_" +  str(file_num-1)
                    temp_file = "chr_of_" + str(chr_num) + "_temp"
                    Popen('cat ' + prev_filename + ' ' + filename + ' > ' + output_dir  +  temp_file,shell=True).wait()
                    Popen('mv ' + output_dir + temp_file + ' ' + prev_filename,shell=True).wait()
                    count = 0
                    mole_end_max = mole_end
                    filename = output_dir + "chr" + str(chr_num) + "_" +  str(file_num)
                    fw = open(filename,"w")
                    fw.writelines(line)
                    count += 1

                else:
                    count = 0
                    mole_end_max = mole_end
                    file_num += 1
                    filename = output_dir + "chr" + str(chr_num) + "_" +  str(file_num)
                    fw = open(filename,"w")
                    fw.writelines(line)
                    count += 1

            else:
                if mole_end > mole_end_max:
                    mole_end_max = mole_end
                
                fw.writelines(line)
                count += 1

    print("Raw clustering finished...")
    #print(file_num)
    return file_num


def Get_h5_for_mole_with_one_hetero_var(h5_file,chr_num,output_dir):
    f = open(h5_file,"r")
    filename = output_dir + "chr" + str(chr_num) + "_one_var"
    fw = open(filename,"w")
    for line in f:
        data = line.rsplit()
        if len(data) == 8:
            fw.writelines(line)

    print("Extract molecule with one heterozygous variant finished...")


def Generate_hyplotype_for_one_cluster(merge_var,merge_var_complement):
    use_var_dict = defaultdict(int)
    one_cluster_var_dict = defaultdict(int)
    count_all = 0
    count_conflict = 0
    _times = 3
    for var_pos, all_hp in merge_var.items():
        all_hp_complement = merge_var_complement[var_pos]
        if all_hp_complement == []:
            if len(set(all_hp)) == 1 and len(all_hp) >= 2:
                one_cluster_var_dict[var_pos] = all_hp[0] + 1
                count_all += 1
            else:
                if len(all_hp) >= 2:
                    count1 = all_hp.count(1)
                    count0 = all_hp.count(0)
                    if count1 >= _times*count0:
                        one_cluster_var_dict[var_pos] = 2 # 1
                    elif count0 >= _times*count1:
                        one_cluster_var_dict[var_pos] = 1 # 0
                    else:
                        count_conflict += 1
                    count_all += 1
                else:
                    count_conflict += 1
        else:
            count_all += 1
            use_var_dict[var_pos] = 1
            count1 = all_hp.count(1)
            count0 = all_hp.count(0)
            if count1 > _times*count0:
                use_hp = 2
            elif count0 > _times*count1:
                use_hp = 1
            else:
                use_hp = "nan"

            count1_complement = all_hp_complement.count(1)
            count0_complement = all_hp_complement.count(0)
            if count1_complement >= _times*count0_complement:
                use_hp_complement = 2
            elif count0_complement >= _times*count1_complement:
                use_hp_complement = 1
            else:
                use_hp_complement = "nan"
            if use_hp == "nan" and use_hp_complement == "nan":
                count_conflict += 1
            else:
                if use_hp == use_hp_complement:
                    count_conflict += 1
                else:
                    if use_hp == "nan":
                        if count1 > count0 and use_hp_complement == 1:
                            one_cluster_var_dict[var_pos] = get_complement_hp(use_hp_complement)
                        elif count0 > count1 and use_hp_complement == 2:
                            one_cluster_var_dict[var_pos] = get_complement_hp(use_hp_complement)
                        else:
                            count_conflict += 1
                    elif use_hp_complement == "nan":
                        if count1_complement > count0_complement and use_hp == 1:
                            one_cluster_var_dict[var_pos] = use_hp
                        elif count0_complement > count1_complement and use_hp == 2:
                            one_cluster_var_dict[var_pos] = use_hp
                        else:
                            count_conflict += 1
                    else:
                        one_cluster_var_dict[var_pos] = use_hp

    # process the uniq variant from complement molecule
    for var_pos, all_hp in merge_var_complement.items():
        if use_var_dict[var_pos] != 1 and all_hp != []:
            if len(set(all_hp)) == 1 and len(all_hp) >= 2 :
                one_cluster_var_dict[var_pos] = get_complement_hp_2(all_hp[0]) + 1
                count_all += 1
            else:
                if len(all_hp) >= 2:
                    count1 = all_hp.count(1)
                    count0 = all_hp.count(0)
                    if count1 >= _times*count0:
                        one_cluster_var_dict[var_pos] = 1 # 0
                    elif count0 >= _times*count1:
                        one_cluster_var_dict[var_pos] = 2 # 1
                    else:
                        count_conflict += 1
                    count_all += 1
                else:
                    count_conflict += 1

    conflict_percent = float(count_conflict)/count_all
    return (one_cluster_var_dict,conflict_percent)


def Generate_hyplotype_for_all_cluster(mole_dict,all_merge_cluster_dict, all_merge_cluster_complement_dict):
    final_cluster = []
    count_total = 0
    conflict_list = []
    for num, one_cluster, in all_merge_cluster_dict.items():
        merge_var = defaultdict(list)
        merge_var_complement = defaultdict(list)
        one_cluster_complement = all_merge_cluster_complement_dict[num]
        count_total += len(one_cluster)
        for mole_num in one_cluster:
            mole = mole_dict[mole_num]
            for var_pos,var_hp in mole.items():
                merge_var[var_pos].append(var_hp)
        for mole_num in one_cluster_complement:
            mole = mole_dict[mole_num]
            for var_pos,var_hp in mole.items():
                merge_var_complement[var_pos].append(var_hp)
        
        #print("------------")
        #print(merge_var)
        #print(merge_var_complement)
        #print("------------")

        one_cluster_var_dict,conflict_percent = Generate_hyplotype_for_one_cluster(merge_var,merge_var_complement)
        conflict_list.append(conflict_percent)
        final_cluster.append(one_cluster_var_dict)
        merge_var = defaultdict(list)

    #print(count_total)
    #print(len(all_merge_cluster_dict))
    #print("average conflict percent: " + str(np.mean(conflict_list)))
    return final_cluster

 
def get_complement_hp(hp):
    if hp == 1:
        return 2
    elif hp == 2:
        return 1


def get_complement_hp_2(hp):
    if hp == 1:
        return 0
    elif hp == 0:
        return 1


def nCr(n,r):
    f = math.factorial
    return float(f(n)) / (f(r) * f(n-r))


def Recursive_Clustering_for_Overlap_Variants(final_cluster,step,h5_filename,overlap_var_threshold):
    cluster_list = []
    curr = 0
    count_nonmerge = 0
    for idx in range(len(final_cluster)):
        var_dict = final_cluster[idx]
        #print(curr)
        curr += 1
        merge_cluster = 0
        if cluster_list != []:
            for num in range(len(cluster_list)):
                one_cluster = cluster_list[num]
                overlap_num = 0
                break_flag = 0
                for locus,hp in var_dict.items():
                    for locus_2,hp_2 in one_cluster.items():
                        if locus == locus_2:
                            overlap_num += 1
                            if overlap_num >= overlap_var_threshold:
                                break_flag = 1
                                break
                    if break_flag == 1:
                        break

                if overlap_num >= overlap_var_threshold:
                    uniq_var = []
                    conflict_var = []
                    nonconflict_var = []
                    for locus,hp in var_dict.items():
                        var_hp_1 = one_cluster[locus]
                        var_hp_2 = var_dict[locus]
                        if (var_hp_1 == 0 and var_hp_2 == 1 ) or (var_hp_1 == 0 and var_hp_2 == 2):
                            uniq_var.append(locus)
                        elif (var_hp_1 == 1 and var_hp_2 == 2) or (var_hp_1 == 2 and var_hp_2 == 1):
                            conflict_var.append(locus)
                        elif (var_hp_1 == 1 and var_hp_2 == 1) or (var_hp_1 == 2 and var_hp_2 == 2):
                            nonconflict_var.append(locus)

                    if len(conflict_var) >= len(nonconflict_var)*5.0:   # default = 3
                        for locus in nonconflict_var:
                            one_cluster.pop(locus)
                        for locus in uniq_var:
                            one_cluster[locus] = get_complement_hp(var_dict[locus])

                    elif len(nonconflict_var) >= len(conflict_var)*5.0:
                        for locus in conflict_var:
                            one_cluster.pop(locus)
                        for locus in uniq_var:
                            one_cluster[locus] = var_dict[locus]

                    else:
                        count_nonmerge += 1
                        merge_cluster = 0
                        continue

                    merge_cluster = 1
                    cluster_list[num] = one_cluster     # update, add 1-15-2018
                    break

        if merge_cluster == 0:
            cluster_list.append(var_dict)


    prev_cluster_num = len(final_cluster)

    if len(cluster_list) == prev_cluster_num:
        print("xin here: ")
        print("Clustering " +  str(step) + ":  " + str(count_nonmerge))
        print(prev_cluster_num,len(cluster_list), prev_cluster_num - len(cluster_list))
        print("Converged...")
        save_file = h5_filename + ".p"
        pickle.dump(cluster_list,open(save_file, "wb"))

        return 0

    step += 1
    if step == 1:
        print("Clustering " +  str(step)+ ":  " + str(count_nonmerge) + "," + str(round(float(count_nonmerge/prev_cluster_num),3)) )
    else:
        print("Clustering " +  str(step)+ ":  " + str(count_nonmerge)  )
    #print(prev_cluster_num,len(cluster_list), prev_cluster_num - len(cluster_list))
    final_cluster = cluster_list
    Recursive_Clustering_for_Overlap_Variants(final_cluster,step,h5_filename,overlap_var_threshold)
                    

def Clustering_for_phaseblock(mole_dict,all_merge_cluster_dict, all_merge_cluster_complement_dict,h5_filename,overlap_var_threshold):
    final_cluster = Generate_hyplotype_for_all_cluster(mole_dict,all_merge_cluster_dict, all_merge_cluster_complement_dict)
    Recursive_Clustering_for_Overlap_Variants(final_cluster,0,h5_filename,overlap_var_threshold)


def Recursive_clustering(all_merge_cluster_dict,all_merge_cluster_complement_dict,all_merge_cluster_dict_start_end, all_merge_cluster_complement_dict_start_end,step,mole_dict,mole_info,h5_filename,overlap_var_threshold,support_threshold):
    support_threshold = int(support_threshold)
    all_merge_cluster_dict_sorted = Sort_All_Clusters_by_Start(all_merge_cluster_dict,all_merge_cluster_dict_start_end)
    all_merge_cluster_dict = {}
    all_merge_cluster_dict = all_merge_cluster_dict_sorted
    cluster_dict = defaultdict(list)
    cluster_complement_dict = defaultdict(list)
    cluster_dict_start_end = defaultdict(list)
    cluster_complement_dict_start_end = defaultdict(list)
    curr = 0
    cluster_num = 0
    total_num = len(all_merge_cluster_dict)
    for _num,mole_num_list in all_merge_cluster_dict.items():
       # print(curr)
        curr += 1
        merge_cluster_flag = 0
        cluster = mole_num_list
        cluster_start = all_merge_cluster_dict_start_end[_num][0]
        cluster_end = all_merge_cluster_dict_start_end[_num][1]

        cluster_complement = all_merge_cluster_complement_dict[_num] 
        cluster_complement_start = all_merge_cluster_complement_dict_start_end[_num][0]
        cluster_complement_end = all_merge_cluster_complement_dict_start_end[_num][1]
        
        if cluster_dict != {}:
            for num, cluster_comp in cluster_dict.items():
                cluster_comp_start = cluster_dict_start_end[num][0]
                cluster_comp_end = cluster_dict_start_end[num][1]
                cluster_comp_complement = cluster_complement_dict[num] 
 
                if cluster_start <= cluster_comp_end:
                    overlap_num = 0
                    for mole_i in cluster:
                        for mole_j in cluster_comp:
                            if mole_i == mole_j:
                                overlap_num += 1
                                if overlap_num >= support_threshold:
                                    merge_cluster_flag = 1
                                    cluster_dict[num] = cluster.union(cluster_comp)
                                    cluster_complement_dict[num] = cluster_complement.union(cluster_comp_complement)
                                    cluster_new_start, cluster_new_end = Get_start_end_for_cluster(cluster_dict[num],mole_info)
                                    cluster_complement_new_start, cluster_complement_new_end = Get_start_end_for_cluster(cluster_complement_dict[num],mole_info)
                                    cluster_dict_start_end[num] = [cluster_new_start, cluster_new_end]
                                    cluster_complement_dict_start_end[num] = [cluster_complement_new_start, cluster_complement_new_end]

                                    break
                        if merge_cluster_flag == 1:
                            break
                    if merge_cluster_flag == 1:
                        break
            if merge_cluster_flag == 0:
                for num, cluster_comp in cluster_complement_dict.items():
                    cluster_comp_start = cluster_complement_dict_start_end[num][0]
                    cluster_comp_end = cluster_complement_dict_start_end[num][1]
                    cluster_comp_complement = cluster_dict[num] 
     
                    if cluster_start <= cluster_comp_end:
                        overlap_num = 0
                        for mole_i in cluster:
                            for mole_j in cluster_comp:
                                if mole_i == mole_j:
                                    overlap_num += 1
                                    if overlap_num >= support_threshold:
                                        merge_cluster_flag = 1
                                        cluster_complement_dict[num] = cluster.union(cluster_comp)
                                        cluster_dict[num] = cluster_complement.union(cluster_comp_complement)
                                        cluster_new_start, cluster_new_end = Get_start_end_for_cluster(cluster_complement_dict[num],mole_info)
                                        cluster_complement_new_start, cluster_complement_new_end = Get_start_end_for_cluster(cluster_dict[num],mole_info)
                                        cluster_complement_dict_start_end[num] = [cluster_new_start, cluster_new_end]
                                        cluster_dict_start_end[num] = [cluster_complement_new_start, cluster_complement_new_end]

                                        break
                            if merge_cluster_flag == 1:
                                break
                        if merge_cluster_flag == 1:
                            break
        
        if merge_cluster_flag == 0:
            cluster_dict[cluster_num] = cluster
            cluster_complement_dict[cluster_num] = cluster_complement
            cluster_dict_start_end[cluster_num] = [cluster_start, cluster_end]
            cluster_complement_dict_start_end[cluster_num] = [cluster_complement_start, cluster_complement_end]
            cluster_num += 1
   
    prev_num = total_num
    all_merge_cluster_dict = cluster_dict    
    all_merge_cluster_dict_start_end = cluster_dict_start_end 
    all_merge_cluster_complement_dict = cluster_complement_dict    
    all_merge_cluster_complement_dict_start_end = cluster_complement_dict_start_end   

    step += 1
    print("Clustering " + str(step) + ":")
    print(total_num,len(cluster_dict),len(cluster_complement_dict), total_num - len(cluster_dict))
    if prev_num == len(all_merge_cluster_dict):
        print("Converged...")
        print(len(all_merge_cluster_dict))
        Clustering_for_phaseblock(mole_dict,all_merge_cluster_dict, all_merge_cluster_complement_dict,h5_filename,overlap_var_threshold)

        return 1

    Recursive_clustering(all_merge_cluster_dict,all_merge_cluster_complement_dict,all_merge_cluster_dict_start_end, all_merge_cluster_complement_dict_start_end,step,mole_dict,mole_info,h5_filename,overlap_var_threshold,int(support_threshold))


def Get_start_end_for_cluster(one_cluster_mole_list,mole_info):
    count = 1
    for mole_num in one_cluster_mole_list:
        mole_start_and_end = mole_info[mole_num]
        mole_start = mole_start_and_end[0]
        mole_end = mole_start_and_end[1]
        if count == 1:
            start_min = mole_start
            end_max = mole_end
        else:
            if start_min > mole_start:
                start_min = mole_start 
            if end_max < mole_end:
                end_max = mole_end

        count += 1

    return (start_min, end_max)


def Sort_All_Clusters_by_Start(all_merge_cluster_dict,all_merge_cluster_dict_start_end):
    all_merge_cluster_dict_start = defaultdict(int)
    for key,val in all_merge_cluster_dict_start_end.items():
        all_merge_cluster_dict_start[key] = val[0]

    all_merge_cluster_dict_start_sorted = sorted(all_merge_cluster_dict_start.items(), key=operator.itemgetter(1))
    sorted_keys = []
    total_num = len(all_merge_cluster_dict_start_sorted)
    for ii in range(total_num):
        key_value = all_merge_cluster_dict_start_sorted[ii]
        key = key_value[0]
        sorted_keys.append(key)

    all_merge_cluster_dict_sorted = defaultdict(list)
    for key in sorted_keys:
        all_merge_cluster_dict_sorted[key] = all_merge_cluster_dict[key]

    return all_merge_cluster_dict_sorted


def process_chr(h5_file,overlap_var_threshold, support_threshold,var_depth_used,xin):
    support_threshold = int(support_threshold)
    f = open(h5_file,"r")
    curr = 0
    merge_dict = defaultdict(list)
    mole_dict = defaultdict()
    mole_info = defaultdict(list)
    for line in f:
        #print(curr)
        curr += 1
        data = line.rsplit()
        #chr_num = int(data[0])
        mole_num = int(data[6])
        mole_start = int(data[1])
        mole_end = int(data[2])
        var_list = data[7:]
        var_dict = defaultdict(int)
        var_hetero_pos = []
        var_hetero_hp = defaultdict(int)
        for var in var_list:
            var_info = var.split(":")
            var_pos = int(var_info[0])
            var_hp = int(var_info[1])
            var_dict[var_pos] = var_hp
            var_hetero_pos.append(var_pos)
            var_hetero_hp[var_pos] = var_hp
        var_hetero_pos_sorted = sorted(var_hetero_pos)

        for i in range(len(var_hetero_pos_sorted)-1):
            pos_1 = var_hetero_pos_sorted[i]
            pos_2 = var_hetero_pos_sorted[i+1]
            hp_1 = var_hetero_hp[pos_1]
            hp_2 = var_hetero_hp[pos_2]
            merge_dict[(pos_1,pos_2,hp_1,hp_2)].append(mole_num)

        """
        for i in range(len(var_hetero_pos_sorted)-1):
            for j in range(i + 1, len(var_hetero_pos_sorted)):
                pos_1 = var_hetero_pos_sorted[i]
                pos_2 = var_hetero_pos_sorted[j]
                hp_1 = var_hetero_hp[pos_1]
                hp_2 = var_hetero_hp[pos_2]
                merge_dict[(pos_1,pos_2,hp_1,hp_2)].append(mole_num)
        """
        mole_dict[mole_num] = var_dict
        mole_info[mole_num] = [mole_start,mole_end]

    count_merge = 0
    count_nomerge = 0
    all_merge_cluster_dict = defaultdict(list)
    all_merge_cluster_dict_start_end = defaultdict(list)
    for key,val in merge_dict.items():
        var_1_pos = key[0]
        var_2_pos = key[1]
        var_1_hp = key[2]
        var_2_hp = key[3]

        val_len = len(val)
        if (var_1_pos,var_2_pos,get_complement_hp_2(var_1_hp),get_complement_hp_2(var_2_hp)) in merge_dict:
            val2_len = len(merge_dict[(var_1_pos,var_2_pos,get_complement_hp_2(var_1_hp),get_complement_hp_2(var_2_hp))]) 
        else:
            val2_len = 0
            
        if (var_1_pos,var_2_pos,get_complement_hp_2(var_1_hp),var_2_hp) in merge_dict:
            val3_len = len(merge_dict[(var_1_pos,var_2_pos,get_complement_hp_2(var_1_hp),var_2_hp)]) 
        else:
            val3_len = 0
            
        if (var_1_pos,var_2_pos,var_1_hp,get_complement_hp_2(var_2_hp)) in merge_dict:
            val4_len = len(merge_dict[(var_1_pos,var_2_pos,var_1_hp,get_complement_hp_2(var_2_hp))]) 
        else:
            val4_len = 0
        
        p1 = 0.95
        p2 = 0.95
        n_total = val_len + val2_len + val3_len + val4_len
        k_false = val3_len + val4_len
        k_true = val_len + val2_len
        n_choose_k_false = int(comb(n_total, k_false, exact=True))
        total_combinations = 0
        for kk in range(n_total):
            n_choose_kk = int(comb(n_total, kk + 1, exact=True))
            total_combinations += n_choose_kk
        try:
            ratio = float(n_choose_k_false)/total_combinations
            final_prob = n_choose_k_false*math.pow((1-p1*p2),k_false)*math.pow((p1*p2),k_true)
        except:
            #print(val_len,val2_len,val3_len,val4_len)
            ratio = 1.0
            final_prob = 0.0

        if float(final_prob)/ratio > 0.9999 and (val_len >= var_depth_used and val2_len >= var_depth_used):
            count_merge += 1
            all_merge_cluster_dict[key] = val
            cluster_start_min, cluster_end_max = Get_start_end_for_cluster(val,mole_info)
            all_merge_cluster_dict_start_end[key] = [cluster_start_min, cluster_end_max]
            #print(val_len,val2_len,val3_len,val4_len)

        else:
            count_nomerge += 1
            #print(val_len,val2_len,val3_len,val4_len)



    #print(count_merge,count_nomerge)

    # sort all the clusters by start
    all_merge_cluster_dict_sorted = Sort_All_Clusters_by_Start(all_merge_cluster_dict,all_merge_cluster_dict_start_end)
    all_merge_cluster_dict = {}
    all_merge_cluster_dict = all_merge_cluster_dict_sorted
    # start first clustering
    cluster_dict = defaultdict(list)
    cluster_complement_dict = defaultdict(list)
    cluster_dict_start_end = defaultdict(list)
    cluster_complement_dict_start_end = defaultdict(list)
    curr = 0
    cluster_num = 0
    use_cluster_dict = defaultdict(int)
    total_num = len(all_merge_cluster_dict)
    
    for key,mole_num_list in all_merge_cluster_dict.items():
        if use_cluster_dict[key] != 1:
           # print(curr)
            curr += 1
            merge_cluster_flag = 0
            cluster = set(mole_num_list)
            cluster_start = all_merge_cluster_dict_start_end[key][0]
            cluster_end = all_merge_cluster_dict_start_end[key][1]
            var_1_pos = key[0]
            var_2_pos = key[1]
            var_1_hp = key[2]
            var_2_hp = key[3]
            key_complement = (var_1_pos,var_2_pos,get_complement_hp_2(var_1_hp),get_complement_hp_2(var_2_hp))  
            #if key_complement not in all_merge_cluster_dict:   # add
                #continue     # add
            cluster_complement = set(all_merge_cluster_dict[key_complement]) 
            cluster_complement_start = all_merge_cluster_dict_start_end[key_complement][0]
            cluster_complement_end = all_merge_cluster_dict_start_end[key_complement][1]
            use_cluster_dict[key_complement] = 1
            
            if cluster_dict != {}:
                for num, cluster_comp in cluster_dict.items():
                    cluster_comp_start = cluster_dict_start_end[num][0]
                    cluster_comp_end = cluster_dict_start_end[num][1]
                    cluster_comp_complement = cluster_complement_dict[num] 
     
                    if cluster_start <= cluster_comp_end:
                        overlap_num = 0
                        for mole_i in cluster:
                            for mole_j in cluster_comp:
                                if mole_i == mole_j:
                                    overlap_num += 1
                                    if overlap_num >= support_threshold:
                                        merge_cluster_flag = 1
                                        cluster_dict[num] = cluster.union(cluster_comp)
                                        cluster_complement_dict[num] = cluster_complement.union(cluster_comp_complement)
                                        cluster_new_start, cluster_new_end = Get_start_end_for_cluster(cluster_dict[num],mole_info)
                                        cluster_complement_new_start, cluster_complement_new_end = Get_start_end_for_cluster(cluster_complement_dict[num],mole_info)
                                        cluster_dict_start_end[num] = [cluster_new_start, cluster_new_end]
                                        cluster_complement_dict_start_end[num] = [cluster_complement_new_start, cluster_complement_new_end]

                                        break
                            if merge_cluster_flag == 1:
                                break
                        if merge_cluster_flag == 1:
                            break
                if merge_cluster_flag == 0:
                    for num, cluster_comp in cluster_complement_dict.items():
                        cluster_comp_start = cluster_complement_dict_start_end[num][0]
                        cluster_comp_end = cluster_complement_dict_start_end[num][1]
                        cluster_comp_complement = cluster_dict[num] 
         
                        if cluster_start <= cluster_comp_end:
                            overlap_num = 0
                            for mole_i in cluster:
                                for mole_j in cluster_comp:
                                    if mole_i == mole_j:
                                        overlap_num += 1
                                        if overlap_num >= support_threshold:
                                            merge_cluster_flag = 1
                                            cluster_complement_dict[num] = cluster.union(cluster_comp)
                                            cluster_dict[num] = cluster_complement.union(cluster_comp_complement)
                                            cluster_new_start, cluster_new_end = Get_start_end_for_cluster(cluster_complement_dict[num],mole_info)
                                            cluster_complement_new_start, cluster_complement_new_end = Get_start_end_for_cluster(cluster_dict[num],mole_info)
                                            cluster_complement_dict_start_end[num] = [cluster_new_start, cluster_new_end]
                                            cluster_dict_start_end[num] = [cluster_complement_new_start, cluster_complement_new_end]

                                            break
                                if merge_cluster_flag == 1:
                                    break
                            if merge_cluster_flag == 1:
                                break


            
            if merge_cluster_flag == 0:
                cluster_dict[cluster_num] = cluster
                cluster_complement_dict[cluster_num] = cluster_complement
                cluster_dict_start_end[cluster_num] = [cluster_start, cluster_end]
                cluster_complement_dict_start_end[cluster_num] = [cluster_complement_start, cluster_complement_end]
                cluster_num += 1
            
    print("Clustering " + str(1) + ":")
    #print(total_num,len(cluster_dict),len(cluster_complement_dict), total_num - len(cluster_dict))

    all_merge_cluster_dict = cluster_dict    
    all_merge_cluster_dict_start_end = cluster_dict_start_end 
    all_merge_cluster_complement_dict = cluster_complement_dict    
    all_merge_cluster_complement_dict_start_end = cluster_complement_dict_start_end   
    Recursive_clustering(all_merge_cluster_dict,all_merge_cluster_complement_dict,all_merge_cluster_dict_start_end, all_merge_cluster_complement_dict_start_end,1,mole_dict,mole_info,h5_file,overlap_var_threshold,support_threshold)


def Phase_start(output_dir,h5_dir,sample_name,chr_start,chr_end,overlap_var_threshold,support_threshold,xin):
    """ Algorithm to assign RAW phase blocks"""
    for chr_num in range(chr_start,chr_end + 1):
        rm_file = output_dir + "chr" + str(chr_num) + "_*"
        del_file_list = glob.glob(rm_file)
        if del_file_list != []:
            Popen("rm " + rm_file,shell=True).wait()
        rm_file = output_dir + "chr" + str(chr_num) + ".phased_final*"
        del_file_list = glob.glob(rm_file)
        print(del_file_list)
        if del_file_list != []:
            Popen("rm " + rm_file,shell=True).wait()

    print("using overlap variants threshold: " + str(overlap_var_threshold)) 
    print("using molecules support threshold: " + str(support_threshold)) 
    print("here: " + output_dir)
    #hetero_var_dict = pickle.load(open(h5_dir + "variant_dict_heterozygous.p","rb"))
    f = open(h5_dir + "median_depth_for_var.txt","r")
    for line in f :
        data = line.rsplit()
        var_depth_median = float(data[0])
    f.close()
    var_depth_used = var_depth_median*0.2
    for chr_num in range(chr_start,chr_end+1):
        print("processing " + str(chr_num) + "...")
        mole_h5_file_origin = h5_dir + sample_name + "_chr" + str(chr_num) + "_sorted.h5"
        total_filenum = split_h5(mole_h5_file_origin,chr_num,output_dir)
        for i in range(total_filenum):
            cur_filenum = i + 1
            cur_filename = output_dir + "chr" + str(chr_num) + "_" + str(cur_filenum)
            #print(cur_filename)
            process_chr(cur_filename,overlap_var_threshold,support_threshold,var_depth_used,"xin")

        """ Assign FINAL phase blocks """
        metric_phase_percent = []
        metric_corr_percent = []
        for i in range(total_filenum):
            cur_filenum = i + 1
            cur_filename = output_dir + "chr" + str(chr_num) + "_" + str(cur_filenum)
            phase_block_file = cur_filename + ".p"
            output_file_raw = cur_filename + ".phased.raw"
            output_file = cur_filename + ".phased"
            print(cur_filename)
            read_phase_block_file(phase_block_file,cur_filename,output_file_raw,metric_phase_percent, metric_corr_percent)
            Finalize_phase_block(output_file_raw, output_file)

        print("------Final Results for chr" + str(chr_num) + "-------")
        #print("phased percent: " + str(np.mean(metric_phase_percent)) + ", " +str(np.median(metric_phase_percent)))
        #print("correct percent: " + str(np.mean(metric_corr_percent)) + ", " +str(np.median(metric_corr_percent)))

        for file_ in glob.glob(output_dir + "chr" + str(chr_num) + "_*.phased"):
            Popen("cat " + file_ + ">> " +  output_dir + "chr" + str(chr_num) + ".phased_final_1",shell=True).wait()

        

        """ Assign phase block for molecule with only one heterzygous variant"""
        phase_file_final_1 = output_dir + "chr" + str(chr_num) + ".phased_final_1"

        Get_h5_for_mole_with_one_hetero_var(mole_h5_file_origin,chr_num,output_dir)

        h5_file_for_mole_with_one_var = output_dir + "chr" + str(chr_num) + "_one_var"
        phase_file_for_mole_with_one_var = h5_file_for_mole_with_one_var + ".phased"
        nonphase_dict = defaultdict(list) 

        var_prob_dict = Calculate_GenotypeProb_for_variants(phase_file_final_1)
        nonphase_dict = Assign_phase_block_for_mole_with_one_var(var_prob_dict,h5_file_for_mole_with_one_var,phase_file_for_mole_with_one_var,nonphase_dict)



        """ Impute the phase block for the rest very few molecules """
        phase_file_final_2 = output_dir +  "chr" + str(chr_num) + ".phased_final_2"

        nonphase_dict = Impute_phase_block(var_prob_dict,phase_file_final_1,phase_file_final_2,nonphase_dict)
        

        h5_file_for_mole_with_one_var = output_dir + "chr" + str(chr_num) + "_one_var"

        phase_file_for_mole_with_one_var_2 = h5_file_for_mole_with_one_var + ".phased_2"
        phase_file_final_3 = output_dir +  "chr" + str(chr_num) + ".phased_final_3"

        Impute_nonphase_variant(nonphase_dict,phase_file_for_mole_with_one_var,phase_file_final_2,phase_file_for_mole_with_one_var_2,phase_file_final_3)

        
        phase_file_final_total = output_dir +  "chr" + str(chr_num) + ".phased_final"
        #### cancatenate all phased files together ####
        Popen("cat " + phase_file_for_mole_with_one_var_2 + " " + phase_file_final_3 + " > " + phase_file_final_total,shell=True).wait()


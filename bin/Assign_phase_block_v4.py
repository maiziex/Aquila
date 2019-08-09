import pickle
from collections import defaultdict
import numpy as np


def read_phase_block_file(pickle_file,h5_file,output_file,metric_phase_percent,metric_corr_percent,hetero_var_dict):
    phase_block = pickle.load(open(pickle_file,"rb"))
    phase_block_2 = [] 
    total_num = len(phase_block)

    for _num in range(total_num):
        var_dict_2 = defaultdict(int)
        var_dict = phase_block[_num]
        for pos, hp in var_dict.items():
            if hp == 1 or hp == 2:
                var_dict_2[pos] = hp
                
        phase_block_2.append(var_dict_2)

    uniq_flag,convert_dict = Assign_phase_block(phase_block_2,h5_file,output_file,metric_phase_percent,metric_corr_percent)
    Refine_phase_block_for_output(uniq_flag,convert_dict,phase_block_2,h5_file,hetero_var_dict)


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


def get_complement_hp_3(hp):
    if hp == "hp1":
        return "hp2"
    elif hp == "hp2":
        return "hp1"


def Assign_phase_block(phase_block_new,h5_file,output_file,metric_phase_percent,metric_corr_percent):
    convert_dict = defaultdict(int)
    uniq_flag = defaultdict(int)
    fw = open(output_file,"w")
    total_blocks = len(phase_block_new)
    prob_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    set_HP_tag = 1001   # use an artificial HP_tag first
    for _num in range(total_blocks):
        one_block = phase_block_new[_num]
        PS_tag = set_HP_tag
        set_HP_tag += 1
        convert_dict[PS_tag] = _num
        HP_tag = 'hp1'
        HP_complement_tag = 'hp2'
        for var_pos, var_hp in one_block.items():
            prob_dict[var_pos][var_hp][(PS_tag,HP_tag)] = 0.99
            prob_dict[var_pos][get_complement_hp(var_hp)][(PS_tag,HP_complement_tag)] = 0.99

    count_phased = 0
    count_phased_error = 0
    #print(prob_dict)    
    f = open(h5_file, "r")
    curr = 0
    final_correct_percent = []
    for line in f:
        data = line.rsplit()
       # print(curr)
        curr += 1
        data = line.rsplit()
        if data[0] == "X":
            chr_num = 23
        else:
            chr_num = int(data[0])
        mole_num = int(data[6])
        mole_start = int(data[1])
        mole_end = int(data[2])
        var_list = data[7:]
        var_dict = defaultdict(int)
        for var in var_list:
            var_info = var.split(":")
            var_pos = int(var_info[0])
            var_hp = int(var_info[1])
            var_dict[var_pos] = var_hp + 1     # add 1 for operation

        mole_prob = defaultdict(int)
        for var_pos, var_hp in var_dict.items():
            prob_info = prob_dict[var_pos]
            if prob_info != {}:
                PS_info = prob_info[var_hp]
                for PS_tag_all, prob in PS_info.items():
                    mole_prob[PS_tag_all] += 1

        if len(mole_prob) >= 1:
            max_num = 0
            for PS_tag, num_of_prob in mole_prob.items():
                if max_num < num_of_prob:
                    use_PS_tag = PS_tag
                    max_num = num_of_prob

            min_num = mole_prob[(use_PS_tag[0],get_complement_hp_3(use_PS_tag[1]))]
            if max_num >= 3*min_num:   # may need update this threshold
                final_correct_percent.append(float(max_num)/(max_num + min_num))
                #print(curr, len(var_dict))
                #print(mole_prob)
                #print("------------------------------------------")
                #print(mole_prob[use_PS_tag])
                #print(min_num)
                #print(use_PS_tag)
                count_phased += 1
                #print(" ")
                fw.writelines("\t".join(data) + "\t" + "PS:" + str(use_PS_tag[0]) + "\t" + use_PS_tag[1] + "\n")
                if len(mole_prob) > 2:
                    for pos_,hp_ in var_dict.items():
                        _info = prob_dict[pos_][hp_]
                        if len(_info) > 1:
                            uniq_flag[pos_] = use_PS_tag


            else:
                count_phased_error += 1
                fw.writelines(line)

        else:
            fw.writelines(line)

    #print("---result---")
    metric_phase_percent.append(float(count_phased)/curr)
    metric_corr_percent.append(np.mean(final_correct_percent))
    #print(count_phased,count_phased_error,curr,float(count_phased)/curr)
    #print(np.mean(final_correct_percent))
    #print(len(phase_block_new))
    #return (metric_phase_percent,metric_corr_percent)
    return (uniq_flag,convert_dict)


def Refine_phase_block_for_output(uniq_flag,convert_dict,phase_block,h5_file,hetero_var_dict):
    output_file = h5_file + "_final.p"
    count_total = 0
    var_count = defaultdict(int)
    print(len(uniq_flag))
    num_of_blocks = len(phase_block)
    for _num in range(num_of_blocks):
        block_ = phase_block[_num]
        for _pos, _hp in block_.items():
            var_count[_pos] += 1
    for _pos, _count in var_count.items():
        if _count > 1:
            count_total += 1

    print(count_total)
    phase_block_eval = []
    # use the uniq_flag for overlap variants
    # for the final phase blocks, refine some overlap variants when maximize the haplotype likelihood of all molecules 
    chr_num_raw = h5_file.split("/chr")[1].split("_")[0]
    if chr_num_raw == "23":
        chr_num = "X"
    else:
        chr_num = int(chr_num_raw)
    for _num in range(num_of_blocks):
        block_new = defaultdict(int)
        block_ = phase_block[_num]
        for _pos, _hp in block_.items():
            if _pos in uniq_flag.keys():
                use_flag = uniq_flag[_pos]
                if convert_dict[use_flag[0]] == _num:
                    block_new[_pos] = _hp
            else:
                block_new[_pos] = _hp

        if block_new != {}:
            phase_block_eval.append(block_new)

    # use the refined phase blocks to do evaluation with gold standard
    # extract gold standard haplotype
    num_of_blocks_new = len(phase_block_eval)
    total_val = 0
    count_nofound = 0
    count_found = 0
    phase_block_eval_new = []
    count_wrong = 0
    for jj in range(num_of_blocks_new):
        block = phase_block_eval[jj]
        new_block = defaultdict(list)
        for _pos,_val in block.items():
            if (chr_num,_pos) in hetero_var_dict:
                new_block[_pos] = [_val,hetero_var_dict[(chr_num,_pos)][0],hetero_var_dict[(chr_num,_pos)][1]]
            else:
                count_wrong += 1
        phase_block_eval_new.append(new_block)  
        total_val += len(block)
    print(total_val)
    pickle.dump(phase_block_eval_new,open(output_file,"wb"))
    

def Finalize_phase_block(output_file_raw, output_file):
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
            HP_tag = data[-1]
            data[-2] = "PS:" + str(PS_flag_dict[PS_flag])
            fw.writelines("\t".join(data) + "\n")
        else:
            fw.writelines(line)

    f.close()
    fw.close()

    return PS_flag_dict


def Calculate_GenotypeProb_for_variants(phase_block_file):
    prob_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    f = open(phase_block_file,"r")
    curr = 0
    for line in f:
        #print(curr)
        curr += 1
        data = line.rsplit()
        flag = data[-1]
        if flag == "hp1" or flag == "hp2":
            var_list = data[7:-2]
            PS_tag = int(data[-2].split(":")[1])
            HP_tag = data[-1]
            HP_complement_tag = get_complement_hp_3(HP_tag)
            for one_var in var_list:
                var_info = one_var.split(":")
                var_pos = int(var_info[0])
                var_hp = int(var_info[1]) + 1    # add 1 for operation
                prob_dict[var_pos][var_hp][(PS_tag,HP_tag)] += 1
                prob_dict[var_pos][get_complement_hp(var_hp)][(PS_tag,HP_complement_tag)] += 1

    f.close()
    return prob_dict


def Assign_phase_block_for_mole_with_one_var(prob_dict,h5_one_var,output_file,nonphase_dict):
    fw = open(output_file,"w")
    f = open(h5_one_var,"r")
    for line in f:
        data = line.rsplit()
        var = data[-1]
        var_info = var.split(":")
        var_pos = int(var_info[0])
        var_hp = int(var_info[1]) + 1   # add 1 for operation

        prob_info = prob_dict[var_pos]
        if prob_info != {}:
            PS_info = prob_info[var_hp]
            curr = 0
            for PS_tag_all, count_ in PS_info.items():
                if curr == 0:
                    count_max = count_
                    PS_HP_use = PS_tag_all
                else:
                    if count_ > count_max:
                        PS_HP_use = PS_tag_all

                curr += 1
            PS_tag_use = PS_HP_use[0]
            HP_tag_use = PS_HP_use[1]
            fw.writelines("\t".join(data) + "\t" + "PS:" + str(PS_tag_use) + "\t" + HP_tag_use + "\n")
        else:
            fw.writelines(line)
            nonphase_dict[data[6]] = data[7:]


    fw.close()
    f.close()
    return nonphase_dict


def Impute_phase_block(prob_dict,phase_file,output_file,nonphase_dict):
    f = open(phase_file,"r")
    fw = open(output_file,"w")
    for line in f:
        data = line.rsplit()
        flag = data[-1]
        if flag != "hp1" and flag != "hp2":
            var_list = data[7:]
            mole_prob = defaultdict(int)
            key_use = []
            for one_var in var_list:
                var_pos = int(one_var.split(":")[0])
                var_hp = int(one_var.split(":")[1]) + 1 # add 1 for operation
                prob_info = prob_dict[var_pos]
                if prob_info != {}:
                    PS_info = prob_info[var_hp]
                    curr_ = 0
                    for PS_tag_all, count_ in PS_info.items():
                        if curr_ == 0:
                            count_max = count_
                            PS_HP_use = PS_tag_all
                        else:
                            if count_ > count_max:
                                PS_HP_use = PS_tag_all

                        curr_ += 1

                    PS_tag_use = PS_HP_use[0]
                    HP_tag_use = PS_HP_use[1]

                    mole_prob[(PS_tag_use,HP_tag_use)] += 1

                    count_val = 0 
                    for key, val in mole_prob.items():
                        if count_val == 0:
                            val_max = val
                            key_use = key
                        else:
                            if val > val_max:
                                val_max = val
                                key_use = key
                        count_val += 1
            if key_use == []:
                fw.writelines(line)
                nonphase_dict[data[6]] = data[7:]
            else:
                fw.writelines("\t".join(data) + "\t" + "PS:" + str(key_use[0]) + "\t" + key_use[1] + "\n")
        else:
            fw.writelines(line)

    f.close()
    fw.close()
    return nonphase_dict


def write_phase_block_into_h5(phase_block_final,phase_file,output_file):
    f = open(phase_file,"r")
    fw = open(output_file,"w")
    for line in f:
        data = line.rsplit()
        flag = data[-1]
        if flag != "hp1" and flag != "hp2":
            mole_num = int(data[6])
            phase_info = phase_block_final[mole_num]
            use_PS_tag = phase_info[0]
            use_HP_tag = phase_info[1]
            fw.writelines("\t".join(data) + "\t" + "PS:" + str(use_PS_tag) + "\t" + use_HP_tag + "\n")
        else:
            fw.writelines(line)
    f.close()
    fw.close()
    

def Impute_nonphase_variant(nonphase_dict,phase_file_1,phase_file_2,output_file_1,output_file_2):
    ###nonphase_dict = pickle.load(open("nonphase_dict.p", "rb"))
    merge_dict = defaultdict(list)
    mole_var_pos_min_dict = defaultdict(int)
    for key, var_list in nonphase_dict.items():
        mole_num = int(key)
        count_ = 0
        for one_var in var_list:
            var_pos = int(one_var.split(":")[0])
            var_hp = int(one_var.split(":")[1])
            merge_dict[(var_pos, var_hp)].append(mole_num)
            if count_ == 0:
                var_pos_min = var_pos
            else:
                if var_pos_min > var_pos:
                    var_pos_min = var_pos
            count_ += 1
        mole_var_pos_min_dict[mole_num] = var_pos_min

    phase_block = defaultdict(lambda: defaultdict(int))
    count_phase = defaultdict(int)
    use_key = defaultdict(int)
    for key, mole_list in merge_dict.items():
        if use_key[key] == 0:
            set_PS_tag = key[0]
            key_complement = (key[0],get_complement_hp_2(key[1]))
            for mole_num in mole_list:
                phase_block[mole_num][(set_PS_tag,"hp1")] += 1
                count_phase[(set_PS_tag,"hp1")] += 1
            if key_complement in merge_dict:
                use_key[key_complement] = 1
                mole_list_complement = merge_dict[key_complement]
                for mole_num in mole_list_complement:
                    phase_block[mole_num][(set_PS_tag,"hp2")] += 1
                    count_phase[(set_PS_tag,"hp2")] += 1


    phase_block_use = defaultdict(list)
    for mole_num, phase_info_all in phase_block.items():
        use_phase_info = ["",""]
        curr_ = 0
        for phase_info,val in phase_info_all.items():
            count_1 = count_phase[phase_info]
            #phase_info_complement = (phase_info[0],get_complement_hp_3(phase_info[1]))
            #count_2 = count_phase[phase_info_complement]
            #count_diff = abs(count_1 - count_2)
            count_diff = count_1
            if curr_ == 0:
                count_diff_max = count_diff
                use_phase_info = phase_info
            else:
                if count_diff_max < count_diff:
                    count_diff_max = count_diff
                    use_phase_info = phase_info

            curr_ += 1
        use_PS_tag = use_phase_info[0]
        use_HP_tag = use_phase_info[1]
        phase_block_use[mole_num] = [use_PS_tag, use_HP_tag] 

    save_min_pos = defaultdict(list)
    save_min_pos_2 = defaultdict(int)
    for mole_num, val in phase_block_use.items():
        min_pos = mole_var_pos_min_dict[mole_num]
        save_min_pos[(val[0],val[1])].append(min_pos)

    for key, val in save_min_pos.items():
        save_min_pos_2[key] = np.min(val)

    phase_block_final = defaultdict(list)
    for mole_num, val in phase_block_use.items():
        use_HP_tag = val[1]
        use_PS_tag = save_min_pos_2[(val[0],val[1])]
        phase_block_final[mole_num]= [use_PS_tag, use_HP_tag]

    
    write_phase_block_into_h5(phase_block_final, phase_file_1,output_file_1)
    write_phase_block_into_h5(phase_block_final, phase_file_2,output_file_2)



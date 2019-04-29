import pdb
#pdb.set_trace()
from collections import defaultdict
import pickle
import os
from argparse import ArgumentParser

def save_pickle_file(dict1,fname):
    for value in dict1:
        dict1[value] = dict(dict1[value])
    my_dict = dict(dict1)
    with open(fname, "wb") as f:
        pickle.dump(my_dict,f) 


def Extract_mole_num_from_phased_mole(phased_file,PS_flag_dict_cut_file,chr_num):
    phase_dict = defaultdict(list)
    PS_flag_dict_cut = pickle.load(open(PS_flag_dict_cut_file,"rb"))
    f = open(phased_file,"r")
    for line in f:
        data = line.rsplit()
        PS_flag_info = data[-2]
        PS_flag = int(PS_flag_info.split(":")[1])
        HP_flag = data[-1]
        mole_num = int(data[6])
        PS_flag_end = PS_flag_dict_cut[PS_flag]
        phase_dict[(PS_flag,PS_flag_end,HP_flag)].append(mole_num)

    #pickle.dump(phase_dict,open("phase_dict_chr" + str(chr_num) + ".p","wb"))
    print("done~")
    return phase_dict

def check_qname_in_PS(block_info,qname_pos_list):
    block_start = block_info[0]
    block_end = block_info[1]
    flag_in = 0
    for pos in qname_pos_list:
        if pos <= block_end and pos >= block_start:
            flag_in = 1
            break
    return flag_in 


def Extract_qname(phased_dict_mole_num,mole_qname_dict,qname_pos_dict,barcoded_fastq_file,chr_num,output_dir):
    phased_dict_qname_2 = defaultdict(lambda: defaultdict(int))
    qname_total_dict = defaultdict(int)
    curr = 0
    flag_in = 0
    for key, mole_num_list in phased_dict_mole_num.items():
        curr += 1
        for mole_num in mole_num_list:
            qname_list = mole_qname_dict[mole_num]
            for qname in qname_list:
                flag_in = check_qname_in_PS(key,qname_pos_dict[qname])
                if flag_in == 1:
                    phased_dict_qname_2[qname][key] = 1
    #save_pickle_file(phased_dict_qname_2,"phased_dict_qname_2.p")
    
    f = open(barcoded_fastq_file,"r")
    count = 0
    flag = 0
    for line in f:
        data = line.rsplit()
        if count%8 == 0:
            qname_curr = data[0][1:]
            if phased_dict_qname_2[qname_curr] != {}:
                flag = 1
                PS_flag_info = phased_dict_qname_2[qname_curr]
                for _PS_HP, _val in PS_flag_info.items():
                    _PS_flag = _PS_HP[0]
                    _PS_flag_end = _PS_HP[1]
                    _HP_flag = _PS_HP[2]
                file_curr = output_dir +"fastq_by_" + str(_PS_flag) + "_" + str(_PS_flag_end) + "_" + _HP_flag + ".fastq"
                if os.path.isfile(file_curr):
                    fw_curr = open(file_curr,"a")
                    fw_curr.write(line)
                else:
                    fw_curr = open(file_curr,"w")
                    fw_curr.write(line)
                del phased_dict_qname_2[qname_curr]

            else:
                del phased_dict_qname_2[qname_curr]

        elif count%8 == 7:
            if flag == 1:
                flag = 0
                fw_curr.write(line)
                fw_curr.close()

        else:
            if flag == 1:
                fw_curr.write(line)

        count += 1

    print("finished extracting...")



def Extract_start(output_dir,chr_num,phased_h5_file,PS_flag_dict_cut_file,mole_qname_dict_file,qname_pos_dict_file,chr_fastq,xin):
    phased_dict_mole_num = Extract_mole_num_from_phased_mole(phased_h5_file,PS_flag_dict_cut_file,chr_num)
    mole_qname_dict = pickle.load(open(mole_qname_dict_file,"rb"))
    qname_pos_dict = pickle.load(open(qname_pos_dict_file,"rb"))
    Extract_qname(phased_dict_mole_num,mole_qname_dict,qname_pos_dict,chr_fastq,chr_num,output_dir)

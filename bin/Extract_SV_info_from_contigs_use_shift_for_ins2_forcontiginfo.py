import pdb
#pdb.set_trace()
from collections import defaultdict
import pickle
import numpy as np
from argparse import ArgumentParser
import os
import sys
parser = ArgumentParser(description="Variants Calling:")
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--var_size','-v',type=int,help="variant size", default=1)
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs", default='results_/')

args = parser.parse_args()


def Extract_SV_info(hap_SV_file,output_file,chr_num,v_size):
    f = open(hap_SV_file,"r")
    contig_dict = defaultdict(list)
    curr = 0
    for line in f:
        #print(curr)
        curr += 1
        data = line.rsplit()
        if data[0] == "V":
            #print(data)
            mapq = int(data[5])
            repeat = int(data[4])
            cur_chr_num = data[1]
            if cur_chr_num == "chrX":
                cur_chr_num = "chr23"
            if mapq >= 20 and repeat == 1 and cur_chr_num == "chr" + str(chr_num):
                cur_chr_num = data[1]
                ref_start = int(data[2])
                ref_end = int(data[3])
                ref = data[6]
                alt = data[7]
                contig_num = data[8]
                contig_start = int(data[9])
                contig_end = int(data[10])
                strand_dir = data[-1]
                contig_dict[(cur_chr_num,ref_start)].append([cur_chr_num,ref_start, ref_end, ref, alt, contig_num,contig_start,contig_end,strand_dir])

    contig_dict_filter = defaultdict(list)
    count = 0
    for key, val in contig_dict.items():
        if len(val) > 1:
            #print(">1:")
            #print(val)
            count += 1
        else:
            contig_dict_filter[key] = val[0]


    print(len(contig_dict))
    print(len(contig_dict_filter))

    
    SV_dict = defaultdict(list)
    for key, val in contig_dict_filter.items():
        ref = val[3] 
        alt = val[4]
        len_var = abs(len(ref)-len(alt))
        if (ref == "-" and len(alt) >= v_size) or (alt == "-" and len(ref) >= v_size):
            SV_dict[key] = val
        
    pickle.dump(SV_dict, open(output_file,"wb"))

    print(len(SV_dict))
    return SV_dict


def compare_two_haploid_SV(SV_dict_contig_1_file, SV_dict_contig_2_file,use_chr_num,ref_dict_contig_1_file,ref_dict_contig_2_file,out_dir):
    ins_homo_dict = defaultdict(int)
    ins_hetero_dict = defaultdict(int)
    ins_homo_sv = defaultdict(list)
    ins_homo_sv_compound = defaultdict(list)
    ins_hetero_sv = defaultdict(list)
    SV_dict_contig_1 = pickle.load(open(SV_dict_contig_1_file,"rb"))
    SV_dict_contig_2 = pickle.load(open(SV_dict_contig_2_file,"rb"))
    ref_dict_contig_1 = pickle.load(open(ref_dict_contig_1_file,"rb"))
    ref_dict_contig_2 = pickle.load(open(ref_dict_contig_2_file,"rb"))
    count_total_share = 0
    homo_key_dict = defaultdict(int)
    count_same_wrong = 0
    for key, val in SV_dict_contig_1.items():
        flag = 0
        if key in SV_dict_contig_2:
            #print(key)
            homo_key_dict[key] = 1
            val_2 = SV_dict_contig_2[key]
            chr_num = val[0]
            start_1 = val[1]
            end_1 = val[2]
            start_2 = val_2[1]
            end_2 = val_2[2]
            ref_1 = val[3]
            ref_2 = val_2[3]
            alt_1 = val[4]
            alt_2 = val_2[4]
            PS_start_1 = int(val[5].split("_PS")[-1].split("_hp1")[0].split(":")[0])
            PS_end_1 = int(val[5].split("_PS")[-1].split("_hp1")[0].split(":")[1])
            PS_start_2 = int(val_2[5].split("_PS")[-1].split("_hp2")[0].split(":")[0])
            PS_end_2 = int(val_2[5].split("_PS")[-1].split("_hp2")[0].split(":")[1])
            if start_1 == start_2 and end_1 == end_2 and ref_1 == ref_2 and alt_1 == alt_2:
                if ref_1 == "-": # ins
                    if PS_end_1 < PS_start_2 or PS_end_2 < PS_start_1:
                        count_same_wrong += 1
                    else:
                        sv_len_1 = len(alt_1)
                        contig_num_1 = val[-4]
                        contig_start_1 = val[-3]
                        contig_end_1 = val[-2]
                        contig_strand_dir_1 = val[-1]
                        contig_num_2 = val_2[-4]
                        contig_start_2 = val_2[-3]
                        contig_end_2 = val_2[-2]
                        contig_strand_dir_2 = val_2[-1]
                        ins_homo_dict[key] = sv_len_1
                        #ins_homo_sv[key] = [chr_num,start_1,end_1,val[5],val[6],val[7],val_2[5],val_2[6],val_2[7]]
                        ins_homo_sv[key] = [chr_num,start_1,end_1,ref_1,alt_1,contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]

            else:
                if ref_1 == "-": # ins
                    if PS_end_1 < PS_start_2 or PS_end_2 < PS_start_1:
                        count_same_wrong += 1
                    else:
                        sv_len_1 = len(alt_1)
                        sv_len_2 = len(alt_2)
                        contig_num_1 = val[-4]
                        contig_start_1 = val[-3]
                        contig_end_1 = val[-2]
                        contig_strand_dir_1 = val[-1]
                        contig_num_2 = val_2[-4]
                        contig_start_2 = val_2[-3]
                        contig_end_2 = val_2[-2]
                        contig_strand_dir_2 = val_2[-1]
                        ins_homo_dict[key] = np.min([sv_len_1,sv_len_2])
                        #if ref_1 == "-" and alt_2 == "-":
                         #   print(ref_1,alt_1)
                        #ins_homo_sv[key] = [chr_num,start_1,end_1,val[5],val[6],val[7],val_2[5],val_2[6],val_2[7]]
                        ins_homo_sv_compound[key] = [chr_num,[start_1,start_2],[end_1,end_2],[ref_1,ref_2],[alt_1,alt_2],contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]

        else:
            chr_num = val[0]
            start_1 = val[1]
            end_1 = val[2]
            ref_1 = val[3]
            alt_1 = val[4]
            sv_len = abs(len(ref_1) - len(alt_1)) + 1
            for _step in range(start_1 - 10, start_1 + 10):
                key_shift = (chr_num,_step)
                if key_shift in SV_dict_contig_2:
                    val_2 = SV_dict_contig_2[key_shift]
                    start_2 = val_2[1]
                    end_2 = val_2[2]
                    ref_2 = val_2[3]
                    alt_2 = val_2[4]
                    if abs(end_1 - end_2) <= 10:
                        flag = 1
                        homo_key_dict[key_shift] = 1
                        if ref_1 == "-": # ins
                            sv_len_1 = len(alt_1)
                            sv_len_2 = len(alt_2)
                            contig_num_1 = val[-4]
                            contig_start_1 = val[-3]
                            contig_end_1 = val[-2]
                            contig_strand_dir_1 = val[-1]
                            contig_num_2 = val_2[-4]
                            contig_start_2 = val_2[-3]
                            contig_end_2 = val_2[-2]
                            contig_strand_dir_2 = val_2[-1]
                            ins_homo_dict[key_shift] = np.min([sv_len_1,sv_len_2])
                            #ins_homo_sv[key_shift] = [chr_num,start_1,end_1,val[5],val[6],val[7],val_2[5],val_2[6],val_2[7]]
                            ins_homo_sv_compound[key_shift] = [chr_num,[start_1,start_2],[end_1,end_2],[ref_1,ref_2],[alt_1,alt_2],contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]

                        break
            if flag == 0:
                if ref_1 == "-": # ins
                    count_cov = 0
                    for _step in range(start_1 - 20, start_1 + 20):
                        if _step in ref_dict_contig_2:
                            count_cov += 1
                    if count_cov == 40:
                        sv_len_1 = len(alt_1)
                        ins_hetero_dict[key] = sv_len_1
                        contig_start_1 = val[-3]
                        contig_end_1 = val[-2]
                        contig_num_1 = val[-4]
                        contig_strand_dir_1 = val[-1]
                        #ins_hetero_sv[key] = [chr_num,start_1,end_1,val[5],val[6],val[7],0,1,2]
                        ins_hetero_sv[key] = [chr_num,start_1,end_1,ref_1,alt_1,contig_num_1,contig_start_1,contig_end_1,contig_strand_dir_1,'0_PS0:0_hp0',0,0,0]



    # processing the uniq in contig_2
    for key_2, val_2 in SV_dict_contig_2.items():
        if homo_key_dict[key_2] == 1:
            count_total_share += 1
        else:
            start_2 = val_2[1]
            end_2 = val_2[2]
            ref_2 = val_2[3]
            alt_2 = val_2[4]
            if ref_2 == "-": # ins
                count_cov = 0
                for _step in range(start_1 - 20, start_1 + 20):
                    if _step in ref_dict_contig_1:
                        count_cov += 1
                if count_cov == 40:
                    sv_len_2 = len(alt_2)
                    ins_hetero_dict[key_2] = sv_len_2
                    contig_num_2 = val_2[-4]
                    contig_start_2 = val_2[-3]
                    contig_end_2 = val_2[-2]
                    contig_strand_dir_2 = val_2[-1]
                    #ins_hetero_sv[key_2] = [val_2[0],start_2,end_2,0,1,2,val_2[5],val_2[6],val_2[7]]
                    #ins_hetero_sv[key_2] = [val_2[0],start_2,end_2,ref_2,alt_2]
                    ins_hetero_sv[key_2] = [val_2[0],start_2,end_2,ref_2,alt_2,'0_PS0:0_hp0',0,0,0,contig_num_2,contig_start_2,contig_end_2,contig_strand_dir_2]


    print("done")
    pickle.dump(ins_homo_sv, open(out_dir + "ins_homo_sv_chr" + str(use_chr_num) + ".p","wb"))
    pickle.dump(ins_homo_sv_compound, open(out_dir + "ins_homo_sv_compound_chr" + str(use_chr_num) + ".p","wb"))
    pickle.dump(ins_hetero_sv, open(out_dir + "ins_hetero_sv_chr" + str(use_chr_num) + ".p","wb"))
    return count_same_wrong


def Write_vcf_for_SV(fw,ins_homo_file,ins_homo_file_compound,ins_hetero_file,count_id):
    ins_homo = pickle.load(open(ins_homo_file,"rb"))
    ins_homo_compound = pickle.load(open(ins_homo_file_compound,"rb"))
    ins_hetero = pickle.load(open(ins_hetero_file,"rb"))
    for key, val in ins_homo.items():
        #print(key)
        chr_num = str(val[0])
        start_ = val[1]
        end_ = val[2]
        ref = val[3]
        alt = val[4]
        GT = "1/1"
        contig_num_1 = val[5]
        contig_start_1 = val[6]
        contig_end_1 = val[7]
        contig_strand_dir_1 = val[8]
        contig_num_2 = val[9]
        contig_start_2 = val[10]
        contig_end_2 = val[11]
        contig_strand_dir_2 = val[12]
        new_contig_num_1 = contig_num_1.replace(":","_")
        new_contig_num_2 = contig_num_2.replace(":","_")
        contig_info = str(new_contig_num_1) + "_" + str(contig_start_1) + "_" + str(contig_end_1) + "_" + contig_strand_dir_1 + "_" +  str(new_contig_num_2) + "_" + str(contig_start_2) + "_" + str(contig_end_2) + "_" + contig_strand_dir_2 
        fw.writelines(chr_num + "\t" + str(start_) + "\t" + "event" + str(count_id) + "\t" + ref + "\t" + alt + "\t" + "." + "\t" + "PASS" + "\t" + "SVTYPE=INS" + "\t" + "GT:Contig" + "\t"  + "1/1"  + ":" + contig_info +  "\n")
        count_id += 1

    for key, val in ins_homo_compound.items():
        #print(key)
        chr_num = str(val[0])
        start_1 = val[1][0]
        start_2 = val[1][1]
        end_1 = val[2][0]
        end_2 = val[2][1]
        ref_1 = val[3][0]
        ref_2 = val[3][1]
        alt_1 = val[4][0]
        alt_2 = val[4][1]
        alt = val[4]
        GT = "1/2"
        contig_num_1 = val[5]
        contig_start_1 = val[6]
        contig_end_1 = val[7]
        contig_strand_dir_1 = val[8]
        contig_num_2 = val[9]
        contig_start_2 = val[10]
        contig_end_2 = val[11]
        contig_strand_dir_2 = val[12]
        new_contig_num_1 = contig_num_1.replace(":","_")
        new_contig_num_2 = contig_num_2.replace(":","_")
        contig_info = str(new_contig_num_1) + "_" + str(contig_start_1) + "_" + str(contig_end_1) + "_" + contig_strand_dir_1 + "_" +  str(new_contig_num_2) + "_" + str(contig_start_2) + "_" + str(contig_end_2) + "_" + contig_strand_dir_2 
        #fw.writelines(chr_num + "\t" + str(start_1) + "," + str(start_2) + "\t" + "event" + str(count_id) + "\t" + ref_1 + "," + ref_2 + "\t" + alt_1 + "," + alt_2 + "\t" + "." + "\t" + "PASS" + "\t" + "SVTYPE=INS" + "\t" + "GT:Contig" + "\t"  + "1/2"  + ":" + contig_info +  "\n")
        fw.writelines(chr_num + "\t" + str(start_1)  + "\t" + "event" + str(count_id) + "\t" + ref_1 + "\t" + alt_1 + "\t" + "." + "\t" + "PASS" + "\t" + "SVTYPE=INS" + "\t" + "GT:Contig" + "\t"  + "1/2"  + ":" + contig_info +  "\n")
        count_id += 1


    for key, val in ins_hetero.items():
        #print(key)
        chr_num = str(val[0])
        start_ = val[1]
        end_ = val[2]
        ref = val[3]
        alt = val[4]
        GT = "0/1"
        contig_num_1 = val[5]
        contig_start_1 = val[6]
        contig_end_1 = val[7]
        contig_strand_dir_1 = val[8]
        contig_num_2 = val[9]
        contig_start_2 = val[10]
        contig_end_2 = val[11]
        contig_strand_dir_2 = val[12]
        new_contig_num_1 = contig_num_1.replace(":","_")
        new_contig_num_2 = contig_num_2.replace(":","_")
        contig_info = str(new_contig_num_1) + "_" + str(contig_start_1) + "_" + str(contig_end_1) + "_" + str(contig_strand_dir_1) + "_" +  str(new_contig_num_2) + "_" + str(contig_start_2) + "_" + str(contig_end_2) + "_" + str(contig_strand_dir_2) 
        fw.writelines(chr_num + "\t" + str(start_) + "\t" + "event" + str(count_id) + "\t" + ref + "\t" + alt + "\t" + "." + "\t" + "PASS" + "\t" + "SVTYPE=INS" + "\t" + "GT:Contig" + "\t"  + "0/1"  + ":"   + contig_info +  "\n")
        count_id += 1

 
    return (fw,count_id)



if __name__ == "__main__":
    out_dir = args.out_dir + "/" 
    if os.path.exists(out_dir):
        print("using existing output folder: " + out_dir)
    else:
        os.makedirs(out_dir)
    chr_start = args.chr_start
    chr_end = args.chr_end
    v_size = args.var_size

    total_same_wrong = 0
    count_id = 1
    for chr_num in range(chr_start,chr_end + 1):
        final_SV_vcf = out_dir + "Aquila_INS_chr" + str(chr_num)  + ".vcf"
        fw = open(final_SV_vcf,"w")
        SV_dict_contig_1 = Extract_SV_info(out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp1.var.txt",out_dir + "ins_dict_contig_1_chr" + str(chr_num) + ".p", chr_num,v_size)
        SV_dict_contig_2 = Extract_SV_info(out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp2.var.txt",out_dir + "ins_dict_contig_2_chr" + str(chr_num) + ".p", chr_num,v_size)
        count_same_wrong = compare_two_haploid_SV(out_dir + "ins_dict_contig_1_chr" + str(chr_num) + ".p", out_dir + "ins_dict_contig_2_chr" + str(chr_num) + ".p", chr_num, out_dir + "ref_dict_contig_1_chr" + str(chr_num) + ".p", out_dir + "ref_dict_contig_2_chr" + str(chr_num) + ".p",out_dir )
        total_same_wrong += count_same_wrong


        fw,count_id = Write_vcf_for_SV(fw,out_dir + "ins_homo_sv_chr" + str(chr_num) + ".p",out_dir + "ins_homo_sv_compound_chr" + str(chr_num) + ".p" , out_dir + "ins_hetero_sv_chr" +str(chr_num) + ".p",count_id)

        fw.close()
    print("total_same_wrong:")
    print(total_same_wrong)

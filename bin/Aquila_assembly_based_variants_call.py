import pdb
#pdb.set_trace()
import os
from argparse import ArgumentParser
from subprocess import Popen, PIPE
from multiprocessing import Pool,cpu_count,active_children,Manager
import time
import sys
parser = ArgumentParser(description="Run depth all:")
parser.add_argument('--chr_start','-start',type=int,help="chromosome start from", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by", default=23)
parser.add_argument('--num_of_threads','-t',type=int,help="number of threads", default=1)
parser.add_argument('--assembly_dir','-i_dir', help="Directory for inputs",required=True)
parser.add_argument('--out_dir','-o_dir', help="Directory to store outputs", default='Variants_results/')
parser.add_argument('--ref_file','-r', help="reference fasta file",required=True)

args = parser.parse_args()
script_path = os.path.dirname(os.path.abspath( __file__ ))
code_path = script_path + "/" 
__author__ = "Xin Zhou@Stanford"

def Split_supercontig_by_haplotype(fasta_file,xin):
    f = open(fasta_file,"r")
    output_prefix = fasta_file.split(".fasta")[0]
    hp1_file = output_prefix + "_hp1.fasta"
    hp2_file = output_prefix + "_hp2.fasta"
    fw_hp1 = open(hp1_file,"w")
    fw_hp2 = open(hp2_file,"w")
    count = 0
    for line in f:
        data = line.rsplit()
        if count%2 == 0:
            hp_flag = data[0].split("_")[-1]
            if hp_flag == "hp1":
                fw_hp1.writelines(line)
            elif hp_flag == "hp2":
                fw_hp2.writelines(line)
        elif count%2 == 1:
            if hp_flag == "hp1":
                fw_hp1.writelines(line)
            elif hp_flag == "hp2":
                fw_hp2.writelines(line)
        count += 1
    f.close()
    fw_hp1.close()
    fw_hp2.close()
    print("done")


def Split_haplotype(chr_start,chr_end,in_dir):
    pool = Pool(processes=chr_end - chr_start +1)
    for chr_num in range(chr_start,chr_end + 1):
        input_file = in_dir + "Aquila_Contig_chr" + str(chr_num) + ".fasta"
        pool.apply_async(Split_supercontig_by_haplotype,(input_file,"xin"))
    pool.close()
    while len(active_children()) > 1:
        time.sleep(0.5)
    pool.join()
    print("all done~")


def get_paf(ref_file,hap_file,hap_paf,xin):
    try:
        use_cmd =  "minimap2 -cx asm5 -t8 --cs " + ref_file + " " + hap_file  + " > " + hap_paf
    except:
        use_cmd =  code_path + "minimap2/" + "minimap2 -cx asm5 -t8 --cs " + ref_file + " " + hap_file  + " > " + hap_paf
    Popen(use_cmd,shell=True).wait()


def sort_paf(hap_paf,hap_paf_sorted,xin):
    use_cmd = "sort -k6,6 -k8,8n " + hap_paf + " > "  + hap_paf_sorted  
    Popen(use_cmd,shell=True).wait()


def get_var(hap_paf_sorted,hap_var_txt,xin):
    use_cmd = code_path + "/k8-0.2.4/k8-Linux "  + code_path + "paftools/" +  "/paftools.js " + " call -l 1 -L 1 -q 20 " +  hap_paf_sorted + " > " +  hap_var_txt 
    Popen(use_cmd,shell=True).wait()
    

def assembly_based_variants_call_paf(chr_start,chr_end,ref_file):
    for chr_num in range(chr_start,chr_end + 1):
        hap1_file = in_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp1.fasta"
        hap2_file = in_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp2.fasta"
        hap1_paf = out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp1.paf"
        hap2_paf = out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp2.paf"
        pool = Pool(processes=2)
        pool.apply_async(get_paf,(ref_file,hap1_file,hap1_paf,"xin"))
        pool.apply_async(get_paf,(ref_file,hap2_file,hap2_paf,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()

    print("all done~")


def assembly_based_variants_call_sort(chr_start,chr_end):
    for chr_num in range(chr_start,chr_end + 1):
        hap1_paf = out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp1.paf"
        hap2_paf = out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp2.paf"
        hap1_paf_sorted = out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp1.paf.sorted" 
        hap2_paf_sorted = out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp2.paf.sorted" 
        pool = Pool(processes=2)
        pool.apply_async(sort_paf,(hap1_paf,hap1_paf_sorted,"xin"))
        pool.apply_async(sort_paf,(hap2_paf,hap2_paf_sorted,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()

    print("all done~")
 

def assembly_based_variants_call(chr_start,chr_end):
    for chr_num in range(chr_start,chr_end + 1):
        hap1_paf_sorted = out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp1.paf.sorted"
        hap2_paf_sorted = out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp2.paf.sorted"
        hap1_var_txt = out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp1.var.txt"
        hap2_var_txt = out_dir + "Aquila_Contig_chr" + str(chr_num) + "_hp2.var.txt"
        pool = Pool(processes=2)
        pool.apply_async(get_var,(hap1_paf_sorted,hap1_var_txt,"xin"))
        pool.apply_async(get_var,(hap2_paf_sorted,hap2_var_txt,"xin"))
        pool.close()
        while len(active_children()) > 1:
            time.sleep(0.5)
        pool.join()

    print("all done~")


def Run_del_or_ins(all_cmd,xin):
    Popen(all_cmd,shell=True).wait()


def Call_SNV_info_from_contigs(chr_start,chr_end,out_dir,num_of_threads):
    total_num = chr_end - chr_start +1
    pool = Pool(num_of_threads)
    count = 1
    for chr_num in range(chr_start,chr_end + 1):
        all_cmd = "python3 " + code_path + "Extract_SNV_info_from_contigs_forcontiginfo_forall.py  --chr_start " + str(chr_num) + " --chr_end " + str(chr_num) + " --out_dir " + out_dir 
        count += 1
        pool.apply_async(Run_del_or_ins,(all_cmd,"xin"))
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()

            if (count - 1) == total_num:
                print("finished chr" + str(chr_num))
            else:
                pool = Pool(num_of_threads)
    curr_file = out_dir + "Aquila_SNPs.vcf"
    exists = os.path.exists(curr_file)
    if exists:
        Popen("rm -rf " + curr_file)
    for chr_num in range(chr_start,chr_end + 1):
        one_vcf = out_dir + "Aquila_SNPs_chr" + str(chr_num) + ".vcf"
        cat_cmd = "cat " + one_vcf  + " >> " + out_dir + "Aquila_SNPs.vcf"
        Popen(cat_cmd,shell=True).wait()
    print("all done~")


def Call_SV_del_from_contigs(chr_start,chr_end,out_dir,num_of_threads):
    total_num = chr_end - chr_start +1
    pool = Pool(num_of_threads)
    count = 1
    for chr_num in range(chr_start,chr_end + 1):
        all_cmd = "python3 " + code_path + "Extract_SV_info_from_contigs_use_overlap_for_del_forcontiginfo.py --chr_start " + str(chr_num) + " --chr_end " + str(chr_num) + " --out_dir " + out_dir 
        count += 1
        pool.apply_async(Run_del_or_ins,(all_cmd,"xin"))
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()

            if (count - 1) == total_num:
                print("finished chr" + str(chr_num))
            else:
                pool = Pool(num_of_threads)
    curr_file = out_dir + "Aquila_DEL.vcf"
    exists = os.path.exists(curr_file)
    if exists:
        Popen("rm -rf " + curr_file)
    for chr_num in range(chr_start,chr_end + 1):
        one_vcf = out_dir + "Aquila_DEL_chr" + str(chr_num) + ".vcf"
        cat_cmd = "cat " + one_vcf  + " >> " + out_dir + "Aquila_DEL.vcf"
        Popen(cat_cmd,shell=True).wait()

    print("all done~")


def Call_SV_ins_from_contigs(chr_start,chr_end,out_dir,num_of_threads):
    total_num = chr_end - chr_start +1
    pool = Pool(num_of_threads)
    count = 1
    for chr_num in range(chr_start,chr_end + 1):
        all_cmd = "python3 " + code_path + "Extract_SV_info_from_contigs_use_shift_for_ins2_forcontiginfo.py --chr_start " + str(chr_num) + " --chr_end " + str(chr_num) + " --out_dir " + out_dir 
        count += 1
        pool.apply_async(Run_del_or_ins,(all_cmd,"xin"))
        if (count - 1)%num_of_threads == 0 or (count - 1) == total_num:
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()

            if (count - 1) == total_num:
                print("finished chr" + str(chr_num))
            else:
                pool = Pool(num_of_threads)
    curr_file = out_dir + "Aquila_INS.vcf"
    exists = os.path.exists(curr_file)
    if exists:
        Popen("rm -rf " + curr_file)
    for chr_num in range(chr_start,chr_end + 1):
        one_vcf = out_dir + "Aquila_INS_chr" + str(chr_num) + ".vcf"
        cat_cmd = "cat " + one_vcf  + " >> " + out_dir + "Aquila_INS.vcf"
        Popen(cat_cmd,shell=True).wait()

    print("all done~")



def Merge_all_variants(out_dir):
    snp_vcf = out_dir + "Aquila_SNPs.vcf"
    del_vcf = out_dir + "Aquila_DEL.vcf"
    ins_vcf = out_dir + "Aquila_INS.vcf"
    final_vcf = out_dir + "Aquila_final.vcf"
    final_sort_vcf = out_dir + "Aquila_final_sorted.vcf"
    """
    header_vcf = out_dir + "header.vcf"
    fw = open(out_dir + "header.vcf","w")
    fw.writelines("##fileformat=VCFv4.2\n")
    fw.writelines("##source=Aquila\n")
    fw.writelines("##contig=<ID=chr1,length=248956422>\n")
    fw.writelines("##contig=<ID=chr2,length=242193529>\n")
    fw.writelines("##contig=<ID=chr3,length=198295559>\n")
    fw.writelines("##contig=<ID=chr4,length=190214555>\n")
    fw.writelines("##contig=<ID=chr5,length=181538259>\n")
    fw.writelines("##contig=<ID=chr6,length=170805979>\n")
    fw.writelines("##contig=<ID=chr7,length=159345973>\n")
    fw.writelines("##contig=<ID=chr8,length=145138636>\n")
    fw.writelines("##contig=<ID=chr9,length=138394717>\n")
    fw.writelines("##contig=<ID=chr10,length=133797422>\n")
    fw.writelines("##contig=<ID=chr11,length=135086622>\n")
    fw.writelines("##contig=<ID=chr12,length=133275309>\n")
    fw.writelines("##contig=<ID=chr13,length=114364328>\n")
    fw.writelines("##contig=<ID=chr14,length=107043718>\n")
    fw.writelines("##contig=<ID=chr15,length=101991189>\n")
    fw.writelines("##contig=<ID=chr16,length=90338345>\n")
    fw.writelines("##contig=<ID=chr17,length=83257441>\n")
    fw.writelines("##contig=<ID=chr18,length=80373285>\n")
    fw.writelines("##contig=<ID=chr19,length=58617616>\n")
    fw.writelines("##contig=<ID=chr20,length=64444167>\n")
    fw.writelines("##contig=<ID=chr21,length=46709983>\n")
    fw.writelines("##contig=<ID=chr22,length=50818468>\n")
    fw.writelines("##contig=<ID=chrX,length=156040895>\n")
    fw.writelines("##contig=<ID=chrY,length=57227415>\n")
    fw.writelines("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV:DEL=Deletion, INS=Insertion, SNP=snps\">\n")
    fw.writelines("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    fw.writelines("##FORMAT=<ID=Contig,Number=1,Type=String,Description=\"Contig information: Contig name, contig start, contig end\">\n")
    fw.writelines("#CHROM" + "\t" + "POS" + "\t" + "ID" + "\t" + "REF" + "\t" + "ALT" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT" + "\t"  + sample_name + "\n")
    fw.close()
    """
    ####use_cmd = "cat " + header_vcf + " " + snp_vcf + " "  + del_vcf + " " + ins_vcf + " > " + final_vcf 
    use_cmd = "cat " + snp_vcf + " "  + del_vcf + " " + ins_vcf + " > " + final_vcf 
    sort_cmd = "cat " + final_vcf + " | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"LC_ALL=C sort -k1,1 -k2,2n\"}' > " + final_sort_vcf
    Popen(use_cmd,shell=True).wait()
    Popen(sort_cmd,shell=True).wait()
    print("all done~")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_assembly_based_variants_call.py -h",shell=True).wait()
    else:
        out_dir = args.out_dir + "/"
        if os.path.exists(out_dir):
            print("using existing output folder: " + out_dir)
        else:
            os.makedirs(out_dir)
        chr_start = args.chr_start
        chr_end = args.chr_end
        num_of_threads = args.num_of_threads
        in_dir = args.assembly_dir + "/"
        ref_file = args.ref_file
        Split_haplotype(chr_start,chr_end,in_dir)
        assembly_based_variants_call_paf(chr_start,chr_end,ref_file)
        assembly_based_variants_call_sort(chr_start,chr_end)
        assembly_based_variants_call(chr_start,chr_end)
        Call_SNV_info_from_contigs(chr_start,chr_end,out_dir,num_of_threads)
        Call_SV_del_from_contigs(chr_start,chr_end,out_dir,num_of_threads)
        Call_SV_ins_from_contigs(chr_start,chr_end,out_dir,num_of_threads)
        Merge_all_variants(out_dir)



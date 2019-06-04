import pdb
#pdb.set_trace()
from collections import defaultdict
from argparse import ArgumentParser
import pickle
from subprocess import Popen
import sys
import os
parser = ArgumentParser(description="Get Uniqness_map folder to run Aquila:")
parser.add_argument('--umap_bed_gz_file','-i',help="Required parameter, k100.umap.bed.gz by hoffmanMappability, k100 is recommended, you could also use k24,k36,k50",required=True)
parser.add_argument('--chr_start','-start',type=int,help="chromosome number start from, default = 1 ", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome number end by, default = 23 ", default=23)
parser.add_argument('--out_dir','-o_dir', help="Directory to store uniqness map files, default = ./Uniqness_map", default='Uniqness_map')

args = parser.parse_args()


def Get_uniqness(input_file,output_file):
    uniq_map = defaultdict(int)
    f = open(input_file,"r")
    for line in f:
        data = line.rsplit()
        _start = int(data[1])
        _end = int(data[2])
        block_len = _end - _start
        if block_len >= 500:
            use_start = _start + 10
            use_end = _end - 10
            for step in range(use_start, use_end+1):
                uniq_map[step] = 1
    pickle.dump(uniq_map,open(output_file,"wb"))




if __name__ == "__main__":
    if len(sys.argv) == 1:
        Popen("python3 " + "Aquila_sortbam.py -h",shell=True).wait()
    else:
        umap_bed_gz = args.umap_bed_gz_file
        output_file_prefix = "uniq_map"
        out_dir = args.out_dir + "/"
        chr_start = args.chr_start
        chr_end = args.chr_end
        if os.path.exists(out_dir):
            print("using existing output folder: " + out_dir)
        else:
            os.makedirs(out_dir)
        for chr_num in range(chr_start,chr_end + 1):
            uniq_file = out_dir +  "chr" + str(chr_num) + "_uniq.bed"
            if chr_num == 23:
                run_cmd = "zcat " + umap_bed_gz + "| grep -w \"chrX\""  + " > " + uniq_file
            else:
                run_cmd = "zcat " + umap_bed_gz + "| grep -w \"chr\"" + str(chr_num) + " > " + uniq_file
            Popen(run_cmd,shell=True).wait()
            output_file = out_dir + output_file_prefix + "_chr" + str(chr_num) + ".p"
            Get_uniqness(uniq_file,output_file)



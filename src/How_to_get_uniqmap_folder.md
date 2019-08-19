The below pipeline to generate K100 Umap is implemented based on <a href="https://bismap.hoffmanlab.org/">hoffmanMappability</a>.  
Credit by Yichen Liu (liuyichen@std.uestc.edu.cn)
1. Use Jupyter notebook "Aquila_Umap.ipynb"    :octocat: <br />
Configure the second block in Aquila_Umap.ipynb to meet your requirements, then run the blocks one by one.
```
#============================================================================================
fa_folder = "path/to/fafolder/"     #The dir of the folder which contains your fasta files
fa_name = "sample.fa"               #The fasta file you want to process
out_dir = "out/put/dir/"            #The output folder which Aquila_Umap will output to
start = 1                           #The start chromosome
end = 23                            #The end chromosome.
kmer_len = 100                      #Length of the kmers generated
mapq_thres = 20                     #MAPQ filter threshold (20 recommended)
bowtie_thread = 20                  #The number of parallel search threads bowtie2 would use 
#============================================================================================
```
Before running Aquila_Umap, make sure you have already installed <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a> and other Aquila dependencies. 

NOTICE: Some of the files (especially fastq and sam files) generated during the process would take up a lot of storage (approximately 5~80G each, depending on the size of the chromosome). Although the script will remove them, be sure your disk has enough space remained. 300G avaliable storage is recommended for running on 23 Human chromosomes.

For testing, give start and end same chromosome number to only process one chromosome (Human chromosome 21 is recommended).

2. Use "Aquila_Umap.py" :octocat: <br />

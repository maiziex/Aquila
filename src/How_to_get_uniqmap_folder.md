The below pipeline to generate K100 Umap is implemented based on <a href="https://bismap.hoffmanlab.org/">hoffmanMappability</a>.  
Credit by Yichen Liu (liuyichen@std.uestc.edu.cn)

### 1. Use Jupyter notebook "<a href="https://github.com/maiziex/Aquila/blob/master/bin/Aquila_Umap.ipynb">Aquila/bin/Aquila_Umap.ipynb</a>"    :octocat: <br />

Configure the second block in Aquila_Umap.ipynb to meet your requirements, then run the whole notebook.
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
Before running Aquila_Umap, make sure you have already installed <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a> and other Aquila dependencies (samtools, python3(+pysam)). 


#### A practical example

Download Rhesus macaque (Macaca mulatta, <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a>) genome fasta file from  <a href="http://xinzhouneuroscience.org/wp-content/uploads/2019/08/macaca_mulatta.fa">HERE</a>.

Your folder structure should be as follows :
```
Rhesus_macaque
    |-fasta
        |-macaca_mulatta.fa
```
"macaca_mulatta.fa" is the fasta file you have just downloaded.

Then edit the second block and run the whole notebook.
```
fa_folder = "path/to/Rhesus_macaque/fasta/"
fa_name = "macaca_mulatta.fa"
out_dir = "path/to/Rhesus_macaque/output/"
start = (your start chromosome)
end = (your end chromosome)
kmer_len = 100 
mapq_thres = 20
bowtie_thread = 20 
```

When finished, you will get: 
```
Rhesus_macaque
|-fasta
|   |-macaca_mulatta.fa
|
|-output
    |-Uniqness_map
        |-uniq_map_chrxxx.p
        |-uniq_map_chrxxx.p
        |-uniq_map_chrxxx.p
        ...
```
The final results are stored in Uniqness_map folder.

### 2. Use "Aquila_Umap.py" :octocat: <br />

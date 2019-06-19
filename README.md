# :milky_way: Aquila :eagle: 


## Download Aquila v1.1:
```
git clone https://github.com/maiziex/Aquila.git
```

## Dependencies:
Aquila utilizes <a href="https://www.python.org/downloads/">Python3</a>, pysam, <a href="http://samtools.sourceforge.net/">SAMtools</a>, and <a href="https://github.com/lh3/minimap2">minimap2</a>. To be able to execute the above programs by typing their name on the command line, the program executables must be in one of the directories listed in the PATH environment variable (".bashrc"). <br />
Or you could just run "./install.sh" to check their availability and install them if not, but make sure you have installed "python3", "conda" and "wget" first. 

## Install:
```
cd Aquila
chmod +x install.sh
./install.sh
```

## source folder:
After running "./install.sh", a folder "source" would be download, it includes human GRCh38 reference fasta file, or you could also just download it by yourself from the corresponding official websites. 

## Running The Code:
Put the "Aquila/bin" in the ".bashrc" file, and source the ".bashrc" file <br />
Or just use the fullpath of "**Aquila_step1.py**" and "**Aquila_step2.py**"


### Step 1: 
```
Aquila/bin/Aquila_step1.py --bam_file possorted_bam.bam --vcf_file S12878_freebayes.vcf --sample_name S12878 --out_dir Assembly_results_S12878 --uniq_map_dir Aquila/Uniqness_map
```
#### *Required parameters
##### --bam_file: "possorted_bam.bam" is bam file generated from barcode-awere aligner like "Lonranger align". How to get bam file, you can also check <a href="https://github.com/maiziex/Aquila/blob/master/src/How_to_get_bam_and_vcf.md">here</a>.

##### --vcf_file: "S12878_freebayes.vcf" is VCF file generated from variant caller like "FreeBayes". How to get vcf file, you can also check <a href="https://github.com/maiziex/Aquila/blob/master/src/How_to_get_bam_and_vcf.md">here</a>. 

#####  --sample_name: "S12878" are the sample name you can define. 

#####  --uniq_map_dir: "Aquila/Uniqness_map" is the uniqness file for GRCh38 you can download by "./install.sh".

#### *Optional parameters
##### --mbq_threshold, default = 13, It's phred-scaled quality score for the assertion made in ALT.
##### --boundary: default = 50000 (50kb), It is the boundary for long fragments with the same barcode. 
##### --out_dir, default = ./Asssembly_results. You can define your own folder, for example "Assembly_results_S12878". 

##### --block_threshold, default = 200000 (200kb)
 
##### --block_len_use, default = 100000 (100kb)

##### --num_threads, default = 8. It's recommended not to change this setting unless large memory node could be used (2*memory capacity(it suggests for assembly below)), then try to use "--num_threads 12". 

##### --num_threads_for_bwa_mem, default = 20. This setting is evoked for "bwa mem".

##### --chr_start, --chr_end: if you only want to assembly some chromosomes or only one chromosome. For example: use "--chr_start 1 --chr_end 5"  will assemble chromsomes 1,2,3,4,5. Use "--chr_start 2 --chr_end 2" will only assemlby chromosome 2. 

To use the above option "--chr_start, --chr_end", it is recommended to run the below command first to save more time later. 
```
python Aquila/bin/Aquila_step0_sortbam.py --bam_file possorted_bam.bam --out_dir Results_S12878 --num_threads_for_bwa_mem 20 
```
<!--   -->
#### Memory/Time Usage For Step 1
##### Running Step 1 for chromosomes parallelly on multiple(23) nodes

Coverage | Memory| Time for chr1 on a single node | 
--- | --- | --- | 
60X | 100GB |   13:12:35 |
90X | 150GB | 1-01:08:38 |

##### Running Step 1 for WGS on a single node with large memory
Coverage | Memory| Time for WGS on a single node  | 
--- | --- | --- | 
60X | 450GB | 2-14:27:37 |
90X | 500GB | 3-12:34:12 |





### Step 2: 
```
Aquila/bin/Aquila_step2.py --out_dir Assembly_results_S12878 --num_threads 30 --reference Aquila/source/ref.fa
```
#### *Required parameters
#####  --reference: "Aquila/source/ref.fa" is the reference fasta file you can download by "./install".

#### *Optional parameters
#####  --out_dir, default = ./Asssembly_results, make sure it's the same as "--out_dir" from ***Step1*** if you want to define your own output directory name.

#####  --num_threads, default = 30, this determines the number of files assembled simultaneously by SPAdes.  

#####  --num_threads_spades, default = 5, this is the "-t" for SPAdes. 

##### --block_len_use, default = 100000 (100kb)

##### --chr_start, --chr_end: if you only want to assembly some chromosomes or only one chromosome. For example: use "--chr_start 1 --chr_end 2" 

<!-- -->
#### Memory/Time Usage For Step 2
##### Running Step 2 for chromosomes parallelly on multiple nodes
Coverage| Memory| Time for chr1 on a single node | --num_threads | --num_threads_spades|
--- | --- | --- | ---|---|
60X| 100GB | 09:50:43 |30 | 20|
90X| 100GB | 13:25:08 |40 | 20|

##### Running Step 2 for WGS on a single node with large memory
Coverage| Memory| Time for WGS on a single node  | --num_threads | --num_threads_spades|
 ---| --- | --- | ---|---|
60X| 100GB | 3-12:16:27 |40 | 20|
90X| 100GB | 4-15:00:00 |40 | 20|
 




### Clean Data
##### If your hard drive storage is limited, it is suggested to quily clean some data by running "Aquila_clean.py". Or you can keep them for some analysis. 
```
Aquila/bin/Aquila_clean.py --out_dir Assembly_results_S12878 
```

## Final Output:
##### Assembly_Results_S12878/Assembly_Contigs_files: Aquila_contig.fasta and Aquila_Contig_chr*.fasta 

## Final Output Format:
Aquila outputs an overall contig file “Aquila_Contig_chr*.fasta” for each chromosome, and one contig file for each haplotype: “Aquila_Contig_chr*_hp1.fasta” and “Aquila_Contig_chr*_hp2.fasta”. For each contig, the header, for an instance, “>36_PS39049620:39149620_hp1” includes contig number “36”, phase block start coordinate “39049620”, phase block end coordinate “39149620”, and haplotype number “1”. Within the same phase block, the haplotype number “hp1” and “hp2” are arbitrary for maternal and paternal haplotypes. For some contigs from large phase blocks, the headers are much longer and complex, for an instance, “>56432_PS176969599:181582362_hp1_ merge177969599:178064599_hp1-177869599:177969599_hp1”. “56” denotes contig number, “176969599” denotes the start coordinate of the final big phase block, “181582362” denotes the end coordinate of the final big phase block, and “hp1” denotes the haplotype “1”. “177969599:178064599_hp1” and “177869599:177969599_hp1” mean that this contig is concatenated from minicontigs in small chunk (start coordinate: 177969599, end coordinate: 178064599, and haplotype: 1) and small chunk (start coordinate: 177869599, end coordinate: 177969599, and haplotype: 1). 

## Assembly Based Variants Calling and Phasing:
##### For example, you can use "Assemlby_results_S12878" as input directory to generate a VCF file which includes SNPs, small Indels and SVs, and the phased profile of all of them. 
##### Please check check <a href="https://github.com/maiziex/Aquila/blob/master/Assembly_based_variants_call/README.md/">Assembly_based_variants_call_and_phasing</a> for details. 

## Aquila assembly for other version of human referece: hg19
##### 1. Download hg19 reference from <a href="https://support.10xgenomics.com/genome-exome/software/downloads/latest">10x Genomics website</a>
##### 2. Download hg19 "Uniqness_map" folder by wget using the link
```
wget http://xinzhouneuroscience.org/wp-content/uploads/2019/06/Uniqness_map_hg19.tar.gz 
```
##### If you want to run Aquila for other diploid species, to generate "Uniqness_map" for Aquila, check the details of  <a href="https://bismap.hoffmanlab.org/">hoffmanMappability</a> to get the corresponding "k100.umap.bed.gz", then run "Aquila/bin/Get_uniqnessmap_for_Aquila.py" to get the final "Uniqness_map" folder to run Aquila.


# Assembly for multiple libraries:

### Step 1: 
```
Aquila/bin/Aquila_step1_multilibs.py --bam_file_list ./S24385_Lysis_2/Longranger_align_bam/S24385_lysis_2/outs/possorted_bam.bam,./S24385_Lysis_2H/Longranger_align_bam/S24385_lysis_2H/outs/possorted_bam.bam --vcf_file_list ./S24385_lysis_2/Freebayes_results/S24385_lysis_2_grch38_ref_freebayes.vcf,./S24385_lysis_2H/Freebayes_results/S24385_lysis_2H_grch38_ref_freebayes.vcf --sample_name_list S24385_lysis_2,S24385_lysis_2H --out_dir Assembly_results_merged --uniq_map_dir Aquila/Uniqness_map
```
#### *Required parameters
##### --bam_file: "possorted_bam.bam" is bam file generated from barcode-awere aligner like "Lonranger align". Each bam file is seperately by comma (",").

##### --vcf_file: "S12878_freebayes.vcf" is VCF file generated from variant caller like "FreeBayes". Each VCF file is seperately by comma (",").


#####  --sample_name: S24385_lysis_2,S24385_lysis_2H are the sample names you can define. Each sample name is seperately by comma (",").

#####  --uniq_map_dir: "Aquila/Uniqness_map" is the uniqness file you can download by "./install.sh".

#### *Optional parameters
#####  --out_dir, default = ./Asssembly_results 

##### --block_threshold, default = 200000 (200kb)
 
##### --block_len_use, default = 100000 (100kb)

##### --num_threads, default = 8. It's recommended not to change this setting unless large memory node could be used (2*memory capacity(it suggests for assembly below)), then try to use "--num_threads 12". 

##### --num_threads_for_bwa_mem, default = 20. This setting is evoked for "bwa mem".


##### --chr_start, --chr_end: if you only want to assembly some chromosomes or only one chromosome. For example: use "--chr_start 1 --chr_end 5"  will assemble chromsomes 1,2,3,4,5. Use "--chr_start 2 --chr_end 2" will only assemlby chromosome 2. 
To use the above option "--chr_start, --chr_end", it is recommended to run the below command first to save more time later. 
```
python Aquila/bin/Aquila_step0_sortbam_hybrid.py --bam_file_list ./S24385_Lysis_2/Longranger_align_bam/S24385_lysis_2/outs/possorted_bam.bam,./S24385_Lysis_2H/Longranger_align_bam/S24385_lysis_2H/outs/possorted_bam.bam --out_dir Assembly_results_merged --num_threads_for_bwa_mem 10 --sample_name_list S24385_lysis_2,S24385_lysis_2H 
```

### Step 2: (The same as single library assembly)
```
Aquila/bin/Aquila_step2.py --out_dir Assembly_results_merged --num_threads 30 --reference Aquila/source/ref.fa
```
#### *Required parameters
#####  --reference: "Aquila/source/ref.fa" is the reference fasta file you can download by "./install".

#### *Optional parameters
#####  --out_dir, default = ./Asssembly_results, make sure it's the same as "--out_dir" from step1 if you want to define your own output directory name.

#####  --num_threads, default = 20 

##### --block_len_use, default = 100000 (100kb)

##### --chr_start, --chr_end: if you only want to assembly some chromosomes or only one chromosome. 


## Cite Aquila:
#### Aquila: diploid personal genome assembly based on linked-reads enabling comprehensive variant detection.
##### <a href="https://www.biorxiv.org/content/10.1101/660605v1">bioRxiv link</a>

## Troubleshooting:
##### Please submit issues on the github page for <a href="https://github.com/maiziex/Aquila/issues">Aquila</a>.


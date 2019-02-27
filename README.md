# XinMagic


## Download:
```
git clone https://github.com/maiziex/XinMagic.git
```

## Dependencies:
XinMagic utilizes <a href="https://www.python.org/downloads/">Python3</a>, <a href="http://samtools.sourceforge.net/">SAMtools</a>, and <a href="https://github.com/lh3/minimap2">minimap2</a>. To be able to execute the above programs by typing their name on the command line, the program executables must be in one of the directories listed in the PATH environment variable (".bashrc"). <br />
Or you could just run "./install.sh" to check their availability and install them if not, but make sure you have installed "python3", "conda" and "wget" first. 

## Install:
```
cd XinMagic
chmod +x install.sh
./install.sh
```

## source folder:
After running "./install.sh", a folder "source" would be download, it includes human GRCh38 reference fasta file, or you could also just download it by yourself from the corresponding official websites. 

## Running The Code:
Put the "XinMagic/bin" in the ".bashrc" file, and source the ".bashrc" file <br />
Or just use the fullpath of "XinMagic_Step1.py" and "XinMagic_Step2.py"


### Step 1: 
```
XinMagic/XinMagic_step1.py --bam_file possorted_bam.bam --vcf_file S12878_freebayes.vcf --sample_name S12878 --out_dir Results_S12878 --uniq_map_dir XinMagic/Uniqness_map
```
#### *Required parameters
##### --bam_file: "possorted_bam.bam" is bam file generated from barcode-awere aligner like "Lonranger align". How to get bam file, you can also check <a href="https://github.com/maiziex/XinMagic/blob/master/src/How_to_get_bam_and_vcf.md">here</a>.

##### --vcf_file: "S12878_freebayes.vcf" is VCF file generated from variant caller like "FreeBayes". How to get vcf file, you can also check <a href="https://github.com/maiziex/XinMagic/blob/master/src/How_to_get_bam_and_vcf.md">here</a>. 

#####  --sample_name: "S12878" are the sample name you can define. 

#####  --uniq_map_dir: "XinMagic/Uniqness_map" is the uniqness file you can download by "./install.sh".

#### *Optional parameters
#####  --out_dir, default = ./Asssembly_results 

##### --block_threshold, default = 200000 (200kb)
 
##### --block_len_use, default = 100000 (100kb)

##### --chr_start, --chr_end: if you only want to assembly some chromosomes or only one chromosome. For example: use "--chr_start 1 --chr_end 5"  will assemble chromsomes 1,2,3,4,5. Use "--chr_start 2 --chr_end 2" will only assemlby chromosome 2. 

To use the above option "--chr_start, --chr_end", it is recommended to run the below command first to save more time later. 
```
python XinMagic/XinMagic_sortbam.py --bam_file possorted_bam.bam --out_dir Results_S12878 --num_threads 20 
```

### Step 2: 
```
XinMagic/XinMagic_step2.py --out_dir Results_S12878 --num_threads 30 --reference XinMagic/source/ref.fa
```
#### *Required parameters
#####  --reference: "XinMagic/source/ref.fa" is the reference fasta file you can download by "./install".

#### *Optional parameters
#####  --out_dir, default = ./Asssembly_results, make sure it's the same as "--out_dir" from step1 if you want to define your own output directory name.

#####  --num_threads, default = 20 

##### --block_len_use, default = 100000 (100kb)

##### --chr_start, --chr_end: if you only want to assembly some chromosomes or only one chromosome. For example: use "--chr_start 1 --chr_end 2" 


## Final Output:
##### Results/Assembly_Contigs_files: XinMagic_contig.fasta and XinMagic_Contig_chr*.fasta 
##### Results/H5_for_molecules: 
##### Results/results_phased_probmodel:



## Final Output Format:
XinMagic outputs an overall contig file “XinMagic_Contig_chr*.fasta” for each chromosome, and one contig file for each haplotype: “XinMagic_Contig_chr*_hp1.fasta” and “XinMagic_Contig_chr*_hp2.fasta”. For each contig, the header, for an instance, “>36_PS39049620:39149620_hp1” includes contig number “36”, phase block start coordinate “39049620”, phase block end coordinate “39149620”, and haplotype number “1”. Within the same phase block, the haplotype number “hp1” and “hp2” are arbitrary for maternal and paternal haplotypes. For some contigs from large phase blocks, the headers are much longer and complex, for an instance, “>56432_PS176969599:181582362_hp1_ merge177969599:178064599_hp1-177869599:177969599_hp1”. “56” denotes contig number, “176969599” denotes the start coordinate of the final big phase block, “181582362” denotes the end coordinate of the final big phase block, and “hp1” denotes the haplotype “1”. “177969599:178064599_hp1” and “177869599:177969599_hp1” mean that this contig is concatenated from minicontigs in small chunk (start coordinate: 177969599, end coordinate: 178064599, and haplotype: 1) and small chunk (start coordinate: 177869599, end coordinate: 177969599, and haplotype: 1). 

## Assembly Based Variants Calling:
##### You can use "Results/Assembly_Contigs_files" as input directory to generate a VCF file which includes SNPs and SVs. 
##### Please check check <a href="https://github.com/maiziex/XinMagic/blob/master/Assembly_based_variants_call/README.md/">Assembly_based_variants_call</a> for detail. 


# Assembly for multiple libraries:

### Step 1: 
```
XinMagic/XinMagic_step1_hybrid.py --bam_file ./S24385_Lysis_2/Longranger_align_bam/S24385_lysis_2/outs/possorted_bam.bam,./S24385_Lysis_2H/Longranger_align_bam/S24385_lysis_2H/outs/possorted_bam.bam --vcf_file ./S24385_lysis_2/Freebayes_results/S24385_lysis_2_grch38_ref_freebayes.vcf,./S24385_lysis_2H/Freebayes_results/S24385_lysis_2H_grch38_ref_freebayes.vcf --sample_name S24385_lysis_2,S24385_lysis_2H --out_dir Results --uniq_map_dir XinMagic/Uniqness_map
```
#### *Required parameters
##### --bam_file: "possorted_bam.bam" is bam file generated from barcode-awere aligner like "Lonranger align". Each bam file is seperately by comma (",").

##### --vcf_file: "S12878_freebayes.vcf" is VCF file generated from variant caller like "FreeBayes". Each VCF file is seperately by comma (",").


#####  --sample_name: S24385_lysis_2,S24385_lysis_2H are the sample names you can define. Each sample name is seperately by comma (",").

#####  --uniq_map_dir: "XinMagic/Uniqness_map" is the uniqness file you can download by "./install.sh".

#### *Optional parameters
#####  --out_dir, default = ./Asssembly_results 

##### --block_threshold, default = 200000 (200kb)
 
##### --block_len_use, default = 100000 (100kb)

##### --chr_start, --chr_end: if you only want to assembly some chromosomes or only one chromosome. For example: use "--chr_start 1 --chr_end 2" 

### Step 2: (The same as single library assembly)
```
XinMagic/XinMagic_step2.py --out_dir Results --num_threads 30 --reference XinMagic/source/ref.fa
```
#### *Required parameters
#####  --reference: "XinMagic/source/ref.fa" is the reference fasta file you can download by "./install".

#### *Optional parameters
#####  --out_dir, default = ./Asssembly_results, make sure it's the same as "--out_dir" from step1 if you want to define your own output directory name.

#####  --num_threads, default = 20 

##### --block_len_use, default = 100000 (100kb)

##### --chr_start, --chr_end: if you only want to assembly some chromosomes or only one chromosome. For example: use "--chr_start 1 --chr_end 5"  will assemble chromsomes 1,2,3,4,5. Use "--chr_start 2 --chr_end 2" will only assemlby chromosome 2. 
To use the above option "--chr_start, --chr_end", it is recommended to run the below command first to save more time later. 
```
python XinMagic/XinMagic_sortbam.py --bam_file possorted_bam.bam --out_dir Results_S12878 --num_threads 20 
```


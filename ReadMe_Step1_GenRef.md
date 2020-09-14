
### Step 1: 
```
Aquila/bin/Aquila_step1_GenRef.py --bam_file possorted_bam.bam --vcf_file ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  --sample_name S12878 --out_dir Assembly_results_S12878 --uniq_map_dir Aquila/Uniqness_map_hg19 --chr_start 21 --chr_end 21
```
#### *Required parameters
##### --bam_file: "possorted_bam.bam" is bam file generated from barcode-awere aligner like "Lonranger align". How to get bam file, you can also check <a href="https://github.com/maiziex/Aquila/blob/master/src/How_to_get_bam_and_vcf.md">here</a>.

##### --vcf_file: "ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" is VCF file of chr21 from 1000 Genomes Phase3. 

#### *Optional parameters are the same as "Aquila_step1.py"

### Step 2: use "Aquila_step2.py" 

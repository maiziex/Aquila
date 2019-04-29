# Assembly-based variants calling 


## Running The Code:
Put the "Aquila/bin" in the ".bashrc" file, and source the ".bashrc" file <br />
Or use the fullpath of "Aquila_assembly_based_variants_call.py" and "Aquila_phasing_all_variants.py"


```
Aquila/Aquilla_assembly_based_variants_call.py --in_dir Assembly_results_S12878/Assembly_Contigs_files --out_dir Variants_results_S12878 --ref_file refdata-GRCh38-2.1.0/fasta/genome.fa --sample_name S12878
```
#### *Required parameters
##### --in_dir: "Assembly_results_S12878/Assembly_Contigs_files" is the output folder from assembly.

##### --ref_file: "refdata-GRCh38-2.1.0/fasta/genome.fa" is the reference fasta file from GRCh38.

#####  --sample_name: "S12878" are the sample name you can define. 


#### *Optional parameters
#####  --out_dir, default = ./Variants_results

##### --chr_start, --chr_end: if you only want to call variants from some chromosomes or only one chromosome. For example: use "--chr_start 1 --chr_end 5"  will call variants from chromsomes 1,2,3,4,5. Use "--chr_start 2 --chr_end 2" will only call variants from chromosome 2. 

# Phasing all types of variants  

```
Aquila/Aquila_phasing_all_variants.py --assembled_vcf Variants_results_S12878/Aquila_final_sorted.vcf.gz --vcf_file S12878_freebayes.vcf --assemble_dir Assembly_results_S12878
```
#### *Required parameters
##### --assemble_dir: "Assembly_results_S12878" is the output folder from assembly.



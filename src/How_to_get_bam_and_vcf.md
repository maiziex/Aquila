Inputs is 10X `fastqs` files, to get `bam` file: <br />

To check how to download and use "longranger align" or download human reference genome, go to <a href="https://support.10xgenomics.com/genome-exome/software/downloads/latest">10X GENOMICS Website</a>.
```
longranger align --id=S12878_lysis --reference=./refdata-GRCh38-2.1.0 --fastqs=./S12878_Lysis/Final_fastqs
```



Input is `bam` file, to get `vcf` file" : <br />
"possorted_bam.bam" is generated from the above "longranger align"

```
freebayes -f ./refdata-GRCh38-2.1.0/fasta/genome.fa possorted_bam.bam > S12878_grch38_ref_freebayes.vcf 
```

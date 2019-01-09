Inputs are 10X fastqs files, to get `bam` files: <br />

```
longranger align --id=S12878_lysis --reference=./refdata-GRCh38-2.1.0 --fastqs=./S12878_Lysis/Final_fastqs
```

To check how to download and use "longranger align" or download human reference genome, go to <a href="https://support.10xgenomics.com/genome-exome/software/downloads/latest">10X GENOMICS Website</a>.



```
freebayes -f ./refdata-GRCh38-2.1.0/fasta/genome.fa possorted_bam.bam > S12878_grch38_ref_freebayes.vcf 
```
possorted_bam.bam: the output bam file from "longranger align".


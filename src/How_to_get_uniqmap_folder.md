1. Inputs is 10X `fastqs` files, to get `bam` file: :octocat: <br />

To check how to download and use "longranger align" or download human reference genome in details, go to <a href="https://support.10xgenomics.com/genome-exome/software/downloads/latest">10X GENOMICS Website</a>.
```
longranger align --id=S12878_lysis --reference=./refdata-GRCh38-2.1.0 --fastqs=./S12878_Lysis/Final_fastqs
```



2. Input is `bam` file, to get `vcf` file" : :octocat: <br />
"possorted_bam.bam" is generated from the above "longranger align". <br />
To check how to download and use "freebayes" in details, go to <a href="https://github.com/ekg/freebayes">FreeBayes Website</a> 

```
freebayes -f ./refdata-GRCh38-2.1.0/fasta/genome.fa possorted_bam.bam > S12878_grch38_ref_freebayes.vcf 
```

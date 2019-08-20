for i in {1..23} 
do
    cp submit_step1.sbatch temp.sbatch
    echo Aquila/bin/Aquila_step1.py --bam_file possorted_bam.bam --vcf_file S12878_freebayes.vcf --sample_name S12878 --out_dir Assembly_results_S12878 --uniq_map_dir Aquila/Uniqness_map --chr_start "$i" --chr_end "$i" >> temp.sbatch
    #echo Aquila/bin/Aquila_step1_multilibs.py --bam_file_list ./S24385_Lysis_2/Longranger_align_bam/S24385_lysis_2/outs/possorted_bam.bam,./S24385_Lysis_2H/Longranger_align_bam/S24385_lysis_2H/outs/possorted_bam.bam --vcf_file_list ./S24385_lysis_2/Freebayes_results/S24385_lysis_2_grch38_ref_freebayes.vcf,./S24385_lysis_2H/Freebayes_results/S24385_lysis_2H_grch38_ref_freebayes.vcf --sample_name_list S24385_lysis_2,S24385_lysis_2H --out_dir Assembly_results_merged --uniq_map_dir Aquila/Uniqness_map --chr_start "$i" --chr_end "$i" >> temp.sbatch
    sed -e 's/--job-name=SV/--job-name='$i'_Aquila/g ; s/--output=SV.out/--output='$i'.out/g ; s/--error=SV.err/--error='$i'.err/g'`` temp.sbatch > submit_step1_"$i".sbatch 
    rm temp.sbatch
    sbatch submit_step1_"$i".sbatch
done

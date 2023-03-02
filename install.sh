# change permission
cd bin
chmod +x *.py
chmod -R 777 k8-0.2.4
tar -xzf SPAdes-3.13.0-Linux.tar.gz
rm SPAdes-3.13.0-Linux.tar.gz
cd ..


# download the reference file (GRCh38)
wget https://zenodo.org/record/7689958/files/source.tar.gz
tar -xvf source.tar.gz
rm source.tar.gz


wget https://zenodo.org/record/7689958/files/Uniqness_map_hg38.tar.gz
tar -xvf Uniqness_map_hg38.tar.gz
rm Uniqness_map_hg38.tar.gz


if ! [ -x "$(command -v samtools)" ];
then
    echo 'Error: samtools is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda samtools
else
    echo 'using existing samtools...'
fi

if ! [ -x "$(command -v minimap2)" ];
then
    echo 'Error: minimap2 is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda minimap2
else
    echo 'using existing minimap2...'
fi



echo 'You have installed Aquila dependencies and downloaded the source files successfully!'
 




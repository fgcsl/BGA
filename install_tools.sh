#!/bin/bash

qc="mamba create -n qc -c bioconda fastqc multiqc bbmap -y"
mappers="mamba create -n mappers -c bioconda bwa bowtie2 bcftools samtools picard -y"
bam2fastq="mamba create -n bam2fastq -c yuxiang bam2fastq -y"
emboss="mamba create -n emboss -c bioconda emboss -y"
assembler="mamba create -n assembler -c bioconda spades -y"
scaffolder="mamba create -n scaffolder -c conda-forge -c bioconda medusa mummer biopython -y"
mauve="mamba create -n mauve -c bioconda mauve -y"
filler="mamba create -n filler -c bioconda soapdenovo2-gapcloser -y"
shovill="mamba create -n shovill -c bioconda shovill -y"
annotator="mamba create -n annotator -c bioconda prokka -y"
pangenome="mamba create -n pangenome -c bioconda roary -y"
busco="mamba create -n busco -c bioconda busco=5.4* -y"
bgc="mamba create -n bgc -c bioconda antismash -y"
dotplot="mamba create -n dotplot -c bioconda dgenies -y"
qiime2="mamba env create -n qiime2 --file qiime2-2023.5-py38-linux-conda.yml"

#Download qiime2 yml file 

while [ ! -f qiime2-2023.5-py38-linux-conda.yml ]
do
	wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-linux-conda.yml
done 

username=`whoami`

mambapath="/home/$username/mambaforge/envs"

for i in $(ls $mambapath); do echo $i; done > installed_envs.txt 

if [ -f install_chk.txt  ]
then

	mamba_envs=$(cat install_chk.txt)
        grepResult=$(grep -v -f installed_envs.txt install_chk.txt)
        greplen=$(echo "$grepResult" | wc -l)
	
	if [ "$grepResult" == "" ]
	then
		echo "All tools were already installed..."
	fi

	while [[ $greplen -gt 0 ]]; do
        	for result in $grepResult; do
                	for env in $mamba_envs; do
                        	if [[ "$env" == "$result" ]]; then
					echo "Installing $env ....."
                                	eval "$"$(eval echo \$env)
				#	echo $env
        	                fi;
                	done;
        	done;
        
		for i in $(ls $mambapath); do echo $i; done > installed_envs.txt 
		grepResult=$(grep -v -f installed_envs.txt install_chk.txt)	
		greplen=$(grep -v -f installed_envs.txt install_chk.txt| wc -l)	
		
		if [ "$grepResult" == "" ]
        	then
                	echo "Installation Successfully Completed......"
        	fi
	done
else
	echo "install_chk.txt file does not exist"
fi

#!/usr/bin/bash

qc="mamba create -n qc -c bioconda fastqc multiqc bbmap -y"

mappers="mamba create -n mappers -c bioconda bwa bowtie2 bcftools samtools -y"

bam2fastq="mamba create -n bam2fastq -c yuxiang bam2fastq -y"

emboss="mamba create -n emboss -c bioconda emboss -y"

assembler="mamba create -n assembler -c bioconda spades -y"

scaffolder="mamba create -n scaffolder -c conda-forge -c bioconda medusa mummer biopython -y"

mauve="mamba create -n mauve -c bioconda mauve -y"

filler="mamba create -n filler -c bioconda soapdenovo2-gapcloser -y"
 
pilon="mamba create -n pilon -c bioconda pilon -y"

annotator="mamba create -n annotator -c bioconda prokka -y"

pangenome="mamba create -n pangenome -c bioconda roary -y"

#qiime2="wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-linux-conda.yml"
#qiime2="mamba env create -n qiime2 --file qiime2-2023.5-py38-linux-conda.yml"

mambapath=$(locate mambaforge/envs |head -1)

for i in $(ls $mambapath); do echo $i; done > installed_envs.txt 

if [ -f install_chk.txt  ]
then

	mamba_envs=$(cat install_chk.txt)
        grepResult=$(grep -v -f installed_envs.txt install_chk.txt)
        greplen=$(echo "$grepResult" | wc -l)
	
	if [ "$grepResult" == "" ]
	then
		echo "All tools and pacakges are already installed..."
	fi

	while [[ $greplen -gt 0 ]]; do
        	for result in $grepResult; do
                	for env in $mamba_envs; do
                        	if [[ "$env" == "$result" ]]; then
					echo "Installing $env....."
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
	echo "install_chk.txt file not exist"
fi

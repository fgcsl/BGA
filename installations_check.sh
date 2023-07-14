#!/bin/bash

echo '----------------- FastQC ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n qc fastqc --version
echo -e '\n----------------- MultiQC ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n qc multiqc --version
echo -e '\n----------------- BBDuk ----------------------'
echo '*** Expected BBMap version info as output from a java command ***'
mamba run -n qc bbduk.sh -Xmx4g --version
echo -e '\n----------------- BWA ----------------------'
echo '*** Expected multi line output with an ERROR ***'
mamba run -n mappers bwa
echo -e '\n----------------- bowtie2 ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n mappers bowtie2 --version | head -n 1
echo -e '\n----------------- bcftools ----------------------'
echo '*** Expected multi line output with version info ***'
mamba run -n mappers bcftools --version
echo -e '\n----------------- bam2fastq ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n bam2fastq bam2fastq --version
echo -e '\n----------------- EMBOSS ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n emboss seqret --version
echo -e '\n----------------- SPAdes ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n assembler spades.py --version
echo -e '\n----------------- MeDuSa ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n scaffolder medusa -h | head -n 1
echo -e '\n----------------- Mauve ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n mauve progressiveMauve --version
echo -e '\n----------------- GapCloser ----------------------'
echo '=== Expected 1st line - Version: 1.12. ERROR expected ==='
mamba run -n filler GapCloser -h
echo -e '\n----------------- Prokka ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n annotator prokka --version
echo -e '\n----------------- Roary ----------------------'
echo '*** Expected Single line output of version info with a Warning ***'
mamba run -n pangenome roary -w
echo -e '\n----------------- BUSCO ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n busco busco -v
echo -e '\n----------------- CheckM ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n checkm checkm -h | head -n 2
echo -e '\n----------------- QIIME2 ----------------------'
echo '*** Expected two line output of version info ***'
mamba run -n qiime2 qiime --version
echo -e '\n----------------- D-Genies ----------------------'
echo '*** Expected multi line output ***'
mamba run -n dotplot dgenies -h
echo -e '\n----------------- Shovill ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n shovill shovill --version
echo -e '\n----------------- Pilon ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n shovill pilon --version
echo -e '\n----------------- antiSMASH ----------------------'
echo '*** Expected Single line output of version info ***'
mamba run -n bgc antismash --version


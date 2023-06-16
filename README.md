# BGA
Bacterial Genome Analysis Piplene
# Bacterial Genome analysis Pipeline -  DBT-CMI 

### Requirements (QC and Assembly)
Installing Anaconda

Check other lighter and faster versions: [Miniconda](https://docs.conda.io/en/latest/miniconda.html), [Miniforge](https://github.com/conda-forge/miniforge), [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) & [Minimamba](https://mamba.readthedocs.io/en/latest/installation.html)

1. Navigate to Desktop 
    ```
    $ cd Desktop
    ```
2. Make a directory 'bgap' and a sub-directory 'resources' & navigate to it
    ```
    $ mkdir bgap
    $ mkdir bgap/resources
    $ cd bgap/resources
    ```
3. Download Anaconda script
    ```
    $ wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
    ```
4. Change permissions and install
    ```
    $ chmod +x Anaconda3-2022.10-Linux-x86_64.sh
    $ Anaconda3-2022.10-Linux-x86_64.sh
    ```
5. Follow onscreen instructions. Install in the default directory. On restarting the Terminal, (base) prefix will be shown. Check the installed environments
    ```
    $ conda info --envs
    ```
6. Add Bioconda channel
    ```
    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict
    ```
7. Update conda (optional step)
    ```
    $ conda update -n base -c defaults conda
    ```
###  Installing sra-tools

1. Create a new environment sra
    ```
    $ conda create -n sra -y
    ```
2. Deactivate base environment and activate sra environment
    ```
    $ conda deactivate
    $ conda activate sra
    ```
3. Install the latest version of sra-tools from anaconda.org
    ```
    $ conda install -c bioconda sra-tools=3.0
    ```
4. Check the installation
    ```
    $ srapath --version
    ```
### Downloading a dataset of raw reads from NCBI SRA

1. While in the resources directory, download Helicobacter pylori raw reads. Check description of the dataset.
    ```
    $ prefetch SRR22388518 -p
    ```
    > **Note:** If downloaded successfully, the following message will appear
        ```
            2023-03-09T14:56:10 prefetch.3.0.3: Current preference is set to retrieve SRA Normalized Format files with full base quality scores.
            2023-03-09T14:56:11 prefetch.3.0.3: 1) Downloading 'SRR22388518'...
            2023-03-09T14:56:11 prefetch.3.0.3: SRA Normalized Format file is being retrieved, if this is different from your preference, it may be due to current file availability.
            2023-03-09T14:56:11 prefetch.3.0.3:  Downloading via HTTPS...
            |-------------------------------------------------- 100%   
            2023-03-09T15:13:42 prefetch.3.0.3:  HTTPS download succeed
            2023-03-09T15:13:43 prefetch.3.0.3:  'SRR22388518' is valid
            2023-03-09T15:13:43 prefetch.3.0.3: 1) 'SRR22388518' was downloaded successfully
            2023-03-09T15:13:43 prefetch.3.0.3: 'SRR22388518' has 0 unresolved dependencies
        ```
        
2. Convert to FASTQ files and delete the downloaded directory
    ```
    $ fasterq-dump SRR22388518/SRR22388518.sra
    $ rm -r SRR22388518
    ```
3. Move FASTQ files to a new directory and rename.
    ```
    $ mkdir ../reads
    $ mv *fastq ../reads/
    $ cd ../reads
    $ mv SRR22388518_1.fastq a45_R1.fastq 
    $ mv SRR22388518_2.fastq a45_R2.fastq
    ```
4. We now have the required raw reads

### Installing QC packages

1. Deactivate sra environment and create qc environment
    ```
    $ conda deactivate
    $ conda create -n qc -y
    ```          

2. Activate qc and install FastQC, BBDuk and Trimmomatic(Optional) - Always check the versions
    ```
    $ conda activate qc
    $ conda install -c bioconda fastqc trimmomatic bbmap
    ```
3. Check installations
    ```
    $ fastqc
    $ trimmomatic -version
    $ bbduk.sh -h
    ```
    
### Installing Mapping Tools

1. Deactivate qc environment and create mappers environment
    ```
    $ conda deactivate
    $ conda create -n mappers -y
    ```
2. Activate mappers and install bwa, bcftools, and samtools
    ```
    $ conda activate mappers
    $ conda install -c bioconda bwa bcftools samtools
    ```
3. Check installations
    ```
    $ bwa
    $ bcftools -h
    $ samtools --help
    ```           

4. Probable Error samtools: error while loading shared libraries: libncurses.so.5: cannot open shared object file: No such file or directory
Fedora users:
    ```
    $ sudo dnf install ncurses-libs
    $ sudo dnf install ncurses-devel
    $ sudo dnf install ncurses-compat-libs
     ```             
    Ubuntu users:
    ```
    $ sudo apt-get install libncurses5
    $ sudo apt-get install ia32-libs
    ```           

5. Install [bam2fastq](). It requires an older version of samtools. Create new environment
    ```
    $ conda deactivate
    $ conda create -n bam2fastq -y
    $ conda activate bam2fastq
    $ conda install -c yuxiang bam2fastq
    ```
6. Install [emboss]() in a new environment
    ```
    $ conda deactivate
    $ conda create -n emboss -y
    $ conda activate emboss
    $ conda install -c bioconda emboss
    ```
### Download Reference Genome from NCBI for Mapping
[Helicobacter pylori A45](https://www.ncbi.nlm.nih.gov/genome/169?genome_assembly_id=901025)
### Installing Assembling Tools
1. Deactivate mappers environment and create assemblers environment
    ```
    $ conda deactivate
    $ conda create -n assemblers -y
    ```            
2. Activate assemblers. Install SPAdes & Barrnap
    ```
    $ conda activate assemblers
    $ conda install -c bioconda spades=3.15* barrnap
    ```
3. Check installations
    ```
    $ spades.py -h
    $ barrnap --help
    ```

######################## QC ################################
    
# Raw Reads QC
## (Bacterial Genome Analysis Piplene)
### K-mers & De Bruijn Graphs

##### k-mers are substrings of length k. [Example](https://en.wikipedia.org/wiki/K-mer#Introduction)
##### [Why k-mers are required]()

### De Bruijn graph a directed graph representing overlaps between sequences of symbols
[Example]()

### FastQC

A Quality Control Tool. Works on FASTQ, SAM and BAM files

1. Navigate to bgap directory and activate qc
    ```
    $ cd Desktop/bgap
    $ conda deactivate
    $ conda activate qc
    ```

2. Open FastQC GUI. Analyze and save the reports
    ```
    $ fastqc
    ```
3. See the Basic statistics, Per base quality, Sequence length distribution, Overrepresented sequences and Adapter content sections

### BBDuk

1. Run bbduk. Copied adapters file
    ```
    $ mkdir bb_out
    $ cd bb_out
    $ bbduk.sh in1=../reads/a45_R1.fastq in2=../reads/a45_R2.fastq out1=a45_R1.fastq out2=a45_R2.fastq ref=adapters.fa k=23 mink=7 ktrim=r hdist=1 qtrim=r trimq=20 minlen=100 tpe tbo
    ```

2. Explanation:
3. Result
    >Input:                  	1741880 reads 		436749684 bases.
    QTrimmed:               	1522186 reads (87.39%) 	103642774 bases (23.73%)
    KTrimmed:               	376743 reads (21.63%) 	13100754 bases (3.00%)
    Trimmed by overlap:     	8692 reads (0.50%) 	88862 bases (0.02%)
    Total Removed:          	170588 reads (9.79%) 	116832390 bases (26.75%)
    Result:                 	1571292 reads (90.21%) 	319917294 bases (73.25%)

4. Open FastQC GUI. Analyze and save the reports
    ```
    $ fastqc
    ```
### Trimmomatic

1. Run trimmomatic. Using BBDuk adapters file
    ```
    $ mkdir trim_out
    $ cd trim_out
    $ trimmomatic PE -phred33 ../reads/a45_R1.fastq ../reads/a45_R2.fastq a45_R1_paired.fq.gz a45_R1_unpaired.fq.gz a45_R2_paired.fq.gz a45_R2_unpaired.fq.gz ILLUMINACLIP:../adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:100
    ```
2. Input Read Pairs: 870940 Both Surviving: 599798 (68.87%) Forward Only Surviving: 160897 (18.47%) Reverse Only Surviving: 28249 (3.24%) Dropped: 81996 (9.41%)
3. Open FastQC GUI. Analyze and save the reports
    ```
    $ fastqc
    ```

######################### mapping assembly ############################

# Mapping and Assembly
## (Bacterial Genome Analysis Pipeline)

After filtering the raw reads, you can choose either of the following methods depending on the availability of reference genome or intra-species variations.
Basically, if you have a reference genome and do not expect much variation from it, then the reads are mapped to the reference. Else, de novo assembly is preferred. 

## Mapping to a Reference
1. Index the reference sequence
    ```
    $ mkdir mapping
    $ cd mapping
    $ conda deactivate
    $ conda activate mappers
    $ cp ../resources/NZ_CP053256.1_A45_Chr.fasta ./Ref_A45_chr.fasta
    $ bwa index -a is Ref_A45_chr.fasta
    ```
2. Align Reads separately
    ```
    $ cp ../bb_out/*.fastq ./
    $ bwa aln -t 12 Ref_A45_chr.fasta a45_R1.fastq > a45_R1.sai
    $ bwa aln -t 12 Ref_A45_chr.fasta a45_R2.fastq > a45_R2.sai
    ```
3. Create SAM file and convert it to BAM
    ```
    $  bwa sampe Ref_A45_chr.fasta a45_R1.sai a45_R2.sai a45_R1.fastq a45_R2.fastq > a45_aln.sam
    $ samtools view -S a45_aln.sam -b -o a45_aln.bam
    ```
4.  Sort and index BAM file
    ```
    $ samtools sort a45_aln.bam -o a45_sorted.bam
    $ samtools index a45_sorted.bam
    ```
5. Creating a consensus
    ```
    $ samtools mpileup -uf Ref_A45_chr.fasta a45_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > a45_consensus.fq
    ```
6. Convert to fasta & change the identifier
    ```
    $ conda deactivate
    $ conda activate emboss
    $ seqret -osformat fasta a45_consensus.fq -out2 a45_consensus.fa
    ```
7. Exporting unmapped reads (Optional)
    ```
    $ conda deactivate
    $ conda activate bam2fastq
    $ bam2fastq --no-aligned --force --strict -o a45_unmapped#.fq a45_sorted.bam
    ```
8. Check the rRNA Genes (Optional)
    ```
    $ conda deactivate
    $ conda activate barrnap
    $ barrnap -o cons_rrna.fa < ../mappers/a45_consensus.fa > cons_rrna.gff
    ```            

# de novo Assembly

 ### SPAdes
 1. Spades
    ```
    $ mkdir spades
    $ spades.py --careful --pe1-1 bb_out/a45_R1.fastq --pe1-2 bb_out/a45_R2.fastq -o spades/ --cov-cutoff auto -t 12
    ```
    > Started at 11:37 pm. Ended at 11:43 pm.
    With 4 threads: 11:49 pm to 11.59 pm
    With 6 threads: 9m:54s


2. Results
Check for insert size in the log file & number of contigs/scaffolds in the fasta file.
    ```
    $ grep -i 'insert' spades.log
    $ grep -i '>' scaffolds.fasta -c
    ```

3. Barrnap (Optional)
    ```
    $ mkdir barrnap_out
    $ cd barrnap_out
    $ barrnap -o spades_rrna.fa < ../spades/scaffolds.fasta > spades_rrna.gff
    ```
======================
## Alternatives
### IDBA
```
$ fq2fa --merge --filter ../bb_out/a45_R1.fastq ../bb_out/a45_R2.fastq a45_reads.fa
$ idba_ud -r a45_reads.fa -o ./
```
### Velvet
```
$ VelvetOptimiser.pl -s 79 -e 159 -f '-shortPaired -fastq -separate ../bb_out/a45_R1.fastq ../bb_out/a45_R2.fastq' -t 12 -d ./ -v
```
### Unicycler
```
$ unicycler -1 ../bb_out/a45_R1.fastq -2 ../bb_out/a45_R2.fastq -o ./ -t 12
```
======================

########################### Deno assembly #############################

# Steps after De-novo Assembly

### Contig Management

>Mamba installation START ---
Yesterday while installing BUSCO, there was an error and could not proceed. It seems the anaconda package manager is compatible to deal with envirnments with lots of packages.
https://stackoverflow.com/questions/72743734/condas-solving-environment-takes-forever
The alternative is Mambaforge. Can be downloaded from here...
https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh

>Before installing mamba, you have to uninstall the anaconda installation.
    1. Go to Home directory and delete the complete anaconda3 directory.
    2. Remove .condarc file (Press Ctrl+H to see hidden files) from Home directory
    3. In Home directory you will find .bashrc also. Open it with a text editor and delete the lines starting with
    # >>> conda initialize >>>
    and ending with
    # <<< conda initialize <<<
    These lines are usually at the end of the file.
    4. After deleting the lines, save it.
    5. Download the mambaforge file from the above link and install in the same way as Anaconda
    6. Go to the downloaded folder.
    7. Right Click and select "Open in Terminal" option
    $ chmod +x Mambaforge-Linux-x86_64.sh
    $ ./Mambaforge-Linux-x86_64.sh
    8. Follow on-screen instructions
    9. While using mamba, it is similar to conda, only thing is we have to replace conda with mamba
    Eg. Instead of conda deactivate, you should use mamba deactivate. Similarly, mamba create -n envname
    10. Note that adding channels is done with conda (not with mamba). So, the following steps to add bioconda channel will be same...
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict
    11. Apart from these all commands will use mamba
    12. Since the anaconda3 directory was deleted, all packages should be installed again. 
    Don't worry. As the name suggests, Mamba is faster than Conda.
    --- Mamba Installation END ---


### MeDuSa
1. Note If there are any colons in the header of fasta
    ```
    $ sed -i 's/:/_/g' scaffolds.fasta
    here, "scaffolds.fasta" spades result
    (If there are any colons in the header of fasta)
    ```
2. Make a new directory - medusa_out
    ```
    $ mkdir medusa_out
    ```
3. Inside medusa_out directory create a new folder Ref then download chromosome and plasmid reference data and Merge Chromosome and Plasmid and save it as full genome. If reference have 2 choromosome only then you can merge and save it full genome or if 3-4 plasmid then also do same merege all and make it full genome.
    ```
    $ cd medusa_out
    $ mkdir Ref
    $ cat Ref_A45_chr.fasta Ref_A45_p.fasta > medusa_out/Ref/Ref_A45_full.fasta
    Ref_A45_chr.fasta Ref_A45_plasmid.fasta reference file yoy can download it from NCBI (https://www.ncbi.nlm.nih.gov/genome/169?genome_assembly_id=901025)
    ```
4. Copy the scaffolds file to current directory (scaffolds.fasta from the spades directory)
Create new environment and install medusa
    ```
    $ mamba deactivate
    $ mamba create -n medusa -y
    $ mamba activate medusa
    $ mamba install -c conda-forge -c bioconda medusa
    $ mamba install -c conda-forge mummer
    $ mamba install -c conda-forge biopython
    ```
5. Run the medusa command
    ```
    $ cd  medusa_out/Ref
    $ medusa -d -f ../Ref/ -i ../Ref/scaffolds.fasta -random 10 -w2 -v
    ```
> If you face cPickle error
open python file and Change cPickle to pickle in Home/mambaforge/envs/medusa/share/medusa-1.6-2/script/netcon_mummer.py

### Mauve

>Mauve 
>Before proceeding with mauve you have to check how many no. of scaffolds are there in  medusa_out/"scaffolds.fastaScaffold.fasta" file , if there is only one or two you can skip Mauve and proceed with gapcloser to close the no of N from genome(scaffolds.fastaScaffold.fasta file).

1. Create mauve env.
    ```
    $ mamba deactivate
    $ mamba create -n mauve -y
    $ mamba activate mauve
    ```
2. Install Mauve
    ```
    $ mamba install -c bioconda mauve
    $ Mauve
    ```
    > This has Graphical UI
    Do Progressive alignemnt
    Then Reorder the contigs
    > tools > move contigs > choose output folder > add sequence(medusa_out/Ref_A45_full.fasta this is combined file plasmid and chromosome) and again add sequence(medusa_out/scaffolds.fastaScaffold.fasta) > start

### GapCloser
1. Create filler env.
    ```
    $ mamba deactivate 
    $ mamba create -n fillers -y
    $ mamba activate fillers
    ```
2. Install soapdenovo2-gapcloser
    
    >gapcloser required config file as a input with some detais of your sequence like raw reads length
    max_read_lenth select from bb_out result html file and avg_ins (insert size) should be select from spades.log you gan grep "Insert" from spades.log file
    >max_rd_len=250
    [LIB]
    name=a45
    avg_ins=452
    reverse_seq=0
    asm_flags=4
    rank=1
    pair_num_cutoff=3
    map_len=32
    q1=bb_out/a45_R1.fastq
    q2=bb_out/a45_R2.fastq
    ```
    $ mamba install -c bioconda soapdenovo2-gapcloser
    >You have to make config file as save is as "BSE6-1_GC.config"
    $ GapCloser -a medusa_out/Ref/scaffolds.fastaScaffold.fasta -b BSE6-1_GC.config -o BSE6-1_GC0.fasta -t 12
    ```

### Pilon
1. Create pilon env. and install
    ```
    $ mamba deactivate
    $ mamba create -n pilon -y
    $ mamba activate pilon
    $ mamba install -c bioconda pilon
    ```
2. Before running Pilon we have to index the fasta file and reads
    ```
    $ mkdir pilon_out
    $ cd pilon_out
    $ mamba deactivate
    $ mamba activate mappers
    $ mamba install -c bioconda bowtie2
    $ cp ../filler/a45_GC.fasta ./
    ```
# Steps for Reference based genome analysis: Index the genome
1. Indexing the genome with bowtie2
    ```
    $ bowtie2-build a45_GC.fasta a45
    Align reads to genome (5 mins with 12 cores)
    $ bowtie2 -x a45 -1 ../bb_out/a45_R1.fastq -2 ../bb_out/a45_R2.fastq -S reads_on_assembly.sam -p 12
    ```
2. Convert SAM to BAM, sort and index
    ```
    $ samtools view reads_on_assembly.sam -b -o reads_on_assembly.bam
    $ samtools sort reads_on_assembly.bam -o reads_on_assembly_sorted.bam
    $ samtools index reads_on_assembly_sorted.bam
    ```
3. Activate pilon
    ```
    $ mamba deactivate
    $ mamba activate pilon
    $ pilon --genome a45_GC.fasta --frags reads_on_assembly_sorted.bam
    ```
If there is a memory error open the following file
> /home/anwesh/mambaforge/envs/pilon/share/pilon-1.24-0
And increase the max memory option from
default_jvm_mem_opts = ['-Xms512m', '-Xmx1g'] to whatever your RAM has (I increase it from 1g to 4g)
default_jvm_mem_opts = ['-Xms512m', '-Xmx4g']

### BUSCO
1. Create busco env and install
    ```
    $ mamba create -n busco -y
    $ mamba activate busco
    $ mamba install -c conda-forge -c bioconda busco=5.4*
    ```
2. Run BUSCO
    ```
    $ busco -m genome -i ../pilon_out/pilon.fasta -o a45 --auto-lineage-prok -c 10
    ```
### CheckM
1. Create CheckM env and install
    ```
    $ mamba deactivate
    $ mamba create -n checkm python=3.9
    $ mamba activate checkm
    $ mamba install numpy matplotlib pysam
    $ mamba install hmmer prodigal pplacer
    $ pip install checkm-genome
    ```
2. Check the installation
    ```
    $ checkm
    ```
    ***Note: Copy spades.fasta pilon.fasta and a45_GC.fasta into a new directory - genomes
Create new directory checkm_out and navigate into it***

3. Run checkm (3 mins)
    ```
    $ checkm lineage_wf -x fasta ../genomes/ ./ -r -f ./results.txt -t 12
    ```
4. Install assembly-stats
    ```
    $ mamba install -c bioconda assembly-stats
    $ assembly-stats ../genomes/*.fasta > assembly_stats.txt
    ```
Split the assembly into Chr and plasmid

RAST
Online tool


#######################genome annotation###############################

# Genome Annotation 
### Prokka
1. Create env and install Prokka
    ```
    $ conda install -c conda-forge -c bioconda -c defaults prokka
    $ mamba create -n prokka -y
    $ mamba activate prokka
    $ mamba install -c conda-forge -c bioconda prokka
    ```
    >Download GBK files as we have a reference genome.
    https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP053256
2. make a directory and run prokka
    ```
    $ mkdir prokka_out
    $ prokka --outdir a45 --prefix a45 genomes/pilon.fasta
    ```
    > pilon.fasta is a genome sequence
    Download reference genomoes and run prokka on all of them

### Roary
1. Create Roary env. and install
   >Run roary if you have more then one genome, roary takes multiple gff file as input of prokka.
   if you have only one genome then skip PAN genome (roray) and procced to next step
    
    ```
    $ mamba create -n roary -y
    $ mamba activate roary
    $ mamba install -c bioconda roary
    ```

3. Make a new directory raory and navigate into it
Run roary ( mins)
    ```
    $ roary ../prokka_out/*/*.gff -e -n -r -v -f tenRefs
    Error: not found File::Find::Rule
    $ cpan File::Find::Rule
    ```
    >FriPan (https://github.com/drpowell/FriPan)
    Change the server module in server.sh, as suggested in this site
    https://stackoverflow.com/questions/17351016/set-up-python-simplehttpserver-on-windows
    Copy roray output to FriPan root and rename it to filename.roary

### Mafft 
1. Create Mafft env. and install
    ```
    $ mamba deactivate
    $ mamba create -n phylogeny -y
    $ mamba activate phylogeny
    $ mamba install -c bioconda mafft
    Run MAFFT
    $ mafft --maxiterate 100 --reorder --thread 10 16S_a45-Ref-Out.fasta > 16S_a45-Ref-Out_aln.fasta
    ```
    >RaxML GUI - https://github.com/AntonelliLab/raxmlGUI/releases/latest/download/raxmlGUI-2.0.10.AppImage
    raxmlHPC-PTHREADS-SSE3 -T 10 -f a -x 288426 -p 288426 -N 100 -m GTRGAMMA -O -o H_acinonychis -n 16S -s 16S_BSE-Ref-Out_aln_modified.fasta 

### TYGS

### Dotplot with D-Genies

### antiSMASH

### Circos
1. Create Circos env. and install
    ```
    $ mamba deactivate
    $ mamba create -n circos -y
    $ mamba activate circos
    $ mamba install -c bioconda circos
    ```
### KBase



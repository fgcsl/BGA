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

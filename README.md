# BGA
Bacterial Genome Analysis Piplene

# Bacterial Genome analysis Pipeline -  DBT-CMI 

### Requirements Tools & pacakges
- Mamba
    You can install mamba by running "install_mamba.sh" script
    ```
    $ bash install_mamba.sh
    ```
- Stand alone tools
  ```
  $ bash install_tools.sh
  $ bash install_checkm.sh
  ```
  - 
### Requirements (QC and Assembly)

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

Bacterial genome analysis - Reference based
===================================================

* * *

Raw Reads QC
------------

FastQC
------

1.  Navigate to bgap directory and activate qc
    ```              
    $ cd Desktop/bgap
    $ mamba deactivate
    $ mamba activate qc
    ```     
    
3.  Open FastQC GUI. Analyze and save the reports
    ```
    (qc)$ fastqc
    ```   
    
4.  See the Basic statistics, Per base quality, Sequence length distribution, Overrepresented sequences and Adapter content sections

### BBDuk

1.  Run bbduk. Copied adapters file
    ```
    (qc)$ mkdir bb_out
    (qc)$ cd bb_out
    (qc)$ bbduk.sh in1=../reads/a45_R1.fastq in2=../reads/a45_R2.fastq out1=a45_R1.fastq out2=a45_R2.fastq ref=adapters.fa k=23 mink=7 ktrim=r hdist=1 qtrim=r trimq=20 minlen=100 tpe tbo
    ```    
    
2.  Result
    
    > Input: 1741880 reads 436749684 bases.  
    > QTrimmed: 1522186 reads (87.39%) 103642774 bases (23.73%)  
    > KTrimmed: 376743 reads (21.63%) 13100754 bases (3.00%)  
    > Trimmed by overlap: 8692 reads (0.50%) 88862 bases (0.02%)  
    > Total Removed: 170588 reads (9.79%) 116832390 bases (26.75%)  
    > Result: 1571292 reads (90.21%) 319917294 bases (73.25%)
    
3.  Open FastQC GUI. Analyze and save the reports
    ``` 
    (qc)$ fastqc
    ```     
    

### Trimmomatic

1.  Run trimmomatic. Using BBDuk adapters file
    ```   
    (qc)$ mkdir trim_out
    (qc)$ cd trim_out
    (qc)$ trimmomatic PE -phred33 ../reads/a45_R1.fastq ../reads/a45_R2.fastq a45_R1_paired.fq.gz a45_R1_unpaired.fq.gz a45_R2_paired.fq.gz a45_R2_unpaired.fq.gz ILLUMINACLIP:../adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:100
    ```      
    
3.  Input Read Pairs: 870940 Both Surviving: 599798 (68.87%) Forward Only Surviving: 160897 (18.47%) Reverse Only Surviving: 28249 (3.24%) Dropped: 81996 (9.41%)
4.  Open FastQC GUI. Analyze and save the reports
     ```   
     (qc)$ fastqc
     ```    
    

Mapping
=======

After filtering the raw reads, you can choose either of the following methods depending on the availability of reference genome or intra-species variations.  
Basically, if you have a reference genome and do not expect much variation from it, then the reads are mapped to the reference. Else, de novo assembly is preferred.

### Mapping to a Reference

1.  Index the reference sequence
    ```
    (qc)$ mkdir mapping
    (qc)$ cd mapping
    (qc)$ mamba deactivate
    (qc)$ mamba activate mappers
    (qc)$ cp ../resources/NZ_CP053256.1_A45_Chr.fasta ./Ref_A45_chr.fasta
    (qc)$ bwa index -a is Ref_A45_chr.fasta
    ```  
    
2.  Align Reads separately
    ```
    (qc)$ cp ../bb_out/*.fastq ./
    (qc)$ bwa aln -t 12 Ref_A45_chr.fasta a45_R1.fastq > a45_R1.sai
    (qc)$ bwa aln -t 12 Ref_A45_chr.fasta a45_R2.fastq > a45_R2.sai
    ``` 
    
3.  Create SAM file and convert it to BAM
    ```
    (qc)$ bwa sampe Ref_A45_chr.fasta a45_R1.sai a45_R2.sai a45_R1.fastq a45_R2.fastq > a45_aln.sam
    (qc)$ samtools view -S a45_aln.sam -b -o a45_aln.bam
    ```
    
4.  Sort and index BAM file
    ```
    (qc)$ samtools sort a45_aln.bam -o a45_sorted.bam
    (qc)$ samtools index a45_sorted.bam
     ``` 
    
5. Creating a consensus
    ```
    (qc)$ samtools mpileup -uf Ref_A45_chr.fasta a45_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > a45_consensus.fq
    ``` 
    
6. Convert to fasta & change the identifier
    ```
    (qc)$ mamba deactivate
    (base)$ mamba activate emboss
    (emboss)$ seqret -osformat fasta a45_consensus.fq -out2 a45_consensus.fa
    ```  
    
7. Exporting unmapped reads (Optional)
    ``` 
    (emboss)$ mamba deactivate
    (base)$ mamba activate bam2fastq
    (bam2fastq)$ bam2fastq --no-aligned --force --strict -o a45_unmapped#.fq a45_sorted.bam
    ```

8.  Check the rRNA Genes (Optional)
    ```
    (bam2fastq)$ mamba deactivate
    (base)$ mamba activate barrnap
    (barrnap)$ barrnap -o cons_rrna.fa < ../mappers/a45_consensus.fa > cons_rrna.gff
    ```     
    
9.  Close the gap
    ```
    (bam2fastq)$ mamba deactivate
    (base)$ mamba activate filler
    (filler)$ barrnap -o cons_rrna.fa < ../mapping/a45_consensus.fa > cons_rrna.gff
    ```

# de novo Assembly

Bacterial Genome analysis Pipeline.md

de novo Assembly (DAY-2)
========================

### SPAdes

1.  Spades  
    ```
    (base)$ mkdir denovo
    (base)$ cd denovo
    (base)$ mkdir spades
    (base)$ mamba activate assembler
    (assembler)$ spades.py --careful --pe1-1 bb_out/a45_R1.fastq --pe1-2 bb_out/a45_R2.fastq -o spades/ --cov-cutoff auto -t 12
    ```    
    > Started at 11:37 pm. Ended at 11:43 pm.  
    > With 4 threads: 11:49 pm to 11.59 pm  
    > With 6 threads: 9m:54s
    
3.  Results
   
    Check for insert size in the log file & number of contigs/scaffolds in the fasta file.
    ```
    (assembler)$ grep -i 'insert' spades.log
    (assembler)$ grep -i '>' scaffolds.fasta -c
    ```  
    
5.  Barrnap (Optional)
    ```
    $ mkdir barrnap_out
    $ cd barrnap_out
    $ barrnap -o spades_rrna.fa < ../spades/scaffolds.fasta > spades_rrna.gff
    ```    
   
===================

Alternatives
------------

1.  ### IDBA
    ```
    $ fq2fa --merge --filter ../bb_out/a45_R1.fastq ../bb_out/a45_R2.fastq a45_reads.fa
    $ idba_ud -r a45_reads.fa -o ./
    ```   
    
2.  ### Velvet
    ```
    $ VelvetOptimiser.pl -s 79 -e 159 -f '-shortPaired -fastq -separate ../bb_out/a45_R1.fastq ../bb_out/a45_R2.fastq' -t 12 -d ./ -v
    ```    
    
3.  ### Unicycler
    ```
    $ unicycler -1 ../bb_out/a45_R1.fastq -2 ../bb_out/a45_R2.fastq -o ./ -t 12
    ```
    
======================

########################### Deno assembly #############################

# Steps after De-novo Assembly

### Contig Management
### Contig Management

### MeDuSa

1.  Note If there are any colons in the header of fasta
    ```
    (base)$ sed -i 's/:/_/g' scaffolds.fasta
    ```  
    
2.  Make a new directory - medusa\_out and Copy the scaffolds file in medusa\_out directory
    ```
    (base)$ mkdir medusa_out
    (base)$ cp spades/scaffolds.fasta medusa_out
    ```  
    
3.  Inside medusa\_out directory create a new folder Ref and download chromosome & plasmid reference data from ([https://www.ncbi.nlm.nih.gov/genome/169?genome\_assembly\_id=901025](https://www.ncbi.nlm.nih.gov/genome/169?genome_assembly_id=901025)) and Merge both files, save it as full genome.  
    **Note: keep only merged files inside Ref directory**
    ```
    (base)$ cd medusa_out
    (base)$ mkdir Ref
    (base)$ cat Ref_A45_chr.fasta Ref_A45_p.fasta > Ref/Ref_A45_full.fasta
    ``` 
    
4. Copy the scaffolds file to current directory (scaffolds.fasta from the spades directory)  
    Create new environment and install medusa
    ```
    (base)$ mamba deactivate
    $ mamba activate scaffolder
    ```   
    
5. Run the medusa command
    ```
    (scaffolder)$ medusa -d -f Ref/ -i scaffolds.fasta -random 10 -w2 -v
    ```  
> If you face cPickle error  
> open python file and Change cPickle to pickle in Home/mambaforge/envs/medusa/share/medusa-1.6-2/script/netcon\_mummer.py

### Mauve

> Mauve  
> Before proceeding with mauve you have to check how many no. of “N” “n” (gaps) are there in “medusa\_out/scaffolds.fastaScaffold.fasta” file.

2.  Activate mauve env.
   
    ```
    (scaffolder)$ cd ../
    (scaffolder)$ mamba deactivate
    (scaffolder)$ mamba activate mauve
    ```
    
4.  Run Mauve
    
    ```
    (mauve)$ mkdir mauve_out
    (mauve)$ Mauve
    (mauve)$ cd mauve_out
    ```
    
    > This has Graphical UI  
    > Do Progressive alignemnt  
    > Then Reorder the contigs  
    > tools > move contigs > choose output folder (choose mauve\_out folder where output will save) > add sequence(medusa\_out/Ref\_A45\_full.fasta this is combined file plasmid and chromosome) and again add sequence(medusa\_out/scaffolds.fastaScaffold.fasta) > start
    

### GapCloser

1.  Activate filler env.
    ```
    (mauve)$ cd ../
    (mauve)$ mamba deactivate 
    (base)$ mamba activate fillers
    ```  
    
2.  Run filler (soapdenovo2-gapcloser)
    
    > to run gapcloser you need to create config file with some sequence parameters like raw reads length  
    > max\_read\_lenth (select from bb\_out result html file) and avg\_ins (insert size: select from spades.log you can grep “Insert” from spades.log file)
    
    a45\_GC.config file
    ```
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
    ```
    (fillers)$ mkdir filler_out
    (fillers)$ GapCloser -a mauve_out/alignment2/scaffolds.fastaScaffold.fasta -b a45_GC.config -o filler_out/a45_GC.fasta -t 6
    ``` 
    
    OR (If you are skipping mauve then you can use medusa\_out as input)
    

### Use Pilon if there is still gaps (N) in Gapcloser result

1.  Before running Pilon we have to index the fasta file and reads
    ```
    (fillers)$ mkdir pilon_out
    (fillers)$ cd pilon_out
    (fillers)$ mamba deactivate
    (base)$ mamba activate mappers
    ```  
    
2.  Index the genome
    ```
    (mappers)$ bowtie2-build ../filler/a45_GC.fasta a45
    ```   
    
3.  Align reads to genome (5 mins with 12 cores)
    ```
    (mappers)$ bowtie2 -x a45 -1 ../bb_out/a45_R1.fastq -2 ../bb_out/a45_R2.fastq -S reads_on_assembly.sam -p 12
    ###### Convert SAM to BAM, sort and index
    (mappers)$ samtools view reads_on_assembly.sam -b -o reads_on_assembly.bam
    (mappers)$ samtools sort reads_on_assembly.bam -o reads_on_assembly_sorted.bam
    (mappers)$ samtools index reads_on_assembly_sorted.bam
    ``` 
    
4.  Activate pilon
    ```
    (mappers)$ mamba deactivate
    (mappers)$ mamba activate pilon
    (pilon)$ pilon --genome ../filler_out/a45_GC.fasta --frags reads_on_assembly_sorted.bam
    ```   
    
    If there is a memory error open the following file
    
    > /home/anwesh/mambaforge/envs/pilon/share/pilon-1.24-0  
    > And increase the max memory option from
        default_jvm_mem_opts = ['-Xms512m', '-Xmx1g'] to whatever your RAM has (I increase it from 1g to 4g)
        default_jvm_mem_opts = ['-Xms512m', '-Xmx4g']
        
    

### BUSCO (for qc)

1.  Activate busco env.
    ```  
    (pilon)$ cd ../
    (pilon)$mamba deactivate
    (base)$ mamba activate busco
    (busco)$ mamba install -c conda-forge -c bioconda busco=5.4*
    (busco)$ mkdir busco_qc
    ```  
    
2.  Run BUSCO
    ```
    (busco)$ busco -m genome -i filler_out/a45_GC.fasta -o busco_qc/a45 --auto-lineage-prok -c 10
    OR
    (busco)$ busco -m genome -i pilon_out/pilon.fasta -o a45 --auto-lineage-prok -c 10
    ```   
    

### CheckM

1.  Activate Checkm
    ```
    (busco)$ mamba deactivate
    (busco)$ mamba activate checkm
    (checkm)$ mkdir checkm_out
    (checkm)$ cd checkm_out
    ```    
    

**_Note: Copy spades.fasta pilon.fasta or gapcloser.fasta and a45\_GC.fasta into a new directory - genomes  
Create new directory checkm\_out and navigate into it_**
    ```
    (checkm)$ mkdir genomes
    (checkm)$ cp ../filler_out/a45_GC.fasta genomes/
    (checkm)$ cp ../spades/scaffolds.fasta genomes/
    (checkm)$ cp ../medusa_out/scaffolds.fastaScaffold.fasta genomes/
     ```   

4. Run checkm
   ```
   (checkm)$ checkm lineage_wf -x fasta ../genomes/ ./ -r -f ./results.txt -t 12
   ```  
    
    > If you face FileNotFoundError: \[Errno 2\] No such file or directory: ‘/home/dbt-cmi/.checkm/hmms/phylo.hmm’ then download reference data ([https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm](https://github.com/Ecogenomics/CheckM/wiki/Installation#how-to-install-checkm)) and extract it in path “/home/dbt-cmi/.checkm/”  
    > Read Installation and download Required reference data ([https://data.ace.uq.edu.au/public/CheckM\_databases](https://data.ace.uq.edu.au/public/CheckM_databases))
    
5.  Run and Check asssembly stats
    ```
    (checkm)$ assembly-stats ../genomes/*.fasta > assembly_stats.txt
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



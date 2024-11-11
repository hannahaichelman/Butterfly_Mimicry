# This workflow includes information on creating the assembly and annotation of Limenitis lorquini genome
# Author: Hannah Aichelman (hannahaichelman @ gmail.com)

# back-up of raw genome files can be found in the BU SCC directory below - i deleted them from my working directory because of space limitations
/projectnb/mullenl/evan/sequence_data/genomes/limenitis_lorquini_genome_raw_sequence/Raw_Data

# Working here in my directory for these analyses:
pwd
/projectnb/mullenl/hannah/limenitis_genomes

# ---------- Set Up Conda Environment ---------- #

module load miniconda/23.1.0
# prompted to run this:
setup_scc_condarc.sh
#A .condarc file has been written to your home directory:  ~/.condarc
#Conda environments and packages will be stored in the directory:
#     /projectnb/mullenl/haich/.conda/

conda --version
# conda 23.1.0

conda config --add channels defaults
# Warning: 'defaults' already in 'channels' list, moving to the top
conda config --add channels bioconda
conda config --add channels conda-forge

# create new conda environment for Genome Assembly
conda create --name Genome

# got a prompt to update conda
conda update -n base -c defaults conda

# activate conda environment
conda activate Genome

# install packages necessary for genome assembly
# for example, install the program "hifiasm"
conda install hifiasm
# worked!

# Install BBTools on conda Genome environment
conda install -c agbiome bbtools

# deactivate conda environment
conda deactivate

# ---------- Filter out Remnant PacBio Adapters with HifiAdaptFilt ---------- #

# Notes on Parameters used to filter adapters using HiFiAdapterFilt
# Used default Parameters
# Removing reads containing adapters a minimum of 44 bp in length and 97% match.
	# From GitHub README
		# [ -l minimum Length of adapter match to remove. Default=44 ] 
		# [ -m minimum percent Match of adapter to remove. Default=97 ]

pwd
/projectnb/mullenl/hannah/limenitis_genomes/software

git clone https://github.com/sheinasim/HiFiAdapterFilt.git # cloned the whole repository in the pwd

# clones the "HiFiAdapterFilt" directory with all its scripts/files as well as its subdirectory "DB"
# didn't need to do chmod at all, seems to have worked just fine, later I add this "HiFiAdapterFilt" directory and its subdirectory "DB" to the path (temporarily?)

# make sure HiFiAdapterFilt dependencies are installed in conda environment

# activate Genome environment
conda activate Genome

# install dependency packages
conda install -c bioconda bamtools
conda install blast
conda deactivate

# running HifiAdapterFilt with the original bam files

pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/

# made job file before entering conda environment
(Genome)[haich@scc1 raw_data]$ cat HifiAdaptFiltpac 
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N hifiadaptfilt # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -pe omp 10

/projectnb/mullenl/hannah/limenitis_genomes/software/HiFiAdapterFilt/pbadapterfilt.sh m64190e_210221_233100.hifi_reads.bam

# ORIGINAL HIFI FILES:                         
m64190e_210221_233100.hifi_reads.fasta.gz
m64190e_210221_233100.hifi_reads.fastq.gz
m64190e_210221_233100.hifi_reads.bam      
m64190e_210221_233100.hifi_reads.bam.pbi

# OUTPUT FILES:
hifiadaptfilt.o799014
hifiadaptfilt.e799014
m64190e_210221_233100.hifi_reads.blocklist
m64190e_210221_233100.hifi_reads.contaminant.blastout
m64190e_210221_233100.hifi_reads.filt.fastq.gz
m64190e_210221_233100.hifi_reads.pbi.stats
m64190e_210221_233100.hifi_reads.stats

# make directory for filtered output - moved all "OUTPUT FILES" above into the directory below
mkdir hifi_filt
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt

# ---------- Pre-assembly stats with Jellyfish and Genomescope ---------- #
# install necessary programs in conda environment
module load miniconda
conda activate Genome
conda install kmer-jellyfish
conda deactivate

# make job file to gunzip filtered file and run jellyfish to count k-mers then produce histogram
# jellyfish details: https://github.com/gmarcais/Jellyfish/tree/master/doc
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt

nano gunzip_jellyfish
#!/bin/bash
#$ -V # inherit the submission environment
#$ -N gunzip_jellyfish # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -pe omp 10

gunzip m64190e_210221_233100.hifi_reads.filt.fastq.gz
jellyfish count -C -m 21 -s 1000000000 -t 10 m64190e_210221_233100.hifi_reads.filt.fastq -o hifi.filt.reads.jf
jellyfish histo -t 10 hifi.filt.reads.jf > hifi.filt.reads.histo

# submit job
conda activate Genome
qsub gunzip_jellyfish

# upload .histo file into the genomescope GUI: http://qb.cshl.edu/genomescope/
# results: http://genomescope.org/analysis.php?code=QY8HjcEyHuISBWtw6uxi
# interpreting the results: https://bioinformaticsworkbook.org/dataAnalysis/GenomeAssembly/genomescope.html#gsc.tab=0


# ---------- Genome assembly with hifiasm ---------- #

# https://hifiasm.readthedocs.io/en/latest/pa-assembly.html
# https://github.com/chhylp123/hifiasm/issues/61
# https://hifiasm.readthedocs.io/en/latest/interpreting-output.html 

# make file for job submission
nano hifiasm
#!/bin/bash
#$ -V # inherit the submission environment
#$ -N hifiasm # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -pe omp 16

hifiasm -o lim_hifiasm_genome.asm -t 16 m64190e_210221_233100.hifi_reads.filt.fastq
# -o = prefix of output files
# -t = # of CPUs
# 1.fq.gz = input files

# activate conda environment and submit
conda activate Genome
qsub hifiasm

# OUTPUT files:
lim_hifiasm_genome.asm.bp.hap1.p_ctg.gfa
lim_hifiasm_genome.asm.bp.hap1.p_ctg.lowQ.bed
lim_hifiasm_genome.asm.bp.hap1.p_ctg.noseq.gfa
lim_hifiasm_genome.asm.bp.hap2.p_ctg.gfa
lim_hifiasm_genome.asm.bp.hap2.p_ctg.lowQ.bed
lim_hifiasm_genome.asm.bp.hap2.p_ctg.noseq.gfa
lim_hifiasm_genome.asm.bp.p_ctg.gfa
lim_hifiasm_genome.asm.bp.p_ctg.lowQ.bed
lim_hifiasm_genome.asm.bp.p_ctg.noseq.gfa
lim_hifiasm_genome.asm.bp.p_utg.gfa
lim_hifiasm_genome.asm.bp.p_utg.lowQ.bed
lim_hifiasm_genome.asm.bp.p_utg.noseq.gfa
lim_hifiasm_genome.asm.bp.r_utg.gfa
lim_hifiasm_genome.asm.bp.r_utg.lowQ.bed
lim_hifiasm_genome.asm.bp.r_utg.noseq.gfa
lim_hifiasm_genome.asm.ec.bin
lim_hifiasm_genome.asm.ovlp.reverse.bin
lim_hifiasm_genome.asm.ovlp.source.bin

# 3 assembly files:
lim_hifiasm_genome.asm.bp.hap1.p_ctg.gfa
lim_hifiasm_genome.asm.bp.hap2.p_ctg.gfa
lim_hifiasm_genome.asm.bp.p_ctg.gfa

# convert .gfa files to .fasta files using "gfatools"
conda install -c ohmeta gfatools

gfatools gfa2fa lim_hifiasm_genome.asm.bp.hap1.p_ctg.gfa > lim_hifiasm_genome.asm.bp.hap1.p_ctg.fa
gfatools gfa2fa lim_hifiasm_genome.asm.bp.hap2.p_ctg.gfa > lim_hifiasm_genome.asm.bp.hap2.p_ctg.fa
gfatools gfa2fa lim_hifiasm_genome.asm.bp.p_ctg.gfa > lim_hifiasm_genome.asm.bp.p_ctg.fa

# new output files:
lim_hifiasm_genome.asm.bp.hap1.p_ctg.fa  
lim_hifiasm_genome.asm.bp.hap2.p_ctg.fa  
lim_hifiasm_genome.asm.bp.p_ctg.fa

# ---------- Evaluating genome quality ---------- #

conda install -c bioconda htseq 

# Import stats script from this github: https://github.com/thh32/Assembly-stats
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/

git clone https://github.com/thh32/Assembly-stats.git

python Assembly-stats/Assembly_stats.py -i lim_hifiasm_genome.asm.bp.hap1.p_ctg.fa -o lim_hifiasm_genome.asm.bp.hap1.p_ctg_stats
There are 283 contigs/reads in this file
Total number of contigs is: 283
Total number of bases in all contigs/reads is: 373582490
The halfway value is: 186791245
N50 is: 12332538
17550572 is the size of the largest contig/read
10775 is the smallest contig/read

python Assembly-stats/Assembly_stats.py -i lim_hifiasm_genome.asm.bp.hap2.p_ctg.fa -o lim_hifiasm_genome.asm.bp.hap2.p_ctg_stats
There are 133 contigs/reads in this file
Total number of contigs is: 133
Total number of bases in all contigs/reads is: 362572347
The halfway value is: 181286173
N50 is: 11551674
17552827 is the size of the largest contig/read
15850 is the smallest contig/read

python Assembly-stats/Assembly_stats.py -i lim_hifiasm_genome.asm.bp.p_ctg.fa -o lim_hifiasm_genome.asm.bp.p_ctg_stats
There are 251 contigs/reads in this file
Total number of contigs is: 251
Total number of bases in all contigs/reads is: 375736066
The halfway value is: 187868033
N50 is: 12744042
17552827 is the size of the largest contig/read
9574 is the smallest contig/read

# Get more stats / verify stats with stats.sh script that I utilized above
# run this as a job
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt

qsub assembly_stats_hifi.qsub
# output is assembly_stats_hifi.qlog:
==========================================================
Starting on : Thu Sep 14 16:03:49 EDT 2023
Running on node : scc-mb5
Current directory : /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt
Current job ID : 875283
Current job name : assembly_stats_hifi
==========================================================
A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.3329	0.1672	0.1671	0.3328	0.0000	0.0000	0.0000	0.3343	0.1465

Main genome scaffold total:         	283
Main genome contig total:           	283
Main genome scaffold sequence total:	373.582 MB
Main genome contig sequence total:  	373.582 MB  	0.000% gap
Main genome scaffold N/L50:         	14/12.629 MB
Main genome contig N/L50:           	14/12.629 MB
Main genome scaffold N/L90:         	31/4.235 MB
Main genome contig N/L90:           	31/4.235 MB
Max scaffold length:                	17.551 MB
Max contig length:                  	17.551 MB
Number of scaffolds > 50 KB:        	106
% main genome in scaffolds > 50 KB: 	98.70%


Minimum 	Number        	Number        	Total         	Total         	Scaffold
Scaffold	of            	of            	Scaffold      	Contig        	Contig  
Length  	Scaffolds     	Contigs       	Length        	Length        	Coverage
--------	--------------	--------------	--------------	--------------	--------
    All 	           283	           283	   373,582,490	   373,582,490	 100.00%
  10 KB 	           283	           283	   373,582,490	   373,582,490	 100.00%
  25 KB 	           210	           210	   372,120,989	   372,120,989	 100.00%
  50 KB 	           106	           106	   368,720,699	   368,720,699	 100.00%
 100 KB 	            77	            77	   366,805,521	   366,805,521	 100.00%
 250 KB 	            56	            56	   363,402,937	   363,402,937	 100.00%
 500 KB 	            46	            46	   359,823,133	   359,823,133	 100.00%
   1 MB 	            40	            40	   355,794,270	   355,794,270	 100.00%
 2.5 MB 	            33	            33	   341,709,866	   341,709,866	 100.00%
   5 MB 	            29	            29	   327,204,508	   327,204,508	 100.00%
  10 MB 	            21	            21	   270,814,303	   270,814,303	 100.00%

A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.3353	0.1647	0.1647	0.3353	0.0000	0.0000	0.0000	0.3294	0.1044

Main genome scaffold total:         	133
Main genome contig total:           	133
Main genome scaffold sequence total:	362.572 MB
Main genome contig sequence total:  	362.572 MB  	0.000% gap
Main genome scaffold N/L50:         	14/11.639 MB
Main genome contig N/L50:           	14/11.639 MB
Main genome scaffold N/L90:         	28/7.151 MB
Main genome contig N/L90:           	28/7.151 MB
Max scaffold length:                	17.553 MB
Max contig length:                  	17.553 MB
Number of scaffolds > 50 KB:        	66
% main genome in scaffolds > 50 KB: 	99.45%


Minimum 	Number        	Number        	Total         	Total         	Scaffold
Scaffold	of            	of            	Scaffold      	Contig        	Contig  
Length  	Scaffolds     	Contigs       	Length        	Length        	Coverage
--------	--------------	--------------	--------------	--------------	--------
    All 	           133	           133	   362,572,347	   362,572,347	 100.00%
  10 KB 	           133	           133	   362,572,347	   362,572,347	 100.00%
  25 KB 	           111	           111	   362,122,677	   362,122,677	 100.00%
  50 KB 	            66	            66	   360,583,416	   360,583,416	 100.00%
 100 KB 	            48	            48	   359,362,407	   359,362,407	 100.00%
 250 KB 	            41	            41	   358,120,498	   358,120,498	 100.00%
 500 KB 	            36	            36	   356,323,992	   356,323,992	 100.00%
   1 MB 	            35	            35	   355,534,668	   355,534,668	 100.00%
 2.5 MB 	            33	            33	   352,691,277	   352,691,277	 100.00%
   5 MB 	            32	            32	   348,960,151	   348,960,151	 100.00%
  10 MB 	            21	            21	   268,235,216	   268,235,216	 100.00%

A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.3329	0.1672	0.1673	0.3326	0.0000	0.0000	0.0000	0.3345	0.1529

Main genome scaffold total:         	251
Main genome contig total:           	251
Main genome scaffold sequence total:	375.736 MB
Main genome contig sequence total:  	375.736 MB  	0.000% gap
Main genome scaffold N/L50:         	14/12.838 MB
Main genome contig N/L50:           	14/12.838 MB
Main genome scaffold N/L90:         	28/7.327 MB
Main genome contig N/L90:           	28/7.327 MB
Max scaffold length:                	17.553 MB
Max contig length:                  	17.553 MB
Number of scaffolds > 50 KB:        	95
% main genome in scaffolds > 50 KB: 	98.81%


Minimum 	Number        	Number        	Total         	Total         	Scaffold
Scaffold	of            	of            	Scaffold      	Contig        	Contig  
Length  	Scaffolds     	Contigs       	Length        	Length        	Coverage
--------	--------------	--------------	--------------	--------------	--------
    All 	           251	           251	   375,736,066	   375,736,066	 100.00%
   5 KB 	           251	           251	   375,736,066	   375,736,066	 100.00%
  10 KB 	           250	           250	   375,726,492	   375,726,492	 100.00%
  25 KB 	           195	           195	   374,586,004	   374,586,004	 100.00%
  50 KB 	            95	            95	   371,278,513	   371,278,513	 100.00%
 100 KB 	            76	            76	   369,945,567	   369,945,567	 100.00%
 250 KB 	            51	            51	   366,033,390	   366,033,390	 100.00%
 500 KB 	            38	            38	   361,585,998	   361,585,998	 100.00%
   1 MB 	            31	            31	   356,865,229	   356,865,229	 100.00%
 2.5 MB 	            30	            30	   355,254,860	   355,254,860	 100.00%
   5 MB 	            30	            30	   355,254,860	   355,254,860	 100.00%
  10 MB 	            24	            24	   309,940,456	   309,940,456	 100.00%





# ---------- Purge duplicate haplotigs ---------- #
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt

# output files from hifiasm assembly:
lim_hifiasm_genome.asm.bp.hap1.p_ctg.fa  
lim_hifiasm_genome.asm.bp.hap2.p_ctg.fa  
lim_hifiasm_genome.asm.bp.p_ctg.fa

# OUTPUT FILES after filtering hifi adapters
hifiadaptfilt.o799014
hifiadaptfilt.e799014
m64190e_210221_233100.hifi_reads.blocklist
m64190e_210221_233100.hifi_reads.contaminant.blastout
m64190e_210221_233100.hifi_reads.filt.fastq.gz
m64190e_210221_233100.hifi_reads.pbi.stats
m64190e_210221_233100.hifi_reads.stats

# Install software for this step
pwd
/projectnb/mullenl/hannah/limenitis_genomes/software

git clone https://github.com/dfguan/purge_dups.git
cd purge_dups/src && make

git clone https://github.com/dfguan/runner.git
cd runner && python3 setup.py install --user

git clone https://github.com/dfguan/KMC.git 
cd KMC && make -j 16

conda activate Genome

# install in conda environment
conda install -c bioconda purge_haplotigs 

# Preparation for purge_dups - 
# Run minimap2 to align pacbio data and generate paf files
# minimap github page for reference: https://github.com/lh3/minimap2

nano minimap
#!/bin/bash
#$ -V # inherit the submission environment
#$ -N minimap # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -pe omp 10

minimap2 -ax map-hifi lim_hifiasm_genome.asm.bp.p_ctg.fa m64190e_210221_233100.hifi_reads.filt.fastq --secondary=no > aln.sam

# output from minimap job = aln.sam

samtools view -S -b aln.sam > output.bam
# output file = output.bam 

# make and submit file sortbam

cat sortbam
#!/bin/bash
#$ -V # inherit the submission environment
#$ -N sortbam # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -pe omp 16

samtools sort output.bam -m 1G -o aligned.bam -T tmp.ali

qsub sortbam

# create indexed files 
samtools index aligned.bam
samtools faidx lim_hifiasm_genome.asm.bp.p_ctg.fa

# Purge haplotigs Step 1:
purge_haplotigs hist -b aligned.bam -g lim_hifiasm_genome.asm.bp.p_ctg.fa

# output files:
aligned.bam.histogram.png
aligned.bam.gencov

# output instructions from purge haplotigs step 1:
#Check your histogram to observe where your haploid and diploid peaks are
#and choose your low, midpoint, and high cutoffs (check the example histogram png 
#in the readme). You will need 'aligned.bam.gencov' and the cutoffs for the next 
#step: 'purge_haplotigs cov'
# more instructions on the bitbucket site for purge haplotigs here:
#https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Tutorial

# Purge haplotigs Step 2:
# Choose some read-depth cutoffs for low read depth (-l), the midpoint between the 'haploid' and 'diploid' peaks (-m), and high (-h) read-depth cutoff 
purge_haplotigs contigcov -i aligned.bam.gencov -l 12 -m 28 -h 190

# output file = coverage_stats.csv, a contig coverage stats csv file with suspect contigs flagged for further analysis or removal

# Purge haplotigs Step 3:
purge_haplotigs purge -g lim_hifiasm_genome.asm.bp.p_ctg.fa -c coverage_stats.csv


#Now try rerunning, this time to produce dotplots of the final alignments (need to use the -d flag and pass the bam file with the -b flag).
purge_haplotigs purge -g lim_hifiasm_genome.asm.bp.p_ctg.fa -c coverage_stats.csv -d -b aligned.bam

# output curated.fasta is the new curated haploid assembly. More details on output here: https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Tutorial

# index and get curated.fasta ready for mapping
module load bowtie2
module load samtools

bowtie2-build curated.fasta curated.fasta

samtools faidx curated.fasta

# look at number of scaffolds in the new assembly:
qsub assembly_stats_curatedfasta.qsub
cat assembly_stats_curatedfasta.qlog
A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.3359	0.1642	0.1643	0.3356	0.0000	0.0000	0.0000	0.3285	0.1023

Main genome scaffold total:         	64
Main genome contig total:           	64
Main genome scaffold sequence total:	360.605 MB
Main genome contig sequence total:  	360.605 MB  	0.000% gap
Main genome scaffold N/L50:         	13/12.951 MB
Main genome contig N/L50:           	13/12.951 MB
Main genome scaffold N/L90:         	26/8.405 MB
Main genome contig N/L90:           	26/8.405 MB
Max scaffold length:                	17.553 MB
Max contig length:                  	17.553 MB
Number of scaffolds > 50 KB:        	41
% main genome in scaffolds > 50 KB: 	99.82%


Minimum 	Number        	Number        	Total         	Total         	Scaffold
Scaffold	of            	of            	Scaffold      	Contig        	Contig  
Length  	Scaffolds     	Contigs       	Length        	Length        	Coverage
--------	--------------	--------------	--------------	--------------	--------
    All 	            64	            64	   360,605,465	   360,605,465	 100.00%
  10 KB 	            64	            64	   360,605,465	   360,605,465	 100.00%
  25 KB 	            57	            57	   360,466,761	   360,466,761	 100.00%
  50 KB 	            41	            41	   359,946,023	   359,946,023	 100.00%
 100 KB 	            40	            40	   359,886,034	   359,886,034	 100.00%
 250 KB 	            35	            35	   359,099,917	   359,099,917	 100.00%
 500 KB 	            33	            33	   358,462,699	   358,462,699	 100.00%
   1 MB 	            31	            31	   356,865,229	   356,865,229	 100.00%
 2.5 MB 	            30	            30	   355,254,860	   355,254,860	 100.00%
   5 MB 	            30	            30	   355,254,860	   355,254,860	 100.00%
  10 MB 	            24	            24	   309,940,456	   309,940,456	 100.00%

# run BUSCO on this new curated haploid assembly 

# copied needed files from here: /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/Limenitis_IPA_noparam_assembly/BUSCO
# to here: /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/BUSCO
# updated qusb file, busco_genome_assembly_haplotig.qsub:

[haich@scc1 BUSCO]$ cat busco_genome_assembly.qsub
#!/bin/bash -l

#This script was written by Jasmine Alqassar 2023 to run BUSCO v3.0.02 

#$ -P mullenl
#$ -N BUSCO_curated # job name, anything you want
#$ -m bea
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o BUSCO_curated.qlog
#$ -l h_rt=100:00:00
#$ -pe omp 28
#$ -l mem_per_core=13G

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "=========================================================="


module load blast+/2.12.0
module load hmmer/3.3.2
module load python2
module load augustus
module load busco/3.0.2

### format of BUSCO command: busco -i <sequence file> -l <lineage> -o <output name> -m <mode> <other options>
#### For some reason you have to use this python script to run BUSCO, you can't just use the busco command

export AUGUSTUS_CONFIG_PATH="/projectnb/mullenl/hannah/software/augustus_config"
export BUSCO_CONFIG_FILE="/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/BUSCO/config.ini"
python /share/pkg.7/busco/3.0.2/install/scripts/run_BUSCO.py -i /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta -l /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/BUSCO/lepidoptera_odb10 -o BUSCO_Output -m geno -c $NSLOTS


# output:
# BUSCO version is: 3.0.2 
# The lineage dataset is: lepidoptera_odb10 (Creation date: 2020-08-05, number of species: 16, number of BUSCOs: 5286)
# To reproduce this run: python /share/pkg.7/busco/3.0.2/install/scripts/run_BUSCO.py -i /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta -o BUSCO_Output -l /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/BUSCO/lepidoptera_odb10/ -m genome -c 28 -sp fly
#
# Summarized benchmarking in BUSCO notation for file /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta
# BUSCO was run in mode: genome

#	C:98.5%[S:97.7%,D:0.8%],F:0.5%,M:1.0%,n:5286

#	5203	Complete BUSCOs (C)
#	5162	Complete and single-copy BUSCOs (S)
#	41	Complete and duplicated BUSCOs (D)
#	26	Fragmented BUSCOs (F)
#	57	Missing BUSCOs (M)
#	5286	Total BUSCO groups searched


# ---------- Run RepeatModeler / RepeatMasker ---------- #
# This step finds, classifies, and quantifies transposable elements (TEs) in an unannotated genome assembly

# Ue RepeatModeler2 for TE assembly and classification (https://github.com/Dfam-consortium/RepeatModeler, https://stab.st-andrews.ac.uk/wiki/index.php/Repeatmodeler)
	# identifies and classifies repetitive elements de novo and generates a library of consensus repeat sequences for the genome
# At the heart of RepeatModeler are three de-novo repeat finding programs ( RECON, RepeatScout and LtrHarvest/Ltr_retriever ) which employ complementary 
	# computational methods for identifying repeat element boundaries and family relationships from sequence data.
# programs installed on BU SCC already

# we want to use the masked version of the genome for all downstream analysis, because lowercase masking is usually taken into account
# during variant calling or mapping quality calculations (but need to check how ANGSD handles this). It is also possible to manually exclude 
# results from masked regions using bedtools if needed

# Info from Nick VanKuren:
# I usually run RepeatModeler, then combine the output with the RepeatMasker "arthropoda" library. 
# I use that combined fasta file as the custom library for the RepeatMasker run. Sometimes the RepeatModeler output can have some short fragments, 
# so I like to include the curated arthropoda dataset to help make up for that.


# Running RepeatModeler/Masker here:
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/RepeatModeler

# Moved curated.fasta file into the directory above

# First build database and run repeat modeler. Can do this without a qsub

module load perl/5.34.0
module load repeatmodeler/2.0.3

BuildDatabase -name lim_draft curated.fasta

# this made 7 files: lim_draft.nhr  lim_draft.nin  lim_draft.nnd  lim_draft.nni  lim_draft.nog  lim_draft.nsq  lim_draft.translation

#  now run RepeatModeler2 using repeatmodeler.qsub
#!/bin/bash -l
#$ -P mullenl
#$ -N run_repeatmodeler # job name, anything you want
#$ -m bea
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o run_repeatmodeler.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 16

module load perl/5.34.0
module load repeatmodeler/2.0.3

RepeatModeler -pa 16 -engine ncbi -database lim_draft 2>&1 | tee lim_draft.repeatmodeler.log

# directory outputted = RM_471226.ThuJan111048382024
# also output lim_draft-families.fa


# Run repeatmasker.qsub -  Will take ~4-12 hours. ${S}.fa.masked is the final masked fasta.
# ${S}-families.fa is the species-specific repeats identified by RMer.

#!/bin/bash -l
#$ -P mullenl
#$ -N run_repeatmasker # job name, anything you want
#$ -m bea
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o run_repeatmasker.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 16

module load python3/3.10.12
module load repeatmasker/4.1.5

RepeatMasker -pa 16 -gff -engine ncbi -lib lim_draft-families.fa curated.fasta

# output:
# curated.fasta.out 
# curated.fasta.out.gff: gff version of .out for visualisation (IGV...)
# curated.fasta.masked: masked version of the genome (.fasta), with TE residues marked as "X"
# curated.fasta.tbl: summary table; this is great to have a quick estimate of the repeat content. The repeat content is broken down by type of elements. 

# Karina had the output file below, but I did not get a .align file out of my run
# pdiv_genome_curated_host_draft.fa.align: alignment details for the retained RM hits. We will use it in the parsing step

# moved the output from this 'hard masking' to a sub-directory 'hard_masked': /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/RepeatModeler/hard_masked
# going to do a soft-masking too, as suggested by Nick VanKuren and which seems more appropriate for annotation
# more info here: https://star-protocols.cell.com/protocols/2992
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/RepeatModeler/soft_masked

cat repeatmodeler.qsub
#!/bin/bash -l
#$ -P mullenl
#$ -N run_repeatmasker_soft # job name, anything you want
#$ -m bea
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o run_repeatmasker_soft.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 16

module load python3/3.10.12
module load repeatmasker/4.1.5

RepeatMasker -xsmall -gff -pa 16 -engine ncbi -dir soft_masked -lib lim_draft-families.fa curated.fasta


###
# trying to remove repeat regions from vcf file - using a piece of code from a github walk-through I found (https://star-protocols.cell.com/protocols/2992)
# ended up just leaving these repeat regions in the final assembly, leaving here in case it is helpful code.
pwd
/projectnb/mullenl/hannah/limenitis_genomes/software

git clone https://github.com/emmiolkkonen/star_protocols_saimaa.git

# create a tab file of repetitive regions - according to this: https://star-protocols.cell.com/protocols/2992
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/RepeatModeler/soft_masked

rmgff=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/RepeatModeler/soft_masked/curated.fasta.out.gff

awk '!/#/{OFS="∖t";print $1,$4,$5}' $rmgff > repeats.tab

/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/RepeatModeler/hard_masked

rmgff=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/RepeatModeler/hard_masked/curated.fasta.out.gff

awk '!/#/{OFS="∖t";print $1,$4,$5}' $rmgff > repeats.tab

# make bed file 
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' curated.fasta.masked.fai > curated.fasta.masked.bed

# make bed file of repeat regions, the gff output from repeatmasker
module load bedops
gff2bed < curated.fasta.out.gff > curated.fasta.out.bed



# ---------- Annotation with Liftoff ---------- #
# more info here: https://github.com/agshumate/Liftoff
# install into its own conda environment

module load miniconda/23.1.0
conda create --name Genome_Annotation
conda activate Genome_Annotation
conda install -c bioconda liftoff
conda list
# all looks good

# required inputs for Liftoff include the target fasta genome to lift genes to (my reference)
# and a reference fasta genome to lift genes from (using Limenitis camilla): https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_905147385.1/
# using ncbi datasets to download genome data: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli

# now download genome 
datasets download genome accession GCA_905147385.1 --include gff3,rna,cds,protein,genome,seq-report
# this didn't download the annotation files, so I found it easier to download from the tree of life website instead:
# https://projects.ensembl.org/darwin-tree-of-life/

# moved Limenitis camilla files into their own folder on the scc:
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/Limenitis_camilla

conda activate Genome_Annotation
qsub run_liftoff.qsub
cat run_liftoff.qsub
#!/bin/bash -l
#$ -P mullenl
#$ -N liftoff # job name, anything you want
#$ -m bea
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o liftoff.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 16

liftoff -g /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/Limenitis_camilla/Limenitis_camilla-GCA_905147385.1-2021_11-genes.gff3 -o curated.fasta.masked.gff3 /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/RepeatModeler/soft_masked/curated.fasta.masked /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/Limenitis_camilla/Limenitis_camilla-GCA_905147385.1-softmasked.fa


# this annotation didn't work so well, so trying again with the Heliconius melpomene v2.5 genome, downloaded from here:
http://download.lepbase.org/v4/provider/
# also going to try annotation on the not masked version of the genome

cat run_liftoff_hmel.qsub
#!/bin/bash -l
#$ -P mullenl
#$ -N liftoff_hmel2.5 # job name, anything you want
#$ -m bea
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o liftoff_hmel2.5.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 16

liftoff -g /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/Hmel2.5/Hmel2.5.gff3 -o curated.fasta.hmel.gff3 /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/Hmel2.5/Hmel2.5.fa

# This liftoff with Heliconius melpomene worked well, using this for the final annotation.


# ---------- Decorating liftoff annotation with eggnogmapper ---------- #

# It looks like there was a partial CDS in the Liftoff annotation, which apparently happens sometimes with Liftoff . 
# This caused eggnog-mapper to fail to extract proteins from my gff following the method eggnog suggests on their website. 
# Nick VanKuren suggests extracting the proteins myself from the gff3 and the fasta file and then using the protein fasta as input to decorate the gff3 with eggnog
# Will do this with gffread

pwd
/projectnb/mullenl/hannah/scripts

git clone https://github.com/gpertea/gffread
cd gffread
make release

# gffread -x option: write a fasta file with spliced CDS for each GFF transcript
# gffread -y option: write a protein fasta file with the translation of CDS for each record
# for the protein fasta, -V option discards mRNAs with CDS that have an in-frame stop codon - this was messing up eggnog mapper with '.' in the middle of the cds downstream

# with Limenitis camilla liftoff annotation:
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation

original_gff=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/curated.fasta.masked.gff3
input_fasta=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/RepeatModeler/soft_masked/curated.fasta.masked


/projectnb/mullenl/hannah/scripts/gffread/gffread -g ${input_fasta} -x gffreadOutput.fa ${original_gff}
/projectnb/mullenl/hannah/scripts/gffread/gffread -g ${input_fasta} -y gffreadProteinOutput.fa -V ${original_gff}


# with Heliconius melpomene liftoff annotation
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation

original_gff=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/curated.fasta.hmel.gff3
input_fasta=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta

/projectnb/mullenl/hannah/scripts/gffread/gffread -g ${input_fasta} -x gffreadHmelOutput.fa ${original_gff}
/projectnb/mullenl/hannah/scripts/gffread/gffread -g ${input_fasta} -y gffreadHmelProteinOutput.fa -V ${original_gff}


# now we can try using this protein input in eggnog-mapper to decorate our gff3
# Note: i deleted the output from the eggnog decorator script above, but saved the script so i can see what i did previously
module load miniconda
conda activate eggnog-mapper

cat run_eggnog_protein.qsub
#!/bin/bash -l
#$ -P mullenl
#$ -N run_eggnog_protein # job name, anything you want
#$ -m bea
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o run_eggnog_protein.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 10

#Limenitis camilla
#input_protein_fasta=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/gffreadProteinOutput.fa
#out_dir=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/eggnog_output_proteins/
#gff=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/curated.fasta.masked.gff3

#emapper.py --itype proteins -i ${input_protein_fasta} -o eggnog_protein_decorate_liftoff --output_dir ${out_dir} --decorate_gff ${gff} --cpu 10 --override --dmnd_ignore_warnings

#Heliconius melpomene
input_protein_fasta=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/gffreadHmelProteinOutput.fa
out_dir=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/eggnog_output_Hmel_proteins/
gff=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/curated.fasta.hmel.gff3

emapper.py --itype proteins -i ${input_protein_fasta} -o eggnog_Hmel_protein_decorate_liftoff --output_dir ${out_dir} --decorate_gff ${gff} --cpu 10 --override --dmnd_ignore_warnings


cat run_eggnog_genes.qsub
#!/bin/bash -l
#$ -P mullenl
#$ -N run_eggnog_genes # job name, anything you want
#$ -m bea
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o run_eggnog_genes.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 10

#Limenitis camilla
#input_gene_fasta=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/gffreadOutput.fa
#out_dir=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/eggnog_output_genes/
#gff=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/curated.fasta.masked.gff3


#emapper.py -i ${input_gene_fasta} --itype CDS --translate -o eggnog_result_decorate_liftoff --output_dir ${out_dir} --decorate_gff ${gff} --cpu 10 --override --dmnd_ignore_warnings 

#Heliconius melpomene
input_gene_fasta=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/gffreadHmelOutput.fa
out_dir=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/eggnog_output_Hmel_genes/
gff=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/curated.fasta.hmel.gff3


emapper.py -i ${input_gene_fasta} --itype CDS --translate -o eggnog_Hmel_result_decorate_liftoff --output_dir ${out_dir} --decorate_gff ${gff} --cpu 10 --override --dmnd_ignore_warnings 

# Both of these work! Need to assess in IGV but the outputs look much more promising than last time.


# moved all Limenitis camilla outputs into their own folder:
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/eggnog_decoration_limenitiscamilla

# moved all Heliconius melpomene outputs into their own folder:
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/Annotation/eggnog_decoration_heliconiusmelpomene



# ---------- Run MitoHiFi ---------- #
# Identify mitochondrial reads from assembled contigs, more info here: https://github.com/marcelauliano/MitoHiFi?tab=readme-ov-file

# Create singularity container - info on how to do this on SCC here: https://www.bu.edu/tech/support/research/software-and-programming/containers/building/#import
# Login to singularity node
ssh scc-i01

pwd
/usr3/graduate/haich

# fetch the pull command from docker hub for mitohifi
singularity pull docker://ghcr.io/marcelauliano/mitohifi:master

# this created an output called mitohifi_master.sif. moved this container here:
mv mitohifi_master.sif /projectnb/mullenl/hannah/software/

# this is how we interact with the singularity container:
singularity exec /projectnb/mullenl/hannah/software/mitohifi_master.sif mitohifi.py -h

# download mitochondrial genome from closely related species
singularity exec /projectnb/mullenl/hannah/software/mitohifi_master.sif findMitoReference.py --species "Limenitis elwesi" --outfolder ./mitogenome_ref/ --min_length 14000


# had to move curated.fasta file to working directory for the job to run, but the below job worked when I ran it here: /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/MitoHiFi
# also request interactive node first 
qrsh -l h_rt=12:00:00 -pe omp 4

my_fasta=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/MitoHiFi/curated.fasta
relatedmito_fasta=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/MitoHiFi/mitogenome_ref/NC_039867.1.fasta
relatedmitogb=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/MitoHiFi/mitogenome_ref/NC_039867.1.gb

singularity exec /projectnb/mullenl/hannah/software/mitohifi_master.sif mitohifi.py -c ${my_fasta} -f ${relatedmito_fasta} -g ${relatedmitogb} -t 4 -o 5 

# Main output files:
final_mitogenome.fasta - the final mitochondria circularized and rotated to start at tRNA-Phe
final_mitogenome.gb - the final mitochondria annotated in GenBank format.
final_mitogenome.coverage.png - the sequencing coverage throughout the final mitogenome
final_mitogenome.annotation.png - the predicted genes throughout the final mitogenome
contigs_annotations.png - annotation plots for all potential contigs
coverage_plot.png - reads coverage plot of filtered reads mapped to all potential contigs
contigs_stats.tsv - containing the statistics of your assembled mitos such as the number of genes, size, whether it was circularized or not, if the sequence has frameshifts and etc...
shared_genes.tsv - show comparison of annotation between close-related mitogenome and all potential contigs assembled


# ---------- Satsuma2 ---------- #

#### Trying chromosembler (Satsuma2) to scaffold the assembly using Synteny
# have not completed this as of November 2024, not sure if it's necessary?
module load miniconda
conda create --name Satsuma
conda activate Satsuma
conda install bioconda::satsuma2

# to map scaffolds onto chromosome coordinates via synteny run Chromosemble (https://satsuma.sourceforge.net/), options:
Chromosemble -t <reference> -q <your_scaffolds> -o <output_dir>

pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/satsuma2

Chromosemble -q /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta -t 







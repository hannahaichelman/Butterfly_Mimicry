# ---------- Sample Summary ---------- #

## Summary of files used for this GWAS
# This information based on talking with Sean on 6/23/23:

# Genomic cline fastq files (32 individuals, each with R1 and R2 reads across two lanes):
# Used Illumina Nextera adapters, sequenced on Illumina HiSeq2500
# Adapter trimming used for R1 and R2 reads: CTGTCTCTTATACACATCT
# Found here:
pwd
/projectnb/mullenl/evan/sequence_data/genomic_cline_sequence/assembly

# Confirmed that genomic cline samples were sequenced on HiSeq2500, the header that designates sequencer starts with @SN which corresponds to Hiseq2500 machines (https://www.biostars.org/p/198143/)
[haich@scc1 genomic_cline]$ head -20 LIB035248_GEN00127679_S1_R1.fastq
@SN363:1241:HKGC3BCX2:1:1107:1481:2179 1:N:0:TAGGCATG+CGTCTAAT
CACTACTACAAGAATAAAAAAAAATTGATCTGCCATCACCCATAAAACTAAAGTTAGGTGGATATTTCATTTTTTTTAATTAGCAGAGGTGGTTTTTAACTTAATTAACACTTACCGAATTTAAAAGAATTTCTCAACCAACACAAAAAA
+
DDDDDIIIIIIIIIIIIIIIIIIIHIIIIIIIIIIIIIIIIIIIIIIIIIHIIHIIIIHHIHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIGHIHIIIIIIIHHIIIIIIIIIIHIIIIIIIHIIIIIIIIHIIIIIIII

# There are additional whole genomes that Dana Cohen originally sequenced that we could also use for the GWAS
# This sequencing was across 3 lanes of the Illumina HiSeq 2500
# 
pwd
/projectnb/mullenl/danacohen

[haich@scc1 dana_cohen_files]$ head RIH_2088_R1.fq
@D6FSQ5P1:553:HCGWYADXX:1:1101:2137:2194 1:N:0:TAAGGCGATAGATCGC
TATTAAAATAAAATGTCTAAAAGGCCTCAATGAAGTACCTGTTTCATTCTTCATCAGGAAAATCAGCTTGCATTCTTGCCACCATTCCTAACGGGTGGGCCCTATAGCTGTCTCTTATACACATCTCCGAGCCCACGAGACTATGGCGAA
+
??@BDB???F?DDF<2<?>FH>>AFFFF9A@EAEEGEGCDEFH@GCGHIFGE?9B?<3B@9D:B3=@<=8=@CG>:@@;=@7(.=C:).;;>A:AB,3;(95<9>:+4::>:>CC::::@3>>@CCCA59.5581+)5&00:4(++(0.&


# Limenitis mapping brood (QTL) fastq files (96 individuals, each with R1 and R2 reads):
# Sequenced on Illumina NovaSeq S4 PE 150 Cycle - feel pretty sure about this 
# Found here:
pwd
/projectnb/mullenl/evan/sequence_data/limenitis_resequencing

# Confirmed sequenced on NovaSeq, matches what is described at the biostars response linked above
[haich@scc1 limenitis_resequencing]$ head 15-1_S1_R1_001.fastq
@A00572:341:HLWJTDSX3:2:1101:1108:1000 1:N:0:TTGTTCCTTG+TGTTGTTCGG
TNCGCAAATTGTGCCACCCCCGTTAACCCCCATAGGGGTTGATTTTCAAAAAACGCTTAGATATGTATCTCATAATTTGTATTGCAGAGCATTCTTTCAGTGTCACGCAAAATATTTCAATAGAATTAAATTTCCCATACAAACTTAGCCA
+
F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFFFFFFFFFFFFFFFFFFFFF,FFF,FF:FF

# Evan has an example of adapter trimming for all of the 60 GWAS files together:
[haich@scc1 all_individuals]$ pwd
/projectnb/mullenl/evan/sequence_data/GWAS/GWAS_prep/all_individuals

Command line parameters: -a CTGTCTCTTATACACATCT ./Run_2-LIB013265_GEN00032321_AGGCAGAA-CTCTCTAT_L001_R1.fastq.gz -o /projectnb/mullenl/evan/sequence_data/Evan_GWAS/trimmed/Run_2-LIB013265_GEN00032321_AGGCAGAA-CTCTCTAT_L001_R1_trimmed.fastq.gz


# take 2 samples, gunzip them, look at # of reads in r1 and r2 before and after trimming

# ---------- Organize Files ---------- #
# Moved files from Evan's folder and Dana Cohen's folder into my own mullenlab folder. 
pwd
/projectnb/mullenl/hannah/limenitis_genomes

# Three sub-folders are in this directory that contain samples for GWAS: 
# one for the genomic cline samples (/projectnb/mullenl/evan/sequence_data/genomic_cline_sequence/assembly) and 
# one for the limenitis_resequencing files(/projectnb/mullenl/evan/sequence_data/limenitis_resequencing)
# one for the extra files from dana cohen, already concatenated (/projectnb/mullenl/hannah/limenitis_genomes/dana_cohen_files)

# Concatenate the genomic cline files (combine L001 and L002 files)
# made job file called cat_files
pwd
/projectnb/mullenl/hannah/limenitis_genomes/genomic_cline

head cat_files
cat LIB035248_GEN00127679_S1_L001_R1.fastq.gz LIB035248_GEN00127679_S1_L002_R1.fastq.gz > LIB035248_GEN00127679_S1_R1.fastq.gz
cat LIB035248_GEN00127679_S1_L001_R2.fastq.gz LIB035248_GEN00127679_S1_L002_R2.fastq.gz > LIB035248_GEN00127679_S1_R2.fastq.gz
cat LIB035248_GEN00127680_S2_L001_R1.fastq.gz LIB035248_GEN00127680_S2_L002_R1.fastq.gz > LIB035248_GEN00127680_S2_R1.fastq.gz
cat LIB035248_GEN00127680_S2_L001_R2.fastq.gz LIB035248_GEN00127680_S2_L002_R2.fastq.gz > LIB035248_GEN00127680_S2_R2.fastq.gz
cat LIB035248_GEN00127681_S3_L001_R1.fastq.gz LIB035248_GEN00127681_S3_L002_R1.fastq.gz > LIB035248_GEN00127681_S3_R1.fastq.gz

../scc6_qsub_launcher.py -N concatenate_files -P davieslab -M hannahaichelman@gmail.com -j y -h_rt 24:00:00 -jobsfile cat_files
qsub concatenate_files_array.qsub

# The concatenate job worked, moved the original (un-concatenated) files into a folder called 'original_files'
pwd
/projectnb/mullenl/hannah/limenitis_genomes/genomic_cline/original_files

# Dana cohen files need to be unzipped first, they came as a .bz2 compressed file and it seems to be interfering with the trim job
cat bzip_files
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N bzip # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

bzip2 -dk *.bz2

# ---------- Trim Files ---------- #
# using fastp to detect adapters and remove them

# Based on this paper (https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13559) and associated Github (https://github.com/therkildsen-lab/batch-effect) 
# we need to trim poly-G tails using a sliding window approach from samples sequenced on new machines (NovaSeq/NextSeq).
# This helped fix batch effects in the paper referenced above, especially when the samples are sequenced on machines using different chemistry
# Here, we have HiSeq data, which is a four-channel system, and also NovaSeq data, which is a two-channel system. The two-channel system in particular
# is susceptible to poly-G tails, which need to be removed. 
# I can see these poly-G tails on the limenitis resequencing data, but not in the other HiSeq data.

# So, really just need to do this for the new limenitis_resequencing files, which were on the NovaSeq:
pwd
/projectnb/mullenl/hannah/limenitis_genomes/limenitis_resequencing

# First, detect and remove adapters
>trim_test
for i in *R1_001.fastq.gz; do 
echo -e "fastp --detect_adapter_for_pe -q 20 \
-i $i -I ${i/R1_001.fastq.gz/}R2_001.fastq.gz -o ${i/_R1_001.fastq.gz/}_trimtest_R1.fastq.gz -O ${i/_R1_001.fastq.gz/}_trimtest_R2.fastq.gz" >>trim_test;
done
# just submitted this as a normal qsub job, with -pe omp 8
# the trim_test.e9095204 output shows that I was trimming the wrong adapter - it actually has TruSeq adapters...lol. 
# still seeing lots of poly-G ends of the reads though, so now move to do that trimming per the Therkildsen paper. 

# sliding window approach detailed in Therkildsen paper:
>trim2
for i in *trimtest_R1.fastq.gz; do 
echo -e "fastp --trim_poly_g -L -A --cut_right \
-i $i -I ${i/trimtest_R1.fastq.gz/}trimtest_R2.fastq.gz -o ${i/trimtest_R1.fastq.gz/}sliding_window_trim_R1.fastq.gz -O ${i/trimtest_R1.fastq.gz/}sliding_window_trim_R2.fastq.gz" >>trim2;
done
module load fastp
../scc6_qsub_launcher.py -N sliding_window_trim -P mullenl -M hannahaichelman@gmail.com -m as -omp 4 -j y -h_rt 24:00:00 -jobsfile trim2
qsub sliding_window_trim_array.qsub

# Moved all outputs into sub-folders called 'sliding_window_trim_output' 
/projectnb/mullenl/hannah/limenitis_genomes/limenitis_resequencing/sliding_window_trim_output

# deleted the intermediate *trimtest* files

# Now trim the genomic cline samples. Have to slightly adapt the code to make the trim job b/c of different file name scheme
pwd
/projectnb/mullenl/hannah/limenitis_genomes/genomic_cline

>trim_test
for i in *R1.fastq; do 
echo -e "fastp --detect_adapter_for_pe -q 20 \
-i $i -I ${i/R1.fastq/}R2.fastq -o ${i/R1.fastq/}trimtest_R1.fastq -O ${i/R1.fastq/}trimtest_R2.fastq" >>trim_test;
done

qsub trim_test # after adding header including -pe omp 8


# Now on to the dana cohen files:
pwd
/projectnb/mullenl/hannah/limenitis_genomes/dana_cohen_files

>trim_test
for i in *R1.fq; do 
echo -e "fastp --detect_adapter_for_pe -q 20 \
-i $i -I ${i/R1.fq/}R2.fq -o ${i/R1.fq/}trimtest_R1.fq -O ${i/R1.fq/}trimtest_R2.fq" >>trim_test;
done

qsub trim_test # after adding header including -pe omp 8




# ---------- Mapping with bwa mem ---------- #

# after new genome assembly
# mapping reads to my purge_dups assembly (curated.fasta)
# decided to stick with bwa-mem after testing bowtie2, both worked but had higher mapping efficiencies with bwa. 

# first, index reference file for bwa
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt

module load bwa
bwa index curated.fasta

# map dana cohen resequencing reads
pwd
/projectnb/mullenl/hannah/limenitis_genomes/dana_cohen_files
export GENOME_BASENAME=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta
export NSLOTS=16

>map_curatedfasta_dana
for i in *trimtest_R1.fq; do
echo -e "bwa mem $GENOME_BASENAME $i ${i/_R1.fq}_R2.fq >$(basename "$i" _R1.fq).sam" >>map_curatedfasta_dana;
done

../scc6_qsub_launcher.py -N map_curatedfasta_dana -P mullenl -M hannahaichelman@gmail.com -m a -omp $NSLOTS -j y -h_rt 72:00:00 -jobsfile map_curatedfasta_dana
qsub map_curatedfasta_dana_array.qsub


# map genomic cline reads
pwd
/projectnb/mullenl/hannah/limenitis_genomes/genomic_cline

>map_curatedfasta_gc
for i in *trimtest_R1.fastq; do
echo -e "bwa mem $GENOME_BASENAME $i ${i/_R1.fastq}_R2.fastq >$(basename "$i" _R1.fastq).sam" >>map_curatedfasta_gc;
done

../scc6_qsub_launcher.py -N map_curatedfasta_gc -P mullenl -M hannahaichelman@gmail.com -m a -omp $NSLOTS -j y -h_rt 72:00:00 -jobsfile map_curatedfasta_gc
qsub map_curatedfasta_gc_array.qsub

# map limenitis resequencing reads

>map_curatedfasta_lim
for i in *sliding_window_trim_R1.fastq.gz; do
echo -e "bwa mem $GENOME_BASENAME $i ${i/_R1.fastq.gz}_R2.fastq.gz >$(basename "$i" _R1.fastq.gz).sam" >>map_curatedfasta_lim;
done

../scc6_qsub_launcher.py -N map_curatedfasta_lim -P mullenl -M hannahaichelman@gmail.com -m a -omp $NSLOTS -j y -h_rt 72:00:00 -jobsfile map_curatedfasta_lim
qsub map_curatedfasta_lim_array.qsub


# ---------- Make BAM fiiles ---------- #
# samtools -F 4 excludes unmapped reads
# samtools -q 20 skips alignments with MAPQ less than 20

export NSLOTS=16
module load samtools

# ran the following steps in each of the 3 folders
pwd
/projectnb/mullenl/hannah/limenitis_genomes/dana_cohen_files & genomic_cline & limenitis_resequencing

>sam2bam
for i in *.sam; do 
echo -e "samtools view -bS -F 4 $i > ${i/.sam/}.bam \
&& samtools view -h -q 20 ${i/.sam/}.bam | samtools view -buS - | samtools sort -o ${i/.sam/}_minq20_sorted.bam \
" >>sam2bam; done
../../scc6_qsub_launcher.py -N sam2bam_lim -P mullenl -M hannahaichelman@gmail.com -m a -omp $NSLOTS -j y -h_rt 24:00:00 -jobsfile sam2bam
qsub sam2bam_lim_array.qsub


# add in a step to remove pcr duplicates from BAM files - James suggested this to me on April 8th 2024
# because i'd already started calling SNPs, did this in the following folder: /projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe
# more info on this MarkDuplicates here: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard
# and here: https://hpc.nih.gov/training/gatk_tutorial/markdup.html

>removedups
for i in *minq20_sorted.bam; do
echo -e "picard MarkDuplicates I=$i O=${i/.bam}.dedup.bam M=${i/.bam}.dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true TMP_DIR=`pwd`/tmp" >> removedups; done

# added header to removedups qsub job

#!/bin/bash -l
#$ -P mullenl
#$ -N removedups # job name, anything you want
#$ -m bea
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o removedups.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 10
module load picard


## now calculate missing data and depth of coverage info from BAM files
# running with the dedup'd files
ls *sorted.dedup.bam > bamlist

export BAMLIST=bamlist
>finddepth_dedup
for SAMPLEBAM in `cat $BAMLIST`; do 
echo -e "samtools depth -aa $SAMPLEBAM | cut -f 3 | gzip > $SAMPLEBAM.depth.gz \
" >> finddepth_dedup; done
# added header below and qsub 'finddepth' job

#!/bin/bash -l
#$ -P mullenl
#$ -N finddepth_dedup # job name, anything you want
#$ -m bea
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o finddepth_dedup.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 4 
module load samtools

# create list of depth files
ls *depth.gz > depthlist

# Run the R scripts in each sub-folder to calculate depth of coverage. The files are:
#/projectnb/mullenl/hannah/limenitis_genomes/limenitis_resequencing/calcdepth_lim.R
#/projectnb/mullenl/hannah/limenitis_genomes/dana_cohen_files/calcdepth_dc.R
#/projectnb/mullenl/hannah/limenitis_genomes/genomic_cline/calcdepth_gc.R

# I use the csv's output from these R scripts to figure out who to toss and who to keep 

# Made a new folder where I moved all bam files from each of the 3 sub-folders with coverage > 0.5X (mean_depth > 0.5):

pwd
/projectnb/mullenl/hannah/limenitis_genomes/snp_calling

# individuals that were removed: - these stay the same whether mapping to the masked genome or not, and whether using bwamem or bowtie2
RIH_2399_trimtest_minq20_sorted.bam.depth.gz 
F2_HSC_7-19_S82_sliding_window_trim_minq20_sorted.bam.depth.gz 
BCB4-25_S62_sliding_window_trim_minq20_sorted.bam.depth.gz 
BCB4-16_S55_sliding_window_trim_minq20_sorted.bam.depth.gz 
BCB4-13_S53_sliding_window_trim_minq20_sorted.bam.depth.gz 
BCB2-1_S37_sliding_window_trim_minq20_sorted.bam.depth.gz 

# when running only wild caught samples only removed:
RIH_2399_trimtest_minq20_sorted.bam.depth.gz 

# have a total of 142 individuals at this step
# leaves 51 individuals  for wild caught (dana cohen + genomic cline), including 32 genomic cline and 19 dana cohen samples
# leaves 91 samples for mapping brood (limenitis resequencing)

# index bam files to .bai files - this helped investigate if there were issues with the bam files
>indexbam
for i in *.bam; do 
echo -e "samtools index $i ${i/.bam}.bam.bai" >>indexbam; done
../scc6_qsub_launcher.py -N indexbam -P mullenl -M hannahaichelman@gmail.com -m a -omp 8 -j y -h_rt 24:00:00 -jobsfile indexbam
qsub indexbam_array.qsub


# ---------- Run ANGSD ---------- #
pwd
/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe

ls *.dedup.bam > bams


# when I tried re-running angsd with only wild-caught individuals, we ran into a bunch of trouble. where we would get 0 SNPs regardless of changing quality filters. 
# after messing with this for what felt like forever, James Fifer found a thread that indicated snpStat filters, and when I removed those I get SNPs again. 
# Here is the thread: https://github.com/ANGSD/angsd/issues/581

# with the updates regarding the snpstat filters, now run angsd with all, wild-caught individuals, and the mapping brood individuals separately:
# excluding setMinDepth after conversation with James

cat angsd_all
#!/bin/bash
#$ -V # inherit the submission environmentm
#$ -cwd # start job in submission directory
#$ -N angsd_allsamps  # job name, anything you want
#$ -P mullenl
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -j y # Join standard output and error to a single file
#$ -o angsd_allsamps.qlog
#$ -m be
#$ -pe omp 28

module load angsd
GENOME_FASTA=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta # for original un-masked genome

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -skipTriallelic 1 -minInd 114 -ref $GENOME_FASTA -minMaf 0.05"
TODO="-doMajorMinor 1 -dobcf 1 -doMaf 1 -doCounts 1 -doHWE 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 4 -doPost 1 -doGlf 2"

angsd -b bams -GL 1 -P 28 $FILTERS $TODO -out AllLoci_allsamps

NSITES=`zcat AllLoci_allsamps.mafs.gz | wc -l`
echo $NSITES
15775104


cat angsd_wild
#!/bin/bash
#$ -V # inherit the submission environmentm
#$ -cwd # start job in submission directory
#$ -N angsd_wild  # job name, anything you want
#$ -P mullenl
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -j y # Join standard output and error to a single file
#$ -o angsd_wild.qlog
#$ -m be
#$ -pe omp 28

module load angsd
GENOME_FASTA=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta # for original un-masked genome

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -skipTriallelic 1 -minInd 41 -ref $GENOME_FASTA -minMaf 0.05"
TODO="-doMajorMinor 1 -dobcf 1 -doMaf 1 -doCounts 1 -doHWE 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 4 -doPost 1 -doGlf 2"

angsd -b bams_wild -GL 1 -P 28 $FILTERS $TODO -out AllLoci_wild

# before removing maxhetfreq filter:
NSITES=`zcat AllLoci_wild.mafs.gz | wc -l`
echo $NSITES
9921954

cat angsd_mappingbrood
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N angsd_mappingbrood  # job name, anything you want
#$ -P mullenl
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -j y # Join standard output and error to a single file
#$ -o angsd_mappingbrood.qlog
#$ -m be
#$ -pe omp 28

module load angsd
GENOME_FASTA=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta # for original un-masked genome
#GENOME_FASTA=/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/RepeatModeler/hard_masked/curated.fasta.masked # for hard masked genome

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -skipTriallelic 1 -minInd 73 -ref $GENOME_FASTA -minMaf 0.05"
TODO="-doMajorMinor 1 -dobcf 1 -doMaf 1 -doCounts 1 -doHWE 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 4 -doPost 1 -doGlf 2"

angsd -b bams_mappingbrood -GL 1 -P 28 $FILTERS $TODO -out AllLoci_mappingbrood

# before removing maxhetfreq filter:
NSITES=`zcat AllLoci_mappingbrood.mafs.gz | wc -l`
echo $NSITES
16342443



### Here is some helpful troubleshooting I did in trying to figure out the 0 SNP problem:

# identify missing data from samples and by sites using plink2: https://www.cog-genomics.org/plink/1.9/basic_stats#missing
pwd
/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/allsamples_plus_subsetfromallsamples

cat missing_data_bwa
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N missing_data_bwa  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -pe omp 16
#$ -j y # Join standard output and error to a single file
#$ -o missing_data_bwa.qlog

module load plink2

plink2 --allow-extra-chr --double-id --vcf AllLoci_bwa_annotated_output.vcf --missing

# did the same for the bowtie2 mapped files. outputs two files:
#--missing: Sample missing data report written to plink2.smiss .
#--missing: Variant missing data report written to plink2.vmiss .
# re-named to include mapping info and copied to local computer for plotting


# ---------- Make and Annotate VCF ---------- #

#convert bcf to vcf
# put this into a job:
pwd
/projectnb/mullenl/hannah/limenitis_genomes/snp_calling

cat bcf2vcf
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N bcf2vcf  # job name, anything you want
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -pe omp 8

module load htslib/1.16
module load bcftools/1.16

bcftools view AllLoci_allsamps.bcf > AllLoci_allsamps.vcf
bcftools view AllLoci_mappingbrood.bcf > AllLoci_mappingbrood.vcf
bcftools view AllLoci_wild.bcf > AllLoci_wild.vcf


# annotate vcf to give SNPs names
cat annotate_vcf_all
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N annotate_vcf_all  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

module load htslib/1.16
module load bcftools/1.9

#need to annotate vcf file and give SNP names otherwise downstream SNP IDs are all "."
awk '/^[^#]/ {print $1 "\t" $2}' AllLoci_allsamps.vcf > AllLoci_allsamps_chrpos.tmp1
awk '/^[^#]/ {print $1"."$2}' AllLoci_allsamps.vcf > AllLoci_allsamps_chrpos.tmp2
paste AllLoci_allsamps_chrpos.tmp1 AllLoci_allsamps_chrpos.tmp2 > AllLoci_allsamps_chrpos
sed -i '1i\CHROM\tPOS\tANNOTATION_DATA' AllLoci_allsamps_chrpos

bgzip -c AllLoci_allsamps.vcf > AllLoci_allsamps.vcf.gz
bcftools annotate --set-id '%CHROM\_%POS' AllLoci_allsamps.vcf.gz > AllLoci_allsamps_annotated_output.vcf

cat annotate_vcf_wild
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N annotate_vcf_wild  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

module load htslib/1.16
module load bcftools/1.9

#need to annotate vcf file and give SNP names otherwise downstream SNP IDs are all "."
awk '/^[^#]/ {print $1 "\t" $2}' AllLoci_wild.vcf > AllLoci_wild_chrpos.tmp1
awk '/^[^#]/ {print $1"."$2}' AllLoci_wild.vcf > AllLoci_wild_chrpos.tmp2
paste AllLoci_wild_chrpos.tmp1 AllLoci_wild_chrpos.tmp2 > AllLoci_wild_chrpos
sed -i '1i\CHROM\tPOS\tANNOTATION_DATA' AllLoci_wild_chrpos

bgzip -c AllLoci_wild.vcf > AllLoci_wild.vcf.gz
bcftools annotate --set-id '%CHROM\_%POS' AllLoci_wild.vcf.gz > AllLoci_wild_annotated_output.vcf

cat annotate_vcf_mapping
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N annotate_vcf_mapping  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be

module load htslib/1.16
module load bcftools/1.9

awk '/^[^#]/ {print $1 "\t" $2}' AllLoci_mappingbrood.vcf > AllLoci_mappingbrood_chrpos.tmp1
awk '/^[^#]/ {print $1"."$2}' AllLoci_mappingbrood.vcf > AllLoci_mappingbrood_chrpos.tmp2
paste AllLoci_mappingbrood_chrpos.tmp1 AllLoci_mappingbrood_chrpos.tmp2 > AllLoci_mappingbrood_chrpos
sed -i '1i\CHROM\tPOS\tANNOTATION_DATA' AllLoci_mappingbrood_chrpos

bgzip -c AllLoci_mappingbrood.vcf > AllLoci_mappingbrood.vcf.gz
bcftools annotate --set-id '%CHROM\_%POS' AllLoci_mappingbrood.vcf.gz > AllLoci_mappingbrood_annotated_output.vcf


# ---------- Do GWAS using GEMMA ---------- #
pwd
/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe

# First LD Prune SNPs
# Use plink_LD job to create a pruned vcf from the loci in approx LD
# James Fifer is removing sites that are not in linkage equilibrium (only to create the relatedness matrix!) because sites with high LD might cause spurious signals of relatedness.
# Example from James: If individual A has some gene variant X at locus 1 and individual B has some gene variant Y at locus 1 and locus 1 is in high LD with 100 other loci 
# you would expect that these two individuals also show divergence at those 100 other loci, but for these 100 other loci its not necessarily a reflection of how non 
# related they are its just a reflection of your loci being being skewed towards those with high LD. In other words if you dont prune out loci in high LD you might 
# accidentally underestimate (or overestimate if you took this same example but now both individuals have gene variant Y) relatedness 

# we only need to remove sites not in LD to generate the relatedness matrix, we can keep all the sites for the gwas with gemma.

cat plink_LD
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N pink_LD  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o plink_LD.qlog


module load plink2
plink --double-id --allow-extra-chr --geno 0.10 --indep-pairwise 200 20 .2 --maf 0.05 --out pruned_mappingbrood --vcf AllLoci_mappingbrood_annotated_output.vcf
plink --double-id --allow-extra-chr --geno 0.10 --indep-pairwise 200 20 .2 --maf 0.05 --out pruned_wild --vcf AllLoci_wild_annotated_output.vcf
plink --double-id --allow-extra-chr --geno 0.10 --indep-pairwise 200 20 .2 --maf 0.05 --out pruned_allsamps --vcf AllLoci_allsamps_annotated_output.vcf

# plink_LD outputs pruned.prune.in and pruned.prune.out, which are a list of sites in linkage equilibrium (pruned.prune.in) and out of linkage equilibrium (pruned.prune.out)

# now prune the vcf file to include only sites that are in linkage equilibrium
cat prune_vcf
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N prune_vcf  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o prune_vcf.qlog

module load vcftools
vcftools --vcf AllLoci_mappingbrood_annotated_output.vcf --snps pruned_mappingbrood.prune.in --recode --out pruned_mappingbrood
vcftools --vcf AllLoci_wild_annotated_output.vcf --snps pruned_wild.prune.in --recode --out pruned_wild
vcftools --vcf AllLoci_allsamps_annotated_output.vcf --snps pruned_allsamps.prune.in --recode --out pruned_allsamps

# prune_vcf outputs pruned.recode.vcf (vcf file with the sites not in linkage equilibrium removed)


# now create the mean genotype files needed for gemma (BIMBAM format)
# the pruned file will be used to make the relatedness matrix, but the un-pruned file will be used as input for our gwas in gemma
cat make_bimbam
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N make_bimbam  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o make_bimbam.qlog


module load qctool/2.0.1 
qctool -g pruned_wild.recode.vcf -ofiletype bimbam_dosage -og AllLoci_wild.pruned.bimbam.geno
qctool -g AllLoci_wild_annotated_output.vcf -ofiletype bimbam_dosage -og AllLoci_wild.bimbam.geno

qctool -g pruned_mappingbrood.recode.vcf -ofiletype bimbam_dosage -og AllLoci_mappingbrood.pruned.bimbam.geno
qctool -g AllLoci_mappingbrood_annotated_output.vcf -ofiletype bimbam_dosage -og AllLoci_mappingbrood.bimbam.geno

qctool -g pruned_allsamps.recode.vcf -ofiletype bimbam_dosage -og AllLoci_allsamps.pruned.bimbam.geno
qctool -g AllLoci_allsamps_annotated_output.vcf -ofiletype bimbam_dosage -og AllLoci_allsamps.bimbam.geno



# now create the phenotype file (BIMBAM format gemma needs) from trait file - this code just takes the trait values and gets rid of header
# since we have a binary trait, 0 indicates control (black wing) and 1 indicates case (orange wing)
cut -f7 -d ',' GWAS_Phenotypes_LowCoverageRemoved.csv | sed '1d' > AllLoci_allsamps.gemma.phenotype
cut -f7 -d ',' GWAS_Phenotypes_LowCoverageRemoved_WildCaughtOnly.csv | sed '1d' > AllLoci_wild.gemma.phenotype
cut -f7 -d ',' GWAS_Phenotypes_LowCoverageRemoved_MappingBrood.csv | sed '1d' > AllLoci_mappingbrood.gemma.phenotype


# Run GEMMA to create relatedness file
# Helpful notes on gemma: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://romainvilloutreix.alwaysdata.net/romainvilloutreix/wp-content/uploads/2017/01/gwas_gemma-2017-01-17.pdf
# Note from gemma manual: 
# in our previous experience based on limited examples, we typically find the centered genotype (-gk 1) matrix provides better control for population structure in lower organisms
# note this relatedness matrix should be made with the pruned genotype file
module load miniconda
conda activate Genome

cat make_relatedness
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N make_relatedness  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o make_relatedness.qlog

gemma -p AllLoci_allsamps.gemma.phenotype -g AllLoci_allsamps.pruned.bimbam.geno -gk 1 -o AllLoci_allsamps.pruned.relatedness -maf 0.05 -hwe 1e-7 -miss 0.10
gemma -p AllLoci_wild.gemma.phenotype -g AllLoci_wild.pruned.bimbam.geno -gk 1 -o AllLoci_wild.pruned.relatedness -maf 0.05 -hwe 1e-7 -miss 0.10
gemma -p AllLoci_mappingbrood.gemma.phenotype -g AllLoci_mappingbrood.pruned.bimbam.geno -gk 1 -o AllLoci_mappingbrood.pruned.relatedness -maf 0.05 -hwe 1e-7 -miss 0.10


# the make_relatedness job automatically creates a folder called 'output' and puts the output in there. 
# the output includes prefix.log.txt file with details about the running parameters, along with a prefix.cXX.txt file that contains an nxn matrix of estimated relatedness

# now run lmm in gemma
# -a (optional) specifies BIMBAM SNP annotation file name

#first have to make BIMBAM SNP annotation file
awk -F'[: ]' '{print $1 ",", $3 ",", $2}' AllLoci_allsamps.bimbam.geno > AllLoci_allsamps.bimbam.snps
awk -F'[: ]' '{print $1 ",", $3 ",", $2}' AllLoci_wild.bimbam.geno > AllLoci_wild.bimbam.snps
awk -F'[: ]' '{print $1 ",", $3 ",", $2}' AllLoci_mappingbrood.bimbam.geno > AllLoci_mappingbrood.bimbam.snps


#running gemma to do a linear model association test
#-lmm [num] option specifies which frequentist test to use, i.e. “-lmm1” performs Wald test,“-lmm2” performs likelihood ratio test,“-lmm3” performs score test, and 
#“-lmm4” performs all three tests
module load miniconda
conda activate Genome
cat run_gemma
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N run_gemma  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o run_gemma.qlog

gemma -p AllLoci_allsamps.gemma.phenotype -g AllLoci_allsamps.bimbam.geno -a AllLoci_allsamps.bimbam.snps -k ./output/AllLoci_allsamps.pruned.relatedness.cXX.txt -lmm 4 -o Butterfly_allsamps.lmm -maf 0.05 -miss 0.10
gemma -p AllLoci_wild.gemma.phenotype -g AllLoci_wild.bimbam.geno -a AllLoci_wild.bimbam.snps -k ./output/AllLoci_wild.pruned.relatedness.cXX.txt -lmm 4 -o Butterfly_wild.lmm -maf 0.05 -miss 0.10
gemma -p AllLoci_mappingbrood.gemma.phenotype -g AllLoci_mappingbrood.bimbam.geno -a AllLoci_mappingbrood.bimbam.snps -k ./output/AllLoci_mappingbrood.pruned.relatedness.cXX.txt -lmm 4 -o Butterfly_mappingbrood.lmm -maf 0.05 -miss 0.10



# use the plot_gemma_results.R script to make a manhattan plot with these results


# info on why we can use relatedness matrix and not include more info to control for population structure: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4102448/
# STAR protocol that had useful information for masking of the genome: https://star-protocols.cell.com/protocols/2992

# ---------- Additional Downstream Analyses ---------- #

#------------------------------Fst per-SNP
# helpful info here: https://speciationgenomics.github.io/per_site_Fst/

mkdir genome_scan
pwd
/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/genome_scan

VCF_wild=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_wild_annotated_output.vcf
VCF_mapping=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_mappingbrood_annotated_output.vcf

# made population files with bam files associated with lorquini (lorquini_wild), weidemeyerii (weidemeyerii_wild), and hybrids
# Sean says we care more about Fst between orange and black phenotypes, so also made those files and re-ran that comparison. 

cat fst_pop
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N fst  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o fst.qlog

module load vcftools

VCF_wild=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_wild_annotated_output.vcf
VCF_mapping=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_mappingbrood_annotated_output.vcf
VCF_all=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_allsamps_annotated_output.vcf

vcftools --vcf ${VCF_wild} --weir-fst-pop orange_wild2 --weir-fst-pop black_wild2 --out ./orange_black_wild
vcftools --vcf ${VCF_mapping} --weir-fst-pop orange_mapping2 --weir-fst-pop black_mapping2 --out ./orange_black_mapping
vcftools --vcf ${VCF_all} --weir-fst-pop orange_allsamps2 --weir-fst-pop black_allsamps2 --out ./orange_black_allsamps

#vcftools --vcf ${VCF_wild} --weir-fst-pop lorquini_wild --weir-fst-pop hybrid_wild --out ./lorquini_hybrid_wild
#vcftools --vcf ${VCF_wild} --weir-fst-pop weidemeyerii_wild --weir-fst-pop hybrid_wild --out ./weidemeyerii_hybrid_wild

## the fst.qlog file also contains genome-wide Fst estimates (weighted and unweighted)
# wild caught:
# Weir and Cockerham mean Fst estimate: 0.10546
# Weir and Cockerham weighted Fst estimate: 0.16452

#------------------------------Fst sliding window
pwd
/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/genome_scan

# some info here: https://www.biostars.org/p/72498/#72660
# we can use vcftools and the following inputs to calculate Fst in sliding windows at distinct window steps
--fst-window-size <integer>
--fst-window-step <integer>

# sean uses 500-bp windows in 50-bp steps in his 2020 Mol. Biol. Evol. paper. will try that first
cat fst_slidingwindow3
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N fst_slidingwindow  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o fst_slidingwindow.qlog

module load vcftools

VCF_wild=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_wild_annotated_output.vcf
VCF_mapping=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_mappingbrood_annotated_output.vcf
VCF_all=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_allsamps_annotated_output.vcf


vcftools --vcf ${VCF_wild} --weir-fst-pop orange_wild2 --weir-fst-pop black_wild2 --fst-window-size 5000 --out ./orange_black_wild_sliding_nostep
vcftools --vcf ${VCF_mapping} --weir-fst-pop orange_mapping2 --weir-fst-pop black_mapping2 --fst-window-size 5000 --out ./orange_black_mapping_sliding_nostep
vcftools --vcf ${VCF_all} --weir-fst-pop orange_allsamps2 --weir-fst-pop black_allsamps2 --fst-window-size 5000 --out ./orange_black_allsamps_sliding_nostep


# optix is ~1500 bp long, so going to increase the window size 

cat fst_slidingwindow2
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N fst_slidingwindow2  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o fst_slidingwindow2.qlog

module load vcftools

VCF_wild=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_wild_annotated_output.vcf
VCF_mapping=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_mappingbrood_annotated_output.vcf
VCF_all=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_allsamps_annotated_output.vcf


vcftools --vcf ${VCF_wild} --weir-fst-pop orange_wild2 --weir-fst-pop black_wild2 --fst-window-size 1500 --fst-window-step 150 --out ./orange_black_wild_sliding
vcftools --vcf ${VCF_mapping} --weir-fst-pop orange_mapping2 --weir-fst-pop black_mapping2 --fst-window-size 1500 --fst-window-step 150 --out ./orange_black_mapping_sliding
vcftools --vcf ${VCF_all} --weir-fst-pop orange_allsamps2 --weir-fst-pop black_allsamps2 --fst-window-size 1500 --fst-window-step 150 --out ./orange_black_allsamps_sliding

# now try the combinations that Nick VanKuren suggested after looking at the above plots
cat fst_slidingwindow3
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N fst_slidingwindow3  # job name, anything you want
#$ -l h_rt=12:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o fst_slidingwindow3.qlog

module load vcftools

VCF_wild=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_wild_annotated_output.vcf
VCF_mapping=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_mappingbrood_annotated_output.vcf
VCF_all=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_allsamps_annotated_output.vcf

# focus on wild caught individuals and run multiple iterations to try and zoom in on peak around optix. first with no overlap
vcftools --vcf ${VCF_wild} --weir-fst-pop orange_wild2 --weir-fst-pop black_wild2 --fst-window-size 5000 --out ./orange_black_wild_5kb_nostep
vcftools --vcf ${VCF_wild} --weir-fst-pop orange_wild2 --weir-fst-pop black_wild2 --fst-window-size 10000 --out ./orange_black_wild_10kb_nostep
vcftools --vcf ${VCF_wild} --weir-fst-pop orange_wild2 --weir-fst-pop black_wild2 --fst-window-size 25000 --out ./orange_black_wild_25kb_nostep

# now with overlap
vcftools --vcf ${VCF_wild} --weir-fst-pop orange_wild2 --weir-fst-pop black_wild2 --fst-window-size 5000 --fst-window-step 1000 --out ./orange_black_wild_5kb_1kbstep
vcftools --vcf ${VCF_wild} --weir-fst-pop orange_wild2 --weir-fst-pop black_wild2 --fst-window-size 10000 --fst-window-step 2000 --out ./orange_black_wild_10kb_2kbstep
vcftools --vcf ${VCF_wild} --weir-fst-pop orange_wild2 --weir-fst-pop black_wild2 --fst-window-size 25000 --fst-window-step 5000 --out ./orange_black_wild_25kb_5kbstep



# in order to plot the Fst in IGV, need to convert to bedgraph file (https://genome.ucsc.edu/goldenPath/help/bedgraph.html)
# first, cut the columns we need from the vcftools output 
cut -f1,2,3,5 orange_black_wild_25kb_5kbstep.windowed.weir.fst > orange_black_wild_25kb_5kbstep.windowed.bed

# add a header to this file "track type=bedGraph" and then save it as a .bedgraph file. This will load into IGV!

# this is how you would add an arbitrary name to a new column that corresponds to line number, but turns out we didn't need this column like Nick suspected
awk 'BEGIN{OFS="\t"; n = 1} {print $1, $2, $3, "window" n, $4;n=n+1}' orange_black_wild_25kb_5kbstep.windowed.bed > orange_black_wild_25kb_5kbstep.windowed.named.bed


#------------------------------Explore Fixed SNPs

# can use tabix to extract SNPs at specific positions from a vcf file, but first have to bgzip and then index it:

pwd
/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe

# compress vcf's
cat zip_vcfs
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N zip_vcfs  # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o zip_vcfs.qlog

module load htslib

bgzip *.vcf

# index vcf's
cat index_vdfs
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N index_vcfs  # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o index_vcfs.qlog

tabix -p vcf AllLoci_wild_annotated_output.vcf.gz


# then to extract specific SNPs:
tabix AllLoci_wild_annotated_output.vcf.gz "ptg000014l:8402414-8402414"

# here's a look at our fixed SNPs:
#CHROM  	POS 		ID  		   REF ALT QUAL FILTER  INFO    						FORMAT      [other samples...]
ptg000014l	7700503	ptg000014l_7700503	T	A	29	PASS	NS=42;DP=280;AF=0.480855	GT:DP:GL:PL:GP	0/0:4:0,-1.20414,-14.2184:0,12,142:0,9.3638,142.85	0/0:11:0,-3.31137,-29.2688:0,33,293:0,30.4362,293.353	0/0:12:0,-3.61241,-30.5963:0,36,306:0,33.4465,306.629	0/0:2:0,-0.602068,-7.29076:0,6,73:0,3.34313,73.5731	0/0:7:0,-2.10724,-20.8623:0,21,209:0,18.3948,209.288	0/0:4:0,-1.20414,-13.7559:0,12,138:0,9.3638,138.224	0/0:7:0,-2.10724,-20.9915:0,21,210:0,18.3948,210.58	0/0:8:0,-2.40827,-23.3775:0,24,234:0,21.4052,234.44	0/0:4:0,-1.20414,-12.9652:0,12,130:0,9.3638,130.318	0/0:2:0,-0.602068,-7.59997:0,6,76:0,3.34313,76.6652	0/0:4:0,-1.20414,-14.6106:0,12,146:0,9.3638,146.771	0/0:1:0,-0.301034,-4:0,3,40:0,0.33279,40.6655	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	1/1:3:-11.2617,-0.903101,0:113,9,0:111.952,5.68796,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	1/1:2:-8,-0.602068,0:80,6,0:79.3345,2.67762,0	1/1:9:-25.7632,-2.7093,0:258,27,0:256.966,23.75,0	1/1:3:-11.1782,-0.903101,0:112,9,0:111.117,5.68796,0	1/1:4:-14.0514,-1.20414,0:141,12,0:139.848,8.6983,0	1/1:2:-7.02374,-0.602068,0:70,6,0:69.5719,2.67762,0	1/1:4:-13.5315,-1.20414,0:135,12,0:134.65,8.6983,0	1/1:2:-7.49995,-0.602068,0:75,6,0:74.334,2.67762,0	1/1:3:-9.98936,-0.903101,0:100,9,0:99.2281,5.68796,0	1/1:2:-7.29994,-0.602068,0:73,6,0:72.3339,2.67762,0	1/1:7:-20.3722,-2.10724,0:204,21,0:203.056,17.7293,0	0/0:12:0,-3.61241,-28.3513:0,36,284:0,33.4465,284.179	0/0:7:0,-2.10724,-20.6086:0,21,206:0,18.3948,206.752	0/0:19:0,-5.71964,-34.0673:0,57,341:0,54.5189,341.338	0/0:12:0,-3.61241,-29.1674:0,36,292:0,33.4465,292.339	0/0:8:0,-2.40827,-23.2133:0,24,232:0,21.4052,232.799	0/0:13:0,-3.91344,-29.0751:0,39,291:0,36.4568,291.417	0/0:20:0,-6.02068,-36.8652:0,60,369:0,57.5292,369.318	0/0:5:0,-1.50517,-15.8998:0,15,159:0,12.3741,159.663	1/1:2:-7.59995,-0.602068,0:76,6,0:75.334,2.67762,0	1/1:8:-22.875,-2.40827,0:229,24,0:228.084,20.7397,0	1/1:14:-32.8124,-4.21447,0:328,42,0:327.459,38.8017,0	1/1:4:-12.451,-1.20414,0:125,12,0:123.845,8.6983,0	1/1:7:-19.2862,-2.10724,0:193,21,0:192.196,17.7293,0	1/1:6:-15.6197,-1.8062,0:156,18,0:155.532,14.719,0	1/1:2:-7.22373,-0.602068,0:72,6,0:71.5718,2.67762,0	1/1:8:-22.5306,-2.40827,0:225,24,0:224.64,20.7397,0	1/1:5:-14.3134,-1.50517,0:143,15,0:142.469,11.7086,0	1/1:9:-26.2694,-2.7093,0:263,27,0:262.028,23.75,0	0/0:12:0,-3.61241,-28.3496:0,36,283:0,33.4465,284.162
ptg000014l	7700505	ptg000014l_7700505	T	C	29	PASS	NS=42;DP=282;AF=0.483208	GT:DP:GL:PL:GP	0/0:4:0,-1.20414,-13.2993:0,12,133:0,9.32287,133.577	0/0:11:0,-3.31137,-29.6509:0,33,297:0,30.3952,297.092	0/0:12:0,-3.61241,-29.9855:0,36,300:0,33.4056,300.439	0/0:2:0,-0.602068,-7.29076:0,6,73:0,3.30219,73.4913	0/0:7:0,-2.10724,-20.9509:0,21,210:0,18.3539,210.093	0/0:3:0,-0.903101,-10.6568:0,9,107:0,6.31253,107.152	0/0:7:0,-2.10724,-21.2953:0,21,213:0,18.3539,213.536	0/0:9:0,-2.7093,-24.7058:0,27,247:0,24.3746,247.642	0/0:3:0,-0.903101,-11.1782:0,9,112:0,6.31253,112.366	0/0:2:0,-0.602068,-7.59997:0,6,76:0,3.30219,76.5833	0/0:4:0,-1.20414,-14.6106:0,12,146:0,9.32287,146.689	0/0:1:0,-0.301034,-4:0,3,40:0,0.291856,40.5836	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	1/1:3:-11.2617,-0.903101,0:113,9,0:112.034,5.7289,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	1/1:2:-8,-0.602068,0:80,6,0:79.4164,2.71856,0	1/1:9:-26.5499,-2.7093,0:265,27,0:264.915,23.7909,0	1/1:3:-11.2617,-0.903101,0:113,9,0:112.034,5.7289,0	1/1:4:-14.469,-1.20414,0:145,12,0:144.106,8.73924,0	1/1:2:-7.29076,-0.602068,0:73,6,0:72.324,2.71856,0	1/1:3:-9.36081,-0.903101,0:94,9,0:93.0245,5.7289,0	1/1:2:-8,-0.602068,0:80,6,0:79.4164,2.71856,0	1/1:3:-9.98936,-0.903101,0:100,9,0:99.31,5.7289,0	1/1:2:-6.39987,-0.602068,0:64,6,0:63.4151,2.71856,0	1/1:8:-21.8001,-2.40827,0:218,24,0:217.417,20.7806,0	0/0:13:0,-3.91344,-29.9303:0,39,299:0,36.4159,299.887	0/0:7:0,-2.10724,-20.6922:0,21,207:0,18.3539,207.506	0/0:20:0,-6.02068,-34.4624:0,60,345:0,57.4883,345.208	0/0:12:0,-3.61241,-29.611:0,36,296:0,33.4056,296.694	0/0:8:0,-2.40827,-23.4594:0,24,235:0,21.3642,235.178	0/0:13:0,-3.91344,-28.9312:0,39,289:0,36.4159,289.895	0/0:20:0,-6.02068,-36.1836:0,60,362:0,57.4883,362.419	0/0:4:0,-1.20414,-13.2315:0,12,132:0,9.32287,132.898	1/1:2:-7.59995,-0.602068,0:76,6,0:75.4159,2.71856,0	1/1:9:-23.6081,-2.7093,0:236,27,0:235.498,23.7909,0	1/1:14:-32.5563,-4.21447,0:326,42,0:324.979,38.8426,0	1/1:4:-12.3812,-1.20414,0:124,12,0:123.228,8.73924,0	1/1:7:-19.0794,-2.10724,0:191,21,0:190.21,17.7703,0	1/1:6:-15.6849,-1.8062,0:157,18,0:156.265,14.7599,0	1/1:2:-7.22373,-0.602068,0:72,6,0:71.6537,2.71856,0	1/1:8:-22.3406,-2.40827,0:223,24,0:222.823,20.7806,0	1/1:6:-18.0446,-1.8062,0:180,18,0:179.862,14.7599,0	1/1:9:-25.9806,-2.7093,0:260,27,0:259.222,23.7909,0	0/0:12:0,-3.61241,-28.9069:0,36,289:0,33.4056,289.653
ptg000014l	8389176	ptg000014l_8389176	T	A	29	PASS	NS=41;DP=260;AF=0.489608	GT:DP:GL:PL:GP	0/0:8:0,-2.40827,-23.3596:0,24,234:0,21.253,233.957	0/0:5:0,-1.50517,-17.4081:0,15,174:0,12.2219,174.442	0/0:2:0,-0.602068,-7.89999:0,6,79:0,3.19093,79.3611	0/0:7:0,-2.10724,-21.9233:0,21,219:0,18.2426,219.594	0/0:6:0,-1.8062,-18.6536:0,18,187:0,15.2323,186.897	0/0:5:0,-1.50517,-16.33:0,15,163:0,12.2219,163.661	0/0:1:0,-0.301034,-3.7:0,3,37:0,0.180593,37.3611	0/0:5:0,-1.50517,-16.1807:0,15,162:0,12.2219,162.168	0/0:1:0,-0.301034,-4:0,3,40:0,0.180593,40.3611	./.:0:0,0,0:0,0,0:0,0,0	0/0:2:0,-0.602068,-7.59998:0,6,76:0,3.19093,76.361	0/0:4:0,-1.20414,-14.527:0,12,145:0,9.21161,145.632	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:1:0,-0.301034,-4:0,3,40:0,0.180593,40.3611	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,./.:0:0,0,0:0,0,0:0,0,0	1/1:3:-11.0947,-0.903101,0:111,9,0:110.586,5.84016,0	1/1:8:-23.8466,-2.40827,0:238,24,0:238.105,20.8919,0	1/1:2:-7.20725,-0.602068,0:72,6,0:71.7114,2.82982,0	1/1:4:-14.3855,-1.20414,0:144,12,0:143.493,8.8505,0	1/1:2:-6.22262,-0.602068,0:62,6,0:61.8651,2.82982,0	1/1:5:-16.0323,-1.50517,0:160,15,0:159.961,11.8608,0	1/1:3:-11.3404,-0.903101,0:113,9,0:113.043,5.84016,0	1/1:4:-14.6106,-1.20414,0:146,12,0:145.745,8.8505,0	1/1:3:-11.1782,-0.903101,0:112,9,0:111.421,5.84016,0	1/1:3:-10.8441,-0.903101,0:108,9,0:108.08,5.84016,0	0/0:14:0,-4.21447,-31.8804:0,42,319:0,39.315,319.165	0/0:5:0,-1.50517,-13.0824:0,15,131:0,12.2219,131.185	0/0:10:0,-3.01034,-25.5967:0,30,256:0,27.2736,256.328	0/0:10:0,-3.01034,-25.7472:0,30,257:0,27.2736,257.833	0/0:11:0,-3.31137,-27.1459:0,33,271:0,30.284,271.82	0/0:22:0,-6.62274,-36.8812:0,66,369:0,63.3977,369.173	0/0:5:0,-1.50517,-16.5735:0,15,166:0,12.2219,166.096	0/0:20:0,-6.02068,-37.0779:0,60,371:0,57.377,371.14	1/1:10:-26.3399,-3.01034,0:263,30,0:263.038,26.9125,0	1/1:5:-16.2659,-1.50517,0:163,15,0:162.297,11.8608,0	1/1:11:-27.6147,-3.31137,0:276,33,0:275.785,29.9229,0	1/1:4:-13.248,-1.20414,0:132,12,0:132.119,8.8505,0	1/1:9:-24.447,-2.7093,0:244,27,0:244.109,23.9022,0	1/1:2:-7.47427,-0.602068,0:75,6,0:74.3816,2.82982,0	1/1:4:-13.2667,-1.20414,0:133,12,0:132.305,8.8505,0	1/1:10:-24.4312,-3.01034,0:244,30,0:243.951,26.9125,0	1/1:12:-27.9349,-3.61241,0:279,36,0:278.988,32.9332,0	1/1:7:-20.0683,-2.10724,0:201,21,0:200.322,17.8815,0	0/0:5:0,-1.50517,-12.1485:0,15,121:0,12.2219,121.846
ptg000014l	8390379	ptg000014l_8390379	T	A	29	PASS	NS=41;DP=326;AF=0.480011	GT:DP:GL:PL:GP	0/0:7:0,-2.10724,-21.6669:0,21,217:0,18.4095,217.364	0/0:5:0,-1.50517,-17.1639:0,15,172:0,12.3888,172.334	0/0:4:0,-1.20414,-13.9226:0,12,139:0,9.37848,139.921	0/0:12:0,-3.61241,-30.7666:0,36,308:0,33.4612,308.36	0/0:2:0,-0.602068,-8:0,6,80:0,3.3578,80.6948	0/0:6:0,-1.8062,-18.9818:0,18,190:0,15.3992,190.512	0/0:1:0,-0.301034,-4:0,3,40:0,0.34746,40.6948	0/0:7:0,-2.10724,-20.9502:0,21,210:0,18.4095,210.196	0/0:4:0,-1.20414,-14.6106:0,12,146:0,9.37848,146.8	0/0:4:0,-1.20414,-13.6693:0,12,137:0,9.37848,137.388	0/0:1:0,-0.301034,-4:0,3,40:0,0.34746,40.6948	0/0:2:0,-0.602068,-8:0,6,80:0,3.3578,80.6948	0/0:2:0,-0.602068,-8:0,6,80:0,3.3578,80.6948	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	1/1:6:-19.1688,-1.8062,0:192,18,0:190.993,14.7043,0	1/1:9:-25.8326,-2.7093,0:258,27,0:257.631,23.7353,0	1/1:4:-14.3855,-1.20414,0:144,12,0:143.16,8.68363,0	1/1:3:-11.0782,-0.903101,0:111,9,0:110.087,5.67329,0	1/1:4:-14.6106,-1.20414,0:146,12,0:145.411,8.68363,0	1/1:2:-7.29076,-0.602068,0:73,6,0:72.2128,2.66296,0	1/1:3:-11.2404,-0.903101,0:112,9,0:111.709,5.67329,0	1/1:4:-13.7008,-1.20414,0:137,12,0:136.313,8.68363,0	1/1:4:-14.6106,-1.20414,0:146,12,0:145.411,8.68363,0	0/0:18:0,-5.41861,-33.6199:0,54,336:0,51.5232,336.894	0/0:6:0,-1.8062,-17.7315:0,18,177:0,15.3992,178.01	0/0:14:0,-4.21447,-30.6246:0,42,306:0,39.4819,306.941	0/0:12:0,-3.61241,-28.264:0,36,283:0,33.4612,283.334	0/0:7:0,-2.10724,-20.7211:0,21,207:0,18.4095,207.906	0/0:28:0,-8.42895,-36.3435:0,84,363:0,81.6266,364.13	0/0:7:0,-2.10724,-21.5985:0,21,216:0,18.4095,216.679	0/0:2:0,-0.602068,-7.69996:0,6,77:0,3.3578,77.6945	1/1:8:-22.8036,-2.40827,0:228,24,0:227.341,20.725,0	1/1:9:-23.2578,-2.7093,0:233,27,0:231.884,23.7353,0	1/1:12:-28.7255,-3.61241,0:287,36,0:286.56,32.7663,0	1/1:24:-38.6291,-7.22481,0:386,72,0:385.596,68.8904,0	1/1:12:-28.621,-3.61241,0:286,36,0:285.515,32.7663,0	1/1:11:-27.7297,-3.31137,0:277,33,0:276.602,29.756,0	1/1:6:-19.2655,-1.8062,0:193,18,0:191.96,14.7043,0	1/1:11:-26.6202,-3.31137,0:266,33,0:265.508,29.756,0	1/1:15:-32.9844,-4.51551,0:330,45,0:329.15,41.7974,0	1/1:11:-27.521,-3.31137,0:275,33,0:274.516,29.756,0	0/0:17:0,-5.11758,-34.4813:0,51,345:0,48.5129,345.507
ptg000014l	8402414	ptg000014l_8402414	A	C	29	PASS	NS=42;DP=469;AF=0.49842	GT:DP:GL:PL:GP	1/1:11:-29.338,-3.31137,0:293,33,0:293.325,30.076,0	1/1:5:-16.4121,-1.50517,0:164,15,0:164.066,12.0139,0	1/1:9:-25.8702,-2.7093,0:259,27,0:258.647,24.0553,0	1/1:15:-34.0297,-4.51551,0:340,45,0:340.242,42.1173,0	1/1:4:-14.3019,-1.20414,0:143,12,0:142.964,9.00361,0	1/1:13:-29.9746,-3.91344,0:300,39,0:299.691,36.0966,0	1/1:3:-11.3404,-0.903101,0:113,9,0:113.349,5.99327,0	1/1:7:-21.9656,-2.10724,0:220,21,0:219.601,18.0346,0	1/1:4:-14.4435,-1.20414,0:144,12,0:144.38,9.00361,0	1/1:2:-7.09993,-0.602068,0:71,6,0:70.9444,2.98293,0	1/1:3:-11.2404,-0.903101,0:112,9,0:112.349,5.99327,0	1/1:2:-8,-0.602068,0:80,6,0:79.9451,2.98293,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:2:0,-0.602068,-8:0,6,80:0,3.03782,80.0549	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:5:0,-1.50517,-17.2628:0,15,173:0,12.0688,172.683	0/0:11:0,-3.31137,-29.3661:0,33,294:0,30.1309,293.716	0/0:17:0,-5.11758,-35.7764:0,51,358:0,48.1929,357.819	0/0:7:0,-2.10724,-21.53:0,21,215:0,18.0895,215.355	0/0:8:0,-2.40827,-24.2719:0,24,243:0,21.0999,242.774	0/0:2:0,-0.602068,-7.89999:0,6,79:0,3.03782,79.0548	0/0:4:0,-1.20414,-13.9259:0,12,139:0,9.0585,139.314	0/0:5:0,-1.50517,-14.8119:0,15,148:0,12.0688,148.174	0/0:4:0,-1.20414,-14.6106:0,12,146:0,9.0585,146.161	0/0:14:0,-4.21447,-32.9783:0,42,330:0,39.1619,329.838	1/1:26:-40.2181,-7.82688,0:402,78,0:402.126,75.231,0	1/1:24:-40.0582,-7.22481,0:401,72,0:400.527,69.2104,0	1/1:26:-39.4459,-7.82688,0:394,78,0:394.404,75.231,0	1/1:16:-33.6725,-4.81654,0:337,48,0:336.67,45.1277,0	1/1:20:-37.2319,-6.02068,0:372,60,0:372.264,57.169,0	1/1:41:-42.8932,-12.3424,0:429,123,0:428.878,120.386,0	1/1:19:-36.1552,-5.71964,0:362,57,0:361.497,54.1587,0	1/1:24:-39.3601,-7.22481,0:394,72,0:393.546,69.2104,0	0/0:5:0,-1.50517,-16.1412:0,15,161:0,12.0688,161.467	0/0:6:0,-1.8062,-17.9187:0,18,179:0,15.0792,179.242	0/0:17:0,-5.11758,-34.1598:0,51,342:0,48.1929,341.653	0/0:8:0,-2.40827,-22.438:0,24,224:0,21.0999,224.435	0/0:21:0,-6.32171,-37.0854:0,63,371:0,60.2342,370.909	0/0:5:0,-1.50517,-15.4371:0,15,154:0,12.0688,154.426	0/0:5:0,-1.50517,-16.2095:0,15,162:0,12.0688,162.15	0/0:10:0,-3.01034,-25.4832:0,30,255:0,27.1205,254.887	0/0:11:0,-3.31137,-25.923:0,33,259:0,30.1309,259.285	0/0:13:0,-3.91344,-29.4436:0,39,294:0,36.1515,294.491	1/1:15:-31.089,-4.51551,0:311,45,0:310.835,42.1173,0
ptg000014l	8404283	ptg000014l_8404283	T	G	29	PASS	NS=42;DP=400;AF=0.490394	GT:DP:GL:PL:GP	0/0:9:0,-2.7093,-24.9506:0,27,250:0,24.2496,249.84	0/0:6:0,-1.8062,-17.9397:0,18,179:0,15.2186,179.731	0/0:12:0,-3.61241,-30.5936:0,36,306:0,33.2807,306.27	0/0:5:0,-1.50517,-15.7399:0,15,157:0,12.2083,157.733	0/0:7:0,-2.10724,-21.6953:0,21,217:0,18.229,217.286	0/0:7:0,-2.10724,-21.6305:0,21,216:0,18.229,216.639	0/0:5:0,-1.50517,-17.193:0,15,172:0,12.2083,172.264	0/0:11:0,-3.31137,-29.3834:0,33,294:0,30.2703,294.168	0/0:2:0,-0.602068,-8:0,6,80:0,3.17728,80.3338	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:2:0,-0.602068,-8:0,6,80:0,3.17728,80.3338	0/0:2:0,-0.602068,-8:0,6,80:0,3.17728,80.3338	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:2:0,-0.602068,-8:0,6,80:0,3.17728,80.3338	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:2:0,-0.602068,-7.39994:0,6,74:0,3.17728,74.3332	./.:0:0,0,0:0,0,0:0,0,0	1/1:4:-13.4663,-1.20414,0:135,12,0:134.33,8.86415,0	1/1:10:-27.1819,-3.01034,0:272,30,0:271.485,26.9262,0	1/1:7:-21.9156,-2.10724,0:219,21,0:218.822,17.8952,0	1/1:7:-20.8349,-2.10724,0:208,21,0:208.015,17.8952,0	1/1:10:-26.9417,-3.01034,0:269,30,0:269.083,26.9262,0	1/1:11:-28.2956,-3.31137,0:283,33,0:282.622,29.9365,0	1/1:9:-26.433,-2.7093,0:264,27,0:263.996,23.9158,0	1/1:7:-20.7146,-2.10724,0:207,21,0:206.812,17.8952,0	1/1:14:-33.3374,-4.21447,0:333,42,0:333.04,38.9675,0	1/1:2:-7.89999,-0.602068,0:79,6,0:78.6661,2.84348,0	0/0:8:0,-2.40827,-22.1678:0,24,222:0,21.2393,222.012	0/0:13:0,-3.91344,-29.5117:0,39,295:0,36.291,295.451	0/0:22:0,-6.62274,-35.738:0,66,357:0,63.384,357.713	0/0:11:0,-3.31137,-27.1542:0,33,272:0,30.2703,271.875	0/0:13:0,-3.91344,-30.4794:0,39,305:0,36.291,305.128	0/0:30:0,-9.03101,-39.7402:0,90,397:0,87.4667,397.735	0/0:19:0,-5.71964,-34.9735:0,57,350:0,54.353,350.069	0/0:8:0,-2.40827,-22.1234:0,24,221:0,21.2393,221.567	1/1:8:-22.8406,-2.40827,0:228,24,0:228.072,20.9055,0	1/1:6:-17.5364,-1.8062,0:175,18,0:175.03,14.8848,0	1/1:17:-34.0574,-5.11758,0:341,51,0:340.24,47.9986,0	1/1:19:-32.3064,-1.61932,0:323,16,0:322.73,13.016,0	1/1:10:-26.2387,-3.01034,0:262,30,0:262.053,26.9262,0	1/1:10:-25.9825,-3.01034,0:260,30,0:259.491,26.9262,0	1/1:6:-20.0445,-1.8062,0:200,18,0:200.111,14.8848,0	1/1:19:-37.0301,-5.71964,0:370,57,0:369.967,54.0192,0	1/1:10:-25.8355,-3.01034,0:258,30,0:258.021,26.9262,0	1/1:11:-27.9071,-3.31137,0:279,33,0:278.738,29.9365,0	0/0:7:0,-2.10724,-20.9897:0,21,210:0,18.229,210.23
ptg000014l	8458650	ptg000014l_8458650	C	T	29	PASS	NS=42;DP=644;AF=0.483942	GT:DP:GL:PL:GP	0/0:15:0,-4.51551,-33.3021:0,45,333:0,42.4238,333.579	0/0:20:0,-6.02068,-38.4236:0,60,384:0,57.4755,384.794	0/0:17:0,-5.11758,-35.7098:0,51,357:0,48.4445,357.656	0/0:23:0,-6.92378,-38.9693:0,69,390:0,66.5065,390.251	0/0:19:0,-5.71964,-37.0733:0,57,371:0,54.4652,371.291	0/0:14:0,-4.21447,-32.2165:0,42,322:0,39.4135,322.723	0/0:4:0,-1.20414,-12.7146:0,12,127:0,9.31011,127.704	0/0:13:0,-3.91344,-30.9972:0,39,310:0,36.4032,310.53	0/0:5:0,-1.50517,-17.1095:0,15,171:0,12.3204,171.653	./.:0:0,0,0:0,0,0:0,0,0	0/0:6:0,-1.8062,-19.7464:0,18,197:0,15.3308,198.022	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:4:0,-1.20414,-13.2667:0,12,133:0,9.31011,133.225	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:1:0,-0.301034,-3.9:0,3,39:0,0.279096,39.5581	./.:0:0,0,0:0,0,0:0,0,0	0/0:2:0,-0.602068,-8:0,6,80:0,3.28943,80.5581	1/1:14:-33.0656,-4.21447,0:331,42,0:330.098,38.8554,0	1/1:14:-32.8133,-4.21447,0:328,42,0:327.575,38.8554,0	1/1:12:-30.274,-3.61241,0:303,36,0:302.182,32.8347,0	1/1:18:-36.5398,-5.41861,0:365,54,0:364.839,50.8967,0	1/1:16:-34.6301,-4.81654,0:346,48,0:345.743,44.8761,0	1/1:7:-20.8689,-2.10724,0:209,21,0:208.131,17.783,0	1/1:2:-7.49995,-0.602068,0:75,6,0:74.4414,2.73132,0	1/1:8:-23.4803,-2.40827,0:235,24,0:234.245,20.7934,0	1/1:5:-17.2628,-1.50517,0:173,15,0:172.07,11.7623,0	1/1:10:-27.3534,-3.01034,0:274,30,0:272.976,26.814,0	0/0:23:0,-6.92378,-36.8763:0,69,369:0,66.5065,369.321	0/0:11:0,-3.31137,-26.8698:0,33,269:0,30.3825,269.256	0/0:29:0,-8.72998,-40.0403:0,87,400:0,84.5686,400.961	0/0:37:0,-11.1383,-42.0461:0,111,420:0,108.651,421.019	0/0:36:0,-10.8372,-42.1559:0,108,422:0,105.641,422.117	0/0:35:0,-10.5362,-40.9559:0,105,410:0,102.631,410.117	0/0:27:0,-8.12791,-39.3484:0,81,393:0,78.5479,394.042	0/0:18:0,-5.41861,-33.1239:0,54,331:0,51.4548,331.797	1/1:8:-23.0449,-2.40827,0:230,24,0:229.891,20.7934,0	1/1:14:-31.9052,-4.21447,0:319,42,0:318.494,38.8554,0	1/1:21:-37.6952,-6.32171,0:377,63,0:376.393,59.9277,0	1/1:18:-35.7618,-5.41861,0:358,54,0:357.06,50.8967,0	1/1:21:-36.8877,-6.32171,0:369,63,0:368.319,59.9277,0	1/1:11:-28.0708,-3.31137,0:281,33,0:280.15,29.8244,0	1/1:14:-31.5592,-4.21447,0:316,42,0:315.034,38.8554,0	1/1:19:-36.4977,-5.71964,0:365,57,0:364.419,53.9071,0	1/1:20:-38.2467,-6.02068,0:382,60,0:381.908,56.9174,0	1/1:20:-37.7944,-6.02068,0:378,60,0:377.386,56.9174,0	0/0:13:0,-3.91344,-29.3272:0,39,293:0,36.4032,293.83
ptg000014l	8460881	ptg000014l_8460881	A	C	29	PASS	NS=44;DP=425;AF=0.473105	GT:DP:GL:PL:GP	0/0:10:0,-3.01034,-27.778:0,30,278:0,27.5608,278.715	0/0:15:0,-4.51551,-33.6986:0,45,337:0,42.6124,337.922	0/0:5:0,-1.50517,-15.3803:0,15,154:0,12.5091,154.738	0/0:7:0,-2.10724,-21.2813:0,21,213:0,18.5297,213.748	0/0:15:0,-4.51551,-33.7145:0,45,337:0,42.6124,338.08	0/0:4:0,-1.20414,-14.0092:0,12,140:0,9.49872,141.028	0/0:5:0,-1.50517,-17.1639:0,15,172:0,12.5091,172.575	0/0:5:0,-1.50517,-17.0308:0,15,170:0,12.5091,171.243	./.:0:0,0,0:0,0,0:0,0,0	0/0:2:0,-0.602068,-8:0,6,80:0,3.47805,80.9353	./.:0:0,0,0:0,0,0:0,0,0	0/0:12:0,-3.61241,-30.4117:0,36,304:0,33.5814,305.052	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:8:0,-2.40827,-23.8378:0,24,238:0,21.5401,239.313	0/0:1:0,-0.301034,-4:0,3,40:0,0.467709,40.9353	./.:0:0,0,0:0,0,0:0,0,0	0/0:1:0,-0.301034,-4:0,3,40:0,0.467709,40.9353	0/0:2:0,-0.602068,-8:0,6,80:0,3.47805,80.9353	0/0:1:0,-0.301034,-3.9:0,3,39:0,0.467709,39.9353	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	1/1:9:-24.3217,-2.7093,0:243,27,0:242.282,23.6151,0	1/1:11:-28.3523,-3.31137,0:284,33,0:282.588,29.6357,0	1/1:3:-11.3404,-0.903101,0:113,9,0:112.469,5.55304,0	1/1:12:-30.165,-3.61241,0:302,36,0:300.714,32.6461,0	1/1:6:-19.2063,-1.8062,0:192,18,0:191.127,14.5841,0	1/1:9:-25.4837,-2.7093,0:255,27,0:253.902,23.6151,0	1/1:7:-20.831,-2.10724,0:208,21,0:207.374,17.5944,0	1/1:3:-9.91953,-0.903101,0:99,9,0:98.26,5.55304,0	1/1:6:-19.1178,-1.8062,0:191,18,0:190.242,14.5841,0	1/1:5:-14.0017,-1.50517,0:140,15,0:139.081,11.5737,0	0/0:18:0,-5.41861,-33.4495:0,54,334:0,51.6435,335.431	0/0:17:0,-5.11758,-34.7829:0,51,348:0,48.6331,348.765	0/0:16:0,-4.81654,-32.3781:0,48,324:0,45.6228,324.716	0/0:16:0,-4.81654,-33.4259:0,48,334:0,45.6228,335.195	0/0:17:0,-5.11758,-33.2368:0,51,332:0,48.6331,333.303	0/0:25:0,-7.52585,-37.7673:0,75,378:0,72.7158,378.608	0/0:16:0,-4.81654,-33.1853:0,48,332:0,45.6228,332.788	0/0:11:0,-3.31137,-27.0355:0,33,270:0,30.5711,271.29	1/1:10:-21.2404,-3.01034,0:212,30,0:211.469,26.6254,0	1/1:8:-21.615,-2.40827,0:216,24,0:215.214,20.6047,0	1/1:16:-32.9953,-4.81654,0:330,48,0:329.018,44.6874,0	1/1:8:-22.9095,-2.40827,0:229,24,0:228.16,20.6047,0	1/1:10:-25.994,-3.01034,0:260,30,0:259.005,26.6254,0	1/1:4:-13.5926,-1.20414,0:136,12,0:134.991,8.56338,0	1/1:9:-24.7038,-2.7093,0:247,27,0:246.103,23.6151,0	1/1:15:-31.8261,-4.51551,0:318,45,0:317.326,41.6771,0	1/1:21:-35.2995,-6.32171,0:353,63,0:352.06,59.7391,0	1/1:19:-35.1889,-5.71964,0:352,57,0:350.954,53.7185,0	0/0:5:0,-1.50517,-15.0897:0,15,151:0,12.5091,151.832
ptg000014l	8483079	ptg000014l_8483079	C	A	29	PASS	NS=41;DP=398;AF=0.463201	GT:DP:GL:PL:GP	0/0:3:0,-0.903101,-11.0782:0,9,111:0,6.66114,112.063	0/0:5:0,-1.50517,-16.3121:0,15,163:0,12.6818,164.402	0/0:14:0,-4.21447,-32.5147:0,42,325:0,39.7749,326.428	0/0:8:0,-2.40827,-23.6887:0,24,237:0,21.7128,238.167	0/0:8:0,-2.40827,-23.8466:0,24,238:0,21.7128,239.747	0/0:8:0,-2.40827,-23.6887:0,24,237:0,21.7128,238.167	0/0:4:0,-1.20414,-14.0514:0,12,141:0,9.67147,141.795	0/0:9:0,-2.7093,-25.6552:0,27,257:0,24.7232,257.833	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:14:0,-4.21447,-32.2451:0,42,322:0,39.7749,323.732	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:8:0,-2.40827,-24.2276:0,24,242:0,21.7128,243.557	./.:0:0,0,0:0,0,0:0,0,0	0/0:4:0,-1.20414,-14.1184:0,12,141:0,9.67147,142.465	0/0:2:0,-0.602068,-7.09985:0,6,71:0,3.6508,72.2793	0/0:4:0,-1.20414,-13.4336:0,12,134:0,9.67147,135.617	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	1/1:9:-24.0111,-2.7093,0:240,27,0:238.831,23.4423,0	1/1:4:-14.3855,-1.20414,0:144,12,0:142.574,8.39063,0	1/1:5:-16.8863,-1.50517,0:169,15,0:167.582,11.401,0	1/1:9:-25.6306,-2.7093,0:256,27,0:255.026,23.4423,0	1/1:8:-23.8674,-2.40827,0:239,24,0:237.394,20.432,0	1/1:5:-17.074,-1.50517,0:171,15,0:169.46,11.401,0	1/1:9:-25.9197,-2.7093,0:259,27,0:257.916,23.4423,0	1/1:5:-17.2246,-1.50517,0:172,15,0:170.965,11.401,0	1/1:6:-19.771,-1.8062,0:198,18,0:196.429,14.4113,0	0/0:22:0,-6.62274,-36.9161:0,66,369:0,63.8576,370.441	0/0:17:0,-5.11758,-33.5111:0,51,335:0,48.8059,336.392	0/0:13:0,-3.91344,-29.161:0,39,292:0,36.7645,292.89	0/0:15:0,-4.51551,-30.887:0,45,309:0,42.7852,310.151	0/0:14:0,-4.21447,-29.9905:0,42,300:0,39.7749,301.185	0/0:30:0,-9.03101,-38.7248:0,90,387:0,87.9403,388.529	0/0:19:0,-5.71964,-35.4846:0,57,355:0,54.8265,356.127	0/0:9:0,-2.7093,-25.5214:0,27,255:0,24.7232,256.495	1/1:14:-31.0089,-4.21447,0:310,42,0:308.808,38.494,0	1/1:3:-9.7264,-0.903101,0:97,9,0:95.9832,5.38029,0	1/1:13:-30.9198,-3.91344,0:309,39,0:307.917,35.4837,0	1/1:7:-20.4098,-2.10724,0:204,21,0:202.817,17.4216,0	1/1:7:-20.7665,-2.10724,0:208,21,0:206.384,17.4216,0	1/1:11:-28.1328,-3.31137,0:281,33,0:280.047,29.463,0	1/1:13:-30.5301,-3.91344,0:305,39,0:304.02,35.4837,0	1/1:13:-27.3886,-3.91344,0:274,39,0:272.605,35.4837,0	1/1:12:-27.4464,-3.61241,0:274,36,0:273.183,32.4733,0	1/1:3:-9.2773,-0.903101,0:93,9,0:91.4922,5.38029,0	0/0:12:0,-3.61241,-30.503:0,36,305:0,33.7542,306.311
ptg000014l	8512954	ptg000014l_8512954	T	G	29	PASS	NS=44;DP=490;AF=0.494307	GT:DP:GL:PL:GP	0/0:18:0,-5.41861,-36.1502:0,54,362:0,51.2747,361.7	0/0:15:0,-4.51551,-33.9542:0,45,340:0,42.2437,339.74	0/0:17:0,-5.11758,-35.8103:0,51,358:0,48.2643,358.301	0/0:9:0,-2.7093,-24.9195:0,27,249:0,24.1816,249.392	0/0:14:0,-4.21447,-33.4997:0,42,335:0,39.2333,335.195	0/0:15:0,-4.51551,-33.573:0,45,336:0,42.2437,335.928	0/0:6:0,-1.8062,-18.7554:0,18,188:0,15.1506,187.752	0/0:8:0,-2.40827,-23.5035:0,24,235:0,21.1713,235.233	0/0:3:0,-0.903101,-11.2617:0,9,113:0,6.11961,112.815	./.:0:0,0,0:0,0,0:0,0,0	0/0:7:0,-2.10724,-21.7598:0,21,218:0,18.161,217.795	0/0:1:0,-0.301034,-4:0,3,40:0,0.0989358,40.1978	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/1:5:-1.89519,0,-13.0425:19,0,130:21.8633,0,133.534	./.:0:0,0,0:0,0,0:0,0,0	0/0:1:0,-0.301034,-4:0,3,40:0,0.0989358,40.1978	0/1:1:-3.5,-0.301034,0:35,3,0:34.9011,0,0.0988591	./.:0:0,0,0:0,0,0:0,0,0	0/0:2:0,-0.602068,-8:0,6,80:0,3.10927,80.1978	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	1/1:8:-24.3579,-2.40827,0:244,24,0:243.382,20.9735,0	1/1:13:-29.9168,-3.91344,0:299,39,0:298.97,36.0252,0	1/1:7:-20.697,-2.10724,0:207,21,0:206.773,17.9632,0	1/1:5:-16.4234,-1.50517,0:164,15,0:164.036,11.9425,0	1/1:10:-27.2004,-3.01034,0:272,30,0:271.806,26.9942,0	1/1:6:-19.5353,-1.8062,0:195,18,0:195.155,14.9528,0	1/1:5:-16.76,-1.50517,0:168,15,0:167.402,11.9425,0	1/1:11:-29.2552,-3.31137,0:293,33,0:292.354,30.0045,0	1/1:2:-7.20725,-0.602068,0:72,6,0:71.8747,2.91148,0	1/1:10:-27.9998,-3.01034,0:280,30,0:279.8,26.9942,0	0/0:30:0,-9.03101,-39.5819:0,90,396:0,87.3987,396.016	0/0:27:0,-8.12791,-38.5599:0,81,386:0,78.3677,385.797	0/0:18:0,-5.41861,-34.4692:0,54,345:0,51.2747,344.89	0/0:17:0,-5.11758,-33.1743:0,51,332:0,48.2643,331.94	0/0:12:0,-3.61241,-28.8654:0,36,289:0,33.2127,288.852	0/0:19:0,-5.71964,-33.3394:0,57,333:0,54.285,333.592	0/0:17:0,-5.11758,-33.671:0,51,337:0,48.2643,336.907	0/0:9:0,-2.7093,-23.5135:0,27,235:0,24.1816,235.332	1/1:11:-27.4616,-3.31137,0:275,33,0:274.418,30.0045,0	1/1:14:-31.169,-4.21447,0:312,42,0:311.492,39.0355,0	1/1:17:-33.1675,-5.11758,0:332,51,0:331.477,48.0666,0	1/1:8:-21.0904,-2.40827,0:211,24,0:210.707,20.9735,0	1/1:17:-33.8035,-5.11758,0:338,51,0:337.837,48.0666,0	1/1:14:-32.1075,-4.21447,0:321,42,0:320.877,39.0355,0	1/1:6:-17.542,-1.8062,0:175,18,0:175.223,14.9528,0	1/1:15:-32.3644,-4.51551,0:324,45,0:323.446,42.0459,0	1/1:17:-32.3962,-5.11758,0:324,51,0:323.764,48.0666,0	1/1:9:-23.7782,-2.7093,0:238,27,0:237.584,23.9838,0	0/0:14:0,-4.21447,-29.6169:0,42,296:0,39.2333,296.366
ptg000014l	9038849	ptg000014l_9038849	A	G	29	PASS	NS=42;DP=379;AF=0.473757	GT:DP:GL:PL:GP	0/0:4:0,-1.20414,-14.6106:0,12,146:0,9.48736,147.018	0/0:5:0,-1.50517,-16.4192:0,15,164:0,12.4977,165.105	0/0:3:0,-0.903101,-11.2404:0,9,112:0,6.47702,113.317	0/0:8:0,-2.40827,-23.626:0,24,236:0,21.5287,237.173	0/0:10:0,-3.01034,-27.2299:0,30,272:0,27.5494,273.212	0/0:2:0,-0.602068,-8:0,6,80:0,3.46668,80.9126	0/0:4:0,-1.20414,-14.093:0,12,141:0,9.48736,141.842	0/0:3:0,-0.903101,-10.6441:0,9,106:0,6.47702,107.354	0/0:3:0,-0.903101,-11.3404:0,9,113:0,6.47702,114.317	./.:0:0,0,0:0,0,0:0,0,0	0/0:2:0,-0.602068,-7.49995:0,6,75:0,3.46668,75.9121	0/0:6:0,-1.8062,-18.4176:0,18,184:0,15.508,185.088	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	1/1:2:-8,-0.602068,0:80,6,0:79.0874,2.55407,0	0/0:6:0,-1.8062,-19.2466:0,18,192:0,15.508,193.379	0/1:5:-9.28488,0,-5.69506:93,0,57:95.4028,0,60.4172	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	0/0:7:0,-2.10724,-22.0626:0,21,221:0,18.5184,221.539	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,0	./.:0:0,0,0:0,0,0:0,0,1/1:13:-26.6085,-3.61241,0:266,36,0:265.172,32.6575,0	1/1:3:-9.98936,-0.903101,0:100,9,0:98.981,5.56441,0	1/1:8:-23.8984,-2.40827,0:239,24,0:238.071,20.6161,0	1/1:6:-19.414,-1.8062,0:194,18,0:193.227,14.5954,0	1/1:12:-30.6835,-3.61241,0:307,36,0:305.922,32.6575,0	./.:0:0,0,0:0,0,0:0,0,0	1/1:9:-25.9675,-2.7093,0:260,27,0:258.763,23.6264,0	1/1:11:-26.7571,-3.31137,0:268,33,0:266.659,29.6471,0	1/1:9:-25.8694,-2.7093,0:259,27,0:257.782,23.6264,0	0/0:17:0,-5.11758,-32.5719:0,51,326:0,48.6218,326.632	0/0:10:0,-3.01034,-27.0884:0,30,271:0,27.5494,271.797	0/0:21:0,-6.32171,-37.344:0,63,373:0,60.6631,374.353	0/0:9:0,-2.7093,-25.1536:0,27,252:0,24.539,252.449	0/0:15:0,-4.51551,-33.4367:0,45,334:0,42.6011,335.279	0/0:21:0,-6.32171,-36.0351:0,63,360:0,60.6631,361.264	0/0:19:0,-5.71964,-35.7995:0,57,358:0,54.6424,358.907	0/0:5:0,-1.50517,-16.2795:0,15,163:0,12.4977,163.708	1/1:10:-25.1265,-3.01034,0:251,30,0:250.352,26.6368,0	1/1:9:-22.1949,-2.7093,0:222,27,0:221.036,23.6264,0	1/1:4:-12.9831,-1.20414,0:130,12,0:128.918,8.57475,0	1/1:9:-25.0741,-2.7093,0:251,27,0:249.828,23.6264,0	1/1:17:-34.1502,-5.11758,0:342,51,0:340.589,47.7091,0	1/1:9:-23.3927,-2.7093,0:234,27,0:233.014,23.6264,0	1/1:12:-27.2779,-3.61241,0:273,36,0:271.866,32.6575,0	1/1:11:-27.4872,-3.31137,0:275,33,0:273.96,29.6471,0	1/1:18:-35.3173,-5.41861,0:353,54,0:352.26,50.7195,0	1/1:10:-24.6612,-3.01034,0:247,30,0:245.7,26.6368,0	0/0:12:0,-3.61241,-29.2613:0,36,293:0,33.5701,293.526

# all snps are homozygous reference or homozygous alternate

# to extract the header from the vcf file:
tabix -H AllLoci_wild_annotated_output.vcf.gz
# i think i will just put this into excel and summarize that way, once i have the header this should be simpler. 


#------------------------------Looking for interesting genes within the window around Optix

## Doing this analysis after presenting data to Reed lab in September 2024
## subset Heliconius lativitta genome from Lewis 2019 (Hel_final_2016.fa from Martik Chatterjee to see if interesting stuff blasts to our genome
pwd
/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/Heliconius_lativitta

# Here are the coordinates Martik Chatterjee sent me, represent genes of interest in the lativitta genome that we would be interested to see if they align to our genome, 
# specifically the region downstream of optix:
LR1: 168283-170346
LR2: 174119-175527
obs132-131000-133624
obs214: 213873-214663
optix(promoter):49339-49966
u1:52112-52435
u2:53231-53681

# Want to add ~1kb on either side for each of these genes, which are all on chromosome 18_1

LR1: 169283-171346
LR2: 175119-176527
obs132-132000-134624
obs214: 214873-215663
optix(promoter):50339-50966
u1:53112-53435
u2:54231-54681

Hel_LR1.bed
Hel_chr18_1	169283	171346

Hel_LR2.bed
Hel_chr18_1	175119	176527

Hel_obs132.bed
Hel_chr18_1	132000	134624

Hel_obs214.bed
Hel_chr18_1	214873	215663

Hel_optix_promoter.bed
Hel_chr18_1	50339	50966

Hel_u1.bed
Hel_chr18_1	53112	53435

Hel_u2.bed
Hel_chr18_1	54231	54681


cat subset_fasta_Hel.qsub
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N subset_fasta_Hel  # job name, anything you want
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -pe omp 8

module load bedtools
bedtools getfasta -fi Hel_final_2016.fa -bed Hel_LR1.bed -fo Hel_LR1.fasta
bedtools getfasta -fi Hel_final_2016.fa -bed Hel_LR2.bed -fo Hel_LR2.fasta
bedtools getfasta -fi Hel_final_2016.fa -bed Hel_obs132.bed -fo Hel_obs132.fasta
bedtools getfasta -fi Hel_final_2016.fa -bed Hel_obs214.bed -fo Hel_obs214.fasta
bedtools getfasta -fi Hel_final_2016.fa -bed Hel_optix_promoter.bed -fo Hel_optix_promoter.fasta
bedtools getfasta -fi Hel_final_2016.fa -bed Hel_u1.bed -fo Hel_u1.fasta
bedtools getfasta -fi Hel_final_2016.fa -bed Hel_u2.bed -fo Hel_u2.fasta



# Now blast the fasta files containing genes of interest against our L. lorquini genome
# turns out the probe sequences are nucleotide seqs, so just use a regular blastn
cat blast_hel_seqs.qsub
#!/bin/bash -l
#$ -P mullenl
#$ -N blast_hel_seqs # job name, anything you want
#$ -m bea
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o blast_hel_seqs.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 10

module load blast+

blastn -query Hel_LR1.fasta -db /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curatedfastadb -out Hel_LR1_blast.csv -outfmt "7"
blastn -query Hel_obs132.fasta -db /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curatedfastadb -out Hel_obs132_blast.csv -outfmt "7"
blastn -query Hel_optix_promoter.fasta -db /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curatedfastadb -out LR1_blast.csv -outfmt "7"
blastn -query Hel_u2.fasta -db /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curatedfastadb -out Hel_optix_promoter_blast.csv -outfmt "7"
blastn -query Hel_LR2.fasta -db /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curatedfastadb -out Hel_LR2_blast.csv -outfmt "7"
blastn -query Hel_obs214.fasta -db /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curatedfastadb -out Hel_obs214_blast.csv -outfmt "7"
blastn -query Hel_u1.fasta -db /projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curatedfastadb -out Hel_u1_blast.csv -outfmt "7"




#------------------------------NgsAdmix
# use NgsAdmix to make admixture plots with the output from angsd

# first, install:
pwd
/projectnb/mullenl/hannah/limenitis_genomes/software

wget popgen.dk/software/download/NGSadmix/ngsadmix32.cpp 

# all set, downloaded, now configure

g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix

# NgsAdmix for K from 2 to 6 : do not run if the dataset contains clones or genotyping replicates!

# submitting as a job bc i think it will take awhile
cat ngsadmix
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N ngsadmix  # job name, anything you want
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -pe omp 10
module load angsd

for K in `seq 2 6` ;
do
/projectnb/mullenl/hannah/limenitis_genomes/software/NGSadmix -likes AllLoci2.beagle.gz -K $K -P 10 -o AllLoci2_k${K};
done

# Removed ANGSD1/ folder because it was taking up a ton of space



#------------------------------Linkage Disequilibrium Decay
# helpful info in the same tutorial as above for Fst, linked here: https://speciationgenomics.github.io/ld_decay/
mkdir ld_decay
pwd
/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/ld_decay

cat plink_LDscan.qsub
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N LD_scan  # job name, anything you want
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o LD_scan.qlog
#$ -pe omp 16

module load plink

VCF=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_allsamps_annotated_output.vcf

plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5 \
--thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out all_loci_ld


# then run python script on the output to calculate linkage across set distances
ld_decay_calc.py -i all_loci_ld.ld.gz -o all_loci_ld


# and look at just wild caught individuals:
cat plink_LDscan_wild.qsub
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N LD_scan_wild  # job name, anything you want
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o LD_scan_wild.qlog
#$ -pe omp 16

module load plink

VCF=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_wild_annotated_output.vcf.gz

plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5 \
--thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out wild_loci_ld


# then run python script on the output to calculate linkage across set distances
ld_decay_calc.py -i wild_loci_ld.ld.gz -o wild_loci_ld


# update with new thresholds based on advice from Nick VanKuren:
# I think you want to focus in on pairs 0 - 10 kb apart - LD in butterflies decays by 0.5 over just a few hundred bp. 
# This has been true in every butterfly I've looked at. So I would look at smaller bins closer together to get a better idea of the decay rate.
# --ld-window specifies the MAXIMUM number of variants between the two for which you're calculating LD. 
# If that's set to 100 you're probably only going to calculate for pairs less than about 5 kb apart (i.e. most pairs 5 kb apart are going to have more than 100 other variant sites between them). 
# So I wouldn't even set this - the default is infinity in plink2, so no threshold there. I would also set --ld-window-kb to 50 or even 10 because it decays so quickly.
# --ld-window-kb - this is the upper end of the LD window. If it is set to 1000, means that we ignore any two sites more than 1 Mb apart in the genome.

cat plink_LDscan_wild2.qsub
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N LD_scan_wild2  # job name, anything you want
#$ -l h_rt=48:00:00 #maximum run time
#$ -M hannahaichelman@gmail.com #your email
#$ -m be
#$ -j y # Join standard output and error to a single file
#$ -o LD_scan_wild2.qlog
#$ -pe omp 16

module load plink

VCF=/projectnb/mullenl/hannah/limenitis_genomes/snp_calling/bwamap_curatedfasta_pe/AllLoci_wild_annotated_output.vcf.gz

plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5 \
--thin 0.1 -r2 gz --ld-window-kb 50 \
--ld-window-r2 0 \
--make-bed --out wild_loci_ld_50kbwindow

plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5 \
--thin 0.1 -r2 gz --ld-window-kb 10 \
--ld-window-r2 0 \
--make-bed --out wild_loci_ld_10kbwindow


# then run python script on the output to calculate linkage across set distances

ld_decay_calc.py -i wild_loci_ld_50kbwindow.ld.gz -o wild_loci_ld_50kbwindow


ld_decay_calc.py -i wild_loci_ld_10kbwindow.ld.gz -o wild_loci_ld_10kbwindow





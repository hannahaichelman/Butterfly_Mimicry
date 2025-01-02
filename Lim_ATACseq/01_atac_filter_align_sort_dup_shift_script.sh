#!/bin/bash -l
#$ -P mullenl
#$ -N atac_step1 # job name, anything you want
#$ -m bea
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o atac_step1.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 16

##atac seq data analysis: from raw fastq files to aligned, fitered bam and bigwig files

# load modules needed
module load trimmomatic
module load bowtie2
module load samtools
module load picard
module load python3/3.8.6
module load deeptools/3.5.1 # alignmentSieve and bamCoverage

# set path where indexed genome lives
export GENOME_PATH="/projectnb/mullenl/hannah/limenitis_genomes/reference_genomes/limenitis_lorquini_genome_raw_sequence/hifi_filt/curated.fasta"
 
# process from raw files to bigwig files in for loop
for fq in /projectnb/mullenl/hannah/limenitis_atacseq/*R1.fastq.gz;
do
	sample=`basename ${fq%_*}`
	trim_sample=`basename ${fq%???????????????????????????????}` #number of question marks determines how many characters to drop from end of string 
	aligned_sample=`basename ${trim_sample}_align_bowtie2`
	aligned_sorted_sample=`basename ${aligned_sample}_sort`
	aligned_sorted_marked=`basename ${aligned_sorted_sample}_marked`
	aligned_sorted_filtered=`basename ${aligned_sorted_marked}_filtered`
	aligned_sorted_filtered_shifted=`basename ${aligned_sorted_filtered}_shifted`
	aligned_sorted_filtered_shifted_sort=`basename ${aligned_sorted_filtered}_shifted_sort`
	aligned_sorted_filtered_shifted_sort_bw=`basename ${aligned_sorted_filtered}_shifted_sort_rpkm`

	echo "trimming ${trim_sample}"
	#trim with trimmomatic 
	trimmomatic PE /projectnb/mullenl/hannah/limenitis_atacseq/${sample}_R1.fastq.gz /projectnb/mullenl/hannah/limenitis_atacseq/${sample}_R2.fastq.gz -baseout "${trim_sample}.fq.gz" \
	 ILLUMINACLIP:/projectnb/mullenl/hannah/limenitis_atacseq/NexteraPE-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

	echo "aligning ${trim_sample}" 
	#run bowtie2 alignemnt 
	bowtie2 -t -p 8 --local --no-mixed --no-discordant -I 25 -X 800 -x $GENOME_PATH -1 ${trim_sample}_1P.fq.gz \
	 -2 ${trim_sample}_2P.fq.gz | samtools view -bS - > ${aligned_sample}.bam

	echo "sorting ${trim_sample}"
	#sort by coordinate 
	samtools sort -o ${aligned_sorted_sample}.bam -@ 16 ${aligned_sample}.bam

	#create index file
	samtools index ${aligned_sorted_sample}.bam

	echo "marking duplicates ${trim_sample}"
	#mark duplicates 
	picard MarkDuplicates I=${aligned_sorted_sample}.bam O=${aligned_sorted_marked}.bam Remove_Duplicates=false \
	  METRICS_FILE=${aligned_sorted_marked}_dup_metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT

	echo "removing multi-mapped reads ${trim_sample}"  
	#remove multi-mapped reads and duplicates. The flag -F 1548 removes duplicates, reads mapped to multple genomic locations, platform-failed reads, and unmapped reads.
	samtools view -h -b -f 2 -F 1548 -q 30 ${aligned_sorted_marked}.bam | samtools sort -o ${aligned_sorted_filtered}.bam

	samtools index ${aligned_sorted_filtered}.bam

	echo "shifting reads ${trim_sample}"
	#shift reads 
	alignmentSieve --numberOfProcessors max --ATACshift --bam ${aligned_sorted_filtered}.bam -o ${aligned_sorted_filtered_shifted}.bam

	samtools sort ${aligned_sorted_filtered_shifted}.bam -@ 16 -o ${aligned_sorted_filtered_shifted_sort}.bam

	samtools index -M *_shifted_sort.bam

	echo "converting to bw ${trim_sample}"
	#convert bam to bigwig for easier visualization
	bamCoverage --numberOfProcessors 16 --binSize 10 --normalizeUsing RPKM --bam ${aligned_sorted_filtered_shifted_sort}.bam -o ${aligned_sorted_filtered_shifted_sort_bw}.bw
done 
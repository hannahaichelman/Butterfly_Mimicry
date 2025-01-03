#!/bin/bash -l
#$ -P mullenl
#$ -N atac_peakcall # job name, anything you want
#$ -m bea
#$ -M hannahaichelman@gmail.com
#$ -j y # Join standard output and error to a single file
#$ -o atac_peakcall.qlog
#$ -l h_rt=48:00:00
#$ -pe omp 16


#atac peak calling and pooling script


conda activate atac

for bam in /projectnb/mullenl/hannah/limenitis_atacseq/*align_bowtie2_sort_marked_filtered_shifted_sort*.bam
do
	sample_bam=`basename ${bam%????}`
	echo "calling peaks for sample ${sample_bam}"
	#Call peaks
	macs3 hmmratac -f BEDPE -l 10 -u 20 -c 2 -n ${sample_bam}_hhmratac -i ${sample_bam}.bed --outdir macs3/hhmratac/${sample_bam} 2> ${sample_bam}_macs3_hhmratac.log

done
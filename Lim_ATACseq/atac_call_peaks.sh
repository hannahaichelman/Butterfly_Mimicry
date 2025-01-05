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

module load miniconda
conda activate atac


## call peaks 
#rename files so that they have sample descriptor and replicate number (i.e. 20hr_treatment_r1.bam)
# need to do this after we figure out who is who...
#echo "renaming"
#sed 's/"//g' jc_atac_file_rename.txt | while IFS=, read orig new; do mv "$orig" "$new"; done 

for bam in /projectnb/mullenl/hannah/limenitis_atacseq/*filtered_shifted_sort.bam
do
	sample_bam=`basename ${bam%????}`
	echo "calling peaks for sample ${sample_bam}"
	#Convert the bam file to BEDPE
	macs3 randsample -i ${sample_bam}.bam -f BAMPE -p 100 -o ${sample_bam}.bed

	#Call peaks
	macs3 callpeak -f BEDPE --nomodel --nolambda --shift -37 --extsize 73 -g 360605465 -B  --keep-dup all --cutoff-analysis -n ${sample_bam} -t ${sample_bam}.bed --outdir macs3/${sample_bam} 2> ${sample_bam}_macs3.log

	##FriP calculation
	#convert BED (the peaks) to SAF
	echo "calculate FRIP ${sample_bam}"
	awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' macs3/${sample_bam}/${sample_bam}_peaks.narrowPeak > ${sample_bam}_peaks.saf
	#count
	featureCounts -p -a ${sample_bam}_peaks.saf -F SAF -o ${sample_bam}-readCountInPeaks.txt ${sample_bam}.bam 2> ${sample_bam}_featureCounts.log

done

##pool peaks - needs to be done after we figure out who is who...likely in separate script.

while IFS= read -r line; do
	sample_peaks="$line" 
	echo "Pooling peaks for ${sample_peaks}"
	macs3 callpeak -q 0.01 -f BEDPE --nomodel --shift -37 --extsize 73 -B -g 360605465 --keep-dup all --cutoff-analysis -n ${sample_peaks}_pooled -t *${sample_peaks}_r1.bed *${sample_peaks}_r2.bed *${sample_peaks}_r3.bed --outdir macs3/${sample_peaks}_pooled 2> macs3_${sample_peaks}_pooled.log
done < /local/workdir/noah/JC/tissue_culture/results/atac/jc_atac_sample_names.txt

##intersect pooled peaks with individual replicates to create a file with only peaks present in all three replicates 
# also need to do this after we figure out who is who...likely in pool peaks script.

while IFS= read -r line; do
	sample_intersect="$line" 
	echo "Intersecting peaks for ${sample_intersect}"
	#First extract pooled peaks which are in replicate 1
	intersectBed -wo -a macs3/${sample_intersect}_pooled/${sample_intersect}_pooled_peaks.broadPeak -b macs3/*_${sample_intersect}_r1/*_${sample_intersect}_r1_peaks.broadPeak  | awk 'BEGIN {FS="\t" ; OFS = "\t"} {s1=$3-$2 ; s2=$12-$11; if(($19/s1 > 0.5) || ($19/s2 > 0.5)) {print $0}}' | cut -f 1-9 > tmp.bed
	intersectBed -wo -a tmp.bed -b macs3/*${sample_intersect}_r2/*${sample_intersect}_r2_peaks.broadPeak  | awk 'BEGIN {FS="\t" ; OFS = "\t"} {s1=$3-$2 ; s2=$12-$11; if(($19/s1 > 0.5) || ($19/s2 > 0.5)) {print $0}}' | cut -f 1-9 > tmp2.bed
	intersectBed -wo -a tmp2.bed -b macs3/*${sample_intersect}_r3/*${sample_intersect}_r3_peaks.broadPeak  | awk 'BEGIN {FS="\t" ; OFS = "\t"} {s1=$3-$2 ; s2=$12-$11; if(($19/s1 > 0.5) || ($19/s2 > 0.5)) {print $0}}' | cut -f 1-9 > ${sample_intersect}_pooled_replicated_broadPeak.bed
done < /local/workdir/noah/JC/tissue_culture/results/atac/jc_atac_sample_names.txt
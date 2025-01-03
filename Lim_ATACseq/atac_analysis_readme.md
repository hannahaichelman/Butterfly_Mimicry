#ATAC analysis README - adapted from Noah Brady
# other good info here: https://github.com/CebolaLab/ATAC-seq?tab=readme-ov-file#pre-alignment-qc

##Step 1: QC raw fasta files. 

# Use fastqc to check quality of raw reads

##Step 2: trim adapters from raw fasta files, align reads to genome, sort by coordinate, filter multi-mapped and duplicate reads, shift reads, and produce bigwig files.

# Run the following script for step 2: 01_atac_filter_align_sort_dup_shift_script.sh
# pwd: /projectnb/mullenl/hannah/limenitis_atacseq


##Step 4: post alignment QC 

# run  atac_ac_script.R

##Step 5: call peaks and pool peaks

# to run this step, need to create conda environment and install macs3
module load miniconda
conda create --name atac
conda activate atac
conda install -c bioconda macs3
conda install bioconda::subread

# then run atac_call_peaks.qsub 

##Step 6: quantify differential accessability 
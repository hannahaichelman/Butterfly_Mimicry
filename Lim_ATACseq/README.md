# ATAC analysis README - adapted from Noah Brady

Other good info about ATACseq pipelines found at these links:
https://github.com/CebolaLab/ATAC-seq?tab=readme-ov-file#pre-alignment-qc
https://ivanek.github.io/analysisOfGenomicsDataWithR/12_ATACSeq_html.html
https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/atacseq_tutorials.html
https://nf-co.re/atacseq/2.1.2/

## Step 1: QC raw fasta files. 

Use fastqc to check quality of raw reads

## Step 2: trim adapters from raw fasta files, align reads to genome, sort by coordinate, filter multi-mapped and duplicate reads, shift reads, and produce bigwig files.

Run the following script for step 2: 01_atac_filter_align_sort_dup_shift_script.sh
Run in this working directory on BU's scc: /projectnb/mullenl/hannah/limenitis_atacseq


## Step 3: post alignment QC 

Run the atac_ac_script.R as an interactive session on BU's scc: /projectnb/mullenl/hannah/limenitis_atacseq

## Step 4: call peaks and pool peaks

To run this step, need to create conda environment and install macs3 and subread
```
module load miniconda
conda create --name atac
conda activate atac
conda install -c bioconda macs3
conda install bioconda::subread
```

Then, run atac_call_peaks.qsub, which uses the conda environment we just created (wd: /projectnb/mullenl/hannah/limenitis_atacseq)

## Step 5: quantify differential accessability 
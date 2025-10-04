
#### Set-Up ####
# Install the DESeq2 library
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(vegan)
library(EnhancedVolcano)




#### Functions ####
# Function to extract gene names from GFF file
extract_gene_names_from_gff <- function(gff_file_path) {
  gff_data <- read.csv(gff_file_path, skip = 4, header = FALSE, sep = "\t", 
                       stringsAsFactors = FALSE, quote = "")
  gene_entries <- gff_data[gff_data$V3 == "gene", ]
  
  gene_mapping <- data.frame(
    gene_id = character(),
    gene_name = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(gene_entries)) {
    id_match <- regmatches(gene_entries$V9[i], regexpr("ID=gene:[^;]+", gene_entries$V9[i]))
    if (length(id_match) > 0) {
      gene_id <- gsub("ID=gene:", "", id_match)
      name_match <- regmatches(gene_entries$V9[i], regexpr("Name=[^;]+", gene_entries$V9[i]))
      gene_name <- if (length(name_match) > 0) {
        gsub("Name=", "", name_match)
      } else {
        gene_id
      }
      
      gene_mapping <- rbind(gene_mapping, data.frame(
        gene_id = gene_id,
        gene_name = gene_name,
        stringsAsFactors = FALSE
      ))
    }
  }
  return(gene_mapping)
}

# Function to add gene names to DESeq2 results
add_gene_names_to_results <- function(results_obj, gene_mapping) {
  results_df <- as.data.frame(results_obj)
  results_df$gene_id <- rownames(results_df)
  
  results_with_names <- results_df %>%
    left_join(gene_mapping, by = "gene_id") %>%
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))
  
  rownames(results_with_names) <- results_with_names$gene_name
  return(results_with_names)
}

#### Read in and organize data ####
# read in counts file (from Cornell sequencing facility)
# this was the first version of the raw counts that the facility sent, which only included half of the samples
#counts = read.table("Lim_RNASeq/6065D_rawCounts.txt", header = TRUE)

# newest version of the raw counts file, containing all 12 samples
counts = read.table("Lim_RNASeq/6065D-All-Samples_rawCounts.txt", header = TRUE)
head(counts)

# based on the cornell differential expression sheet linked their id's back to the sample type
# d3_fw = LL07	LL09	LL11
# d3_hw = LL08	LL10	LL12
# d6_fw = LL13	LL15	LL17
# d6_hw = LL14	LL16	LL18

# make experimental design matrix
sample = colnames(counts)
expDesign = data.frame(sample) 
expDesign$wing_type = as.factor(c("fw","fw","fw","hw","hw","hw","fw","fw","fw","hw","hw","hw"))
expDesign$age = as.factor(c("d3","d3","d3","d3","d3","d3","d6","d6","d6","d6","d6","d6"))
expDesign$wing_age = as.factor(c("fw_d3","fw_d3","fw_d3","hw_d3","hw_d3","hw_d3","fw_d6","fw_d6","fw_d6","hw_d6","hw_d6","hw_d6"))
expDesign$individual = as.factor(c("1","2","3","1","2","3","4","5","6","4","5","6"))

str(expDesign)


# explore total counts across treatment
totalCounts=colSums(counts)
totalCounts
barplot(totalCounts, col=expDesign$age, ylab="raw counts", main = "total counts")

min(totalCounts) 
max(totalCounts) 


#### Differential Expression ####
# DESeq for differential expression
dds<-DESeqDataSetFromMatrix(countData=counts, colData=expDesign, design=~wing_age) 
dds = DESeq(dds)
results = results(dds)
summary(results)

# write out normalized counts file
head(results)
norm.counts = counts(dds, normalized = TRUE) # these are the counts DESeq uses
#write.csv(norm.counts,"Lim_RNASeq/normalized_counts.csv") 

# perform rlog transformation - useful for various unsupervised clustering analyses moving forward
rlogged = rlogTransformation(dds, blind = TRUE)

# Forewing vs. Hindwing, day 3 expression comparison
# second term here is the "control"
res_wing_d3 <- results(dds, contrast=c("wing_age","fw_d3","hw_d3"))
head(res_wing_d3)

#how many FDR < 10%?
table(res_wing_d3$padj<0.1)
# 0.1=2
summary(res_wing_d3)

# Forewing vs. Hindwing, day 6 expression comparison
res_wing_d6 <- results(dds, contrast=c("wing_age","fw_d6","hw_d6"))
head(res_wing_d6)

#how many FDR < 10%?
table(res_wing_d6$padj<0.1)
# 0.1=3
summary(res_wing_d6)


# Explore day 3 vs. day 6 expression comparison, Forewing
# second term here is the "control"
res_age_fw <- results(dds, contrast=c("wing_age","fw_d6","fw_d3"))
head(res_age_fw)

#how many FDR < 10%?
table(res_age_fw$padj<0.01)
# 0.1=1880
# 0.05=1339
# 0.01=720

summary(res_age_fw)
# out of 11933 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1044, 8.7%
# LFC < 0 (down)     : 836, 7%
# outliers [1]       : 254, 2.1%
# low counts [2]     : 1850, 16%
# (mean count < 2)

# MA Plot
DESeq2::plotMA(res_age_fw, main = "Forewing, Day 6 vs. Day 3")

# Explore day 3 vs. day 6 expression comparison, Hindwing
# second term here is the "control"
res_age_hw <- results(dds, contrast=c("wing_age","hw_d6","hw_d3"))
head(res_age_hw)

#how many FDR < 10%?
table(res_age_hw$padj<0.1)
# 0.1=1989
# 0.05=1419
# 0.01=797

summary(res_age_hw)
# out of 11933 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1127, 9.4%
# LFC < 0 (down)     : 862, 7.2%
# outliers [1]       : 254, 2.1%
# low counts [2]     : 1388, 12%
# (mean count < 1)

# MA Plot
DESeq2::plotMA(res_age_hw, main = "Hindwing, Day 6 vs. Day 3")


#### Make RLD P-val table ####
#get pvals
val=cbind(res_age_fw$pvalue, res_age_fw$padj, res_age_hw$pvalue, res_age_hw$padj)
head(val)
colnames(val)=c("pval.age_fw", "padj.age_fw", "pval.age_hw", "padj.age_hw")
length(val[,1])
table(complete.cases(val))

#make rlog data and pvals table
rld=assay(rlogged)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,val)
head(rldpvals)
dim(rldpvals)
# [1] 12093     16
table(complete.cases(rldpvals))
# FALSE  TRUE 
# 2264  9829 

# write.csv(rldpvals, "Lim_RNASeq/RLDandPVALS.csv", quote=F)

# identify degs:
rldpvals = read.table("Lim_RNASeq/RLDandPVALS.csv", header = TRUE, sep = ",")
str(rldpvals)

rldpvals_sig = rldpvals %>%
  filter(padj.age_fw < 0.1)
rldpvals_sig

sig.genes = rldpvals_sig %>%
  select(X) %>%
  rename(sig_gene_name = X)
sig.genes

#sig_gene_names <- filter(gff3, grepl(paste(sig.genes$sig_gene_name, collapse='|'), gff3$V9))
#write.csv(sig_gene_names, "Lim_RNASeq/sig_degs_gene_names.csv")

#### Extract gene names from gff ####
gff3 = read.csv("Lim_RNASeq/eggnog_Hmel_result_decorate_liftoff.emapper.decorated.gff", skip = 4, header = F, sep = "\t")
View(gff3)

# Extract gene names from GFF file
gene_mapping <- extract_gene_names_from_gff("Lim_RNASeq/eggnog_Hmel_result_decorate_liftoff.emapper.decorated.gff")


# compare with gff3 file I provided to Cornell to see where the genes are


#### PCA ####
# First create a PCA data frame and calculate the variance estimated by PC1 and PC2
pcadata = DESeq2::plotPCA(rlogged, intgroup = c("wing_age","individual"), returnData = TRUE)
percentVar = round(100 * attr(pcadata, "percentVar"))

# PCA with prcomp to calculate percent variance by hand
ntop=500
# calculate the variance for each gene
rv <- rowVars(assay(rlogged))
# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
pca_500 <- prcomp(t(assay(rlogged)[select,]))
# the contribution to the total variance for each component
percentVar_prcomp <- (pca_500$sdev^2 / sum(pca_500$sdev^2 ))*100

# or just do it for all genes instead of top 500
pca = prcomp(t(assay(rlogged)), center = TRUE, scale. = FALSE)
summary_pca = summary(pca)
summary_pca


# Using adonis2 from library(vegan) we can see if there are any significant differences based on wing type
adonis2(formula = pca$x ~ wing_age + individual, data = pcadata, method = 'eu')
#           Df SumOfSqs      R2      F Pr(>F)    
# Model     7    93972 0.95877 13.288  0.001 ***
# Residual  4     4041 0.04123                  
# Total    11    98013 1.00000                  

# Plot PCA
cols_wing_age = c("hw_d3"="#bdbdbd", "fw_d3" = "#fdae6b","hw_d6" = "#636363", "fw_d6" = "#e6550d")

# Plot PC1 and PC2
pca_12 = DESeq2::plotPCA(rlogged, returnData = TRUE, intgroup = c("wing_age","individual")) %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = wing_age, shape = individual), size = 3, stroke = 1) +
  scale_colour_manual(values = cols_wing_age,
                      name = "Wing Type + Age") +
  scale_shape_manual(values = c(21, 22, 23, 7, 8, 9),
                     name = "Individual") +
  #stat_ellipse(aes(color=rlogged$wing_type), type = "t", linetype = 2, lwd = 1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #labs(title = "~ treatment, p<0.001")+
  theme_bw()
pca_12
ggsave(pca_12, file = "Lim_RNASeq/allsample_pca_12.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


#### Volcano Plot ####
# forewing, d6 vs. d3
EnhancedVolcano(
  res_age_fw,  
  lab = rownames(res_age_fw),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.01,
  FCcutoff = 2,
  title = 'Forewing: day 6 vs. day 3',
  #drawConnectors = TRUE,
  #widthConnectors = 0.75
)

# hindwing, d6 vs. d3
EnhancedVolcano(
  res_age_hw,  
  lab = rownames(res_age_hw),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.01,
  FCcutoff = 2,
  title = 'Hindwing: day 6 vs. day 3',
  #drawConnectors = TRUE,
  #widthConnectors = 0.75
)

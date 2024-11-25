
#### Set-Up ####
# Install the DESeq2 library
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(vegan)

#### Read in and organize data ####
# read in counts file (from Cornell sequencing facility)
counts = read.table("Lim_RNASeq/6065D_rawCounts.txt", header = TRUE)
head(counts)

# make experimental design matrix
sample = colnames(counts)
expDesign = data.frame(sample) 
expDesign$wing_type = as.factor(c("fw","fw","fw","hw","hw","hw"))
expDesign$individual = as.factor(c("1","2","3","1","2","3"))
str(expDesign)


# explore total counts across treatment
totalCounts=colSums(counts)
totalCounts
barplot(totalCounts, col=expDesign$wing_type, ylab="raw counts", main = "total counts")

min(totalCounts) 
max(totalCounts) 


#### Differential Expression ####
# DESeq for differential expression
dds<-DESeqDataSetFromMatrix(countData=counts, colData=expDesign, design=~wing_type+individual) #can only test for the main effects of treatment
dds = DESeq(dds)
results = results(dds)
summary(results)

# write out normalized counts file
head(results)
norm.counts = counts(dds, normalized = TRUE) # these are the counts DESeq uses
# write.csv(norm.counts,"Lim_RNASeq/normalized_counts.csv") #these are all counts, not considering treatment comparisons

# perform rlog transformation - useful for various unsupervised clustering analyses moving forward
rlogged = rlogTransformation(dds, blind = TRUE)

# Explore Forewing vs. Hindwing expression comparison
# second term here is the "control"
res_wing <- results(dds, contrast=c("wing_type","fw","hw"))
head(res_wing)

#how many FDR < 10%?
table(res_wing$padj<0.05)
# 0.1=38
# 0.05=31
# 0.01=13

summary(res_wing)
# out of 11249 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 10, 0.089%
# LFC < 0 (down)     : 28, 0.25%
# outliers [1]       : 0, 0%
# low counts [2]     : 1086, 9.7%
# (mean count < 1)


# MA Plot
DESeq2::plotMA(res_wing, main = "Forewing vs. Hindwing")


# Make RLD P-val table
#get pvals
val=cbind(res_wing$pvalue, res_wing$padj)
head(val)
colnames(val)=c("pval.fw", "padj.fw")
length(val[,1])
table(complete.cases(val))

#make rlog data and pvals table
rld=assay(rlogged)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,val)
head(rldpvals)
dim(rldpvals)
# [1] 12093     8
table(complete.cases(rldpvals))
# FALSE  TRUE 
# 1930 10163 

# write.csv(rldpvals, "Lim_RNASeq/RLDandPVALS.csv", quote=F)

# identify degs:
rldpvals = read.table("Lim_RNASeq/RLDandPVALS.csv", header = TRUE, sep = ",")
str(rldpvals)

rldpvals_sig = rldpvals %>%
  filter(padj.fw < 0.1)
rldpvals_sig


#### PCA ####
# First create a PCA data frame and calculate the variance estimated by PC1 and PC2
pcadata = DESeq2::plotPCA(rlogged, intgroup = c("wing_type","individual"), returnData = TRUE)
percentVar = round(100 * attr(pcadata, "percentVar"))
# PCA with prcomp
pca = prcomp(t(assay(rlogged)), center = TRUE, scale. = FALSE)
summary_pca = summary(pca)
summary_pca

#percentVar <- pca$sdev^2/sum(pca$sdev^2)


# Using adonis from lipcadata# Using adonis from library(vegan) we can see if there are any significant differences based on wing type
adonis2(pca$x ~ wing_type + individual, data = pcadata, method = 'eu')
#          Df SumOfSqs      R2      F Pr(>F)
# Model     1    750.6 0.02624 0.1078    0.7
# Residual  4  27858.1 0.97376              
# Total     5  28608.7 1.00000              

# Plot PCA
cols_wing = c("hw" = "#636363", "fw" = "#fd8d3c")

# Plot PC1 and PC2
pca_12 = DESeq2::plotPCA(rlogged, returnData = TRUE, intgroup = c("wing_type","individual")) %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = wing_type, shape = individual), size = 3, stroke = 1) +
  scale_colour_manual(values = cols_wing,
                      name = "Wing Type") +
  scale_shape_manual(values = c(21, 22, 23),
                     name = "Individual") +
  #stat_ellipse(aes(color=rlogged$wing_type), type = "t", linetype = 2, lwd = 1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #labs(title = "~ treatment, p<0.001")+
  theme_bw()
pca_12
ggsave(pca_12, file = "Lim_RNASeq/pca_12.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# Plot PC1 and PC3 (separates wing types better)
pca_13 = DESeq2::plotPCA(rlogged, returnData = TRUE, intgroup = c("wing_type","individual"), pcs = c(1,3)) %>% 
  ggplot(aes(x = PC1, y = PC3)) +
  geom_point(aes(colour = wing_type, shape = individual), size = 3, stroke = 1) +
  scale_colour_manual(values = cols_wing,
                      name = "Wing Type") +
  scale_shape_manual(values = c(21, 22, 23),
                     name = "Individual") +
  #stat_ellipse(aes(color=rlogged$wing_type), type = "t", linetype = 2, lwd = 1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
  #labs(title = "~ treatment, p<0.001")+
  theme_bw()
pca_13
ggsave(pca_13, file = "Lim_RNASeq/pca_13.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

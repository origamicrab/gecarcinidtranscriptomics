#deRNAseq with deseq2

#install relevant packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximport")
#BiocManager::install("DESeq2")
#BiocManager::install('EnhancedVolcano')
#install.packages("tidyverse")
#install.packages("readr") #increases tximport's speed
#install.packages("Rcpp")
#install.packages("pheatmap")
#install.packages("ggfortify")
#install.packages('gplots')
#install.packages('RColorBrewer')

#load libraries
library(tximport)
library(Rcpp)
library(DESeq2)
library(tidyverse)
library(readr)

#read in sample index file
samples <- read.csv("tc_samples_GL_acute.csv")

#import gene level counts into R using tximport
txi <- tximport(files = samples$quant_file, type = "salmon", txOut = TRUE)

#summarize the new file
summary(txi)

#make column headers more informative
colnames(txi$counts) <- samples$sample

#make sure they make sense now
head(txi$counts)

#transform our txi object into something deseq can work with. error messages about factor conversion are fine!
dds <- DESeqDataSetFromTximport(txi = txi, colData = samples, design = ~time_point)

#make sure the baseline is the first level in the treatment factor
dds$time_point <- relevel(dds$time_point, "T0")

#how many reads are there per contig? Save this to import into python after dropping all rows that have only zeroes
#Using python script, select only transcripts that have more than 0 copies in at least three individuals per treatment

#keep <- rowSums(counts(dds)) > 0
#dds_init <- dds[keep,]
#write.csv(counts(dds_init), "TCGL_acute_contig_counts.csv")

#import a list of all transcripts that have more than 0 copies in at least three individuals in one treatment
goodTranscripts <- read.csv("TCGL_kept.csv", header = FALSE)

index <- match(goodTranscripts[,1],rownames(counts(dds)))

dds <- dds[index,]

#perform differential expression
dds <- DESeq(dds)

#output matrix of normalized counts
normcounts <- counts(dds, normalized = TRUE) #save this as a csv/rds matrix
#write.csv(normcounts, "tc_norm_reads_GL_acute.csv")

#examine results
res <- results(dds)
head(res)

#sort and filter our output based on adjusted p-value and log2foldchange
res_sig <- subset(res, padj<.05)
res_lfc <- subset(res_sig, abs(log2FoldChange) > 1) 
head(res_lfc)

#for Fisher table construction, create dataframe of insignificant genes and significant ones with lfc <1 and export
res_insig <- subset(res, padj>.05)
res_siglowlfc <- subset(res_sig, abs(log2FoldChange) < 1)

res_others <- rbind(res_insig, res_siglowlfc)

#QC: save res_lfc to compare to single copy orthologs from WGCNA and see histogram of l2fc values
#write.csv(res_lfc, "lfc_celeste_GL_acute_DEG.csv")
#write.csv(res_others, "lfc_celeste_GL_acute_nonDEG.csv")

#plot volcano distribution of all genes
par(mfrow=c(1,1))
plotMA(res, ylim = c(-15, 15))

#create publication ready volcano plots
library(EnhancedVolcano)

EnhancedVolcano(res, 
                lab = rownames(res), 
                x = 'log2FoldChange', 
                y = 'pvalue', xlim = c(-15,15),
                pCutoff = 10e-4,
                FCcutoff = 2, #higher here than what we used for rest of analysis (|l2fc| < 1)
                title = 'T. celeste Posterior Gill',
                subtitle = 'Baseline (T0) vs. Acute (T1) Expression',
                caption = bquote(~Log[2]~ "fold change cutoff: 2, p-value cutoff: 10e-4"),
                labSize = .1, #I don't want gene labels because these aren't informative yet without annotation. Made them v. small
                legendLabSize = 16,
                legendIconSize = 5.0) + coord_flip() # will show plot horizontally instead of vertically

#plotDispEsts( dds, ylim = c(1e-6, 1e4))

#what gene has the most significant difference in gene expression?
plotCounts(dds, gene=which.min(res$padj), intgroup = "time_point")


hist( res$pvalue, breaks=20, col="grey" )

#regularized log transformation of our count data
rld <- rlog(dds)
head(assay(rld))

par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )

#variance stabilized transformation of our count data
vsd <- vst(dds)

#plot. Calculate sample distances and plot them
sample_dists <- assay(vsd) %>%
  t() %>%
  dist() %>%
  as.matrix() 

head(sample_dists)

#from tutorial
#distances between samples
sampleDists <- dist(t(assay(rld)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$time_point,
                                     rld$sample, sep="-" )
colnames(sampleDistMatrix) <- NULL

library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

#calculate MDS values from the distance matrix
mdsData <- data.frame(cmdscale(sample_dists))
mds <- cbind(mdsData, as.data.frame(colData(vsd))) # combine with sample data
head(mds)

library(ggfortify)
library(cluster)

#plot PCA with ggplot2 with 95% confidence intervals
ggplot(mds, aes(X1, X2, color = time_point, shape = time_point)) + 
  geom_point(size = 4) +
  theme_classic(base_size = 20) + 
  theme(legend.position= "none") +
  theme(panel.grid.major = element_line(colour = "lightgrey", size = rel(0.75)),
        panel.grid.minor = element_line(colour = "lightgrey", size = rel(0.25))) +
  labs(y= 'PC2', x = 'PC1') +
  geom_polygon(stat = "ellipse", aes(fill = time_point), alpha = 0.3)

#generate a heatmap of the data
library(pheatmap)

#select 25 significant genes with the highest log2fold change
genes <- order(res_lfc$log2FoldChange, decreasing = TRUE) [1:25]

#make a data.frame that contains info about samples that will appear in the heatmap
annot_col <- samples %>%
  column_to_rownames('sample') %>%
  select(time_point) %>%
  as.data.frame()

head(annot_col)

#for creating count table by expression change, how many genes are up or down regulated?
up <- subset(res_lfc, log2FoldChange > 0)
down <- subset(res_lfc, log2FoldChange < 0)

#plot the heatmap -- there is a difference between using vsd and rld here. Plot with both and see which is best
pheatmap(assay(rld)[genes, ], cluster_rows = TRUE, show_rownames = TRUE, cluster_cols = FALSE, annotation_col = annot_col, scale="row")

#to generate heatmap of candidate genes only from supplement
#pull out transcripts of interest
cand_genes <- c('TRINITY_DN20685_c0_g2_i1', 'TRINITY_DN22631_c4_g1_i2', 'TRINITY_DN23469_c3_g1_i2', 'TRINITY_DN23608_c2_g1_i2', 
           'TRINITY_DN23619_c5_g4_i1')

#make a data.frame that contains info about samples that will appear in the heatmap
annot_col <- samples %>%
  column_to_rownames('sample') %>%
  select(time_point) %>%
  as.data.frame()

head(annot_col)

#plot the heatmap
anno.colors <- list(time_point = c(T0 = "white", T1 = "lightgray", T3 = "darkgray", T4 = "black")) #choose nicer colors for timepoints
pheatmap(assay(vsd)[cand_genes, ], cluster_rows = TRUE, show_colnames = FALSE, show_rownames = TRUE, cluster_cols = FALSE, 
         annotation_colors = anno.colors, annotation_col = annot_col, scale="row", 
         color = brewer.pal(8, "BrBG")) #palette to match the other heatmap in paper

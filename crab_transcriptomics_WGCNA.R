library(tximport)
library(dplyr)
library(ggplot2)
library(DESeq2)

#Read in quant data by species

###### CELESTE #########
#read in all count files
tc_samples <- read.csv("tc_samples_ALL.csv")
txi_celeste <- tximport(files = tc_samples$quant_file, type = "salmon", txOut = TRUE)
#summarize the new file
summary(txi_celeste)
#make column headers more informative
colnames(txi_celeste$abundance) <- tc_samples$sample
#make sure they make sense now
head(txi_celeste$abundance)

###### MAGNA #########
#read in all count files
tm_samples <- read.csv("tm_samples_ALL.csv")
txi_magna <- tximport(files = tm_samples$quant_file, type = "salmon", txOut = TRUE)
#summarize the new file
summary(txi_magna)
#make column headers more informative
colnames(txi_magna$abundance) <- tm_samples$sample
#make sure they make sense now
head(txi_magna$abundance)

###### NATALIS #########
#read in all count files
gn_samples <- read.csv("gn_samples_ALL.csv")
txi_natalis <- tximport(files = gn_samples$quant_file, type = "salmon", txOut = TRUE)
#summarize the new file
summary(txi_natalis)
#make column headers more informative
colnames(txi_natalis$abundance) <- gn_samples$sample
#make sure they make sense now
head(txi_natalis$abundance)

##Read in ortholog group table
ortho.counts <- read.delim("Orthogroups.GeneCount.csv")

#how many are present in only one copy in all three species?
singles <- which(ortho.counts$celeste_pep==1 & ortho.counts$magna_pep==1 & ortho.counts$natalis_pep==1)
length(singles) # These are 1:1:1 orthologs. There are 2185 of them.

#how many are present in one or more copies in all three species?
nozeros <- which(ortho.counts$celeste_pep>0 & ortho.counts$magna_pep>0 & ortho.counts$natalis_pep>0)
length(nozeros) # These are all groups with at least one in each species. There are 9516 of them.

##Read in orthogroups
orthos <- read.delim("Orthogroups.csv") #Had to find/replace "__" with "::" in Excel before reading in
#write.csv(orthos, "ortho_checkFeb12022.csv")
#single.orthos <- orthos[singles,]
#write.csv(single.orthos, "single_orthos.csv") #export this output and make sure all transcript names match salmon outputs, drop index row, and rename OG column

##combine one-to-one orthologs into one frame
#read in edited singles csv
single.orthos <- read.csv("single_orthos.csv")
singles.celeste <- txi_celeste$abundance[match(single.orthos$celeste_pep,rownames(txi_celeste$abundance)),]
singles.magna <- txi_magna$abundance[match(single.orthos$magna_pep,rownames(txi_magna$abundance)),]
singles.natalis <- txi_natalis$abundance[match(single.orthos$natalis_pep,rownames(txi_natalis$abundance)),]
all.singles <- cbind(singles.celeste,singles.magna,singles.natalis)
rownames(all.singles) <- single.orthos[,1]
#write.csv(all.singles, "all.singles_norm_ALL.csv")

#load WGCNA library
library(WGCNA)

##Read in gene expression and metadata
ge <- read.csv("all.singles_norm_ALL.csv", header = T) #all.singles normalized count output from DESeq2
rownames(ge) <- ge[,1]
ge[,1] <- NULL
datExpr <- t(ge)
#QC: do rowSums for each sample to see if there are different numbers of reads between tissues
#write.csv(datExpr, "datExpr_Mar14_2022.csv")


##Read in metadata
meta <- read.csv("all_samples.csv", header = T)
meta$baseline <- as.factor(meta$baseline)
meta$acute <- as.factor(meta$acute)
meta$extreme <- as.factor(meta$extreme)
meta$recovery <- as.factor(meta$recovery)
meta$celeste <- as.factor(meta$celeste)
meta$magna <- as.factor(meta$magna) 
meta$natalis <- as.factor(meta$natalis) 
meta$antennal_gland <- as.factor(meta$antennal_gland)
meta$gill <- as.factor(meta$gill)
ordermeta <- meta[match(rownames(datExpr),meta$sample),]
ordermeta <- tail(ordermeta, -1) #drop row of NA
datTraits <- ordermeta[,c("celeste", "magna", "natalis", "baseline", "acute", "extreme", "recovery", "antennal_gland", "gill")]
rownames(datTraits) <- ordermeta$sample

##Filter for genes with zero/very low expression
low.thresh=0.1
#up.thresh=10000
means <- colMeans(datExpr)
keep <- which(means>low.thresh)
datExpr <- datExpr[,keep]

sampleTree = hclust(dist(datExpr), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)  #DC17 is an outlier

# Plot a line to show the cut
abline(h = 375, col = "red");

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 375, minSize = 10)
table(clust) # clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

## Soft Thresholding 
options(stringsAsFactors = FALSE)

#choose a set of soft thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

#call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

#scale free topology fit index as a function of the soft thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
#this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
#mean connectivity as a function of the soft thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#(I ended up choosing 20 -- lowest power for which the scale-free topology fit index reaches 0.9)

softPower = 20;
adjacency = adjacency(datExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

#we like large modules, so we set a min module size of 30
minModuleSize = 30;

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#merge modules with similar expression
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicMods)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#choose .25, corresponds to correlation of .75
MEDissThres = 0.25

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicMods, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
# moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

## Relating modules to traits 
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
# MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

dim(textMatrix) = dim(moduleTraitCor)
par(mfrow=c(1,1),mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#these colors are horrific and painful to look at, so I ended up remaking this figure in Python for the publication

#save corr info to make prettier corr matrix 
write.csv(moduleTraitCor, file = "module_trait_corr_matrix.csv")
write.csv(moduleTraitPvalue, file = 'module_trait_p_value.csv')

# Define variable weight containing the species column of datTrait FOR EACH SPP
#celeste_first
celeste_MM <- as.data.frame(datTraits$celeste)
names(celeste_MM) = 'celeste'

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, celeste_MM, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(celeste_MM), sep="");
names(GSPvalue) = paste("p.GS.", names(celeste_MM), sep="");

# Create the starting data frame
geneInfo0 = data.frame(moduleGenes = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for species
modOrder = order(-abs(cor(MEs, celeste_MM, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleGenes, -abs(geneInfo0$GS.celeste));
geneInfo = geneInfo0[geneOrder, ]

#save gene info table. Has OG name by its respective module membership for each module. 
#We used this table to try to figure out module functions based on genes with greates abs. module membership values
#make sure to save this before moving to next species as it will overwrite these variables with new info
write.csv(geneInfo, file = "geneInfo_celeste.csv")





####   MAGNA
magna_MM <- as.data.frame(datTraits$magna)
names(magna_MM) = 'magna'

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, magna_MM, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(magna_MM), sep="");
names(GSPvalue) = paste("p.GS.", names(magna_MM), sep="");

# Create the starting data frame
geneInfo0 = data.frame(moduleGenes = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for species
modOrder = order(-abs(cor(MEs, magna_MM, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleGenes, -abs(geneInfo0$GS.magna));
geneInfo = geneInfo0[geneOrder, ]

#save gene info table
write.csv(geneInfo, file = "geneInfo_magna.csv")







#### NATALIS
natalis_MM <- as.data.frame(datTraits$natalis)
names(natalis_MM) = 'natalis'

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, natalis_MM, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(natalis_MM), sep="");
names(GSPvalue) = paste("p.GS.", names(natalis_MM), sep="");

# Create the starting data frame
geneInfo0 = data.frame(moduleGenes = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for species
modOrder = order(-abs(cor(MEs, natalis_MM, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleGenes, -abs(geneInfo0$GS.natalis));
geneInfo = geneInfo0[geneOrder, ]

#save gene info table
write.csv(geneInfo, file = "geneInfo_natalis.csv")


###### BASELINE


baseline_MM <- as.data.frame(datTraits$baseline)
names(baseline_MM) = 'baseline'

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, baseline_MM, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(baseline_MM), sep="");
names(GSPvalue) = paste("p.GS.", names(baseline_MM), sep="");

# Create the starting data frame
geneInfo0 = data.frame(moduleGenes = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for species
modOrder = order(-abs(cor(MEs, baseline_MM, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleGenes, -abs(geneInfo0$GS.baseline));
geneInfo = geneInfo0[geneOrder, ]

#save gene info table
write.csv(geneInfo, file = "geneInfo_baseline.csv")


##### ACUTE

acute_MM <- as.data.frame(datTraits$acute)
names(acute_MM) = 'acute'

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, acute_MM, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(acute_MM), sep="");
names(GSPvalue) = paste("p.GS.", names(acute_MM), sep="");

# Create the starting data frame
geneInfo0 = data.frame(moduleGenes = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for species
modOrder = order(-abs(cor(MEs, acute_MM, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleGenes, -abs(geneInfo0$GS.acute));
geneInfo = geneInfo0[geneOrder, ]

#save gene info table
write.csv(geneInfo, file = "geneInfo_acute.csv")




###### EXTREME
extreme_MM <- as.data.frame(datTraits$extreme)
names(extreme_MM) = 'extreme'

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, extreme_MM, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(extreme_MM), sep="");
names(GSPvalue) = paste("p.GS.", names(extreme_MM), sep="");

# Create the starting data frame
geneInfo0 = data.frame(moduleGenes = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for species
modOrder = order(-abs(cor(MEs, extreme_MM, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleGenes, -abs(geneInfo0$GS.extreme));
geneInfo = geneInfo0[geneOrder, ]

#save gene info table
write.csv(geneInfo, file = "geneInfo_extreme.csv")







###### RECOVERY

recovery_MM <- as.data.frame(datTraits$recovery)
names(recovery_MM) = 'recovery'

# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, recovery_MM, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(recovery_MM), sep="");
names(GSPvalue) = paste("p.GS.", names(recovery_MM), sep="");

# Create the starting data frame
geneInfo0 = data.frame(moduleGenes = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for species
modOrder = order(-abs(cor(MEs, recovery_MM, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleGenes, -abs(geneInfo0$GS.recovery));
geneInfo = geneInfo0[geneOrder, ]

#save gene info table
write.csv(geneInfo, file = "geneInfo_recovery.csv")

# save modules to create eigengene expression plots
# edit in Excel to create a table with columns (sample name, species, time point, tissue type, and 
# [module numbers with eigengene expression values in the cells])
write.csv(MEs, file = "MEs_ALL.csv")

#for each module, use linear model to see p-values for various combinations of spp and time points and interaction terms
Lm <- lm(MEs[,3]~datTraits$baseline*datTraits$celeste)
summary(Lm) #take summary to see p-value for each factor and interaction between factors
#make a table of these for p-value x module x trait and interaction

#choose hubgene in each module
hubs = chooseTopHubInEachModule(datExpr, colnames(MEs))
hubs

#pull out ME OGs for enrichment analysis. Can export these for post-annotation downstream analysis.
ME0 <- names(datExpr[1,])[moduleColors=="0"]
ME1 <- names(datExpr[1,])[moduleColors=="1"]
ME2 <- names(datExpr[1,])[moduleColors=="2"]
ME3 <- names(datExpr[1,])[moduleColors=="3"]
ME4 <- names(datExpr[1,])[moduleColors=="4"]
ME5 <- names(datExpr[1,])[moduleColors=="5"]
ME6 <- names(datExpr[1,])[moduleColors=="6"]
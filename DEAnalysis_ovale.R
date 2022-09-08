library(ggplot2)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(edgeR)
library(pcaMethods)
library(mixOmics)
library(gprofiler2)
library(GO.db)

##Read in counts, remove globin genes
human_counts.txt <- read.delim("/Volumes/projects-t3/SerreDLab-3/kieran.tebben/mali/ovale/count_tables/human_counts.txt.gz", comment.char="#")
globin_genes <- c("HBB", "HBD", "HBA1", "HBA2", "HBE1", "HBEGF", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
human_counts.txt <- human_counts.txt[!(human_counts.txt$Geneid %in% globin_genes),]
#counts <- human_counts.txt[,c(7,8,10:15)]
counts <- human_counts.txt
colnames(counts) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length", "109_10_29_09", "109_2_10_11", "357_11_4_10", "357_3_26_11", "357_5_26_11", "357_8_30_10", "371_12_28_10", "371_3_12_11", "371_9_9_11")
rownames(counts) <- human_counts.txt$Geneid
counts <- counts[,-c(9)]

##Make DGEList
dgList <- DGEList(counts=counts[,c(7:14)], genes=counts[c(1:6)])

##Filter for only genes expressed at at least 10 cpm in >50% of samples
keep <- rowSums(cpm(dgList)>10) >= 4
table(keep)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]
summary(cpm(dgList))

##Normalization with TMM 
dgList <- calcNormFactors(dgList, method="TMM")
#dgList_cpm <- cpm(dgList)

##Variance partition
library(variancePartition)
gene_count <- data.frame(cpm(dgList$counts))

###Read in gene metadata for each sample
all <- read.delim("po_allvars.txt")
all <- all[all$Sample %in% colnames(counts),]
all <- all[ order(match(all$Sample, colnames(counts))), ]
all$Sex <- all$gender
all$Individual <- factor(all$pid)
all$Species <- c("pf", "po", "po", "po", "pf", "pf", "po", "pf")
all$Parasitemia <- log(all$parasitemia)
all$Age <- all$age_infection

###Define formula
form <- ~ (1|Species) + (1|Individual) + Parasitemia + Age

###Create model
varPart <- fitExtractVarPartModel(gene_count, form, all)
df <- data.frame(varPart)
df$gene <- rownames(df)
df_long <- pivot_longer(df, !gene, names_to = "variable", values_to = "Percent")
df_long$Percent <- 100 * df_long$Percent

###Sort variables by median fraction of variance explained
vp <- sortCols( varPart )
plotVarPart(vp, text.size = 5)

ggplot(df_long) + geom_violin(aes(x = variable, y = Percent))

##Data exploration - MDS
plotMDS(dgList, gene.selection = "common")

##PCA 
counts_tmm <- cpm(dgList)
d_t <- as.data.frame(t(counts_tmm))
d_t$species <- c("pf", "po", "po", "po", "pf", "pf", "po", "pf")
d_t$pid <- c("A1", "A2", "B1", "B2", "B3", "C1", "C2", "C3")
d_t$sex <- c("F", "F", "F", "F", "F", "M", "M", "M")
d_t$age <- c("<5", "<5", ">5", ">5", ">5", ">5", ">5", ">5")

plot(tune.pca(d_t, ncomp = 8, center = TRUE, scale = FALSE))

result <- pca(d_t, ncomp = 3, center = TRUE, scale = FALSE)
result

plotIndiv(result, 
          #Change ground variable to color by different category
          group = d_t$age, 
          style = "ggplot2", 
          col.per.group = c("blue", "red"), 
          legend = T, cex = 7, 
          ind.names = d_t$pid,
          size.legend = 20,
          size.legend.title = 20,
          size.xlabel = 20,
          size.ylabel = 20,
          size.axis = 20,
          title = "PCA")

##PCA of only Po samples - all genes, not orthologs##
counts_po <- counts[,c(8,9,10,13)]
dgList <- DGEList(counts=counts_po, genes=rownames(counts_po))

#Filter for only genes expressed at at least 10 cpm
keep <- rowSums(cpm(dgList)>10) >= 4
table(keep)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]
summary(cpm(dgList))

dgList <- calcNormFactors(dgList, method="TMM")

#PCA of ovale samples only
counts_tmm <- cpm(dgList)
d_t <- as.data.frame(t(counts_tmm))
d_t$species <- c("P. ovale curtisi", "P. ovale curtisi", "P. ovale curtisi", "P. ovale wallikeri")
d_t$pid <- c("A2", "B1", "B2", "C2")

plot(tune.pca(d_t, ncomp = 4, center = TRUE, scale = FALSE))

result <- pca(d_t, ncomp = 3, center = TRUE, scale = FALSE)
result

plotIndiv(result, 
          #Change group to color by different variable
          group = d_t$species, 
          style = "ggplot2", 
          ind.names = d_t$pid,
          col.per.group = c("blue", "red"), 
          legend = T, cex = 7, 
          size.legend = 20,
          size.legend.title = 20,
          size.xlabel = 20,
          size.ylabel = 20,
          size.axis = 20,
          title = "PCA")

##Design Matrix
d <- dgList[,c(1,2,3,5,7,8)]

pid <- factor(c("109", "109", "357", "357", "371", "371"))
species <- factor(c("pf", "po", "po", "pf", "po", "pf"))
parasitemia <- c(16950, 2025, 11375, 600, 4825, 11950)

design <- model.matrix(~pid+species)
rownames(design) <- colnames(d)
design

##Estimate dispersion 
y <- estimateDisp(d, design, robust=TRUE) #Estimates common, trended and tagwise together
#y <- estimateGLMCommonDisp(dgList, design, verbose=TRUE)
#y <- estimateGLMTrendedDisp(y, design)
#y <- estimateGLMTagwiseDisp(y, design)

##Differential expression analysis
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
summary(decideTests(lrt))

results <- as.data.frame(lrt$table)
results$FDR <- p.adjust(results$PValue, method="BH")


results_up <- subset(results, logFC > 0)
results_up <- subset(results_up, FDR < 0.1)
results_up$logpval <- -log10(results_up$PValue)
#gp_up <- gost(row.names(results_up), organism = "hsapiens")

results_down <- subset(results, logFC < 0)
results_down <- subset(results_down, FDR < 0.1)

##Volcano plot
results$deexpressed <- ifelse(results$FDR <= 0.1 & results$logFC >= 0, "UP", 
                              ifelse(results$FDR <= 0.1 & results$logFC < 0, "DOWN", "NO"))

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results, aes(x = logFC, y = -log10(PValue), col = deexpressed)) + 
  geom_point() + scale_colour_manual(values = mycolors, label = c("Down", "Up", "Not significant"), name = "Key") + theme_minimal() + theme(text = element_text(size = 25))

################ PLASMODIUM #############
##Read in Pf data and subset to only orthologous genes and primary pf samples
pf_counts.txt <- read.delim("/Volumes/projects-t3/SerreDLab-3/kieran.tebben/mali/ovale/count_tables/pf_counts.txt.gz", comment.char="#")
counts_pf <- pf_counts.txt[,c(7,8,10:15)]
colnames(counts_pf) <- c("109_10_29_09", "109_2_10_11", "357_3_26_11", "357_5_26_11", "357_8_30_10", "371_12_28_10", "371_3_12_11", "371_9_9_11")
rownames(counts_pf) <- pf_counts.txt$Geneid

orthologs <- read.delim("pf_po_orthologs.txt")
rownames(orthologs) <- orthologs$Pf_geneID

counts_pf <- merge(orthologs, counts_pf, by = 0)
counts_pf <- counts_pf[,c(2,4,8,9,11)]
rownames(counts_pf) <- counts_pf$Pf_geneID
counts_pf <- counts_pf[,-c(1)]

##Read in Po data and subset to only orthologous genes and primary po infections
po_counts.txt <- read.delim("/Volumes/projects-t3/SerreDLab-3/kieran.tebben/mali/ovale/count_tables/po_counts.txt.gz", comment.char="#")
counts_po <- po_counts.txt[,c(7,8,10:15)]
colnames(counts_po) <- c("109_10_29_09", "109_2_10_11", "357_3_26_11", "357_5_26_11", "357_8_30_10", "371_12_28_10", "371_3_12_11", "371_9_9_11")
rownames(counts_po) <- po_counts.txt$Geneid

rownames(orthologs) <- orthologs$Po_geneID
counts_po <- merge(orthologs, counts_po, by = 0)
counts_po <- counts_po[,c(2,5,6,7,10)]
rownames(counts_po) <- counts_po$Pf_geneID
counts_po <- counts_po[,-c(1)]

counts <- merge(counts_po, counts_pf, by = 0)
rownames(counts) <- counts$Row.names
counts <- counts[,-c(1)]

##Make DGEList
dgList <- DGEList(counts=counts, genes=rownames(counts))

##Filter for only genes expressed at at least 10 cpm
keep <- rowSums(cpm(dgList)>10) >= 4
table(keep)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]
summary(cpm(dgList))

dgList <- calcNormFactors(dgList, method="TMM")

#dgList_cpm <- cpm(dgList)

plotMDS(dgList, gene.selection="common", col = as.numeric(dgList$samples$group), labels = dgList$samples$group)

###PCA
counts_tmm <- cpm(dgList)
d_t <- as.data.frame(t(counts_tmm))
d_t$species <- c("po", "po", "po", "po", "pf", "pf", "pf", "pf")
d_t$pid <- c("A2", "B1", "B2", "C2", "A1", "B3", "C1", "C3")

plot(tune.pca(d_t, ncomp = 8, center = TRUE, scale = FALSE))

result <- pca(d_t, ncomp = 3, center = TRUE, scale = FALSE)
result

plotIndiv(result, 
          group = d_t$species, 
          style = "ggplot2", 
          ind.names = d_t$pid,
          col.per.group = c("blue", "red"), 
          legend = T, cex = 7, 
          size.legend = 20,
          size.legend.title = 20,
          size.xlabel = 20,
          size.ylabel = 20,
          size.axis = 20,
          title = "PCA")

plotIndiv(result, 
          comp = c(1,3),
          group = d_t$species, 
          style = "ggplot2", 
          col.per.group = c("blue", "red"), 
          legend = T, cex = 7, 
          ind.names = rownames(d_t),
          size.legend = 20,
          size.legend.title = 20,
          size.xlabel = 20,
          size.ylabel = 20,
          size.axis = 20)

##PCA of only Po samples - all genes, not orthologs##
po_counts.txt <- read.delim("/Volumes/projects-t3/SerreDLab-3/kieran.tebben/mali/ovale/count_tables/po_counts.txt.gz", comment.char="#")
counts_po <- po_counts.txt[,c(7,8,10:15)]
colnames(counts_po) <- c("109_10_29_09", "109_2_10_11", "357_3_26_11", "357_5_26_11", "357_8_30_10", "371_12_28_10", "371_3_12_11", "371_9_9_11")
rownames(counts_po) <- po_counts.txt$Geneid
counts_po <- counts_po[,c(2,3,4,7)]

dgList <- DGEList(counts=counts_po, genes=rownames(counts_po))

#Filter for only genes expressed at at least 10 cpm
keep <- rowSums(cpm(dgList)>10) >= 4
table(keep)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]
summary(cpm(dgList))

dgList <- calcNormFactors(dgList, method="TMM")

#PCA
counts_tmm <- cpm(dgList)
d_t <- as.data.frame(t(counts_tmm))
d_t$species <- c("P. ovale curtisi", "P. ovale curtisi", "P. ovale curtisi", "P. ovale wallikeri")
d_t$pid <- c("A2", "B1", "B2", "C2")

plot(tune.pca(d_t, ncomp = 3, center = TRUE, scale = FALSE))

result <- pca(d_t, ncomp = 3, center = TRUE, scale = FALSE)
result

plotIndiv(result, 
          group = d_t$species, 
          style = "ggplot2", 
          ind.names = d_t$pid,
          col.per.group = c("blue", "red"), 
          legend = T, cex = 7, 
          size.legend = 20,
          size.legend.title = 20,
          size.xlabel = 20,
          size.ylabel = 20,
          size.axis = 20,
          title = "PCA")


##Variance partition
gene_count <- cpm(dgList)

###Read in gene metadata for each sample
all <- read.delim("po_allvars.txt")
all <- all[all$Sample %in% colnames(counts),]
all <- all[ order(match(all$Sample, colnames(counts))), ]
all$Sex <- all$gender
all$Individual <- factor(all$pid)
all$Species <- c("po", "po", "po", "po", "pf", "pf", "pf", "pf")
all$Parasitemia <- log(all$parasitemia)
all$Age <- all$age_infection

###Define formula
form <- ~ (1|Species) + (1|Individual) + Parasitemia + Age

###Create model
varPart <- fitExtractVarPartModel(gene_count, form, all)

###Sort variables by median fraction of variance explained
vp <- sortCols( varPart )
plotVarPart(vp, text.size = 5)

###Variance partition with just Pf samples###
gene_count_pf <- cpm(dgList)
gene_count_pf <- gene_count_pf[,c(5:8)]

all <- read.delim("po_allvars.txt")
all <- all[all$Sample %in% colnames(counts),]
all <- all[ order(match(all$Sample, colnames(counts))), ]
all$Sex <- all$gender
all$Individual <- factor(all$pid)
all$Species <- c("po", "po", "po", "po", "pf", "pf", "pf", "pf")
all$Parasitemia <- log(all$parasitemia)
all$Age <- all$age_infection
all$Stage <- all$Trophs

##See if Proportion of trophs correlates with PC2 of Pf samples
pf <- subset(all, Species == "pf")
ggscatter(pf, x = "Trophs", y = "PC2", 
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.coef.size = 10, cor.method = "pearson",
          xlab = "Proportion of Trophozoites", ylab = "PC2", size = 5) + theme(text = element_text(size = 30))

###Sort variables by median fraction of variance explained
vp <- sortCols( varPart )
plotVarPart(vp, text.size = 5)

##Differential expression - no covariates (compare species)
species <- as.factor(c("po", "po", "po", "po", "pf", "pf", "pf", "pf"))

design <- model.matrix(~species)
design

##Estimate dispersion and calculate differential expression 
y <- estimateDisp(dgList, design, robust = TRUE)

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
summary(decideTests(lrt))

results <- as.data.frame(lrt$table)
results$FDR <- p.adjust(results$PValue, method="BH")

results_up <- subset(results, logFC > 0)
results_up <- subset(results_up, FDR < 0.1)

results_down <- subset(results, logFC < 0)
results_down <- subset(results_down, FDR < 0.1)

##Volcano plot
results$deexpressed <- ifelse(results$FDR <= 0.1 & results$logFC >= 0, "UP", 
                              ifelse(results$FDR <= 0.1 & results$logFC < 0, "DOWN", "NO"))

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results, aes(x = logFC, y = -log10(PValue), col = deexpressed)) + 
  geom_point() + scale_colour_manual(values = mycolors) + theme_minimal()

##Differential expression - adjusted for stage composition
#Make dgList
dgList <- DGEList(counts=counts, genes=rownames(counts))

##Filter for only genes expressed at at least 10 cpm in >4/8 samples
keep <- rowSums(cpm(dgList)>10) >= 4
table(keep)
dgList <- dgList[keep, , keep.lib.sizes=FALSE]
summary(cpm(dgList))

dgList <- calcNormFactors(dgList, method="TMM")

#Design matrix including stage composition
#Correct for ratio of rings to trophs? or correct for either rings, trophs and schizonts
ciber <- read.csv("plasmodium_deconvolution.csv")

group <- as.factor(c("po", "po", "po", "po", "pf", "pf", "pf", "pf"))
rings <- c(42.08, 54.53, 52.75, 78.05, 5.66, 15.95, 32.11, 33.47)
trophs <- c(40.95, 43.27, 26.40, 9.72, 78.75, 65.83, 45.37, 21.89)
schiz <- c(9.26, 1.15, 4.96, 1.94, 0.00, 9.08, 3.72, 19.59)
male <- c(1.54, 0.63, 4.13, 4.73, 8.57, 5.29, 14.79, 19.13)
female <- c(6.17, 0.43, 11.76, 5.55, 7.01, 3.85, 4.01, 5.92)

#Remove covariates from design matrix to run unadjusted model
design <- model.matrix(~group+rings+trophs+schiz+female)
design

##Estimate dispersion and calculate differential expression 
y <- estimateDisp(dgList, design, robust = TRUE)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
summary(decideTests(lrt))

results <- as.data.frame(lrt$table)
results$FDR <- p.adjust(results$PValue, method="BH")


results_up <- subset(results, logFC > 0)
results_up <- subset(results_up, FDR < 0.1)

results_down <- subset(results, logFC < 0)
results_down <- subset(results_down, FDR < 0.1)

##Plot results
results$deexpressed <- ifelse(results$FDR <= 0.1 & results$logFC >= 0, "UP", 
                              ifelse(results$FDR <= 0.1 & results$logFC < 0, "DOWN", "NO"))

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(results, aes(x = logFC, y = -log10(PValue), col = deexpressed, label = rownames(results))) + 
  geom_point() + scale_colour_manual(values = mycolors, label = c("Down", "Up", "Not DE"), name = "Key") + theme_minimal() + theme(text = element_text(size = 25))

# ##DE analysis adjusting for PC2 instead of stage composition 
# #Design matrix including PC2
# results <- data.frame(result$variates$X)
# loadings <- data.frame(result$loadings$X)
# 
# group <- group <- as.factor(c("po", "po", "po", "po", "pf", "pf", "pf", "pf"))
# PCs <- data.frame(cbind(results$PC1, results$PC2, results$PC3))
# colnames(PCs) <- c("PC1", "PC2", "PC3")
# 
# ggplot(PCs) + geom_point(aes(x = PC1, y = PC3))
# 
# design <- model.matrix(~PC2+group)
# rownames(design) <- rownames(dgList$samples)
# design
# 
# y <- estimateDisp(dgList, design, robust = TRUE)
# 
# fit <- glmFit(y, design)
# lrt <- glmLRT(fit)
# topTags(lrt)
# summary(decideTests(lrt))
# 
# results <- as.data.frame(lrt$table)
# results$FDR <- p.adjust(results$PValue, method="BH")
# 
# 
# results_up <- subset(results, logFC > 0)
# results_up <- subset(results_up, FDR < 0.1)
# 
# results_down <- subset(results, logFC < 0)
# results_down <- subset(results_down, FDR < 0.1)
# 
# ##Plot results
# results$deexpressed <- ifelse(results$FDR <= 0.1 & results$logFC >= 0, "UP", 
#                               ifelse(results$FDR <= 0.1 & results$logFC < 0, "DOWN", "NO"))
# 
# mycolors <- c("blue", "red", "black")
# names(mycolors) <- c("DOWN", "UP", "NO")
# 
# ggplot(results, aes(x = logFC, y = -log10(PValue), col = deexpressed)) + 
#   geom_point() + scale_colour_manual(values = mycolors) + theme_minimal()

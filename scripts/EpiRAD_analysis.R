##################################################################
## This script reads in a text file derived from the "Base Counts"
## .vcf output from ipyrad. Base counts are used for analysis of 
## EpiRADseq data.

# Read in data file
data <- read.delim("data3-2.txt", header=FALSE)
# Since the base counts were split into four columns, these need 
# to be summed
data2 <- t(sapply(seq(4,ncol(data), by=4), function(i) {
     indx <- i:(i+3)
     rowSums(data[indx[indx <= ncol(data)]])}))

#The resulting file needs to be transposed and turned into a dataframe
data3 <- as.data.frame(t(data2))
#Add column with locus number (CHROM from .vcf file)
locus <- data[,1]
row.names(data3) <- locus
#Add header names
header <- read.delim("header_data3.txt", header=FALSE)
names <- as.vector(t(header))
colnames(data3) <- names
#Select samples of interest (some have very low sample sizes)
data4 <- data3[,c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                  24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,
                  41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56)]

#Remove ddr rows that have any zeros. The premise here is that zeros 
#in the EpiRAD dataset are informative because they may reflect 
#methylation, but they could also relfect true absence of the locus
#in the library. Here the ddRAD library serves to standarize the EpiRAD
#library. Any zeros in the ddRAD libary are treated as absence of the
#locus, thereby leaving zeros in the EpiRAD library only where the 
#locus was counted in the ddRAD library.

data5 <- data4[apply(data4[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,
                             31,33,35,37,39,41,43,45,47,49)],1,
                     function(z) !any(z==0)),] 


#################################################################
# Now use edgeR package to standardize count data by library size

library("edgeR")
#read in the file to edgeR
counts <- DGEList(counts=data5)
counts$samples
#TMM normalization (corrects for library size)
counts2 <- calcNormFactors(counts)
counts2$samples
#extract normalized counts
counts2_cpm <- cpm(counts2, normalized.lib.sizes=TRUE, log=TRUE)

##Plots to show ddRAD vs EpiRAD library (before normalization)
par(mfrow = c(5, 5))
par(mar = c(2, 2 ,2 ,2), oma = c(4, 4, 0.5, 0.5))

for (i in seq(1,49, by = 2)){
  plot(data5[,i], data5[,i+1], main = colnames(data5[i]))
}

#plot normalized counts
par(mfrow = c(5, 5))
par(mar = c(2, 2, 2, 2), oma = c(4, 4, 0.5, 0.5)) 

for (i in seq(1,49, by = 2)){
  plot(counts2_cpm[,i], counts2_cpm[,i+1], main = colnames(counts2_cpm[i]))
}


##################################################################
#Using lm to get residuals

models <- list()
for (i in seq(1,49, by = 2)){
  models[[colnames(counts2_cpm)[i]]] <- lm(counts2_cpm[,i+1] ~ counts2_cpm[,i])
}

residuals <- lapply(models, '[[', 2)
resid_all <- as.data.frame(residuals)  

#plot residuals
par(mfrow = c(5, 5))
par(mar = c(2,2, 2, 2), oma = c(4, 4, 0.5, 0.5)) 

for (i in 1:25){
  plot(resid_all[,i])
}

#Plot to compare raw data to residuals
par(mfrow = c(2, 1))
par(mar = c(4, 4.5, 2, 1), oma = c(1, 1, 0, 0))
plot(data5[,13], data5[,14], xlab = "ddRAD read counts", ylab = "EpiRAD read counts", col = "blue")
plot(resid_all[,7], ylab = "Residual", col = "blue")

########################################################
#Read in sample info
sinfo <- read.table("sample_info.txt", colClasses = 'character', header = FALSE)
#without sample 101 and 112
sinfo2 <- sinfo[c(2:10,12:27),]
#transpose
tsinfo <- t(sinfo2)
#create vectors for morphotype and diameter (note whether sample 101
#was included or not)
type <- tsinfo[3,]
diam <- tsinfo[5,]

########################################################
#MDS
resid_t <- t(resid_all)
d <- dist(resid_t) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

#This adds a vector of color values based on the diam values
data_seq = seq(min(as.numeric(diam)), max(as.numeric(diam)), length=25)
col_pal = colorRampPalette(c('blue', 'green', 'red'))(25+1)
cols = col_pal[ cut(as.numeric(diam), data_seq, include.lowest=T) ]


layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric  MDS",  type="n")
points(x,y, col = cols, pch = 19)

legend_image <- as.raster(matrix(sort(cols, decreasing = TRUE), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'branch diameter')
text(x=1.5, y = seq(0,1,l=5), labels = seq(6,26,l=5))
rasterImage(legend_image, 0, 0, 1,1)


########################################################
#PCA
pca <- prcomp(resid_t)
summary(pca)
eig <- (pca$sdev)^2

#plot
biplot(pca)

# Correlation between variables and principal components
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}
# Variable correlation/coordinates
loadings <- pca$rotation
sdev <- pca$sdev
var.coord <- var.cor <- t(apply(loadings, 1, var_cor_func, sdev))
head(var.coord[, 1:4])

# Plot the correlation circle
a <- seq(0, 2*pi, length = 100)
plot( cos(a), sin(a), type = 'l', col="gray",
      xlab = "PC1",  ylab = "PC2")
abline(h = 0, v = 0, lty = 2)
# Add active variables
arrows(0, 0, var.coord[, 1], var.coord[, 2], 
       length = 0.1, angle = 15, code = 2)
# Add labels
text(var.coord, labels=rownames(var.coord), cex = 1, adj=1)

#################################################################
#Make binary dataset based on residuals <=-1

resid_all_binary <- ifelse(resid_all<=-1, 1, 0)

#proportion of methylated cutsites
prop_methyl <- colSums(resid_all_binary) / nrow(resid_all_binary)
barplot(prop_methyl)
dens <- density(prop_methyl)
plot(dens)

#Get only rows that are differentially methylated
resid1 <- resid_all_binary[rowSums(resid_all_binary) < 25, ]
resid2 <- resid1[rowSums(resid1) >= 1, ]

#MDS
resid_t_binary <- t(resid2)
d <- dist(resid_t_binary) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]

#This adds a column of color values
# based on the y values
data_seq = seq(min(as.numeric(diam)), max(as.numeric(diam)), length=25)
col_pal = colorRampPalette(c('blue', 'green', 'red'))(25+1)
cols = col_pal[ cut(as.numeric(diam), data_seq, include.lowest=T) ]
text(x, y, labels = row.names(resid_t_binary), 
     col = cols, cex=.7)
points(x,y, col = cols, pch = 19)


layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
         main="Metric  MDS",  type="n")
points(x,y, col = cols, pch = 19)

legend_image <- as.raster(matrix(sort(cols, decreasing = TRUE), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'branch diameter')
text(x=1.5, y = seq(0,1,l=5), labels = seq(6,26,l=5))
rasterImage(legend_image, 0, 0, 1,1)

#heatmap

heatmap(resid_t_binary, scale = "none")

####################################################################
#DAPC using adegenet

library("adegenet")
library("ade4")

#Find optimal number of clusters irrespective of species id
#In this case, best to retain all PCs
groups <- find.clusters(resid_t_binary, max.n.clust=24)

#Cross validation to determine number of PCs to retain
xval <- xvalDapc(resid_t_binary, type, n.pca.max = 24, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)

#show max number of PCs to retain
xval[2:6]

#perform dapc using groups defined above group (groups$grp)
dapc1 <- dapc(resid_t_binary, grp = type)
scatter(dapc1)

compoplot(dapc1, posi="bottomright",
          txt.leg=paste("Cluster", 1:2), lab="",
          ncol=1, xlab="individuals", col=funky(2), 
          label.inds = TRUE)

#perform dapc using coral type as group (grp)
dapc1 <- dapc(a6, grp = type[2:27])
scatter(dapc1)



#look at loadings
set.seed(4)
contrib <- loadingplot(dapc1$var.contr, axis=1, 
                       thres=.015, lab.jitter=1)


####################################################################
#DESeq2 analysis of differentially methylated loci
#NOTE: this does not seem to return useful information (it is skewed
#towards high count data, though maybe this could be dealt with)

library(DESeq2)
# Specify which columns are in which groups
sinfo <- read.delim("data3_sample_info.txt", header=FALSE)
deseq2.colData <- data.frame(condition=factor(sinfo[2:27,3]), 
                             type=factor(rep("paired-end", 26)))
rownames(deseq2.colData) <- colnames(data5)
deseq2.dds <- DESeqDataSetFromMatrix(countData = data5,
                                     colData = deseq2.colData, 
                                     design = ~ condition)
# Run Analysis
deseq2.dds <- DESeq(deseq2.dds)
deseq2.res <- results(deseq2.dds)
deseq2.res <- deseq2.res[order(rownames(deseq2.res)), ]

head(deseq2.res)

# Count number of hits with adjusted p-value less then 0.05
dim(deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ])

tmp <- deseq2.res
# The main plot
plot(tmp$baseMean, tmp$log2FoldChange, pch=20, cex=0.45, ylim=c(-3, 3), log="x", col="darkgray",
     main="Diff. methylyated loci  (pval <= 0.05)",
     xlab="mean of normalized counts",
     ylab="Log2 Fold Change")

tmp.sig <- deseq2.res[!is.na(deseq2.res$padj) & deseq2.res$padj <= 0.05, ]
points(tmp.sig$baseMean, tmp.sig$log2FoldChange, pch=20, cex=0.45, col="chartreuse")
# 2 FC lines
abline(h=c(-1,1), col="blue")

#Get normalized counts and plot PCA
rld <- rlog(deseq2.dds)
plotPCA(rld, intgroup=c("condition", "type"))

counts2_cpm <- counts(deseq2.dds, normalized = TRUE)

sig_counts <- deseq_normcounts[c("11617", "445", "6418","6943", "791", "958"), ]


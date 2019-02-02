############################################################################
## This script analyzes ddRADseq and EpiRADseq data in Porites astreoides


#Set directory
setwd("~/Documents/Projects/PoritesRADseq/P.ast-transplant/analyses/ipyrad_analysis/data2_outfiles") 

#Load libraries
library(ape, quietly = TRUE)
library(reshape2, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(edgeR, quietly = TRUE)
library(adegenet, quietly = TRUE)
library(ade4, quietly = TRUE)
library(hierfstat, quietly = TRUE)
library(vegan, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(superheat, quietly = TRUE)
library(limma, quietly = TRUE)


## Read in a .geno file from the ipyrad output and extract 
## ddRAD data without missing values. This is the best filetype for MDS.

geno1 <- read.table("data2.u.geno", colClasses = 'character', header = FALSE)
geno2 <- read.fwf("data2.u.geno", widths=rep(1, max(nchar(geno1$V1))), colClasses = 'numeric', header=FALSE)
header <- read.delim("header.txt", header=FALSE)
names <- t(header)
names2 <-as.vector(names)
colnames(geno2) <- names2

#########################################################################
#Matrix with ddr and epi loci for comparison of SNP genotyping error
#:44,59:70,73:76,83:90
geno6 <- geno2[,c(37:38)] #subset of samples of interest
geno7 <- geno6[!rowSums(geno6 == 9) >= 1,] #this removes missing values for a complete dataset
geno4 <- t(geno7)
#create distanc ematrix
epidd_dist <- dist.gene(geno4, method = "percent", pairwise.deletion = FALSE,
                        variance = FALSE)
epidd_dist2 <- as.matrix(epidd_dist)
epidd_dist3 <- epidd_dist2[c(seq(from =1, to = nrow(epidd_dist2), by= 2)), 
                           c(seq(from =2, to = ncol(epidd_dist2), by= 2))]

melted <- melt(epidd_dist3, na.rm = TRUE)

ggplot(data = melted, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "blue", limit = c(0,0.25), space = "Lab", 
                       name="SNP 
Mismatches") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9.5, hjust = 1))+
          labs(x= "EpiRADseq samples", y = "ddRADseq samples")+
  coord_fixed()

#mean genotyping error (mean of diagonal) 
mean(diag(epidd_dist3)) #0.013
sd(diag(epidd_dist3)) #0.008

#mean resampling error
resamp1 <- epidd_dist3[row(epidd_dist3) == (col(epidd_dist3) - 1)][c(TRUE, FALSE)]
resamp2 <- epidd_dist3[row(epidd_dist3) == (col(epidd_dist3) + 1)][c(TRUE, FALSE)]
mean(c(resamp1,resamp2)) #0.017
sd(c(resamp1,resamp2)) #0.008

#mean genetic distance for non replicates and resamples
mean(melted$value[melted$value >0.1]) #0.162
sd(melted$value[melted$value >0.1]) #0.015

##############################################################################
## Next we read in a text file derived from the "Read Depth"
## .vcf output from ipyrad (v.0.5.15). Read counts are used for analysis of 
## EpiRADseq data.

#Read in data file. The file "data3-2.txt" was generated from the notebook "VCF_readcounts.ipynb"
Epidata <- read.delim("out.DP.FORMAT", header=TRUE)

#Use aggregate to get means for duplicate records of CHROM (locus ID)
Epidata2 <- aggregate(.~CHROM, data=Epidata, mean)
rownames(Epidata2) <- Epidata2$CHROM

Epidata3 <- Epidata2[,c(39:46,61:72,75:78,85:92)]



#Remove ddr rows that have any zeros. The premise here is that zeros 
#in the EpiRAD dataset are informative because they may reflect 
#methylation, but they could also relfect true absence of the locus
#in the library. Here the ddRAD library serves to standarize the EpiRAD
#library. Any zeros in the ddRAD libary are treated as absence of the
#locus, thereby leaving zeros in the EpiRAD library only where the 
#locus was counted in the ddRAD library.

Epidata4 <- Epidata3[apply(Epidata3[c(seq(1, 32, by = 2))],1,
                     function(z) !any(z<=15)),] #increased from z==0


#################################################################
# Now use edgeR package to standardize EpiRAD count data by library size

#read in the file for edgeR
counts <- DGEList(counts=Epidata4)
counts$samples
#TMM normalization (corrects for library size)
counts2 <- calcNormFactors(counts)
counts2$samples
#extract normalized counts
counts2_cpm <- cpm(counts2, normalized.lib.sizes=TRUE, log=TRUE)

##Plots to show ddRAD vs EpiRAD library (before normalization)
par(mfrow = c(5, 5))
par(mar = c(2, 2 ,2 ,2), oma = c(4, 4, 0.5, 0.5))

for (i in seq(1,ncol(Epidata4), by = 2)){
  plot(Epidata4[,i], Epidata4[,i+1], main = colnames(Epidata4[i]), col = "blue")
}


#plot normalized counts
par(mfrow = c(5, 5))
par(mar = c(2, 2, 2, 2), oma = c(4, 4, 0.5, 0.5)) 

for (i in seq(1,ncol(counts2_cpm), by = 2)){
  plot(counts2_cpm[,i], counts2_cpm[,i+1], main = colnames(counts2_cpm[i]), col = "blue")
}

 
##################################################################
#Using lm to get residuals

models <- list()
for (i in seq(1,ncol(counts2_cpm), by = 2)){
  models[[colnames(counts2_cpm)[i]]] <- lm(counts2_cpm[,i+1] ~ counts2_cpm[,i])
}

residuals <- lapply(models, '[[', 2)
resid_all <- as.data.frame(residuals)  

#plot residuals
par(mfrow = c(5, 5))
par(mar = c(2,2, 2, 2), oma = c(4, 4, 0.5, 0.5)) 

for (i in 1:ncol(resid_all)){
  plot(resid_all[,i], col = "blue", ylim = c(-10, 4))
}

#Use k-means clustering to obtain methylated and unmethylated loci
set.seed(1234)
clusters <- lapply(residuals, kmeans, 2)
clusters2 <- lapply(clusters, '[[', 1)
clusters3 <- as.data.frame(clusters2)
clusters3[,c(2,4,5,8,10,12:14)] <-ifelse(clusters3[,c(2,4,5,8,10,12:14)] == 2, 1, 2)
clusters4 <- as.matrix(clusters3)

#plot residuals with k-means colors
par(mfrow = c(5, 5))
par(mar = c(2,2, 2, 2), oma = c(4, 4, 0.5, 0.5)) 

for (i in 1:ncol(resid_all)){
  plot(resid_all[,i], col = clusters4[,i], ylim = c(-10, 4))
}


heatmap(clusters4, scale = "none")

 
#################################################################
#Make binary dataset of EpiRAD data 

resid_all_binary <- clusters4-1

#proportion of methylated cutsites
prop_methyl <- colSums(resid_all_binary) / nrow(resid_all_binary)
dens <- density(prop_methyl)
par(mfrow = c(1, 1))
par(mar = c(5, 4, 4, 2))
plot(dens, xlab = "Methylation", ylab = "Density", main = "")

#mean and sd of methylation
mean(prop_methyl) #0.181
sd(prop_methyl) #0.01

#Get rows that are differentially methylated
resid0 <- length(which(rowSums(resid_all_binary) == 0))/(length(resid_all_binary[,1])) #loci consitutively unmethylated
resid1 <- length(which(rowSums(resid_all_binary) == 1))/(length(resid_all_binary[,1])) #loci with 1 instance of methylation
resid2 <- length(which(rowSums(resid_all_binary) == 2))/(length(resid_all_binary[,1])) #loci with 2 instances of methylation
resid3 <- length(which(rowSums(resid_all_binary) == 3))/(length(resid_all_binary[,1])) #loci with 3 instances of methylation
resid4 <- length(which(rowSums(resid_all_binary) == 4))/(length(resid_all_binary[,1])) 
resid5 <- length(which(rowSums(resid_all_binary) == 5))/(length(resid_all_binary[,1])) 
resid6 <- length(which(rowSums(resid_all_binary) == 6))/(length(resid_all_binary[,1])) 
resid7 <- length(which(rowSums(resid_all_binary) == 7))/(length(resid_all_binary[,1])) 
resid8 <- length(which(rowSums(resid_all_binary) == 8))/(length(resid_all_binary[,1])) 
resid9 <- length(which(rowSums(resid_all_binary) == 9))/(length(resid_all_binary[,1])) 
resid10 <- length(which(rowSums(resid_all_binary) == 10))/(length(resid_all_binary[,1])) 
resid11 <- length(which(rowSums(resid_all_binary) == 11))/(length(resid_all_binary[,1])) 
resid12 <- length(which(rowSums(resid_all_binary) == 12))/(length(resid_all_binary[,1])) 
resid13 <- length(which(rowSums(resid_all_binary) == 13))/(length(resid_all_binary[,1])) 
resid14 <- length(which(rowSums(resid_all_binary) == 14))/(length(resid_all_binary[,1])) 
resid15 <- length(which(rowSums(resid_all_binary) == 15))/(length(resid_all_binary[,1]))
resid16 <- length(which(rowSums(resid_all_binary) == 16))/(length(resid_all_binary[,1])) 
obs <- seq(0,16,1)
meth_obs <- rbind(resid0, resid1, resid2, resid3, resid4, resid5, resid6, resid7, resid8, resid9,
                  resid10, resid11, resid12, resid13, resid14, resid15, resid16)
methobs <- cbind(meth_obs,obs)


#################################################################
#Box plots / data summary

#distance matrix of methylation differences (to get between colony per year differences)
temp <- t(resid_all_binary)
temp2 <- dist.gene(temp, method = "percent", pairwise.deletion = FALSE,
                   variance = FALSE)
dist1 <- as.matrix(temp2)
dist2 <- epidd_dist2[c(seq(from =1, to = nrow(epidd_dist2), by= 2)), 
                     c(seq(from =2, to = ncol(epidd_dist2), by= 2))]

#differences between colonies in 2015
temp <- t(resid_all_binary[,c(seq(1,16, by = 2))])
temp2 <- dist.gene(temp, method = "percent", pairwise.deletion = FALSE,
                   variance = FALSE)
dist3 <- as.matrix(temp2)
get_upper_tri <- function(dist3){
  dist3[lower.tri(dist3)]<- NA
  return(dist3)
}
upper_tri <- get_upper_tri(dist3)
dist4 <- melt(upper_tri, na.rm = TRUE)
dist5 <- dist4[!(dist4$value == 0) >= 1,]

#differences between colonies in 2016
temp <- t(resid_all_binary[,c(seq(2,16, by = 2))])
temp2 <- dist.gene(temp, method = "percent", pairwise.deletion = FALSE,
                   variance = FALSE)
dist6 <- as.matrix(temp2)
get_upper_tri <- function(dist6){
  dist6[lower.tri(dist6)]<- NA
  return(dist6)
}
upper_tri <- get_upper_tri(dist6)
dist7 <- melt(upper_tri, na.rm = TRUE)
dist8 <- dist7[!(dist7$value == 0) >= 1,]

#create dataset with 2015 and 2016 differences with year factor
year <- as.factor(rep(c("2015", "2016"), each = 28))
per_change <- rbind(dist5, dist8)
per_change2 <- cbind(per_change$value, year)

#percent of methylated CpGs with year factor
year2 <- as.factor(rep(c("2015", "2016"), 8))
prop_methyl2 <- cbind(prop_methyl,year2)


############################################################################
#heatmap of EpiRAD data 

#Matrix for heatmap
heat.mat <- cbind(prop_methyl,t(resid_all_binary))
names <- c("pa10-15", "pa10-16", "pa11-15", "pa11-16", "pa2-15", "pa2-16", "pa3-15",
           "pa3-16", "pa5-15", "pa5-16", "pa6-15", "pa6-16", "pa8-15", "pa8-16",
           "pa9-15", "pa9-16")
heat.mat2 <- t(heat.mat)
colnames(heat.mat2) <- names
heat.mat2 <- heat.mat2[,sort(colnames(heat.mat2))]


superheat(heat.mat2[2:650, c(1:16)], bottom.label.text.size = 5,
          scale = FALSE, row.dendrogram = TRUE, col.dendrogram = FALSE, 
          pretty.order.cols = FALSE, pretty.order.rows = TRUE,
          heat.pal = c("#6baed6", "#08519c"), yt = heat.mat2[1,c(1:16)],
          legend = FALSE, grid.vline.col = "white", left.label = "none",
          grid.hline.size = 0.1,bottom.label.text.angle = 90,
          yt.axis.name = "Methylation", yt.plot.type = "bar")

par(mfrow = c(2, 1))
par(mar = c(4, 4.5, 2, 1), oma = c(1, 1, 0, 0))

#Percent CpG methylation by year
boxplot(prop_methyl2[,1]*100 ~ year2, col= c("#efedf5", "#bcbddc"),  
        main = "B", cex = 2, ylab = "Percent")
#Percentage of loci by incidences of methylation
barplot(methobs[,1]*100, names.arg = methobs[,2], col = "#756bb1", 
        main = "C", ylab = "Percent", xlab = "No. of samples")

#test for homogeneity of variances differences between years
bartlett.test(prop_methyl ~ year2, data = prop_methyl2) #variances not sig. different

#paired t-test
t.test(prop_methyl ~ year2, data = prop_methyl2, var.equal = TRUE, paired = TRUE) #p-value = 0.5685


###########################################################################
#MDS of EpiRAD data

names <- c("pa10-15", "pa10-16", "pa11-15", "pa11-16", "pa2-15", "pa2-16", "pa3-15",
           "pa3-16", "pa5-15", "pa5-16", "pa6-15", "pa6-16", "pa8-15", "pa8-16",
           "pa9-15", "pa9-16")
colnames(resid_all_binary) <- names
resid_t_binary <- t(resid_all_binary)
# euclidean distances between the rows
epidist <- dist(resid_t_binary) 
epifit <- cmdscale(epidist,eig=TRUE, k=2)
epix <- epifit$points[,1]
epiy <- epifit$points[,2]

layout(matrix(c(1, 1, 2,
                1, 1, 3), nrow=2, byrow=TRUE))
layout.show(n=3)

#MDS plot
plot(epix, epiy, xlab="Coordinate 1", ylab="Coordinate 2", type = 'n', main = "A")
for (i in seq(1,17, by = 2)){
  arrows(epix[i], epiy[i], epix[i+1], epiy[i+1], length = 0.1, col = "gray")
}

names2 <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)
palette(brewer.pal(n = 8, name = "Set2"))
text(epix, epiy, labels = row.names(resid_t_binary), cex=.7, font = 2)

#Percent pairwise difference between colonies
boxplot(per_change2[,1]*100 ~ year, col= c("#efedf5", "#bcbddc"), 
        main = "B", cex = 2, ylab = "Percent")
#Percent change in 'clones' after transplantation
boxplot(diag(dist2)*100, col= "#756bb1", 
        main = "C", cex = 2, ylab = "Percent")

#test for homogeneity of variances differences between years
bartlett.test(V1 ~ year, data = per_change2) #variances not sig. different

#paired t-test
t.test(V1 ~ year, data = per_change2, var.equal = TRUE, paired = TRUE) #p-value = 0.03892
mean(per_change2[c(1:28),1]) #2015 0.0379705
sd(per_change2[c(1:28),1]) #2015 0.005844253
mean(per_change2[c(29:56),1]) #2016 0.03565926
sd(per_change2[c(29:56),1]) #2016 0.007706606
mean(diag(dist2))

#################################################################
# Venn diagram looking at shared loci across individuals

#just get diff methylated loci 
diffmeth <- resid_all_binary[which(rowSums(resid_all_binary) >= 1 ),]
diffmeth2 <- diffmeth[which(rowSums(diffmeth) < 16 ),]
first <- seq(1,ncol(diffmeth2),2)  
second <- seq(2,ncol(diffmeth2),2)  
diffmeth3 <- abs(diffmeth2[, first] - diffmeth2[, second])
diffmeth4 <- which(rowSums(diffmeth3) >=2)
diffmeth5 <- as.numeric(names(diffmeth4))-1

#################################################################
#Was the degree of methylation change associated with inshore /offshore differences?

epidist2 <- as.matrix(epidist)

epidist3 <- epidist2[row(epidist2) == (col(epidist2) - 1)]

inshoreoffshore <- as.data.frame(cbind(epidist3[c(1,3,9,11,5,7,13,15)],
                                       c("Inshore", "Inshore", "Inshore", "Inshore",
                                         "Offshore", "Offshore", "Offshore", "Offshore")))

inshoreoffshore$V1 <- as.numeric(as.character(inshoreoffshore$V1))

boxplot(inshoreoffshore$V1 ~ inshoreoffshore$V2)
bartlett.test(inshoreoffshore$V1 ~ inshoreoffshore$V2) #variances not sig. different
t.test(inshoreoffshore$V1 ~ inshoreoffshore$V2, var.equal = TRUE) #df = 6, p-value = 0.6959

#################################################################
#Compare pairwise genetic distance with pairwise epigenetic distance

snpdist <- t(geno7[,c(seq(1,32, by = 2))])

snpdist2 <- dist.gene(snpdist, method = "percent", pairwise.deletion = FALSE,
                        variance = FALSE)
snpdist3 <- as.matrix(snpdist2)

# Get lower triangle of the  matrix
get_lower_tri<-function(snpdist3){
  snpdist3[upper.tri(snpdist3)] <- NA
  return(snpdist3)
}
# Get upper triangle of the  matrix
get_upper_tri <- function(snpdist3){
  snpdist3[lower.tri(snpdist3)]<- NA
  return(snpdist3)
}

#get just the upper triangle of the matrix
upper_tri <- get_upper_tri(snpdist3)
snpdist4 <- melt(upper_tri, na.rm = TRUE)
snpdist5 <- snpdist4[!(snpdist4$value == 0) >= 1,]

#Now the same thing for methylation data

resid_diff <- t(resid_t_binary)
methdist <- t(resid_diff)
methdist2 <- dist.gene(methdist, method = "percent", pairwise.deletion = FALSE,
                        variance = FALSE)
methdist3 <- as.matrix(methdist2)

# Get lower triangle of the  matrix
get_lower_tri<-function(methdist3){
  methdist3[upper.tri(methdist3)] <- NA
  return(methdist3)
}
# Get upper triangle of the  matrix
get_upper_tri <- function(methdist3){
  methdist3[lower.tri(methdist3)]<- NA
  return(methdist3)
}

#get just the upper triangle of the matrix
upper_tri <- get_upper_tri(methdist3)
methdist4 <- melt(upper_tri, na.rm = TRUE)
methdist5 <- methdist4[!(methdist4$value == 0) >= 1,]

# #linear regression of the snp and meth data
# epi_snp_lm <- lm(snpdist5[,3] ~ methdist5[,3])
# summary(epi_snp_lm)
# 
# #linear regression without outliers
# snpdist6 <- snpdist5[which(snpdist5[,3]>0.1),]
# methdist6 <- methdist5[which(methdist5[,3]>0.159),]
# epi_snp_lm_no_out <- lm(snpdist6[,3] ~ methdist6[,3])
# summary(epi_snp_lm_no_out)

dev.off()
plot(snpdist5[,3], methdist5[,3], ylim = c(0.02, 0.07), xlim = c(0.12,0.20),col = "blue", 
     xlab = "Genetic distance", ylab = "Epigenetic distance")
# abline(epi_snp_lm, col = "orange")
# abline(epi_snp_lm_no_out, col = "green")


#######################################################################################################
# #Select samples of interest (some have very low sample sizes)
# # This is for independent data sets
# 
# Epidata_10 <- Epidata2[,c(39:42)]
# Epidata_11 <- Epidata2[,c(43:46)]
# Epidata_02 <- Epidata2[,c(61:64)]
# Epidata_03 <- Epidata2[,c(65:68)]
# Epidata_05 <- Epidata2[,c(69:72)]
# Epidata_06 <- Epidata2[,c(75:78)]
# Epidata_08 <- Epidata2[,c(85:88)]
# Epidata_09 <- Epidata2[,c(89:92)]
# 
# 
# #Remove ddr rows that have any zeros. The premise here is that zeros 
# #in the EpiRAD dataset are informative because they may reflect 
# #methylation, but they could also relfect true absence of the locus
# #in the library. Here the ddRAD library serves to standarize the EpiRAD
# #library. Any zeros in the ddRAD libary are treated as absence of the
# #locus, thereby leaving zeros in the EpiRAD library only where the 
# #locus was counted in the ddRAD library.
# 
# Epidata_10_1 <- Epidata_10[apply(Epidata_10[c(seq(1, ncol(Epidata_10), by = 2))],1,
#                                  function(z) !any(z <= 15)), ] #increased from z==0
# Epidata_11_1 <- Epidata_11[apply(Epidata_11[c(seq(1, ncol(Epidata_11), by = 2))],1,
#                                  function(z) !any(z <= 15)), ] #increased from z==0
# Epidata_02_1 <- Epidata_02[apply(Epidata_02[c(seq(1, ncol(Epidata_02), by = 2))],1,
#                                  function(z) !any(z <= 15)), ] #increased from z==0
# Epidata_03_1 <- Epidata_03[apply(Epidata_03[c(seq(1, ncol(Epidata_03), by = 2))],1,
#                                  function(z) !any(z <= 15)), ] #increased from z==0
# Epidata_05_1 <- Epidata_05[apply(Epidata_05[c(seq(1, ncol(Epidata_05), by = 2))],1,
#                                  function(z) !any(z <= 15)), ] #increased from z==0
# Epidata_06_1 <- Epidata_06[apply(Epidata_06[c(seq(1, ncol(Epidata_06), by = 2))],1,
#                                  function(z) !any(z <= 15)), ] #increased from z==0
# Epidata_08_1 <- Epidata_08[apply(Epidata_08[c(seq(1, ncol(Epidata_08), by = 2))],1,
#                                  function(z) !any(z <= 15)), ] #increased from z==0
# Epidata_09_1 <- Epidata_09[apply(Epidata_09[c(seq(1, ncol(Epidata_09), by = 2))],1,
#                                  function(z) !any(z <= 15)), ] #increased from z==0
# 
# 
# Epidata4 <- list(Epidata_10_1, Epidata_11_1, Epidata_02_1, Epidata_03_1, Epidata_05_1, 
#                  Epidata_06_1, Epidata_08_1, Epidata_09_1)
# 
# ##Plots to show ddRAD vs EpiRAD library (before normalization)
# par(mfrow = c(5, 5))
# par(mar = c(2, 2 ,2 ,2), oma = c(4, 4, 0.5, 0.5))
# 
# for (i in seq(1,ncol(Epidata4[[1]]), by = 2)){
#   for (j in 1:8){
#     plot(Epidata4[[j]][,i], Epidata4[[j]][,i+1], main = colnames(Epidata4[[j]][i]), col = "blue")
#   }
# }
# 
# # remove_outliers <- function(x, na.rm = TRUE, ...) {
# #   qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
# #   H <- 1.5 * IQR(x, na.rm = na.rm)
# #   y <- x
# #   y[x < (qnt[1] - H)] <- NA
# #   y[x > (qnt[2] + H)] <- NA
# #   y
# # }
# 
# 
# #################################################################
# # Now use edgeR package to standardize EpiRAD count data by library size
# 
# #read in the file for edgeR
# counts <- lapply(Epidata4, DGEList)
# counts
# #TMM normalization (corrects for library size)
# counts2 <- lapply(counts, calcNormFactors)
# counts2
# #extract normalized counts
# counts2_cpm <- lapply(counts2, cpm, normalized.lib.sizes=TRUE, log=TRUE)
# 
# 
# #plot normalized counts
# par(mfrow = c(5, 5))
# par(mar = c(2, 2 ,2 ,2), oma = c(4, 4, 0.5, 0.5))
# 
# for (i in seq(1,ncol(counts2_cpm[[1]]), by = 2)){
#   for (j in 1:8){
#     plot(counts2_cpm[[j]][,i], counts2_cpm[[j]][,i+1], main = colnames(counts2_cpm[[j]][i]), col = "blue")
#   }
# }
# 
# 
# ##################################################################
# #Using lm to get residuals
# 
# #models <- list()  NOT WORKING
# #for (i in seq(1,4, by = 2)){
# #  for (j in 1:4){
# #    models <- lapply(counts2_cpm, lm, formula = counts2_cpm[[j]][,i] ~ counts2_cpm[[j]][,i+1])
# #  }
# #}
# 
# 
# model_10_15 <- lm(counts2_cpm[[1]][,2] ~ counts2_cpm[[1]][,1])
# model_10_16 <- lm(counts2_cpm[[1]][,4] ~ counts2_cpm[[1]][,3])
# model_11_15 <- lm(counts2_cpm[[2]][,2] ~ counts2_cpm[[2]][,1])
# model_11_16 <- lm(counts2_cpm[[2]][,4] ~ counts2_cpm[[2]][,3])
# model_02_15 <- lm(counts2_cpm[[3]][,2] ~ counts2_cpm[[3]][,1])
# model_02_16 <- lm(counts2_cpm[[3]][,4] ~ counts2_cpm[[3]][,3])
# model_03_15 <- lm(counts2_cpm[[4]][,2] ~ counts2_cpm[[4]][,1])
# model_03_16 <- lm(counts2_cpm[[4]][,4] ~ counts2_cpm[[4]][,3])
# model_05_15 <- lm(counts2_cpm[[5]][,2] ~ counts2_cpm[[5]][,1])
# model_05_16 <- lm(counts2_cpm[[5]][,4] ~ counts2_cpm[[5]][,3])
# model_06_15 <- lm(counts2_cpm[[6]][,2] ~ counts2_cpm[[6]][,1])
# model_06_16 <- lm(counts2_cpm[[6]][,4] ~ counts2_cpm[[6]][,3])
# model_08_15 <- lm(counts2_cpm[[7]][,2] ~ counts2_cpm[[7]][,1])
# model_08_16 <- lm(counts2_cpm[[7]][,4] ~ counts2_cpm[[7]][,3])
# model_09_15 <- lm(counts2_cpm[[8]][,2] ~ counts2_cpm[[8]][,1])
# model_09_16 <- lm(counts2_cpm[[8]][,4] ~ counts2_cpm[[8]][,3])
# 
# 
# models <- list(model_10_15, model_10_16, model_11_15, model_11_16, model_02_15, 
#                model_02_16, model_03_15, model_03_16, model_05_15, model_05_16, 
#                model_06_15, model_06_16, model_08_15, model_08_16, model_09_15, model_09_16)
# 
# 
# residuals <- lapply(models, '[[', 2)
# 
# #plot residuals
# par(mfrow = c(5, 5))
# par(mar = c(2, 2 ,2 ,2), oma = c(4, 4, 0.5, 0.5))
# 
# for (i in 1:length(residuals)){
#   plot(residuals[[i]], col = "blue", ylim = c(-10, 4))
# }
# 
# #k-means clustering to separate methylated and unmethylated loci
# set.seed(1234)
# clusters <- lapply(residuals, kmeans, 2)
# clusters2 <- lapply(clusters, '[[', 1)
# clusters3 <- lapply(clusters2[c(3,4,5,6,7,11,12,13,14,15)], function(x) ifelse(x == 2, 1, 2))
# clusters4 <- list(clusters2[[1]], clusters2[[2]], clusters3[[1]], clusters3[[2]], 
#                   clusters3[[3]], clusters3[[4]], clusters3[[5]], clusters2[[8]], 
#                   clusters2[[9]], clusters2[[10]], clusters3[[6]], clusters3[[7]],
#                   clusters3[[8]], clusters3[[9]], clusters3[[10]], clusters2[[16]])
# 
# #plot residuals with k-means colors
# par(mfrow = c(5, 5))
# par(mar = c(2,2, 2, 2), oma = c(4, 4, 0.5, 0.5)) 
# 
# for (i in 1:length(clusters4)){
#   plot(residuals[[i]], col = clusters4[[i]], ylim = c(-10, 4))
# }
# 
# #################################################################
# #Make binary dataset of EpiRAD data based on residuals <=-1
# #All methylated loci converted to 1, nonmethylated to zero
# #Select samples of interest (some have very low sample sizes)
# 
# Sym_resid_all_binary <- lapply(clusters4, function(x) x-1)
# 
# #proportion of methylated cutsites
# Sym_prop_methyl <- lapply(Sym_resid_all_binary, function(x) sum(x) / length(x))
# 
# #Get differences between 2015 and 2016 methylation for each replicate
# #0 is no change, 1 is demethylation, -1 is methylation
# # diff <- list()
# # for (i in seq(1,length(unlist(Sym_resid_all_binary)), by = 2)){
# #   diff[[i]] <- Sym_resid_all_binary[[i]] - Sym_resid_all_binary[[i+1]]
# # }
# # 
# # past10diff <- as.data.frame(unlist(diff[[1]]))
# # past11diff <- as.data.frame(unlist(diff[[3]]))
# # past02diff <- as.data.frame(unlist(diff[[5]]))
# # past03diff <- as.data.frame(unlist(diff[[7]]))
# # past05diff <- as.data.frame(unlist(diff[[9]]))
# # past06diff <- as.data.frame(unlist(diff[[11]]))
# # past08diff <- as.data.frame(unlist(diff[[13]]))
# # past09diff <- as.data.frame(unlist(diff[[15]]))
# # 
# # temp <- merge(past10diff, past11diff, by = 0, all.x = TRUE, all.y = TRUE)
# # temp <- merge(temp, past02diff, by.x = 1, by.y = 0, all.x = TRUE, all.y = TRUE)
# # temp <- merge(temp, past03diff, by.x = 1, by.y = 0, all.x = TRUE, all.y = TRUE)
# # temp <- merge(temp, past05diff, by.x = 1, by.y = 0, all.x = TRUE, all.y = TRUE)
# # temp <- merge(temp, past06diff, by.x = 1, by.y = 0, all.x = TRUE, all.y = TRUE)
# # temp <- merge(temp, past08diff, by.x = 1, by.y = 0, all.x = TRUE, all.y = TRUE)
# # temp <- merge(temp, past09diff, by.x = 1, by.y = 0, all.x = TRUE, all.y = TRUE)
# # row.names(temp) <- temp$Row.names
# # difftable <- temp[,c(2:9)]
# # diff1 <- difftable[rowSums(difftable, na.rm = TRUE) < 8, ]
# # diff2 <- diff1[abs(rowSums(difftable, na.rm = TRUE)) >= 3, ]
# # 
# # #get seq numbers, subtracting one from each number to match `.loci` numbers
# # seqs <- as.numeric(rownames(diff2))-1
# # #write to file
# # write.table(seqs, file = "seqs_test.txt", row.names = FALSE, col.names = FALSE)
# # #then in bash use grep -B 1 -wFf seqs_test.txt data2.loci > diffseqs.txt
# # #awk -v n=3 'NR%n==1' diffseqs.txt | awk '{print $2}' > diffseqs2.txt

##############################

# #apply density function to residuals
# resid_dens <- lapply(residuals, density)
# 
# #extract x,y coordinates from density list
# resid_dens_df <- as.data.frame(lapply(resid_dens, function(x) {x[c(1,2)]}))  
# 
# #plot density curves
# par(mfrow = c(5, 5))
# par(mar = c(2,2, 2, 2), oma = c(4, 4, 0.5, 0.5)) 
# 
# for (i in seq(1,34, by = 2)){
#   plot(resid_dens_df[,i], resid_dens_df[,i+1], main = colnames(resid_dens_df[i]), col = "blue")
# }
# 
# #run mixture model analysis on residuals. mu (mean) was taken from average of model run
# #with unspecified mu
# mixmdls <- list()
# for (i in 1:17){
#   mixmdls[[colnames(resid_all)[i]]] <- normalmixEM(resid_all[,i], mu = c(-8, 1.76))
# }
# 
# #plot density curves from normalmixEM
# par(mfrow = c(5, 5))
# par(mar = c(2,2, 2, 2), oma = c(4, 4, 0.5, 0.5)) 
# for (i in 1:17){
#   plot(mixmdls[[i]], whichplots = 2)
# }
# 
# #get posterior probabilities from normalmixEM
# probs <- as.data.frame(lapply(mixmdls, function(x) {x[6]})) 
# probs1 <- as.matrix(probs[seq(1, ncol(probs), 2)])
# names <- c("pa10-15", "pa10-16", "pa11-15", "pa11-16", "pa2-15", "pa2-16", "pa3-15",
#            "pa3-16", "pa5-15", "pa5-16", "pa5h-16", "pa6-15", "pa6-16", "pa8-15", "pa8-16",
#            "pa9-15", "pa9-16")
# colnames(probs1) <- names
# rownames(probs1) <- rownames(resid_all)
# 
# cols1 <- rbPal(10)[as.numeric(cut(probs1[,1],breaks = 10))]
# cols <- rainbow(length(probs1[,3]))[order(order(probs1[,3]))]
# plot(resid_all[,17], col = clusters4[,17])
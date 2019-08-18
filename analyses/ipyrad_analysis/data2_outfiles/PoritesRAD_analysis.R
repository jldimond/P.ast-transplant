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
library(gridExtra, quietly = TRUE)


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
geno6 <- geno2[,c(37:46,59:76,83:90)] #subset of samples of interest
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
  scale_fill_gradient2(low = "white", high = "blue", limit = c(0,0.20), space = "Lab", 
                       name="SNP 
Mismatches") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9.5, hjust = 1))+
          labs(x= "EpiRADseq samples", y = "ddRADseq samples")+
  coord_fixed()

#mean genotyping error (mean of diagonal) 
mean(diag(epidd_dist3)) #0.012
sd(diag(epidd_dist3)) #0.007

#mean resampling error
melted2 <- melted[which(melted$value < 0.1),]
mean(melted2$value[c(2:3,6:8,10:12,15:16,19:20,23:25,27:29,32:33,36:37,40:41)]) #0.017
sd(melted2$value[c(2:3,6:8,10:12,15:16,19:20,23:25,27:29,32:33,36:37,40:41)]) #0.009

#mean genetic distance for non replicates and resamples
mean(melted$value[melted$value >0.1]) #0.163
sd(melted$value[melted$value >0.1]) #0.014


##############################################################################
## Next we read in a text file derived from the "Read Depth"
## .vcf output from ipyrad (v.0.5.15). Read counts are used for analysis of 
## EpiRADseq data.

#Read in data file. The file "data3-2.txt" was generated from the notebook "VCF_readcounts.ipynb"
Epidata <- read.delim("out.DP.FORMAT", header=TRUE)

#Use aggregate to get means for duplicate records of CHROM (locus ID)
Epidata2 <- aggregate(.~CHROM, data=Epidata, mean)
rownames(Epidata2) <- Epidata2$CHROM

Epidata3 <- Epidata2[,c(39:48,61:78,85:92)]



#Remove ddr rows that have any zeros. The premise here is that zeros 
#in the EpiRAD dataset are informative because they may reflect 
#methylation, but they could also relfect true absence of the locus
#in the library. Here the ddRAD library serves to standarize the EpiRAD
#library. Any zeros in the ddRAD libary are treated as absence of the
#locus, thereby leaving zeros in the EpiRAD library only where the 
#locus was counted in the ddRAD library.

Epidata4 <- Epidata3[apply(Epidata3[c(seq(1, 36, by = 2))],1,
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
clusters3[,c(1:5,7:8,11:13,18)] <-ifelse(clusters3[,c(1:5,7:8,11:13,18)] == 2, 1, 2)
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
mean(prop_methyl) #0.186
sd(prop_methyl) #0.009

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
resid17 <- length(which(rowSums(resid_all_binary) == 17))/(length(resid_all_binary[,1])) 
resid18 <- length(which(rowSums(resid_all_binary) == 18))/(length(resid_all_binary[,1])) 

obs <- seq(0,18,1)
meth_obs <- rbind(resid0, resid1, resid2, resid3, resid4, resid5, resid6, resid7, resid8, resid9,
                  resid10, resid11, resid12, resid13, resid14, resid15, resid16, resid17, resid18)
methobs <- cbind(meth_obs,obs)


#################################################################
#Box plots / data summary

#distance matrix of methylation differences (to get between colony per year differences)
temp <- t(resid_all_binary)
temp2 <- dist.gene(temp, method = "percent", pairwise.deletion = FALSE,
                   variance = FALSE)
dist1 <- as.matrix(temp2)

#differences between colonies in 2015
temp <- t(resid_all_binary[,c(1,3,6,8,10,13,15,17)])
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
temp <- t(resid_all_binary[,c(2,4,7,9,11,14,16,18)])
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
year <- rep(c("2015", "2016"), each = 28)
per_change <- rbind(dist5, dist8)
per_change2 <- as.data.frame(cbind(per_change$value, year))
per_change2$V1 <- as.numeric(as.character(per_change2$V1))

#percent of methylated CpGs with year factor
year2 <- rep(c("2015", "2016"), 8)
prop_methyl2 <- as.data.frame(cbind(prop_methyl[c(1:4,6:11,13:18)],year2))
prop_methyl2$V1 <- as.numeric(as.character(prop_methyl2$V1))

############################################################################
#heatmap of EpiRAD data 

#Matrix for heatmap
heat.mat <- cbind(prop_methyl[c(1:4,6:11,13:18)],t(resid_all_binary[,c(1:4,6:11,13:18)]))
names <- c("pa10-15", "pa10-16", "pa11-15", "pa11-16", "pa2-15", "pa2-16", "pa3-15",
           "pa3-16", "pa5-15", "pa5-16", "pa6-15", "pa6-16", "pa8-15", "pa8-16",
           "pa9-15", "pa9-16")
heat.mat2 <- t(heat.mat)
colnames(heat.mat2) <- names
heat.mat2 <- heat.mat2[,sort(colnames(heat.mat2))]


heat <- superheat(heat.mat2[2:630, c(1:16)], bottom.label.text.size = 5,
          scale = FALSE, row.dendrogram = TRUE, col.dendrogram = FALSE, 
          pretty.order.cols = FALSE, pretty.order.rows = TRUE,
          heat.pal = c("#6baed6", "#08519c"), yt = heat.mat2[1,c(1:16)],
          legend = FALSE, grid.vline.col = "white", left.label = "none",
          grid.hline.size = 0.1,bottom.label.text.angle = 90,
          yt.axis.name = "Methylation", yt.plot.type = "bar")

par(mfrow = c(2, 1))
par(mar = c(4, 4.5, 2, 1), oma = c(1, 1, 0, 0))


###########################################################################
#MDS of EpiRAD data

names <- c("pa10-15", "pa10-16", "pa11-15", "pa11-16", "pa11h-16","pa2-15", "pa2-16", "pa3-15",
           "pa3-16", "pa5-15", "pa5-16", "pa5h-16", "pa6-15", "pa6-16", "pa8-15", "pa8-16",
           "pa9-15", "pa9-16")
colnames(resid_all_binary) <- names
resid_t_binary <- t(resid_all_binary)
# euclidean distances between the rows
epidist <- dist(resid_t_binary) 
epifit <- cmdscale(epidist,eig=TRUE, k=2)
epix <- epifit$points[,1]
epiy <- epifit$points[,2]

#MDS plot
p <- plot(epix, epiy, xlab="Coordinate 1", ylab="Coordinate 2", type = 'n')
for (i in seq(1,4, by = 2)){
  arrows(epix[i], epiy[i], epix[i+1], epiy[i+1], length = 0, col = "dark gray")
}
for (i in seq(6,11, by = 2)){
  arrows(epix[i], epiy[i], epix[i+1], epiy[i+1], length = 0, col = "dark gray")
}
for (i in seq(13,18, by = 2)){
  arrows(epix[i], epiy[i], epix[i+1], epiy[i+1], length = 0, col = "dark gray")
}
arrows(epix[3], epiy[3], epix[5], epiy[5], length = 0, col = "dark gray")
arrows(epix[10], epiy[10], epix[12], epiy[12], length = 0, col = "dark gray")

names2 <- c(1,1,2,2,2,3,3,4,4,5,5,5,6,6,7,7,8,8)
shapes <- c(1,2,1,2,3,1,2,1,2,1,2,3,1,2,1,2,1,2)
palette(brewer.pal(n = 8, name = "Set2"))
points(epix, epiy, col = names2, pch = shapes, cex = 1.5, lwd = 2)
legend(1.4,-1.5, legend = c("2015", "2016", "Control"), pch = (1:3), cex = 1.0, pt.lwd = 1.5)
mtext("A", side = 3, adj = 0)

#######################################################
#summary plots and analyses

#Methylated CpGs by year
#test for homogeneity of variances differences between years
bartlett.test(V1 ~ year2, data = prop_methyl2) #variances sig. different p-value = 0.4579

#paired t-test
t.test(V1 ~ year2, data = prop_methyl2, var.equal = FALSE, paired = TRUE) #p-value = 0.5111

#Pairwise methylation difference
#test for homogeneity of variances differences between years
bartlett.test(V1 ~ year, data = per_change2) #variances sig. different p-value = 0.4471

#paired t-test
t.test(V1 ~ year, data = per_change2, var.equal = FALSE, paired = TRUE) #p-value = 0.0005982
mean(per_change2[c(1:28),1]) #2015 0.03957529
sd(per_change2[c(1:28),1]) #2015 0.005906561
mean(per_change2[c(29:56),1]) #2016 0.03577107
sd(per_change2[c(29:56),1]) #2016 0.006848356

#Percent CpG methylation by year
p1 <- ggplot(prop_methyl2, aes(x=year2, y=V1*100, group=year2, fill=year2)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#e5f5f9", "#99d8c9")) +
  labs(x ="Year", y = "Methylated CpGs (%)", title = "B") +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none") +
  annotate(geom="text", x=1.6, y=20.1, label="p = 0.511",
           color="black") 

#Percent pairwise difference between colonies
p2 <- ggplot(per_change2, aes(x=year, y=V1*100, group=year, fill=year)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#e5f5f9", "#99d8c9")) +
  labs(x ="Year", y = "Pairwise methylation difference (%)", title = "C") +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.position = "none") +
  annotate(geom="text", x=1.6, y=4.88, label="p < 0.001",
           color="black")

p3 <- grid.arrange(p1, p2, ncol=2)
plot(p3)

#################################################################
# Venn diagram looking at shared loci across individuals

#just get diff methylated loci 
diffmeth <- resid_all_binary[,c(1:4,6:11,13:18)][which(rowSums(resid_all_binary[,c(1:4,6:11,13:18)]) >= 1 ),]
diffmeth2 <- diffmeth[which(rowSums(diffmeth) < 16 ),]
first <- seq(1,ncol(diffmeth2),2)  
second <- seq(2,ncol(diffmeth2),2)  
diffmeth3 <- abs(diffmeth2[, first] - diffmeth2[, second])
#differentially methylated loci
diffmeth4 <- which(rowSums(diffmeth3) >=2)
#need locus names to match .loci file (subtract 1 from name)
diffmeth5 <- as.numeric(names(diffmeth4))-1
#table of locus IDs
write.table(diffmeth5, "diffmeth5.txt", sep="\t", row.names=F)



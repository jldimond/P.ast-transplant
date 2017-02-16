## ---
## output: github_document
## ---

############################################################################
## This script analyzes ddRADseq and EpiRADseq data in branching Porites spp.


#Set directory
setwd("~/Documents/Projects/PoritesRADseq/Branching-Porites/analyses/ipyrad_analysis/data3_outfiles")


## First we read in a .str file from the ipyrad output and subset the
## data to get only ddRAD data and extract only SNPs without missing data.


# Read in data file
ustr <- read.delim("data3.u.str", header=FALSE)
ddata <- ustr[,colSums(is.na(ustr))<nrow(ustr)]
#Vector of  unlinkedSNP IDs corresponding to data3.snps.map
unlinkedsnps <- read.delim("loc_id4.txt", header=FALSE)
unlinkedsnps2 <- t(unlinkedsnps)
ddata2 <- ddata[,2:11823]
colnames(ddata2) <- unlinkedsnps2
ddata3 <- cbind(ddata[,1],ddata2)
ddata4 <- t(ddata3)
ddata5 <- as.data.frame(ddata4)
#Extract only ddRAD data
ddata6 <- ddata5[,c(1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,42,47,48,51,52,55,56,59,60,63,64,69,70,73,74,77,78,81,82,85,86,89,90,93,94,97,98,101,102,105,106,109,110)]
#Remove any SNPs with missing data (-9 is the NA value)
ddata7 <- ddata6[!rowSums(ddata6 == -9) >= 1,]
ddata8 <-t(ddata7)
rownames(ddata8) <- ddata8[,1]
ddata9 <- ddata8[,2:1114]
write.table(ddata9, file = "data3-2.str", row.names = TRUE, col.names = TRUE, quote = FALSE)

## Next we read in a .geno file from the ipyrad output and extract 
## ddRAD data without missing values. This is the best filetype for MDS.

geno1 <- read.table("data3.u.geno", colClasses = 'character', header = FALSE)
geno2 <- read.fwf("data3.u.geno", widths=rep(1, max(nchar(geno1$V1))), colClasses = 'numeric', header=FALSE)
header <- read.delim("header_data3.txt", header=FALSE)
names <- t(header)
names2 <-as.vector(names)
colnames(geno2) <- names2

#Select samples of interest (some have very low sample sizes)

geno3 <- geno2[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56)]

#Matrix with only ddr loci **if including sample 101

geno4 <- geno3[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52)]

#Get rid of rows with any NAs (9)

geno5 <- geno4[!rowSums(geno4 == 9) >= 1,]
head(geno5)


## Next we read in a text file derived from the "Base Counts"
## .vcf output from ipyrad. Base counts are used for analysis of 
## EpiRADseq data.


#Read in data file
Epidata <- read.delim("data3-2.txt", header=FALSE)
#Since the base counts were split into four columns for each base, these need 
#to be summed
Epidata2 <- t(sapply(seq(4,ncol(Epidata), by=4), function(i) {
  indx <- i:(i+3)
  rowSums(Epidata[indx[indx <= ncol(Epidata)]])}))

#The resulting file needs to be transposed and turned into a dataframe
Epidata3 <- as.data.frame(t(Epidata2))
#Add column with locus number (CHROM from .vcf file)
locus <- Epidata[,1]
row.names(Epidata3) <- locus
#Add header names
header <- read.delim("header_data3.txt", header=FALSE)
names <- as.vector(t(header))
colnames(Epidata3) <- names
#Select samples of interest (some have very low sample sizes)
Epidata4 <- Epidata3[,c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                  24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,
                  41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56)]

#Remove ddr rows that have any zeros. The premise here is that zeros 
#in the EpiRAD dataset are informative because they may reflect 
#methylation, but they could also relfect true absence of the locus
#in the library. Here the ddRAD library serves to standarize the EpiRAD
#library. Any zeros in the ddRAD libary are treated as absence of the
#locus, thereby leaving zeros in the EpiRAD library only where the 
#locus was counted in the ddRAD library.

Epidata5 <- Epidata4[apply(Epidata4[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,
                             31,33,35,37,39,41,43,45,47,49)],1,
                     function(z) !any(z==0)),] 


#################################################################
# Now use edgeR package to standardize EpiRAD count data by library size

library("edgeR")
#read in the file to edgeR
counts <- DGEList(counts=Epidata5)
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
  plot(Epidata5[,i], Epidata5[,i+1], main = colnames(Epidata5[i]), col = "blue")
}

#plot normalized counts
par(mfrow = c(5, 5))
par(mar = c(2, 2, 2, 2), oma = c(4, 4, 0.5, 0.5)) 

for (i in seq(1,49, by = 2)){
  plot(counts2_cpm[,i], counts2_cpm[,i+1], main = colnames(counts2_cpm[i]), col = "blue")
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
  plot(resid_all[,i], col = "blue")
}

#Plot to compare raw data to residuals
par(mfrow = c(2, 1))
par(mar = c(4, 4.5, 2, 1), oma = c(1, 1, 0, 0))
plot(Epidata5[,13], Epidata5[,14], xlab = "ddRAD read counts", ylab = "EpiRAD read counts", col = "blue")
plot(resid_all[,7], ylab = "Residual", col = "blue")
mtext('A', side=3, line=-1.6, at = 0.2, outer=TRUE)
mtext('B', side=3, line=-16.8, at = 0.2, outer=TRUE)

#################################################################
#Make binary dataset of EpiRAD data based on residuals <=-1

#All methylated loci converted to 1, nonmethylated to zero
resid_all_binary <- ifelse(resid_all<=-1, 1, 0)

#proportion of methylated cutsites
prop_methyl <- colSums(resid_all_binary) / nrow(resid_all_binary)
barplot(prop_methyl)
dens <- density(prop_methyl)
plot(dens)

#Get only rows that are differentially methylated
resid1 <- resid_all_binary[rowSums(resid_all_binary) < 25, ]
resid2 <- resid1[rowSums(resid1) >= 1, ]

########################################################
#Read in sample info (sample #, depth, symbiont type, diameter)
sinfo <- read.table("sample_info.txt", colClasses = 'character', header = FALSE)
#transpose
tsinfo <- t(sinfo)
#create vectors for diameter (note whether sample 101
#was included or not)
diam <- tsinfo[4,]

###########################################################################
#MDS of ddRAD and EpiRAD data
#First ddRAD

geno6 <- t(geno5)
ddist <- dist(geno6) # euclidean distances between the rows
ddfit <- cmdscale(ddist,eig=TRUE, k=2)
ddx <- ddfit$points[,1]
ddy <- ddfit$points[,2]

#Now EpiRAD

resid_t_binary <- t(resid_all_binary)
epidist <- dist(resid_t_binary) # euclidean distances between the rows
epifit <- cmdscale(epidist,eig=TRUE, k=2)
epix <- epifit$points[,1]
epiy <- epifit$points[,2]

#Plot both MDS plots

par(mfrow = c(2, 1))
par(mar = c(4, 4.5, 2, 1), oma = c(1, 1, 0, 0))
plot(ddx, ddy, xlab="Coordinate 1", ylab="Coordinate 2", col = "blue")
plot(epix, epiy, xlab="Coordinate 1", ylab="Coordinate 2", col = "blue")
mtext('A', side=3, line=-1.6, at = 0.22, outer=TRUE)
mtext('B', side=3, line=-20, at = 0.22, outer=TRUE)

####################################################################
#DAPC (discriminant analysis of principal components) of SNPs using adegenet

library("adegenet")
library("ade4")

#Read in unlinked SNP file created at the top of this script.
#Note: must manually delete 1st entry in 1st row of "data3-2.str"
genind1 <- read.structure("data3-2.str", n.ind = 27, n.loc = 1113, 
                          onerowperind = FALSE, col.lab = 1, 
                          NA.char = "-9", ask = FALSE, 
                          row.marknames = 1, quiet = FALSE)

#Find optimal number of clusters irrespective of species id
#In this case, best to retain all PCs
groups <- find.clusters(genind1, max.n.clust=10, n.pca = 24,
                        choose.n.clust = FALSE, criterion = "min")


#Cross validation to determine number of PCs to retain
xval <- xvalDapc(genind1@tab, groups$grp, n.pca.max = 25, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)
dev.off() 
#show max number of PCs to retain
xval[2:6]

#perform dapc using groups defined above group (groups$grp). Note 
#that n.pca and n.da can be left blank and the program will query 
#which values to choose.
dapc1 <- dapc(genind1, pop = groups$grp, n.pca=9, n.da = 2)
scatter(dapc1, posi.da = "bottomleft", scree.pca = TRUE, posi.pca = "bottomright")

#Plot cluster vs. BIC and DAPC together
par(mfrow = c(2, 1))
par(mar = c(6, 4.5, 2, 2), oma = c(1, 1, 1, 1))
plot(groups$Kstat, xlab = "Groups (K)", ylab = "BIC", pch = 16, 
     xaxp = c(1, 9, 4))
lines(groups$Kstat)
cols <- c("red", "orange", "purple")
scatter(dapc1, #label.inds = list(air = 0.1, pch = 0.5),
        posi.da = "bottomleft", scree.pca = TRUE, posi.pca = "bottomright",
        cell=0, cstar=0, clab=0, cex=3, solid=.4, bg="white", 
        leg=TRUE, posi.leg="topleft", col = cols)

#look at loadings of individual loci
set.seed(4)
contrib <- loadingplot(dapc1$var.contr, axis=1, 
                       threshold= quantile(dapc1$var.contr,0.90), lab.jitter=1)

###########################################################
#Fst
library("hierfstat")

fst <- pairwise.fst(genind1, pop = groups$grp)
fst

####################################################################
#DAPC (discriminant analysis of principal components) of Epi-loci using adegenet

#This uses the binary methylation file generated above (must be transposed)
resid_t_binary <- t(resid_all_binary)

#Find optimal number of clusters irrespective of species id
#In this case, best to retain all PCs
Epigroups <- find.clusters(resid_t_binary, max.n.clust=10, n.pca = 24,
                           choose.n.clust = TRUE, criterion = "min")

Epidapc1 <- dapc(resid_t_binary, grp = Epigroups$grp, n.pca=7, n.da = 2)
scatter(Epidapc1)

###########################################################
#boxplot comparing branch diameter among groups from SNP dapc

diam2 <- as.numeric(sinfo[,4])
boxplot(diam2 ~ dapc1$assign, xlab = "Group", ylab =  "Diameter (mm)")
text(1.376804, 22.84152, "a", cex = 1)
text(1.990675, 13.63027, "b", cex = 1)
text(2.985303, 20.00729, "b", cex = 1)

#ANOVA
#Run lm on diameter by group
model <- lm(diam2 ~ dapc1$assign)
#Check model (qq plot, etc)
par(mfrow=c(2,2))
plot(model)
library(lmtest)
#Breush Pagan Test for heteroscadisticity
bpt <- bptest(model)
print(bpt)
#Run ANOVA
aov <- anova(model)
print(summary(aov))
#pairwise t test with bonferonni adjustment
ttest <- pairwise.t.test(diam2, dapc1$assign, p.adj = "bonf")
print(ttest)

######################################################################
#Multiple regression: test linear model of depth, symbiont type, and
#branch diameter vs. dapc of SNPs using a single DA from the model above

#dapc with single DF
dapc1 <- dapc(genind1, pop = groups$grp, n.pca=9, n.da=1)
#Select # pcs = 9, # df = 1
scatter(dapc1)
#convert DF coord to vector
dapc1_da1 <- dapc1$ind.coord
#replace unknown symbiont type NAs to "U"
sinfo[is.na(sinfo)] <- "U"
#fit model with water depth (sinfo$V2), branch diameter (sinfo$V5), 
# and symbiont type (sinfo$V4) 
fit <- lm(na.omit(dapc1_da1 ~ as.numeric(sinfo$V2) + as.factor(sinfo$V4) + as.numeric(sinfo$V5)))
#Relative importance of different variables in model
library(relaimpo)
relimp <- calc.relimp(fit,type=c("lmg","last","first"),
                      rela=TRUE)
print(relimp)

# Bootstrap Measures of Relative Importance (1000 samples) 
boot <- boot.relimp(fit, b = 1000, type = c("lmg", "last", "first"), 
                    rank = TRUE, diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result

############################################################################
#heatmap of EpiRAD data 

#transpose differentially expressed dataset
resid_t_diff <- t(resid2)

#color pallet matrix for SNP groups and branch diam groups (in between whitespace)
groupvec <- as.character(groups$grp)
SNP_Groups <- groupvec[c(2:10,12:27)]
SNP_Groups <- replace(SNP_Groups, which(SNP_Groups == 1), "red")
SNP_Groups <- replace(SNP_Groups, which(SNP_Groups == 2), "orange")
SNP_Groups <- replace(SNP_Groups, which(SNP_Groups == 3), "purple")
col_pal = colorRampPalette(c('light gray', 'black'))(25+1)
Diameter = col_pal[ cut(as.numeric(sinfo2$V5), data_seq, include.lowest=T) ]
white <- colorRampPalette(colors= "#ffffff")
whitespace <- white(25)
myCols = cbind(SNP_Groups, whitespace, Diameter)

library("heatmap.plus")

heatmap.plus(resid_t_diff, scale = "none", labRow = sinfo2$V1, labCol = FALSE,
        RowSideColors = myCols, col = c("#6baed6", "#08519c"))

############################################################################
#Determine loci that are both weighting the SNP grouping and 
#are differentially methylated

loci.90 <- unique(gsub("\\..*","",contrib$var.names))
loci.Epi <- colnames(resid_t_diff)
loci.both <- intersect(loci.90,loci.Epi)

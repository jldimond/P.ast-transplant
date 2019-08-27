# This script analyzes the following AQUA MODIS satellite data:
# AQUA MODIS SST 11u nighttime 4km monthly 
# AQUA MODIS Chlorophyll concentration, OCx algorithm 4km monthly
# AQUA MODIS Diffuse attenuation coefficient at 490 nm, KD2 algorithm 4km monthly
# Downloaded from https://oceancolor.gsfc.nasa.gov/cgi/l3 and accessed 8/24/2018
# Example citation for diff dataset:
# NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group; (2014): Sea-viewing Wide Field-of-view Sensor (SeaWiFS) Ocean Color Data, NASA OB.DAAC. http://doi.org/10.5067/ORBVIEW-2/SEAWIFS_OC.2014.0. Accessed on 2016/02/29.
# Monthly data, range = Nov 2015 - Oct 2016
# Belize basemap from http://www.biodiversity.bz/Belize_Basemap.zip
# UN Environment coral reef atlas from http://datadownload.web-production.linode.unep-wcmc.org/download/qesFJJYG/WCMC008_CoralReefs2010_v3.zip

library(raster)
library(ncdf4)
library(sp)
library(rgdal)
library(dplyr)
library(plyr)
library(reshape2)
library(ggplot2)
library(cowplot)

setwd("/Users/jd/Documents/Projects/PoritesRADseq/P.ast-transplant/analyses/satellite-data-analysis/")

# list files
sst = list.files("./sst", pattern="sst_4km.nc", full.names=TRUE) 

# create ratser stack
sst_stack <- stack(sst)

# crop stack to new extent
sst_stack_crop <- crop(sst_stack, extent(-88.4,-87.9,16.55,17.125))

# get mean and range of sst over 13 months of data
sst_avg <- mean(sst_stack_crop)
sst_range <- max(sst_stack_crop)-min(sst_stack_crop)

#now for chla data
# list files
chla = list.files("./chla", pattern="chl_ocx_4km.nc", full.names=TRUE) 

# create raster stack
chla_stack <- stack(chla)

# crop stack to new extent
chla_stack_crop <- crop(chla_stack, extent(-88.4,-87.9,16.55,17.125))

# get mean chla over 13 months of data
chla_avg <- mean(chla_stack_crop, na.rm = TRUE)

# now for the attenuation data
# list files
atten = list.files("./attenuation", pattern="Kd_490_4km.nc", full.names=TRUE) 

# create ratser stack
atten_stack <- stack(atten)

# crop stack to new extent
atten_stack_crop <- crop(atten_stack, extent(-88.4,-87.9,16.55,17.125))

# get mean chla over 13 months of data
atten_avg <- mean(atten_stack_crop, na.rm = TRUE)

# open Belize basemap
belize <- readOGR(dsn = "./Belize_Basemap", layer = "Belize_Basemap")

# transform projection to the projection of the MODIS data
# "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
belize_trans <- spTransform(belize, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
belize_trans_clip <- crop(belize_trans, extent(-88.42,-87.9,16.54,17.125))

plot(belize_trans, col = "gray")

# open coral reef basemap
coralreef <- readOGR(dsn = "./14_001_WCMC008_CoralReefs2010_v3/01_Data", layer = "WCMC008_CoralReef2010_Py_v3")

# transform projection to the projection of the MODIS data
# "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
coralreef_trans <- spTransform(coralreef, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coralreef_trans_clip <- crop(coralreef_trans, extent(-88.42,-87.9,16.54,17.125))

plot(coralreef_trans_clip, col = "gray")

# get coordinates of coral specimens 
corals <- read.delim("./specimens.txt", header = FALSE)
names(corals) <- c("Latitude","Longitude","IDnumber")

# extract sst data from rasters using coral site data

corals2 <- cbind(corals[,2], corals[,1], corals[,3])
sst_corals <- extract(sst_stack_crop, corals2[,1:2])
sst_corals <- as.data.frame(sst_corals)
rownames(sst_corals) <- corals$IDnumber
sst_corals <- cbind(rep(c("Offshore", "Offshore", "Inshore", "Inshore"), 2), sst_corals)
colnames(sst_corals) <- c("Location","11","12","1","2","3","4","5","6","7","8","9","10")

melted <- melt(sst_corals, id.vars= "Location")

sst_summary <- ddply(melted, "Location", summarise,
                     Mean = mean(value), sd = sd(value),
                     sem = sd(value)/sqrt(length(value)),
                     Range = max(value)-min(value), na.rm = TRUE)

sst_mean <- apply(sst_corals[,2:13], 1, mean)
sst_mean2 <- cbind(sst_corals$Location, sst_mean)

#test for homogeneity of variances differences between locations
bartlett.test(sst_mean ~ V1, data = sst_mean2) #variances not sig. different

#paired t-test
t.test(sst_mean ~ V1, data = sst_mean2, var.equal = TRUE, paired = FALSE) #df = 6, p-value = 0.0018

sst_range <- apply(sst_corals[,2:13], 1, function(x) max(x)-min(x))
sst_range2 <- cbind(sst_corals$Location, sst_range)

# extract chla data from rasters using coral site data

chla_corals <- extract(chla_stack_crop, corals2[,1:2])
rownames(chla_corals) <- corals$IDnumber
offshore_chl <- chla_corals[c(1:2,5:6),]
inshore_chl <- chla_corals[c(3:4,7:8),]
offshore_chl_mean <- mean(offshore_chl, na.rm=TRUE)
offshore_chl_se <- sd(offshore_chl, na.rm=TRUE)/sqrt(44)
inshore_chl_mean <- mean(inshore_chl, na.rm=TRUE)
inshore_chl_se <- sd(inshore_chl, na.rm=TRUE)/sqrt(42)
chla_summary1 <- as.numeric(rbind(inshore_chl_mean,offshore_chl_mean))
chla_summary2 <- as.numeric(rbind(inshore_chl_se,offshore_chl_se))
chla_summary <- cbind(chla_summary1, chla_summary2, c("Inshore", "Offshore"))
colnames(chla_summary) <- c("Mean", "SEM", "Location")
chla_summary <- as.data.frame(chla_summary)
chla_summary$Mean <- as.numeric(as.character(chla_summary$Mean))
chla_summary$SEM <- as.numeric(as.character(chla_summary$SEM))

inshore_chl2 <- apply(inshore_chl, 1, mean, na.rm = TRUE)
offshore_chl2 <- apply(offshore_chl, 1, mean, na.rm = TRUE)
chl_mean2 <- as.data.frame(cbind(c(rep("Inshore",4),rep("Offshore",4)),c(inshore_chl2, offshore_chl2)))
chl_mean2$V2 <- as.numeric(as.character(chl_mean2$V2))


# extract atten data from rasters using coral site data

atten_corals <- extract(atten_stack_crop, corals2[,1:2])
rownames(atten_corals) <- corals$IDnumber
offshore_atten <- atten_corals[c(1:2,5:6),]
inshore_atten <- atten_corals[c(3:4,7:8),]
offshore_atten_mean <- mean(offshore_atten, na.rm=TRUE)
offshore_atten_se <- sd(offshore_atten, na.rm=TRUE)/sqrt(44)
inshore_atten_mean <- mean(inshore_atten, na.rm=TRUE)
inshore_atten_se <- sd(inshore_atten, na.rm=TRUE)/sqrt(42)
atten_summary1 <- as.numeric(rbind(inshore_atten_mean,offshore_atten_mean))
atten_summary2 <- as.numeric(rbind(inshore_atten_se,offshore_atten_se))
atten_summary <- cbind(atten_summary1, atten_summary2, c("Inshore", "Offshore"))
colnames(atten_summary) <- c("Mean", "SEM", "Location")
atten_summary <- as.data.frame(atten_summary)
atten_summary$Mean <- as.numeric(as.character(atten_summary$Mean))
atten_summary$SEM <- as.numeric(as.character(atten_summary$SEM))

inshore_atten2 <- apply(inshore_atten, 1, mean, na.rm = TRUE)
offshore_atten2 <- apply(offshore_atten, 1, mean, na.rm = TRUE)
atten_mean2 <- as.data.frame(cbind(c(rep("Inshore",4),rep("Offshore",4)),c(inshore_atten2, offshore_atten2)))
atten_mean2$V2 <- as.numeric(as.character(atten_mean2$V2))

#dataframe with all means
merge_sat <- merge(sst_mean2, sst_range2, by="row.names")
merge_sat <- merge(merge_sat, atten_mean2, by.x="Row.names", by.y="row.names")
merge_sat <- merge(merge_sat, chl_mean2, by.x="Row.names", by.y="row.names")
merge_sat <- merge_sat[,c(1,3,5,7,9)]
names(merge_sat) <- c("Colony", "SST_Mean", "SST_Range", "kd490", "Chl a")
merge_sat2 <- cbind(rep(1, 8), merge_sat)
names(merge_sat2) <- c("i", "Colony", "SST_Mean", "SST_Range", "kd490", "Chl_a")
merge_sat2$i <- as.factor(merge_sat2$i)
#add row for common garden
cg <-data.frame("1",1,28.40292, 3.395, 0.04608333, 0.2462378)
names(cg)<-colnames(merge_sat2)
merge_sat2 <- rbind(cg, merge_sat2)
CG2 <- subset(merge_sat2, Colony == "CG")

plot1 <-ggplot(merge_sat2, aes(x=i, y=SST_Mean)) +
  geom_dotplot(binaxis='y', stackdir='center') +
  labs(x="", y= "SST Mean (°C)") +
  geom_point(data=merge_sat2[1, ], aes(x=i, y=SST_Mean), colour="red") +
  theme_gray() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 

plot2 <-ggplot(merge_sat2, aes(x=i, y=SST_Range)) +
  geom_dotplot(binaxis='y', stackdir='center') +
  labs(x="", y= "SST Range (°C)") +
  geom_point(data=merge_sat2[1, ], aes(x=i, y=SST_Range), colour="red") +
  theme_gray() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 

plot3 <-ggplot(merge_sat2, aes(x=i, y=kd490)) +
  geom_dotplot(binaxis='y', stackdir='center') +
  ylab(bquote('kd490 ('*m^-1*')')) +
  xlab("") +
  geom_point(data=merge_sat2[1, ], aes(x=i, y=kd490), colour="red") +
  theme_gray() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 
  
plot4 <-ggplot(merge_sat2, aes(x=i, y=Chl_a)) +
  geom_dotplot(binaxis='y', stackdir='center') +
  ylab(bquote('Chl a ('*mg~ m^-3*')')) +
  xlab("") +
  geom_point(data=merge_sat2[1, ], aes(x=i, y=Chl_a), colour="red") +
  theme_gray() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 
  
plot5 <- plot_grid(plot1, plot2, plot4, nrow=1)


# turn rasters into dataframes in order to plot with ggplot
sst_avg_df <- data.frame(rasterToPoints(sst_avg)) 
names(sst_avg_df) <- c("Longitude","Latitude","SST_Mean")

sst_range_df <- data.frame(rasterToPoints(sst_range)) 
names(sst_range_df) <- c("Longitude","Latitude","SST_Range")

chla_avg_df <- data.frame(rasterToPoints(chla_avg)) 
names(chla_avg_df) <- c("Longitude","Latitude","Chl_a")

atten_avg_df <- data.frame(rasterToPoints(atten_avg)) 
names(atten_avg_df) <- c("Longitude","Latitude","kd490")

p <- ggplot(sst_avg_df)
p <- p + geom_tile(aes(x=Longitude,y=Latitude,fill=SST_Mean)) 
p <- p + scale_fill_gradientn(colors=c('#313695','#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43','#d73027','#a50026'), 
                              name="SST Mean (°C)")
p <- p + coord_equal() + xlab("Longitude") + ylab("Latitude")
p <- p + geom_polygon(data=belize_trans_clip, aes(x=long, y=lat, group=group))
p <- p + geom_polygon(data=coralreef_trans_clip, aes(x=long, y=lat, group=group), color = "gray", fill= "gray")
p <- p + geom_text(aes(x=Longitude,y=Latitude, label=IDnumber),data=corals)
p <- p + theme_minimal()
p <- p + annotate("segment", x = -88, xend = -88.08, y = 16.8, yend = 16.8, colour = "black", size=1, arrow=arrow(length=unit(0.30,"cm"), type = "closed"))

p2 <- ggplot(sst_summary, aes(x=Location, y=Mean, color=Location, group = Location))
p2 <- p2 + geom_errorbar(aes(ymin=Mean-sem, ymax=Mean+sem), width=.1, color = c("#525252", "#969696"))
p2 <- p2 + geom_col(fill = c("#525252", "#969696"), color= c("#525252", "#969696"))
p2 <- p2 + theme_bw(base_size = 9)
p2 <- p2 + theme(axis.line = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x=element_blank())
p2 <- p2 + theme(aspect.ratio=1)
p2 <- p2 + ylab("SST (°C)")
p2 <- p2 + theme(legend.position="none")

p3 <- ggplot(sst_range_df)
p3 <- p3 + geom_tile(aes(x=Longitude,y=Latitude,fill=SST_Range)) 
p3 <- p3 + scale_fill_gradientn(colors=c('#313695','#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43','#d73027','#a50026'), 
                                name="SST Range (°C)")
p3 <- p3 + coord_equal() + xlab("Longitude") + ylab("Latitude")
p3 <- p3 + geom_polygon(data=belize_trans_clip, aes(x=long, y=lat, group=group))
p3 <- p3 + geom_polygon(data=coralreef_trans_clip, aes(x=long, y=lat, group=group), color = "gray", fill= "gray")
p3 <- p3 + geom_text(aes(x=Longitude,y=Latitude, label=IDnumber),data=corals)
p3 <- p3 + theme_minimal()
p3 <- p3 + annotate("segment", x = -88, xend = -88.08, y = 16.8, yend = 16.8, colour = "black", size=1, arrow=arrow(length=unit(0.30,"cm"), type = "closed"))


p4 <- ggplot(sst_summary, aes(x=Location, y=Range, colour=Location, group = Location))
p4 <- p4 + geom_col(fill = c("#525252", "#969696"), color= c("#525252", "#969696"))
p4 <- p4 + theme_bw(base_size = 9)
p4 <- p4 + theme(axis.line = element_line(color = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 axis.title.x=element_blank())
p4 <- p4 + theme(aspect.ratio=1)
p4 <- p4 + ylab("SST Range (°C)")
p4 <- p4 + theme(legend.position="none")

p5 <- ggplot(chla_avg_df)
p5 <- p5 + geom_tile(aes(x=Longitude,y=Latitude,fill=Chl_a)) 
p5 <- p5 + scale_fill_gradientn(colors=c('#f7fcfd','#e5f5f9','#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45', '#006d2c','#00441b'), 
                                name=(bquote('Chl a ('*mg~ m^-3*')')))
p5 <- p5 + coord_equal() + xlab("Longitude") + ylab("Latitude")
p5 <- p5 + geom_polygon(data=belize_trans_clip, aes(x=long, y=lat, group=group))
p5 <- p5 + geom_polygon(data=coralreef_trans_clip, aes(x=long, y=lat, group=group), color = "gray", fill= "gray")
p5 <- p5 + geom_text(aes(x=Longitude,y=Latitude, label=IDnumber),data=corals)
p5 <- p5 + theme_minimal() 
p5 <- p5 + annotate("segment", x = -88, xend = -88.08, y = 16.8, yend = 16.8, colour = "black", size=1, arrow=arrow(length=unit(0.30,"cm"), type = "closed"))


p6 <- ggplot(chla_summary, aes(x=Location, y=Mean, colour=Location, group = Location))
p6 <- p6 + geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.1, color = c("#525252", "#969696"))
p6 <- p6 + geom_col(fill = c("#525252", "#969696"), color= c("#525252", "#969696"))
p6 <- p6 + theme_bw(base_size = 9)
p6 <- p6 + theme(axis.line = element_line(color = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 axis.title.x=element_blank())
p6 <- p6 + theme(aspect.ratio=1)
p6 <- p6 + ylab(bquote('Chl a ('*mg~ m^-3*')'))
p6 <- p6 + theme(legend.position="none")

p7 <- ggplot(atten_avg_df)
p7 <- p7 + geom_tile(aes(x=Longitude,y=Latitude,fill=kd490)) 
p7 <- p7 + scale_fill_gradientn(colors=c('#313695','#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43','#d73027','#a50026'), 
                                name=(bquote('kd490 ('*m^-1*')')))
p7 <- p7 + coord_equal() + xlab("Longitude") + ylab("Latitude")
p7 <- p7 + geom_polygon(data=belize_trans_clip, aes(x=long, y=lat, group=group))
p7 <- p7 + geom_polygon(data=coralreef_trans_clip, aes(x=long, y=lat, group=group), color = "gray", fill= "gray")
p7 <- p7 + geom_text(aes(x=Longitude,y=Latitude, label=IDnumber),data=corals)
p7 <- p7 + theme_minimal()

p8 <- ggplot(atten_summary, aes(x=Location, y=Mean, colour=Location, group = Location))
p8 <- p8 + geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=.1, color = c("#525252", "#969696"))
p8 <- p8 + geom_col(fill = c("#525252", "#969696"), color= c("#525252", "#969696"))
p8 <- p8 + theme_bw(base_size = 9)
p8 <- p8 + theme(axis.line = element_line(color = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 axis.title.x=element_blank())
p8 <- p8 + theme(aspect.ratio=1)
p8 <- p8 + ylab(bquote('kd490 ('*m^-1*')'))
p8 <- p8 + scale_y_continuous(breaks = seq(0, 0.1, 0.1), limits = c(0, 0.13))
p8 <- p8 + theme(legend.position="none")

p9 <- ggdraw() 
p9 <- p9 + draw_plot(p + theme(legend.position = "top", legend.box = "horizontal"), 0, 0, 1, 1)
p9 <- p9 + draw_plot(p2, 0.31, 0.47, 0.29, 0.29) 

p10 <- ggdraw() 
p10 <- p10 + draw_plot(p3 + theme(legend.position = "top", legend.box = "horizontal"), 0, 0, 1, 1)
p10 <- p10 + draw_plot(p4, 0.30, 0.47, 0.29, 0.29)

p11 <- ggdraw() 
p11 <- p11 + draw_plot(p5 + theme(legend.position = "top", legend.box = "horizontal"), 0, 0, 1, 1)
p11 <- p11 + draw_plot(p6, 0.31, 0.475, 0.29, 0.29)
  
p12 <- ggdraw() 
p12 <- p12 + draw_plot(p7 + theme(legend.position = "top", legend.box = "horizontal"), 0, 0, 1, 1)
p12 <- p12 + draw_plot(p8, 0.31, 0.475, 0.29, 0.29)

plot_grid(p9, p10, p11, p12)
plot_grid(p, p3, p5, plot5, nrow = 2, ncol = 2, labels="AUTO")

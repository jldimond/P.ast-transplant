#Set working directory
setwd("~/Documents/Projects/PoritesRADseq/P.ast-transplant/analyses/SymbiontGenoytping/")
#Get list of filenames to import with *.blast suffix
filenames <- list.files(pattern="*.blast", full.names=TRUE)
#import files as a list of dataframes
sym <- lapply(filenames, read.delim, sep = "\t", skip =1, header = F)
#extract sample names from filenames
sampnames <- gsub('^.*./\\s*|\\s*_.*$', '', filenames)
#add sample names to list
sym <- Map(cbind, sym, sample = sampnames)
#determine clades from blast output
clade <- lapply(sym, function(x) ifelse(grepl("sp. A|clade A", x$V3, ignore.case = T), "A", 
                       ifelse(grepl("sp. B|clade B", x$V3, ignore.case = T), "B",
                              ifelse(grepl("sp. C|clade C", x$V3, ignore.case = T), "C",
                                     ifelse(grepl("sp. D|clade D", x$V3, ignore.case = T), "D", NA)))))
#add clades to list
sym2 <- mapply(cbind, sym, clade = clade, SIMPLIFY = FALSE)
#remove instances where one read mapped to more than one clade
sym3 <- lapply(sym2, function(x) x[!duplicated(c(x$V1, x$clade)),])
#rename list elements with sample names
names(sym3) <- sampnames
#bind rows
library(dplyr)
sym4 <- bind_rows(sym3, .id="df")
#remove all instances of 18S rDNA, since this gives ambiguous
#matches to many clades
sym5 <- sym4[!grepl("18S", sym4$V3),]
#omit NA values
sym6 <- na.omit(sym5)
library(ggplot2)
#summarize data
symsum <- sym6 %>% 
  group_by(sample,clade) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))
#add vector for species
symsum$species <- c(rep("Branching Porites", 32), rep("Porites astreoides", 72))
library(ggstance)
#plot it
ggplot(subset(symsum,sample %in% c("pa10-15", "pa10-16","pa11-15", "pa11-16",
                                   "pa2-15", "pa2-16", "pa3-15", "pa3-16", "pa5-15",
                                   "pa5-16", "pa6-15", "pa6-16", "pa8-15", "pa8-16",
                                   "pa9-15", "pa9-16",
                                   "103", "104", "105", "124", "125", "126",
                                   "128", "129", "130", "131"))) + 
aes(x = perc, y = sample, fill = factor(clade)) +
geom_barh(stat="identity", width = 0.9) +
scale_x_continuous(limits = c(0,1), expand = c(0, 0))+
scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"))+
labs(x = "Proportion", fill = "Clade") +
theme_gray()+
theme(axis.text.y  = element_blank(), axis.ticks.y=element_blank(), 
axis.title.y = element_blank(), legend.position="top")+
facet_grid(species ~ ., scales = "free_y", space = "free_y")
    
  
#plot before and after transplantation
  ggplot(subset(symsum,sample %in% c("pa10-15", "pa11-15", "pa12-15", "pa17-15",
                                   "pa18-15", "pa2-15", "pa3-15", "pa5-15",
                                   "pa6-15", "pa7-15", "pa8-15", "pa9-15",
                                   "pa10-16", "pa11-16", "pa12-16", "pa17-16",
                                   "pa18-16", "pa2-16", "pa3-16", "pa5-16",
                                   "pa6-16", "pa7-16", "pa8-16", "pa9-16"))) + 
    aes(x = factor(sample), y = perc*100, fill = factor(clade)) +
    geom_bar(stat="identity", width = 0.7) +
    labs(x = "Sample", y = "percent", fill = "clade") +
    theme_minimal(base_size = 14)




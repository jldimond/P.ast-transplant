# Branching _Porites_ project repository

This is a repository to accompany the manuscript ["Genetic And epigenetic insight into morphospecies in a reef coral"](http://biorxiv.org/content/early/2017/03/22/119156).


The study involves analysis of [ddRAD-seq](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037135) and [EpiRAD-seq](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12435/abstract) data for 27 samples of branching _Porites_ spp. corals collected in Belize. There is taxonomic uncertainty whether Caribbean branching _Porites_ spp. comprise 3 different species or a single polymorphic species. The species are distinguished primarily by branch thickness.

![_Porites porites_](./images/Screen%20Shot%202016-11-03%20at%206.56.27%20PM.png)

The goal is to evaluate (1) if there is any genetic structuring of these individuals using the ddRAD data and (2) if there is any epigenetic structuring of these individuals using the EpiRAD data. 


# Methods overview

## Library preparation

A starting set of 48 samples of _Porites_ spp. were prepared for ddRADseq and EpiRADseq. Only 30 of these samples were the branching _Porites_ spp. analyzed in this study; the rest were _Porites astreoides_ destined for analysis in a separate study. A separate ddRAD and EpiRAD library was created for each sample, meaning that a total of 96 (48 ddRAD and 48 EpiRAD) samples were run on a 96 well plate. These were split into 12 pools that were run on a single Illumina HiSeq 4000 lane (100 bp paired-end reads). The [data](https://github.com/jldimond/Branching-Porites/tree/master/data) directory contains sample metadata and barcode files associated with each library. 


## Assembly 

Sequences were assembled using ipyrad v0.3.41 (Eaton 2014). We used the ‘denovo - reference’ assembly method with the Symbiodinium minutum (clade B; (Shoguchi et al. 2013)) and Symbiodinium kawagutii (clade F; (Lin et al. 2015) ) draft genomes used as reference to subtract symbiont reads from the de novo assembly. The [assembly workflow](https://github.com/jldimond/Branching-Porites/blob/master/notebooks/ipyrad_assembly.ipynb) details the assembly in a Jupyter notebook. Further specifics of the ipyrad parameters used for assembly can be seen in the aforementioned workflow or directly in the [parameters file](https://github.com/jldimond/Branching-Porites/blob/master/analyses/ipyrad_analysis/params-data3.txt). Information on these parameters and how they are used in ipyrad can be found [here](http://ipyrad.readthedocs.io/).


## Analysis 

Each sample has both a ddRAD-seq and EpiRAD-seq library associated with it. The only difference between the libraries is that the ddRAD library was generated with a methylation-insensitive common cutter (_MspI_) while the EpiRAD library was generated with a methylation-sensitive common cutter (_HpaII_). Both restriction enzymes recognize the same cut site, CCGG, but if this site is methylated, _HpaII_ will not cut and the locus will not be present in the EpiRAD library. Both libraries had the same rare cutter (_PstI_). Thus, both libraries should contain similar sets of loci unless the locus is methylated. 

In essence the analysis involves (1) analysis of single nucleotide polymorphisms (SNPs) in the ddRAD data and (2) analysis of read counts in the EpiRAD data. ipyrad generates many [output files](https://github.com/jldimond/Branching-Porites/tree/master/analyses/ipyrad_analysis/data3_outfiles) that can be readily used for the ddRAD SNP anlysis, but the EpiRAD analysis requires a bit more work to extract workable data. First, read count data are extracted from the .vcf format output from iPyrad. This file is too large to be included in the online repository, but the workflow can be found [here](https://github.com/jldimond/Branching-Porites/blob/master/notebooks/VCF_readcounts.ipynb). The remainder of the analysis is performed in R. The R notebook is located [here](https://github.com/jldimond/Branching-Porites/blob/master/analyses/ipyrad_analysis/data3_outfiles/PoritesRAD_analysis.md).


## Software

This study made use of the following software:

ipyrad v0.3.41

SnpEff(SnpSift) v4.2

R v3.1.3

RStudio v1.0.136

#### R packages:

ape v3.4

edgeR v3.8.6

adegenet v2.0.1

ade4 v1.7-4

hierfstat v0.04-22

lmtest v0.9-34

relaimpo v2.2-2


## Directory structure


`analyses/` - Files resulting from ipyrad assembly and R analyses.

`data/` -  Metadata and links to methods and raw data.

`images/` - Coral photos, gel images, figures.

`notebooks/` - Jupyter notebooks.

`scripts/` - Scripts such as R scripts used for analyses.

![_Porites astreoides_](./images/titleimage.png)

This is a repository to accompany a manuscript in preparation titled ["Response of DNA methylation to experimental transplantation in the reef coral Porites astreoides"](https://docs.google.com/document/d/1iN7kf_chLfLjr0Nh5mjNlU0gexdgP9ThgfcxaZgm1Gw/edit).


The study involves analysis of [ddRAD-seq](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037135) and [EpiRAD-seq](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12435/abstract) data for samples of _Porites astreoides_ corals in Belize. Several coral colonies were transplanted to a common garden in November 2015 and sampled again the following year.
The goal is to evaluate whether DNA methylation levels in these corals change in response to simulated environmental change. 


# Methods overview

One-year common garden transplantation experiment on Belize Barrier Reef
   colonies moved to common garden from 10-20 km away
ddRADseq coupled with EpiRADseq (ddRAD variant) to assess methylation
   EpiRADseq uses methylation-sensitive restriction enzyme
   Methylation estimated by read counts 
   Reads low/absent in EpiRADseq library = methylated
De novo assembly using ipyrad (no P. astreoides genome available yet)
   Symbiodinium Clade A genome used to remove symbiont reads

## Library preparation

ddRADseq and EpiRADseq library preparation protocols are described in detail [here](http://onsnetwork.org/jdimond/2016/08/). Most of the data used for this study was generated in a single set of ddRADseq and EpiRADseq libraries (JD002_A-L) produced in January 2017; technical details of these libraries prior to sequencing are [here](http://onsnetwork.org/jdimond/2017/02/14/rad-library-prep/). A separate ddRAD and EpiRAD library was created for each sample. These were split into 12 pools (A-L) that were run on a single Illumina HiSeq 4000 lane (100 bp paired-end reads). The [data](https://github.com/jldimond/Branching-Porites/tree/master/data) directory contains sample metadata and barcode files associated with each library. 


## Assembly 

Sequences were assembled using ipyrad v.0.5.15 (Eaton 2014). We used the ‘denovo - reference’ assembly method with the Symbiodinium microadriaticum (clade A; [Aranda et al. 2016](https://www.nature.com/articles/srep39734)) genome used as reference to subtract symbiont reads from the de novo assembly. The [assembly workflow](https://github.com/jldimond/Branching-Porites/blob/master/notebooks/ipyrad_assembly.ipynb) details the assembly in a Jupyter notebook. Further specifics of the ipyrad parameters used for assembly can be seen in the aforementioned workflow or directly in the [parameters file](https://github.com/jldimond/Branching-Porites/blob/master/analyses/ipyrad_analysis/params-data3.txt). Information on these parameters and how they are used in ipyrad can be found [here](http://ipyrad.readthedocs.io/).


## Analysis 

Each sample has both a ddRAD-seq and EpiRAD-seq library associated with it. The only difference between the libraries is that the ddRAD library was generated with a methylation-insensitive common cutter (_MspI_) while the EpiRAD library was generated with a methylation-sensitive common cutter (_HpaII_). Both restriction enzymes recognize the same cut site, CCGG, but if this site is methylated, _HpaII_ will not cut and the locus will not be present in the EpiRAD library. Both libraries had the same rare cutter (_PstI_). Thus, both libraries should contain similar sets of loci unless the locus is methylated. 

ipyrad generates many [output files](https://github.com/jldimond/Branching-Porites/tree/master/analyses/ipyrad_analysis/data3_outfiles) that can be readily used for the ddRAD SNP anlysis, but the EpiRAD analysis requires a bit more work to extract workable data. First, read count data are extracted from the .vcf format output from iPyrad. This file is too large to be included in the online repository, but the workflow can be found [here](https://github.com/jldimond/Branching-Porites/blob/master/notebooks/VCF_readcounts.ipynb). The remainder of the analysis is performed in R. The R notebook is located [here](https://github.com/jldimond/Branching-Porites/blob/master/analyses/ipyrad_analysis/data3_outfiles/PoritesRAD_analysis.md).


## Software

This study made use of the following software:

ipyrad v.0.5.15

R v3.1.3

RStudio v1.0.136


#### R packages:

ape v3.4

edgeR v3.8.6

adegenet v2.0.1

ade4 v1.7-4

hierfstat v0.04-22

iheatmapr v0.2.5

reshape2 v1.4.2

ggplot2 v2.2.1


## Directory structure


`analyses/` - Files resulting from ipyrad assembly and R analyses.

`data/` -  Metadata and links to methods and raw data.

`images/` - Coral photos, gel images, figures.

`notebooks/` - Jupyter notebooks.

`scripts/` - Scripts such as R scripts used for analyses.

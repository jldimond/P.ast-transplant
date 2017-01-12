#FISH546 2016 project repository

This is a repository for my 2016 FISH 546 project, involving analysis of _Porites_ spp. ddRAD-seq and EpiRAD-seq data.

I have recently obtained [ddRAD-seq](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0037135) and [EpiRAD-seq](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12435/abstract) data for 48 samples of _Porites_ spp. corals collected in Belize. Most of these samples consist of branching _Porites_ spp., for which there is taxonomic uncertainty whether these comprise 3 different species or a single polymorphic species. The species are identitied primarily by branch thickness.

![_Porites porites_](./images/Screen%20Shot%202016-11-03%20at%206.56.27%20PM.png)

My goal is to evaluate (1) if there is any genetic structuring of these individuals using the ddRAD data and (2) if there is any epigenetic structuring of these individuals using the EpiRAD data. I will be using iPyrad to obtain the clusters of loci, and will be filtering out symbiont sequences by using the "denovo-reference" option to obtain sequences that do not match the _Symbiodinium minutum_ and _Symbiodinium kawagutii_ genomes.

##Assembly overview

My first attempt at an assembly used the `reference` option in iPyrad, with a _Porites astreoides_ transcriptome as reference. This assembly returned a very small number of number of loci, so I decided to consider other options. My second attempt used the `denovo - reference` assembly method with the _Symbiodinium minutum_ draft genome as reference. This returned many more loci, but in the end I opted to use both the _Symbiodinium minutum_ (clade B) and _Symbiodinium kawagutii_ (clade F) draft genomes to maximize the number of symbiont reads subtracted from the assembly. These genomes should enable removal of any highly conserved symbiont sequences, but what about potential sequences that are specific to the _Symbiodinium_ clade A and C symbionts found in my samples*? My rationale here is that by refining the final RAD assembly to include only loci present in all samples (i.e., no missing data), this should exclude less conserved sequences specific to either _Symbiodinium_ clades A or C. This is the approach currently adopted in the [assembly workflow](https://github.com/jldimond/jldimond-fish546-2016/blob/master/notebooks/ipyrad_assembly.ipynb). Further specifics of the parameters used for assembly can be seen in the aforementioned  workflow or directly in the [parameters file](https://github.com/jldimond/jldimond-fish546-2016/blob/master/analyses/ipyrad_analysis/params-data3.txt). Information on these parameters and how they are used in iPyrad can be found [here](http://ipyrad.readthedocs.io/).

*symbionts were genotyped via cp23S Sanger sequencing earlier this year

##Analysis overview

Each sample has both a ddRAD-seq and EpiRAD-seq library associated with it. The only difference between the libraries is that the ddRAD library was generated with a methylation-insensitive common cutter (_MspI_) while the EpiRAD library was generated with a methylation-sensitive common cutter (_HpaII_). Both restriction enzymes recognize the same cut site, CCGG, but if this site is methylated, (_HpaII_) will not cut and the locus will not be present in the EpiRAD library. Both libraries had the same rare cutter (_PstI_). Thus, both libraries should contain similar sets of loci unless the locus is methylated. 

In essence the analysis involves (1) analysis of single nucleotide polymorphisms (SNPs) in the ddRAD data and (2) analysis of read counts in the EpiRAD data. iPyrad provides many useful [output files](https://github.com/jldimond/jldimond-fish546-2016/tree/master/analyses/ipyrad_analysis/data3_outfiles) that can be readily used for the ddRAD SNP anlysis, but the EpiRAD analysis requires a bit more work to extract workable data. First, read count data are extracted from the .vcf format output from iPyrad. This is a large file that is not included in the online repository, but the workflow can be found [here](https://github.com/jldimond/jldimond-fish546-2016/blob/master/notebooks/VCF_readcounts.ipynb). Next, data are imported into R for analysis. Samples with lots of missing data are excluded, then rows in ddRAD libraries with any zeros are excluded. The premise here is that zeros in the EpiRAD dataset are informative because they may reflect methylation, but they could also reflect true absence of the locus in the library. Here the ddRAD library serves to standarize the EpiRAD library. Any zeros in the ddRAD libary are treated as absence of the locus, thereby leaving zeros in the EpiRAD library only where the locus was counted in the ddRAD library. Next, these data are standardized by library size using TMM normalization in the edgeR package. Currently, I am looking at approximately 1400 loci shared across all samples, with approximately 15% of these methylated. The R script for the EpiRAD analysis is located [here](https://github.com/jldimond/jldimond-fish546-2016/blob/master/scripts/EpiRAD_analysis.R), and the script for the ddRAD analysis is [here](https://github.com/jldimond/jldimond-fish546-2016/blob/master/scripts/ddRAD_analysis.R).

From here, the analysis is still evolving. I am experimenting with using the [residuals](https://github.com/jldimond/jldimond-fish546-2016/blob/master/analyses/residuals.pdf) from a linear model of ddRAD libraries vs. EpiRAD libraries as the working data for the EpiRAD analysis. Moving forward, I will implement principal components analysis to both assess dispersion of samples and evaluate differentially methylated loci. A heat map of differentially methylated loci would be nice. 

For the ddRAD SNP analysis, I will also use MDS/PCA approaches, as well as phylogenetic tree approaches. I am looking into methods used for species delimitation.

##Directory structure

The directory structure is as follows:

`analyses/` - Files resulting from analyses.

`data/` -  Information and links to methods and raw data.

`notebooks/` - Jupyter notebooks.

`scripts/` - Scripts such as R scripts used for analyses.

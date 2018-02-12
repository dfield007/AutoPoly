# AutoPoly 

This is the repository for the AutoPoly, an R based package that performs population assignment in autotetraploids and autohexaploids using genetic marker information. The analysis uses allele frequency information from reference genotypes sampled from candidate (source) populations, and then assigns a set of genotyped (or phenotyped) individuals of unknown origin (i.e. offspring) to their most likely: 

* joint maternal and paternal population, when maternal genotype is unknown (e.g. seed dispersal), or 
* paternal population, given the known maternal genotype of each offspring (e.g. pollen dispersal).

# Quick guide

To install, one approach is from the command line:

  `R CMD INSTALL AutoPoly_0.2.0.tar.gz`

Or you can start R and download from R cran (if updated):

  `install.packages('AutoPoly')`

Next start R and run package:

  `library(AutoPoly)`

Now your set to go. Please look at the Vignette for examples on how to run data sets through the package.

# News/updates

12/2/2018 - Updated for R v3.4.3, AutoPoly_0.2.0.tar.gz

04/10/2017 - Manuscript published online in Heredity - Population assignment in autopolyploids
https://www.nature.com/articles/hdy201751




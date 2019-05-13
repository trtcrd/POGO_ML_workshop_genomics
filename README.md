# POGO_ML_workshop_genomics

***Genomics tutorial for the POGO 2019 Workshop on Machine Learning and Artificial Intelligence in Biological Oceanographic Observations, Tristan Cordier and Anders Lanzén, Ostend, Belgium, 20-22th May 2019.***


## Part 1 - Prediction of environmental impact using metabarcoding

This tutorial will cover training and cross-validation of a random forest classifier, to predict environmental impact from fish farming, based on metabarcoding (amplicon sequencing) of microbial benthic communities. AMBI (AZTI Biotic Index) and other measured indices based on visual observation of macrobenthic species composition is used as a "ground truth". The tutorial is done completely in R, by connecting to RStudio running in an instance launched at [AWS](https://eu-west-3.console.aws.amazon.com/ec2).

__The R script [genomics_tutorial.R](genomics_tutorial.R) covers Part 1 of the Tutorial.__

You can read more about the analysed dataset in [Cordier et al. 2017](dx.doi.org/10.1021/acs.est.7b01518) and [Cordier et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12926).


## Part 2 – Trait prediction of MAGs

The first part of the tutorial focused on metabarcoding data, i.e. profiling of the __community structure__ using rRNA as a __taxonomic marker__. As you saw, this can be a powerful method for predicting functional traits of the community as a whole, for example its adaptation to impacts from fish farming, specifically eutrophication. There are specific tools for predicting the functional traits of communities based only on taxonomic structure, such as [PICRUST](https://picrust.github.io/picrust/) and [Tax4Fun](http://tax4fun.gobics.de/). These use ancestral state reconstruction and other supervised methods (though not strictly ML) on training data from sequenced genomes, to predict functional gene content from rRNA data. However, the most accurate is to profile traits by sequencing communities and genomes directly.

__Metagenomics__,  shotgun sequencing of all environmental DNA in a sample, is still too expensive and complex for routine tasks like monitoring of aquaculture sites, but is fast becoming a more popular tool in marine ecology (along with metatranscriptomics, single cell genomics, proteomics and other molecular high throughput methods). However, the large datasets generated are complex and demanding to analyse manually, making them ideal targets for machine learning. Thanks to dropping sequencing costs, public metagenome data is now accumulating faster than our ability to study it in depth.

Instead of focusing on the whole community, we will "zoom in" on a number of __Metagenome Assembled Genomes (MAGs)__ from the [__Tara Oceans Expeditions__](http://oceans.taraexpeditions.org). We will focus on only 15 MAGs from nitrogen-fixing bacterioplankton, out of nearly a thousand MAGs, painstakingly [curated from automatic binning](http://merenlab.org/data/tara-oceans-mags/). Binning, by the way, also is another good example of an application for ML. Have a look at [the paper by Delmont et al.](dx.doi.org/10.1038/s41564-018-0176-9) to learn more about these intersting microbes. 

We will use two different tools for __predicting and analysing functional capabilities__, or traits, of the organisms represented by the MAGs:  __[Traitar](https://github.com/hzi-bifo/traitar/blob/master/INSTALL.md)__ published by [Weimann et al, in 2016](https://msystems.asm.org/content/1/6/e00101-16), and [a recently published method by Farrell et al](https://www.biorxiv.org/content/10.1101/307157v1).

To start the tutorial, connect to your AWS instance using an SSH enabled terminal (unless you are in Linux or OSX, you can use PuTTY, Cygwin or similar on Windows). You can see the DNS or IP in the same way as for Part 1. (for full login details right-click on the instance in the AWS EC2 console and select Connect, which will display something like `ssh -i my_aws_key.pem ubuntu@ec2-35-181-51-240.eu-west-3.compute.amazonaws.com`).


Traitar uses __[support vector machines (SVMs)](https://en.wikipedia.org/wiki/Support_vector_machine)__, which are a class of algorithms using linear N-dimensional hyperplanes based on supervised training. They can also use a [kernel function](https://en.wikipedia.org/wiki/Kernel_method) or "data transformation" to be turned into non-linear methods. 



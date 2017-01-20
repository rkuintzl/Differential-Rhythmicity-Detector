# LLC_Analysis (DRD?) :
## [Full package name -- (differential rhythmicity detector)?]

LICENSE?

Authors:
Rachael C. Kuintzle (rkuintzl@caltech.edu), California Institute of Technology
David A. Hendrix (David.Hendrix@oregonstate.edu), Oregon State University

This document is intended to give a minimal description of the DRD pipeline, a collection of Perl scripts used to identify genes with differential rhythmicity based on round-the-clock RNA expression data from two different conditions (in this case, young and old age), as described in
[REFERENCE]

Initial rhythmicity assessment for individual genes in a given condition is accomplished upstream of this pipeline using the package ARS.R:
[REFERENCE]

We are sharing this code in the spirit of experimental transparency and open science.
The scripts in this collection are not generalized; they are tailored specifically
for our data using custom file formats that are detailed below.
However, you are free to adapt this code for your own applications, and to send any comments/questions to David.Hendrix@oregonstate.edu. Thanks!

## Introduction

Currently, [DRD] is composed of 3 modules [RENAME!!!]:

        * parseArser_get_LLCscores_basedOnZ-scores.pl   :   Compute S_DRs (differential rhythmicity scores)
        * fitGaussian_DAH_script.pl    :   Fit distribution of S_DRs to a Gaussian curve; obtain mean and variance parameters
        * parseArser_get_LLCscore_p-values_and_FDR : [include only the "print_everything" one?]  Assign p-values and Benjamini-Hochberg corrected p-values (q-values) to each S_DR

To get help on each module, you can type :

        compute_DR_scores.pl --help

## Input files

The formats used to describe genes, transcripts, exon is **.GTF** and **.FASTA** for genome file.

Basically, FEELnc users should have the following minimal input files:

        - Infile.GTF          (-i,--infile)   : input GTF file (e.g cufflinks/stringtie transcripts.GTF)
        - ref_annotation.GTF  (-a,--mRNAfile) : GTF annotation file*
        - ref_genome.FASTA    (-g,--genome)   : genome FASTA file or directory with individual chrom FASTA files


    -------------------------
## Installation and requirements

### Requirements

The following software and libraries must be installed on your machine:

- [Perl5+](https://www.perl.org/) : tested with version 5.18.2
 * [Bioperl](http://www.bioperl.org/wiki/Main_Page)  : tested with version BioPerl-1.6.924
- R [Rscript](http://cran.r-project.org): tested with version 3.1.0.

### Installation?????

To run these scripts, you need simply to download the files, and execute them in the command line.

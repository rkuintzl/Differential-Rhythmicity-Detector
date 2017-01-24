# DRD :
## Differential Rhythmicity Detector

Authors:
Rachael C. Kuintzle (rkuintzl@caltech.edu), California Institute of Technology
David A. Hendrix (David.Hendrix@oregonstate.edu), Oregon State University

This document is intended to give a minimal description of the DRD pipeline, a collection of Perl scripts used to identify genes with differential rhythmicity based on round-the-clock RNA expression data from two different conditions (in this case, young and old age), as described in:
Kuintzle, R. C. et al. Circadian deep sequencing reveals stress-response genes that adopt robust rhythmic expression during aging. Nat. Commun. 8, 14529 doi: 10.1038/ncomms14529 (2017).

We are sharing this code in the spirit of experimental transparency and open science.
The scripts in this collection are tailored to run downstream of the R implementation of ARSER contained in the R library [MetaCycle](https://cran.r-project.org/web/packages/MetaCycle/MetaCycle.pdf), which detects rhythmic patterns in time-series data. You are free to adapt our code for your own applications, and to send any comments/questions to David.Hendrix@oregonstate.edu. Thanks!


## Introduction

DRD is composed of 3 Perl scripts:

        * compute_DR_scores.pl   :   Compute S_DRs (differential rhythmicity scores)
        * DR_scores_Gaussian_fit.pl    :   Fit histogram of S_DRs to a Gaussian distribution; obtain mean and variance parameters
        * compute_DR_p-values_and_FDR.pl  :  Assign p-values and Benjamini-Hochberg corrected p-values (q-values) to each S_DR -- larger S_DRs are more significant

Additionally, a module called Histogram.pm is provided. This module is used by the script DR_scores_Gaussian_fit.pl.

To get help on each script, you can type :

        compute_DR_scores.pl --help


## Input files

The input files used in our first script are modified versions of the output file produced by the R implementation of ARSER. ARSER generates tab-delimited text files with user-provided, unique gene or transcript IDs as the first column. To prepare an ARSER output file for use in our DRD pipeline, merely add a column containing unique gene symbols to the beginning of the file. Alternatively, you can simply duplicate the ID column if you do not wish to provide gene symbols.

[MetaCycle](https://cran.r-project.org/web/packages/MetaCycle/MetaCycle.pdf) contains the R implementation of the program [ARSER](https://github.com/cauyrd/ARSER). The function, called "meta2d", is available as part of the R library "MetaCycle". When this function is executed, "ARS" should be specified as the parameter "cycMethod".

-------------------------

## Dependencies

The following software and libraries must be installed on your machine:

- [Perl5+](https://www.perl.org/) : tested with version [5.18.2] DAVE?
 * [Math::CDF](http://search.cpan.org/~callahan/Math-CDF-0.1/CDF.pm)  :
 * [Statistics::Basic](http://search.cpan.org/~jettero/Statistics-Basic-1.6611/lib/Statistics/Basic.pod)
 * [Getopt::Long](http://search.cpan.org/~jv/Getopt-Long-2.49.1/lib/Getopt/Long.pm)

Additionally, the input files with rhythmicity parameters for each gene are modified slightly from those produced by the R implementation of ARSER, as described above. Consult the [MetaCycle manual](https://cran.r-project.org/web/packages/MetaCycle/MetaCycle.pdf) for information on its installation and dependencies.


## Usage

To run these scripts, you need only to download the Perl files and execute them in the command line. Usage and example commands are shown for each script below.


### compute_DR_scores.pl

perl ./compute_DR_scores.pl <ARSER output for condition 1> <ARSER output for condition 2> [options]

Example usage:
perl ./compute_DR_scores.pl namesAdded_young_dm6_genes_ARSER_result.txt namesAdded_old_dm6_genes_ARSER_result.txt -m 1

The output file produced by this script will be used in subsequent scripts in the pipeline.


### DR_scores_Gaussian_fit.pl

perl ./DR_scores_Gaussian_fit.pl -i <output file from compute_DR_scores.pl> -b <bin size> -o <outfile prefix> -c <column in text file to use for histogram building> -R

Example usage:
perl ./DR_scores_Gaussian_fit.pl -i DR_scores_medianFpkmCutoff1.txt -o DR_scores_binSize0.1 

Note that this script uses the Histogram.pm module that is also provided. It should work if kept in the same directory as the script; otherwise, put it into your PERL5LIB path. 

Output: Consider primarily the *.fitParams file, which reports the Gaussian mean and variance needed to run the third and final script. There are also histogram files produced, as well as information on the Gaussian fit. 


### compute_DR_p-values_and_FDR.pl

perl ./compute_DR_p-values_and_FDR.pl <output file from compute_DR_scores.pl> -m <DR mean (Gaussian)> -v <DR variance (Gaussian)> [options]

Note: The mean and variance parameters are reported by the script from step 2: DR_scores_Gaussian_fit.pl, as part of the *.fitParams file. 

Example usage:
perl ./compute_DR_p-values_and_FDR.pl DR_scores_medianFpkmCutoff1.txt -m -0.0270973 -v 0.859165 -f 0.01

Note: The default FDR is 0.05 when no -f flag is provided.

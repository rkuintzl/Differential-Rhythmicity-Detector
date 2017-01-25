#!/usr/bin/perl -w

# Copyright 2016 David Hendrix

# This file is part of the DRD script collection.

#    DRD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    DRD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DRD.  If not, see <http://www.gnu.org/licenses/>.
#

# Author:
# David A. Hendrix (David.Hendrix@oregonstate.edu), Oregon State University

# This script is one in a series of scripts for characterizing differential rhythmicity in
# gene expression analysis as described in
# Kuintzle, R. C. et al. Circadian deep sequencing reveals stress-response genes that adopt robust rhythmic expression during aging. Nat. Commun. 8, 14529 doi: 10.1038/ncomms14529 (2017).

# Please use DR_scores_Gaussian_fit.pl -h to see the help screen for further instructions on running this script.

use Getopt::Long;
use Histogram;
use strict;
$|=1;

my $usage = "
Usage:\n\n
        $0 -i <DR scores file> [options]\n

Options:\n
        -i --infile      the input file of data points to analyze -- the output file produced by the script compute_DR_scores.pl.\n
        -o --out         output file base.\n
        -b --binsize     the size of bins for the histogram. (default=0.1)\n
        -l --lower       a lower bound on the fitting.
        -u --upper       an upper bound on the fitting.
        -m --max         an upper-bound on what is read in. filter out data greater than max.\n\n";

my $data = new Histogram;

# set defaults
$data->col(11);
$data->bin(0.1);
$data->randomize(1);

Getopt::Long::Configure("no_ignore_case");
# read in the options
GetOptions( 'infile=s' => \$data->{_in},
	    'out=s' => \$data->{_out},
	    'binsize=f' => \$data->{_bin},
	    'pseudocount=f' => \$data->{_pseudoCount},
	    'lower=f' => \$data->{_fitMin},
	    'uppper=f' => \$data->{_fitMax},
	    'floor' => \$data->{_lower},   # an lowerbound(filter) on what is read in.
	    'max=f' => \$data->{_upper},   # an upperbound(filter) on what is read in.
 	    'help' => \$data->{_help},
	    'raw' => \$data->{_raw},
	    'Log' => \$data->{_logTransform},
	    'Randomize' => \$data->{_randomize},
	    'Cummulative' => \$data->{_Cummulative}
    );

if($data->help()) {
    die $usage;
}

# read the data, one row/line at a time.
if($data->in()) {
    open(INFILE,$data->in()) or die "could not open ",$data->in(),"\n";
    while(<INFILE>) {
	unless(/^\#/) {
	    chomp;
	    my @row = split;
	    $data->tally($row[$data->{_col}-1]);
	}
    }
} else {
    die $usage;
}
my $histoFile = $data->printHistogram();
my $fitParams = $data->fitGaussianDistribution($histoFile);
$data->printGaussianFit($fitParams);
$data->printGaussianFitParameters($fitParams);

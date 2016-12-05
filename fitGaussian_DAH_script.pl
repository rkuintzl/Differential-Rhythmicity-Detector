#!/usr/bin/perl -w
use Getopt::Long;
use Histogram;
use strict;
$|=1;


my $usage = "
Usage:\n\n
        $0 -i <infile> [options]\n
        or\n
        cat infile | perl $0 [options]\n\n
Options:\n
        -i --infile      the input file of data points to analyze.\n
        -o --out         output file base.\n
        -b --binsize     the size of bins for the histogram.\n
        -c --column      the column to evaluate the histogram of.\n
        -l --lower       a lower bound on the fitting.
        -u --upper       an upper bound on the fitting.
        -m --max         an upper-bound on what is read in. filter out data greater than max.\n
        -r --raw         print the raw (not normalized) histogram. [1=true,0=false=default]\n
        -L --Log         log transform the input data. [1=true,0=false=default]\n
        -R --Randomize   randomize the initial value for the mean in the fit.
        -C --Cummulative print the cummulative distribution instead of histogram\n
        -p --pseudocount pseudo count to add when taking log-scale\n\n";


my $data = new Histogram;

Getopt::Long::Configure("no_ignore_case");
# read in the options
GetOptions( 'infile=s' => \$data->{_in},
	    'out=s' => \$data->{_out},
	    'column=i' => \$data->{_col},
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
    while(<STDIN>) {
	unless(/^\#/) {
	    chomp;
	    my @row = split;
	    $data->tally($row[$data->{_col}-1]);
	}
    }
}
my $histoFile = $data->printHistogram();
my $fitParams = $data->fitGaussianDistribution($histoFile);
$data->printGaussianFit($fitParams);
$data->printGaussianFitParameters($fitParams);

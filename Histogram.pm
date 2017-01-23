# Class Histogram
package Histogram;
use Graphics::GnuplotIF qw(GnuplotIF);
use Math::CDF;
use strict;
return 1;
$|=1;

sub new {
    # object constructor
    # create the initial structure to hold the data.
    my $class = shift;
    my $self = {};
    # parameters
    $self->{_in} = "";
    $self->{_out} = "";
    $self->{_col} = 1;
    $self->{_bin} = 1;
    $self->{_lower} = "";
    $self->{_upper} = "";
    $self->{_fitMin} = "";
    $self->{_fitMax} = "";
    $self->{_unnormalized} = 0;
    $self->{_logTransform} = 0;
    $self->{_randomize} = 0;
    $self->{_help} = "";
    # measured values
    $self->{_max} = -1e1000;
    $self->{_min} = 1e1000;
    $self->{_tot} = 0;
    $self->{_sum} = 0;
    $self->{_ss} = 0;
    $self->{_histo} = {};
    bless($self,$class);
    return $self;
}

sub tally {
    # this method adds one value to the "self".
    my $self = shift;
    my($value) = @_;
    if(defined($value)) {
	# include, by default. don't if out of range.
	my $INCLUDE = 1;
	if($self->upper()) {
	    if($self->upper() < $value) {
		$INCLUDE = 0;
	    }
	}
	if($self->lower()) {
	    if($value < $self->lower()) {
		$INCLUDE = 0;
	    }
	}
	if($INCLUDE) {
	    $self->updateHisto($value);
	    $self->updateMax($value);
	    $self->updateMin($value);
	    $self->increment($value);
	}
    }
}

sub printHistogram {
    my $self = shift;
    my $pseudoCount = $self->pseudoCount() ? $self->pseudoCount : 0;
    # print the histogram.    
    my $suffix = $self->Cummulative() ? ".cum.hist" : ".hist";
    my $outFile = $self->in() ? $self->in().".c".$self->col().$suffix : "histogram.txt";
    $outFile = $self->out() ? $self->out().$suffix : $outFile;
    open(OUTFILE,">$outFile") or die "could not open $outFile for writing.\n";
    # compute the binIndex of the max and min values.
    my $minX = int sprintf("%.0f",$self->min()/$self->bin()-0.49999999);
    my $maxX = int sprintf("%.0f",$self->max()/$self->bin()-0.49999999);
    my $sum = 0;
    foreach my $x ($minX..$maxX) {
        $self->{_histo}->{$x} = $self->{_histo}->{$x} ? $self->{_histo}->{$x} : 0;
	#my $Xval = $self->logTransform() ? 10**($x*$self->bin())-$pseudoCount : $x*$self->bin();
	my $Xval = $x*$self->bin();
	unless($self->unnormalized()) { $self->{_histo}->{$x} /= $self->tot(); }	
	if($self->Cummulative()) {
	    $sum += $self->{_histo}->{$x} ? $self->{_histo}->{$x} : 0.000;
	    printf(OUTFILE "%.3e\t%.3e\n",$Xval,$sum);
	} else {
	    printf(OUTFILE "%.3e\t%.3e\n",$Xval,$self->{_histo}->{$x} ? $self->{_histo}->{$x} : 0.000);
	}
    }
    close(OUTFILE);
    # print the data summary file.
    my $summary = $self->in() ? $self->in().".c".$self->col().".summary" : "summary.txt";
    $summary = $self->out() ? $self->out().".summary" : $summary;
    open(SUMMARY,">$summary") or die "could not open $summary for writing.\n";
    printf(SUMMARY "mean=%.3f\n",$self->mean());
    printf(SUMMARY "var=%.3f\n",$self->var());
    printf(SUMMARY "min=%.3f\n",$self->min());
    printf(SUMMARY "max=%.3f\n",$self->max());
    close(SUMMARY);
    return $outFile;
}

sub fitGaussian {
    my $self = shift;
    # print the histogram. This will normalize it if needed.
    my $histoFile = $self->printHistogram();
    my $mean = $self->mean();
    my $var = $self->var();
    my $height = $self->height();
    my $fitParameters = $self->fitGaussianDistribution($histoFile);
    $self->printGaussianFit($fitParameters);
    $self->printGaussianFitParameters($fitParameters);
    return $fitParameters;
}

sub fitGaussianDistribution {
    my $self = shift;
    my($histoFile) = @_;    
    # set the initial conditions for the fit.
    my $mean = $self->mean();
    if($self->randomize()) {
	$mean = 0.5*(rand() - 1.0);
    }
    my $var = $self->var();
    my $height = $self->height();
    #print "entered fitDistribution with $histoFile $mean $var $height\n";
    # define the commands for fitting the distribution.
    my $fitLog = $histoFile.".fit.log";
    my $range1 = $self->fitMin();
    my $range2 = $self->fitMax();
    my $funcDef =  "f(x) = height*exp(-((x-mean)**2)/(2*var))";
    my $init = "mean=$mean; var=$var; height=$height;";
    my $command = "fit [$range1:$range2] f(x) \'$histoFile\' via mean,var,height";    
    # run gnuplot.
    my $gnuplot  = Graphics::GnuplotIF->new; 
    $gnuplot->gnuplot_cmd($funcDef);
    sleep(5);
    $gnuplot->gnuplot_cmd($init);
    sleep(5);
    $gnuplot->gnuplot_cmd($command);
    sleep(5); # sleep to ensure the fitting is complete.
    system("mv fit.log $fitLog");
    return $self->readFitParameters($histoFile);
}

sub readFitParameters {
    my $self = shift;
    my $histoFile = shift;
    # by construction, the fit log is defined in terms of the histoFile above.
    my $fitLog = $histoFile.".fit.log";
    #print "attempting to open $fitLog\n";
    my %parameters;
    open(FL,$fitLog) or die "could not open $fitLog\n";
    while(<FL>) {
        chomp;
        if(/FIT:    data read from \'$histoFile\'/) {
            my $READ=1;
            while($READ) {
                $_ = <FL>;
                chomp;
                if(/(\S+)\s+\=\s+(\S+)\s+\+\/\-\s+(.*)\s+\((.*)\%\)/) {
                    my $name = $1;
                    my $value = $2;
                    #print "read in: $name $value\n";
                    $parameters{$name} = $value;
                } elsif(/correlation/) {
                    $READ=0;
                } elsif(/BREAK/) {
                    die "reached breakpoint:\n$_\n";
                }       
            }
        }
    }
    die "Could not read in parameters from $fitLog\n" unless(keys(%parameters));
    return \%parameters;
}

sub printGaussianFit {
    my $self = shift;
    my($fitParameters) = @_;
    my $outFile = $self->in() ? $self->in().".c".$self->col().".Gfit" : "gaussianFit.txt";
    $outFile = $self->out() ? $self->out().".Gfit" : $outFile;
    open(FD,">".$outFile);
    my $height = $fitParameters->{height};
    my $mean = $fitParameters->{mean};
    my $var = $fitParameters->{var};
    my $sum = 0;
    my $minX = int sprintf("%.0f",$self->min()/$self->bin()-0.49999999);
    my $maxX = int sprintf("%.0f",$self->max()/$self->bin()-0.49999999);
    # x is over the indices of the histo, not actual values.
    foreach my $x ($minX..$maxX) {
	# X (capital) is the actual value.
	my $X = $x*$self->bin();
	$self->{_histo}->{$x} = $self->{_histo}->{$x} ? $self->{_histo}->{$x} : 0;
	# use $X (rather than $x here) $X is the actual x value
        my $f = $height*exp((-($X-$mean)**2)/(2*$var));
        my $d = $self->{_histo}->{$x};
        my $diff = $d - $f;
        $sum += $diff;
        print FD "$X\t$f\t$d\t$diff\t$sum\n";
    }
    close(FD);    
}

sub printGaussianFitParameters {
    my $self = shift;
    my($fitParameters) = @_;
    my $outFile = $self->in() ? $self->in().".c".$self->col().".fitParams" : "gaussianFitParams.txt";
    $outFile = $self->out() ? $self->out().".fitParams" : $outFile;
    open(FP,">".$outFile);
    my $height = $fitParameters->{height};
    my $mean = $fitParameters->{mean};
    my $var = $fitParameters->{var};
    print FP "mean = $mean\nvar = $var\nstdDev = ", sqrt($var), "\nheight = $height\n";
}

sub normalPValue {
    my $self = shift;
    my($x,$mean,$stdDev) = @_;
    my $z = ($x - $mean)/$stdDev;
    return 1.0 - Math::CDF::pnorm($z);
}

sub binomialPValue {
    # all scaling has already been done here.
    my $self = shift;
    my($expCounts,$ctlCounts) = @_;
    if(my $cdf = Math::CDF::pbinom($expCounts-1,int($expCounts+$ctlCounts),0.5)) {
	return 1 - $cdf;
    } else {
        #print STDERR "...Could not get cdf for $expCounts $ctlCounts\n";
    }
    return 1.0;
}

###################################
##### BASIC RETRIEVAL METHODS #####
###################################

sub valueToBinIndex {
    my $self = shift;
    my($value) = @_;
    my $index = int sprintf("%.0f",$value/$self->bin()-0.49999999);
    return $index;
}

sub binIndexToValue {
    my $self = shift;
    my($binIndex) = @_;
    my $value = $binIndex*$self->bin();
    return $value;
}

sub unnormalized {
    my $self = shift;
    if(@_) {
	$self->{_unnormalized} = shift;
    }
    return $self->{_unnormalized};
}

sub logTransform {
    my $self = shift;
    if(@_) {
	$self->{_logTransform} = shift;
    }
    return $self->{_logTransform};
}

sub randomize {
    my $self = shift;
    if(@_) {
	$self->{_randomize} = shift;
    }
    return $self->{_randomize};
}

sub help {
    my $self = shift;
    if(@_) {
	$self->{_help} = shift;
    }
    return $self->{_help};
}

sub in {
    my $self = shift;
    if(@_) {
	$self->{_in} = shift;
    }
    return $self->{_in};
}

sub out {
    my $self = shift;
    if(@_) {
	$self->{_out} = shift;
    }
    return $self->{_out};
}

sub Cummulative {
    my $self = shift;
    if(@_) {
	$self->{_Cummulative} = shift;
    }
    return $self->{_Cummulative};
}

sub col {
    my $self = shift;
    if(@_) {
	$self->{_col} = shift;
    }
    return $self->{_col};
}

sub bin {
    my $self = shift;
    if(@_) {
	$self->{_bin} = shift;
    }
    return $self->{_bin};
}

sub pseudoCount {
    my $self = shift;
    if(@_) {
	$self->{_pseudoCount} = shift;
    }
    return $self->{_pseudoCount};
}

sub upper {
    my $self = shift;
    if(@_) {
	$self->{_upper} = shift;
    }
    return $self->{_upper};
}

sub lower {
    my $self = shift;
    if(@_) {
	$self->{_lower} = shift;
    }
    return $self->{_lower};
}

sub max {
    my $self = shift;
    if(@_) {
	$self->{_max} = shift;
    }
    return $self->{_max};
}

sub updateMax {
    my $self = shift;
    my $pseudoCount = $self->pseudoCount() ? $self->pseudoCount : 0;
    if(@_) {
	my $value = shift;	
	$value = log($value+$pseudoCount)/log(10) if($self->logTransform());
	$self->max($value) if($value > $self->max());
    }
    return $self->max();
}

sub min {
    my $self = shift;
   if(@_) {
	$self->{_min} = shift;
    }
    return $self->{_min};
}

sub updateMin {
    my $self = shift;
    my $pseudoCount = $self->pseudoCount() ? $self->pseudoCount : 0;
    if(@_) {
	my $value = shift;	
	$value = log($value+$pseudoCount)/log(10) if($self->logTransform());
	$self->min($value) if($value < $self->min());
    }
    return $self->min();
}

sub fitMin {
    # a lower bound used in fitting.
    my $self = shift;
    if(@_) {
	$self->{_fitMin} = shift;
    }
    return $self->{_fitMin};
}

sub fitMax {
    # an upper bound used in fitting.
    my $self = shift;
    if(@_) {
	$self->{_fitMax} = shift;
    }
    return $self->{_fitMax};
}

sub tot {
    my $self = shift;
    if(@_) {
	$self->{_tot} = shift;
    }
    return $self->{_tot};
}

sub sum {
    my $self = shift;
    if(@_) {
	$self->{_sum} = shift;
    }
    return $self->{_sum};
}

sub ss {
    my $self = shift;
    if(@_) {
	$self->{_ss} = shift;
    }
    return $self->{_ss};
}

sub increment {
    my $self = shift;
    if(@_) {
	my $value = shift;
	$self->{_tot}++;
	$self->{_sum} += $value;
	$self->{_ss} += $value*$value;
    }
}

sub histo {
    # the hash reference to the histogram.
    my $self = shift;
    if(@_) {
	$self->{_histo} = shift;
    }
    return $self->{_histo};
}

sub updateHisto {
    my $self = shift;
    my $pseudoCount = $self->pseudoCount() ? $self->pseudoCount : 0;
    if(@_) {
	my($value) = @_;
	$value = log($value+$pseudoCount)/log(10) if($self->logTransform());
	my $x = int sprintf("%.0f",$value/$self->bin()-0.49999999);
	$self->{_histo}->{$x}++;
    }
}

sub mean {
    my $self = shift;
    if($self->tot()) {
	return $self->sum() / $self->tot();
    }
    return 0.0;
}

sub var {
    my $self = shift;
    if($self->tot()) {
	return ($self->ss() - $self->sum()*$self->sum()/$self->tot())/($self->tot() - 1);	
    } 
    return 0.0;
}

sub height {
    my $self = shift;
    my $maxX = 0;
    my $height = 0;
    foreach my $x (keys %{$self->{_histo}}) {
	if($height < $self->{_histo}->{$x}) {
	    $maxX = $x;
	    $height = $self->{_histo}->{$x};
	}
    }
    return $height;
}

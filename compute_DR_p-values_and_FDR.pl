#!/usr/bin/perl -w

# Copyright 2016 Rachael Kuintzle

# This file is part of the DRD script collection.

#    DRD is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    DRD is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with DRD.  If not, see <http://www.gnu.org/licenses/>.  

# Author:
# Rachael C. Kuintzle (rkuintzl@caltech.edu), California Institute of Technology                                                                                
# This script is one in a series of scripts for characterizing differential rhythmicity in
# gene expression analysis as described in:
# Kuintzle, R. C. et al. Circadian deep sequencing reveals stress-response genes that adopt robust rhythmic expression during aging. Nat. Commun. 8, 14529 doi: 10.1038/ncomms14529 (2017).

use Getopt::Long;
use Math::CDF;
use strict;

$|=1;

# This script scores rhythmicity in young and old as the -log of the Arser p-value. 
# Thus, rhythmicity fold changes are the differences in -log(p-value) between two conditions. 
# It calculates p-values from z-scores using the Perl module Math::CDF
# and uses the Benjamini-Hochberg procedure to adjust the p-value and produce q-values.

my $usage = "
Usage:\n\n
      $0 <DR score file> [options]\n\n
      -f --fdr           a false discovery rate (FDR) for the BH correction (default=0.05)
      -m --mean          a value for the mean of the Gaussian fit.
      -v --variance      a value for the variance of the Gaussian fit. 
      -o --outBase       a suffix for the output file (default=LLC_stats)
      -h --help          print this help message
";


my $outBase = 'LLC_stats';
my $FDR = 0.05; # Example: 0.05 for 5% false discovery rate (FDR)
my $mean = 0.0;
my $variance = 1.0; 
my $help = 0;

my $SDRfile = $ARGV[0] or die $usage; # Output file of compute_DR_scores.pl

Getopt::Long::Configure("no_ignore_case");
# read in the options
GetOptions( 'fdr=f' => \$FDR,
	    'mean=f' => \$mean,
	    'variance=f' => \$variance,
	    'outBase=s' => \$outBase,
 	    'help' => \$help
    );

if($help) {
    die $usage;
}

my $stdDev = sqrt($variance);
my $LLC_outFile = $outBase . '_FDR' . $FDR . '_mean' . $mean . '_stdDev' . $stdDev . '.txt';

my($rhythmVals) = readSDRfile($SDRfile);
my($LLCstats,$ELCstats) = getPVals($rhythmVals,$mean,$stdDev); # ELC stats is not used in this version
my($sigInfo_LLCs) = BH_test($LLCstats,$FDR);
printSDRinfo($sigInfo_LLCs,$LLC_outFile);

sub getPVals {
    my($rhythmVals,$mean,$stdDev) = @_;
    my @LLCstats;
    my @ELCstats;
    foreach my $geneId (keys %{$rhythmVals}) {
	my($name,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$DR_score) = @{$rhythmVals->{$geneId}};
	my $zScore = ($DR_score - $mean)/$stdDev;
	my $LLCpVal = computePValue($zScore);
	my $ELCpVal = 1-$LLCpVal;
	push(@LLCstats,[$LLCpVal,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$DR_score,$geneId,$name]);
	push(@ELCstats,[$ELCpVal,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$DR_score,$geneId,$name]);
    }
    my @sortedLLCstats = sort {$a->[0] <=> $b->[0]} @LLCstats;
    my @sortedELCstats = sort {$a->[0] <=> $b->[0]} @ELCstats;
    return(\@sortedLLCstats,\@sortedELCstats);
}

sub printSDRinfo {
    my($sigInfo,$outFile) = @_;
    open(OUT,">$outFile") or die "Could not open $outFile for writing.\n";
    print OUT "#GeneID\tSymbol\tARSER Young p-value\tARSER Old p-value\tZp (periodicity dif z-score)\tYoung Max-Min Exp\tOld Max-Min Exp\tZr (robustness dif z-score)\tDR Score\tZ-score\tP-value\tQ-value\tSignificant?\tRank\tOld Exp FC\n";
    foreach my $gene (@{$sigInfo}) {
	my($pValue,$rank,$qValue,$zScore,$yPval,$oPval,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$yExpFC,$oExpFC,$DR_score,$sigStatus,$geneId,$name) = @{$gene};
	print OUT "$geneId\t$name\t$yPval\t$oPval\t$normRhythmDif\t$yExpDif\t$oExpDif\t$normLogExpDifFC\t$DR_score\t$zScore\t$pValue\t$qValue\t$sigStatus\t$rank\t$oExpFC\n";
    }
    close(OUT);
}

sub computePValue {
    my($zScore) = @_;
    return 1 - Math::CDF::pnorm($zScore);
}

sub BH_test {
    my($pInfo,$FDR) = @_;
    my @rankInfo;
    my @sigInfo;
    my $i = 0;
    my $n = scalar(@{$pInfo});
#    print "Number of genes in set: $n\n";
    my $k = 0;
    foreach my $gene (@{$pInfo}) {
	my($p,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$DR_score,$geneId,$name) = @{$gene};
        $i++;
	# q-tilde, not real q-value
        my $q = ($p * $n)/$i; 
        if($p <= $FDR * ($i/$n)) {
	    # $k will store largest rank that satistifies inequality
            $k = $i;
        }
        push(@rankInfo,[$p,$i,$q,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$DR_score,$geneId,$name]); # where $i is rank of $p
    }
    my $N = @rankInfo;
    my @list;
    my %qValues;
    for(my $r=$N-1; $r>=0; $r--) {
        my($p,$i,$q,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$DR_score,$geneId,$name) = @{$rankInfo[$r]};
        push(@list,$q);
        @list = sort {$a <=> $b} @list;
        $qValues{$geneId} = $list[0];
    }
    foreach my $gene (@rankInfo) {
        my($p,$i,$q,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$DR_score,$geneId,$name) = @{$gene};
        my $sigStatus;
        if($i <= $k) {
            $sigStatus = 'yes';
        } else {
            $sigStatus = 'no';
        }
	# the real q-value
        my $qValue = $qValues{$geneId};
	my $rank = $i;
        push(@sigInfo,[$p,$rank,$qValue,$zScore,$yPval,$oPval,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$yExpFC,$oExpFC,$DR_score,$sigStatus,$geneId,$name]);
    }
    return(\@sigInfo);
}

sub readSDRfile {
    my($file) = @_;
    my %info;
    open(FILE,$file) or die "Could not open $file\n";
    while(<FILE>) {
        chomp;
	unless(/^\#/) {
	    my($id,$name,$yExpFC,$oExpFC,$yPval,$oPval,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$DR_score) = split(/\t/);
	    $info{$id} = [$name,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$DR_score];
	}
    }
    close(FILE);
    return(\%info);
}


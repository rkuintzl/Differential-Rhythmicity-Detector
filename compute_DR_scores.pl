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

# This script is one in a series of scripts for characterizing differential rhythmicity in:
# gene expression analysis as described in
# Kuintzle, R. C. et al. Circadian deep sequencing reveals stress-response genes that adopt robust rhythmic expression during aging. Nat. Commun. 8, 14529 doi: 10.1038/ncomms14529 (2017).

use strict;
use Getopt::Long;
use Statistics::Basic;
$|=1;

# Imposed restrictions:
#     Median expression value >= 1 FPKM in young OR old.
#     Young and old p-values cannot both be 1--these genes are excluded.
#     Each gene must have at least one non-zero FPKM value in young AND old
#     (to avoid using a pseudo-count when computing the fold change)
#
# Differential rhythmicity score (S_DR) is defined as: (Zp + Zr)/2,
# where Zp is the periodicity difference z-score, and Zr is the robustness difference z-score.
# For more details, see our publication referenced above.


my $usage = "
Usage:\n\n
      $0 <ARSER young output> <ARSER old output> [options]\n\n
      -m --minFpkm       a cutoff for the median expression (default=1).
      -o --outBase       a suffix for the output file (default=DR_scores_medianFpkmCutoff)
      -h --help          print this help message
";

my $youngFile = $ARGV[0] or die $usage;
my $oldFile = $ARGV[1] or die $usage;
my $minFpkm = 1; # This is a cutoff for median expression.
my $outBase = "";
my $help = 0;

Getopt::Long::Configure("no_ignore_case");
# read in the options
GetOptions( 'minFpkm=f' => \$minFpkm,
	    'outBase=s' => \$outBase,
 	    'help' => \$help
    );

if($help) {
    die $usage;
}

my $outFile = "DR_scores_medianFpkmCutoff" . $minFpkm . '.txt';
$outFile = $outBase . "_" . $outFile if($outBase);

my($yVals,$oVals,$scores,$stats) = getRhythmInfo($youngFile,$oldFile);
printRhythmInfo($yVals,$oVals,$scores,$stats);

sub getRhythmInfo {
    my($youngFile,$oldFile) = @_;
    my %scores;
    my @rhythmDifs;
    my @logExpDifFCs;
    my $expressions = [];
    my $yVals = {};
    my $oVals = {};
    ($yVals,$expressions) = readFile($youngFile,$expressions);
    ($oVals,$expressions) = readFile($oldFile,$expressions);
    my $expPsCount = nonZeroMin($expressions);
    foreach my $id (keys %{$oVals}) { 
	my($yPval,$yAmp,$yExpDif,$yMin,$yMed,$yMax,$yRhythmScore,$yName) = @{$yVals->{$id}};
        my($oPval,$oAmp,$oExpDif,$oMin,$oMed,$oMax,$oRhythmScore,$oName) = @{$oVals->{$id}};
	$yVals->{$id}->[2] = $yExpDif;
	$oVals->{$id}->[2] = $oExpDif;
	my $yFC = ($yMax + $expPsCount)/($yMin + $expPsCount);
	my $oFC = ($oMax + $expPsCount)/($oMin + $expPsCount);
	unless($yPval == 1 && $oPval == 1) {
	    if($yMed >= $minFpkm || $oMed >= $minFpkm) {
		if($yMax > 0 && $oMax > 0) {
		    my $rhythmDif = $oRhythmScore-$yRhythmScore;
		    my $logExpDifFC = log($oExpDif/$yExpDif)/log(2);
		    $scores{$id} = [$rhythmDif,$logExpDifFC,$yFC,$oFC];
		    push(@rhythmDifs,$rhythmDif);
		    push(@logExpDifFCs,$logExpDifFC);
		}
	    }
	}
    }
    # Compute stats. r = robustness (expression differences). p = periodicity differences.
    my $rMean = mean(\@logExpDifFCs);
    my $rStdDev = Statistics::Basic::stddev(\@logExpDifFCs);
    my $pMean = mean(\@rhythmDifs);
    my $pStdDev = Statistics::Basic::stddev(\@rhythmDifs);
    my @stats = ($rMean,$rStdDev,$pMean,$pStdDev);
    return($yVals,$oVals,\%scores,\@stats);
}

sub printRhythmInfo {
    my($yVals,$oVals,$scores,$stats) = @_;
    my($rMean,$rStdDev,$pMean,$pStdDev) = @{$stats};
    open(OUT,">$outFile") or die "Could not open $outFile for writing.\n";
    print OUT "#Gene ID\tSymbol\tYoung Exp Max/Min\tOld Exp Max/Min\tARS Young p-val\tARS Old p-val\tZp (Periodicity Diff Z-score)\tYoung Exp Dif\tOld Exp Dif\tZr (Log2 Expression Diff FC Z-score)\tDR Score\n";
    foreach my $id (keys %{$oVals}) { 
	my($yPval,$yAmp,$yExpDif,$yMin,$yMed,$yMax,$yRhythmScore,$yName) = @{$yVals->{$id}};
        my($oPval,$oAmp,$oExpDif,$oMin,$oMed,$oMax,$oRhythmScore,$oName) = @{$oVals->{$id}};
	unless($yPval == 1 && $oPval == 1) {
	    if($yMed >= $minFpkm || $oMed >= $minFpkm) {
		if($yMax > 0 && $oMax > 0) { # This requirement is imposed so that a pseudo-count is not required.
		    my($rhythmDif,$logExpDifFC,$yExpFC,$oExpFC) = @{$scores->{$id}};
		    # Compute Zr (robustness z-score):
		    my $Zr = ($logExpDifFC - $rMean)/$rStdDev;
		    # Compute Zp (periodicity z-score):
		    my $Zp = ($rhythmDif - $pMean)/$pStdDev;
		    my $DR_score = ($Zr + $Zp) / sqrt(2);
		    print OUT "$id\t$yName\t$yExpFC\t$oExpFC\t$yPval\t$oPval\t$Zp\t$yExpDif\t$oExpDif\t$Zr\t$DR_score\n";
		}
	    }
	}
    }
    close(OUT);
}

sub readFile {
    my($file,$expVals) = @_;
    my %sigInfo;
    open(FILE,$file) or die "Could not open $file\n";
    while(<FILE>) {
        chomp;
	unless(/^Gene/) {
	    my @terms = split(/\t/);
	    my @expressions = @terms[14..25];
	    push(@{$expVals},@expressions);
	    my($name,$id,$filterType,$arMethod,$periodNumber,$period,$amplitude,$phase,$mean,$Rsquared,$adjustedRsquared,$coefVar,$pValue,$qValue) = split(/\t/); # $qValue is fdr_bh
	    $amplitude = 0 if($amplitude eq 'NA');
	    my $maxExp = max(\@expressions);
	    my $minExp = min(\@expressions);
	    my $medExp = median(\@expressions);
	    my $expDif = $maxExp-$minExp;
	    my $rhythmScore = -log($pValue);
	    $sigInfo{$id} = [$pValue,$amplitude,$expDif,$minExp,$medExp,$maxExp,$rhythmScore,$name];
	}
    }
    close(FILE);
    return(\%sigInfo,$expVals);
}

sub median {
    my($array) = @_;
    my @array = sort {$a <=> $b} @{$array};
    return $array[int(@array/2)] if(@array % 2);
    return ($array[int(@array/2)-1] + $array[int(@array/2)])/2;
}

sub max {
    my($array) = @_;
    my $max = 0;
    foreach my $item (@{$array}) {
        if($item>$max) {
            $max = $item;
        }
    }
    return $max;
}

sub min {
    my($array) = @_;
    my $min = 'Inf';
    foreach my $item (@{$array}) {
        if($item<$min) {
            $min = $item;
        }
    }
    return $min;
}

sub nonZeroMin {
    my($array) = @_;
    my $min = 'Inf';
    foreach my $item (@{$array}) {
        if($item<$min) {
	    unless($item == 0) {
		$min = $item;
	    }
	}
    }
    return $min;
}

sub mean { 
    my($array) = @_;
    my $count = scalar(@{$array});
    return sum($array)/$count;
}

sub sum {
    my($array) = @_;
    my $sum = 0;
    foreach my $item (@{$array}) {
        $sum += $item;
    }
    return $sum;
}


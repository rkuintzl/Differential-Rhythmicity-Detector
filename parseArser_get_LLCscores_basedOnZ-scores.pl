#!/usr/bin/perl -w

# Created on [2015]
#
# Authors:
# Rachael C. Kuintzle (rkuintzl@caltech.edu), California Institute of Technology
# David A. Hendrix (David.Hendrix@oregonstate.edu), Oregon State University
#
# This script is one in a series of scripts for characterizing differential rhythmicity in
# gene expression analysis as described in
# Kuintzle, R. C. et al. Circadian deep sequencing reveals stress-response genes that adopt robust rhythmic expression during aging. Nat. Commun. 8, 14529 doi: 10.1038/ncomms14529 (2017).

use strict;
use Scalar::Util qw(looks_like_number);
use Statistics::Basic;
use Math::CDF;

$|=1;

# Restrictions:
#     MEDIAN expression value >= 1 FPKM in young OR old.
#     Young and old p-values cannot both be 1--these genes are excluded.
#     Each gene must have at least one non-zero FPKM value in young AND old
#     (to avoid using a pseudo-count when computing the fold change)

# Differential rhythmicity score (S_DR) is defined as: (Zp + Zr)/2,
# where Zp is the periodicity difference z-score, and Zr is the robustness difference z-score.
# For more details, see our publication referenced above.
#
# Please use [./fitGaussian_DAH_script.pl -h] to see the help screen for further instructions on running this script.

my $usage = "Usage:\n$0 <ARS Cycle young output> <ARS Cycle old output>\n";

my $youngFile = $ARGV[0] or die $usage;
my $oldFile = $ARGV[1] or die $usage;

my $minFpkm = 1; # This is a cutoff for median expression.

my $outFile = 'Arser_oldVsYoung_SDR_scores_medianFpkmCutoff' . $minFpkm . '.txt';

my($yVals,$oVals,$scores,$stats) = getRhythmInfo($youngFile,$oldFile);
printRhythmInfo($yVals,$oVals,$scores,$stats);

sub getRhythmInfo {
    my($youngFile,$oldFile) = @_;
    my %scores;
    my @rhythmDifs;
    my @logExpDifFCs;
    my $expDiffs = [];
    my $expressions = [];
    my($yVals,$expDiffs,$expressions) = readFile($youngFile,$expDiffs,$expressions);
    my($oVals,$expDiffs,$expressions) = readFile($oldFile,$expDiffs,$expressions);
    my $pseudoCount = nonZeroMin($expDiffs);
    my $expPsCount = nonZeroMin($expressions);
    print "Genes with zero expression in young or old across all time points:\n";
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
		} else {
		    print "$yName\t$id\n";
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
    my($file,$expDiffs,$expVals) = @_;
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
	    push(@{$expDiffs},$expDif);
	    my $rhythmScore = -log($pValue);
	    $sigInfo{$id} = [$pValue,$amplitude,$expDif,$minExp,$medExp,$maxExp,$rhythmScore,$name];
	}
    }
    close(FILE);
    return(\%sigInfo,$expDiffs,$expVals);
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

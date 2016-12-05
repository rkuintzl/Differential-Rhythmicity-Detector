#!/usr/bin/perl -w
use strict;
use Math::CDF;

$|=1;

# This script scores rhythmicity in young and old as the suprisal, or the -log of the Arser p-value. 
# Thus, rhythmicity fold changes are the difference of the suprisal. 
# It calculates p-values from z-scores using Math::CDF
# and adjust the p-value to produce q-values as the findAntiCorrelatedGenes script does.

# This script does not pre-filter on fold change or p-value.

my $usage = "Usage:\n$0 <Arser LLC score file> <FDR> <mean> <stddev> <Arser p-value cutoff (old flies)> <max/min requirement (1.5)>\n";

my $LLCfile = $ARGV[0] or die $usage; # Created by ~/Scripts/parseArser_get_LLCscores_basedOnZscores, etc.
my $FDR = $ARGV[1] or die $usage;
# Determined from a normal fit to the 
# LLC score distribution 
# by Dave's script ~/Scripts/fitGaussian_DAH_script.pl
my $mean = $ARGV[2] or die $usage; 
my $stdDev = $ARGV[3] or die $usage; 
my $cutoff = $ARGV[4] or die $usage;
my $FCcutoff = $ARGV[5] or die $usage;

my $LLC_outFile = 'LLC_stats_FDR'. $FDR . '_ARS.Old_pValCutoff' . $cutoff . '_minFC' . $FCcutoff . '_mean' . $mean . '_stdDev' . $stdDev . '.txt';
my $ELC_outFile = 'ELC_stats_FDR'. $FDR . '_ARS.Old_pValCutoff' . $cutoff . '_minFC' . $FCcutoff . '_mean' . $mean . '_stdDev' . $stdDev . '.txt';

my($rhythmVals) = readLLCfile($LLCfile);
my($LLCstats,$ELCstats) = getPVals($rhythmVals,$mean,$stdDev);
my($sigInfo_LLCs) = BH_test($LLCstats,$FDR);
my($sigInfo_ELCs) = BH_test($ELCstats,$FDR);
printLLCfile($sigInfo_LLCs,$FCcutoff,$cutoff,$LLC_outFile);
printELCfile($sigInfo_ELCs,$FCcutoff,$cutoff,$ELC_outFile);

sub getPVals {
    my($rhythmVals,$mean,$stdDev) = @_;
    my @LLCstats;
    my @ELCstats;
    foreach my $geneId (keys %{$rhythmVals}) {
	my($name,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$LLCscore) = @{$rhythmVals->{$geneId}};
	my $zScore = ($LLCscore - $mean)/$stdDev;
	my $LLCpVal = computePValue($zScore);
	my $ELCpVal = 1-$LLCpVal;
	push(@LLCstats,[$LLCpVal,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$LLCscore,$geneId,$name]);
	push(@ELCstats,[$ELCpVal,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$LLCscore,$geneId,$name]);
    }
    my @sortedLLCstats = sort {$a->[0] <=> $b->[0]} @LLCstats;
    my @sortedELCstats = sort {$a->[0] <=> $b->[0]} @ELCstats;
    return(\@sortedLLCstats,\@sortedELCstats);
}

sub printLLCfile {
    my($sigInfo,$FCcutoff,$cutoff,$outFile) = @_;
    open(OUT,">$outFile") or die "Could not open $outFile for writing.\n";
    print OUT "#GeneID\tSymbol\tArser Young p-value\tArser Old p-value\tZp (periodicity dif z-score)\tYoung Max-Min Exp\tOld Max-Min Exp\tZr (robustness dif z-score)\tLLC Score\tZ-score\tP-value\tQ-value\tSignificant?\tRank\tOld Exp FC\n";
    foreach my $gene (@{$sigInfo}) {
	my($pValue,$rank,$qValue,$zScore,$yPval,$oPval,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$yExpFC,$oExpFC,$LLCscore,$sigStatus,$geneId,$name) = @{$gene};
	if($oPval < $cutoff && $oExpFC >= $FCcutoff) { 
	    print OUT "$geneId\t$name\t$yPval\t$oPval\t$normRhythmDif\t$yExpDif\t$oExpDif\t$normLogExpDifFC\t$LLCscore\t$zScore\t$pValue\t$qValue\t$sigStatus\t$rank\t$oExpFC\n";
	}
    }
    close(OUT);
}

sub printELCfile {
    my($sigInfo,$FCcutoff,$cutoff,$outFile) = @_;
    open(OUT,">$outFile") or die "Could not open $outFile for writing.\n";
    print OUT "#GeneID\tSymbol\tArser Young p-value\tArser Old p-value\tZp (periodicity dif z-score)\tYoung Max-Min Exp\tOld Max-Min Exp\tZr (robustness dif z-score)\tLLC Score\tZ-score\tP-value\tQ-value\tSignificant?\tRank\tOld Exp FC\n";
    foreach my $gene (@{$sigInfo}) {
        my($pValue,$rank,$qValue,$zScore,$yPval,$oPval,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$yExpFC,$oExpFC,$LLCscore,$sigStatus,$geneId,$name) = @{$gene};
        if($yPval < $cutoff && $yExpFC >= $FCcutoff) {
            print OUT "$geneId\t$name\t$yPval\t$oPval\t$normRhythmDif\t$yExpDif\t$oExpDif\t$normLogExpDifFC\t$LLCscore\t$zScore\t$pValue\t$qValue\t$sigStatus\t$rank\t$oExpFC\n";
        }
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
    print "Number of genes in set: $n\n";
    my $k = 0;
    foreach my $gene (@{$pInfo}) {
	my($p,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$LLCscore,$geneId,$name) = @{$gene};
        $i++;
	# q-tilde, not real q-value
        my $q = ($p * $n)/$i;
        if($p <= $FDR * ($i/$n)) {
	    # $k will store largest rank that satistifies inequality
            $k = $i;
        }
        push(@rankInfo,[$p,$i,$q,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$LLCscore,$geneId,$name]); # where $i is rank of $p
    }
    print "k = $k\n";
    my $N = @rankInfo;
    my @list;
    my %qValues;
    for(my $r=$N-1; $r>=0; $r--) {
        my($p,$i,$q,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$LLCscore,$geneId,$name) = @{$rankInfo[$r]};
        push(@list,$q);
        @list = sort {$a <=> $b} @list;
        $qValues{$geneId} = $list[0];
    }
    foreach my $gene (@rankInfo) {
        my($p,$i,$q,$zScore,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$LLCscore,$geneId,$name) = @{$gene};
        my $sigStatus;
        if($i <= $k) {
            $sigStatus = 'yes';
        } else {
            $sigStatus = 'no';
        }
	# the real q-value
        my $qValue = $qValues{$geneId};
	my $rank = $i;
        push(@sigInfo,[$p,$rank,$qValue,$zScore,$yPval,$oPval,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$yExpFC,$oExpFC,$LLCscore,$sigStatus,$geneId,$name]);
    }
    return(\@sigInfo);
}

sub readLLCfile {
    my($file) = @_;
    my %info;
    open(FILE,$file) or die "Could not open $file\n";
    while(<FILE>) {
        chomp;
	unless(/^\#/) {
	    my($id,$name,$yExpFC,$oExpFC,$yPval,$oPval,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$LLCscore) = split(/\t/);
	    $info{$id} = [$name,$yPval,$oPval,$yExpFC,$oExpFC,$normRhythmDif,$yExpDif,$oExpDif,$normLogExpDifFC,$LLCscore];
	}
    }
    close(FILE);
    return(\%info);
}


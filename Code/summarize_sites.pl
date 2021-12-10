 #!/usr/bin/env perl
 
 ### Author : Mathieu Flamand - Duke University
 ### Version 1.21
 ### date of last modification : 2021-12-3
 
use v5.26;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

# This script will find sites in multiple replicates and output the number of replicates and average edit

my $min=1;
my $ratio=1.5;
my $minMut=2;
my $add_flag;
my $minEdit = 0;
my $minCov = 0;
my $maxEdit = 1;

GetOptions ("r|MinRep:i"=>\$min,
            "fold=f"=>\$ratio,
            "repOnly"=>\$add_flag,
            "mut:i"=>\$minMut,
            "edit|percent|p:f"=>\$minEdit,
            "MaxEdit:f"=>\$maxEdit,
			"cov|coverage:i"=>\$minCov,
			)or die "no option were given";

my %h_site;

unless (@ARGV){
say "please provide at least 2 files to summarize";
say "Options:";
say "\t--MinRep  \# minimum number of replicates to keep sites";
say "\t--fold 1.5 \# minimum fold over control to keep sites";
say "\t--repOnly \# use when merging the same file against different controls.";
say "\t--mut \# minimum number of mutation in each replicate to keep sites.";
say "\t--edit \# minimum editing rate to keep sites (as of fraction of 1, eg 0.05.";
say "\t--coverage \# minimum coverage to keep site";
exit(1);

}

while(<>)
{
    next if $_ =~ /^\#/;
    chomp;
	my ($chr, $bp_m1, $bp, $ID, $dart_edit_ratio_out, $strand, $control_edit_ratio_out, $control_total, $dart_edit_ratio_out1, $dart_total, $experiment,$score) = split("\t",$_) ;
    my ($gene,$type,$n_mut,$foldover)=(split(/\|/, $ID))[0,1,3,4];
    
    $n_mut=(split("=", $n_mut))[1];

    $foldover = $ratio if ($foldover eq "NA");

    next if ($foldover < $ratio or $n_mut < $minMut);
	next if ($dart_total < $minCov);
    next if ($dart_edit_ratio_out < $minEdit );
    next if ($dart_edit_ratio_out > $maxEdit );

    if (! defined $h_site{$chr}{$bp_m1})
    {
        $h_site{$chr}{$bp_m1}={'gene'=>$gene, 'strand'=>$strand, 'experiment'=>$experiment,'type'=>$type};
        $h_site{$chr}{$bp_m1}{'n_mut'}=$n_mut;
        $h_site{$chr}{$bp_m1}{'n'}=1;
        # $h_site{$chr}{$bp_m1}{'fold'}=$foldover;
        $h_site{$chr}{$bp_m1}{econtrol}=$control_edit_ratio_out;
        $h_site{$chr}{$bp_m1}{edart}=$dart_edit_ratio_out;
        $h_site{$chr}{$bp_m1}{control_total}=$control_total;
        $h_site{$chr}{$bp_m1}{dart_total}=$dart_total;
    }
    else
    {
        $h_site{$chr}{$bp_m1}{'n_mut'}+=$n_mut;
        # $h_site{$chr}{$bp_m1}{'fold'}+=$foldover;
        $h_site{$chr}{$bp_m1}{'n'}++;
        $h_site{$chr}{$bp_m1}{edart}+=$dart_edit_ratio_out;
        $h_site{$chr}{$bp_m1}{dart_total}+=$dart_total;
        $h_site{$chr}{$bp_m1}{econtrol}+=$control_edit_ratio_out;
        $h_site{$chr}{$bp_m1}{control_total}+=$control_total;
    }
}

foreach my $chr (sort keys %h_site)
{
    foreach my $pos (sort {$a <=> $b} keys %{$h_site{$chr}})
    {
		my $hr_site = $h_site{$chr}{$pos}; #hash_ref
        
        my $num = $hr_site->{'n'};
        next unless ($num >= $min);

        # my $foldover=$hr_site->{'fold'}/$d_num;

        my $dart_edit_ratio_out= sprintf("%.4f", $hr_site->{edart}/$num);
        my $control_edit_ratio_out= sprintf("%.4f", $hr_site->{econtrol}/$num);
        
        my $new_ratio= "NA";
        if ($hr_site->{econtrol} > 0) {$new_ratio = sprintf("%.4f", $hr_site->{edart}/$hr_site->{econtrol});}
        my $gene = $hr_site->{gene};
        my $experiment= $hr_site->{experiment};
        my $n_mut= $hr_site->{'n_mut'};
        $n_mut = $n_mut / $num if $add_flag;
        my $type= $hr_site->{type};
        my $strand= $hr_site->{strand};
        my $control_total= $hr_site->{control_total};
        my $dart_total=$hr_site->{dart_total};
        $dart_total/=$num if $add_flag;
        my  $ID = $gene . "|". $type ."|" . $experiment . "|mut=" . $n_mut . "|". $new_ratio ."|rep=". $num;
        
        say join("\t", $chr, $pos, $pos+1, $ID, $dart_edit_ratio_out, $strand, $control_edit_ratio_out, $control_total, $dart_edit_ratio_out, $dart_total, $experiment,$num);

    }
}
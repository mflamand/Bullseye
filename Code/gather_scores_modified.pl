#!/bin/env perl

#V1.1; Mathieu Flamand, Université Laval

use v5.26;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
#use Data::Dumper

my $out_score;
my $coverage;
my $mutations;
my $outfile;

GetOptions ("score"=>\$out_score,
			"coverage|cov"=>\$coverage,
			"mutations|mut"=>\$mutations,
			"outfile|o:s"=>\$outfile,
		) or die "Error parsing command line options\n";

unless ($out_score or $coverage or $mutations)
{
	say "Please use at least one of the following options: score, coverage, mutations";
	exit(0);
}

if (! $outfile and (($out_score and $coverage) and ($out_score and $mutations) or ($coverage and $mutations) ))
{
	say "Please provide an output file name with : -outfile if more than one outfile is requested";
	exit(0);
}

unless (@ARGV) {
    die "Please provide one or more input BED files\n";
}

my $score_fh;
my $cov_fh;
my $mut_fh;

if ($outfile)
{
	$outfile =~ /(.+)\.[^.]+$/;
	my $out_stem = $1; # removes extension from outfile
	my $out_dir = dirname($outfile);
	mkdir($out_dir) unless (-d $out_dir);
	### print header of output file

	if ($out_score)
	{
		my $out_file= "$out_stem".".score.txt";
		open($score_fh,">", $out_file) or die "Can't open $out_file for writing: $!\n";
	}
	if ($coverage)
	{
		my $out_file= "$out_stem".".coverage.txt";
		open($cov_fh,">", $out_file) or die "Can't open $out_file for writing: $!\n";
	}
	if ($mutations)
	{
		my $out_file= "$out_stem".".mut.txt";
		open($mut_fh,">", $out_file) or die "Can't open $out_file for writing: $!\n";
	}
}
else
{ ## if no output file given
	$score_fh = \*STDOUT;
	$cov_fh = \*STDOUT;
	$mut_fh = \*STDOUT;
}

my %hash;
my @header=("\#chr", "start", "end", "cluster","strand");
my @index;

foreach my $file (@ARGV)
{
	my $name=basename($file, ".bed");
	push @header, $name;
	push @index, $name;

	open(my $fh, "<", $file);
	
	while(<$fh>)
	{
		chomp;
		next if ($_=~/^[\#\n]/);
		my ($chr,$start,$end,$cluster,$score,$strand, $n_C, $n_E, $e, $cov, $type) = (split(/\t/,$_))[0..10];
		my $pos = join("_",$chr, $start,"$end");

		my $mut= int($cov * $e + 0.5);

		if (not defined $hash{$pos})
		{
			$hash{$pos}={'CHR'=>$chr, 'START'=>$start, 'END'=>$end, 'STRAND'=>$strand, 'CLUSTER'=>$cluster};
		}
		$score = 'NA' if $cov == 0; # in case there is not enough coverage
		$hash{$pos}{$name}{'SCORE'}= $score ;
		$hash{$pos}{$name}{'COV'}= $cov ;
		$hash{$pos}{$name}{'MUT'}= $mut ;
	}
	close $fh;
}

if ($out_score)
{
	say $score_fh join("\t", @header);
}
if ($coverage)
{
	say $cov_fh join("\t", @header);
}
if ($mutations)
{
	say $mut_fh join("\t", @header);
}

foreach my $pos (sort { $hash{$a}{'CHR'} cmp $hash{$b}{'CHR'} || $hash{$a}{'START'} <=> $hash{$b}{'START'} } keys %hash)
{
	# my $hash_ref = %{$hash{$cluster}};
	my $chr = $hash{$pos}{'CHR'};
	my $start = $hash{$pos}{'START'};
	my $end = $hash{$pos}{'END'};
	my $strand = $hash{$pos}{'STRAND'};
	my $cluster =   $hash{$pos}{'CLUSTER'};

	my @score;
    my @cov;
    my @mut;

	foreach my $sample (@index) {

		if (not defined $hash{$pos}{$sample})
		{
			$hash{$pos}{$sample}={'SCORE'=>'NA', 'COV'=>0, 'MUT'=>0};
		}
		
		push @score, $hash{$pos}{$sample}{SCORE};
        push @cov,   $hash{$pos}{$sample}{COV};
        push @mut,   $hash{$pos}{$sample}{MUT};
    }

	if ($out_score)
	{
		say $score_fh join("\t", $chr, $start,$end,$cluster,$strand,@score);
	}
	if ($coverage)
	{
		say $cov_fh join("\t", $chr, $start,$end,$cluster,$strand,@cov);
	}
	if ($mutations)
	{
		say $mut_fh join("\t", $chr, $start,$end,$cluster,$strand,@mut);
	}

}

#print Dumper(\%hash);


close $score_fh if $score_fh;
close $cov_fh if $score_fh;
close $mut_fh if $mut_fh;
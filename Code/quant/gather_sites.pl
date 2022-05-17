#!/bin/env perl

use v5.26;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $out_score;
my $coverage;
my $mutations;
my $outfile;

GetOptions ("score"=>\$out_score,
			"coverage|cov"=>\$coverage,
			"mutations|mut"=>\$mutations,
			"outfile|o:s"=>\$outfile,
		) or error_out();

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

{

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
	foreach my $file (@ARGV)
	{
		my $name=basename($file, ".bed");
		push @header, $name;
		
		open(my $fh, "<", $file);
		
		while(<$fh>)
		{
			chomp;
			next if ($_=~/^[\#\n]/);
			my ($chr,$start,$end,$cluster,$score,$strand, $n_C, $n_E, $e, $nc, $cov) = (split(/\t/,$_))[0..10];
			my $pos = join("_",$chr, $start,"$end");
			if (not defined $hash{$pos})
			{
				$hash{$pos}={'CHR'=>$chr, 'START'=>$start, 'END'=>$end, 'STRAND'=>$strand, 'CLUSTER'=>$cluster};
			}
			$score = 'NA' if $cov == 0; # in case there is not enough coverage
			push( @{$hash{$pos}{'SCORE'}}, $score );
			push( @{$hash{$pos}{'COV'}}, $cov );
			push( @{$hash{$pos}{'MUT'}}, $e );
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

		if ($out_score)
		{
			my @score= @{$hash{$pos}{'SCORE'}};
			say $score_fh join("\t", $chr, $start,$end,$cluster,$strand,@score);
		}
		if ($coverage)
		{
			my @cov= @{$hash{$pos}{'COV'}};
			say $cov_fh join("\t", $chr, $start,$end,$cluster,$strand,@cov);
		}
		if ($mutations)
		{
			my @mut= @{$hash{$pos}{'MUT'}};
			say $mut_fh join("\t", $chr, $start,$end,$cluster,$strand,@mut);
		}

	}
}
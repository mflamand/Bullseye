#!/usr/bin/env perl

 ### Author : Mathieu Flamand - Duke University
 ### version : 0.2 
 ### date of last modification : 2022-10-13

 ### This short program will merge matrix files 

use strict;
use warnings;
use v5.26;
use Getopt::Long;

my ($file,$outfile);

GetOptions ("i|input:s"=>\$file,
			"o|output:s"=>\$outfile,
		) or die ;
my $fh;
if ($file =~ /\.gz$/){
	open($fh, "zcat $file |") or die "can't open file $file for reading\n";
}else{
	open ($fh, "<", $file) or die "can't open file $file for reading\n";
}
open (my $fout, ">", $outfile) or die "can't open file $outfile for writing\n";

my $count_hash = {};
my $prev_chr='';
my $prev_coord='';
while(my $line = <$fh>)
{ 
	chomp $line;
	my($chr,$coord, $A, $T, $C, $G,$N, $total) = (split(/\t/,$line))[0,1,2,3,4,5,6,7];
	
	$count_hash->{$chr}->{$coord}->{A} += $A;
	$count_hash->{$chr}->{$coord}->{T} += $T;
	$count_hash->{$chr}->{$coord}->{C} += $C;
	$count_hash->{$chr}->{$coord}->{G} += $G;
	$count_hash->{$chr}->{$coord}->{N} += $N;
	$count_hash->{$chr}->{$coord}->{total} += $total;

	
	if ($prev_chr ne $chr and $prev_chr ne '')
	{
		print_line($prev_chr, $prev_coord, $count_hash)
	}
	
	if ($prev_coord ne $coord and $prev_coord ne '' and $prev_chr eq $chr)
	{
		print_line($chr, $prev_coord, $count_hash);
	}
					
			
	$prev_chr = $chr;
	$prev_coord = $coord;

}

foreach my $pos (keys %{$count_hash->{$prev_chr}})
{
	print_line($prev_chr, $pos, $count_hash);
}

close($fh);
close($fout);

system("bgzip -f $outfile") == 0 or die "error compressing final file with bgzip";
system("tabix -b 2 -e 2 $outfile.gz") == 0 or die "error indexing final file with tabix";



sub print_line
{
	my $chr = shift @_;
	my $pos = shift @_;
	my $hash_ref = shift @_;
	my $A = $hash_ref->{$chr}->{$pos}->{A};
	my $T = $hash_ref->{$chr}->{$pos}->{T};
	my $C = $hash_ref->{$chr}->{$pos}->{C};
	my $G = $hash_ref->{$chr}->{$pos}->{G};
	my $N = $hash_ref->{$chr}->{$pos}->{N};
	my $total = $A + $T + $C + $G + $N;
	my  @outline = ($chr,$pos,$A,$T,$C,$G,$N,$total);
	say $fout join("\t",@outline);		
	delete $hash_ref->{$chr}->{$pos};
}


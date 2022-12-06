 #!/usr/bin/env perl
 
 ### Author : Mathieu Flamand - Duke University
 ### version : 0.2
 ### date of creation: 2022-12-06 
 
### This scripts converts a gtf file to a refFlat file for use in Bullseye
### The GTF file needs to have the following fields: gene_name transcript_id transcript_version

$| = 1;

use v5.26;
use strict;
use warnings;
use Getopt::Long;

my $USAGE = "$0 --gtf ensembl.gtf --out ensembl.refFlat \n";

my ($input,$output);

GetOptions ("i|gtf=s"=>\$input,
			"o|out=s"=>\$output,
) or error_out();


my $in_fh;
if ($input =~ /.gz$/){
   	open($in_fh, "zcat $input |") or die "Can't zcat $input for reading: $!\n";
}else{
    open($in_fh, "<", $input) or die "Can't read $input: $!\n";
}

open(my $out_fh, ">",$output) or die "Can't open file $output for writing: $1\n";

my %genes;
my $index=1;
my $split_start_flag = 0;
my $split_stop_flag = 0;
my $last_transcript = '';

while (<$in_fh>){
    next if $_ =~ /^\#/;
    my ($chr,$def,$type,$start,$end,$score,$strand,$phase,$group) = split(/\t/, $_);
    $group =~ /gene_name\s\"(.+?)\"/;
    $start--; #offset
    my $name = $1;

    if ($type eq "gene"){
            $genes{$chr}{$name}{'index'}=$index;
            $index++;
    }else{
        $group =~ /transcript_id\s\"(.+?)\"/;
        my $transcript = $1;
        $group =~ /transcript_version\s\"(\d+?)\"/;
        my $transcript_version = $1;
        my $transcript_id = $transcript.".".$transcript_version; 
        if($last_transcript ne $transcript_id and $last_transcript ne ''){
            $split_start_flag = 0;
            $split_stop_flag = 0;
        }

       

        if ($type eq "transcript"){
            if($strand eq "+"){
                    $genes{$chr}{$name}{$transcript_id}={strand=>$strand,TxStart=>$start,TxEnd=>$end,CDS_start=>$end,CDS_end=>$end,exon_count=>0};
                }else{
                    $genes{$chr}{$name}{$transcript_id}={strand=>$strand,TxStart=>$start,TxEnd=>$end,CDS_start=>$end,CDS_end=>$end,exon_count=>0};
                }
        }
        if ($type eq "exon"){
            $genes{$chr}{$name}{$transcript_id}{'exon_count'}++;
            push ( @{$genes{$chr}{$name}{$transcript_id}{'start_list'}},$start );
            push ( @{$genes{$chr}{$name}{$transcript_id}{'end_list'}},$end );
        }
        if ($type eq "start_codon"){
            if($split_start_flag==0){
                if($strand eq "+"){
                    $genes{$chr}{$name}{$transcript_id}{'CDS_start'}=$start;
                }else{
                    $genes{$chr}{$name}{$transcript_id}{'CDS_end'}=$end;
                }
            }else{
                if($strand eq "+" and $start >= $genes{$chr}{$name}{$transcript_id}{'CDS_start'}){
                    $genes{$chr}{$name}{$transcript_id}{'CDS_start'}=$start;
                }elsif($strand eq "-" and $end >= $genes{$chr}{$name}{$transcript_id}{'CDS_end'}){
                    $genes{$chr}{$name}{$transcript_id}{'CDS_end'}=$end;
                }
            }
            $split_start_flag=1;
        }
        if ($type eq "stop_codon"){
            if($split_stop_flag==0){
                if($strand eq "+"){
                    $genes{$chr}{$name}{$transcript_id}{'CDS_end'}=$end;
                }else{
                    $genes{$chr}{$name}{$transcript_id}{'CDS_start'}=$start;
                }
            }else{
                if($strand eq "+"){
                    $genes{$chr}{$name}{$transcript_id}{'CDS_end'}=$end if $end >= $genes{$chr}{$name}{$transcript_id}{'CDS_end'};
                }else{
                    $genes{$chr}{$name}{$transcript_id}{'CDS_start'}=$start if $start >= $genes{$chr}{$name}{$transcript_id}{'CDS_start'};
                }
            }
            $split_stop_flag=1;
        }
        $last_transcript=$transcript_id;
    }
    
}

close($in_fh);

foreach my $chr (sort {$a cmp $b}  keys %genes){
    foreach my $gene (sort { $genes{$chr}{$a}{'index'} <=> $genes{$chr}{$b}{'index'} } keys %{$genes{$chr}}){
        delete $genes{$chr}{$gene}{'index'};
        foreach my $transcript (keys %{$genes{$chr}{$gene}}){
            my $hash_ref=$genes{$chr}{$gene}{$transcript};
            my $strand = $hash_ref->{'strand'};
            my $TxStart = $hash_ref->{'TxStart'};
            my $TxEnd = $hash_ref->{'TxEnd'};
            my $CDS_start = $hash_ref->{'CDS_start'};
            my $CDS_end = $hash_ref->{'CDS_end'};
            my $exon_count = $hash_ref->{'exon_count'};
            my @start_list = @{$hash_ref->{'start_list'}};
            my @end_list = @{$hash_ref->{'end_list'}};

            # case when start codon is define, but not stop codon #

            if ( $strand eq "-" and ($CDS_start == $TxEnd and $CDS_end != $TxEnd) ){
                    $CDS_start = $TxStart;
            }
            if ( $strand eq "+" and ($CDS_end == $TxEnd and $CDS_start != $TxEnd) ){
                    $CDS_end = $TxEnd;
            }
            

            if ($strand eq "-"){
                @start_list=reverse(@start_list);
                @end_list=reverse(@end_list);
            }
            my $start_list =join(",", @start_list).",";
            my $end_list =join(",", @end_list).",";
            say $out_fh join("\t", $gene,$transcript,$chr,$strand,$TxStart,$TxEnd,$CDS_start,$CDS_end,$exon_count,$start_list,$end_list);
        }
    }
}

close($out_fh);

my $cmd = "LC_ALL=C sort -f -k3,3 -k5,5n -o $output $output";
system($cmd) == 0 or die "Failed to sort file $output:$!";


sub error_out {
	say STDERR "\n\n$USAGE";

print STDERR <<'USAGE::BLOCK';

complete list of options: 

	--gtf	: Input GTF file
	--out	: output genePred file (refFlat)
	
USAGE::BLOCK

	exit (1);

}

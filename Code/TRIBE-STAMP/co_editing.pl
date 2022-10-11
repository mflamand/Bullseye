#!/usr/bin/env perl
 
### Author : Mathieu Flamand - Duke University
### version : 1.2
### date of last modification : 2021-12-8
### This program looks for, and quantifies co-editing events occuring within the same reads. 

use v5.26;
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Data::Dumper;
use MCE::Loop;
use File::Basename;
use Text::NSP::Measures::2D::Fisher::twotailed;
use Array::IntSpan; 

my $USAGE = "$0 --bed1 sn1.bed --bed2 sn2.bed --bam file.bam --outfile output_file.tsv\n";


my ($bed1,$bed2,$bam,$remove_duplicate,$remove_MM,$verbose,$filter_noedit,$OUTFILE,$annotationfile,$help,$return_pairs)='';
my $n_iter =0;
my $distance = 10; # default option
my $min_cov= 10; # default option
my $ncpu= $ENV{'SLURM_CPUS_PER_TASK'} // 1; # This if for running on slurm systems, number of cores can be otherwise selected with the --cpu option
my $min_QUAL=20; #min quality of 20 (1/100)

#Get option from commmand line
GetOptions ("bed1=s"=>\$bed1,
			"bed2=s"=>\$bed2,
            "bam:s"=>\$bam,
			"d|distance:i"=>\$distance,
            "cpu:i"=>\$ncpu,
			"removeDup"=>\$remove_duplicate,
			"removeMultiMapped"=>\$remove_MM,
			"o|outfile=s"=>\$OUTFILE,
			"minCov=i"=>\$min_cov,
			"FilterNeg"=>\$filter_noedit,
			"refFlat=s"=>\$annotationfile,
			"PHRED:i"=>\$min_QUAL,
			"v|verbose"=>\$verbose,
			"h|help"=>\$help,
			"ReturnPairsOnly"=>\$return_pairs,
			"PermutationTest=i"=>\$n_iter,
		) or die("Error in command line arguments, please use --help option for usage information\n");

if ( $help ) { error_out(); }

unless(-f $bed1){say "Error, bed file $bed1 could not be found";error_out();}
unless(-f $bed2){say "Error, bed file $bed2 could not be found";error_out();}
unless(!$return_pairs and -f $bam){say "Error, bam file could not be found";error_out();}

if ($verbose and ! $n_iter){
	say STDERR "comparing pairs in $bed1 and $bed2, using Fisher's exact test";
}
elsif($verbose and $n_iter){
	say STDERR "comparing pairs in $bed1 and $bed2, using $n_iter iteration of the Chi-squared test";
}

$min_QUAL+=33; #For illumina read quality

my $out_fh;
if ($OUTFILE){
	my $out_dir = dirname($OUTFILE);
	mkdir($out_dir) unless (-d $out_dir); # create output directory if it does not exists
	open($out_fh,">", $OUTFILE) or croak "Can't open $OUTFILE for writing: $!\n";
}
else{
	$out_fh = \*STDOUT; #MCE compatible standard output
}

### identify the reverse complement for negative strand
my %complement = ("A"=>"T","T"=>"A","C"=>"G","G"=>"C"); 

my $introns = {};# build hash for intron positions
# This is for the optional mapping of reads spanning introns (only one intron supported for now).
# A RefFlat file is read to determine introns positions
if ($annotationfile){
	open(my $annotation_handle,"<", $annotationfile) or die "Can't open $annotationfile: $!\n";   
	foreach my $line (<$annotation_handle>){  
		chomp $line;
		next if ($line=~/^[\#\n]/);
		next if ($line=~/^name/); 
		my($chr,$strand,$exon_count,$sString,$eString) = (split(/\t/,$line))[2,3,8,9,10];
		my @starts = split(/\,/,$sString); 
		my @ends = split(/\,/,$eString);

		my $intron_count = $exon_count-1;

		if ($intron_count>0){
			my @intronStart;
			my @intronEnds;

			#Get start and ends of introns from exons ends and start
			for (my $i = 0; $i <= $#ends-1;$i++){
				my $j = $i +1;
				my $intron_start=$ends[$i]+1;
				my $intron_end= $starts[$j]-1;
				push @intronStart, $intron_start;
				push @intronEnds, $intron_end;
		
			}
			for (my $i = 0; $i <= $#intronStart;$i++){    
				my $legnth = $intronEnds[$i]-$intronStart[$i];
				$introns->{$chr}->{$strand}->{start}->{$intronStart[$i]}=$legnth;
				$introns->{$chr}->{$strand}->{end}->{$intronEnds[$i]}=$legnth;
					
			}
		}
	}
}

#Store sites of first Bed file (outputed from Bullseye find_edit_site.pl)
my $sites={};
open (my $dart_fh, "<", "$bed1") or croak "can't open file $bed1";		
while(my $line = <$dart_fh>){
    	chomp $line;
		my ($chr,$pos,$strand,$cov) = (split(/\t/, $line))[0,2,5,9];
		my ($type,$n_mut) = (split(/\|/, $line))[2,3];
		$n_mut=~ s/mut\=(\d+)$/$1)/;
        $sites->{$chr}->{$strand}->{$pos}=$type;
}

#open second bed file and identify all pairs to measure
open (my $trybe_fh, "<", "$bed2") or croak "can't open file $bed2";

my @pairs;
my %dup_hash;

my $i=1;
while(my $line = <$trybe_fh>){
    	chomp $line;
		my($chr,$pos,$name,$strand,$cov) = (split(/\t/, $line))[0,2,3,5,9];
        my ($type2,$n_mut);
        ($name, $type2) = (split(/\|/,$name))[0,2,3];
		
		my $start = $pos - $distance;
        my $end = $pos + $distance;

        foreach my $nuc ($start..$end){
            my $pair;
			## check for sites spanning one intron and allow pairing if a refFlat is provided
			## This is for the case when introns are on the left of probed sites

			if ($nuc < $pos and defined $introns->{$chr}->{$strand}->{end}->{$nuc}){  # if the position map to an intron end, it means that there was an intron upstream.
                
				my $length=$introns->{$chr}->{$strand}->{end}->{$nuc};
				my $start2=$start - $length; #go back and search before intron
				my $end2=$nuc-$length;             

				foreach my $nuc2 ($start2..$end2){
                    if (defined $sites->{$chr}->{$strand}->{$nuc2}){
                        my $type = $sites->{$chr}->{$strand}->{$nuc2};
						my $sitedistance = abs($nuc2-$pos)-$length;
                        $pair = "$chr:$nuc2:$type:$pos:$type2:$name:$strand:$sitedistance";
						next if defined $dup_hash{"$chr:$pos:$nuc2:$strand:$sitedistance"}; # avoid using the pair if it already exists 
                        next if defined $dup_hash{"$chr:$nuc2:$pos:$strand:$sitedistance"};
                        $dup_hash{"$chr:$nuc2:$pos:$strand:$sitedistance"}=1;

                        push @pairs, $pair;
                    }
                }
            }

			## This is for the case when introns are on the right of probed sites
			if ( $nuc > $pos and defined $introns->{$chr}->{$strand}->{start}->{$nuc}){  #case when nucleotide matched intron boundary
               
			   	my $length=$introns->{$chr}->{$strand}->{start}->{$nuc};
			   	my $start2=$nuc+$length;
                my $end2=$end+$length;
                
				foreach my $nuc2 ($start2..$end2){
                    if (defined $sites->{$chr}->{$strand}->{$nuc2}) {
                        my $type = $sites->{$chr}->{$strand}->{$nuc2};
						my $sitedistance = abs($nuc2-$pos)-$length;
                        $pair = "$chr:$nuc2:$type:$pos:$type2:$name:$strand:$sitedistance";
                        next if defined $dup_hash{"$chr:$pos:$nuc2:$strand:$sitedistance"};
                        next if defined $dup_hash{"$chr:$nuc2:$pos:$strand:$sitedistance"}; 
                        $dup_hash{"$chr:$nuc2:$pos:$strand:$sitedistance"}=1;

                        push @pairs, $pair;
                    }
                }
            }
            # with no intron
            if (defined $sites->{$chr}->{$strand}->{$nuc}) {
                next if $pos == $nuc; 
                my $type = $sites->{$chr}->{$strand}->{$nuc};
				my $sitedistance = abs($nuc-$pos);
				$pair = "$chr:$nuc:$type:$pos:$type2:$name:$strand:$sitedistance";
				next if defined $dup_hash{"$chr:$pos:$nuc:$strand:$sitedistance"}; 
				next if defined $dup_hash{"$chr:$nuc:$pos:$strand:$sitedistance"};  
				$dup_hash{"$chr:$nuc:$pos:$strand:$sitedistance"}=1;
	            push @pairs, $pair;
            }
        }
}

say STDERR $#pairs+1 ." pairs will be compared for analysis" if $verbose;
undef $sites;
undef %dup_hash;
	
if ($return_pairs) {
	say $out_fh join("\t", ('#chr', 'site1', 'site2','name','distance','Strand'));
	foreach (@pairs){
		my ($chr,$nuc,$type1,$pos,$type2,$name,$strand,$sitedistance)=split(":",$_);
		my $start = $nuc < $pos ? $nuc : $pos ;
		my $end = $nuc < $pos ? $pos : $nuc ;
		say $out_fh join("\t",$chr,$start,$end, $name, $sitedistance, $strand)
	}
	exit(0);
}

# Prepare output file header	
print $out_fh join("\t", ('#chr', 'pos1', 'pos2','Strand', 'Total1','ratio1','Total2','ratio2', 'SUM', 'NEITHER', 'SITE1', 'SITE2', 'BOTH', 'RATIO', 'Fisher_P'));
print $out_fh ("\tRand_P\t") if $n_iter; 
say $out_fh join("\t",('Distance','GeneID'));

MCE::Loop->init(
	chunk_size => 1,
	max_workers => $ncpu, ## use number of specified threads for processing 
);

mce_loop
{
    my ($mce, $chunk_ref, $chunk_id) = @_;
   	my ($chr,$nuc,$type1,$pos,$type2,$name,$strand,$sitedistance)=split(":",$_);

    my $outline = parse_bam($bam,$chr,$nuc,$type1,$pos,$type2,$strand);
		
	MCE->next if( not defined $outline or (split("\t",$outline))[6] < $min_cov);
	
	$outline =join("\t", $outline, $sitedistance, $name);
	MCE->say($out_fh,$outline);
	
}@pairs; 

sub parse_bam {
    # parse_bam($bam,$chr,$nuc,$pos,$strand);
    my $file = shift @_;
    my $chr = shift @_;
    my $site1 = shift @_;
	my $type1 = shift @_;
    my $site2 = shift @_;
    my $type2 = shift @_;
	my $strand =shift @_;

    my ($left,$right) = (0,0);
    my %data;
	my $mut_flag;
	my ($noedit_flag,$flag1,$flag2,$both_flag)=0;
	my $interval;

    $interval = $chr.":".$site1."-".$site2 if ($site1 < $site2);
    $interval = $chr.":".$site2."-".$site1 if ($site1 > $site2);
    
	return unless ($file and $interval);
	
	open( my $fin, "samtools view $file $interval|") or die "Can't read bam file"; # read indexed bam file by chromosome

  	while (my $line = <$fin>){
		chomp $line;
	
    	$mut_flag=0; #6 = no mutation, 4 = C2U mutation, 5 = A2I mutation, 3 = both mutation 
		my($flag,$start,$matchinfo,$seq,$phred) = (split(/\t/,$line))[1,3,5,9,10];
		### check for PCR/optical duplicates and remove if needed
		if ($remove_duplicate or $remove_MM){
			$flag = sprintf ("%012b", $flag);
			my ($dup,$secondary_alignement) = (split (//, $flag))[1,3];
			if ($remove_duplicate){
				next if $dup eq 1;
			}
			if ($remove_MM){
				next if $secondary_alignement eq 1;
			}
		}

        next if ($site1 < $site2 and $start > $site1);
		next if ($site1 > $site2 and $start > $site2);

		($matchinfo, $seq,$phred) = parse_cigar($matchinfo, $seq,$phred);

        my @matches= $matchinfo =~ /(\d+M)/g;
        my @skip = $matchinfo =~ /(\d+N)/g;
        my @datapts = split('',$seq); # split sequence in each nucleotide
        my @qual = unpack("C*", $phred); # split quality 
       
        foreach my $b ("A","T","C","G","N")	{  # initialize hash to 0
            $data{$site1}{$b} += 0;
            $data{$site2}{$b} += 0;
        }
        	
        #deconstruct reads in smaller parts removing introns ("N" in cigar)
        my $offset = 0;
        my $mapped_array = Array::IntSpan->new();	

        for (my $i = 0; $i <= $#matches; $i++ ){
            $matches[$i] =~ /(\d+)M/;
            my $length= $1;
            my $end = $start+$length-1;
            if ($i > 0){
                  $skip[$i-1] =~ /(\d+)N/; 
                  $offset += $1;
            }
            $mapped_array->set_range($start,$end,1); #define array for each contiguous region of read
            $offset += $length;
        }
	
        next unless (defined $mapped_array->lookup($site1) and defined $mapped_array->lookup($site2)); #if both sites are not found in read, skip it. 
	    $offset = 0;
        for (my $i = 0; $i <= $#matches; $i++ ){
            my $match = $matches[$i];
            $match =~ /(\d+)M/;
            my $length = $1;
            if ($i > 0){
                 $skip[$i-1] =~ /(\d+)N/; $offset += $1;
            }
            my @sequence = splice(@datapts,0,$length );

            for (my $j = 0; $j <= $#sequence; $j++){
                my $base = $sequence[$j];
                next if $base eq '-';	
                my $bp = $start + $offset + $j;
		
                next unless ($bp == $site1 or $bp == $site2);
                $mut_flag += return_flag($strand,$bp,$base,$site1,$site2,$type1,$type2);
	
			   	$data{$bp}{$base}++;
			
            }
            $offset += $#sequence + 1;
        }

		$mut_flag == 6 ? $noedit_flag+=1 : $noedit_flag+=0;
		$mut_flag == 4 ? $flag1+=1 : $flag1+=0;
		$mut_flag == 5 ? $flag2+=1 : $flag2+=0;
		$mut_flag == 3 ? $both_flag+=1 : $both_flag+=0;

	}

	return unless (defined $both_flag and defined $flag1 and defined $flag2 and defined $noedit_flag);
	return if $filter_noedit and ($both_flag + $flag1 == 0 or $both_flag + $flag2 == 0);
	return if ($both_flag+$flag2+$flag1+$noedit_flag <1);

    my ($type1_non_edit,$type1_edit) = return_base($strand, $type1);
    my ($type2_non_edit,$type2_edit) = return_base($strand, $type2);
      
    my $A = $data{$site1}{"A"};
    my $T = $data{$site1}{"T"};
    my $C = $data{$site1}{"C"};
    my $G = $data{$site1}{"G"};
    my $N = $data{$site1}{"N"};
    my $total1 = $A + $T + $C + $G ;

    # my $total1 = $both_flag+$flag2+$flag1+$noedit_flag;
	return if $total1 < $min_cov;
    my $edit_ratio1 = sprintf( "%.6f", ($data{$site1}{$type1_edit}/$total1));

    $A = $data{$site2}{"A"};
    $T = $data{$site2}{"T"};
    $C = $data{$site2}{"C"};
    $G = $data{$site2}{"G"};
	$N = $data{$site2}{"N"};
   	my $total2 = $A + $T + $C + $G ;
	# my $total2 = $both_flag+$flag2+$flag1+$noedit_flag;
	
	return if $total2 < $min_cov;

    my $edit_ratio2 = sprintf( "%.6f", ($data{$site2}{$type2_edit}/$total2));
    my $both_ratio = sprintf( "%.6f", ($both_flag/($both_flag+$flag2+$flag1+$noedit_flag)));
	my $Fischer = 1;
	# get two-sided Fisher exact test result
	$Fischer = calculateStatistic( n11=>$both_flag,
                                n1p=>$both_flag+$flag2,
                                np1=>$both_flag+$flag1,
                                npp=>$both_flag+$flag2+$flag1+$noedit_flag);
    
	# if using the randomization test, perform n iteration of chi square test
	my $p_val_iter=0;
    if($n_iter and (($data{$site1}{$type1_edit}+$data{$site1}{$type1_non_edit})==0 or ($data{$site2}{$type2_edit}+$data{$site2}{$type2_non_edit})==0 )){
		$p_val_iter=1;
	}elsif($n_iter){
		my $chi_res=chi_squared($both_flag,$flag1,$flag2,$noedit_flag);
		my $freq_site1= ($data{$site1}{$type1_edit}/($data{$site1}{$type1_edit}+$data{$site1}{$type1_non_edit})); #don't trim
		my $freq_site2= ($data{$site2}{$type2_edit}/($data{$site2}{$type2_edit}+$data{$site2}{$type2_non_edit}));
		# my $both_freq= ($both_flag/($both_flag+$flag2+$flag1+$noedit_flag));
		my $both_freq= $freq_site1*$freq_site2; #expected frequency instead of observed

		for (1..$n_iter){
			my ($a,$b,$c,$d)=(0,0,0,0);
		
			for(1..$total1)	{
				my $val=rand();
				if ($val <= $both_freq)	{
					$a+=1;
				}elsif($val <= $both_freq+$freq_site1){
					$b+=1;
				}elsif($val <= $both_freq+$freq_site1+$freq_site2){
					$c+=1;
				}else{
					$d+=1;
				}
			}
			my $iter_res=chi_squared($a,$b,$c,$d);
			$p_val_iter+=1 if $iter_res>= $chi_res ;

		}
		$p_val_iter = $p_val_iter/$n_iter; # probability that randomized values are more extreme than observed values 
	}

	my $sum = $noedit_flag+$flag1+$flag2+$both_flag;
	my @outline = ($chr,$site1,$site2,$strand,$total1,$edit_ratio1,$total2,$edit_ratio2, $sum, $noedit_flag, $flag1, $flag2, $both_flag, $both_ratio, $Fischer);
	push @outline, $p_val_iter if $n_iter;
    my $outline = join("\t",@outline);
    return $outline;
    
}

sub return_base {
    #return_base($strand,$type)
    my $strand = shift @_;
    my $type = shift @_;
    my ($base1,$base2) = (split("2", $type))[0,1];
   
    $base1 = "T" if $base1 eq "U";
    $base1 = "G" if $base1 eq "I";
    $base2 = "T" if $base2 eq "U";
    $base2 = "G" if $base2 eq "I";
    $base1 = $complement{$base1} if ($strand eq "-");
    $base2 = $complement{$base2} if ($strand eq "-");
    return ($base1,$base2);
}

sub return_flag {
    #return_flag($strand,$bp,$base,$site1,$site2,$type1,$type2)
	
    my $strand = shift @_;
	my $bp = shift @_;
	my $base =shift @_;
	my $site1 =shift @_;
	my $site2 = shift @_;
    my $type1= shift @_;
    my $type2 = shift @_;

    my ($type1_non_edit,$type1_edit) = (return_base($strand, $type1));
    my ($type2_non_edit,$type2_edit) = (return_base($strand, $type2));
   
	my $num=0;

	if ($bp == $site1 and $base eq $type1_edit){
		$num+=1;
	}
	elsif($bp == $site2 and $base eq $type2_edit){
		$num+=2;
	}
	elsif($bp == $site1 and $base eq $type1_non_edit){
		$num+=3;
	}
	elsif($bp == $site2 and $base eq $type2_non_edit){
		$num+=3;
	}
	else{
		$num+=6; # in this case, we have a different mutation, this read will not be counted for stats.
	}
    return $num;
}

sub chi_squared {
     my ($a,$b,$c,$d) = @_;
     return 0 if($b+$d == 0 or $a + $b==0 or $a + $c ==0 or $c+$d==0);
     my $n= $a + $b + $c + $d;
     return (($n*($a*$d - $b*$c)**2) / (($a + $b)*($c + $d)*($a + $c)*($b + $d)));
}

sub error_out {
	say STDERR "\n\n$USAGE";

print STDERR <<'USAGE::BLOCK';

co_editing.pl identifies and quantifies co-editing sites in bam files using two bed files containing the editing sites. This program is usually run after Find_edit_site.pl and requieres information about editing type in the bed file name field

complete list of options: 

	--bed1		: first bed file containing edited positions (s)
	--bed2		: second bed file containing edited positions (s)
	--bam		: bam file to be looked at(s)
	--distance	: (or -d) maximal distance between sites for co-editing analysis(i)
	--removeDup : do not count duplicate reads (recommended if files contain the proper sam flag) 
	--minCov	: minimum coverage for in for pairs of sites to be considered (i)
	--outfile	: Output bed file, defaults to STDOUT (s)
	--verbose	: display extra log information
	--refFlat	: refFlat information file to consider location of introns for and mapping of pairs across intron boundaries. Experimental, for now only looks across a single intron. (s)
	--FilterNeg	: do not print pair of sites if either of them have no editing
	--PHRED		: Do not count bases below this PHRED threshold, defaults to 20 (1/100 miscall), +33 will be added for final filtering
	--PermutationTest n : Performs a Randomization test (comparing Chi squared values) instead of a Fisher exact test for calculating probabilities of observed data over expected data. Provide number of iteration to be performed (i)  
	--help		

USAGE::BLOCK

	exit (1);

}


sub parse_cigar{
    my $matchinfo = shift @_;
    my $seq = shift @_;
	my $phred = shift @_; 

    my @datapts = split('',$seq); # split sequence in each nucleotide
	my @qual =split('',$phred);

    if($matchinfo =~ /^(\d+)S(.+?)(\d+)S$/){   ### if cigar begins and ends with soft clipped, remove softclip 
        $matchinfo = $2;
        splice(@datapts, 0, $1); 
        splice(@datapts, -$3);# remove the number of soft clipped nucleotide
		splice(@qual, 0, $1); 
        splice(@qual, -$3);
    }

    if($matchinfo =~ /(.+[M|N|I|D])(\d+)S$/) {  ### if cigar ends with soft clipped, remove it 
        $matchinfo = $1;
        splice(@datapts,-$2); 
		splice(@qual,-$2); 
    }

    if($matchinfo =~ /^(\d+)S(.+[M|N|I|D])$/) { ### if cigar begins with soft clipped, remove the clipped portion  
        $matchinfo = $2;
        splice(@datapts, 0, $1);  # remove the number of soft clipped nucleotide
		 splice(@qual, 0, $1);
    }
    $matchinfo=join_M($matchinfo);

    while ($matchinfo =~ /I/) {
        if($matchinfo =~ /(^.+?)(\d+)I(.+[M|N|I|D])$/) {  ### if cigar begins with soft clipped, remove the clipped portion  
            my $cigar_start = $1;
            my $in_size = $2;
            my $cigar_end = $3;

            my @positions= $cigar_start =~ /(\d+)[M|I]/g;
            my $offset = 0;
            foreach(@positions) {$offset += $_;}         
            
            $matchinfo = $cigar_start . $cigar_end;
            splice(@datapts, $offset, $in_size);
			splice(@qual, $offset, $in_size);
            $matchinfo=join_M($matchinfo);
        }
    }

    while ($matchinfo =~ /D/) {
        if($matchinfo =~ /(^.+?)(\d+)D(.+[M|N])$/) {
            my $cigar_start = $1;
            my $cigar_end = $3;
            my $del_length= $2;

            my @missing = map{"N"} 1..$del_length; # add 'N' for each missing base
			my @missing_phred = map{"!"} 1..$del_length; # add '-' for each missing base
            
            my @positions= $cigar_start =~ /(\d+)M/g; 
            my $offset = 0;
            foreach(@positions) {$offset += $_;} 
         
            $matchinfo = $cigar_start.$del_length.'M'. $cigar_end;  # add number of deleted bases to previous match to account for added "-"
            
            splice(@datapts, $offset,0, @missing); 
			splice(@qual, $offset,0, @missing_phred);   
            $matchinfo=join_M($matchinfo);            
        }
    }
    
    my $out = join("",@datapts);
	my $phred_out = join("",@qual);
    return($matchinfo, $out,$phred_out);
}

sub join_M {
    my $cigar = shift @_;

    while($cigar =~ /(^.*?)(\d+)M(\d+)M(.*$)/){
        my $start=$1;
        my $offset = $2+$3;
        my $end= $4;
        $cigar = $start . $offset .'M' .$end;
    }
    return $cigar;
}


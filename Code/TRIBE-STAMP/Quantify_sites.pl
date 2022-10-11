#!/usr/bin/env perl
 
 ### Author : Mathieu Flamand - Duke University
 ### version : 1.1
 ### date of last modification : 2021-11-03
 
### This programs identifies editing sites by comparing a DART/TRYBE matrix to a control matrix file or to the genomic sequence. 
### For strand information, a refFlat file is provided, sites found within annotated features will be  idenitified according to provided settings.

use v5.26;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Carp;
use MCE::Loop;
use Math::CDF qw(:all);
use Data::Dumper; ### for debugging - to be removed
use Array::IntSpan; 

$|=1;

my $USAGE = "$0 --ClusterFile clusters.bed --EditedMatrix DART.matrix.gz --outfile output_file.out\n";

###set option of program using Getopt::long
my ($ClusterFile, $tablename, $OUTFILE, $verbose, $help,$barcode_option,$extUTR)='';
my $scratch = "./SCRATCH";

##########################################################
### Read environmental variable for slurm task manager ###
##########################################################
my $memory = $ENV{'SLURM_MEM_PER_NODE'} // 4096; #4G default for Duke nodes 
my $ncpu= $ENV{'SLURM_CPUS_PER_TASK'} // 1;

###settings for selection of sites
my $experiment = 'C2T'; # by default it will look for C-to-U transition
my $low_thresh = 2.5;
my $high_thresh =  95;
my $Edited_mincovthresh = 10;

###set option of program using Getopt::long
GetOptions ("ClusterFile:s"=>\$ClusterFile,
			"d|EditedMatrix=s"=>\$tablename,
			"e|editType=s"=>\$experiment,
			"minEdit=f"=>\$low_thresh,
			"maxEdit=f"=>\$high_thresh,
			"o|outfile=s"=>\$OUTFILE,
			"v|verbose"=>\$verbose,
			"h|help"=>\$help,
			"EditedMinCoverage|cov=i"=>\$Edited_mincovthresh,
			"SingleCell|bc|barcode"=>\$barcode_option,
			"cpu=i"=>\$ncpu,
		) or error_out();



if ( $help ) { error_out(); } 

unless($ClusterFile && $tablename)
{
	say STDERR "Error. Please provide a Bed file, an edited matrix and a genome fasta file";
	error_out();
}

### check for all files before proceeding ###
if($ClusterFile and ! -f $ClusterFile){say "Error. Bed file $ClusterFile could not be found";exit();}
unless(-f $tablename){say "Error, could not find edited matrix $tablename";exit();}
###		

### check which type of chr in annotation file
my $cluster_stem ="";
if ($ClusterFile){$cluster_stem = check_chr($ClusterFile, 0);}

my $table_stem = check_chr($tablename, 0);
if ($table_stem =~ /^chr/i and ! $cluster_stem =~ /^chr/i){say "Error, annotation file and matrix files do not have the same chromosome annotation. Cannot match UCSC or Ensembl style."; exit();}

### Parse options ###

unless ($experiment =~ /[ACUTG]2[ACTGUI]/i)
{
	say "Error. the type of experiment provided did not meet requierements. Please use the following format :C2U, A2I ...\n";
	error_out();
}

#process edit threshold data
$low_thresh /= 100;
$high_thresh /= 100;

# define which bases will be looked at for conversion from the options, than adjust to DNA bases.
my $noneditbase = (split("2",$experiment))[0]; 
my $editbase = (split("2",$experiment))[1];
$noneditbase = "T" if $noneditbase eq "U";
$editbase = "T" if $editbase eq "U";
$editbase = "G" if $editbase eq "I";

### identify the reverse complement for negative strand
my %complement = ("A"=>"T","T"=>"A","C"=>"G","G"=>"C"); 
my $noneditbaseREV = $complement{$noneditbase};
my $editbaseREV = $complement{$editbase};

if ($verbose)
{
	say STDERR "Processing files $ClusterFile and $tablename";
	say STDERR "Counting $noneditbase - to - $editbase transitions";
	# say STDERR "Detecting bases with $low_thresh to $high_thresh edit";
	# say STDERR "Miminum enrichement of editing of control is $WToverKO";
	say "processing using $ncpu cpu cores" if $verbose; 
}

## manage different options for output files and formats
my $out_fh;
my $bed_fh;

if ($OUTFILE)
{
	
	$OUTFILE =~ /(.+)\.[^.]+$/;
	my $out_stem = $1; # removes extension from outfile
	my $out_dir = dirname($OUTFILE);
	mkdir($out_dir) unless (-d $out_dir);
	### print header of output file
	my $out_file= "$out_stem".".bed";
	open($out_fh,">", $out_file) or croak "Can't open $out_file for writing: $!\n";
	# say $out_fh join("\t", qq(\#chr start end gene score strand Possible_sites Edited_sites Non_edited Edited Cov));
}
else{ ## if not output file given
	$out_fh = \*STDOUT;
	# say join("\t", qq(\#chr start end gene score strand Possible_sites Edited_sites Non_edited Edited Cov));

}

my $hr_genes = {};
open(my $annotation_handle,"sort -k1,1 -k2,2n $ClusterFile |") or die "Can't open $ClusterFile: $!\n";  

foreach my $line (<$annotation_handle>)
{  
    chomp $line;
    next if ($line=~/^[\#\n]/);
    my ($chr,$start,$end,$ID,$strand) = (split(/\t/,$line))[0,1,2,3,5];
    #adjust by 1
    $hr_genes->{$chr}->{$end} = {END=>$end+1,ID=>$ID,STRAND=>$strand}; 

}
close $annotation_handle;

my $index_list = {};

if (-f $tablename and $tablename =~ /\.gz$/ and -e "$tablename.tbi")
{
	open(my $index, "tabix --list-chroms $tablename |") or die "Failed to read index of $tablename";
	while(<$index>){chomp;$index_list->{$_}=1;}
}
else
{
	say "The Input matrix needs to be bgzip compressed and indexed with Tabix";exit(0);
}

my @chr_list = sort keys %{$hr_genes};

MCE::Loop->init(
	chunk_size => 1,
	max_workers => $ncpu, ## use number of specified threads for processing 
	loop_timeout => 1800, # defaults to 60minutes max per chromosome
	on_post_exit => sub 
	{
		my ($mce, $e) = @_;
		say "done with chromosome: $e->{msg}" if $verbose;
		unless ($e->{status} eq 42){
			say "worker died, Something went wrong, perhaps check if ran out of memory?";
			MCE->abort();}
		$mce->restart_worker; # restart worker after each run to release memory 
	},
);

mce_loop
{
	my ($mce, $chunk_ref, $chunk_id)=@_;
	my $chr = $_;
	
	#open filehandles for contol and dart files	and return if they don't exist. avoids extra processing before reading %gene
		my $dart_fh;
		return unless defined $index_list->{$chr};
        return unless defined $hr_genes->{$chr};
		open ( $dart_fh, "tabix $tablename $chr|") or say "can't read file $tablename with tabix" and return;

	    my $pos_array_for = Array::IntSpan->new(); # initialize new ranges
        my $pos_array_rev = Array::IntSpan->new();	

        # prepare look up ranges for each type of mRNA features
        my ($left,$right);
        
        my $last_dart_line='';
		my ($dart_line, $dart_flag) = read_matrix_line($dart_fh, $last_dart_line); # read next line are return last_line and trigger eof flag.
        my ($potential_sites,$edited_sites, $sum_edit,$sum_nonedit,$sum_score) = map {0}(0 .. 4);
        
        while($dart_line)
        {
            last if $dart_flag;
            my ($dart, $ecoord) = load_line($dart_line);
            
            my @outline;
            
            if (defined $hr_genes->{$chr}->{$ecoord})
			{
				my $strand = $hr_genes->{$chr}->{$ecoord}->{STRAND};
				my $ID = $hr_genes->{$chr}->{$ecoord}->{ID};
			
				@outline = find_sites($chr, $ID, $ecoord, $strand, $dart );
    
                if (@outline)
                {
                    my $result = join("\t", @outline);
                    MCE->say ($out_fh, $result);
                }
                @outline = ();
			}

            ($dart_line, $dart_flag) = read_matrix_line($dart_fh, $last_dart_line); # read next line
            $last_dart_line = $dart_line;
        }
	
		MCE->exit(42, "$chr");
}\@chr_list;

close($out_fh);



if ($OUTFILE)
{

$OUTFILE =~ /(.+)\.[^.]+$/;
my $path_name = $1;
my $out_file= "$path_name".".bed";
my $temp = $out_file . ".tmp";

open(my $in_fh,"<", $out_file);
open(my $fh_temp, ">", $temp);


	foreach my $line (<$in_fh>)
	{  
		chomp $line;
		my ($chr,$start,$end,$ID,$strand) = (split(/\t/,$line))[0,1,2,3,5];
	
		if (defined $hr_genes->{$chr}->{$end}){
			say $fh_temp $line; 
			delete $hr_genes->{$chr}->{$end};
		}
	
	}
	close($in_fh);
	foreach my $chr (keys %{$hr_genes})
	{
		foreach my $missing (keys %{$hr_genes->{$chr}})
		{
			my $IDs = $hr_genes->{$chr}->{$missing}->{ID};
			my $strand = $hr_genes->{$chr}->{$missing}->{STRAND};
			my @ret = ($chr, $missing-1, $missing, $IDs, 0, $strand, 1, 0, 0, 0, 0 );
			say $fh_temp join("\t",@ret); 
		}
	}
	close($fh_temp);

	say "sorting $out_file" if $verbose;
	system("LC_ALL=C sort -k1,1 -k2,2n -u $temp > $out_file");
	system("rm $temp");
	# system("mv $temp $out_file");

}
	

### subroutines ###

#Sub for exiting upon wrong option
sub error_out
{
	say STDERR "\n\n$USAGE";

print STDERR <<'USAGE::BLOCK';

Find_edit_site identifies editing sites in RNAseq dataset by comaparing a matrix of coverage to another or to the sequence of a provided genome. The detection will be done only at positions identified in a reference sequence file (refFlat).
Sites can be filtered through a minimal coverage, minimum and maximum editing percentage and on a threshold of editing in the control matrix. 
Filterering based on a score based on the probability of edinting based on the beta distribution of nucleotides in the edited matrix over the control matrix is also available with the --score option.

complete list of options: 

	--ClusterFile	: Cluster bed file containing cluster region to analyze (s)
	--EditedMatrix		: DART RNA count matrix (s)
	--genome		: instead of using a control file, DART dataset will be compared to sequence in fasta file. Provide path to folder containing multiple fasta or to a single fasta file (s)
	--editType		: Conversion to quantify in the following format : C2U (default) or A2I (s)
	--minEdit		: minimal editing percentage for sites to be considered default to 10%. please write as a number without % symbol (f)
	--maxEdit		: maximum editing percentage for sites to be considered default to 80%. please write as a number without % symbol (f)
	--editFoldThreshold	: Threshold edited/control mutation rate, defaults to 1: by defaults no filtering is done. number can be increased for filtering here, or data can be filtered later. (f)
	--MinEditSites		: minimum number of edit sites, defaults to 2 (i)
	--ControlMinCoverage	: minimum coverage for in control file for sites to be considered(i)
	--EditedMinCoverage	: minimum coverage for in control file for sites to be considered: defaults to 10 (i)
	--outfile		: Output bed file, defaults to STDOUT (s)
	--verbose		: display extra log information
	
USAGE::BLOCK

	exit (1);

}

sub find_sites
{
	#my ($chr, $ID, $bp, $strand, $dart) = @_;
	
	my $chr = shift @_;
	my $ID = shift @_;
	my $bp = shift @_;
	my $strand = shift @_;
	my $dart = shift @_;

    ## define values for dart dataset
	my $dart_non_edit = ($strand eq '+') ? $dart->{$bp}->{$noneditbase} : $dart->{$bp}->{$noneditbaseREV}; 
	my $dart_edit = ($strand eq '+') ? $dart->{$bp}->{$editbase} : $dart->{$bp}->{$editbaseREV}; 
	
	return if ($dart_non_edit + $dart_edit == 0); # prevents possible divisions by 0

	#edit for output
	my $dart_total = $dart->{$bp}->{total} - $dart->{$bp}->{N};
	return if  $dart_total < $Edited_mincovthresh;
	my $dart_edit_ratio = $dart_edit / $dart_total;
	my $dart_edit_ratio_out = sprintf( "%.6f", $dart_edit_ratio);
	my $edited_sites=0; 
    if ($dart_edit_ratio > $low_thresh){$edited_sites++;}
    
    my @outline = ($chr, $bp-1, $bp, $ID, $dart_edit_ratio_out, $strand, 1 ,$edited_sites, $dart_edit,$dart_non_edit,$dart_total) ;
	return @outline;
}


#subroutine to read next matrix line
sub read_matrix_line 
{
	my $fh = shift @_;
	my $last_line = shift @_;
	
	if ($fh and my $line = <$fh>)
	{
		chomp $line;
		return ($line,"0");
	}
	else
	{
		return ($last_line,"1"); # This will return the last line of again if eof is reached
	}
}

sub load_line
{
	my $line = shift @_;
	my %hash;
	my($coord, $A, $T, $C, $G,$N, $total) = (split(/\t/, $line))[1,2,3,4,5,6,7];

	$hash{$coord}{A} = $A;
	$hash{$coord}{T} = $T;
	$hash{$coord}{C} = $C;
	$hash{$coord}{G} = $G;
	$hash{$coord}{N} = $N;
	$hash{$coord}{total} = $total;
	return (\%hash,$coord);
}


sub check_chr 
{
	my ($file, $index) = @_;
	my $fh;
	if ($file =~ /\.gz$/)
	{
		open($fh, "zcat $file | head -5 |") or die "Cannot open ($file) for reading: $!";
	} 
	else
	{
		open($fh, '<', $file) or die "Cannot open ($file) for reading: $!";
	}

	my $stem ='';
	while (defined (my $line = <$fh>)) 
	{
		my @fields = split("\t", $line);
		next if $line=~/^[\#\n]/;
		$stem = 'chr' if $fields[$index] =~ /^chr/;
		last if $stem;
	}
	close  $fh;
	return $stem;
}


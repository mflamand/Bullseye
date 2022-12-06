#!/usr/bin/env perl
 
 ### Author : Mathieu Flamand - Duke University
 ### version : 1.5.1
 ### date of last modification : 2022-12-6
 
### This programs identifies editing sites by comparing a DART/TRYBE matrix to a control matrix file or to the genomic sequence. 
### For strand information, a refFlat file is provided, sites found within annotated features will be  idenitified according to provided settings.

use v5.26;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Carp;
use MCE::Loop;
use Bio::DB::Fasta;
use Math::CDF qw(:all);
use Array::IntSpan; 

$|=1;

my $USAGE = "$0 --annotationFile refFlat --EditedMatrix DART_RNA_matrix --controlMatrix control_RNA_matrix --outfile output_file\n";

###set option of program using Getopt::long
my ($annotationfile,$known_sites, $tablename, $ControlTableName, $OUTFILE, $verbose, $help, $bed_option,$barcode_option, $intron_flag, $genome,$extUTR,$dbFasta,$fallback,$score,$bed_flag,$skip_chr_check)='';
my @bedfiles;


### Read environmental variable for slurm task manager ###
# if you are not running with SLURM, CPU and available memory can be set with the --cpu #core and --mem 4096 (in Megabytes)
my $memory = $ENV{'SLURM_MEM_PER_NODE'} // 4096; #4G default for Duke nodes 
my $ncpu= $ENV{'SLURM_CPUS_PER_TASK'} // 1;

###settings for selection of sites
my $experiment = 'C2T'; # by default it will look for C-to-U transition
my $low_thresh = 10;
my $high_thresh =  80;
my $WToverKO = 1.5; #1.5 fold over control default
my $min_edit_site = 2;
my $Edited_mincovthresh = 10;
my $Control_mincovthresh = 10;
my $max_bkg; 
my $bkg_ratio; 
my $stranded;

###set option of program using Getopt::long
GetOptions ("a|annotationFile:s"=>\$annotationfile,
			"d|EditedMatrix=s"=>\$tablename,
			"c|controlMatrix=s"=>\$ControlTableName,
			"e|editType=s"=>\$experiment,
			"minEdit=f"=>\$low_thresh,
			"maxEdit=f"=>\$high_thresh,
			"t|editFoldThreshold=f"=>\$WToverKO,
			"o|outfile=s"=>\$OUTFILE,
			"v|verbose"=>\$verbose,
			"h|help"=>\$help,
			"EditedMinCoverage=i"=>\$Edited_mincovthresh,
			"ControlMinCoverage=i"=>\$Control_mincovthresh,
			"bed6"=>\$bed_option,
			"SingleCell|bc|barcode"=>\$barcode_option,
			"cpu=i"=>\$ncpu,
			"mem=i"=>\$memory,
			"me|MinEditSites=i"=>\$min_edit_site,
			"intron"=>\$intron_flag,
			"genome=s"=>\$genome,
			"extUTR:5000"=>\$extUTR,
			"fallback=s"=>\$fallback,
			"score:f"=>\$score,
			"filterBed=s"=>\@bedfiles,
			"printFilteredSites"=>\$bed_flag,
			"KnownSites:s"=>\$known_sites,
			"MaxBckg:s"=>\$max_bkg,
			"BckgRatio:s"=>\$bkg_ratio,
			"stranded"=>\$stranded,
			"skip-chr-check"=>\$skip_chr_check,
		) or error_out();


# Error checking before running
if ( $help ) { error_out(); } 
if ($annotationfile and $known_sites){
	say STDERR "Error, please choose either a refflat annotation file (-a) or bed file with known sites (--bedfile), not both";
	error_out();
}
if($ControlTableName and $genome){
	say STDERR "Error, please choose either a control matrix file (-c) or a fasta file (-g), not both";
	error_out();
}
unless(($annotationfile or $known_sites) && $tablename && ($ControlTableName or $genome)){
	say STDERR "Error. Please provide at least an annotation file, an edited matrix and a control matrix or genome file";
	error_out();
}
if($fallback and $genome){
	say STDERR "\n fallback and genome options cannot be used at the same time. plase select one or the other\n";
	error_out();
}
if($annotationfile and ! -f $annotationfile){say "Error. Annotation file $annotationfile could not be found";exit();}
if($known_sites and ! -f $known_sites){say "Error. Annotation file $known_sites could not be found";exit();}
unless(-f $tablename){say "Error, could not find edited matrix $tablename";exit();}
if(! $genome and ! -f $ControlTableName){say "Error, could not find control matrix $ControlTableName";exit();}
if($genome and !(-f $genome or -d $genome)){say "Error, could not find genome file $genome";exit();}
if($fallback and !(-f $fallback or -d $fallback)){say "Error, could not find genome file $fallback";exit();}

unless ($experiment =~ /[ACUTG]2[ACTGUI]/i){
	say "Error. the type of experiment provided did not meet requierements. Please use the following format :C2U, A2I ...\n";
	error_out();
}

### check which type of chr in annotation file
my $annotation_stem ="";
if ($annotationfile){$annotation_stem = check_chr($annotationfile, 2);}
elsif($known_sites){$annotation_stem = check_chr($known_sites, 0);}

my $table_stem = check_chr($tablename, 0);
if (! $skip_chr_check and $table_stem =~ /^chr/i and ! $annotation_stem =~ /^chr/i){say "Error, annotation file and matrix files do not have the same chromosome annotation. Cannot match UCSC or Ensembl style."; exit();}

### Parse options ###
#initialize db if needed
my $genome_stem='';
if ($genome) {
	$dbFasta = Bio::DB::Fasta->new($genome);
	my @IDs = $dbFasta->get_all_primary_ids;
	if( $IDs[0] =~ /^chr/){$genome_stem = 'chr';}
}
if ($fallback) {
	$dbFasta = Bio::DB::Fasta->new($fallback);
	my @IDs = $dbFasta->get_all_primary_ids;
	if( $IDs[0] =~/^chr/){$genome_stem = 'chr';}
}
if($genome or $fallback){
	if (! $skip_chr_check and $table_stem =~ /^chr/i and $genome_stem ne 'chr'){say "Error, genome file and matrix files do not have the same chromosome annotation. Cannot match UCSC or Ensembl style."; exit();}
}

say "processing using $ncpu cpu cores" if $verbose; 

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

if ($verbose){
	say STDERR "Processing files $annotationfile, $ControlTableName and $tablename" unless ($genome or $fallback or $known_sites);
	say STDERR "Processing files $annotationfile, $genome and $tablename" if $genome and not $fallback;
	say STDERR "Processing files $annotationfile, $ControlTableName, $tablename and $fallback" if $fallback and not $genome;
	say STDERR "Detecting $noneditbase - to - $editbase transitions";
	say STDERR "Minimum coverages of $Control_mincovthresh for control matrix and $Edited_mincovthresh for edited matrix";
	say STDERR "Detecting bases with $low_thresh to $high_thresh edit";
	say STDERR "Miminum enrichement of editing of control is $WToverKO";
	say STDERR "3'UTR will be extented by $extUTR nt" if $extUTR;
}

# select default option for scoring
if( length $score && $score == 0){
	$score = 1; #minimum ratio for comparing Dart to Control
}
say "Scoring will be done with following ratio : $score" if ($verbose and $score);


## manage different options for output files and formats
my $out_fh;
if ($OUTFILE){
	$OUTFILE =~ /(.+)\.[^.]+$/;
	my $path_name = $1;
	
	my $out_file= "$path_name".".bed";

	open($out_fh,">", $out_file) or croak "Can't open $out_file for writing: $!\n";
		
	### print header of output file
	say $out_fh join("\t", qq(\#chr start end gene score strand control_ratio control_total dart_ratio dart_total conversion)) unless ($score);
	say $out_fh join("\t", qq(\#chr start end gene score strand control_ratio control_total dart_ratio dart_total conversion score)) if ($score);

}
else{ ## if not output file given
	$out_fh = \*STDOUT;
}

### Start processing ###
### Here we are reading an optional bed files which contains known SNPs or regions to be excluded
### several file can be provided, all regions will be added to the same hash table, strand is not considered for filtering these regions
my $excluded_sites={};
if (@bedfiles){
	say "sites in files: @bedfiles will not be considered in analysis";# if $verbose;
	foreach my $files (@bedfiles){
	my $bed_stem = check_chr($files,0);
	if (! $skip_chr_check and $bed_stem =~ /^chr/i and ! $annotation_stem =~ /^chr/i){say "Error, annotation file and bed files do not have the same chromosome annotation. Cannot match UCSC or Ensembl style."; exit();}
	
	open(my $fh, "<",  $files);
	while(<$fh>){
			next if $_=~/^[\#\n]/;
			my($chr,$start,$end,$strand) = (split("\t",$_))[0,1,2,5];
			foreach my $pos ($start+1..$end){
				$excluded_sites->{$chr}->{$pos}=1;
			}
		}
	}
}

## open annotation filehandle and store in %genes
## there a 2 options: 
## 1) with knownsites option: find sites in the region found in the provided bed file. 
## 2) look in all exons found in a refFlat file. Use --intron to detect in introns. --extUTR to extend 3'UTR for coding RNAs only.
## Sites will be annotated from this gtf file. Priority is given as: 1(low): extended 3'UTR 2 : Introns 3: ncRNA 4: 5'UTR, 5: 3'UTR 6(high) : CDS	

my $genes = {};
if($known_sites){
	open(my $annotation_handle,"<", $known_sites) or die "Can't open $known_sites: $!\n";  

	foreach my $line (<$annotation_handle>){  
		chomp $line;
		next if ($line=~/^[\#\n]/);
		my ($chr,$start,$end,$ID,$strand) = (split(/\t/,$line))[0,1,2,3,5];
		
		my $gene = (split(/\|/, $ID))[0];
	
		# say join("\t", $chr,$start,$end,$strand, $gene);
	
		my $length = $end - $start;
		my $type="Known_site";

		$genes->{'6'}->{$chr}->{$gene}->{length} = $length;
		$genes->{'6'}->{$chr}->{$gene}->{$start} = {END=>$end,ID=>$gene,STRAND=>$strand,TYPE=>$type}; 

	}
	close $annotation_handle;
}
else{
	open(my $annotation_handle,"<", $annotationfile) or die "Can't open $annotationfile: $!\n";   
	print STDERR "loading annotation file $annotationfile... " if $verbose;

	foreach my $line (<$annotation_handle>)	{  
		chomp $line;
		next if ($line=~/^[\#\n]/);
		next if ($line=~/^name/); 
		my($gene,$trxid,$chr,$strand,$txStart,$txEnd,$cdsStart,$cdsEnd,$sString,$eString) = (split(/\t/,$line))[0,1,2,3,4,5,6,7,9,10];
		my @starts = split(/\,/,$sString); 
		my @ends = split(/\,/,$eString);
		my $length = $txEnd - $txStart;

		my $contain_introns = $#ends > 1 ? 1 : 0;

		my @intronStart;
		my @intronEnds;
		
		if ($intron_flag and $contain_introns)	{
			my @intronStart;
			my @intronEnds;

			for (my $i = 0; $i <= $#ends-1;$i++){
				my $j = $i +1;
				my $intron_start=$ends[$i]+1;
				my $intron_end= $starts[$j]-1;
				next if $intron_start >= $intron_end;
				push @intronStart, $intron_start;
				push @intronEnds, $intron_end;
			}

			for (my $i = 0; $i <= $#intronStart;$i++){    
				my $type = "intron"; #default
				$genes->{'2'}->{$chr}->{$trxid}->{length} = $length;
				$genes->{'2'}->{$chr}->{$trxid}->{$intronStart[$i]} = {END=>$intronEnds[$i],ID=>$gene,STRAND=>$strand, TYPE=>$type}; 
			}
		}

		unless($cdsStart == $cdsEnd) ## case when ncRNA
		{
			
			if ($extUTR){  # if extUTR option was used, add to 3'UTR of coding genes. coding only to avoid small sno/sn/U RNAs etc may be missing lncRNAs.
			
				if ($strand eq "+")	{
					my $type = "3UTRe"; 
					my $new_end = $ends[-1] + $extUTR;
					
					$genes->{'1'}->{$chr}->{$trxid}->{length} = $length;
					$genes->{'1'}->{$chr}->{$trxid}->{$ends[-1]} = {END=>$new_end,ID=>$gene,STRAND=>$strand, TYPE=>$type}; 
				}
				else{
					my $type = "3UTRe"; 
					my $new_start = $starts[0] - $extUTR;

					$genes->{'1'}->{$chr}->{$trxid}->{length} = $length;
					$genes->{'1'}->{$chr}->{$trxid}->{$new_start} = {END=>$ends[-1],ID=>$gene,STRAND=>$strand, TYPE=>$type}; 
				}
			}
			
			my $flag_5prime=0;
			my $flag_3prime=0; 

			for (my $i = 0; $i <= $#starts;$i++){    # define 5'UTR and 3'UTR regions in coding RNAs. For annotation purposes, creates new exon boundary
			
				$flag_5prime = 1 if $ends[$i] == $cdsStart;
				$flag_3prime = 1 if $ends[$i] == $cdsEnd;
				
				if ($ends[$i] > $cdsStart and $flag_5prime ==0){
					push @starts, $cdsStart;
					my $new_end = $cdsStart-1;
					push @ends, $new_end;
					$flag_5prime=1;
				}
				if ($ends[$i] > $cdsEnd and $flag_3prime ==0){
					my $new_start = $cdsEnd+1;
					push @starts, $new_start;
					push @ends,$cdsEnd;
					$flag_3prime=1;
				}
			}
		}
		@starts = sort { $a <=> $b } @starts;
		@ends = sort { $a <=> $b } @ends;

		for (my $i = 0; $i <= $#starts;$i++){   
			my $type;
			if($cdsStart == $cdsEnd){
				$type="ncRNA";
			}
			else{
				$type = "CDS"; #default
				
				if($starts[$i] < $cdsStart){ $type = $strand eq "+" ? "5UTR" : "3UTR" ;}
				if($starts[$i] > $cdsEnd){ $type = $strand eq "+" ? "3UTR" : "5UTR" ;}
			}
			if ($type eq "ncRNA"){
				$genes->{'3'}->{$chr}->{$trxid}->{length} = $length;
				$genes->{'3'}->{$chr}->{$trxid}->{$starts[$i]} = {END=>$ends[$i],ID=>$gene,STRAND=>$strand, TYPE=>$type};
			}
			elsif($type eq "5UTR"){
				$genes->{'4'}->{$chr}->{$trxid}->{length} = $length;
				$genes->{'4'}->{$chr}->{$trxid}->{$starts[$i]} = {END=>$ends[$i],ID=>$gene,STRAND=>$strand, TYPE=>$type};
			}
			elsif($type eq "3UTR"){
				$genes->{'5'}->{$chr}->{$trxid}->{length} = $length;
				$genes->{'5'}->{$chr}->{$trxid}->{$starts[$i]} = {END=>$ends[$i],ID=>$gene,STRAND=>$strand, TYPE=>$type};
			}
			elsif($type eq "CDS"){
				$genes->{'6'}->{$chr}->{$trxid}->{length} = $length;
				$genes->{'6'}->{$chr}->{$trxid}->{$starts[$i]} = {END=>$ends[$i],ID=>$gene,STRAND=>$strand, TYPE=>$type};
			}
		}
	}
	close $annotation_handle;
}

## When using a control dataset, we verify that the matrix file exist and is indexed. Then we extract the name of all chromosomes.   
my $index_list = {};
unless ($genome)
{
	if (-f $ControlTableName and $ControlTableName =~ /\.gz$/ and -e "$ControlTableName.tbi"){
		open(my $index, "tabix --list-chroms $ControlTableName |") or die "Failed to read index of $ControlTableName";
		while(<$index>){chomp;$index_list->{'control'}->{$_}=1;}
	}
	else{
		say "Please make sure the matrix files are bgzip compressed and indexed with Tabix";
		error_out();
	}
}

## we do the same for the actual matrix
if (-f $tablename and $tablename =~ /\.gz$/ and -e "$tablename.tbi"){
	open(my $index, "tabix --list-chroms $tablename |") or die "Failed to read index of $tablename";
	while(<$index>){chomp;$index_list->{'DART'}->{$_}=1;}
}
else{
	say "Please make sure the matrix files are bgzip compressed and indexed with Tabix";
	error_out();
}

# we will look in all chromosomes that have a least 1 ORF in the refFlat file or is in the knownsite bed file
my @chr_list = sort keys %{$genes->{'6'}};

MCE::Loop->init(
	chunk_size => 1,
	max_workers => $ncpu, ## use number of specified threads for processing 
	loop_timeout => 1800, # defaults to 60minutes max per chromosome
	on_post_exit => sub{
		my ($mce, $e) = @_;
		say "done with chromosome: $e->{msg}" if $verbose;
		unless ($e->{status} eq 42){
			say "worker died, Something went wrong, perhaps check if ran out of memory?";
			MCE->abort();}
		$mce->restart_worker; # restart worker after each run to release memory 
	},
);

## This anonymous sub will do the heavy lifting, splitting the work for each chromosome
## 1) we define genomic ranges containing features of interest on each strand of the chromosome using Array::IntSpan
## 2) We read the matrix file in parallel with a control matrix file and/or a handle to an indexed fasta file
## 3) We identitify to editing sites with the find_site() sub. These sites a printed to the output file
## 4) the sites found in the exluded file(s) are returned from the sub.

my @filtered_sites = mce_loop
{
	my ($mce, $chunk_ref, $chunk_id)=@_;
	my $chr = $_;
	
	my @bad_sites;
	my ($control_fh, $dart_fh);
		
	unless ($genome){
		return unless defined $index_list->{'control'}->{$chr};
		open ( $control_fh, "tabix $ControlTableName $chr|") or say "can't read file $ControlTableName with tabix" and return;
	}
	
	return unless defined $index_list->{'DART'}->{$chr};
	open ( $dart_fh, "tabix $tablename $chr|") or say "can't read file $tablename with tabix" and return;
	
	my $pos_array_for = Array::IntSpan->new(); # initialize new ranges
	my $pos_array_rev = Array::IntSpan->new();	

	# prepare look up ranges for each type of mRNA features
	foreach my $type_order (sort {$a <=> $b}  keys %{$genes}){
		foreach my $transcript ( sort { $genes->{$type_order}->{$chr}->{$a}->{length} <=> $genes->{$type_order}->{$chr}->{$b}->{length} } keys %{$genes->{$type_order}->{$chr}}){
			# my ($strand,$type,$ID);
			my ($left,$right);
		
			delete $genes->{$type_order}->{$chr}->{$transcript}->{length}; #remove length from hash once sorted in order to get only starts

			foreach my $start (sort {$a <=> $b} keys %{$genes->{$type_order}->{$chr}->{$transcript}}) { 
				my $end = $genes->{$type_order}->{$chr}->{$transcript}->{$start}->{END};
				my $ID = $genes->{$type_order}->{$chr}->{$transcript}->{$start}->{ID}; ## gene name
				my $strand = $genes->{$type_order}->{$chr}->{$transcript}->{$start}->{STRAND};
				my $type = $genes->{$type_order}->{$chr}->{$transcript}->{$start}->{TYPE};
				$ID = "$ID|$type";

				if ($start <= $end){
					$pos_array_for->set_range($start,$end,$ID) if $strand eq "+";
					$pos_array_rev->set_range($start,$end,$ID) if $strand eq "-";
				}
				else{
					$pos_array_for->set_range($end,$start,$ID) if $strand eq "+";
					$pos_array_rev->set_range($end,$start,$ID) if $strand eq "-";
				}
			} 
		}
	}
	
	my $last_control_line='';
	my $last_dart_line='';
	my $barcode='';

	# read first lines
	my ($control_line, $control_flag) = read_matrix_line($control_fh, $last_control_line) unless ($genome);
	my ($dart_line, $dart_flag) = read_matrix_line($dart_fh, $last_dart_line);

	my $chrom_seq;
	#load chromosome sequence in memory if needed. This will increase memory usage, but should decrease access time?
	if ($genome or $fallback){
		$chrom_seq=$dbFasta->seq($chr);
	}

	# process both edited matrix anx control matrix in parallel (unless genome option is used)
	while($dart_line and ($control_line or $genome)){
		last if $dart_flag;
		$barcode = (split(/\t/, $dart_line))[8] if $barcode_option;

		#load current lines in hash_ref
		my ($dart, $ecoord) = load_line($dart_line);
		my ($control, $ccoord) = load_line($control_line) unless ($genome);
		if ($genome){
			$control ={};
			$ccoord = $ecoord;
		}
		# in case fall back option is used. when value is lower in absent in control matrix, fetch genomic coordinate
		if($ecoord < $ccoord and $fallback and not $dart_flag or $control_flag)	{
			$control ={};
			$ccoord = $ecoord;
		}
		#case when edit coordinate is lower than control - meaning no coverage in control for this position
		if($ecoord < $ccoord and not $dart_flag or $control_flag){
			($dart_line, $dart_flag) = read_matrix_line($dart_fh, $last_dart_line);
		}
		#case when control coordinates are lower than edit - no coverage at given position - read next line
		elsif(($ecoord > $ccoord and not $control_flag or $dart_flag) and not $genome){
			($control_line, $control_flag) = read_matrix_line($control_fh, $last_control_line); 
		}else{
			my @outline;
			my $dart_strand;
			my $control_strand;
			if(defined $dart->{$ecoord}->{strand}){
				if (defined $pos_array_for->lookup($ecoord) or defined $pos_array_rev->lookup($ecoord) ){
					$dart_strand=$dart->{$ecoord}->{strand};
					if (not defined $control->{$ecoord}->{strand}){
						$control_strand = $dart_strand;# if no coverage, or genomic assign same strand
					}else{
						$control_strand = $control->{$ecoord}->{strand}; #otherwise get strand
					}
					if ($dart_strand eq $control_strand){ #when both strands are the same
						
						if ($dart_strand eq "-" and defined $pos_array_rev->lookup($ecoord)){
							my $ID = $pos_array_rev->lookup($ecoord);
							@outline = find_sites($chr, $ID, $ecoord, $dart_strand, $control, $dart, \$chrom_seq)
						}elsif ($dart_strand eq "+" and defined $pos_array_for->lookup($ecoord)){
							my $ID = $pos_array_for->lookup($ecoord);
							@outline = find_sites($chr, $ID, $ecoord, $dart_strand, $control, $dart, \$chrom_seq)
						}				
					}elsif($dart_strand eq "+"){ # case when dart is ahead of control # - before + #
						($control_line, $control_flag) = read_matrix_line($control_fh, $last_control_line); 
					}elsif($dart_strand eq "-"){ # this is when dart is behind control
						if(not $fallback){
							($dart_line, $dart_flag) = read_matrix_line($dart_fh, $last_dart_line);
						}elsif(defined $pos_array_rev->lookup($ecoord)){
							my $ID = $pos_array_rev->lookup($ecoord);
							@outline = find_sites($chr, $ID, $ecoord, $dart_strand, $control, $dart, \$chrom_seq)
						}
					}
				}
				if (@outline){
					my $strand = $outline[5];
					push(@outline, $barcode) if $barcode_option; # add barcode when processing them
					my $result = join("\t", @outline);
					if (defined $excluded_sites->{$chr}->{$ecoord})	{push(@bad_sites, $result);}
					else{MCE->say ($out_fh, $result);}
					@outline = ();
				}
				$dart ={};
				$control = {};	
				($dart_line, $dart_flag) = read_matrix_line($dart_fh, $last_dart_line); # read next line
			}else{
				if (defined $pos_array_for->lookup($ecoord)){
					my $strand = "+";
					my $ID = $pos_array_for->lookup($ecoord);
					@outline = find_sites($chr, $ID, $ecoord, $strand, $control, $dart, \$chrom_seq);
				}
			
				if (@outline){
					my $strand = $outline[5];
					push(@outline, $barcode) if $barcode_option; # add barcode when processing them
					my $result = join("\t", @outline);
					if (defined $excluded_sites->{$chr}->{$ecoord})	{push(@bad_sites, $result);}
					else{MCE->say ($out_fh, $result);}
					@outline = ();
				}
				if (defined $pos_array_rev->lookup($ecoord)){
					my $strand = "-";
					my $ID = $pos_array_rev->lookup($ecoord);
					@outline = find_sites($chr, $ID, $ecoord, $strand, $control, $dart, \$chrom_seq);
				}
				if (@outline){
					my $strand = $outline[5];
					push(@outline, $barcode) if $barcode_option;
			
					my $result = join("\t", @outline);

					if (defined $excluded_sites->{$chr}->{$ecoord}){push(@bad_sites, $result);}
					else{MCE->say ($out_fh, $result);}
					@outline = ();
				}
			$dart ={};
			$control = {};			
			
			($dart_line, $dart_flag) = read_matrix_line($dart_fh, $last_dart_line); # read next line

			}

		$last_control_line = $control_line;
		$last_dart_line = $dart_line;

		}
	}
	close ($dart_fh);
	delete $excluded_sites->{$chr} if defined $excluded_sites->{$chr};
	delete $genes->{$chr}; 
	MCE->gather(@bad_sites);
	MCE->exit(42, "$chr");
}\@chr_list;

# sort output files and print excluded_sites if requested
if ($OUTFILE){
	$OUTFILE =~ /(.+)\.[^.]+$/;
	my $path_name = $1;
	my $out_file= "$path_name".".bed";
	my $temp = $out_file . ".tmp";
	
	say "sorting $out_file" if $verbose;
	system("head -1 $out_file >  $temp");
	system("tail -n +2 $out_file | sort -k1,1 -k2,2n -k5,5n -k12,12 >> $temp");
	system("mv $temp $out_file");
	
	if ($bed_option){
		print_bed($out_file);
	}

	if (@bedfiles)	{
		my $excluded_sites = @filtered_sites;
		say STDERR "$excluded_sites were excluded from the analysis";
		if ($bed_flag)	{
			my $out_file= "$path_name".".excluded_sites.bed";
			say "Excluded sites are found in $out_file";
			open (my $fh, ">", $out_file) or die "Can't open $out_file for writing:$!";
			foreach my $lines(@filtered_sites){
				say $fh $lines;
			}
			close $fh;
		} 
	}
}
else{
	if (@bedfiles){
		my $excluded_sites = @filtered_sites;
		say STDERR "$excluded_sites were excluded from the analysis" if $verbose;
	}
}

### subroutines ###

#Sub for exiting upon wrong option
sub error_out {
	say STDERR "\n\n$USAGE";

print STDERR <<'USAGE::BLOCK';

Find_edit_site identifies editing sites in RNAseq dataset by comaparing a matrix of coverage to another or to the sequence of a provided genome. The detection will be done only at positions identified in a reference sequence file (refFlat).
Sites can be filtered through a minimal coverage, minimum and maximum editing percentage and on a threshold of editing in the control matrix. 
Filterering based on a score based on the probability of edinting based on the beta distribution of nucleotides in the edited matrix over the control matrix is also available with the --score option.

complete list of options: 

	--annotationFile	: mm10 FlatRef Annotation file (s)
	--controlMatrix		: control count matrix (wtRNA or mutYTH) (s)
	--EditedMatrix		: DART RNA count matrix (s)
	--genome		: instead of using a control file, DART dataset will be compared to sequence in fasta file. Provide path to folder containing multiple fasta or to a single fasta file (s)
	--fallback		: Optional Path to genome file (fasta). If provided, will be used to compare mutation to in the case where control matrix has no/low coverage (s)
	--editType		: Conversion to quantify in the following format : C2U (default) or A2I (s)
	--minEdit		: minimal editing percentage for sites to be considered default to 10%. please write as a number without % symbol (f)
	--maxEdit		: maximum editing percentage for sites to be considered default to 80%. please write as a number without % symbol (f)
	--editFoldThreshold	: Threshold edited/control mutation rate, defaults to 1: by defaults no filtering is done. number can be increased for filtering here, or data can be filtered later. (f)
	--MinEditSites		: minimum number of edit sites, defaults to 2 (i)
	--ControlMinCoverage	: minimum coverage for in control file for sites to be considered: if using genome and score option, this number will be used for computation. defaults to 10 (i)
	--EditedMinCoverage	: minimum coverage for in control file for sites to be considered: defaults to 10 (i)
	--outfile		: Output bed file, defaults to STDOUT (s)
	--bed6			: output an extra bed6 file, if -o is used, usual file will be generates as .bed, and and extra BED6 will be used with .bed6
	--verbose		: display extra log information
	--inton			: site detection will be done in both exons and introns
	--extUTR		: Extend 3'UTR region by 5kb for site detection
	--score			: score and filter reads based on probability of edit based on beta distribution curve. Sites with 10 fold higher probability in DART/control will be kept (f)
	--filterBed		: bed file of sites to exclude from analysis (SNP or APOBEC background). Multiple files can be provided by using --filterBed multiple times
	--printFilteredSites	: Additionally output fitered sites in a outfile.excluded_sites.bed
	--KnownSites	: Instead of a refFlat annotation file, a bed file containing regions to probed can be defined. For now this needs to be a bed6 file.
	
USAGE::BLOCK

	exit (1);

}

## subroutine to read next matrix line
## we take in the filehandle and the valid last line. We return the next line and the 0 flag. If end of file is reached, we return last line and the 1 FLAG. 
sub read_matrix_line {
	my $fh = shift @_;
	my $last_line = shift @_;
		
	if ($fh and my $line = <$fh>){
		chomp $line;
		return ($line,"0");
	}
	else{
		return ($last_line,"1"); 
	}
}

## subroutine to read matrix line
## We read 1 line, and return a hash reference containing the number of nucleotide at a given position
sub load_line {
	my $line = shift @_;
	my %hash;
	my($coord, $A, $T, $C, $G, $N, $total,$strand) = (split(/\t/, $line))[1,2,3,4,5,6,7,8];
	$hash{$coord}{A} = $A;
	$hash{$coord}{T} = $T;
	$hash{$coord}{C} = $C;
	$hash{$coord}{G} = $G;
	$hash{$coord}{N} = $N;
	$hash{$coord}{total} = $total;
	$hash{$coord}{strand} = $strand;
	return (\%hash,$coord);
}

## Main function to identify edit sites
## This subroutine reads in the current position, strand, data from matrix files and chromosomal sequence if option was used. 
## it finds sites with passing the coverage threshold, mutation thresholds and fold over controls. Optionally return score based on beta distribution
sub find_sites {
	#my ($chr, $ID, $bp, $strand, $control, $dart, $dbFasta) = @_;
	
	my $chr = shift @_;
	my $ID = shift @_;
	my $bp = shift @_;
	my $strand = shift @_;
	my $control = shift @_;
	my $dart = shift @_;
	my $chrom_seq_ref = shift @_;

	#load DART site first, avoid fetching genomic sequencing when using the genome option

	my $dart_non_edit = ($strand eq '+') ? $dart->{$bp}->{$noneditbase} : $dart->{$bp}->{$noneditbaseREV}; 
	my $dart_edit = ($strand eq '+') ? $dart->{$bp}->{$editbase} : $dart->{$bp}->{$editbaseREV}; 
	
	return if ($dart_non_edit + $dart_edit == 0); 
	return unless ($dart_edit >= $min_edit_site);

	my $dart_total = $dart->{$bp}->{total} - $dart->{$bp}->{N}; 
	return if  $dart_total < $Edited_mincovthresh; 
	
	my $dart_edit_ratio = $dart_edit / $dart_total;
	return if ($dart_edit_ratio < $low_thresh or $dart_edit_ratio > $high_thresh); 
	
	my $dart_others = $dart_total - $dart_non_edit - $dart_edit;
	my $other_ratio  = $dart_others/$dart_total;
	if ($max_bkg){
		return if $other_ratio >= $max_bkg; # return if ratio of mutations is higher than max allowed levels of background mutation 
	}
	if ($bkg_ratio){ 
		return if ($dart_edit_ratio <= $other_ratio * $bkg_ratio);
	} 
	
	my $dart_edit_ratio_out = sprintf( "%.4f", $dart_edit_ratio);

	my $dart_edit_conf;
	# calculate confidence based on beta distribution
	if ($score){
	$dart_edit_conf = pbeta($low_thresh, $dart_edit,$dart_non_edit); # get probability that value is higher than minimal threshold using beta distribution
	$dart_edit_conf = 1 unless length $dart_edit_conf; # This is to avoid division by zero when no significance is found
	$dart_edit_conf = 4.94e-324 if $dart_edit_conf == 0; # This is to avoid divisions by zero when value is smaller the double-precision floating-point value 
	}

	#process control 
	my ($control_non_edit,$control_edit, $control_total, $base); #declare variables
	if ($genome or ($fallback and (! defined $control->{$bp}->{'total'} or $control->{$bp}->{'total'} < $Control_mincovthresh))){
		my $nuc = uc(substr ($$chrom_seq_ref, $bp-1, 1));
		$control = getEmptyHash($bp, $nuc);
	}
	
	$control_non_edit = ($strand eq '+') ? $control->{$bp}->{$noneditbase} : $control->{$bp}->{$noneditbaseREV}; 
	$control_edit = ($strand eq '+') ? $control->{$bp}->{$editbase} : $control->{$bp}->{$editbaseREV}; 
	$control_total = $control->{$bp}->{total} - $control->{$bp}->{N};# remove ambiguous nucleotides from total
	return if $control_total < $Control_mincovthresh;

	# check if non edit in control is at least 60%. avoids divisions by 0 in some cases. Should eliminate SNPs hets ~50% and remove highly edited sites
	return if $control_non_edit / $control_total < 0.60; 

	#edit for output
	my $control_edit_ratio = $control_edit / $control_total;
	return if ($control_edit_ratio > 0.0 && $dart_edit_ratio/$control_edit_ratio < $WToverKO ); ## this is hard cutoff for edit DART/mut

	my $control_others = $control_total - $control_non_edit - $control_edit;
	my $other_ratio_ctrl  = $control_others/$control_total;
	if ($max_bkg){
		return if $other_ratio_ctrl >= $max_bkg; # return if ratio of mutations is higher than max allowed levels of background mutation 
	}

	my $control_edit_ratio_out = ($control_edit_ratio != 0) ? sprintf( "%.4f", $control_edit_ratio) : 0  ;

	my $control_edit_conf;
	my $edit_ratio;

	if ($score){
	# calculate confidence based on beta distribution
	$control_edit_conf = pbeta($low_thresh, $control_edit,$control_non_edit);
	$control_edit_conf = 1 unless length $control_edit_conf;
	$control_edit_conf = 4.94e-324 if $control_edit_conf == 0;
	# check control/dart ratio
	$edit_ratio = log10($control_edit_conf/$dart_edit_conf);
	return if $edit_ratio < $score;
	}

	$ID = $ID . "|" . $experiment . "|mut=" . $dart_edit;

	if ($control_edit_ratio_out){
		$ID = $ID . "|" . sprintf( "%.4f", $dart_edit_ratio/$control_edit_ratio);
	} 
	else{
		$ID = $ID . "|NA" ;
	}

	if ($fallback and defined $control->{$bp}->{'tag'}){
		$ID = $ID . "|genome";
	}
	elsif($fallback){
		$ID = $ID . "|control";
	}

	my @outline = ($chr, $bp-1, $bp, $ID, $dart_edit_ratio_out, $strand, $control_edit_ratio_out, $control_total, $dart_edit_ratio_out, $dart_total, $experiment); #prepare output line
	push (@outline, $edit_ratio) if $score;
	return @outline;
}

sub getEmptyHash{

	my $bp = shift @_;
	my $nuc = shift @_;

	my $hash_ref={};

	$hash_ref->{$bp}->{A} += 0;
	$hash_ref->{$bp}->{T} += 0;
	$hash_ref->{$bp}->{C} += 0;
	$hash_ref->{$bp}->{G} += 0;
	$hash_ref->{$bp}->{N} += 0;
	$hash_ref->{$bp}->{$nuc} += $Control_mincovthresh;
	$hash_ref->{$bp}->{total} += $Control_mincovthresh;
	$hash_ref->{$bp}->{'tag'} = 'g';

	return $hash_ref;

}

## helper sub
sub log10 { 
    my $n = shift; 
    return log($n) / log(10); 
} 

## check that all files have the same types of chrmosomes in first column eg.: chr1 vs 1
sub check_chr{
	my ($file, $index) = @_;
	my $fh;
	if ($file =~ /\.gz$/){
		open($fh, "zcat $file | head -5 |") or die "Cannot open ($file) for reading: $!";
	} 
	else{
		open($fh, '<', $file) or die "Cannot open ($file) for reading: $!";
	}

	my $stem ='';
	while (defined (my $line = <$fh>)){
		my @fields = split("\t", $line);
		next if $line=~/^[\#\n]/;
		$stem = 'chr' if $fields[$index] =~ /^chr/;
		last if $stem;
	}
	close  $fh;
	return $stem;
}

## print a bed6 file, used if bed6 is asked as output
sub print_bed{
	my $file = shift @_;

	$file =~ /(.+)\.[^.]+$/;
	my $path_name = $1;
	my $bed_file= "$path_name".".bed6";

	open(my $fh, "<", $file) or croak "Can't open $file for reading";
	open(my $bed_fh, ">" ,"$bed_file") or croak "Can't open $bed_file for writing";

	foreach (<$fh>){
		my @out = (split("\t",$_))[0..5];
		my $barcode = (split("\t",$_))[-1] if $barcode_option;;
		push @out, $barcode  if $barcode_option;
		my $out = join("\t", @out);
		say $bed_fh $out;
	}
}
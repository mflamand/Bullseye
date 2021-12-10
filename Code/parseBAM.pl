#!/usr/bin/env perl

 ### Author : Mathieu Flamand - Duke University
 ### version : 1.22 
 ### date of last modification : 2021-12-6

use strict;
use warnings;
use v5.26;
use Getopt::Long;
use File::Basename;
use Carp;
use MCE::Loop;
use File::Path;
use Array::IntSpan; 

my ($file, $verbose, $outfile, $max_barcode,$read_thresh, $barcode_filter, $remove_duplicate,$remove_MM, $cell_ID_flag,$bc_pattern,$help);
my $mode = "bulk";
my $min_cov = 1;
my $time_out = 7200;
my $memory = $ENV{'SLURM_MEM_PER_NODE'} // 8192; 
my $ncpu= $ENV{'SLURM_CPUS_PER_TASK'} // 1;

my $exclude;

###set option of program using Getopt::long
GetOptions ("i|input:s"=>\$file,
			"v|verbose"=>\$verbose,
			"o|output:s"=>\$outfile,
			"mode:s"=>\$mode, 
			"min|minCoverage:i"=>\$min_cov,
			"nb|MaxBarcodes:i"=>\$max_barcode,
			"rt|readThreshold:i"=>\$read_thresh,
			"Cell_ID_pattern=s"=>\$bc_pattern,  ## to be changed to single cell, keep barcode for compatibility for now
			"cpu:i"=>\$ncpu,
			"filterBarcodes=s"=>\$barcode_filter,
			"removeDuplicates"=>\$remove_duplicate,
            "removeMultiMapped"=>\$remove_MM,
			"h|help"=>\$help,
			"mem=i"=>\$memory,
            "ribosome:s"=>\$exclude,
		) or die "Error in command line arguments, please use --help for information on usage\n";;

my $usage = "$0 -i input.bam (--cpu 1 -o output.matrix -b barcode.txt -min 1 -nb 1000 -rt 50000 --verbose)";

################################################
### Check if all options are correctely used ###
################################################

if ( $help ) { error_out(); } 

unless (-e $file and $file =~/\.bam$/){say "Error, could not find input file"; error_out();}
### Make sure index exists and read file to get list of chromosomes in bam ###
unless (-e "$file.bai"){
	say "creating index file for $file";
	system("samtools index -@ $ncpu $file") and croak "Could not index $file with samtools"; # Exit if command uncesseful (system() return 0 on sucess)
}

if($mode =~ /bulk/i){
	$mode = "bulk";
}
elsif($mode =~ /singlecell/i or $mode eq "sc" or $mode eq "bc" or $mode eq "barcode"){
	$mode="SingleCell";
}
elsif($mode eq "ExtractBarcode"){
	$mode="cell_ID_flag";
}
else{
	say "Error please select one of 'Bulk', 'SingleCell' or 'ExtractBarcode' for the --mode option";
	error_out();
}

my $sorted_tag = qx(samtools view -H $file | grep SO);
unless ($sorted_tag =~ /coordinate/){
	say "could not detect if file is sorted by coordinate, please make file contains a header and is sorted by coordinate";
	error_out();
}

if ($bc_pattern){
	if ($bc_pattern =~ /10x/i){$bc_pattern = 'CB:Z:';}
	elsif ($bc_pattern =~ /SMART/i) {$bc_pattern = 'RG:Z:';}  
	say "using regex $bc_pattern for Cell ID" if $verbose;
}

if ($barcode_filter){
	unless ($bc_pattern){say "\nError barcode option used but no barcode pattern specified, program will exit\n"; error_out();}
}

open(my $index, "samtools idxstats $file |") or croak "cannot run samtools idxstats";
my $chr_list = {};
while(<$index>){
	chomp;
	next if $_ =~ /^\*/; # ignore unmapped reads
	my ($chr_id, $chr_length, $mapped_reads) = (split("\t", $_))[0,1,2];
	$chr_list->{$chr_id}=$chr_length unless $mapped_reads < 1; # only adds if it contains reads
}
close($index);
my @chromosomes= sort keys %{$chr_list};

# define output file
my $fout;
unless($cell_ID_flag){
    if(not $outfile){$fout =  \*STDOUT;}
    elsif($outfile eq '-'){$fout =  \*STDOUT;}
    else{open($fout,">", $outfile) or croak "Can't open file for writing: $!\n";}
}

if ($mode eq "cell_ID_flag"){
	extract_barcodes($file);
    say "Done extracting cell IDs. Exiting now" if $verbose;
	exit(0);
}

my $hr_exclude={};
if ($exclude){
    open(my $bed_in, "<", $exclude) or die "cannot open $exclude for reading";
    while(<$bed_in>){
        chomp;
        next if $_ =~ /^\#/; 
        my ($chr, $start,$end, $strand) = (split("\t", $_))[0,1,2,5];
        unless(defined $hr_exclude->{$chr}){
            $hr_exclude->{$chr}= Array::IntSpan->new();
        }
        $hr_exclude->{$chr} = Array::IntSpan->set_range($start,$end,$strand);
    }
}

# inititalize MCE for faster processing.
MCE::Loop->init(chunk_size => 1,max_workers => $ncpu,loop_timeout => $time_out,
    on_post_exit => sub {
    my ($mce, $e) = @_;
    say "done with chromosome: $e->{msg}" if $verbose;
    unless ($e->{status} eq 1){
        say "worker died, please use more memory or less threads and run again";
        MCE->abort();}
    $mce->restart_worker;},
);

if($mode eq "SingleCell"){
	my %barcode_hash;
	if ($barcode_filter){
		open(my $bc_fh, "<", $barcode_filter) or die "can't open $barcode_filter for reading: $!\n";
		my $line_counter=0;
		while(<$bc_fh>){
			(split, "\t")[0] =~ /:?([A-Z0-9]+)(\-\d)?/i;
			my $barcode = $1;
			last if $read_thresh and (split, "\t")[1] < $read_thresh;
			$barcode_hash{$barcode} = 1;
			$line_counter++;
			if (defined $max_barcode){last if $line_counter == $max_barcode;}
		}
	}

    mce_loop{
		my ($mce, $chunk_ref, $chunk_id)=@_;
		my $chr = $_; 
        my $matrix = {};
        my $counter = 0;
        my $prevchr = "NA";

        open( my $fin, "samtools view $file $chr|") or croak "Can't read bam file $file for chr $chr"; 

        LINE:while (my $line = <$fin>){
            chomp;
            next if $line =~ /^\@/;
            my($flag,$start,$matchinfo,$seq) = (split(/\t/,$line))[1,3,5,9];
            ### check for PCR/optical duplicates and remove if needed ###
            if ($remove_duplicate or $remove_MM){
                $flag = sprintf ("%012b", $flag);
                my ($dup,$secondary_alignement) = (split (//, $flag))[1,3];
                if ($remove_duplicate){
                   next LINE if $dup eq 1;
                }
                if ($remove_MM){
                    next LINE if $secondary_alignement eq 1;
                }
            }


            if ($counter == 10000){
                foreach my $barcode (keys %{$matrix}){
                    foreach my $pos (sort { $a <=> $b } keys %{$matrix->{$barcode}}){
                        if ($pos < $start){
                            my @outline = return_line($pos, $matrix->{$barcode}, $chr);
                            push @outline, $barcode;
                            my $out_line = join("\t", @outline);
                            MCE->say($fout, $out_line ) if ($outline[7] >= $min_cov);
                            delete($matrix->{$barcode}->{$pos});
                        }
                    }
                }
                $counter = 0; #reset counter
            }
            
            #Parse Cigar String to remove soft clips, ins and dels.
            ($matchinfo, $seq) = parse_cigar($matchinfo, $seq,$line);

            my @matches= $matchinfo =~ /(\d+M)/g;
            my @skip = $matchinfo =~ /(\d+N)/g;
            my @datapts = split('',$seq); # split sequence in each nucleotide
            
            my $offset = 0;

            # for each matching region in cigar, splice sequence and adjust starting position with offset
            for (my $i = 0; $i <= $#matches; $i++ ){
                my $match = $matches[$i];
                $match =~ /(\d+)M/;
                my $length = $1;
                if ($i > 0) { $skip[$i-1] =~ /(\d+)N/; $offset += $1;}
                my @sequence = splice(@datapts,0,$length );

                for (my $j = 0; $j <= $#sequence; $j++){
                    my $base = $sequence[$j];
                    next if $base eq '-';	
                    my $bp = $start + $offset + $j;
                    if ($exclude){
                        next if defined  $hr_exclude->{$chr}->lookup($bp) ;
                    }
                    $matrix->{$bp}->{$base}++;
                    foreach my $b ("A","T","C","G","N"){$matrix->{$bp}->{$b}+= 0;} # make sure each entry is non empty
                }
                $offset += $#sequence + 1;
            }
            $counter++;
        }   
        # deal with any leftovers
        foreach my $barcode (keys %{$matrix}){
            foreach my $pos (sort { $a <=> $b } keys %{$matrix->{$barcode}}){
                my @outline = return_line($pos, $matrix->{$barcode}, $chr);
                push @outline, $barcode;
                my $out_line = join("\t", @outline);
                MCE->say($fout, $out_line ) if ($outline[7] >= $min_cov);
                delete($matrix->{$barcode}->{$pos});
            }
        }
        MCE->exit(1, "$chr"); # exit when done with chromosome this will report in STDERR if verbose is used
    }\@chromosomes;
}elsif($mode eq "bulk"){

    mce_loop{
        my ($mce, $chunk_ref, $chunk_id)=@_;
        my $chr = $_; 
        
        my $matrix = {};
        my $counter = 0;
        my $prevchr = "NA";

        open( my $fin, "samtools view $file $chr|") or croak "Can't read bam file $file for chr $chr"; 

        LINE:while (my $line = <$fin>)
        {
            chomp;
            next if $line =~ /^\@/;
            my($flag,$start,$matchinfo,$seq) = (split(/\t/,$line))[1,3,5,9];
        
           if ($remove_duplicate or $remove_MM){
                $flag = sprintf ("%012b", $flag);
                my ($dup,$secondary_alignement) = (split (//, $flag))[1,3];
                if ($remove_duplicate){
                   next LINE if $dup eq 1;
                }
                if ($remove_MM){
                    next LINE if $secondary_alignement eq 1;
                }
            }
        
            if ($counter == 50000){
               foreach my $pos (sort { $a <=> $b } keys %{$matrix}){
                        if ($pos < $start){
                            my @outline = return_line($pos, $matrix, $chr);
                            my $out_line = join("\t", @outline);
                            MCE->say($fout, $out_line ) if ($outline[7] >= $min_cov);
                            delete($matrix->{$pos});
                        }
                    }
                $counter = 0; #reset counter
            }

            ($matchinfo, $seq) = parse_cigar($matchinfo, $seq,$line);
            my @matches= $matchinfo =~ /(\d+M)/g;
            my @skip = $matchinfo =~ /(\d+N)/g;
            my @datapts = split('',$seq); # split sequence in each nucleotide
            my $offset = 0;

            for (my $i = 0; $i <= $#matches; $i++ ){
                my $match = $matches[$i];
                $match =~ /(\d+)M/;
                my $length = $1;
                if ($i > 0) { $skip[$i-1] =~ /(\d+)N/; $offset += $1;}
                my @sequence = splice(@datapts,0,$length );

                for (my $j = 0; $j <= $#sequence; $j++){
                    my $base = $sequence[$j];
                    next if $base eq '-'; # skip if dels
                    my $bp = $start + $offset + $j;
                    if ($exclude){
                        next if defined  $hr_exclude->{$chr}->lookup($bp) ;
                    }
                    $matrix->{$bp}->{$base}++;
                    foreach my $b ("A","T","C","G","N"){$matrix->{$bp}->{$b} += 0;}
                    }
                $offset += $#sequence+1;
            }
            $counter++;
        }
    
        foreach my $pos (sort { $a <=> $b } keys %{$matrix}){
            my @outline = return_line($pos, $matrix, $chr);
            my $out_line = join("\t", @outline);
            MCE->say($fout, $out_line ) if ($outline[7] >= $min_cov);
        }
        MCE->exit(1, "$chr");
    }\@chromosomes;# ($matchinfo, $seq) = parse_cigar($matchinfo, $seq);
} # end of bulk

####################
### File sorting ###
####################

my $outdir = dirname($outfile); # temp files will be the output folder if they are needed
my $sort_memory = int($memory*0.8); # use 80% of available memory for sorting
# fast in place sort with C locale
system("LC_ALL=C sort -k1,1 -k2,2n --parallel=$ncpu -T $outdir -S${sort_memory}M -o $outfile $outfile") == 0 or die "Failed to sort file $file:$!";
# compress with bgzip 
system("bgzip -f -@ $ncpu $outfile") == 0 or die "error compressing final file with bgzip";
# index with tabi x
say "indexing $outfile with tabix" if $verbose;
system("tabix -b 2 -e 2 $outfile.gz") == 0 or die "error indexing final file with tabix";

say STDERR "All done" if $verbose;


# subroutines

sub error_out{
say "\n\n$usage";
print << 'EOF';
	
This program will build a matrix of nucleotide count for every positions mapped in a bam file:

	--mode : one of 'Bulk', 'SingleCell' or 'ExtractBarcodes'. Default to Bulk for processing without considering Barcodes. 
	--input: input bam file used to build matrix. Make sure files are coordinated sorted
	--output: output file name, defaults to STDOUT
	--removeDuplicates: To ignore reads marked as PCR or optical duplicates
    --removeMultiMapped: To ignore multi mapped reads 
	--verbose: display extra information
	--cpu: number of threads to use for processing (only necessary if you are not using slurm)
	--mem: available memory for sorting (M) (only necessary if you are not using slurm)
	--minCoverage: minimum base coverage to output to final file (default = 1)
	--filterBarcodes: only keep the barcodes included in the first column of a provided file. 
				A second column with the number of reads for each barcode can be provided and used for filtering with the following options:
	--MaxBarcodes: (sc) number of barcodes to process (defaults to all). Will stop reading provided barcode file after hitting this number, barcodes should be filered by mapped reads.
	--readThreshold: (sc) minimum number of reads for barcodes to be considered
	--Cell_ID_pattern: pattern for to use for single cell processing. use 10X, SMART or enter a sam tag RegEx for cell identification.
    --exclude: bed file with coordinates for regions to be skipped. (eg. Ribosomal RNAs.) 
	
EOF
    exit(1);
}


sub parse_cigar{
    my $matchinfo = shift @_;
    my $seq = shift @_;
    my $line = shift @_;

    my @datapts = split('',$seq); # split sequence in each nucleotide


    if($matchinfo =~ /^(\d+)S(.+?)(\d+)S$/){   ### if cigar begins and ends with soft clipped, remove softclip 
        $matchinfo = $2;
        splice(@datapts, 0, $1); 
        splice(@datapts, -$3);# remove the number of soft clipped nucleotide
    }

    if($matchinfo =~ /(.+[M|N|I|D])(\d+)S$/) {  ### if cigar ends with soft clipped, remove it 
        $matchinfo = $1;
        splice(@datapts,-$2); 
    }

    if($matchinfo =~ /^(\d+)S(.+[M|N|I|D])$/) { ### if cigar begins with soft clipped, remove the clipped portion  
        $matchinfo = $2;
        splice(@datapts, 0, $1);  # remove the number of soft clipped nucleotide
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
            $matchinfo=join_M($matchinfo);
        }
    }

    while ($matchinfo =~ /D/) {
       if($matchinfo =~ /(^.+?)(\d+)D(.+[M|N])$/) {
       
            my $cigar_start = $1;
            my $cigar_end = $3;
            my $del_length= $2;

            my @missing = map{"-"} 1..$del_length; # add '-' for each missing base
            
            my @positions= $cigar_start =~ /(\d+)M/g; 
            my $offset = 0;
            foreach(@positions) {$offset += $_;} 
         
            $matchinfo = $cigar_start.$del_length.'M'. $cigar_end;  # add number of deleted bases to previous match to account for added "-"
            
            splice(@datapts, $offset,0, @missing);  
            $matchinfo=join_M($matchinfo);        
  
        }
    }

    my $out = join("",@datapts);
    return($matchinfo, $out);
 
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

sub extract_barcodes{
    my $input = shift @_;
    my $stem = basename($input, ".bam");
	my $barcode_file = "barcode.$stem.txt";
    say "Counting barcode and generating list in barcode_file" if $verbose;
    my $linecount = 0; 
    open my $bc_fh, ">", "$barcode_file";
	
    my %barcode_list;

    MCE::Loop->init(chunk_size => 1,max_workers => $ncpu);

    my %barcode_lists_chr = mce_loop{
		my ($mce, $chunk_ref, $chunk_id)=@_;
		my $chr = $_; 
		my %ret;
		
		open( my $fin, "samtools view $input $chr|") or croak "can not open file $input to read:$!\n";; # read indexed bam file by chromosome
	
		while(<$fin>){
            next if $_ =~ /^\@/;
            my $flag = (split(/\t/))[1];
            
            if ($remove_duplicate){
                $flag = sprintf ("%012b", $flag);
                my $dup = (split (//, $flag))[1];
                next if $dup eq 1;
            }
            
            next unless $_ =~ /${bc_pattern}([A-Z0-9]+)(-\d)?\s/i;
            my $barcode = $1;
            $ret{$barcode} += 1;
				
        }
	    MCE->gather($chr, \%ret);
	}\@chromosomes;

    foreach my $chr (@chromosomes){
		my $hash = $barcode_lists_chr{$chr};
		foreach my $barcode ( keys %{$hash}){
			$barcode_list{$barcode} += $hash->{$barcode};
		}
	}
	
	## store file with barcodes
	foreach my $bc (sort { $barcode_list{$b} <=> $barcode_list{$a}  }  keys %barcode_list){
		say $bc_fh "$bc\t$barcode_list{$bc}";
	}
}

sub return_line{
    my $pos = shift @_;
    my $matrix_ref = shift @_;
    my $chr = shift @_;

    my $A = $matrix_ref->{$pos}->{"A"};
    my $T = $matrix_ref->{$pos}->{"T"};
    my $C = $matrix_ref->{$pos}->{"C"};
    my $G = $matrix_ref->{$pos}->{"G"};
    my $N = $matrix_ref->{$pos}->{"N"};
    my $total = $A + $T + $C + $G + $N;
    my @outline = ($chr,$pos,$A,$T,$C,$G,$N,$total);
    
    return @outline; 
}
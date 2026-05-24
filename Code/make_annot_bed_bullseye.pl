#!/usr/bin/perl

# Modified from metaPlotR by Mathieu Flamand, Université Laval, 2026-5-23
# Use Bio::DB::HTS::Faidx instead of Bio::DB::Fasta
# Fetch each exon once, then iterate over the exon sequence string
# Requires a fasta file with .fai index or a bgziped fasta file with index 

use strict;
use warnings;
use Bio::DB::HTS::Faidx;
use Getopt::Long;

my $genomeDir = '';
my $genePred  = '';
GetOptions ('genomeDir=s' => \$genomeDir, 'genePred=s' => \$genePred);

&checkArgs($genomeDir, $genePred) || die "Check command line arguments: $!";

open (my $annot, "<", "$genePred") || die "Cant open genePred file: $!";

my $db = Bio::DB::HTS::Faidx->new($genomeDir);

my @IDs = $db->get_all_sequence_ids();
my %ID_hash = map { $_ => 1 } @IDs;

my $counter = 0;

while (my $line = <$annot>) {
    chomp($line);
    next if ($line =~ /^[\#\n]/);

    my @temp = split(/\t/, $line);
    my ($transcript_ID, $chr, $strand) = @temp[1..3];
    my ($txStart, $txEnd, $cdsStart, $cdsEnd) = @temp[4..7];
    my $exonCount = $temp[8];
    my @exonFrame = split(/,/, $temp[-1]);
    my @exon_starts = split(/,/, $temp[9]);
    my @exon_ends   = split(/,/, $temp[10]);
    my $gene_name   = $temp[12];

    next unless defined $ID_hash{$chr};

    my $gene_size = 0;
    my $num_exons = scalar(@exon_starts);
    for (my $i = 0; $i < $num_exons; $i++) {
        $gene_size += abs($exon_ends[$i] - $exon_starts[$i]);
    }

    $txStart  = $txStart + 1;
    $cdsStart = $cdsStart + 1;
    foreach my $x (@exon_starts) {
        $x = $x + 1;
    }

    my $nucl_j = -1;
    my $count   = 0;

    for (my $i = 0; $i < $exonCount; $i++) {
        my $exon_start = $exon_starts[$i];
        my $exon_end   = $exon_ends[$i];

        my $location = $chr . ":" . $exon_start . "-" . $exon_end;
        my $exon_seq = uc($db->get_sequence_no_length($location));

        next unless defined $exon_seq && $exon_seq ne '';

        my $exon_len = length($exon_seq);

        for (my $offset = 0; $offset < $exon_len; $offset++) {
            my $coord = $exon_start + $offset;

            $count++;
            my $mrna_pos = &determine_mrna_pos($strand, \$gene_size, $count);
            my $feature  = &utrs_or_cds($coord, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd);
            my $nucleotide = substr($exon_seq, $offset, 1);
            my $codon_pos  = &determine_codon_pos($feature, $nucl_j, $strand);

            print $chr, "\t", $coord - 1, "\t", $coord, "\t",
                  join('|', ($transcript_ID, $gene_name, $feature, $mrna_pos)), "\t",
                  $nucleotide, "\t", $strand, "\n";
        }
    }
}

sub determine_codon_pos() {
    my ($feature, $nucl_j, $strand) = @_;
    my $codon_pos = 'NA';

    if ($feature eq 'cds') {
        $nucl_j++;
        $codon_pos = $nucl_j % 3;

        if ($strand eq '-') {
            if ($codon_pos == 0) { $codon_pos = 2; }
            elsif ($codon_pos == 2) { $codon_pos = 0; }
        }
    }
    return $codon_pos;
}

sub determine_mrna_pos() {
    my ($strand, $gene_size_ref, $count) = @_;
    my $mrna_pos = $count;
    if ($strand eq '-') {
        $mrna_pos = $$gene_size_ref;
        $$gene_size_ref -= 1;
    }
    return $mrna_pos;
}

sub checkArgs() {
    my ($genDir, $genPred) = @_;
    my $pass = 1;
    if (! -d $genDir and ! -f $genDir) {
        print "No such directory exists: $genDir \n";
        $pass = 0;
    }
    if (! -f $genPred) {
        print "No such file exists: $genDir \n";
        $pass = 0;
    }
    return ($pass);
}

sub utrs_or_cds () {
    my ($coord, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd) = @_;
    my $feature = 'na';

    if (abs($cdsStart - $cdsEnd) <= 1) {
        $feature = 'ncRNA';
    }
    elsif ($coord >= $txStart && $coord < $cdsStart) {
        if ($strand eq '+') { $feature = '5utr'; }
        else                { $feature = '3utr'; }
    }
    elsif ($coord > $cdsEnd && $coord <= $txEnd) {
        if ($strand eq '+') { $feature = '3utr'; }
        else                { $feature = '5utr'; }
    }
    else {
        $feature = 'cds';
    }

    return $feature;
}
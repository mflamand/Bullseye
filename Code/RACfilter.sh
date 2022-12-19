#!/bin/bash

# use option -r to remove RAC sequence from 7th column useful for metagene analysis

GENOME=/path/to/genome/file.fasta

# option for parsing -r
OPTIND=1 
KEEPRAC=0

while getopts ":r" opt; do
	case "$opt" in
	r)
		REMOVE=1
		;;
	\? ) 
	echo "Usage: RACfilter.sh [-r] file1.bed"
	exit
		;;
	esac
done
shift $((OPTIND-1))

STEM=$(basename "$1" .bed)

if [[ $REMOVE -eq 1 ]]; then
	perl -F'\t' -ane 'if ($_=~/^\#/){next;} if($F[5] =~ /\+/){ $F[1] -= 2; } else { $F[2] += 2; } print join("\t", @F)' $STEM.bed | bedtools getfasta -bedOut -s -fi $GENOME -bed - > $STEM.3nt.sequence.bed
	perl -F'\t' -ane 'if($F[5] =~ /\+/){ $F[1] += 2; } else { $F[2] -= 2; } if( $F[-1] =~ m/[AG]AC/i ){chomp;pop @F; print join ("\t", @F ), "\n"}' $STEM.3nt.sequence.bed > $STEM.RAC.bed
	rm $STEM.3nt.sequence.bed
else
# to keep RAC in 7th column:
	perl -F'\t' -ane 'if ($_=~/^\#/){next;}if($F[5] =~ /\+/){ $F[1] -= 2; } else { $F[2] += 2; } print join("\t", @F)' $STEM.bed | bedtools getfasta -bedOut -s -fi $GENOME -bed - > $STEM.3nt.sequence.bed
	perl -F'\t' -ane 'if($F[5] =~ /\+/){ $F[1] += 2; } else { $F[2] -= 2; } if( $F[-1] =~ m/[AG]AC/i ){ print join "\t", @F }' $STEM.3nt.sequence.bed > $STEM.RAC.bed
	rm $STEM.3nt.sequence.bed
fi

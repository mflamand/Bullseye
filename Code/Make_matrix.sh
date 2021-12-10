#!/bin/bash
#SBATCH --mem=8G
#SBATCH -c 4
#SBATCH -a 1-n
#SBATCH -p scavenger

#adjust number of arrays to correspond to number of bam files to process
#depending on the size of your bam files, you may need to adjust the amount of RAM used

file=$(ls *.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)   ### can change ls to subset the analysis to some bam files (.sorted.bam)
STEM=$(basename "$file" .bam)  

SOFTWARE=/datacommons/meyerlab/DARTseq/Bullseye/

####### list of options for parseBAM.pl #######
# This program will build a matrix of nucleotide count for every positions mapped in a bam file:
	# --mode : one of 'Bulk', 'SingleCell' or 'ExtractBarcodes'. Default to Bulk. 
	# --input: input bam file used to build matrix. Make sure files are coordinated sorted
	# --output: output file name, defaults to STDOUT
	# --removeDuplicates: To ignore reads marked as PCR or optical duplicates
    # --removeMultiMapped: To ignore multi mapped reads 
	# --verbose: display extra information
	# --cpu: number of threads to use for processing (only necessary if you are not using slurm)
	# --mem: available memory for sorting (M) (only necessary if you are not using slurm)
	# --minCoverage: minimum base coverage to output to final file (default = 1)
	# --filterBarcodes: only keep the barcodes included in the first column of a provided file. 
	# 			A second column with the number of reads for each barcode can be provided and used for filtering with the following options:
	# --MaxBarcodes: (sc) number of barcodes to process (defaults to all). Will stop reading provided barcode file after hitting this number, barcodes should be filered by mapped reads.
	# --readThreshold: (sc) minimum number of reads for barcodes to be considered
	# --Cell_ID_pattern: pattern for to use for single cell processing. use 10X, SMART or enter a sam tag RegEx for cell identification.
    # --exclude: bed file with coordinates for regions to be skipped. (eg. Ribosomal RNAs.) 
#######			end of options			#######
	

echo 'processing' $file

echo 'Started at' `date +"%D %H:%M"`

perl $SOFTWARE/parseBAM.pl --input $file --output $STEM.matrix --verbose --removeDuplicates --minCoverage 10

echo 'Ended at' `date +"%D %H:%M"`


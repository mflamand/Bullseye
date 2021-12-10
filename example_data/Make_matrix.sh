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
#--input: input bam file used to build matrix. Make sure files are coordinated sorted
#--output: output file name, defaults to STDOUT
#--filterBarcodes: When processing bulk data, only keep reads matching to specified barcodes
#--removeDuplicates: To ignore reads marked as PCR or optical duplicates
#--verbose: display extra information
#--cpu: number of threads to use for processing (only necessary if you are not using slurm)
#--mem: available memory for sorting (M) (only necessary if you are not using slurm)
#--minCoverage: minimum base coverage to output to final file (default = 1)
#--Nobgzip: This option will skip compression and indexing of final matrix
#######			end of options			#######

echo 'processing' $file

echo 'Started at' `date +"%D %H:%M"`

perl $SOFTWARE/parseBAM.pl --input $file --output $STEM.matrix --verbose --removeDuplicates --minCoverage 10

echo 'Ended at' `date +"%D %H:%M"`


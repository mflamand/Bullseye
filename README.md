# Bullseye
Bullseye is a pipeline for detection of RNA editing sites in DART-seq datasets. 

## Bullseye was initially branched from HyperTRIBE and modified to allow the following:

1. Enable setting detection of various RNA editing (TRIBE / DART). By default C-to-U transitions are detected, other types of editing events can be indicated on the command line. 
2. removed the use of a mySQL server.   
3. Support detection of editing in single cell sequencing datasets.
4. Improve speed through multicore processing.
5. Allows detection of editing site by comparison to genomic sequence or to a control file.

Bullseye is a set of Perl script, but make heavy use of Samtools for reading and indexing files. The requierements are:

- Perl > v.5.26.0 (untested on earlier version, but could work with minimal changes)
- Perl modules: 
  - MCE, multicore engine
  - Bio::DB::Fasta
  - Array::IntSpan
- Samtools (>v.1.10)
- Tabix 
- bedtools 

## Installation of prerequisite softwares
Prerequisites:

Perl > 5.26 Samtools Bedtools Tabix Perl modules MCE, Math::CDF and Bio::DB::Fasta

The easiest way to get everything up and running without root access is through conda.
create the conda environment with the provided yml file and activate the environment

	conda env create -f bullseye.yml

	conda activate Bullseye

#### then install Bio::DB::Fasta
XML::Parser is required for Bio::DB::Fasta, but causes problems with conda/cpanm. I opt to install it manually, providing the expat library path:

	wget http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz 
	tar -xf XML-Parser-2.46.tar.gz
	cd XML-Parser-2.46 
	perl Makefile.PL EXPATLIBPATH=$CONDA_PREFIX/lib EXPATINCPATH=$CONDA_PREFIX/include
	make
	make install

To install the remaining perl packages:

	cpanm Bio::DB::Fasta
	cpanm Text::NSP
	cpanm Array::IntSpan
	cpanm MCE

For the use of Bullseye in the detection of m6A in single cells, see: 

1. Tegowski, M. , Flamand, M.N., Meyer, K.M. scDART-seq reveals distinct m6 A signatures and mRNA methylation heterogeneity in single cells. Mol Cell 82, 868–878, February 17, 2022. https://doi.org/10.1016/j.molcel.2021.12.038 

For information of original HyperTRIBE pipeline, see : 

## HyperTRIBE

https://github.com/rosbashlab/HyperTRIBE

For more details please see:

1. Xu, W., Rahman, R., Rosbash, M. Mechanistic Implications of Enhanced Editing by a HyperTRIBE RNA-binding protein. RNA 24, 173-182 (2018). doi:10.1261/rna.064691.117

2. McMahon, A.C.,  Rahman, R., Jin, H., Shen, J.L., Fieldsend, A., Luo, W., Rosbash, M., TRIBE: Hijacking an RNA-Editing Enzyme to Identify Cell-Specific Targets of RNA-Binding Proteins. Cell 165, 742-753 (2016). doi: 10.1016/j.cell.2016.03.007.


## Running Bullseye

Bullseyes takes coordinate-sorted and indexed bam files as input and outputs a bed like file with a list of editing sites. Two scripts (parseBAM.pl and find_RNA_edit_sites.pl) are used to first parse individual bam files to generate coverage matrices, which are then compared to identify editing sites. 

Typically, raw fastq files are trimmed with flexbar and aligned to the genome with STAR. The use of other aligner should also be supported but was not tested. If PCR and optical duplicates are marked in the sam flag, they can be ignored during parsing 

bash script are provided and provide examples on how to run the perl script on an hpc cluster. for each perl script, use --help to get the complete list of options

1. Running parseBAM.pl:
	This first script will parse the aligned and sorted BAM files to output a tab delimited file with the count of each nucleotides at each position in the genome.
  	A good starting point for bulk DART-seq processing would be to run it as such, removing PCR duplicate reads, and keeping all positions covered by at least 10 reads:
  	
  		perl parseBAM.pl --input file.bam --output output.matrix --cpu 4 --minCoverage 10 --removeDuplicates
	
	Addional options are available and can be listed with the --help option	  
	
	Single cell dataset data can be processed using the same script using the --mode SingleCell option and indicating the type of barcode used with the -Cell_ID_pattern option. For example: 
		
		perl parseBAM.pl --mode SingleCell --Cell_ID_pattern 10X --input file.bam --output output.matrix --cpu 4

	Additionally, we provide a bash script (Make_matrix.sh) in the example data folder. This script provides a usage example of parseBAM.pl under the Slurm workload manager
  	
2. Running Find_edit_site.pl:
  This script compares editing in a matrix file from parseBAM.pl to a control matrix (Mettl3 KO or YTHmut-APOBEC1) or to the genomic sequence to identify m6A sites.
  Only the sites found in transcripts defined in a refFlat annotation file will be identified. RefFlat files can be downloaded from UCSC. an example file is provided in the example directory.
  
		perl Find_edit_site.pl --annotationFile $annotation_file \
		--EditedMatrix dart.matrix.gz \
		--controlMatrix control.matrix.gz \
		--minEdit 5 \ #minimal editing rates
		--maxEdit 90 \ #maximal editing rates
		--editFoldThreshold 1.5 \ # minimal editing ration over control sample
		--MinEditSites 3 \ #minimal number of mutations for detection of site
		--cpu 4 \
		--outfile output.bed \
		--fallback genome.fasta \ # fasta file of genome in case there was no coverage for a given position in the control file
		--verbose
	
	As before, additional options are available and can be listed with the --help option.

	The output file is a bed like file:
	the first 6 columns being the standard BED6 format: chr start end gene score strand
	file also contains additional columns:
		Control_edit_frequency: editing rate in control sample
		control_covorage: coverage in control samples
		dart_edit_frequency: editing rate in DART sample (same as score column)
		dart_coverage: coverage in control samples
		conversion_type: C2U for DART-seq

3. Filtering

	Beyond the core scripts listed above, I provide helper scripts to summarize and filter sites after their identification:
	
	1. summarize_sites.pl merges sites across output of multiple replicates. 

		For example, to merge the sites that were identified in a single sample, but against several control files (3), we can run:
		
			perl summarize_sites.pl \
			--MinRep 3 \	# minimum number of replicates for sites to be kept
			--mut 3 \	# minimum number of mutation for a site to be used in merging
			--repOnly \	# This option will merge files without adding coverage and average editing frequency. This is useful when merging the output from the same matrix against many control file
			sample1.*.bed > sample1.output.bed	# input bed files
			
		To then find common site between biological replicates, we can run:
		
			perl summarize_sites.pl \
			--MinRep 3 \
			--mut 3 \	
			sample1.output.bed sample2.output.bed sample3.output.bed > common_sites.min3rep.bed	
		
	2.  RAC filtering. For identification of m6A sites we may want to remove editing sites that are not found in the canonical RAC motif.
		RACfilter.sh is a bash script wrapper that uses bedtools getFasta and perl to fetch the sequence at each site and keep only those that are found in RAC.
		Before running, the bash script will need to be modified to include the path to a fasta file of the appropriate genome.
		we can run it as such :
			
			#first convert to bed6 if necessary
			perl -anE 'say join("\t", @F[0..5])' file.bed > out.bed
			bash RACfilter.sh -r *.bed # keep RAC sites only

> We provide example data and script for bulk DART-seq data on a single gene (APC).
> This scripts are for use on a SLURM based high performance cluster but can be adapted to be run directly : 

  1. In the directory containing the BAM files (WT_Soma_ctrl1.bam,WT_Soma_ctrl2.bam, Mettl3KO_Soma_ctrl1.bam, Mettl3KO_Soma_ctrl2.bam ), we first parse those files using Make_matrix.sh :

  			sbatch -a 1-4 Make_matrix.sh

  2. For each DART matrix we can find edit sites against each control file with the find_RNA_edit_sites.sh helper script :
		
			#for sample 1 (WT_Soma_Ctrl1_Apc.matrix.gz)
			sbatch -a 1-2 find_RNA_edit_sites.sh WT_Soma_Ctrl1_Apc.matrix.gz Mettl3KO_Soma_Ctrl1_Apc.matrix.gz Mettl3KO_Soma_Ctrl2_Apc.matrix.gz
			#for sample 1 (WT_Soma_Ctrl1_Apc.matrix.gz)
			sbatch -a 1-2 find_RNA_edit_sites.sh WT_Soma_Ctrl2_Apc.matrix.gz Mettl3KO_Soma_Ctrl1_Apc.matrix.gz Mettl3KO_Soma_Ctrl2_Apc.matrix.gz

  3. We can summarize and filter the sites : 

  to merge each replicates against both controls:

  			perl summarize_sites.pl --repOnly WT_Soma_Ctrl1_Apc.Mettl3KO_Soma_Ctrl1_Apc.bed WT_Soma_Ctrl1_Apc.Mettl3KO_Soma_Ctrl2_Apc.bed > WT_Soma_Ctrl1_Apc.bed
  			perl summarize_sites.pl --repOnly WT_Soma_Ctrl2_Apc.Mettl3KO_Soma_Ctrl1_Apc.bed WT_Soma_Ctrl2_Apc.Mettl3KO_Soma_Ctrl2_Apc.bed > WT_Soma_Ctrl2_Apc.bed

  Both replicates can then be merged with : 

  			perl summarize_sites.pl --minRep 2 WT_Soma_Ctrl1_Apc.bed WT_Soma_Ctrl2_Apc.bed > WT_Soma_Ctrl.bed


  The expected output can be found in the expected_output directory

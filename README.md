# Bullseye
 Bullseye is a pipeline for detection of RNA editing sites in DART-seq datasets. 

## Bullseye was initially branched and inspired from rosbashlab/HyperTRIBE and modified. The main differences are as follow:

1. Enable detection of any RNA editing types (TRIBE / DART). By default C-to-U transitions are detected, other types of editing events can be indicated on the command line. 
2. removed the use of a mySQL server. Coverage matrices are instead stored as Tabix indexed files.   
3. Support detection of editing sites in single cell sequencing datasets.
4. Used perl's MCE module for multicore processing.
5. Allows detection of editing site by comparison to genomic sequence or to a control file.

Bullseye is a set of Perl script, but make use of Samtools input files and for and compressing and indexing intermediart files. 

The Prerequisites are:

- Linux/Unix System: tested with CentOS8 and with Ubuntu 20.04 running in the Windows subsytem for Linux (WSL)
- Perl, tested with v.5.26.0 
- Perl modules: 
    - BioPerl
    - MCE, multicore engine
    - Array::IntSpan
    - Math::CDF
    - Text::NSP
    - Bio::DB::Fasta
- Samtools and Tabix (>v.1.10)
- bedtools 

## Installation of prerequisite softwares

The easiest way to get everything up and running without root access is through conda.
create the conda environment with the provided yml file and activate the environment

conda env create -f bullseye.yml

conda activate Bullseye
You will then need to install some perl modules, first start with MCE

cpanm MCE

#### then install Bio::DB::Fasta
XML::Parser is required for Bio::DB::Fasta, but causes problems with conda/cpanm. I opt to install it manually, providing the expat library path:

wget http://www.cpan.org/authors/id/T/TO/TODDR/XML-Parser-2.46.tar.gz 
tar -xf XML-Parser-2.46.tar.gz
cd XML-Parser-2.46 
perl Makefile.PL EXPATLIBPATH=$CONDA_PREFIX/lib EXPATINCPATH=$CONDA_PREFIX/include
make
make install

cpanm Bio::DB::Fasta
cpanm Text::NSP
cpanm Array::IntSpan
cpanm Math::CDF

## Running Bullseye

> Bullseyes takes coordinate-sorted and indexed bam files as input and outputs a bed like file with a list of editing sites. Two scripts (parseBAM.pl and find_RNA_edit_sites.pl) are used to first parse individual bam files to generate coverage matrices, which are then compared to identify editing sites. 

bash script are provided and provide examples on how to run the perl script on an hpc cluster. for each perl script, use --help to get the complete list of options

1. Running parseBAM.pl:

  perl parseBAM.pl --input file.bam 
  options:
  --output output.matrix          # outputs to STDOUT by default
  --cpu 4                         # for multicore processing
  --minCoverage 10                # minimum coverage for reads in output matrix
  --removeDuplicates              # ignore duplicate reads

  See Make_matrix.sh in the example data folder for and of how to use this parseBam.pl with Slurm workload manager

  The output will be a tab separated file containing coverage for each nucleotide at each position in the genome. The file is compressed with bgzip and indexed with Tabix

2. Running Find_edit_site.pl:
  This script compares editing in a matrix file from parseBAM.pl to a control matrix (Mettl3 KO or YTHmut-APOBEC1) or to the genomic sequence to identify m6A sites.
  Only the sites found in transcripts defined in a refFlat annotation file will be identified. RefFlat files can be downloaded from UCSC. an example file can be found here as well.

  perl Find_edit_site.pl --annotationFile $annotation_file \       # RefFlat file
	--EditedMatrix $dart.matrix \                                    # DART-seq matrix
	--controlMatrix $control.matrix \                                # Control matrix / alternatively, use the --genome genome.fasta to find sites against genome
	--minEdit 10 \                                                   # minimum editing percentage
	--maxEdit 80 \                                                   # maximum editing percentage
	--editFoldThreshold 1.5 \                                        # minimal fold editing over control file
	--MinEditSites 2 \                                               # minimal number of mutations to consider editing sites
	--cpu 4 \                                                        # for multicore processing 
	--outfile $BED/$dart.$control.txt \                              # output file, defaults to standard output
	--fallback $GENOME \                                             # optional genome fasta file. If coverage is unsufficient in the control matrix, the genomic sequence will be used instead
	--verbose

The output file is a bed like file (first 6 columns are bed6): #chr start end gene score strand Control_edit_frequency control_covorage dart_edit_frequency dart_coverage conversion_type score

3. Filtering
  Accessory bash scripts help with filtering of site across replicates.
  a. summarize_sites.pl merges sites across output of multiple replicates. Additional filtering can be done at this step, by default, the total coverage at a given site is added and the editing frequency is averaged:

  perl summarize_sites.pl --MinRep 3 \         # minimum number of replicates for sites to be kept
    --mut 2 \                                  # minimum number of mutation for a site to be used in merging
    --fold 1.5 \                               # minimum fold editing over control file for each site to be used in merging
    --repOnly \                                # This option will merge files without adding coverage and average editing frequency.
                                               # This is usefull when merging the output from the same matrix against many control file
    *.bed > output.bed                         # input bed files

  b.  RAC filtering. RACfilter.sh is a wrapper for bedtools getFasta that allows to keep only sites that are in the RAC consensus motif:
  
  perl -anE 'say join("\t", @F[0..5])' file.bed > out.bed   # convert to bed 6

  bash RACfilter.sh -r *.bed # keep RAC sites only

> For example, on a SLURM based cluster we can process example data as such: (here we are looking at sites on a single gene (Apc)):

  in directory with bam files (WT_Soma_ctrl1.bam,WT_Soma_ctrl2.bam, Mettl3KO_Soma_ctrl1.bam, Mettl3KO_Soma_ctrl2.bam ), we can first make matrix with :

   sbatch -a 1-4 Make_matrix.sh

  for each DART matrix we can find edit sites against each control file with :

   sbatch -a 1-2 find_RNA_edit_sites.sh WT_Soma_Ctrl1_Apc.matrix.gz Mettl3KO_Soma_Ctrl1_Apc.matrix.gz Mettl3KO_Soma_Ctrl2_Apc.matrix.gz

   sbatch -a 1-2 find_RNA_edit_sites.sh WT_Soma_Ctrl2_Apc.matrix.gz Mettl3KO_Soma_Ctrl1_Apc.matrix.gz Mettl3KO_Soma_Ctrl2_Apc.matrix.gz

  Sites can then be summarized with : 

  to merge each replicates against controls:

  perl summarize_sites.pl --repOnly WT_Soma_Ctrl1_Apc.Mettl3KO_Soma_Ctrl1_Apc.bed WT_Soma_Ctrl1_Apc.Mettl3KO_Soma_Ctrl2_Apc.bed > WT_Soma_Ctrl1_Apc.bed
  perl summarize_sites.pl --repOnly WT_Soma_Ctrl2_Apc.Mettl3KO_Soma_Ctrl1_Apc.bed WT_Soma_Ctrl2_Apc.Mettl3KO_Soma_Ctrl2_Apc.bed > WT_Soma_Ctrl2_Apc.bed

  Both replicates can then be merged with : 

  perl summarize_sites.pl --minRep 2 WT_Soma_Ctrl1_Apc.bed WT_Soma_Ctrl2_Apc.bed > WT_Soma_Ctrl.bed

  The expected output can be found in the expected_output directory



For more information on the original HyperTRIBE pipeline, you can see : 

## HyperTRIBE

https://github.com/rosbashlab/HyperTRIBE

For more details please see:

1. Xu, W., Rahman, R., Rosbash, M. Mechanistic Implications of Enhanced Editing by a HyperTRIBE RNA-binding protein. RNA 24, 173-182 (2018). doi:10.1261/rna.064691.117

2. McMahon, A.C.,  Rahman, R., Jin, H., Shen, J.L., Fieldsend, A., Luo, W., Rosbash, M., TRIBE: Hijacking an RNA-Editing Enzyme to Identify Cell-Specific Targets of RNA-Binding Proteins. Cell 165, 742-753 (2016). doi: 10.1016/j.cell.2016.03.007.

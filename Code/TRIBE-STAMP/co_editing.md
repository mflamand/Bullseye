# TRIBE-STAMP analysis pipeline

This pipeline is for the analysis of co-editing as described in the submitted manuscript Flamand MN, Ke K et al.

These script expand Bullseye to enable analysis of co-editing by APOBEC and ADAR fusion proteins at the single molecule level. 

For detail on the original Bullseye pipeline, please see: 

https://github.com/mflamand/Bullseye

Bullseyes takes coordinate-sorted and indexed bam files as input and outputs a bed like file with a list of editing sites. Two scripts (parseBAM.pl and find_RNA_edit_sites.pl) are used to first parse individual bam files to generate coverage matrices, which are then compared to identify editing sites.

Typically, raw fastq files are trimmed with flexbar and aligned to the genome with STAR. The use of other aligner should also be supported but was not tested. If PCR and optical duplicates are marked in the sam flag, they can be ignored during parsing

bash script are provided and provide examples on how to run the perl script on an hpc cluster. for each perl script, use --help to get the complete list of options

### Intructions for TRIBE-STAMP

1. parsing BAM file with parseBAM.pl:
	This first script will parse the aligned and sorted BAM files to output a tab delimited file with the count of each nucleotides at each position in the genome.
    The --removeDuplicates and --removeMultiMapped options can be used to ignore reads that are PCR duplicates or mapped to multiple location in the genome

    ```{bash}
            perl parseBAM.pl --input input.bam \
            --output output.matrix \
            --verbose \
            --removeDuplicates \
            --removeMultiMapped \
            --minCoverage 10
    ```

2. Identification of editing sites in the transcriptome with Find_edit_site.pl:
  
	This script compares editing in a matrix file from parseBAM.pl to a control matrix (APOBEC alone or ADAR alone expressing cells),

	Editing site will be identified only in regions specified in an annotation file. For this a GenePrediction file (refFlat) is provided. This allows annotation of each sites to a feature of the transcriptome. Optionally, the flags "--intron" and "--extUTR size" can be added to identify sites in the introns of transcripts and in an extended region after the annotated 3'UTRs (protein coding genes only).

	A refFlat files can be downloaded from UCSC table browser, in the refFlat table of the NCBI RefSeq track.

    For example, to identify A-to-I editing sites we run :

    ```{bash}
    perl Find_edit_site.pl --annotationFile ucsc_genepred.refFlat \
	--EditedMatrix input.matrix.gz \
	--controlMatrix ADAR_alone.matrix.gz \
	--editType A2I \
	--minEdit 1 \
	--maxEdit 99 \
	--editFoldThreshold 3 \
	--EditedMinCoverage 50 \
	--ControlMinCoverage 10 \
	--MinEditSites 3 \
	--outfile output.A2I.bed \
	--fallback hg38.fasta \
	--filterBed dbSNP153_v36.bed \
	--extUTR \
	--verbose
    ```

    and for C-to-U editing sites:

    ```{bash}
    perl Find_edit_site.pl --annotationFile ucsc_genepred.refFlat \
	--EditedMatrix input.matrix.gz \
	--controlMatrix APOBEC1_alone.matrix.gz \
	--editType C2U \
	--minEdit 1 \
	--maxEdit 99 \
	--editFoldThreshold 3 \
	--EditedMinCoverage 50 \
	--ControlMinCoverage 10 \
	--MinEditSites 3 \
	--outfile output.C2U.bed \
	--fallback hg38.fasta \
	--filterBed dbSNP153_v36.bed \
	--extUTR \
	--verbose
    ```

    By using the --filterBed option, we are removing any known SNP sites found in the provided bed file

    --falback is used detect sites when there is insufficient coverage in the TRIBE-STAMP samples but not in the control samples (ADAR or APOBEC1 alone). In this case, the mutations rates are compared to the reference sequence found in the provided fasta file.

3. Filter for high confidence sites

	1. Merges sites for each samples with summarize_sites.pl.
    Each sample is compared to at least 2 negative controls, we only keep sites that are found against all negative controls.

    ```{bash}
    perl summarize_sites.pl \
        --MinRep 2 \
        --repOnly \
        sample1.control1.bed sample1.control2.bed > sample1.output.bed	
    ```
    by using the --repOnly, we are merging the sites without changing summarizing the number of mutation, coverage and editing rates

    We can then identify all sites found in multiple biological replicates using the same script:

    For example to find sites identified in at least 4 replicates in DF1_ADAR expressing cells we could run:

    ```{bash}	
			perl summarize_sites.pl \
			--MinRep 4 
			DF1_ADAR_sample*.output.bed > DF1_ADAR.min4.bed 
	```

    if input files are named DF1_ADAR_sample1.output.bed DF1_ADAR_sample2.output.bed ... etc

4. Having identified the editing sites for each DF fusion protein, we can now look identify pairs of editing sites within a specified windows and look in all reads that cover both sites for co-editing using co_editing.pl

    For example to look at co-editing in cells expressing DF1-ADAR and DF2-APOBEC1, we can run:

    ```{bash}
    perl co_editing.pl --bed1 DF1_C2U.bed \
        --bed2 DF2_APOBEC1.bed \
        --bam DF1_ADAR.DF2_APOBEC.rep1.bam \
        --outfile DF1_ADAR.DF2_APOBEC.rep1.C2U.A2I.tsv \
        --refFlat ucsc_genepred.refFlat
        --minCov 20 \
        --distance 150 \
        --removeDup \
        --removeMultiMapped \
        --PermutationTest 10000 \
        --FilterNeg
    ```

    --removeDup and --removeMultiMapped are as before used to ignore duplicate and multi mapped reads

    --distance is used to indicate the maximum distance between each editing sites. Generally this would be the sequencing length or less

    --refFlat: An annotation file can be provided to allow the detection of co-editing on pairs of sites that are on different exons, but found within the distance window. As of now, only sites separated by a single intron will be probed 

    --PermutationTest N, is to perform a permutation test on the distribution of edited reads (Non-edited, edited at site1 only, edited at site2 only, edited at both sites). This option can be ommited.

    --FilterNeg is used to remove pairs with no editing from the output file. If the option is ommited, editing at all pairs will be reported, even if no editing is detected in this particular biological replicate.

    The output of this script will be a tab separated file with the following information: 

    /chr/site1_position/site2_position/strand/coverage_site1/edit_rate_site1/coverage_site2/edit_rate_site2/reads_covering_both_sites/non_edited_reads/ reads_edited_site1_alone/reads_edited_site2_alone/reads_edited_both_sites/co_editing_ration/Fisher_p_value/Permutation_p_value/Distance/GeneID);


### Intructions for quantification of editing in RIP-seq and polysome-associated RNAs

Steps 1 through 3 are performed as for TRIBE-STAMP, to identify editing sites in cells expressing each fusion protein.

4. Quantify editing at all called sites in each biological replicates with quant_sites.pl:

    From the list of sites, we look in the matrix generated in step 1 and report the coverage, number of mutation and editing rates at each called site

    ```{bash}
    perl Quantify_sites.pl --ClusterFile DF1_ADAR.bed \
    --EditedMatrix DF1_ADAR.DF2_APOBEC1_rep1.matrix.gz \
    --editType A2I \
    --outfile DF1_ADAR.DF2_APOBEC1_rep1.A2I.txt \
    --verbose

   perl Quantify_sites.pl --ClusterFile DF2_APOBEC1.bed \
    --EditedMatrix DF1_ADAR.DF2_APOBEC1_rep1.matrix.gz \
    --editType C2U \
    --outfile DF1_ADAR.DF2_APOBEC1_rep1.A2I.txt \
    --verbose
    ```

    We can then gather the information in a table containing the coverage and the number of mutation at each site using gather_score.pl:

    ```{bash}
     perl gather_score.pl --coverage --mutations --outfile DF1_ADAR.txt DF1_ADAR*.txt # with all files generated in previous step as input

     perl gather_score.pl --coverage --mutations --outfile DF2_APOBEC.txt *DF2_APOBEC1*.txt 
    ```

    The output of this script will be 2 matrix containing a single line for each site, and the coverage or number of mutation for each biological replicate for the files provided. (e.g. DF1_ADAR.coverage.txt and DF1_ADAR.mutations.txt are generated when providing the --coverage and --mutations options)

    These data are then imported in R to fit the mutation rates (mut/cov) in a quasibinomial GLM for each editing sites.
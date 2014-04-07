Ver0.01 - April 2014 
DeNovoGear - paired sample caling with Beta Binomial genotype Likelihoods
Avinash Ramu, Don Conrad, Conrad Lab, WUSTL.

Introduction
Traditionally variant callers use a binomial assumption for allele counts coming from a sequencing experiment. This assumption shows good performance for Whole Genome data but can be improved for Whole Exome data where other biases come into play. We have shown that the beta-binomial model captures the variance in the Whole Exome data better and calling variants based on this assumption controls sensitivity and specificity better. This is a variant caller which compares two samples hence paired-sample and tries to find mutations that are in one sample but not in the other. Currently this has been implemented only for SNV's.

INSTALL
For installation instructions, see INSTALL

Sample Usage
denovogear-bb bb_paired --ped example.ped --normal_bam normal_bam --tumor_bam tumor_bam --ref reference.fa --pp_cutoff 0.1 --rd_cutoff 10 > dng-bbOP.txt 


Arguments
Required arguments.
--ped  Tab delimited file, specifies the normal and tumor sample names. It contains just one line of the format 
      "ExperimentID    TumorID	    NormalID   0  0  0"
--normal_bam  BAM file of the normal sample.
--tumor_bam   BAM file of the tumor sample.
--ref  Reference genome the BAM was mapped against. This is in the FASTA format, an index will be created if it doesn't exist already.

Optional Arguments
--pp_cutoff  posterior probability cutoff for DNG op. It has to be a real number between 0 and 1. For example if PP_cutoff is set to 0.9, all calls with posterior probability >= 0.9 will be printed to screen. I would suggest setting it to somewhere between 0.1 and 0.5 and looking at the number of calls that you get and if that matches your theoretical expectation. default = 0.001
--rd_cutoff  read_depth cutoff for the calls to be made. For example if set to 10, only sites that have a read depth of >=10 in each sample will be considered for calling, the same sites are also used for estimating the alpha and beta parameters of the model. default = 10.
--skip_mpileup skip the mpileup step and only perform calling. [This saves some time, for development usage]


Example Output
DENOVO-PAIR-SNP TUMOR_ID: tumor NORMAL_ID: normal chr: chr1 pos: 47604741 ref: C alt: N maxlike_null: 6.30957e-09 pp_null: 0.864664 tgt_null(normal/tumor): CC/CC maxlike_dnm: 1e-09 pp_dnm: 0.135336 tgt_dnm(normal/tumor): CG/CC READ_DEPTH tumor: 177 normal: 1194 null_snpcode: 1 dnm_snpcode: 3 tumor_Alt_RD: 0 normal_Alt_RD: 2

Explanation of output columns
TUMOR_ID - ID of the tumor sample.
NORMAL_ID - ID of the normal sample.
chr - chromosome 
pos - position
ref - reference Allele
alt - alternate Allele. If set to 'N' it means only the ref allele was observed or more than two alleles observed at the site.
maxlike_null - maximum likelihood of the NULL genotype configurations.
pp_null - posterior probability of the NULL genotype configuration i.e tgt_null.
tgt_null - genotype configuration under the NULL model.
maxlike_dnm - maximum likelihood of the DENOVO genotype configurations.
pp_dnm - posterior probability of the DENOVO genotype configuration, this is the field used to rank the denovo calls, values close to 1 indicate high probability of observing a denovo mutation at this site.
tgt_dnm -  DENOVO genotype configuration.
READ_DEPTH tumor, normal - read depth of tumor and normal samples.
null_snpcode, dnm_snpcode - codes to describe tgt_null and tgt_dnm columns, 0 = hom/hom 1 = hom/het 2 = het/hom 3 = het/het. Depending on the experiment you might choose to filter out calls based on a particular configuration.
tumor_Alt_RD = number of reads with the alternate allele in the tumor sample
normal_Alt_RD = number of reads with the alternate allele in the normal sample

Notes
1. The confidence of the denovo model is given by the pp_dnm field. This is the posterior probabily of observing a denovo mutation at a given site i.e values closer to 1 have a higher probability of being a true positive.
2. This code makes use of SamTools libraries for parsing files. We thank the authors of SamTools for having made their code open-source.

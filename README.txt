Ver1.0 - June 7 2013
DeNovoGear with Beta Binomial genotype Likelihoods
Avinash Ramu, Conrad Lab, WUSTL.

Sample Usage
/home/comp/exlab/aramu/files/DeNovoGear/BetaBinomial/BetaBinomial_paired/build/src/denovogear-bb bb_paired --ped temp.ped --normal_bam /home/comp/tdlab/ahughes/DeNovoGear_Test_Data/PAMBCT_Non_Cancer.bam --tumor_bam /home/comp/tdlab/ahughes/DeNovoGear_Test_Data/PAMBCT_Leukemia.bam --ref /home/comp/exlab/aramu/files/20110915_CEUtrio/WEx/ref/hg19.fa --PP_cutoff 0.1 --RD_cutoff 10 > drew.dngOP.txt 


Arguments
--ped  Tab delimited file, specifies the normal and tumor sample names. It contains just one line of the format 
      "ExperimentID    TumorID	    NormalID   0  0  0"
--normal_bam  BAM file of the normal sample.
--tumor_bam   BAM file of the tumor sample.
--ref  Reference genome the BAM was mapped against.
--PP_cutoff  posterior probability cutoff for DNG op. It has to be a real number between 0 and 1. For example if PP_cutoff is set to 0.9, all calls with posterior probability >= 0.9 will be printed to screen. I would suggest setting it to somewhere between 0.1 and 0.5 and looking at the number of calls that you get and if that matches your theoretical expectation. default = 0.001
--RD_cutoff  read_depth cutoff for the calls to be made. For example if set to 10, only sites that have a read depth of >=10 in each sample will be considered for calling, the same sites are also used for estimating the alpha and beta parameters of the model. default = 10.
--skip_mpileup skip the mpileup step and only perform calling. [for development usage]


Example Output
DENOVO-PAIR-SNP TUMOR_ID: tumor NORMAL_ID: normal chr: chr1 pos: 47604741 ref: C alt: N maxlike_null: 6.30957e-09 pp_null: 0.864664 tgt_null(normal/tumor): CC/CC maxlike_dnm: 1e-09 pp_dnm: 0.135336 tgt_dnm(normal/tumor): CG/CC READ_DEPTH tumor: 177 normal: 1194 MAPPING_QUALITY tumor: -1 normal: -1 null_snpcode: 1 dnm_snpcode: 3 tumor_Alt_RD: 0 normal_Alt_RD: 2

Output Columns
TUMOR_ID - ID of the tumor sample.
NORMAL_ID - ID of the normal sample.
chr - chromosome 
pos - position
ref - ref Allele
alt - alternate Allele. If set to 'N' more than two alleles observed at the site.
maxlike_null - maximum likelihood under the NULL model.
pp_null - posterior probability of the NULL model.
tgt_null - genotype configuration under the NULL model.
maxlike_dnm - maximum likelihood under the DENOVO model.
tgt_dnm -  genotype configuration under the DENOVO model.
READ_DEPTH tumor, normal - read depth of tumor and normal samples.
MAPPING QUALITY tumor, normal - IGNORE for now
null_snpcode, dnm_snpcode - codes to describe tgt_null and tgt_dnm columns, 0 = hom/hom 1 = hom/het 2 = het/hom 3 = het/het
tumor_Alt_RD = number of Reads with the alternate allele in the tumor sample
normal_Alt_RD = number of Reads with the alternate allele in the normal sample


The confidence of the denovo model is given by the pp_dnm field.

This code makes use of SamTools libraries for parsing files.

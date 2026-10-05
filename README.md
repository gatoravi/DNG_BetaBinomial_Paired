# DeNovoGear-BB: paired-sample mutation calling with beta-binomial likelihoods

A variant caller that compares two sequencing samples from the same individual,
such as a tumor and its matched normal, and finds single-nucleotide mutations
present in one sample but not the other. It extends
[DeNovoGear](https://github.com/ultimatesource/denovogear) by replacing the
usual binomial read-count model with a **beta-binomial** model.

*Avinash Ramu and Don Conrad, Conrad Lab, Washington University in St. Louis
(v0.01, April 2014).*

## Why beta-binomial?

Most variant callers assume that reads supporting each allele follow a
binomial distribution. That works well for whole-genome sequencing, but in
whole-exome data, capture efficiency, PCR duplication and other biases add
extra variance (overdispersion). A beta-binomial model captures that
overdispersion, which gives better control of sensitivity and specificity
when calling mutations from exomes.

## How it works

1. Runs `samtools mpileup` on the two BAM files to get per-site read counts.
2. At sites passing the read-depth cutoff, fits the beta-binomial α and β
   parameters by maximum likelihood (`src/BBFit.R`, called from C++ through
   RInside/Rcpp).
3. For every site, compares the most likely "no mutation" genotype
   configuration with the most likely "new mutation" configuration and reports
   the posterior probability of a mutation (`pp_dnm`).

Only SNVs are supported.

## Install

Requires CMake, R (3.0 or later) and the R packages Rcpp, RcppClassic and
RInside. A patched copy of samtools 0.1.18 is bundled.

```bash
git clone https://github.com/gatoravi/DNG_BetaBinomial_Paired.git
cd DNG_BetaBinomial_Paired

# Build the bundled samtools libraries
(cd src/samtools-0.1.18/bcftools && make)
(cd src/samtools-0.1.18 && make lib)

# The caller looks for the R fitting script in ~/src/
mkdir -p ~/src && cp src/BBFit.R ~/src/

mkdir build && cd build
cmake ..
make
make install
```

See [INSTALL](INSTALL) for details.

## Usage

```bash
denovogear-bb bb_paired \
    --ped example.ped \
    --normal_bam normal.bam \
    --tumor_bam tumor.bam \
    --ref reference.fa \
    --pp_cutoff 0.1 --rd_cutoff 10 > dng-bb_output.txt
```

`src/callDNGBB.sh` has another example command.

### Arguments

| Argument | Required | Description |
| --- | --- | --- |
| `--ped` | yes | One-line tab-delimited file naming the samples: `ExperimentID  TumorID  NormalID  0  0  0` |
| `--normal_bam` | yes | BAM file of the normal sample |
| `--tumor_bam` | yes | BAM file of the tumor sample |
| `--ref` | yes | Reference FASTA the BAMs were aligned to (indexed if needed) |
| `--pp_cutoff` | no | Report calls with `pp_dnm` at or above this value (default 0.001; 0.1 to 0.5 is a sensible starting range) |
| `--rd_cutoff` | no | Minimum read depth in each sample, used both for calling and for fitting α and β (default 10) |
| `--skip_mpileup` | no | Skip the mpileup step and only call (for development) |

## Output

One line per candidate site, for example:

```
DENOVO-PAIR-SNP TUMOR_ID: tumor NORMAL_ID: normal chr: chr1 pos: 47604741 ref: C alt: N maxlike_null: 6.30957e-09 pp_null: 0.864664 tgt_null(normal/tumor): CC/CC maxlike_dnm: 1e-09 pp_dnm: 0.135336 tgt_dnm(normal/tumor): CG/CC READ_DEPTH tumor: 177 normal: 1194 null_snpcode: 1 dnm_snpcode: 3 tumor_Alt_RD: 0 normal_Alt_RD: 2
```

| Field | Meaning |
| --- | --- |
| `TUMOR_ID`, `NORMAL_ID` | Sample IDs |
| `chr`, `pos` | Genomic position |
| `ref`, `alt` | Reference and alternate allele (`N` means only the reference allele, or more than two alleles, was seen) |
| `maxlike_null`, `pp_null`, `tgt_null` | Likelihood, posterior probability and genotypes (normal/tumor) of the best no-mutation configuration |
| `maxlike_dnm`, `pp_dnm`, `tgt_dnm` | The same for the best mutation configuration. **`pp_dnm` is the score to rank calls by**; values near 1 indicate a likely true mutation |
| `READ_DEPTH tumor`, `normal` | Read depth in each sample |
| `null_snpcode`, `dnm_snpcode` | Configuration codes: 0 = hom/hom, 1 = hom/het, 2 = het/hom, 3 = het/het |
| `tumor_Alt_RD`, `normal_Alt_RD` | Reads supporting the alternate allele in each sample |

## Related

- [DeNovoGear](https://github.com/ultimatesource/denovogear): the original
  trio and paired-sample de novo mutation caller.
- [denovogear-bb](https://github.com/gatoravi/denovogear-bb): earlier
  development snapshot of this code.

## Acknowledgements

Uses the samtools libraries for reading BAM and BCF files, and the newmat
matrix library.

## License

MIT. See [COPYING](COPYING).

# Population structure, divergence, and introgression in _Heliconius_

This repository contains the main analysis scripts used to characterize population
structure, divergence history, and introgression in the comimetic butterflies
_Heliconius melpomene_ and _Heliconius erato_ and their close relatives.

The repository is intended as a publication-code archive. It is not a fully
automated workflow and is not designed to run from beginning to end without
editing paths, sample lists, software environments, and cluster-specific
settings. Instead, it documents the core commands, parameter choices, and custom
analysis code used in the study, with local machine paths replaced by
placeholders.

## Repository contents

```text
metadata/
  erato_metadata.tsv
  melp_metadata.tsv

scripts/
  0_raw_read_trimming/
  1.0_mito_assembly/
  1.1_wgs_alignment/
  1.2_wgs_genotyping/
  2.1_accessible_genome/
  3.1_vcf_processing/
  3.2_genotyping_outgroups/
  4.1_eigenstrat/
  4.2_admixture/
  4.3_feems/
  4.4_popgen_stats/
  4.5_psi/
  4.6_allele_sharing/
  4.7_smcpp/
  4.8_fstatistics/
  5.1_mitotree/
  5.2_njtree/
  5.3_treemix/
  6.1_bpp/
  7.1_introgressed_haplotypes/
```

The `metadata/` directory contains sample metadata used to organize the two
main datasets. The `scripts/` directory is organized in the approximate order of
the analyses in the paper.

## What is not included

This repository does not include raw sequencing reads, large alignment files,
large VCF/BCF files, intermediate genomic outputs, or final figure panels. Those
files should be obtained from the public repositories and accession numbers
listed in the manuscript data availability statement.

The scripts assume that large input files such as raw reads, reference genomes,
indexed BAM files, VCF files, BED masks, and sample lists already exist in the
user's working environment.

## Analysis overview

| Directory | Purpose |
| --- | --- |
| `0_raw_read_trimming` | Trim paired-end reads with `fastp`. |
| `1.0_mito_assembly` | Assemble mitochondrial genomes and screen assemblies with BLAST. |
| `1.1_wgs_alignment` | Align whole-genome reads and mark duplicates. |
| `1.2_wgs_genotyping` | Generate GVCFs and jointly genotype samples with GATK. |
| `2.1_accessible_genome` | Build accessible-genome masks from depth, soft clipping, and repeats. |
| `3.1_vcf_processing` | Split, subset, filter, and prepare high-confidence VCFs. |
| `3.2_genotyping_outgroups` | Genotype outgroup samples for downstream polarization and comparison. |
| `4.1_eigenstrat` | Convert VCFs, filter sites, run PCA, and prepare EIGENSTRAT files. |
| `4.2_admixture` | Run ADMIXTURE analyses. |
| `4.3_feems` | Test FEEMS parameters and run final FEEMS analyses. |
| `4.4_popgen_stats` | Estimate windowed population-genetic statistics. |
| `4.5_psi` | Polarize alleles and calculate directionality statistics. |
| `4.6_allele_sharing` | Rarefy groups and summarize shared/private derived alleles. |
| `4.7_smcpp` | Convert VCFs and test SMC++ parameter combinations. |
| `4.8_fstatistics` | Precompute f2 blocks, run f4 tests, and search/plot qpGraphs. |
| `5.1_mitotree` | Align mitochondrial assemblies and infer mitochondrial trees. |
| `5.2_njtree` | Build nuclear neighbor-joining trees with leave-one-window-out support. |
| `5.3_treemix` | Convert PLINK files to TreeMix input and run TreeMix. |
| `6.1_bpp` | Generate loci and run BPP topology, MSCM, MSCI, and final model analyses. |
| `7.1_introgressed_haplotypes` | Impute genotypes, run SPrime, and score candidate introgressed haplotypes. |

## How to use these scripts

Each script begins with a short header describing its purpose, expected input,
expected output, and required software. Most user-specific paths are written as
placeholders such as:

```text
<path/to/input.vcf.gz>
<path/to/output_dir>
<path/to/sample_list.txt>
```

To reuse a script, replace the placeholders with paths in your own working
directory, confirm that the sample names match your metadata and sample lists,
and run the script in the relevant software environment.

These scripts are intentionally lean. They preserve the core commands and
analysis parameters used in the study, but they omit local scheduling details,
temporary logging, and machine-specific paths.

## Adapting to another study system

The directory structure can be reused for another population-genomic study. In
that case, the main changes will usually be:

1. Replace sample metadata and sample lists.
2. Replace reference genome paths and chromosome names.
3. Replace population, species, and geography labels.
4. Reconsider thresholds that depend on depth, sample size, or genome assembly.
5. Keep analysis parameters that represent deliberate scientific choices, such
   as replicate numbers, window sizes, missingness thresholds, and model
   settings.

Some scripts in `6.1_bpp/` and `7.1_introgressed_haplotypes/` are more
analysis-specific because they encode model topologies or focal comparison
groups from this study. They can still be reused, but those labels and model
definitions should be checked carefully before applying them to a new system.

## Software

The analyses use a mixture of command-line population-genomic software, Python,
and R. Major tools include:

- `fastp`
- `bwa-mem2`
- `samtools`
- `bcftools`
- `bedtools`
- `GATK`
- `mosdepth`
- `PLINK`
- `KING`
- `EIGENSOFT`
- `ADMIXTURE`
- `FEEMS`
- `pixy`
- `scikit-allel`
- `SMC++`
- `admixtools`
- `MAFFT`
- `trimAl`
- `IQ-TREE2`
- `BEAST`
- `TreeAnnotator`
- `TreeMix`
- `BPP`
- Python packages including `numpy`, `pandas`, and `scikit-allel`
- R packages including `tidyverse`, `admixtools`, `ape`, and related plotting
  packages

Exact software versions should be reported in the manuscript methods or in a
separate `software_versions.md` file before final release.

## Data availability

Raw sequencing reads and large derived genomic files are not stored in this
repository. Public accessions for previously published and newly generated data
should be provided in the manuscript data availability statement and, where
appropriate, in the metadata files in `metadata/`.

Before publication, all newly generated sequencing data should have stable public
accessions or a clearly stated access policy.

## Code availability

The release version of this repository should be archived in a DOI-minting
repository such as Zenodo or Code Ocean before publication. The manuscript code
availability statement should cite the archived version, not only the live
GitHub URL.

Example wording to update before submission:

```text
Code used for the analyses is available at GitHub: <repository_url>. The
version associated with this manuscript is archived at Zenodo: <DOI>.
```

## Citation

If using this repository, cite the associated manuscript and the archived code
release DOI once available.

## License

A license should be added before public release. Without a license, others can
view the code but do not have clear permission to reuse, modify, or redistribute
it. For a publication-code repository, a common choice is the MIT License for
code. If metadata are intended to be reused independently, a data license such
as CC BY 4.0 can also be considered.

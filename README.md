# Population structures and incipient speciation in comimetic _Heliconius_ butterflies

This repository documents the analytical pipeline behind the paper, tentatively
titled "Population histories and incipient speciation in mimetic butterflies",
including scripts used to characterize population structure, demographic
histories, and speciation histories of the comimetic butterflies _Heliconius
melpomene_ and _H. erato_, together with their close relatives, which are
incipient species with continued gene flow.

Analysis scripts are provided in `scripts/`, roughly following the order of the
analyses in the paper, and raw data or derived genomic files can be found at
`<data_repository_or_accessions>`.

These scripts document the parameters used in the main analyses, with directory as place holders, 
To reuse a script, replace these placeholders with the correct paths, check that the sample names match your metadata, and run the script in
the relevant software environment. Paths are written as placeholders, for example:

```text
<path/to/input.vcf.gz>
<path/to/output_dir>
<path/to/sample_list.txt>
```

## Analysis overview

| Directory | Purpose |
| --- | --- |
| `0_raw_read_trimming` | Trim paired-end raw reads with `fastp`. |
| `1.0_mito_assembly` | Assemble mitochondrial genomes and screen assembled contigs against custom BLAST databases. |
| `1.1_wgs_alignment` | Align whole-genome reads to the reference genome and mark duplicate reads. |
| `1.2_wgs_genotyping` | Call per-sample GVCFs, import them into GenomicsDB, and jointly genotype the dataset. |
| `2.1_accessible_genome` | Define accessible genome masks using depth, soft-clipping, and repeat filters. |
| `3.1_vcf_processing` | Split, subset, and filter VCFs into the high-confidence datasets used downstream. |
| `3.2_genotyping_outgroups` | Genotype outgroup samples for polarization and comparative analyses. |
| `4.1_eigenstrat` | Convert VCFs, rename contigs, filter sites, run PCA, and prepare EIGENSTRAT files. |
| `4.2_admixture` | Run ADMIXTURE across parameter settings. |
| `4.3_feems` | Test FEEMS parameters and run final spatial population-structure analyses. |
| `4.4_popgen_stats` | Calculate windowed population-genetic statistics, including FST, DXY, pi, and Tajima's D. |
| `4.5_psi` | Polarize ancestral alleles and calculate directionality statistics. |
| `4.6_allele_sharing` | Rarefy sample groups and summarize private/shared derived allele counts. |
| `4.7_smcpp` | Convert VCFs for SMC++ and test demographic parameter settings. |
| `4.8_fstatistics` | Precompute f2 blocks, calculate f4 statistics, and run qpGraph analyses. |
| `5.1_mitotree` | Align mitochondrial assemblies, trim alignments, and infer mitochondrial trees. |
| `5.2_njtree` | Build nuclear neighbor-joining trees with leave-one-window-out support. |
| `5.3_treemix` | Convert allele-frequency files and run TreeMix. |
| `6.1_bpp` | Generate loci and run BPP topology, MSCM, MSCI, and final parameter-estimation analyses. |
| `7.1_introgressed_haplotypes` | Impute genotypes, run SPrime, and score candidate introgressed haplotypes. |

## Software

The analyses use command-line genomics tools, Python, and R. Major software
includes:

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
- R packages including `tidyverse`, `admixtools`, `ape`, and plotting packages

Exact software versions should be reported with the final archived version of
the code.

## License

This repository is released under the MIT License.

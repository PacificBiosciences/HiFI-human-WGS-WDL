# tertiary.wdl analysis workflow

This is a simple, opinionated subworkflow for tertiary analysis in rare disease research.  It starts with small variants and structural variants in VCF format, filters to remove variants that are common in the population, annotates with functional impact, and then prioritizes based on the predicted impact on the gene and the gene's relevance to the phenotype.  It has been designed for ~30x WGS HiFi for the proband and ~10-30x WGS HiFi for parents (optional).

## Inputs

- Small variants and structural variants are provided to this workflow in VCF format.  If multiple family members have been sequenced, they are provided as a single joint-called VCF per variant type per family.  If only the proband has been sequenced, the VCFs are provided for the proband only.
- We generate a pedigree describing sample relationships and phenotype status, based on the input provided to the entrypoint workflow.  In the case of a singleton, the pedigree is a single row.
- Using the comma-delimited list of HPO terms provided to the entrypoint workflow, we generate a Phenotype Rank (Phrank) lookup table, a simple two column lookup table mapping gene symbols to Phrank score.  Phrank scores are positive real numbers (or null) such that higher scores indicate a gene is more likely to be relevant to the phenotypes.  The Phrank lookup is used to prioritize variants based on the predicted impact on the gene and the gene's relevance to the phenotype.  Phrank scores are not normalized, and providing more phenotypes for a sample will result in a higher maximum Phrank score.
- Reference data is provided by the [`ref_map_file`](./ref_map.md) input.  This workflow is currently only compatible with the GRCh38 human reference.
- Population data, other supplemental data, and allele thresholds are provided by the [`tertiary_map_file`](./tertiary_map.md) input.  We provide a version of this file that uses population data from [gnomAD v4.1](https://gnomad.broadinstitute.org/news/2024-05-gnomad-v4-1-updates/) and [CoLoRSdb](https://colorsdb.org) v1.2.0 [<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14814308.svg" alt="10.5281/zenodo.14814308">](https://zenodo.org/records/14814308). We provide the ability to tweak the allele thresholds, but the default values are recommended, as increasing these will result in much higher resource usage.

## Process

### Small variants

We use [`slivar`](https://github.com/brentp/slivar) and [`bcftools csq`](https://samtools.github.io/bcftools/howtos/csq-calling.html) to filter and annotate small variants, and to identify potential compound heterozygous ("comphet") candidate pairs.  Slivar uses variant annotations stored in "gnotate" databases.  We use the following steps (nb: some steps performed within the same command).

1. Ignore variants with a non-passing `FILTER` value.
2. Ignore variants that are present at > 3% (`slivar_max_af`) in any of the population datasets.
3. Ignore variants with more than 4 homozygous alternate ("homalt") calls (`slivar_max_nhomalt`) in any of the population datasets.  For the purposes of this tool, we count hemizygous ("hemialt") calls on the X chromosome as homalt.
4. To be tagged as a potential "dominant" variant, the site must be high quality[^1] in all relevant samples, present as homref in all unaffected samples, present as homalt or hetalt in all affected samples, and have allele count < 4 (`slivar_max_ac`) in the population datasets.
5. To be tagged as a potential "recessive" variant, the site must be high quality[^1] in all relevant samples, present as homalt or hemi in all affected samples, and present as homref or hetalt in all unaffected samples.
6. To be tagged in comphet analysis, the site must be have GQ > 5 (`slivar_min_gq`) and present as hetalt in all affected samples.
7. All remaining "tagged" variants are annotated with predicted impact using Ensembl GFF3 gene set and `bcftools csq`.  This annotated VCF is provided for downstream analysis.
8. All variants considered for comphet analysis with high potential impacts[^2] are considered in pairs.  If the pair of variants are shown to be _in cis_ according to HiPhase phasing, they are rejected.  The passing pairs are stored in a second VCF for downstream analysis.

We use [`slivar tsv`](https://github.com/brentp/slivar/wiki/tsv:-creating-a-spreadsheet-from-a-filtered-VCF) to produce TSVs from the VCFs generated above.  These TSVs have many of the relevant fields from the VCF, as well as:

- Clinvar annotations for the gene
- gnomAD [loss-of-function tolerance metrics](https://gnomad.broadinstitute.org/downloads#v2-lof-curation-results)
- Phrank scores for the gene

### Structural variants

We use [`svpack`](https://github.com/PacificBiosciences/svpack) to filter and annotate SVs, with the following steps.

1. Remove variants with a non-passing `FILTER` value.
2. Remove variants < 50bp
3. Remove variants that match any existing variant in: gnomAD v4.1 (n=) or CoLoRSdb (n=).  In this case, "match" means that the variant is the same type, the difference in position is <= 100bp, and the difference in size is <= 100bp.
4. Annotate `INFO/BCSQ` with predicted impact using Ensembl GFF3 gene set.
5. Annotate `INFO/homalt` and `INFO/hetalt` with the names of samples in this cohort that have the variant in homozygous or heterozygous form, respectively.

We use [`slivar tsv`](https://github.com/brentp/slivar/wiki/tsv:-creating-a-spreadsheet-from-a-filtered-VCF) to produce a TSV of structural variants that impact genes in affected samples.  This TSV has many of the relevant fields from the VCF, as well as:

- Clinvar annotations for the gene
- gnomAD [loss-of-function tolerance metrics](https://gnomad.broadinstitute.org/downloads#v2-lof-curation-results)
- Phrank scores for the gene

[^1]: High quality is defined as:
  GQ >= 20 (GQ >= 10 for males on chrX)
  DP >= 6
  0.2 <= hetalt AB <= 0.8
  homref AB < 0.02
  homalt AB > 0.98

[^2]: For more description of considered impacts, see [`slivar` documentation](https://github.com/brentp/slivar/wiki/compound-heterozygotes).  We alter the default "skip" list to:
  non_coding_transcript
  intron
  non_coding
  upstream_gene
  downstream_gene
  non_coding_transcript_exon
  NMD_transcript
  5_prime_UTR
  3_prime_UTR

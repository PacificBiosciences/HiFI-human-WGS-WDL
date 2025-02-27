# singleton.wdl inputs and outputs

- [singleton.wdl inputs and outputs](#singletonwdl-inputs-and-outputs)
  - [DAG (simplified)](#dag-simplified)
  - [Inputs](#inputs)
  - [Outputs](#outputs)
    - [Alignments, Coverage, and QC](#alignments-coverage-and-qc)
    - [Small Variants (\<50 bp)](#small-variants-50-bp)
    - [Structural Variants (≥50 bp)](#structural-variants-50-bp)
    - [Copy Number Variants (≥100 kb)](#copy-number-variants-100-kb)
    - [Tandem Repeat Genotyping](#tandem-repeat-genotyping)
    - [Variant Phasing](#variant-phasing)
    - [Variant Calling in Dark Regions](#variant-calling-in-dark-regions)
    - [5mCpG Methylation Calling](#5mcpg-methylation-calling)
    - [PGx Typing](#pgx-typing)
    - [Tertiary Analysis](#tertiary-analysis)

## DAG (simplified)

```mermaid
---
title: singleton.wdl
---
flowchart TD
  subgraph "`**Upstream of Phasing**`"
    subgraph "per-movie"
      ubam[/"HiFi uBAM"/] --> pbmm2_align["pbmm2 align"]
      pbmm2_align --> sawfish_discover["Sawfish discover"]
    end
    pbmm2_align --> merge_read_stats["merge read statistics"]
    pbmm2_align --> samtools_merge["samtools merge"]
    samtools_merge --> mosdepth["mosdepth"]
    samtools_merge --> paraphase["Paraphase"]
    samtools_merge --> hificnv["HiFiCNV"]
    samtools_merge --> trgt["TRGT"]
    samtools_merge --> trgt_dropouts["TR coverage dropouts"]
    samtools_merge --> deepvariant["DeepVariant"]
    sawfish_discover --> sawfish_call["Sawfish call"]
  end
  subgraph "`**Phasing and Downstream**`"
    deepvariant --> hiphase["HiPhase"]
    trgt --> hiphase
    pbsv_call --> hiphase
    hiphase --> bcftools_roh["bcftools roh"]
    hiphase --> bcftools_stats["bcftools stats\n(small variants)"]
    hiphase --> sv_stats["SV stats"]
    hiphase --> cpg_pileup["5mCpG pileup"]
    hiphase --> starphase["StarPhase"]
    hiphase --> pharmcat["PharmCat"]
    starphase --> pharmcat
  end
  subgraph "`**Tertiary Analysis**`"
    hiphase --> slivar_small_variants["slivar small variants"]
    hiphase --> svpack["svpack filter and annotate"]
    svpack --> slivar_svpack["slivar svpack tsv"]
  end
```

## Inputs

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| String | sample_id | Unique identifier for the sample | Alphanumeric characters, periods, dashes, and underscores are allowed. |
| String? | sex | Sample sex<br/>`["MALE", "FEMALE"]` | Used by HiFiCNV and TRGT for genotyping. Allosome karyotype will default to XX unless sex is specified as `"MALE"`. |
| Array\[File\] | hifi_reads | Array of paths to HiFi reads in unaligned BAM format. |  |
| File | [ref_map_file](./ref_map.md) | TSV containing reference genome file paths; must match backend |  |
| String? | phenotypes | Comma-delimited list of HPO terms. | [Human Phenotype Ontology (HPO) phenotypes](https://hpo.jax.org/app/) associated with the cohort.<br/><br/>If omitted, tertiary analysis will be skipped. |
| File? | [tertiary_map_file](./tertiary_map.md) | TSV containing tertiary analysis file paths and thresholds; must match backend | `AF`/`AC`/`nhomalt` thresholds can be modified, but this will affect performance.<br/><br/>If omitted, tertiary analysis will be skipped. |
| Boolean | gpu | Use GPU when possible<br/><br/>Default: `false` | [GPU support](./gpu.md#gpu-support) |
| String | backend | Backend where the workflow will be executed<br/><br/>`["GCP", "Azure", "AWS-AGC", "AWS-HealthOmics", "HPC"]` |  |
| String? | zones | Zones where compute will take place; required if backend is set to 'AWS' or 'GCP'. | [Determining available zones in GCP](./backends/gcp.md#determining-available-zones) |
| String? | gpuType | GPU type to use; required if gpu is set to `true` for cloud backends; must match backend  | [Available GPU types](./gpu.md#gpu-types) |
| String? | container_registry | Container registry where workflow images are hosted.<br/><br/>Default: `"quay.io/pacbio"` | If omitted, [PacBio's public Quay.io registry](https://quay.io/organization/pacbio) will be used.<br/><br/>Custom container_registry must be set if backend is set to 'AWS-HealthOmics'. |
| Boolean | preemptible | Where possible, run tasks preemptibly<br/><br/>`[true, false]`<br/><br/>Default: `true` | If set to `true`, run tasks preemptibly where possible. If set to `false`, on-demand VMs will be used for every task. Ignored if backend is set to HPC. |

## Outputs

### Alignments, Coverage, and QC

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| String | workflow_name | Workflow name |  |
| String | workflow_version | Workflow version |  |
| File | stats_file | Table of summary statistics |  |
| File | bam_stats | BAM stats | Per-read length and read-quality |
| File | read_length_plot | Read length plot |  |
| File? | read_quality_plot | Read quality plot |  |
| File | merged_haplotagged_bam | Merged, haplotagged alignments | Includes unmapped reads |
| File | merged_haplotagged_bam_index |  |  |
| File | mosdepth_summary | Summary of aligned read depth. |  |
| File | mosdepth_region_bed | Median aligned read depth by 500bp windows. |  |
| File | mosdepth_region_bed_index |  |  |
| File | mosdepth_depth_distribution_plot |  |  |
| File | mapq_distribution_plot | Distribution of mapping quality per alignment | |
| File | mg_distribution_plot | Distribution of gap-compressed identity score per alignment | |
| String | stat_num_reads | Number of reads |  |
| String | stat_read_length_mean | Mean read length |  |
| String | stat_read_length_median | Median read length |  |
| String | stat_read_quality_mean | Mean read quality |  |
| String | stat_read_quality_median | Median read quality |  |
| String | stat_mapped_read_count | Count of reads mapped to reference |  |
| String | stat_mapped_percent | Percent of reads mapped to reference |  |
| String | inferred_sex | Inferred sex | Sex is inferred based on relative depth of chrY alignments. |
| String | stat_mean_depth | Mean depth | |

### Small Variants (<50 bp)

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | phased_small_variant_vcf | Phased small variant VCF |  |
| File | phased_small_variant_vcf_index |  |  |
| File | small_variant_gvcf | Small variant GVCF | Can be used for joint-calling. |
| File | small_variant_gvcf_index |  |  |
| File | small_variant_stats | Small variant stats | Generated by `bcftools stats`. |
| String | stat_small_variant_SNV_count | SNV count | (PASS variants) |
| String | stat_small_variant_INDEL_count | INDEL count | (PASS variants) |
| String | stat_small_variant_TSTV_ratio | Ts/Tv ratio | (PASS variants) |
| String | stat_small_variant_HETHOM_ratio | Het/Hom ratio | (PASS variants) |
| File | snv_distribution_plot | Distribution of SNVs by REF, ALT | |
| File | indel_distribution_plot | Distribution of indels by size | |

### Structural Variants (≥50 bp)

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | phased_sv_vcf | Phased structural variant VCF |  |
| File | phased_sv_vcf_index | Index for phased structural variant VCF |  |
| String | stat_sv_DUP_count | Structural variant DUP count | (PASS variants) |
| String | stat_sv_DEL_count | Structural variant DEL count | (PASS variants) |
| String | stat_sv_INS_count | Structural variant INS count | (PASS variants) |
| String | stat_sv_INV_count | Structural variant INV count | (PASS variants) |
| String | stat_sv_INVBND_count | Structural variant INVBND count | (PASS variants) |
| String | stat_sv_BND_count | Structural variant BND count | (PASS variants) |
| File | bcftools_roh_out | ROH calling |  `bcftools roh` |
| File | bcftools_roh_bed | Generated from above, without filtering |  |

### Copy Number Variants (≥100 kb)

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | cnv_vcf | CNV VCF |  |
| File | cnv_vcf_index | Index for CNV VCF |  |
| File | cnv_copynum_bedgraph | CNV copy number BEDGraph |  |
| File | cnv_depth_bw | CNV depth BigWig |  |
| File | cnv_maf_bw | CNV MAF BigWig |  |
| String | stat_cnv_DUP_count | Count of DUP events | (for PASS variants) |
| String | stat_cnv_DEL_count | Count of DEL events | (PASS variants) |
| String | stat_cnv_DUP_sum | Sum of DUP bp | (PASS variants) |
| String | stat_cnv_DEL_sum | Sum of DEL bp | (PASS variants) |

### Tandem Repeat Genotyping

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | phased_trgt_vcf | Phased TRGT VCF |  |
| File | phased_trgt_vcf_index |  |  |
| File | trgt_spanning_reads | TRGT spanning reads |  |
| File | trgt_spanning_reads_index |  |  |
| File | trgt_coverage_dropouts | TRGT coverage dropouts |  |
| String | stat_trgt_genotyped_count | Count of genotyped sites |  |
| String | stat_trgt_uncalled_count | Count of ungenotyped sites |  |

### Variant Phasing

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | phase_stats | Phasing stats |  |
| File | phase_blocks | Phase blocks |  |
| File | phase_haplotags | Per-read haplotag assignment |  |
| String | stat_phased_basepairs | Count of bp within phase blocks |  |
| String | stat_phase_block_ng50 | Phase block NG50 |  |

### Variant Calling in Dark Regions

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | paraphase_output_json | Paraphase output JSON |  |
| File | paraphase_realigned_bam | Paraphase realigned BAM |  |
| File | paraphase_realigned_bam_index |  |  |
| File? | paraphase_vcfs | Paraphase VCFs | Compressed as `.tar.gz` |

### 5mCpG Methylation Calling

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File? | cpg_hap1_bed | CpG hap1 BED |  |
| File? | cpg_hap1_bed_index |  |  |
| File? | cpg_hap2_bed | CpG hap2 BED |  |
| File? | cpg_hap2_bed_index |  |  |
| File? | cpg_combined_bed | CpG combined BED |  |
| File? | cpg_combined_bed_index |  |  |
| File? | cpg_hap1_bw | CpG hap1 BigWig |  |
| File? | cpg_hap2_bw | CpG hap2 BigWig |  |
| File? | cpg_combined_bw | CpG combined BigWig |  |
| String | stat_cpg_hap1_count | Hap1 CpG count |  |
| String | stat_cpg_hap2_count | Hap2 CpG count |  |
| String | stat_cpg_combined_count | Combined CpG count |  |

### PGx Typing

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File | pbstarphase_json | PBstarPhase JSON | Haplotype calls for PGx loci |
| File? | pharmcat_match_json | PharmCAT match JSON |  |
| File? | pharmcat_phenotype_json | PharmCAT phenotype JSON |  |
| File? | pharmcat_report_html | PharmCAT report HTML |  |
| File? | pharmcat_report_json | PharmCAT report JSON |  |

### Tertiary Analysis

| Type | Name | Description | Notes |
| ---- | ---- | ----------- | ----- |
| File? | pedigree | Pedigree file in PLINK PED [format](https://zzz.bwh.harvard.edu/plink/data.shtml#ped) |  |
| File? | tertiary_small_variant_filtered_vcf | Filtered, annotated small variant VCF |  |
| File? | tertiary_small_variant_filtered_vcf_index |  |  |
| File? | tertiary_small_variant_filtered_tsv | Filtered, annotated small variant calls |  |
| File? | tertiary_small_variant_compound_het_vcf | Filtered, annotated compound heterozygous small variant VCF |  |
| File? | tertiary_small_variant_compound_het_vcf_index |  |  |
| File? | tertiary_small_variant_compound_het_tsv | Filtered, annotated compound heterozygous small variant calls |  |
| File? | tertiary_sv_filtered_vcf | Filtered, annotated structural variant VCF |  |
| File? | tertiary_sv_filtered_vcf_index |  |  |
| File? | tertiary_sv_filtered_tsv | Filtered, annotated structural variant TSV |  |

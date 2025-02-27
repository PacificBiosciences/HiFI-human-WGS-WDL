version 1.0

import "humanwgs_structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "upstream/upstream.wdl" as Upstream
import "downstream/downstream.wdl" as Downstream
import "tertiary/tertiary.wdl" as TertiaryAnalysis
import "wdl-common/wdl/tasks/utilities.wdl" as Utilities


workflow humanwgs_singleton {
  meta {
    description: "PacBio HiFi human whole genome sequencing pipeline for individual samples."
  }

  parameter_meta {
    sample_id: {
      name: "Unique identifier for the sample"
    }
    sex: {
      name: "Sample sex",
      choices: ["MALE", "FEMALE"]
    }
    hifi_reads: {
      name: "Array of paths to HiFi reads in unaligned BAM format."
    }
    ref_map_file: {
      name: "TSV containing reference genome file paths; must match backend"
    }
    pharmcat_min_coverage: {
      name: "Minimum coverage for PharmCAT"
    }
    phenotypes: {
      name: "Comma-delimited list of HPO codes for phenotypes"
    }
    tertiary_map_file: {
      name: "TSV containing tertiary analysis file paths and thresholds; must match backend"
    }
    gpu: {
      name: "Use GPU when possible"
    }
    backend: {
      name: "Backend where the workflow will be executed",
      choices: ["GCP", "Azure", "AWS-HealthOmics", "HPC"]
    }
    zones: {
      name: "Zones where compute will take place; required if backend is set to 'GCP'"
    }
    gpuType: {
      name: "GPU type to use; required if gpu is set to `true` for cloud backends; must match backend"
    }
    container_registry: {
      name: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used. Must be set if backend is set to 'AWS-HealthOmics'",
      default: "quay.io/pacbio"
    }
    preemptible: {
      name: "Where possible, run tasks preemptibly"
    }
    debug_version: {
      name: "Debug version for testing purposes"
    }
  }

  input {
    String sample_id
    String? sex
    Array[File] hifi_reads

    File ref_map_file

    Int pharmcat_min_coverage = 10

    String phenotypes = "HP:0000001"
    File? tertiary_map_file

    Boolean gpu = false

    # Backend configuration
    String backend
    String? zones
    String? gpuType
    String? container_registry

    Boolean preemptible = true

    String? debug_version
  }

  call BackendConfiguration.backend_configuration {
    input:
      backend            = backend,
      zones              = zones,
      gpuType            = gpuType,
      container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

  call Upstream.upstream {
    input:
      sample_id                    = sample_id,
      sex                          = sex,
      hifi_reads                   = hifi_reads,
      ref_map_file                 = ref_map_file,
      single_sample                = true,
      gpu                          = gpu,
      default_runtime_attributes   = default_runtime_attributes
  }

  call Downstream.downstream {
    input:
      sample_id                  = sample_id,
      small_variant_vcf          = upstream.small_variant_vcf,
      small_variant_vcf_index    = upstream.small_variant_vcf_index,
      sv_vcf                     = select_first([upstream.sv_vcf]),
      sv_vcf_index               = select_first([upstream.sv_vcf_index]),
      trgt_vcf                   = upstream.trgt_vcf,
      trgt_vcf_index             = upstream.trgt_vcf_index,
      aligned_bam                = upstream.out_bam,
      aligned_bam_index          = upstream.out_bam_index,
      pharmcat_min_coverage      = pharmcat_min_coverage,
      ref_map_file               = ref_map_file,
      default_runtime_attributes = default_runtime_attributes
  }

  Map[String, Array[String]] stats = {
    'sample_id': [sample_id],
    'num_reads': [upstream.stat_num_reads],
    'read_length_mean': [upstream.stat_read_length_mean],
    'read_length_median': [upstream.stat_read_length_median],
    'read_quality_mean': [upstream.stat_read_quality_mean],
    'read_quality_median': [upstream.stat_read_quality_median],
    'mapped_read_count': [downstream.stat_mapped_read_count],
    'mapped_percent': [downstream.stat_mapped_percent],
    'mean_depth': [upstream.stat_mean_depth],
    'inferred_sex': [upstream.inferred_sex],
    'stat_phased_basepairs': [downstream.stat_phased_basepairs],
    'phase_block_ng50': [downstream.stat_phase_block_ng50],
    'cpg_combined_count': [downstream.stat_combined_cpg_count],
    'cpg_hap1_count': [downstream.stat_hap1_cpg_count],
    'cpg_hap2_count': [downstream.stat_hap2_cpg_count],
    'SNV_count': [downstream.stat_SNV_count],
    'TSTV_ratio': [downstream.stat_TSTV_ratio],
    'HETHOM_ratio': [downstream.stat_HETHOM_ratio],
    'INDEL_count': [downstream.stat_INDEL_count],
    'sv_DUP_count': [downstream.stat_sv_DUP_count],
    'sv_DEL_count': [downstream.stat_sv_DEL_count],
    'sv_INS_count': [downstream.stat_sv_INS_count],
    'sv_INV_count': [downstream.stat_sv_INV_count],
    'sv_INVBND_count': [downstream.stat_sv_INVBND_count],
    'sv_BND_count': [downstream.stat_sv_BND_count],
    'cnv_DUP_count': [upstream.stat_cnv_DUP_count],
    'cnv_DEL_count': [upstream.stat_cnv_DEL_count],
    'cnv_DUP_sum': [upstream.stat_cnv_DUP_sum],
    'cnv_DEL_sum': [upstream.stat_cnv_DEL_sum],
    'trgt_genotyped_count': [upstream.stat_trgt_genotyped_count],
    'trgt_uncalled_count': [upstream.stat_trgt_uncalled_count]
  }

  call Utilities.consolidate_stats {
    input:
      id                 = sample_id,
      stats              = stats,
      runtime_attributes = default_runtime_attributes
  }

  Map[String, String] pedigree_sex = {
    "MALE": "1",
    "FEMALE": "2",
    "": "."
  }

  # write sample metadata similar to pedigree format
  # family_id, sample_id, father_id, mother_id, sex, affected
  Array[String] sample_metadata = [
    sample_id, sample_id,
    ".", ".",
    pedigree_sex[upstream.inferred_sex], "2"
  ]

  if (defined(tertiary_map_file)) {
    call TertiaryAnalysis.tertiary_analysis {
      input:
        sample_metadata            = [sample_metadata],
        phenotypes                 = phenotypes,
        is_trio_kid                = [false],
        is_duo_kid                 = [false],
        small_variant_vcf          = downstream.phased_small_variant_vcf,
        small_variant_vcf_index    = downstream.phased_small_variant_vcf_index,
        sv_vcf                     = downstream.phased_sv_vcf,
        sv_vcf_index               = downstream.phased_sv_vcf_index,
        ref_map_file               = ref_map_file,
        tertiary_map_file          = select_first([tertiary_map_file]),
        default_runtime_attributes = default_runtime_attributes
    }
  }

  output {
    # consolidated stats
    File stats_file = consolidate_stats.output_tsv

    # bam stats
    File   bam_stats                = upstream.read_length_and_quality
    File   read_length_plot         = upstream.read_length_plot
    File?  read_quality_plot        = upstream.read_quality_plot
    String stat_num_reads           = upstream.stat_num_reads
    String stat_read_length_mean    = upstream.stat_read_length_mean
    String stat_read_length_median  = upstream.stat_read_length_median
    String stat_read_quality_mean   = upstream.stat_read_quality_mean
    String stat_read_quality_median = upstream.stat_read_quality_median

    # merged, haplotagged alignments
    File   merged_haplotagged_bam       = downstream.merged_haplotagged_bam
    File   merged_haplotagged_bam_index = downstream.merged_haplotagged_bam_index
    String stat_mapped_read_count       = downstream.stat_mapped_read_count
    String stat_mapped_percent          = downstream.stat_mapped_percent
    File   mapq_distribution_plot       = downstream.mapq_distribution_plot
    File   mg_distribution_plot         = downstream.mg_distribution_plot

    # mosdepth outputs
    File   mosdepth_summary                 = upstream.mosdepth_summary
    File   mosdepth_region_bed              = upstream.mosdepth_region_bed
    File   mosdepth_region_bed_index        = upstream.mosdepth_region_bed_index
    File   mosdepth_depth_distribution_plot = upstream.mosdepth_depth_distribution_plot
    String stat_mean_depth                  = upstream.stat_mean_depth
    String inferred_sex                     = upstream.inferred_sex

    # phasing stats
    File   phase_stats           = downstream.phase_stats
    File   phase_blocks          = downstream.phase_blocks
    File   phase_haplotags       = downstream.phase_haplotags
    String stat_phased_basepairs = downstream.stat_phased_basepairs
    String stat_phase_block_ng50 = downstream.stat_phase_block_ng50

    # cpg_pileup outputs
    File?  cpg_combined_bed        = downstream.cpg_combined_bed
    File?  cpg_combined_bed_index  = downstream.cpg_combined_bed_index
    File?  cpg_hap1_bed            = downstream.cpg_hap1_bed
    File?  cpg_hap1_bed_index      = downstream.cpg_hap1_bed_index
    File?  cpg_hap2_bed            = downstream.cpg_hap2_bed
    File?  cpg_hap2_bed_index      = downstream.cpg_hap2_bed_index
    File?  cpg_combined_bw         = downstream.cpg_combined_bw
    File?  cpg_hap1_bw             = downstream.cpg_hap1_bw
    File?  cpg_hap2_bw             = downstream.cpg_hap2_bw
    String stat_cpg_hap1_count     = downstream.stat_hap1_cpg_count
    String stat_cpg_hap2_count     = downstream.stat_hap2_cpg_count
    String stat_cpg_combined_count = downstream.stat_combined_cpg_count

    # sv outputs
    File phased_sv_vcf       = downstream.phased_sv_vcf
    File phased_sv_vcf_index = downstream.phased_sv_vcf_index

    # sv stats
    String stat_sv_DUP_count    = downstream.stat_sv_DUP_count
    String stat_sv_DEL_count    = downstream.stat_sv_DEL_count
    String stat_sv_INS_count    = downstream.stat_sv_INS_count
    String stat_sv_INV_count    = downstream.stat_sv_INV_count
    String stat_sv_INVBND_count = downstream.stat_sv_INVBND_count
    String stat_sv_BND_count    = downstream.stat_sv_BND_count

    # small variant outputs
    File phased_small_variant_vcf       = downstream.phased_small_variant_vcf
    File phased_small_variant_vcf_index = downstream.phased_small_variant_vcf_index
    File small_variant_gvcf             = upstream.small_variant_gvcf
    File small_variant_gvcf_index       = upstream.small_variant_gvcf_index

    # small variant stats
    File   small_variant_stats             = downstream.small_variant_stats
    File   bcftools_roh_out                = downstream.bcftools_roh_out
    File   bcftools_roh_bed                = downstream.bcftools_roh_bed
    String stat_small_variant_SNV_count    = downstream.stat_SNV_count
    String stat_small_variant_INDEL_count  = downstream.stat_INDEL_count
    String stat_small_variant_TSTV_ratio   = downstream.stat_TSTV_ratio
    String stat_small_variant_HETHOM_ratio = downstream.stat_HETHOM_ratio
    File   snv_distribution_plot           = downstream.snv_distribution_plot
    File   indel_distribution_plot         = downstream.indel_distribution_plot

    # trgt outputs
    File   phased_trgt_vcf           = downstream.phased_trgt_vcf
    File   phased_trgt_vcf_index     = downstream.phased_trgt_vcf_index
    File   trgt_spanning_reads       = upstream.trgt_spanning_reads
    File   trgt_spanning_reads_index = upstream.trgt_spanning_reads_index
    File   trgt_coverage_dropouts    = upstream.trgt_coverage_dropouts
    String stat_trgt_genotyped_count = upstream.stat_trgt_genotyped_count
    String stat_trgt_uncalled_count  = upstream.stat_trgt_uncalled_count

    # paraphase outputs
    File  paraphase_output_json         = upstream.paraphase_output_json
    File  paraphase_realigned_bam       = upstream.paraphase_realigned_bam
    File  paraphase_realigned_bam_index = upstream.paraphase_realigned_bam_index
    File? paraphase_vcfs                = upstream.paraphase_vcfs

    # per sample cnv outputs
    File   cnv_vcf              = upstream.cnv_vcf
    File   cnv_vcf_index        = upstream.cnv_vcf_index
    File   cnv_copynum_bedgraph = upstream.cnv_copynum_bedgraph
    File   cnv_depth_bw         = upstream.cnv_depth_bw
    File   cnv_maf_bw           = upstream.cnv_maf_bw
    String stat_cnv_DUP_count   = upstream.stat_cnv_DUP_count
    String stat_cnv_DEL_count   = upstream.stat_cnv_DEL_count
    String stat_cnv_DUP_sum     = upstream.stat_cnv_DUP_sum
    String stat_cnv_DEL_sum     = upstream.stat_cnv_DEL_sum

    # PGx outputs
    File  pbstarphase_json        = downstream.pbstarphase_json
    File? pharmcat_match_json     = downstream.pharmcat_match_json
    File? pharmcat_phenotype_json = downstream.pharmcat_phenotype_json
    File? pharmcat_report_html    = downstream.pharmcat_report_html
    File? pharmcat_report_json    = downstream.pharmcat_report_json

    # tertiary analysis outputs
    File? tertiary_small_variant_filtered_vcf           = tertiary_analysis.small_variant_filtered_vcf
    File? tertiary_small_variant_filtered_vcf_index     = tertiary_analysis.small_variant_filtered_vcf_index
    File? tertiary_small_variant_filtered_tsv           = tertiary_analysis.small_variant_filtered_tsv
    File? tertiary_small_variant_compound_het_vcf       = tertiary_analysis.small_variant_compound_het_vcf
    File? tertiary_small_variant_compound_het_vcf_index = tertiary_analysis.small_variant_compound_het_vcf_index
    File? tertiary_small_variant_compound_het_tsv       = tertiary_analysis.small_variant_compound_het_tsv
    File? tertiary_sv_filtered_vcf                      = tertiary_analysis.sv_filtered_vcf
    File? tertiary_sv_filtered_vcf_index                = tertiary_analysis.sv_filtered_vcf_index
    File? tertiary_sv_filtered_tsv                      = tertiary_analysis.sv_filtered_tsv

    # workflow metadata
    String workflow_name    = "humanwgs_family"
    String workflow_version = "v3.0.0-alpha1" + if defined(debug_version) then "~{"-" + debug_version}" else ""
  }
}

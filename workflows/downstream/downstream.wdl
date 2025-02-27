version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/tasks/hiphase.wdl" as Hiphase
import "../wdl-common/wdl/tasks/bcftools.wdl" as Bcftools
import "../wdl-common/wdl/tasks/cpg_pileup.wdl" as Cpgpileup
import "../wdl-common/wdl/tasks/pbstarphase.wdl" as Pbstarphase
import "../wdl-common/wdl/workflows/pharmcat/pharmcat.wdl" as Pharmcat

workflow downstream {
  meta {
    description: "Phases small variants, SVs, and TRGTs, haplotags alignments, calls HLA and PGx alleles."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    small_variant_vcf: {
      name: "Small variant VCF"
    }
    small_variant_vcf_index: {
      name: "Small variant VCF index"
    }
    sv_vcf: {
      name: "Structural variant VCF"
    }
    sv_vcf_index: {
      name: "Structural variant VCF index"
    }
    trgt_vcf: {
      name: "TRGT VCF"
    }
    trgt_vcf_index: {
      name: "TRGT VCF index"
    }
    aligned_bam: {
      name: "Aligned BAM"
    }
    aligned_bam_index: {
      name: "Aligned BAI"
    }
    pharmcat_min_coverage: {
      name: "Minimum coverage for PharmCAT"
    }
    ref_map_file: {
      name: "Reference map file"
    }
    default_runtime_attributes: {
      name: "Default runtime attributes"
    }
  }

  input {
    String sample_id

    File small_variant_vcf
    File small_variant_vcf_index
    File sv_vcf
    File sv_vcf_index
    File trgt_vcf
    File trgt_vcf_index

    File aligned_bam
    File aligned_bam_index

    Int pharmcat_min_coverage

    File ref_map_file

    RuntimeAttributes default_runtime_attributes
  }

  Map[String, String] ref_map = read_map(ref_map_file)

  Array[File] hiphase_input_vcfs = [small_variant_vcf, sv_vcf, trgt_vcf]
  Array[File] hiphase_input_vcf_indices = [small_variant_vcf_index, sv_vcf_index, trgt_vcf_index]

  scatter (vcf_index in range(length(hiphase_input_vcfs))) {
    # generate an array of phased VCF names that match the input VCFs
    String phased_vcf_name       = basename(hiphase_input_vcfs[vcf_index], ".vcf.gz") + ".phased.vcf.gz"
    String phased_vcf_index_name = basename(hiphase_input_vcf_indices[vcf_index], ".vcf.gz.tbi") + ".phased.vcf.gz.tbi"
  }

  call Hiphase.hiphase {
    input: 
      sample_id              = sample_id,
      vcfs                   = hiphase_input_vcfs,
      vcf_indices            = hiphase_input_vcf_indices,
      phased_vcf_names       = phased_vcf_name,
      phased_vcf_index_names = phased_vcf_index_name,
      aligned_bam            = aligned_bam,
      aligned_bam_index      = aligned_bam_index,
      ref_name               = ref_map["name"],
      ref_fasta              = ref_map["fasta"],          # !FileCoercion
      ref_index              = ref_map["fasta_index"],    # !FileCoercion
      runtime_attributes     = default_runtime_attributes
  }

  # hiphase.phased_vcfs[0] -> phased small variant VCF
  # hiphase.phased_vcfs[1] -> phased SV VCF
  # hiphase.phased_vcfs[2] -> phased TRGT VCF

  call Bcftools.bcftools_stats_roh_small_variants {
    input:
      sample_id          = sample_id,
      vcf                = hiphase.phased_vcfs[0],
      ref_fasta          = ref_map["fasta"],       # !FileCoercion
      ref_name           = ref_map["name"],
      runtime_attributes = default_runtime_attributes
  }

  call Bcftools.sv_stats {
    input:
      vcf                = hiphase.phased_vcfs[1],
      runtime_attributes = default_runtime_attributes
  }

  call Cpgpileup.cpg_pileup {
    input: 
      haplotagged_bam       = hiphase.haplotagged_bam,
      haplotagged_bam_index = hiphase.haplotagged_bam_index,
      out_prefix            = "~{sample_id}.~{ref_map['name']}",
      ref_fasta             = ref_map["fasta"],                  # !FileCoercion
      ref_index             = ref_map["fasta_index"],            # !FileCoercion
      runtime_attributes    = default_runtime_attributes
  }

  call Pbstarphase.pbstarphase_diplotype {
    input:
      sample_id                           = sample_id,
      phased_small_variant_vcf            = hiphase.phased_vcfs[0],
      phased_small_variant_vcf_index      = hiphase.phased_vcf_indices[0],
      phased_structural_variant_vcf       = hiphase.phased_vcfs[1],
      phased_structural_variant_vcf_index = hiphase.phased_vcf_indices[1],
      aligned_bam                         = hiphase.haplotagged_bam,
      aligned_bam_index                   = hiphase.haplotagged_bam_index,
      ref_fasta                           = ref_map["fasta"],              # !FileCoercion
      ref_index                           = ref_map["fasta_index"],        # !FileCoercion
      runtime_attributes                  = default_runtime_attributes
  }

  call Pharmcat.pharmcat {
    input:
      sample_id                  = sample_id,
      haplotagged_bam            = hiphase.haplotagged_bam,
      haplotagged_bam_index      = hiphase.haplotagged_bam_index,
      phased_vcf                 = hiphase.phased_vcfs[0],
      phased_vcf_index           = hiphase.phased_vcf_indices[0],
      input_tsvs                 = [pbstarphase_diplotype.pharmcat_tsv],
      ref_fasta                  = ref_map["fasta"],                        # !FileCoercion
      ref_index                  = ref_map["fasta_index"],                  # !FileCoercion
      pharmcat_positions         = ref_map["pharmcat_positions_vcf"],       # !FileCoercion
      pharmcat_positions_index   = ref_map["pharmcat_positions_vcf_index"], # !FileCoercion
      pharmcat_min_coverage      = pharmcat_min_coverage,
      default_runtime_attributes = default_runtime_attributes
  }

  output {
    # hiphase outputs
    File   merged_haplotagged_bam         = hiphase.haplotagged_bam
    File   merged_haplotagged_bam_index   = hiphase.haplotagged_bam_index
    File   phased_small_variant_vcf       = hiphase.phased_vcfs[0]
    File   phased_small_variant_vcf_index = hiphase.phased_vcf_indices[0]
    File   phased_sv_vcf                  = hiphase.phased_vcfs[1]
    File   phased_sv_vcf_index            = hiphase.phased_vcf_indices[1]
    File   phased_trgt_vcf                = hiphase.phased_vcfs[2]
    File   phased_trgt_vcf_index          = hiphase.phased_vcf_indices[2]
    File   phase_stats                    = hiphase.phase_stats
    File   phase_blocks                   = hiphase.phase_blocks
    File   phase_haplotags                = hiphase.phase_haplotags
    String stat_phased_basepairs          = hiphase.stat_phased_basepairs
    String stat_phase_block_ng50          = hiphase.stat_phase_block_ng50
    String stat_mapped_read_count         = hiphase.stat_mapped_read_count
    String stat_mapped_percent            = hiphase.stat_mapped_percent
    File   mapq_distribution_plot         = hiphase.mapq_distribution_plot
    File   mg_distribution_plot           = hiphase.mg_distribution_plot

    # small variant stats
    File   small_variant_stats     = bcftools_stats_roh_small_variants.stats
    File   bcftools_roh_out        = bcftools_stats_roh_small_variants.roh_out
    File   bcftools_roh_bed        = bcftools_stats_roh_small_variants.roh_bed
    String stat_SNV_count          = bcftools_stats_roh_small_variants.stat_SNV_count
    String stat_INDEL_count        = bcftools_stats_roh_small_variants.stat_INDEL_count
    String stat_TSTV_ratio         = bcftools_stats_roh_small_variants.stat_TSTV_ratio
    String stat_HETHOM_ratio       = bcftools_stats_roh_small_variants.stat_HETHOM_ratio
    File   snv_distribution_plot   = bcftools_stats_roh_small_variants.snv_distribution_plot
    File   indel_distribution_plot = bcftools_stats_roh_small_variants.indel_distribution_plot

    # sv stats
    String stat_sv_DUP_count    = sv_stats.stat_sv_DUP_count
    String stat_sv_DEL_count    = sv_stats.stat_sv_DEL_count
    String stat_sv_INS_count    = sv_stats.stat_sv_INS_count
    String stat_sv_INV_count    = sv_stats.stat_sv_INV_count
    String stat_sv_INVBND_count = sv_stats.stat_sv_INVBND_count
    String stat_sv_BND_count    = sv_stats.stat_sv_BND_count

    # cpg_pileup outputs
    File?  cpg_combined_bed        = cpg_pileup.combined_bed
    File?  cpg_combined_bed_index  = cpg_pileup.combined_bed_index
    File?  cpg_hap1_bed            = cpg_pileup.hap1_bed
    File?  cpg_hap1_bed_index      = cpg_pileup.hap1_bed_index
    File?  cpg_hap2_bed            = cpg_pileup.hap2_bed
    File?  cpg_hap2_bed_index      = cpg_pileup.hap2_bed_index
    File?  cpg_combined_bw         = cpg_pileup.combined_bw
    File?  cpg_hap1_bw             = cpg_pileup.hap1_bw
    File?  cpg_hap2_bw             = cpg_pileup.hap2_bw
    String stat_hap1_cpg_count     = cpg_pileup.stat_hap1_cpg_count
    String stat_hap2_cpg_count     = cpg_pileup.stat_hap2_cpg_count
    String stat_combined_cpg_count = cpg_pileup.stat_combined_cpg_count

    # pbstarphase outputs
    File pbstarphase_json = pbstarphase_diplotype.out_json

    # pharmcat and pangu outputs
    File? pharmcat_match_json     = pharmcat.pharmcat_match_json
    File? pharmcat_phenotype_json = pharmcat.pharmcat_phenotype_json
    File? pharmcat_report_html    = pharmcat.pharmcat_report_html
    File? pharmcat_report_json    = pharmcat.pharmcat_report_json
  }
}
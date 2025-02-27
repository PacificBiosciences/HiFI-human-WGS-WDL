version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/tasks/glnexus.wdl" as Glnexus
import "../wdl-common/wdl/tasks/sawfish.wdl" as Sawfish
import "../wdl-common/wdl/tasks/bcftools.wdl" as Bcftools

workflow joint {
  meta {
    description: "Tasks for joint-calling variants from a set of samples and splitting the joint calls by sample for parallel phasing."
  }

  parameter_meta {
    family_id: {
      name: "Cohort ID"
    }
    sample_ids: {
      name: "Sample IDs"
    }
    gvcfs: {
      name: "GVCFs"
    }   
    gvcf_indices: {
      name: "GVCF Indices"
    }
    discover_tars: {
      name: "Sawfish discover output tarballs"
    }
    aligned_bams: {
      name: "Aligned BAMs"
    }
    aligned_bam_indices: {
      name: "Aligned BAM Indices"
    }
    ref_map_file: {
      name: "Reference Map File"
    }
    glnexus_mem_gb: {
      name: "GLnexus Memory (GB)"
    }
    default_runtime_attributes: {
      name: "Default Runtime Attribute Struct"
    }
    split_joint_structural_variant_vcfs: {
      name: "Joint-call structural variant VCF, split by sample"
    }
    split_joint_structural_variant_vcf_indices: {
      name: "Joint-call structural variant VCF indices, split by sample"
    }
    split_joint_small_variant_vcfs: {
      name: "Joint-call small variant VCF, split by sample"
    }
    split_joint_small_variant_vcf_indices: {
      name: "Joint-call small variant VCF indices, split by sample"
    }
  }

  input {
    String family_id
    Array[String] sample_ids

    Array[File] gvcfs
    Array[File] gvcf_indices

    Array[File] discover_tars
    Array[File] aligned_bams
    Array[File] aligned_bam_indices

    File ref_map_file

    Int? glnexus_mem_gb

    RuntimeAttributes default_runtime_attributes
  }

  Map[String, String] ref_map = read_map(ref_map_file)

  call Sawfish.sawfish_call {
    input:
      discover_tars       = discover_tars,
      aligned_bams        = aligned_bams,
      aligned_bam_indices = aligned_bam_indices,
      ref_fasta           = ref_map["fasta"],                                      # !FileCoercion
      ref_index           = ref_map["fasta_index"],                                # !FileCoercion
      out_prefix          = "~{family_id}.joint.~{ref_map['name']}.structural_variants",
      runtime_attributes  = default_runtime_attributes
  }

  String sv_vcf_basename = basename(sawfish_call.vcf, ".vcf.gz")

  scatter (sample_id in sample_ids) {
    String split_sv_vcf_name = "~{sample_id}.~{sv_vcf_basename}.vcf.gz"
    String split_sv_vcf_index_name = "~{sample_id}.~{sv_vcf_basename}.vcf.gz.tbi"
  }

  call Bcftools.split_vcf_by_sample as split_sawfish {
    input:
      sample_ids            = sample_ids,
      vcf                   = sawfish_call.vcf,
      vcf_index             = sawfish_call.vcf_index,
      split_vcf_names       = split_sv_vcf_name,
      split_vcf_index_names = split_sv_vcf_index_name,
      runtime_attributes    = default_runtime_attributes
  }

  call Glnexus.glnexus {
    input:
      cohort_id          = family_id + ".joint",
      gvcfs              = gvcfs,
      gvcf_indices       = gvcf_indices,
      ref_name           = ref_map["name"],
      mem_gb             = glnexus_mem_gb,
      runtime_attributes = default_runtime_attributes
  }

  String glnexus_vcf_basename = basename(glnexus.vcf, ".vcf.gz")

  scatter (sample_id in sample_ids) {
    String split_glnexus_vcf_name = "~{sample_id}.~{glnexus_vcf_basename}.vcf.gz"
    String split_glnexus_vcf_index_name = "~{sample_id}.~{glnexus_vcf_basename}.vcf.gz.tbi"
  }

  call Bcftools.split_vcf_by_sample as split_glnexus {
    input:
      sample_ids            = sample_ids,
      vcf                   = glnexus.vcf,
      vcf_index             = glnexus.vcf_index,
      split_vcf_names       = split_glnexus_vcf_name,
      split_vcf_index_names = split_glnexus_vcf_index_name,
      runtime_attributes    = default_runtime_attributes
  }

  output {
    Array[File] split_joint_structural_variant_vcfs        = split_sawfish.split_vcfs
    Array[File] split_joint_structural_variant_vcf_indices = split_sawfish.split_vcf_indices
    Array[File] split_joint_small_variant_vcfs             = split_glnexus.split_vcfs
    Array[File] split_joint_small_variant_vcf_indices      = split_glnexus.split_vcf_indices
  }
}

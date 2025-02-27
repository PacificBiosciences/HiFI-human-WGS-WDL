version 1.0

import "../wdl-common/wdl/structs.wdl"
import "../wdl-common/wdl/tasks/pbmm2.wdl" as Pbmm2
import "../wdl-common/wdl/tasks/merge_bam_stats.wdl" as MergeBamStats
import "../wdl-common/wdl/tasks/sawfish.wdl" as Sawfish
import "../wdl-common/wdl/workflows/deepvariant/deepvariant.wdl" as DeepVariant
import "../wdl-common/wdl/tasks/samtools.wdl" as Samtools
import "../wdl-common/wdl/tasks/mosdepth.wdl" as Mosdepth
import "../wdl-common/wdl/tasks/trgt.wdl" as Trgt
import "../wdl-common/wdl/tasks/paraphase.wdl" as Paraphase
import "../wdl-common/wdl/tasks/hificnv.wdl" as Hificnv

workflow upstream {
  meta {
    description: "Given a set of HiFi reads for a human sample, run steps upstream of phasing."
  }

  parameter_meta {
    sample_id: {
      name: "Sample ID"
    }
    sex: {
      name: "Sample sex",
      choices: ["MALE", "FEMALE"]
    }
    hifi_reads: {
      name: "HiFi reads (BAMs)"
    }
    ref_map_file: {
      name: "TSV containing reference genome information"
    }
    single_sample: {
      name: "Single sample workflow"
    }
    gpu: {
      name: "Use GPU for DeepVariant"
    }
    default_runtime_attributes: {
      name: "Runtime attribute structure"
    }
  }

  input {
    String sample_id
    String? sex
    Array[File] hifi_reads

    File ref_map_file

    Boolean single_sample = false

    Boolean gpu

    RuntimeAttributes default_runtime_attributes
  }

  Map[String, String] ref_map = read_map(ref_map_file)

  scatter (hifi_read_bam in hifi_reads) {
    call Pbmm2.pbmm2_align_wgs as pbmm2_align {
      input:
        sample_id          = sample_id,
        bam                = hifi_read_bam,
        ref_fasta          = ref_map["fasta"],       # !FileCoercion
        ref_index          = ref_map["fasta_index"], # !FileCoercion
        ref_name           = ref_map["name"],
        runtime_attributes = default_runtime_attributes
    }
  }

  call MergeBamStats.merge_bam_stats {
    input:
      sample_id            = sample_id,
      bam_stats            = pbmm2_align.bam_stats,
      runtime_attributes   = default_runtime_attributes
  }

  # merge aligned bams if there are multiple
  if (length(pbmm2_align.aligned_bam) > 1) {
    call Samtools.samtools_merge {
      input:
        bams               = pbmm2_align.aligned_bam,
        out_prefix         = "~{sample_id}.~{ref_map['name']}",
        runtime_attributes = default_runtime_attributes
    }
  }

  # select the merged bam if it exists, otherwise select the first (only) aligned bam
  File aligned_bam_data  = select_first([samtools_merge.merged_bam, pbmm2_align.aligned_bam[0]])
  File aligned_bam_index = select_first([samtools_merge.merged_bam_index, pbmm2_align.aligned_bam_index[0]])

  call Mosdepth.mosdepth {
    input:
      sample_id          = sample_id,
      ref_name           = ref_map["name"],
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      infer_sex          = true,
      runtime_attributes = default_runtime_attributes
  }

  call DeepVariant.deepvariant {
    input:
      sample_id                    = sample_id,
      aligned_bams                 = [aligned_bam_data],
      aligned_bam_indices          = [aligned_bam_index],
      ref_fasta                    = ref_map["fasta"],             # !FileCoercion
      ref_index                    = ref_map["fasta_index"],       # !FileCoercion
      ref_name                     = ref_map["name"],
      gpu                          = gpu,
      default_runtime_attributes   = default_runtime_attributes
  }

  call Sawfish.sawfish_discover {
    input:
      sex                 = select_first([sex, mosdepth.inferred_sex]),
      aligned_bam         = aligned_bam_data,
      aligned_bam_index   = aligned_bam_index,
      ref_fasta           = ref_map["fasta"],                           # !FileCoercion
      ref_index           = ref_map["fasta_index"],                     # !FileCoercion
      out_prefix          = "~{sample_id}.~{ref_map['name']}",
      expected_male_bed   = ref_map["hificnv_expected_bed_male"],       # !FileCoercion
      expected_female_bed = ref_map["hificnv_expected_bed_female"],     # !FileCoercion
      runtime_attributes  = default_runtime_attributes
  }

  call Trgt.trgt {
    input:
      sample_id          = sample_id,
      sex                = select_first([sex, mosdepth.inferred_sex]),
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      ref_fasta          = ref_map["fasta"],                           # !FileCoercion
      ref_index          = ref_map["fasta_index"],                     # !FileCoercion
      trgt_bed           = ref_map["trgt_tandem_repeat_bed"],          # !FileCoercion
      out_prefix         = "~{sample_id}.~{ref_map['name']}",
      runtime_attributes = default_runtime_attributes
  }

  call Trgt.coverage_dropouts {
    input: 
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      trgt_bed           = ref_map["trgt_tandem_repeat_bed"], # !FileCoercion
      out_prefix         = "~{sample_id}.~{ref_map['name']}",
      runtime_attributes = default_runtime_attributes
  }

  call Paraphase.paraphase {
    input:
      aligned_bam        = aligned_bam_data,
      aligned_bam_index  = aligned_bam_index,
      ref_fasta          = ref_map["fasta"],         # !FileCoercion
      ref_index          = ref_map["fasta_index"],   # !FileCoercion
      sample_id          = sample_id,
      runtime_attributes = default_runtime_attributes
  }

  call Hificnv.hificnv {
    input:
      sample_id           = sample_id,
      sex                 = select_first([sex, mosdepth.inferred_sex]),
      aligned_bam         = aligned_bam_data,
      aligned_bam_index   = aligned_bam_index,
      vcf                 = deepvariant.vcf,
      vcf_index           = deepvariant.vcf_index,
      ref_fasta           = ref_map["fasta"],                           # !FileCoercion
      ref_index           = ref_map["fasta_index"],                     # !FileCoercion
      ref_name            = ref_map["name"],
      exclude_bed         = ref_map["hificnv_exclude_bed"],             # !FileCoercion
      exclude_bed_index   = ref_map["hificnv_exclude_bed_index"],       # !FileCoercion
      expected_male_bed   = ref_map["hificnv_expected_bed_male"],       # !FileCoercion
      expected_female_bed = ref_map["hificnv_expected_bed_female"],     # !FileCoercion
      runtime_attributes  = default_runtime_attributes
  }

  if (single_sample) {
    call Sawfish.sawfish_call {
      input: 
        discover_tars       = [sawfish_discover.discover_tar],
        aligned_bams        = [aligned_bam_data],
        aligned_bam_indices = [aligned_bam_index],
        ref_fasta           = ref_map["fasta"],                                      # !FileCoercion
        ref_index           = ref_map["fasta_index"],                                # !FileCoercion
        out_prefix          = "~{sample_id}.~{ref_map['name']}.structural_variants",
        runtime_attributes  = default_runtime_attributes
    }
  }

  output {
    # bam stats
    File   read_length_and_quality  = merge_bam_stats.read_length_and_quality
    File   read_length_plot         = merge_bam_stats.read_length_plot
    File?  read_quality_plot        = merge_bam_stats.read_quality_plot
    String stat_num_reads           = merge_bam_stats.stat_num_reads
    String stat_read_length_mean    = merge_bam_stats.stat_read_length_mean
    String stat_read_length_median  = merge_bam_stats.stat_read_length_median
    String stat_read_quality_mean   = merge_bam_stats.stat_read_quality_mean
    String stat_read_quality_median = merge_bam_stats.stat_read_quality_median

    # alignments
    File out_bam       = aligned_bam_data
    File out_bam_index = aligned_bam_index

    # mosdepth outputs
    File   mosdepth_summary                 = mosdepth.summary
    File   mosdepth_region_bed              = mosdepth.region_bed
    File   mosdepth_region_bed_index        = mosdepth.region_bed_index
    File   mosdepth_depth_distribution_plot = mosdepth.depth_distribution_plot
    String inferred_sex                     = mosdepth.inferred_sex
    String stat_mean_depth                  = mosdepth.stat_mean_depth

    # per sample sv signatures
    File discover_tar = sawfish_discover.discover_tar

    # sawfish outputs for single sample
    File? sv_vcf       = sawfish_call.vcf
    File? sv_vcf_index = sawfish_call.vcf_index

    # small variant outputs
    File small_variant_vcf        = deepvariant.vcf
    File small_variant_vcf_index  = deepvariant.vcf_index
    File small_variant_gvcf       = deepvariant.gvcf
    File small_variant_gvcf_index = deepvariant.gvcf_index

    # trgt outputs
    File   trgt_vcf                  = trgt.vcf
    File   trgt_vcf_index            = trgt.vcf_index
    File   trgt_spanning_reads       = trgt.bam
    File   trgt_spanning_reads_index = trgt.bam_index
    File   trgt_coverage_dropouts    = coverage_dropouts.dropouts
    String stat_trgt_genotyped_count = trgt.stat_genotyped_count
    String stat_trgt_uncalled_count  = trgt.stat_uncalled_count

    # paraphase outputs
    File  paraphase_output_json         = paraphase.out_json
    File  paraphase_realigned_bam       = paraphase.bam
    File  paraphase_realigned_bam_index = paraphase.bam_index
    File? paraphase_vcfs                = paraphase.vcfs_tar

    # per sample hificnv outputs
    File   cnv_vcf              = hificnv.cnv_vcf
    File   cnv_vcf_index        = hificnv.cnv_vcf_index
    File   cnv_copynum_bedgraph = hificnv.copynum_bedgraph
    File   cnv_depth_bw         = hificnv.depth_bw
    File   cnv_maf_bw           = hificnv.maf_bw
    String stat_cnv_DUP_count   = hificnv.stat_DUP_count
    String stat_cnv_DEL_count   = hificnv.stat_DEL_count
    String stat_cnv_DUP_sum     = hificnv.stat_DUP_sum
    String stat_cnv_DEL_sum     = hificnv.stat_DEL_sum
  }
}

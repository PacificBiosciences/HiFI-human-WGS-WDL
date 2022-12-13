version 1.0

import "../common/structs.wdl"
import "../common/tasks/mosdepth.wdl" as Mosdepth

workflow smrtcell_analysis {
	input {
		Sample sample

		ReferenceData reference

		String deepvariant_version
		DeepVariantModel? deepvariant_model

		String container_registry
	}

	Int deepvariant_threads = 64

	scatter (movie_bam in sample.movie_bams) {
		call pbmm2_align {
			input:
				sample_id = sample.sample_id,
				bam = movie_bam.data,
				reference = reference.fasta.data,
				reference_index = reference.fasta.data_index,
				reference_name = reference.name,
				container_registry = container_registry
		}

		call Mosdepth.mosdepth {
			input:
				aligned_bam = pbmm2_align.aligned_bam,
				aligned_bam_index = pbmm2_align.aligned_bam_index,
				container_registry = container_registry
		}

		call pbsv_discover {
			input:
				aligned_bam = pbmm2_align.aligned_bam,
				aligned_bam_index = pbmm2_align.aligned_bam_index,
				reference_tandem_repeat_bed = reference.tandem_repeat_bed,
				container_registry = container_registry
		}

		IndexData aligned_bam = {
			"data": pbmm2_align.aligned_bam,
			"data_index": pbmm2_align.aligned_bam_index
		}
	}

	call deepvariant_make_examples {
		input:
			sample_id = sample.sample_id,
			aligned_bams = pbmm2_align.aligned_bam,
			aligned_bam_indices = pbmm2_align.aligned_bam_index,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			deepvariant_threads = deepvariant_threads,
			deepvariant_version = deepvariant_version
	}

	call deepvariant_call_variants {
		input:
			sample_id = sample.sample_id,
			reference_name = reference.name,
			example_tfrecords = deepvariant_make_examples.example_tfrecords,
			deepvariant_model = deepvariant_model,
			deepvariant_threads = deepvariant_threads,
			deepvariant_version = deepvariant_version
	}

	call deepvariant_postprocess_variants {
		input:
			sample_id = sample.sample_id,
			tfrecord = deepvariant_call_variants.tfrecord,
			nonvariant_site_tfrecords = deepvariant_make_examples.nonvariant_site_tfrecords,
			reference = reference.fasta.data,
			reference_index = reference.fasta.data_index,
			reference_name = reference.name,
			deepvariant_threads = deepvariant_threads,
			deepvariant_version = deepvariant_version
	}

	call bcftools_stats {
		input:
			vcf = deepvariant_postprocess_variants.vcf,
			params = "--apply-filters PASS --samples ~{sample.sample_id}",
			reference = reference.fasta.data,
			container_registry = container_registry
	}

	call bcftools_roh {
		input:
			vcf = deepvariant_postprocess_variants.vcf,
			container_registry = container_registry
	}

	output {
		Array[File] bam_stats = pbmm2_align.bam_stats
		Array[File] read_length_summary = pbmm2_align.read_length_summary
		Array[File] read_quality_summary = pbmm2_align.read_quality_summary
		Array[IndexData] aligned_bams = aligned_bam
		Array[File] aligned_bam_mosdepth_summary = mosdepth.summary
		Array[File] aligned_bam_mosdepth_region_bed = mosdepth.region_bed
		Array[File] svsigs = pbsv_discover.svsig
		IndexData small_variant_vcf = {"data": deepvariant_postprocess_variants.vcf, "data_index": deepvariant_postprocess_variants.vcf_index}
		IndexData small_variant_gvcf = {"data": deepvariant_postprocess_variants.gvcf, "data_index": deepvariant_postprocess_variants.gvcf_index}
		File small_variant_vcf_stats = bcftools_stats.stats
		File small_variant_roh_bed = bcftools_roh.roh_bed
	}

	parameter_meta {
		sample: {help: "Sample information and associated data files"}
		reference: {help: "Reference genome data"}
		deepvariant_version: {help: "Version of deepvariant to use"}
		deepvariant_model: {help: "Optional deepvariant model file to use"}
		container_registry: {help: "Container registry where docker images are hosted"}
	}
}

task pbmm2_align {
	input {
		String sample_id
		File bam

		File reference
		File reference_index
		String reference_name

		String container_registry
	}

	String movie = basename(bam, ".bam")

	Int threads = 24
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		pbmm2 align \
			--num-threads ~{threads} \
			--preset CCS \
			--sample ~{sample_id} \
			--log-level INFO \
			--sort \
			--unmapped \
			-c 0 \
			-y 70 \
			~{reference} \
			~{bam} \
			~{sample_id}.~{movie}.~{reference_name}.aligned.bam

		# movie stats
		extract_read_length_and_qual.py \
			~{bam} \
		> ~{sample_id}.~{movie}.read_length_and_quality.tsv

		awk '{{ b=int($2/1000); b=(b>39?39:b); print 1000*b "\t" $2; }}' \
			~{sample_id}.~{movie}.read_length_and_quality.tsv \
			| sort -k1,1g \
			| datamash -g 1 count 1 sum 2 \
			| awk 'BEGIN {{ for(i=0;i<=39;i++) {{ print 1000*i"\t0\t0"; }} }} {{ print; }}' \
			| sort -k1,1g \
			| datamash -g 1 sum 2 sum 3 \
		> ~{sample_id}.~{movie}.read_length_summary.tsv

		awk '{{ print ($3>50?50:$3) "\t" $2; }}' \
				~{sample_id}.~{movie}.read_length_and_quality.tsv \
			| sort -k1,1g \
			| datamash -g 1 count 1 sum 2 \
			| awk 'BEGIN {{ for(i=0;i<=60;i++) {{ print i"\t0\t0"; }} }} {{ print; }}' \
			| sort -k1,1g \
			| datamash -g 1 sum 2 sum 3 \
		> ~{sample_id}.~{movie}.read_quality_summary.tsv
	>>>

	output {
		File aligned_bam = "~{sample_id}.~{movie}.~{reference_name}.aligned.bam"
		File aligned_bam_index = "~{sample_id}.~{movie}.~{reference_name}.aligned.bam.bai"
		File bam_stats = "~{sample_id}.~{movie}.read_length_and_quality.tsv"
		File read_length_summary = "~{sample_id}.~{movie}.read_length_summary.tsv"
		File read_quality_summary = "~{sample_id}.~{movie}.read_quality_summary.tsv"
	}

	runtime {
		docker: "~{container_registry}/pbmm2:b1a46c6"
		cpu: threads
		memory: "256 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task pbsv_discover {
	input {
		File aligned_bam
		File aligned_bam_index

		File reference_tandem_repeat_bed

		String container_registry
	}

	String prefix = basename(aligned_bam, ".bam")
	Int disk_size = ceil((size(aligned_bam, "GB") + size(reference_tandem_repeat_bed, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		pbsv discover \
			--log-level INFO \
			--hifi \
			--tandem-repeats ~{reference_tandem_repeat_bed} \
			~{aligned_bam} \
			~{prefix}.svsig.gz
	>>>

	output {
		File svsig = "~{prefix}.svsig.gz"
	}

	runtime {
		docker: "~{container_registry}/pbsv:b1a46c6"
		cpu: 4
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task deepvariant_make_examples {
	input {
		String sample_id
		Array[File] aligned_bams
		Array[File] aligned_bam_indices

		File reference
		File reference_index

		Int deepvariant_threads
		String deepvariant_version
	}

	Int disk_size = 500

	command <<<
		set -euo pipefail

		seq 0 ~{deepvariant_threads - 1} \
		| parallel \
			--jobs ~{deepvariant_threads} \
			/opt/deepvariant/bin/make_examples \
				--norealign_reads \
				--vsc_min_fraction_indels 0.12 \
				--pileup_image_width 199 \
				--track_ref_reads \
				--phase_reads \
				--partition_size=25000 \
				--max_reads_per_partition=600 \
				--alt_aligned_pileup=diff_channels \
				--add_hp_channel \
				--sort_by_haplotypes \
				--parse_sam_aux_fields \
				--min_mapping_quality=1 \
				--mode calling \
				--ref ~{reference} \
				--reads ~{sep="," aligned_bams} \
				--examples ~{sample_id}.examples.tfrecord@~{deepvariant_threads}.gz \
				--gvcf ~{sample_id}.gvcf.tfrecord@~{deepvariant_threads}.gz \
				--task {}
	>>>

	output {
		Array[File] example_tfrecords = glob("~{sample_id}.examples.tfrecord*.gz")
		Array[File] nonvariant_site_tfrecords = glob("~{sample_id}.gvcf.tfrecord*.gz")
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: deepvariant_threads
		memory: "256 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task deepvariant_call_variants {
	input {
		String sample_id
		String reference_name
		Array[File] example_tfrecords
		DeepVariantModel? deepvariant_model

		Int deepvariant_threads
		String deepvariant_version
	}

	String deepvariant_model_path = if (defined(deepvariant_model)) then sub(select_first([deepvariant_model]).model.data, "\\.data.*", "") else "/opt/models/pacbio/model.ckpt"
	String example_tfrecord_path = sub(example_tfrecords[0], "/" + basename(example_tfrecords[0]), "")
	Int disk_size = 500

	command <<<
		set -euo pipefail

		/opt/deepvariant/bin/call_variants \
			--outfile ~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz \
			--examples ~{example_tfrecord_path}/~{sample_id}.examples.tfrecord@~{deepvariant_threads}.gz \
			--checkpoint ~{deepvariant_model_path}
	>>>

	output {
		File tfrecord = "~{sample_id}.~{reference_name}.call_variants_output.tfrecord.gz"
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: deepvariant_threads
		memory: "256 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task deepvariant_postprocess_variants {
	input {
		String sample_id
		File tfrecord
		Array[File] nonvariant_site_tfrecords

		File reference
		File reference_index
		String reference_name

		Int deepvariant_threads
		String deepvariant_version
	}

	String nonvariant_site_tfrecord_path = sub(nonvariant_site_tfrecords[0], "/" + basename(nonvariant_site_tfrecords[0]), "")
	Int disk_size = ceil((size(tfrecord, "GB") + size(reference, "GB") + size(nonvariant_site_tfrecords[0], "GB") * length(nonvariant_site_tfrecords)) * 2 + 20)

	command <<<
		set -euo pipefail

		/opt/deepvariant/bin/postprocess_variants \
			--ref ~{reference} \
			--infile ~{tfrecord} \
			--outfile ~{sample_id}.~{reference_name}.deepvariant.vcf.gz \
			--nonvariant_site_tfrecord_path ~{nonvariant_site_tfrecord_path}/~{sample_id}.gvcf.tfrecord@~{deepvariant_threads}.gz \
			--gvcf_outfile ~{sample_id}.~{reference_name}.deepvariant.g.vcf.gz
	>>>

	output {
		File vcf = "~{sample_id}.~{reference_name}.deepvariant.vcf.gz"
		File vcf_index = "~{sample_id}.~{reference_name}.deepvariant.vcf.gz.tbi"
		File gvcf = "~{sample_id}.~{reference_name}.deepvariant.g.vcf.gz"
		File gvcf_index = "~{sample_id}.~{reference_name}.deepvariant.g.vcf.gz.tbi"
		File report = "~{sample_id}.~{reference_name}.deepvariant.visual_report.html"
	}

	runtime {
		docker: "gcr.io/deepvariant-docker/deepvariant:~{deepvariant_version}"
		cpu: 4
		memory: "30 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task bcftools_stats {
	input {
		File vcf
		String? params

		File? reference

		String container_registry
	}

	String vcf_basename = basename(vcf, ".gz")

	Int threads = 4
	Int disk_size = ceil((size(vcf, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		bcftools stats \
			--threads ~{threads - 1} \
			~{params} \
			~{"--fasta-ref " + reference} \
			~{vcf} \
		> ~{vcf_basename}.stats.txt
	>>>

	output {
		File stats = "~{vcf_basename}.stats.txt"
	}

	runtime {
		docker: "~{container_registry}/bcftools:b1a46c6"
		cpu: threads
		memory: "14 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

task bcftools_roh {
	input {
		File vcf

		String container_registry
	}

	String vcf_basename = basename(vcf, ".vcf.gz")

	Int disk_size = ceil(size(vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		echo -e "#chr\tstart\tend\tqual" > ~{vcf_basename}.roh.bed
		bcftools roh \
			--AF-dflt 0.4 \
			~{vcf} \
		| awk -v OFS='\t' '$1=="RG" {{ print $3, $4, $5, $8 }}' \
		>> ~{vcf_basename}.roh.bed
	>>>

	output {
		File roh_bed = "~{vcf_basename}.roh.bed"
	}

	runtime {
		docker: "~{container_registry}/bcftools:b1a46c6"
		cpu: 2
		memory: "4 GB"
		disk: disk_size + " GB"
		preemptible: true
		maxRetries: 3
	}
}

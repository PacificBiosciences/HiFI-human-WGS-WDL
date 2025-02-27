version 1.0

import "../humanwgs_structs.wdl"
import "../wdl-common/wdl/tasks/write_phrank.wdl" as Write_phrank
import "../wdl-common/wdl/tasks/utilities.wdl" as Utilities

workflow tertiary_analysis {
  meta {
    description: "Run tertiary analysis on small and structural variants."
  }

  parameter_meta {
    sample_metadata: {
      name: "PLINK pedigree (PED) formatted lines."
    }
    phenotypes: {
      name: "Comma-delimited list of HPO codes for phenotypes"
    }
    is_trio_kid: {
      name: "Boolean array indicating if the sample is a child with both parents defined"
    }
    is_duo_kid: {
      name: "Boolean array indicating if the sample is a child with only one parent defined"
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
    ref_map_file: {
      name: "Reference map file"
    }
    tertiary_map_file: {
      name: "Tertiary map file"
    }
    default_runtime_attributes: {
      name: "Runtime attribute structure"
    }
    small_variant_filtered_vcf: {
      name: "Filtered and annotated small variant VCF"
    }
    small_variant_filtered_vcf_index: {
      name: "Filtered and annotated small variant VCF index"
    }
    small_variant_filtered_tsv: {
      name: "Filtered and annotated small variant TSV"
    }
    small_variant_compound_het_vcf: {
      name: "Filtered and annotated compound heterozygous small variant VCF"
    }
    small_variant_compound_het_vcf_index: {
      name: "Filtered and annotated compound heterozygous small variant VCF index"
    }
    small_variant_compound_het_tsv: {
      name: "Filtered and annotated compound heterozygous small variant TSV"
    }
    sv_filtered_vcf: {
      name: "Filtered and annotated structural variant VCF"
    }
    sv_filtered_vcf_index: {
      name: "Filtered and annotated structural variant VCF index"
    }
    sv_filtered_tsv: {
      name: "Filtered and annotated structural variant TSV"
    }
  }

  input {
    Array[Array[String]] sample_metadata
    String phenotypes

    Array[Boolean] is_trio_kid  # !UnusedDeclaration
    Array[Boolean] is_duo_kid   # !UnusedDeclaration

    File small_variant_vcf
    File small_variant_vcf_index
    File sv_vcf
    File sv_vcf_index

    File ref_map_file
    File tertiary_map_file

    RuntimeAttributes default_runtime_attributes
  }

  Map[String, String] ref_map      = read_map(ref_map_file)
  Map[String, String] tertiary_map = read_map(tertiary_map_file)

  call Write_phrank.write_phrank {
    input:
      phenotypes         = phenotypes,
      runtime_attributes = default_runtime_attributes
  }

  call Utilities.split_string as split_gnotate_files {
    input:
      concatenated_string = tertiary_map["slivar_gnotate_files"],
      delimiter           = ",",
      runtime_attributes  = default_runtime_attributes
  }

  call Utilities.split_string as split_gnotate_prefixes {
    input:
      concatenated_string = tertiary_map["slivar_gnotate_prefixes"],
      delimiter           = ",",
      runtime_attributes  = default_runtime_attributes
  }

  scatter (gnotate_prefix in split_gnotate_prefixes.array) {
    # These would ideally be within slivar_small_variant, but
    # cromwell doesn't yet support versions of WDL with the suffix function
    # allele frequencies <= max_af in each of the frequency databases
    String slivar_af_expr = "INFO.~{gnotate_prefix}_af <= ~{tertiary_map['slivar_max_af']}"
    # nhomalt <= max_nhomalt in each of the frequency databases
    String slivar_nhomalt_expr = "INFO.~{gnotate_prefix}_nhomalt <= ~{tertiary_map['slivar_max_nhomalt']}"
    # allele counts <= max_ac in each of the frequency databases
    String slivar_ac_expr = "INFO.~{gnotate_prefix}_ac <= ~{tertiary_map['slivar_max_ac']}"
    # info fields for slivar tsv
    Array[String] info_fields = ["~{gnotate_prefix}_af","~{gnotate_prefix}_nhomalt","~{gnotate_prefix}_ac"]
	}

  call slivar_small_variant {
    input:
      vcf                = small_variant_vcf,
      vcf_index          = small_variant_vcf_index,
      sample_metadata    = sample_metadata,
      phrank_lookup      = write_phrank.phrank_lookup,
      reference          = ref_map["fasta"],               # !FileCoercion
      reference_index    = ref_map["fasta_index"],         # !FileCoercion
      gff                = tertiary_map["ensembl_gff"],    # !FileCoercion
      lof_lookup         = tertiary_map["lof_lookup"],     # !FileCoercion
      clinvar_lookup     = tertiary_map["clinvar_lookup"], # !FileCoercion
      slivar_js          = tertiary_map["slivar_js"],      # !FileCoercion
      gnotate_files      = split_gnotate_files.array,      # !FileCoercion
      af_expr            = slivar_af_expr,
      nhomalt_expr       = slivar_nhomalt_expr,
      ac_expr            = slivar_ac_expr,
      info_fields        = flatten(info_fields),
      min_gq             = tertiary_map["slivar_min_gq"],
      runtime_attributes = default_runtime_attributes
  }

  call Utilities.split_string as split_sv_vcfs {
    input:
      concatenated_string = tertiary_map["svpack_pop_vcfs"],
      delimiter           = ",",
      runtime_attributes  = default_runtime_attributes
  }

  call Utilities.split_string as split_sv_vcf_indices {
    input:
      concatenated_string = tertiary_map["svpack_pop_vcf_indices"],
      delimiter           = ",",
      runtime_attributes  = default_runtime_attributes
  }

  call svpack_filter_annotated {
    input:
      sv_vcf                 = sv_vcf,
      sample_metadata        = sample_metadata,
      population_vcfs        = split_sv_vcfs.array,         # !FileCoercion
      population_vcf_indices = split_sv_vcf_indices.array,  # !FileCoercion
      gff                    = tertiary_map["ensembl_gff"], # !FileCoercion
      runtime_attributes     = default_runtime_attributes
  }

  call slivar_svpack_tsv {
    input:
      filtered_vcf       = svpack_filter_annotated.svpack_vcf,
      sample_metadata    = sample_metadata,
      lof_lookup         = tertiary_map["lof_lookup"],         # !FileCoercion
      clinvar_lookup     = tertiary_map["clinvar_lookup"],     # !FileCoercion
      phrank_lookup      = write_phrank.phrank_lookup,
      runtime_attributes = default_runtime_attributes
  }

  output {
    File small_variant_filtered_vcf       = slivar_small_variant.filtered_vcf
    File small_variant_filtered_vcf_index = slivar_small_variant.filtered_vcf_index
    File small_variant_filtered_tsv       = slivar_small_variant.filtered_tsv

    File small_variant_compound_het_vcf       = slivar_small_variant.compound_het_vcf
    File small_variant_compound_het_vcf_index = slivar_small_variant.compound_het_vcf_index
    File small_variant_compound_het_tsv       = slivar_small_variant.compound_het_tsv

    File sv_filtered_vcf       = svpack_filter_annotated.svpack_vcf
    File sv_filtered_vcf_index = svpack_filter_annotated.svpack_vcf_index
    File sv_filtered_tsv       = slivar_svpack_tsv.svpack_tsv
  }
}

task slivar_small_variant {
  meta {
    description: "Filter and annotate small variants with slivar."
  }
  parameter_meta {
    vcf: {
      name: "Small variant VCF"
    }
    vcf_index: {
      name: "Small variant VCF index"
    }
    sample_metadata: {
      name: "PLINK pedigree (PED) formatted lines."
    }
    phrank_lookup: {
      name: "Gene symbol -> Phrank phenotype rank score lookup table"
    }
    reference: {
      name: "Reference genome FASTA"
    }
    reference_index: {
      name: "Reference genome FASTA index"
    }
    gff: {
      name: "Ensembl GFF annotation"
    }
    lof_lookup: {
      name: "Gene symbol -> LoF score lookup table"
    }
    clinvar_lookup: {
      name: "Gene symbol -> ClinVar lookup table"
    }
    slivar_js: {
      name: "Slivar functions JS file"
    }
    gnotate_files: {
      name: "Slivar gnotate files with Allele Frequency (AF), Allele Count (AC), and Number of Homozygotes (nhomalt)"
    }
    af_expr: {
      name: "Allele frequency expressions for slivar"
    }
    nhomalt_expr: {
      name: "nhomalt expressions for slivar"
    }
    ac_expr: {
      name: "Allele count expressions for slivar"
    }
    min_gq: {
      name: "Min genotype quality"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    filtered_vcf: {
      name: "Filtered and annotated small variant VCF"
    }
    filtered_vcf_index: {
      name: "Filtered and annotated small variant VCF index"
    }
    compound_het_vcf: {
      name: "Filtered and annotated compound heterozygous small variant VCF"
    }
    compound_het_vcf_index: {
      name: "Filtered and annotated compound heterozygous small variant VCF index"
    }
    filtered_tsv: {
      name: "Filtered and annotated small variant TSV"
    }
    compound_het_tsv: {
      name: "Filtered and annotated compound heterozygous small variant TSV"
    }
  }

  input {
    File vcf
    File vcf_index

    Array[Array[String]] sample_metadata
    File phrank_lookup

    File reference
    File reference_index

    File gff
    File lof_lookup
    File clinvar_lookup

    File slivar_js
    Array[File] gnotate_files

    Array[String] af_expr
    Array[String] nhomalt_expr
    Array[String] ac_expr
    Array[String] info_fields

    String min_gq

    RuntimeAttributes runtime_attributes
  }

  # First, select only passing variants with AF and nhomalt lower than the specified thresholds
  # The af_expr and nhomalt_expr arrays will be concatenated with this array
  Array[Array[String]] info_expr = [['variant.FILTER=="PASS"'], af_expr, nhomalt_expr]

  # Implicit "high quality" filters from slivar_js are also applied in steps below
  # min_GQ: 20, min_AB: 0.20, min_DP: 6, min_male_X_GQ: 10, min_male_X_DP: 6
  # hom_ref AB < 0.02, hom_alt AB > 0.98, het AB between min_AB and (1-min_AB)

  # Label recessive if all affected samples are HOMALT and all unaffected samples are HETALT or HOMREF
  # Special case of x-linked recessive is also handled, see segregating_recessive_x in slivar docs
  Array[String] family_recessive_expr = ['recessive:fam.every(segregating_recessive)']

  # Label dominant if all affected samples are HETALT and all unaffected samples are HOMREF
  # Special case of x-linked dominant is also handled, see segregating_dominant_x in slivar docs
  # The ac_expr array will be concatenated with this array
  Array[Array[String]] family_dominant_expr = [['dominant:fam.every(segregating_dominant)'], ac_expr]

  # Label comphet_side if the sample is HETALT and the GQ is above the specified threshold
  Array[String] sample_expr = [
    'comphet_side:sample.het',
    'sample.GQ > ~{min_gq}'
  ]

  # Skip these variant types when looking for compound hets
  Array[String] skip_list = [
    'non_coding_transcript',
    'intron',
    'non_coding',
    'upstream_gene',
    'downstream_gene',
    'non_coding_transcript_exon',
    'NMD_transcript',
    '5_prime_UTR',
    '3_prime_UTR'
  ]

  String vcf_basename = basename(vcf, ".vcf.gz")

  Int threads   = 8
  Int mem_gb    = 2 * threads
  Int disk_size = ceil((size(vcf, "GB") + size(reference, "GB") + size(gnotate_files, "GB") + size(gff, "GB") + size(lof_lookup, "GB") + size(clinvar_lookup, "GB") + size(phrank_lookup, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    bcftools --version

    bcftools norm \
      --threads ~{threads - 1} \
      --multiallelics \
      - \
      --output-type b \
      --fasta-ref ~{reference} \
      ~{vcf} \
    | bcftools sort \
      --output-type b \
      --output ~{vcf_basename}.norm.bcf

    bcftools index \
      --threads ~{threads - 1} \
      ~{vcf_basename}.norm.bcf

    # slivar has no version option
    slivar expr 2>&1 | grep -Eo 'slivar version: [0-9.]+ [0-9a-f]+' 

    pslivar \
      --processes ~{threads} \
      --fasta ~{reference} \
      --pass-only \
      --js ~{slivar_js} \
      --info '~{sep=" && " flatten(info_expr)}' \
      --family-expr '~{sep=" && " family_recessive_expr}' \
      --family-expr '~{sep=" && " flatten(family_dominant_expr)}' \
      --sample-expr '~{sep=" && " sample_expr}' \
      ~{sep=" " prefix("--gnotate ", gnotate_files)} \
      --vcf ~{vcf_basename}.norm.bcf \
      --ped ~{write_tsv(sample_metadata)} \
    | bcftools csq \
      --local-csq \
      --samples - \
      --ncsq 40 \
      --gff-annot ~{gff} \
      --fasta-ref ~{reference} \
      - \
      --output-type z \
      --output ~{vcf_basename}.norm.slivar.vcf.gz

    bcftools index \
      --threads ~{threads - 1} \
      --tbi ~{vcf_basename}.norm.slivar.vcf.gz

    slivar \
      compound-hets \
      --skip ~{sep=',' skip_list} \
      --vcf ~{vcf_basename}.norm.slivar.vcf.gz \
      --sample-field comphet_side \
      --ped ~{write_tsv(sample_metadata)} \
      --allow-non-trios \
    | add_comphet_phase.py \
    | bcftools view \
      --output-type z \
      --output ~{vcf_basename}.norm.slivar.compound_hets.vcf.gz

    bcftools index \
      --threads ~{threads - 1} \
      --tbi ~{vcf_basename}.norm.slivar.compound_hets.vcf.gz

    slivar tsv \
      --info-field ~{sep=' --info-field ' info_fields} \
      --sample-field dominant \
      --sample-field recessive \
      --csq-field BCSQ \
      --gene-description ~{lof_lookup} \
      --gene-description ~{clinvar_lookup} \
      --gene-description ~{phrank_lookup} \
      --ped ~{write_tsv(sample_metadata)} \
      --out /dev/stdout \
      ~{vcf_basename}.norm.slivar.vcf.gz \
    | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
    > ~{vcf_basename}.norm.slivar.tsv

    slivar tsv \
      --info-field ~{sep=' --info-field ' info_fields} \
      --sample-field slivar_comphet \
      --info-field slivar_comphet \
      --csq-field BCSQ \
      --gene-description ~{lof_lookup} \
      --gene-description ~{clinvar_lookup} \
      --gene-description ~{phrank_lookup} \
      --ped ~{write_tsv(sample_metadata)} \
      --out /dev/stdout \
      ~{vcf_basename}.norm.slivar.compound_hets.vcf.gz \
    | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
    > ~{vcf_basename}.norm.slivar.compound_hets.tsv
  >>>

  output {
    File filtered_vcf           = "~{vcf_basename}.norm.slivar.vcf.gz"
    File filtered_vcf_index     = "~{vcf_basename}.norm.slivar.vcf.gz.tbi"
    File compound_het_vcf       = "~{vcf_basename}.norm.slivar.compound_hets.vcf.gz"
    File compound_het_vcf_index = "~{vcf_basename}.norm.slivar.compound_hets.vcf.gz.tbi"
    File filtered_tsv           = "~{vcf_basename}.norm.slivar.tsv"
    File compound_het_tsv       = "~{vcf_basename}.norm.slivar.compound_hets.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/slivar@sha256:f71a27f756e2d69ec30949cbea97c54abbafde757562a98ef965f21a28aa8eaa"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
  }
}

task svpack_filter_annotated {
  meta {
    description: "Filter and annotate structural variants with svpack."
  }

  parameter_meta {
    sample_metadata: {
      name: "PLINK pedigree (PED) formatted lines."
    }
    sv_vcf: {
      name: "Structural variant VCF"
    }
    population_vcfs: {
      name: "SV population VCFs"
    }
    population_vcf_indices: {
      name: "SV population VCF indices"
    }
    gff: {
      name: "Ensembl GFF annotation"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    svpack_vcf: {
      name: "Filtered and annotated structural variant VCF"
    }
    svpack_vcf_index: {
      name: "Filtered and annotated structural variant VCF index"
    }
  }

  input {
    File sv_vcf
    Array[Array[String]] sample_metadata

    Array[File] population_vcfs
    Array[File] population_vcf_indices

    File gff

    RuntimeAttributes runtime_attributes
  }

  Int threads   = 2
  Int mem_gb    = 16
  Int disk_size = ceil(size(sv_vcf, "GB") * 2 + 20)

  String out_prefix = basename(sv_vcf, ".vcf.gz")

  command <<<
    echo "svpack version:"
    cat /opt/svpack/.git/HEAD

    affected=$(awk -F'\t' '$6 ~ /2/ {{ print $2 }}' ~{write_tsv(sample_metadata)} | paste -sd',')  # TODO: potentially replace awk

    svpack \
      filter \
      --pass-only \
      --min-svlen 50 \
      ~{sv_vcf} \
    ~{sep=' ' prefix('| svpack match -v - ', population_vcfs)} \
    | svpack \
      consequence \
      - \
      <(zcat ~{gff} || cat ~{gff}) \
    | svpack \
      tagzygosity \
      --samples "${affected}" \
      - \
    > ~{out_prefix}.svpack.vcf

    bgzip --version

    bgzip ~{out_prefix}.svpack.vcf

    tabix --version

    tabix --preset vcf ~{out_prefix}.svpack.vcf.gz
  >>>

  output {
    File svpack_vcf       = "~{out_prefix}.svpack.vcf.gz"
    File svpack_vcf_index = "~{out_prefix}.svpack.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/svpack@sha256:628e9851e425ed8044a907d33de04043d1ef02d4d2b2667cf2e9a389bb011eba"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
  }
}

task slivar_svpack_tsv {
  meta {
    description: "Create spreadsheet-friendly TSV from svpack annotated VCFs."
  }

  parameter_meta {
    filtered_vcf : {
      name: "Filtered and annotated structural variant VCF"
    }
    sample_metadata: {
      name: "PLINK pedigree (PED) formatted lines."
    }
    lof_lookup: {
      name: "Gene symbol -> LoF score lookup table"
    }
    clinvar_lookup: {
      name: "Gene symbol -> ClinVar lookup table"
    }
    phrank_lookup: {
      name: "Gene symbol -> Phrank phenotype rank score lookup table"
    }
    runtime_attributes: {
      name: "Runtime attribute structure"
    }
    svpack_tsv: {
      name: "Filtered and annotated structural variant TSV"
    }
  }

  input {
    File filtered_vcf

    Array[Array[String]] sample_metadata
    File lof_lookup
    File clinvar_lookup
    File phrank_lookup

    RuntimeAttributes runtime_attributes
  }

  Array[String] info_fields = [
    'SVTYPE',
    'SVLEN',
    'MATEID',
    'END'
  ]

  String filtered_vcf_basename = basename(filtered_vcf, ".vcf.gz")

  Int threads   = 2
  Int mem_gb    = 4
  Int disk_size = ceil((size(filtered_vcf, "GB") + size(lof_lookup, "GB") + size(clinvar_lookup, "GB") + size(phrank_lookup, "GB")) * 2 + 20)

  command <<<
    set -euo pipefail

    # slivar has no version option
    slivar expr 2>&1 | grep -Eo 'slivar version: [0-9.]+ [0-9a-f]+'

    slivar tsv \
      --info-field ~{sep=' --info-field ' info_fields} \
      --sample-field hetalt \
      --sample-field homalt \
      --csq-field BCSQ \
      --gene-description ~{lof_lookup} \
      --gene-description ~{clinvar_lookup} \
      --gene-description ~{phrank_lookup} \
      --ped ~{write_tsv(sample_metadata)} \
      --out /dev/stdout \
      ~{filtered_vcf} \
    | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
    > ~{filtered_vcf_basename}.tsv
  >>>

  output {
    File svpack_tsv = "~{filtered_vcf_basename}.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/slivar@sha256:f71a27f756e2d69ec30949cbea97c54abbafde757562a98ef965f21a28aa8eaa"
    cpu: threads
    memory: mem_gb + " GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
  }
}

version 1.0

import "../humanwgs_structs.wdl"
import "../wdl-common/wdl/tasks/split_string.wdl" as Split_string

workflow tertiary_analysis {
  meta {
    description: "Run tertiary analysis on small and structural variants."
  }

  parameter_meta {
    pedigree: {
      name: "PLINK pedigree (PED) format"
    }
    phrank_lookup: {
      name: "Gene symbol -> Phrank phenotype rank score lookup table"
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
    filtered_small_variant_vcf: {
      name: "Filtered and annotated small variant VCF"
    }
    filtered_small_variant_vcf_index: {
      name: "Filtered and annotated small variant VCF index"
    }
    filtered_small_variant_tsv: {
      name: "Filtered and annotated small variant TSV"
    }
    compound_het_small_variant_vcf: {
      name: "Filtered and annotated compound heterozygous small variant VCF"
    }
    compound_het_small_variant_vcf_index: {
      name: "Filtered and annotated compound heterozygous small variant VCF index"
    }
    compound_het_small_variant_tsv: {
      name: "Filtered and annotated compound heterozygous small variant TSV"
    }
    filtered_svpack_vcf: {
      name: "Filtered and annotated structural variant VCF"
    }
    filtered_svpack_vcfdv: {
      name: "Filtered and annotated structural variant VCF index"
    }
    filtered_svpack_tsv: {
      name: "Filtered and annotated structural variant TSV"
    }
  }

  input {
    File pedigree
    File phrank_lookup

    File small_variant_vcf
    File small_variant_vcf_index
    File sv_vcf
    File sv_vcf_index

    File ref_map_file
    File tertiary_map_file

    RuntimeAttributes default_runtime_attributes
  }

  Map[String, String] ref_map       = read_map(ref_map_file)
  Map[String, String] tertiary_map = read_map(tertiary_map_file)

  call slivar_small_variant {
    input:
      vcf                = small_variant_vcf,
      vcf_index          = small_variant_vcf_index,
      pedigree           = pedigree,
      reference          = ref_map["fasta"],                    # !FileCoercion
      reference_index    = ref_map["fasta_index"],              # !FileCoercion
      slivar_js          = tertiary_map["slivar_js"],           # !FileCoercion
      gnomad_af          = tertiary_map["slivar_gnomad_af"],    # !FileCoercion
      hprc_af            = tertiary_map["slivar_hprc_af"],      # !FileCoercion
      gff                = tertiary_map["ensembl_gff"],         # !FileCoercion
      lof_lookup         = tertiary_map["lof_lookup"],          # !FileCoercion
      clinvar_lookup     = tertiary_map["clinvar_lookup"],      # !FileCoercion
      phrank_lookup      = phrank_lookup,
      runtime_attributes = default_runtime_attributes
  }

  call Split_string.split_string as split_sv_vcfs {
    input:
      concatenated_string = tertiary_map["svpack_pop_vcfs"],
      delimiter           = ",",
      runtime_attributes  = default_runtime_attributes
  }

  call Split_string.split_string as split_sv_vcf_indices {
    input:
      concatenated_string = tertiary_map["svpack_pop_vcf_indices"],
      delimiter           = ",",
      runtime_attributes  = default_runtime_attributes
  }

  call svpack_filter_annotated {
    input:
      sv_vcf                 = sv_vcf,
      population_vcfs        = split_sv_vcfs.array,         # !FileCoercion
      population_vcf_indices = split_sv_vcf_indices.array,  # !FileCoercion
      gff                    = tertiary_map["ensembl_gff"], # !FileCoercion
      runtime_attributes     = default_runtime_attributes
  }

  call slivar_svpack_tsv {
    input:
      filtered_vcf       = svpack_filter_annotated.svpack_vcf,
      pedigree           = pedigree,
      lof_lookup         = tertiary_map["lof_lookup"],         # !FileCoercion
      clinvar_lookup     = tertiary_map["clinvar_lookup"],     # !FileCoercion
      phrank_lookup      = phrank_lookup,
      runtime_attributes = default_runtime_attributes
  }

  output {
    File filtered_small_variant_vcf       = slivar_small_variant.filtered_vcf
    File filtered_small_variant_vcf_index = slivar_small_variant.filtered_vcf_index
    File filtered_small_variant_tsv       = slivar_small_variant.filtered_tsv

    File compound_het_small_variant_vcf       = slivar_small_variant.compound_het_vcf
    File compound_het_small_variant_vcf_index = slivar_small_variant.compound_het_vcf_index
    File compound_het_small_variant_tsv       = slivar_small_variant.compound_het_tsv

    File filtered_svpack_vcf       = svpack_filter_annotated.svpack_vcf
    File filtered_svpack_vcf_index = svpack_filter_annotated.svpack_vcf_index
    File filtered_svpack_tsv       = slivar_svpack_tsv.svpack_tsv
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
    pedigree: {
      name: "PLINK pedigree (PED) format"
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
    gnomad_af: {
      name: "gnomAD gnotate file"
    }
    hprc_af: {
      name: "HPRC gnotate"
    }
    max_gnomad_af: {
      name: "Max gnomAD allele frequency"
    }
    max_hprc_af: {
      name: "Max HPRC allele frequency"
    }
    max_gnomad_nhomalt: {
      name: "Max gnomAD count of HOMALT alleles"
    }
    max_hprc_nhomalt: {
      name: "Max HPRC count of HOMALT alleles"
    }
    max_gnomad_ac: {
      name: "Max gnomAD allele count"
    }
    max_hprc_ac: {
      name: "Max HPRC allele count"
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

    File pedigree
    File phrank_lookup

    File reference
    File reference_index

    File gff
    File lof_lookup
    File clinvar_lookup

    File slivar_js
    File gnomad_af
    File hprc_af

    Float max_gnomad_af      = 0.03
    Float max_hprc_af        = 0.03
    Int   max_gnomad_nhomalt = 4
    Int   max_hprc_nhomalt   = 4
    Int   max_gnomad_ac      = 4
    Int   max_hprc_ac        = 4
    Int   min_gq             = 5

    RuntimeAttributes runtime_attributes
  }

  # First, select only passing variants with AF and nhomalt lower than the specified thresholds
  Array[String] info_expr = [
    'variant.FILTER=="PASS"',
    'INFO.gnomad_af <= ~{max_gnomad_af}',
    'INFO.hprc_af <= ~{max_hprc_af}',
    'INFO.gnomad_nhomalt <= ~{max_gnomad_nhomalt}',
    'INFO.hprc_nhomalt <= ~{max_hprc_nhomalt}'
  ]

  # Implicit "high quality" filters are also applied in steps below
  # min_GQ: 20, min_AB: 0.20, min_DP: 6, min_male_X_GQ: 10, min_male_X_DP: 6
  # hom_ref AB < 0.02, hom_alt AB > 0.98, het AB between min_AB and (1-min_AB)

  # Label recessive if all affected samples are HOMALT and all unaffected samples are HETALT or HOMREF
  # Special case of x-linked recessive is also handled, see segregating_recessive_x in slivar docs
  Array[String] family_recessive_expr = [
    'recessive:fam.every(segregating_recessive)'
  ]

  # Label dominant if all affected samples are HETALT and all unaffected samples are HOMREF
  # Special case of x-linked dominant is also handled, see segregating_dominant_x in slivar docs
  Array[String] family_dominant_expr = [
    'dominant:fam.every(segregating_dominant)',
    'INFO.gnomad_ac <= ~{max_gnomad_ac}',
    'INFO.hprc_ac <= ~{max_hprc_ac}'
  ]

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

  # Fields to include in the output TSV
  Array[String] info_fields = [
    'gnomad_af',
    'hprc_af',
    'gnomad_nhomalt',
    'hprc_nhomalt',
    'gnomad_ac',
    'hprc_ac'
  ]

  String vcf_basename = basename(vcf, ".vcf.gz")
  Int    threads      = 8
  Int    disk_size    = ceil((size(vcf, "GB") + size(reference, "GB") + size(gnomad_af, "GB") + size(hprc_af, "GB") + size(gff, "GB") + size(lof_lookup, "GB") + size(clinvar_lookup, "GB") + size(phrank_lookup, "GB")) * 2 + 20)

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
      --info '~{sep=" && " info_expr}' \
      --family-expr '~{sep=" && " family_recessive_expr}' \
      --family-expr '~{sep=" && " family_dominant_expr}' \
      --sample-expr '~{sep=" && " sample_expr}' \
      --gnotate ~{gnomad_af} \
      --gnotate ~{hprc_af} \
      --vcf ~{vcf_basename}.norm.bcf \
      --ped ~{pedigree} \
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
      --ped ~{pedigree} \
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
      --ped ~{pedigree} \
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
      --ped ~{pedigree} \
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
    docker: "~{runtime_attributes.container_registry}/slivar@sha256:0a09289ccb760da310669906c675be02fd16b18bbedc971605a587275e34966c"
    cpu: threads
    memory: "16 GB"
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

    Array[File] population_vcfs
    Array[File] population_vcf_indices

    File gff

    RuntimeAttributes runtime_attributes
  }

  String sv_vcf_basename = basename(sv_vcf, ".vcf.gz")
  Int    disk_size       = ceil(size(sv_vcf, "GB") * 2 + 20)

  command <<<
    set -euo pipefail

    echo "svpack version:"
    cat /opt/svpack/.git/HEAD

    svpack \
      filter \
      --pass-only \
      --min-svlen 50 \
      ~{sv_vcf} \
    ~{sep=' ' prefix('| svpack match -v - ', population_vcfs)} \
    | svpack \
      consequence \
      - \
      ~{gff} \
    | svpack \
      tagzygosity \
      - \
    > ~{sv_vcf_basename}.svpack.vcf

    bgzip --version

    bgzip ~{sv_vcf_basename}.svpack.vcf

    tabix --version

    tabix -p vcf ~{sv_vcf_basename}.svpack.vcf.gz
  >>>

  output {
    File svpack_vcf       = "~{sv_vcf_basename}.svpack.vcf.gz"
    File svpack_vcf_index = "~{sv_vcf_basename}.svpack.vcf.gz.tbi"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/svpack@sha256:a680421cb517e1fa4a3097838719a13a6bd655a5e6980ace1b03af9dd707dd75"
    cpu: 2
    memory: "16 GB"
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
    pedigree: {
      name: "PLINK pedigree (PED) format"
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

    File pedigree
    File lof_lookup
    File clinvar_lookup
    File phrank_lookup

    RuntimeAttributes runtime_attributes
  }

  Array[String] info_fields = [
    'SVTYPE',
    'SVLEN',
    'SVANN',
    'CIPOS',
    'MATEID',
    'END'
  ]

  String filtered_vcf_basename = basename(filtered_vcf, ".vcf.gz")
  Int    disk_size             = ceil((size(filtered_vcf, "GB") + size(lof_lookup, "GB") + size(clinvar_lookup, "GB") + size(phrank_lookup, "GB")) * 2 + 20)

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
      --ped ~{pedigree} \
      --out /dev/stdout \
      ~{filtered_vcf} \
    | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
    > ~{filtered_vcf_basename}.tsv
  >>>

  output {
    File svpack_tsv = "~{filtered_vcf_basename}.tsv"
  }

  runtime {
    docker: "~{runtime_attributes.container_registry}/slivar@sha256:0a09289ccb760da310669906c675be02fd16b18bbedc971605a587275e34966c"
    cpu: 2
    memory: "4 GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    zones: runtime_attributes.zones
  }
}

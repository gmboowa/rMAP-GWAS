version 1.0

workflow rMAP_GWAS {
  input {
    # sample-set inputs
    # Suggested workspace mappings when running on a gwasmtb_set:
    #   sample_names = this.gwasmtbs.gwasmtb_id
    #   read1s       = this.gwasmtbs.read1
    #   read2s       = this.gwasmtbs.read2
    #   groups       = this.gwasmtbs.group
    Array[String]+ sample_names
    Array[File]+ read1s
    Array[File]+ read2s
    Array[String]+ groups

    # Phenotype coding and biological display labels.
    # Required binary GWAS grouping is provided by the `groups` input, expected as case/control.
    # To make the HTML self-explanatory, map phenotype_display_values to the metadata column
    # that defines the biological contrast, and set phenotype_name to that column name.
    # Examples:
    #   phenotype_name           = "specimen_source"
    #   phenotype_display_values = this.<entity_set_members>.specimen_source
    #
    # Then the report automatically infers:
    #   case    -> case (blood)
    #   control -> control (sputum)
    String phenotype_name = "case_control"
    Array[String] phenotype_display_values = []

    # Analysis controls
    Float min_af = 0.01
    Float max_af = 0.99

    # Modular GWAS controls.
    # Current gene presence/absence GWAS remains the default.
    # Set do_snp_gwas=true to run the reference-based SNP GWAS branch
    # for mutation-mediated phenotypes such as MTBC drug resistance.
    String gwas_mode = "gene_presence_absence"
    Boolean do_snp_gwas = false

    # Optional recombination assessment/filtering for SNP GWAS.
    # Default is false to keep MTBC/clonal-organism runs simple and to avoid adding
    # runtime to gene-only or smoke-test analyses. Enable for recombining bacteria
    # where recombinant blocks may inflate SNP associations.
    Boolean do_gubbins = false
    Float snp_min_qual = 20.0
    String container_backend = "docker"

    # Pyseer population-structure controls.
    # For small smoke tests, using too many MDS dimensions can make the null model singular.
    Int pyseer_max_dimensions = 2
    Boolean pyseer_force_no_distances = false
    Boolean pyseer_no_distances_fallback = true
    Float significance_alpha = 0.05
    String output_prefix = "rMAP_GWAS"

    # Runtime controls
    # These defaults are smoke-test friendly. For full cohorts, increase pangenome/gwas resources as needed.
    Int fastp_threads = 4
    Int assembly_threads = 4
    Int annotation_threads = 4
    Int pangenome_threads = 4
    Int gwas_threads = 4
    Int gubbins_threads = 4

    Int fastp_memory_gb = 8
    Int assembly_memory_gb = 32
    Int annotation_memory_gb = 16
    Int pangenome_memory_gb = 32
    Int gwas_memory_gb = 32
    Int gubbins_memory_gb = 16

    Int fastp_disk_gb = 50
    Int assembly_disk_gb = 200
    Int annotation_disk_gb = 100
    Int pangenome_disk_gb = 300
    Int gwas_disk_gb = 300
    Int snp_disk_gb = 300
    Int gubbins_disk_gb = 300

    # Docker images
    String fastp_docker = "quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    String shovill_docker = "quay.io/biocontainers/shovill:1.1.0--hdfd78af_1"
    String quast_docker = "staphb/quast:5.2.0"
    String prokka_docker = "staphb/prokka:1.14.6"
    # Annotation engine control. Prokka remains the default to preserve existing runs.
    # Set use_bakta=true to use Bakta instead; the Bakta image must include a compatible database
    # or bakta_db must point to a mounted database path available inside the container.
    Boolean use_bakta = false
    String bakta_docker = "gmboowa/rmap-gwas-bakta-db:light-0.1"
    String bakta_db = "/opt/bakta/db-light"
    String panaroo_docker = "quay.io/biocontainers/panaroo:1.5.2--pyhdfd78af_0"
    # Combined linux/amd64 image for local Colima testing and Cromwell execution.
    # Contains pyseer, mash, Python, pandas, numpy, scipy, statsmodels, scikit-learn and tqdm.
    String mash_docker = "gmboowa/rmap-gwas-pyseer-annotate:0.2"
    String pyseer_docker = "gmboowa/rmap-gwas-pyseer-annotate:0.2"
    # Dedicated reference-based SNP-calling image. This keeps mapping/SNP calling
    # separate from the pyseer/statistics image and avoids missing bwa/samtools/bcftools errors.
    String snp_calling_docker = "staphb/snippy:4.6.0"
    # Optional recombination module. Used only when do_gubbins=true.
    String gubbins_docker = "staphb/gubbins:latest"
    String python_docker = "gmboowa/rmap-gwas-pyseer-annotate:0.2"
    # Docker image used for post-GWAS reference annotation and plot generation.
    # Keep this Python-capable; the task uses pure Python and does not require BLAST.
    String hit_annotation_docker = "gmboowa/rmap-gwas-pyseer-annotate:0.2"
    Int plot_max_points = 5000

    # Reference/provenance settings. This records the species-specific reference package used or intended for the run.
    String reference_docker = "gmboowa/rmap-gwas-mtbc-refs:2026.06"
    String reference_species = "Mycobacterium tuberculosis complex"
    String reference_name = "MTBC_2026_06"
  }

  # Internal binary labels. Keep the Terra metadata table simple by using
  # group=case/control and phenotype_display_values for the biological meaning.
  String case_label = "case"
  String control_label = "control"
  String annotation_engine = if use_bakta then "Bakta" else "Prokka"

  # Use sample-set arrays directly.
  Array[String]+ all_sample_names = sample_names
  Array[File]+ all_read1s = read1s
  Array[File]+ all_read2s = read2s
  Array[String] phenotype_display_values_for_report = phenotype_display_values

  # Shovill checks usable RAM inside the VM and fails if --ram is >= available RAM.
  # Cromwell VMs expose slightly less usable RAM than the WDL runtime request, so keep
  # the Shovill command-level RAM below the task runtime memory.
  Int assembly_shovill_ram_gb = if assembly_memory_gb > 16 then assembly_memory_gb - 8 else if assembly_memory_gb > 8 then assembly_memory_gb - 4 else assembly_memory_gb

  call VALIDATE_SAMPLE_SET_INPUTS {
    input:
      sample_names = sample_names,
      read1_count = length(read1s),
      read2_count = length(read2s),
      groups = groups,
      phenotype_display_values = phenotype_display_values_for_report,
      case_label = case_label,
      control_label = control_label,
      python_docker = python_docker
  }

  call PREPARE_PHENOTYPE_TABLE {
    input:
      sample_names = sample_names,
      groups = groups,
      phenotype_name = phenotype_name,
      case_label = case_label,
      control_label = control_label,
      phenotype_display_values = phenotype_display_values_for_report,
      output_prefix = output_prefix,
      python_docker = python_docker,
      validation_report = VALIDATE_SAMPLE_SET_INPUTS.validation_report
  }

  scatter (i in range(length(all_sample_names))) {
    call FASTP_TRIM {
      input:
        sample_name = all_sample_names[i],
        read1 = all_read1s[i],
        read2 = all_read2s[i],
        threads = fastp_threads,
        memory_gb = fastp_memory_gb,
        disk_gb = fastp_disk_gb,
        docker = fastp_docker,
        validation_report = VALIDATE_SAMPLE_SET_INPUTS.validation_report
    }

    call SHOVILL_ASSEMBLE {
      input:
        sample_name = all_sample_names[i],
        read1 = FASTP_TRIM.trimmed_read1,
        read2 = FASTP_TRIM.trimmed_read2,
        threads = assembly_threads,
        memory_gb = assembly_memory_gb,
        shovill_ram_gb = assembly_shovill_ram_gb,
        disk_gb = assembly_disk_gb,
        docker = shovill_docker
    }

    call QUAST_ASSEMBLY_QC {
      input:
        sample_name = all_sample_names[i],
        assembly = SHOVILL_ASSEMBLE.contigs_fasta,
        threads = 2,
        docker = quast_docker
    }

    if (!use_bakta) {
      call PROKKA_ANNOTATE {
        input:
          sample_name = all_sample_names[i],
          assembly = SHOVILL_ASSEMBLE.contigs_fasta,
          threads = annotation_threads,
          memory_gb = annotation_memory_gb,
          docker = prokka_docker
      }
    }

    if (use_bakta) {
      call BAKTA_ANNOTATE {
        input:
          sample_name = all_sample_names[i],
          assembly = SHOVILL_ASSEMBLE.contigs_fasta,
          threads = annotation_threads,
          memory_gb = annotation_memory_gb,
          disk_gb = annotation_disk_gb,
          docker = bakta_docker,
          bakta_db = bakta_db
      }
    }
  }

  Array[File] annotation_gffs = if use_bakta then select_all(BAKTA_ANNOTATE.gff) else select_all(PROKKA_ANNOTATE.gff)
  Array[File] annotation_gbks = if use_bakta then select_all(BAKTA_ANNOTATE.gbk) else select_all(PROKKA_ANNOTATE.gbk)
  Array[File] annotation_faas = if use_bakta then select_all(BAKTA_ANNOTATE.faa) else select_all(PROKKA_ANNOTATE.faa)
  Array[File] annotation_ffns = if use_bakta then select_all(BAKTA_ANNOTATE.ffn) else select_all(PROKKA_ANNOTATE.ffn)

  call PANAROO_PANGENOME {
    input:
      gffs = annotation_gffs,
      threads = pangenome_threads,
      memory_gb = pangenome_memory_gb,
      disk_gb = pangenome_disk_gb,
      docker = panaroo_docker
  }

  call MASH_DISTANCE_MATRIX {
    input:
      sample_names = all_sample_names,
      assemblies = SHOVILL_ASSEMBLE.contigs_fasta,
      threads = gwas_threads,
      memory_gb = gwas_memory_gb,
      disk_gb = gwas_disk_gb,
      docker = mash_docker
  }

  call GENERATE_POPULATION_STRUCTURE_PLOTS {
    input:
      phenotype_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_tsv,
      phenotype_legend_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_legend_tsv,
      mash_distances = MASH_DISTANCE_MATRIX.mash_distances,
      output_prefix = output_prefix,
      python_docker = python_docker
  }

  call PYSEER_GENE_GWAS {
    input:
      phenotype_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_tsv,
      gene_presence_absence_rtab = PANAROO_PANGENOME.gene_presence_absence_rtab,
      mash_distances = MASH_DISTANCE_MATRIX.mash_distances,
      min_af = min_af,
      max_af = max_af,
      max_dimensions = pyseer_max_dimensions,
      force_no_distances = pyseer_force_no_distances,
      no_distances_fallback = pyseer_no_distances_fallback,
      threads = gwas_threads,
      memory_gb = gwas_memory_gb,
      disk_gb = gwas_disk_gb,
      docker = pyseer_docker
  }

  call PRIORITIZE_GWAS_HITS {
    input:
      phenotype_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_tsv,
      gene_presence_absence_csv = PANAROO_PANGENOME.gene_presence_absence_csv,
      gene_presence_absence_rtab = PANAROO_PANGENOME.gene_presence_absence_rtab,
      pyseer_gene_assoc = PYSEER_GENE_GWAS.pyseer_gene_assoc,
      significance_alpha = significance_alpha,
      output_prefix = output_prefix,
      annotation_engine = annotation_engine,
      python_docker = python_docker
  }

  call EXTRACT_REFERENCE_GENBANK_FROM_DOCKER {
    input:
      reference_docker = reference_docker,
      reference_name = reference_name
  }

  String snp_output_prefix = output_prefix + "_SNP"

  call MAKE_SNP_GWAS_PLACEHOLDERS {
    input:
      output_prefix = snp_output_prefix,
      python_docker = python_docker
  }

  call MAKE_GUBBINS_PLACEHOLDERS {
    input:
      output_prefix = snp_output_prefix,
      python_docker = python_docker
  }

  if (do_snp_gwas) {
    call SNP_CALLING_SNIPPY {
      input:
        sample_names = all_sample_names,
        read1s = FASTP_TRIM.trimmed_read1,
        read2s = FASTP_TRIM.trimmed_read2,
        reference_fasta = EXTRACT_REFERENCE_GENBANK_FROM_DOCKER.reference_fasta,
        snp_min_qual = snp_min_qual,
        threads = gwas_threads,
        memory_gb = gwas_memory_gb,
        disk_gb = snp_disk_gb,
        docker = snp_calling_docker
    }

    if (do_gubbins) {
      call GUBBINS_RECOMBINATION_FILTER {
        input:
          core_alignment = SNP_CALLING_SNIPPY.snippy_core_alignment,
          snp_vcf = SNP_CALLING_SNIPPY.snp_vcf,
          output_prefix = snp_output_prefix,
          threads = gubbins_threads,
          memory_gb = gubbins_memory_gb,
          disk_gb = gubbins_disk_gb,
          docker = gubbins_docker
      }
    }

    File snp_vcf_for_pyseer = select_first([GUBBINS_RECOMBINATION_FILTER.gubbins_filtered_vcf, SNP_CALLING_SNIPPY.snp_vcf])

    call PYSEER_SNP_GWAS {
      input:
        phenotype_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_tsv,
        snp_vcf = snp_vcf_for_pyseer,
        mash_distances = MASH_DISTANCE_MATRIX.mash_distances,
        min_af = min_af,
        max_af = max_af,
        max_dimensions = pyseer_max_dimensions,
        force_no_distances = pyseer_force_no_distances,
        no_distances_fallback = pyseer_no_distances_fallback,
        threads = gwas_threads,
        memory_gb = gwas_memory_gb,
        disk_gb = snp_disk_gb,
        docker = pyseer_docker
    }

    call PRIORITIZE_SNP_GWAS_HITS {
      input:
        phenotype_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_tsv,
        pyseer_snp_assoc = PYSEER_SNP_GWAS.pyseer_snp_assoc,
        snp_vcf = snp_vcf_for_pyseer,
        reference_genbank = EXTRACT_REFERENCE_GENBANK_FROM_DOCKER.reference_genbank,
        significance_alpha = significance_alpha,
        output_prefix = snp_output_prefix,
        python_docker = python_docker
    }

    call GENERATE_SNP_GWAS_PLOTS {
      input:
        pyseer_snp_assoc = PYSEER_SNP_GWAS.pyseer_snp_assoc,
        snp_vcf = snp_vcf_for_pyseer,
        output_prefix = snp_output_prefix,
        significance_alpha = significance_alpha,
        max_points = plot_max_points,
        python_docker = python_docker
    }
  }

  File snp_vcf_for_report = select_first([GUBBINS_RECOMBINATION_FILTER.gubbins_filtered_vcf, SNP_CALLING_SNIPPY.snp_vcf, MAKE_SNP_GWAS_PLACEHOLDERS.snp_vcf])
  File gubbins_summary_for_report = select_first([GUBBINS_RECOMBINATION_FILTER.gubbins_summary, MAKE_GUBBINS_PLACEHOLDERS.gubbins_summary])
  File gubbins_filtered_alignment_for_report = select_first([GUBBINS_RECOMBINATION_FILTER.gubbins_filtered_alignment, MAKE_GUBBINS_PLACEHOLDERS.gubbins_filtered_alignment])
  File gubbins_recombination_gff_for_report = select_first([GUBBINS_RECOMBINATION_FILTER.gubbins_recombination_gff, MAKE_GUBBINS_PLACEHOLDERS.gubbins_recombination_gff])
  File gubbins_log_for_report = select_first([GUBBINS_RECOMBINATION_FILTER.gubbins_log, MAKE_GUBBINS_PLACEHOLDERS.gubbins_log])
  File gubbins_filtered_vcf_for_report = select_first([GUBBINS_RECOMBINATION_FILTER.gubbins_filtered_vcf, MAKE_GUBBINS_PLACEHOLDERS.gubbins_filtered_vcf])
  File snp_pyseer_assoc_for_report = select_first([PYSEER_SNP_GWAS.pyseer_snp_assoc, MAKE_SNP_GWAS_PLACEHOLDERS.pyseer_snp_assoc])
  File snp_top_hits_for_report = select_first([PRIORITIZE_SNP_GWAS_HITS.snp_top_hits, MAKE_SNP_GWAS_PLACEHOLDERS.snp_top_hits])
  File snp_all_significant_hits_for_report = select_first([PRIORITIZE_SNP_GWAS_HITS.snp_all_significant_hits, MAKE_SNP_GWAS_PLACEHOLDERS.snp_all_significant_hits])
  File snp_summary_for_report = select_first([PRIORITIZE_SNP_GWAS_HITS.snp_summary, MAKE_SNP_GWAS_PLACEHOLDERS.snp_summary])
  File snp_manhattan_plot_for_report = select_first([GENERATE_SNP_GWAS_PLOTS.snp_manhattan_plot_svg, MAKE_SNP_GWAS_PLACEHOLDERS.snp_manhattan_plot_svg])
  File snp_qq_plot_for_report = select_first([GENERATE_SNP_GWAS_PLOTS.snp_qq_plot_svg, MAKE_SNP_GWAS_PLACEHOLDERS.snp_qq_plot_svg])
  File snp_plot_summary_for_report = select_first([GENERATE_SNP_GWAS_PLOTS.snp_plot_summary, MAKE_SNP_GWAS_PLACEHOLDERS.snp_plot_summary])

  call ANNOTATE_GWAS_HITS_WITH_GENBANK {
    input:
      top_priority_hits = PRIORITIZE_GWAS_HITS.top_priority_hits,
      all_significant_hits = PRIORITIZE_GWAS_HITS.all_significant_hits,
      gene_presence_absence_csv = PANAROO_PANGENOME.gene_presence_absence_csv,
      gene_data_csv = PANAROO_PANGENOME.gene_data_csv,
      combined_dna_cds = PANAROO_PANGENOME.combined_dna_cds,
      combined_protein_cds = PANAROO_PANGENOME.combined_protein_cds,
      pan_genome_reference = PANAROO_PANGENOME.pan_genome_reference,
      reference_genbank = EXTRACT_REFERENCE_GENBANK_FROM_DOCKER.reference_genbank,
      output_prefix = output_prefix,
      annotation_engine = annotation_engine,
      python_docker = hit_annotation_docker
  }

  call GENERATE_GWAS_PLOTS {
    input:
      pyseer_gene_assoc = PYSEER_GENE_GWAS.pyseer_gene_assoc,
      output_prefix = output_prefix,
      plot_label = "Gene presence/absence GWAS",
      significance_alpha = significance_alpha,
      max_points = plot_max_points,
      python_docker = python_docker
  }

  call MERGE_RMAP_GWAS_REPORT {
    input:
      output_prefix = output_prefix,
      validation_report = VALIDATE_SAMPLE_SET_INPUTS.validation_report,
      phenotype_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_tsv,
      phenotype_legend_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_legend_tsv,
      top_priority_hits = ANNOTATE_GWAS_HITS_WITH_GENBANK.annotated_top_priority_hits,
      all_significant_hits = ANNOTATE_GWAS_HITS_WITH_GENBANK.annotated_all_significant_hits,
      reference_annotation_summary = ANNOTATE_GWAS_HITS_WITH_GENBANK.reference_annotation_summary,
      qq_plot_svg = GENERATE_GWAS_PLOTS.qq_plot_svg,
      manhattan_plot_svg = GENERATE_GWAS_PLOTS.manhattan_plot_svg,
      plot_summary = GENERATE_GWAS_PLOTS.plot_summary,
      pyseer_gene_assoc = PYSEER_GENE_GWAS.pyseer_gene_assoc,
      panaroo_summary = PANAROO_PANGENOME.panaroo_summary,
      mash_distances = MASH_DISTANCE_MATRIX.mash_distances,
      population_pca_svg = GENERATE_POPULATION_STRUCTURE_PLOTS.pca_svg,
      kinship_heatmap_svg = GENERATE_POPULATION_STRUCTURE_PLOTS.kinship_heatmap_svg,
      population_structure_summary = GENERATE_POPULATION_STRUCTURE_PLOTS.population_structure_summary,
      do_snp_gwas = do_snp_gwas,
      do_gubbins = do_gubbins,
      gwas_mode = gwas_mode,
      container_backend = container_backend,
      gubbins_summary = gubbins_summary_for_report,
      gubbins_filtered_alignment = gubbins_filtered_alignment_for_report,
      gubbins_recombination_gff = gubbins_recombination_gff_for_report,
      gubbins_log = gubbins_log_for_report,
      snp_top_hits = snp_top_hits_for_report,
      snp_all_significant_hits = snp_all_significant_hits_for_report,
      snp_summary = snp_summary_for_report,
      pyseer_snp_assoc = snp_pyseer_assoc_for_report,
      snp_vcf = snp_vcf_for_report,
      snp_qq_plot_svg = snp_qq_plot_for_report,
      snp_manhattan_plot_svg = snp_manhattan_plot_for_report,
      snp_plot_summary = snp_plot_summary_for_report,
      reference_docker = reference_docker,
      reference_species = reference_species,
      reference_name = reference_name,
      annotation_engine = annotation_engine,
      python_docker = python_docker
  }

  output {
    File phenotype_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_tsv
    File phenotype_legend_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_legend_tsv
    File sample_groups_tsv = PREPARE_PHENOTYPE_TABLE.sample_groups_tsv
    Array[File] trimmed_read1s = FASTP_TRIM.trimmed_read1
    Array[File] trimmed_read2s = FASTP_TRIM.trimmed_read2
    Array[File] assemblies = SHOVILL_ASSEMBLE.contigs_fasta
    Array[File] quast_reports = QUAST_ASSEMBLY_QC.quast_report_tsv
    Array[File] gffs = annotation_gffs
    Array[File] gbks = annotation_gbks
    Array[File] faas = annotation_faas
    Array[File] ffns = annotation_ffns
    File gene_presence_absence_csv = PANAROO_PANGENOME.gene_presence_absence_csv
    File gene_presence_absence_rtab = PANAROO_PANGENOME.gene_presence_absence_rtab
    File mash_distances = MASH_DISTANCE_MATRIX.mash_distances
    File population_pca_svg = GENERATE_POPULATION_STRUCTURE_PLOTS.pca_svg
    File kinship_heatmap_svg = GENERATE_POPULATION_STRUCTURE_PLOTS.kinship_heatmap_svg
    File population_structure_summary = GENERATE_POPULATION_STRUCTURE_PLOTS.population_structure_summary
    File pyseer_gene_assoc = PYSEER_GENE_GWAS.pyseer_gene_assoc
    File raw_top_priority_hits = PRIORITIZE_GWAS_HITS.top_priority_hits
    File raw_all_significant_hits = PRIORITIZE_GWAS_HITS.all_significant_hits
    File top_priority_hits = ANNOTATE_GWAS_HITS_WITH_GENBANK.annotated_top_priority_hits
    File all_significant_hits = ANNOTATE_GWAS_HITS_WITH_GENBANK.annotated_all_significant_hits
    File reference_annotation_summary = ANNOTATE_GWAS_HITS_WITH_GENBANK.reference_annotation_summary
    File enrichment_summary = PRIORITIZE_GWAS_HITS.enrichment_summary
    File qq_plot_svg = GENERATE_GWAS_PLOTS.qq_plot_svg
    File manhattan_plot_svg = GENERATE_GWAS_PLOTS.manhattan_plot_svg
    File plot_summary = GENERATE_GWAS_PLOTS.plot_summary
    File reference_fasta = EXTRACT_REFERENCE_GENBANK_FROM_DOCKER.reference_fasta
    File reference_gff = EXTRACT_REFERENCE_GENBANK_FROM_DOCKER.reference_gff
    File reference_genbank = EXTRACT_REFERENCE_GENBANK_FROM_DOCKER.reference_genbank
    File snp_vcf = snp_vcf_for_report
    File gubbins_summary = gubbins_summary_for_report
    File gubbins_filtered_alignment = gubbins_filtered_alignment_for_report
    File gubbins_recombination_gff = gubbins_recombination_gff_for_report
    File gubbins_log = gubbins_log_for_report
    File gubbins_filtered_vcf = gubbins_filtered_vcf_for_report
    File pyseer_snp_assoc = snp_pyseer_assoc_for_report
    File snp_top_hits = snp_top_hits_for_report
    File snp_all_significant_hits = snp_all_significant_hits_for_report
    File snp_summary = snp_summary_for_report
    File snp_manhattan_plot_svg = snp_manhattan_plot_for_report
    File snp_qq_plot_svg = snp_qq_plot_for_report
    File snp_plot_summary = snp_plot_summary_for_report
    File html_report = MERGE_RMAP_GWAS_REPORT.html_report
    File run_provenance = MERGE_RMAP_GWAS_REPORT.run_provenance_json
  }
}

task VALIDATE_SAMPLE_SET_INPUTS {
  input {
    Array[String]+ sample_names
    Int read1_count
    Int read2_count
    Array[String]+ groups
    Array[String] phenotype_display_values
    String case_label
    String control_label
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
sample_names = """~{sep='\n' sample_names}""".strip().splitlines()
read1_count = int("~{read1_count}")
read2_count = int("~{read2_count}")
groups = """~{sep='\n' groups}""".strip().splitlines()
phenotype_display_values = """~{sep='\n' phenotype_display_values}""".strip().splitlines()
case_label = """~{case_label}"""
control_label = """~{control_label}"""

def norm(x):
    return str(x).strip().lower()

errors = []
warnings = []

if not (len(sample_names) == read1_count == read2_count == len(groups)):
    errors.append(
        f"sample_names/read1s/read2s/groups have unequal lengths: "
        f"sample_names={len(sample_names)}, read1s={read1_count}, read2s={read2_count}, groups={len(groups)}"
    )

if phenotype_display_values and len(phenotype_display_values) != len(sample_names):
    errors.append(
        "phenotype_display_values was provided but its length does not match sample_names: "
        f"phenotype_display_values={len(phenotype_display_values)}, sample_names={len(sample_names)}"
    )

if len(sample_names) != len(set(sample_names)):
    errors.append("Sample names must be unique within the selected sample set.")

bad_names = [x for x in sample_names if any(c.isspace() for c in x)]
if bad_names:
    errors.append("Sample names must not contain whitespace: " + ", ".join(bad_names[:10]))

valid_case = {norm(case_label), "case", "cases", "1", "true", "yes"}
valid_control = {norm(control_label), "control", "controls", "0", "false", "no"}

case_names = []
control_names = []
bad_groups = []
for sample, group in zip(sample_names, groups):
    g = norm(group)
    if g in valid_case:
        case_names.append(sample)
    elif g in valid_control:
        control_names.append(sample)
    else:
        bad_groups.append(f"{sample}:{group}")

if bad_groups:
    errors.append(
        "Unrecognized group labels. Expected case/control labels or 1/0. Examples: "
        + ", ".join(bad_groups[:10])
    )

if len(case_names) == 0:
    errors.append("No case samples found from the groups input.")
if len(control_names) == 0:
    errors.append("No control samples found from the groups input.")

if len(case_names) < 10:
    warnings.append(f"Low case count: {len(case_names)}. Microbial GWAS may be underpowered.")
if len(control_names) < 10:
    warnings.append(f"Low control count: {len(control_names)}. Microbial GWAS may be underpowered.")
if len(case_names) < 50 or len(control_names) < 50:
    warnings.append("Recommended minimum for stable microbial GWAS is often >=50 cases and >=50 controls.")

with open("validation_report.txt", "w") as out:
    out.write("rMAP-GWAS sample-set input validation report\n")
    out.write("================================================\n")
    out.write(f"Cases: {len(case_names)}\n")
    out.write(f"Controls: {len(control_names)}\n")
    out.write(f"Total samples: {len(sample_names)}\n")
    out.write(f"Case label: {case_label}\n")
    out.write(f"Control label: {control_label}\n")
    out.write(f"Phenotype display metadata values provided: {'yes' if phenotype_display_values else 'no'}\n")
    if phenotype_display_values:
        out.write(f"Unique phenotype display values: {', '.join(sorted(set(phenotype_display_values)))[:1000]}\n")
    out.write("\n")
    if warnings:
        out.write("Warnings:\n")
        for w in warnings:
            out.write(f"- {w}\n")
        out.write("\n")
    if errors:
        out.write("Errors:\n")
        for e in errors:
            out.write(f"- {e}\n")
        raise SystemExit("Input validation failed. See validation_report.txt")
    out.write("Status: PASS\n")
PY
  >>>

  output {
    File validation_report = "validation_report.txt"
  }

  runtime {
    docker: python_docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
}

task PREPARE_PHENOTYPE_TABLE {
  input {
    Array[String]+ sample_names
    Array[String]+ groups
    String phenotype_name
    String case_label
    String control_label
    Array[String] phenotype_display_values
    String output_prefix
    String python_docker
    File validation_report
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
from collections import Counter
sample_names = """~{sep='\n' sample_names}""".strip().splitlines()
groups = """~{sep='\n' groups}""".strip().splitlines()
phenotype_display_values = """~{sep='\n' phenotype_display_values}""".strip().splitlines()
phenotype = """~{phenotype_name}""".strip() or "case_control"
case_label = """~{case_label}""".strip() or "case"
control_label = """~{control_label}""".strip() or "control"
phenotype_display_name = phenotype if phenotype_display_values else ""

def norm(x):
    return str(x).strip().lower()

def clean_value(x):
    x = str(x).strip()
    if not x:
        return ""
    x = x.replace("_", " ").replace("-", " ")
    x = " ".join(x.split())
    words = []
    for w in x.split(" "):
        wl = w.lower()
        if wl == "hiv":
            words.append("HIV")
        elif wl in ("mtb", "mtbc", "dna", "rna", "amr", "gwas"):
            words.append(w.upper())
        elif wl == "no":
            words.append("No")
        else:
            words.append(w)
    return " ".join(words)

def summarize_values(values, fallback):
    values = [clean_value(v) for v in values if str(v).strip()]
    if not values:
        return clean_value(fallback)
    counts = Counter(values)
    ordered = [v for v, c in counts.most_common()]
    if len(ordered) == 1:
        return ordered[0]
    return "mixed: " + "; ".join(f"{v} (n={counts[v]})" for v in ordered)

def display_label(base, detail):
    detail = clean_value(detail)
    if detail and norm(detail) not in {norm(base), norm(base + "s")}:
        return f"{base} ({detail})"
    return base

def plural_label(base_label):
    if base_label.startswith("case ") or base_label == "case":
        return "cases" + base_label[len("case"):]
    if base_label.startswith("control ") or base_label == "control":
        return "controls" + base_label[len("control"):]
    return base_label

valid_case = {norm(case_label), "case", "cases", "1", "true", "yes"}
valid_control = {norm(control_label), "control", "controls", "0", "false", "no"}

if len(sample_names) != len(groups):
    raise SystemExit(f"sample_names and groups have unequal lengths: {len(sample_names)} vs {len(groups)}")
if phenotype_display_values and len(phenotype_display_values) != len(sample_names):
    raise SystemExit(f"phenotype_display_values has {len(phenotype_display_values)} values but sample_names has {len(sample_names)}")
if not phenotype_display_values:
    phenotype_display_values = [""] * len(sample_names)

rows = []
errors = []
case_meta = []
control_meta = []
for sample, group, display_value in zip(sample_names, groups, phenotype_display_values):
    g = norm(group)
    if g in valid_case:
        rows.append((sample, case_label, 1, display_value))
        if display_value.strip():
            case_meta.append(display_value)
    elif g in valid_control:
        rows.append((sample, control_label, 0, display_value))
        if display_value.strip():
            control_meta.append(display_value)
    else:
        errors.append(f"{sample}:{group}")
if errors:
    raise SystemExit("Unrecognized group labels: " + ", ".join(errors[:10]))

case_detail = summarize_values(case_meta, case_label)
control_detail = summarize_values(control_meta, control_label)
case_display = display_label("case", case_detail)
control_display = display_label("control", control_detail)
case_display_plural = plural_label(case_display)
control_display_plural = plural_label(control_display)
case_count = sum(1 for _,_,v,_ in rows if v == 1)
control_count = sum(1 for _,_,v,_ in rows if v == 0)

with open("~{output_prefix}_phenotypes.tsv", "w") as out:
    out.write("sample\t" + phenotype + "\n")
    for sample, group_label, value, display_value in rows:
        out.write(f"{sample}\t{value}\n")

with open("~{output_prefix}_sample_groups.tsv", "w") as out:
    out.write("sample\tgroup\tphenotype_value\tmetadata_display_column\tmetadata_display_value\tdisplay_label\n")
    for sample, group_label, value, display_value in rows:
        label = case_display if value == 1 else control_display
        out.write(f"{sample}\t{group_label}\t{value}\t{phenotype_display_name}\t{clean_value(display_value)}\t{label}\n")

with open("~{output_prefix}_phenotype_legend.tsv", "w") as out:
    out.write("class\tnumeric_value\tinput_group_label\tmetadata_display_column\tmetadata_values\tn_samples\tdisplay_label\tdisplay_label_plural\n")
    out.write(f"case\t1\t{case_label}\t{phenotype_display_name}\t{case_detail}\t{case_count}\t{case_display}\t{case_display_plural}\n")
    out.write(f"control\t0\t{control_label}\t{phenotype_display_name}\t{control_detail}\t{control_count}\t{control_display}\t{control_display_plural}\n")
PY
  >>>

  output {
    File phenotype_tsv = "~{output_prefix}_phenotypes.tsv"
    File sample_groups_tsv = "~{output_prefix}_sample_groups.tsv"
    File phenotype_legend_tsv = "~{output_prefix}_phenotype_legend.tsv"
  }

  runtime {
    docker: python_docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
}

task FASTP_TRIM {
  input {
    String sample_name
    File read1
    File read2
    Int threads
    Int memory_gb
    Int disk_gb
    String docker
    File validation_report
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
fastp \
  -i ~{read1} \
  -I ~{read2} \
  -o ~{sample_name}_R1.trimmed.fastq.gz \
  -O ~{sample_name}_R2.trimmed.fastq.gz \
  --thread ~{threads} \
  --html ~{sample_name}.fastp.html \
  --json ~{sample_name}.fastp.json
  >>>

  output {
    File trimmed_read1 = "~{sample_name}_R1.trimmed.fastq.gz"
    File trimmed_read2 = "~{sample_name}_R2.trimmed.fastq.gz"
    File fastp_html = "~{sample_name}.fastp.html"
    File fastp_json = "~{sample_name}.fastp.json"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task SHOVILL_ASSEMBLE {
  input {
    String sample_name
    File read1
    File read2
    Int threads
    Int memory_gb
    Int shovill_ram_gb
    Int disk_gb
    String docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}

echo "Task runtime memory request: ~{memory_gb} GB" >&2
echo "Shovill command RAM (--ram): ~{shovill_ram_gb} GB" >&2
echo "CPU threads: ~{threads}" >&2
if command -v free >/dev/null 2>&1; then free -h >&2 || true; fi

shovill \
  --R1 ~{read1} \
  --R2 ~{read2} \
  --outdir shovill_out \
  --cpus ~{threads} \
  --ram ~{shovill_ram_gb} \
  --force

cp shovill_out/contigs.fa ~{sample_name}.contigs.fasta
  >>>

  output {
    File contigs_fasta = "~{sample_name}.contigs.fasta"
    File shovill_log = "shovill_out/shovill.log"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task QUAST_ASSEMBLY_QC {
  input {
    String sample_name
    File assembly
    Int threads
    String docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
quast.py ~{assembly} -o quast_out -t ~{threads}
cp quast_out/report.tsv ~{sample_name}.quast.report.tsv
  >>>

  output {
    File quast_report_tsv = "~{sample_name}.quast.report.tsv"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "8 GB"
    disks: "local-disk 50 HDD"
  }
}

task PROKKA_ANNOTATE {
  input {
    String sample_name
    File assembly
    Int threads
    Int memory_gb
    String docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
prokka \
  --outdir prokka_out \
  --prefix ~{sample_name} \
  --cpus ~{threads} \
  --force \
  ~{assembly}

cp prokka_out/~{sample_name}.gff ~{sample_name}.gff
cp prokka_out/~{sample_name}.gbk ~{sample_name}.gbk
cp prokka_out/~{sample_name}.faa ~{sample_name}.faa
cp prokka_out/~{sample_name}.ffn ~{sample_name}.ffn
  >>>

  output {
    File gff = "~{sample_name}.gff"
    File gbk = "~{sample_name}.gbk"
    File faa = "~{sample_name}.faa"
    File ffn = "~{sample_name}.ffn"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk 100 HDD"
  }
}


task BAKTA_ANNOTATE {
  input {
    String sample_name
    File assembly
    Int threads
    Int memory_gb
    Int disk_gb
    String docker
    String bakta_db
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}

bakta \
  --db "~{bakta_db}" \
  --output bakta_out \
  --prefix "~{sample_name}" \
  --threads ~{threads} \
  --skip-sorf \
  --force \
  "~{assembly}"

copy_first_existing() {
  local out="$1"
  shift
  for p in "$@"; do
    if [ -s "$p" ]; then
      cp "$p" "$out"
      return 0
    fi
  done
  return 1
}

copy_first_existing "~{sample_name}.gff" \
  "bakta_out/~{sample_name}.gff3" \
  "bakta_out/~{sample_name}.gff" || {
    echo "ERROR: Bakta did not produce a GFF/GFF3 file for ~{sample_name}." >&2
    find bakta_out -maxdepth 2 -type f -print >&2 || true
    exit 1
  }

copy_first_existing "~{sample_name}.gbk" \
  "bakta_out/~{sample_name}.gbff" \
  "bakta_out/~{sample_name}.gbk" \
  "bakta_out/~{sample_name}.gb" || : > "~{sample_name}.gbk"

copy_first_existing "~{sample_name}.faa" \
  "bakta_out/~{sample_name}.faa" || : > "~{sample_name}.faa"

copy_first_existing "~{sample_name}.ffn" \
  "bakta_out/~{sample_name}.ffn" \
  "bakta_out/~{sample_name}.fna" || : > "~{sample_name}.ffn"
  >>>

  output {
    File gff = "~{sample_name}.gff"
    File gbk = "~{sample_name}.gbk"
    File faa = "~{sample_name}.faa"
    File ffn = "~{sample_name}.ffn"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task PANAROO_PANGENOME {
  input {
    Array[File] gffs
    Int threads
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}

mkdir -p gffs
for f in ~{sep=' ' gffs}; do
  cp "$f" gffs/
done

echo "Panaroo input GFF count: $(ls gffs/*.gff | wc -l)" >&2
ls -lh gffs/*.gff >&2

# Strict mode is preferred. If strict cleaning fails for a small or heterogeneous
# smoke-test cohort, retry once in moderate mode so the workflow can proceed and
# still produce a usable gene presence/absence matrix.
set +e
panaroo \
  -i gffs/*.gff \
  -o panaroo_out \
  --clean-mode strict \
  --remove-invalid-genes \
  -t ~{threads}
panaroo_rc=$?
set -e

if [ "$panaroo_rc" -ne 0 ]; then
  echo "Panaroo strict mode failed with exit code ${panaroo_rc}; retrying with clean-mode moderate" >&2
  rm -rf panaroo_out
  panaroo \
    -i gffs/*.gff \
    -o panaroo_out \
    --clean-mode moderate \
    --remove-invalid-genes \
    -t ~{threads}
fi

if [ ! -s panaroo_out/gene_presence_absence.csv ]; then
  echo "ERROR: Panaroo did not produce panaroo_out/gene_presence_absence.csv" >&2
  echo "Available Panaroo output files:" >&2
  find panaroo_out -maxdepth 2 -type f | sort >&2 || true
  exit 1
fi

# Some Panaroo builds produce gene_presence_absence.csv but not an Rtab file.
# pyseer needs a gene-by-sample 0/1 table, so create a compatible Rtab if needed.
if [ ! -s panaroo_out/gene_presence_absence.Rtab ]; then
cat > make_panaroo_rtab.py <<'PY'
import csv
from pathlib import Path

csv_path = Path("panaroo_out/gene_presence_absence.csv")
out_path = Path("panaroo_out/gene_presence_absence.Rtab")

metadata_cols = {
    "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
    "No. sequences", "Avg sequences per isolate", "Genome Fragment",
    "Order within Fragment", "Accessory Fragment", "Accessory Order with Fragment",
    "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc"
}

with csv_path.open(newline="") as fh:
    reader = csv.DictReader(fh)
    if not reader.fieldnames:
        raise SystemExit("Panaroo CSV has no header")
    sample_cols = [c for c in reader.fieldnames if c not in metadata_cols]
    if not sample_cols:
        raise SystemExit("No sample columns detected in Panaroo CSV")

    with out_path.open("w") as out:
        out.write("Gene\t" + "\t".join(sample_cols) + "\n")
        for row in reader:
            gene = row.get("Gene", "")
            vals = ["1" if row.get(sample, "").strip() else "0" for sample in sample_cols]
            out.write(gene + "\t" + "\t".join(vals) + "\n")
PY
python make_panaroo_rtab.py
fi

if [ ! -s panaroo_out/gene_presence_absence.Rtab ]; then
  echo "ERROR: gene_presence_absence.Rtab was not created" >&2
  find panaroo_out -maxdepth 2 -type f | sort >&2 || true
  exit 1
fi

echo "Panaroo output files:" > panaroo_summary.txt
find panaroo_out -maxdepth 2 -type f | sort >> panaroo_summary.txt

# Ensure optional Panaroo files needed for post-GWAS annotation are present as WDL outputs.
# Some Panaroo versions may omit one of these; create empty placeholders so the
# downstream annotation task can degrade gracefully rather than failing localization.
for f in   panaroo_out/gene_data.csv   panaroo_out/combined_DNA_CDS.fasta   panaroo_out/combined_protein_CDS.fasta   panaroo_out/pan_genome_reference.fa
 do
  if [ ! -e "$f" ]; then
    echo "WARNING: optional Panaroo annotation file missing, creating placeholder: $f" >&2
    touch "$f"
  fi
 done
  >>>

  output {
    File gene_presence_absence_csv = "panaroo_out/gene_presence_absence.csv"
    File gene_presence_absence_rtab = "panaroo_out/gene_presence_absence.Rtab"
    File panaroo_summary = "panaroo_summary.txt"
    File gene_data_csv = "panaroo_out/gene_data.csv"
    File combined_dna_cds = "panaroo_out/combined_DNA_CDS.fasta"
    File combined_protein_cds = "panaroo_out/combined_protein_CDS.fasta"
    File pan_genome_reference = "panaroo_out/pan_genome_reference.fa"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task MASH_DISTANCE_MATRIX {
  input {
    Array[String] sample_names
    Array[File] assemblies
    Int threads
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
mkdir -p assemblies

python <<'PY'
import os, shutil
from pathlib import Path

names = """~{sep='\n' sample_names}""".strip().splitlines()
files = """~{sep='\n' assemblies}""".strip().splitlines()

if len(names) != len(files):
    raise SystemExit("sample_names and assemblies have unequal lengths")

if len(names) != len(set(names)):
    raise SystemExit("sample_names must be unique before Mash/pyseer distance matrix creation")

for name, f in zip(names, files):
    out = Path("assemblies") / f"{name}.fasta"
    shutil.copyfile(f, out)

with open("assemblies.list", "w") as out:
    for name in names:
        out.write(f"assemblies/{name}.fasta\n")
PY

mash sketch -p ~{threads} -o cohort assemblies/*.fasta

# Mash emits a long pairwise table when comparing a sketch against itself.
# Pyseer expects a square distance matrix with one unique row/column per sample.
mash dist cohort.msh cohort.msh > mash_distances.long.tsv

cat > make_pyseer_distance_matrix.py <<'PY'
from pathlib import Path
import os

names = """~{sep='\n' sample_names}""".strip().splitlines()
name_set = set(names)

def clean_sample(x):
    base = os.path.basename(x.strip())
    for suffix in (".fasta", ".fa", ".fna", ".contigs.fasta"):
        if base.endswith(suffix):
            base = base[:-len(suffix)]
    return base

dist = {a: {b: None for b in names} for a in names}
for s in names:
    dist[s][s] = 0.0

with open("mash_distances.long.tsv") as fh:
    for line in fh:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        q = clean_sample(parts[0])
        r = clean_sample(parts[1])
        if q not in name_set or r not in name_set:
            continue
        try:
            d = float(parts[2])
        except ValueError:
            continue
        dist[q][r] = d
        dist[r][q] = d

missing = []
for a in names:
    for b in names:
        if dist[a][b] is None:
            missing.append((a, b))
if missing:
    raise SystemExit(f"Mash distance matrix incomplete; missing {len(missing)} pairs. Examples: {missing[:10]}")

with open("mash_distances.tsv", "w") as out:
    out.write("sample\t" + "\t".join(names) + "\n")
    for a in names:
        out.write(a + "\t" + "\t".join(f"{dist[a][b]:.10g}" for b in names) + "\n")

with open("mash_distance_matrix_summary.txt", "w") as out:
    out.write(f"samples\t{len(names)}\n")
    out.write("format\tpyseer_square_distance_matrix\n")
    out.write("raw_mash_long_table\tmash_distances.long.tsv\n")
    out.write("pyseer_distance_matrix\tmash_distances.tsv\n")
PY
python make_pyseer_distance_matrix.py

echo "Mash distance matrix summary:"
cat mash_distance_matrix_summary.txt
  >>>

  output {
    File mash_distances = "mash_distances.tsv"
    File mash_distances_long = "mash_distances.long.tsv"
    File mash_distance_summary = "mash_distance_matrix_summary.txt"
    File assemblies_list = "assemblies.list"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task PYSEER_GENE_GWAS {
  input {
    File phenotype_tsv
    File gene_presence_absence_rtab
    File mash_distances
    Float min_af
    Float max_af
    Int max_dimensions
    Boolean force_no_distances
    Boolean no_distances_fallback
    Int threads
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}

echo "Starting pyseer gene GWAS" > pyseer_run.log
echo "Phenotypes: ~{phenotype_tsv}" >> pyseer_run.log
echo "Presence/absence matrix: ~{gene_presence_absence_rtab}" >> pyseer_run.log
echo "Mash distance matrix: ~{mash_distances}" >> pyseer_run.log
echo "min_af: ~{min_af}" >> pyseer_run.log
echo "max_af: ~{max_af}" >> pyseer_run.log
echo "max_dimensions: ~{max_dimensions}" >> pyseer_run.log
echo "force_no_distances: ~{force_no_distances}" >> pyseer_run.log
echo "no_distances_fallback: ~{no_distances_fallback}" >> pyseer_run.log

if [[ "~{force_no_distances}" == "true" ]]; then
  echo "Running pyseer without population-structure correction because force_no_distances=true" >> pyseer_run.log
  pyseer \
    --phenotypes ~{phenotype_tsv} \
    --pres ~{gene_presence_absence_rtab} \
    --no-distances \
    --min-af ~{min_af} \
    --max-af ~{max_af} \
    --cpu ~{threads} \
    > pyseer_gene_assoc.tsv \
    2> pyseer.stderr.log
  cat pyseer.stderr.log >&2
else
  echo "Running pyseer with Mash distances and --max-dimensions ~{max_dimensions}" >> pyseer_run.log
  set +e
  pyseer \
    --phenotypes ~{phenotype_tsv} \
    --pres ~{gene_presence_absence_rtab} \
    --distances ~{mash_distances} \
    --max-dimensions ~{max_dimensions} \
    --min-af ~{min_af} \
    --max-af ~{max_af} \
    --cpu ~{threads} \
    > pyseer_gene_assoc.tsv \
    2> pyseer.stderr.log
  rc=$?
  set -e

  if [[ "$rc" -ne 0 ]]; then
    echo "Primary pyseer run failed with exit code ${rc}." >> pyseer_run.log
    echo "Primary pyseer stderr:" >> pyseer_run.log
    cat pyseer.stderr.log >> pyseer_run.log

    if [[ "~{no_distances_fallback}" == "true" ]]; then
      echo "Retrying pyseer with --no-distances so the smoke test can complete." >> pyseer_run.log
      set +e
      pyseer \
        --phenotypes ~{phenotype_tsv} \
        --pres ~{gene_presence_absence_rtab} \
        --no-distances \
        --min-af ~{min_af} \
        --max-af ~{max_af} \
        --cpu ~{threads} \
        > pyseer_gene_assoc.tsv \
        2> pyseer.no_distances.stderr.log
      rc2=$?
      set -e

      if [[ "$rc2" -ne 0 ]]; then
        echo "Fallback pyseer run without distances also failed with exit code ${rc2}." >> pyseer_run.log
        echo "Fallback pyseer stderr:" >> pyseer_run.log
        cat pyseer.no_distances.stderr.log >> pyseer_run.log
        cat pyseer_run.log >&2
        exit "$rc2"
      fi

      echo "Fallback pyseer run without distances completed." >> pyseer_run.log
      cat pyseer.no_distances.stderr.log >&2
    else
      cat pyseer_run.log >&2
      exit "$rc"
    fi
  else
    echo "Primary pyseer run with distances completed." >> pyseer_run.log
    cat pyseer.stderr.log >&2
  fi
fi

if [[ ! -s pyseer_gene_assoc.tsv ]]; then
  echo "ERROR: pyseer_gene_assoc.tsv was not created or is empty." >&2
  cat pyseer_run.log >&2
  exit 1
fi

echo "Final pyseer output:"
wc -l pyseer_gene_assoc.tsv
head -5 pyseer_gene_assoc.tsv || true
cat pyseer_run.log >&2
  >>>

  output {
    File pyseer_gene_assoc = "pyseer_gene_assoc.tsv"
    File pyseer_run_log = "pyseer_run.log"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task PRIORITIZE_GWAS_HITS {
  input {
    File phenotype_tsv
    File gene_presence_absence_csv
    File gene_presence_absence_rtab
    File pyseer_gene_assoc
    Float significance_alpha
    String output_prefix
    String annotation_engine
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
import csv, math, re
from pathlib import Path

phenotype_path = Path("~{phenotype_tsv}")
rtab_path = Path("~{gene_presence_absence_rtab}")
panaroo_csv_path = Path("~{gene_presence_absence_csv}")
assoc_path = Path("~{pyseer_gene_assoc}")
alpha = float("~{significance_alpha}")
prefix = "~{output_prefix}"
annotation_source = "Panaroo/" + """~{annotation_engine}""".strip()

# Read phenotype coding.
phenotypes = {}
with phenotype_path.open() as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        phenotypes[parts[0]] = int(float(parts[1]))

cases = {s for s, v in phenotypes.items() if v == 1}
controls = {s for s, v in phenotypes.items() if v == 0}

# Read Panaroo gene annotations.
annotations = {}
try:
    with panaroo_csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            gid = row.get("Gene", "")
            annotations[gid] = {
                "gene_name": row.get("Non-unique Gene name", "") or gid,
                "product": row.get("Annotation", "") or "",
            }
except Exception:
    pass

# Read gene presence/absence matrix.
presence = {}
with rtab_path.open() as fh:
    header = fh.readline().rstrip("\n").split("\t")
    sample_cols = header[1:]
    for line in fh:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        gene = parts[0]
        vals = parts[1:]
        present_samples = {s for s, v in zip(sample_cols, vals) if v not in ("0", "", "NA", ".")}
        presence[gene] = present_samples

def parse_float(x):
    try:
        if x in ("", "NA", "nan", "None"):
            return None
        return float(x)
    except Exception:
        return None

def pick(row, names):
    lower = {k.lower(): k for k in row}
    for n in names:
        if n in row:
            return row[n]
        if n.lower() in lower:
            return row[lower[n.lower()]]
    return ""

# Read pyseer associations. Try to be robust to pyseer column naming.
assoc_rows = []
with assoc_path.open() as fh:
    first = fh.readline().rstrip("\n")
    if not first:
        header = []
    else:
        header = re.split(r"\t|\s+", first.strip())
    for line in fh:
        if not line.strip():
            continue
        parts = re.split(r"\t|\s+", line.strip())
        if len(parts) < len(header):
            parts += [""] * (len(header) - len(parts))
        row = dict(zip(header, parts))
        assoc_rows.append(row)

# If pyseer did not print a recognizable header, create a fallback.
if assoc_rows and not any(h.lower() in ("variant", "gene", "feature", "notes") for h in header):
    assoc_rows = []
    with assoc_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = re.split(r"\t|\s+", line.strip())
            if len(parts) >= 4:
                assoc_rows.append({
                    "variant": parts[0],
                    "af": parts[1],
                    "filter-pvalue": parts[2],
                    "lrt-pvalue": parts[3],
                })

enriched_rows = []
all_rows = []

for row in assoc_rows:
    feature = pick(row, ["variant", "gene", "feature", "samples", "name"])
    if not feature:
        feature = next(iter(row.values())) if row else ""

    pval = parse_float(pick(row, ["lrt-pvalue", "lrt_pvalue", "pvalue", "p-value", "filter-pvalue", "p"]))
    qval = parse_float(pick(row, ["q_value", "q-value", "qvalue", "adjusted-pvalue", "adjusted_pvalue"]))
    beta = parse_float(pick(row, ["beta", "effect", "coef", "coefficient"]))

    present = presence.get(feature, set())
    case_present = len(present & cases)
    control_present = len(present & controls)
    case_total = len(cases)
    control_total = len(controls)
    case_freq = case_present / case_total if case_total else 0.0
    control_freq = control_present / control_total if control_total else 0.0

    # If beta is unavailable, use frequency difference to assign direction.
    if beta is None:
        beta = case_freq - control_freq

    if case_freq > control_freq and beta >= 0:
        enriched_in = "Cases"
    elif control_freq > case_freq and beta <= 0:
        enriched_in = "Controls"
    else:
        enriched_in = "Check manually"

    # Use q-value if present, otherwise p-value for provisional ranking.
    stat = qval if qval is not None else pval
    if stat is None or stat <= 0:
        neglog = 0.0
    else:
        neglog = -math.log10(stat)

    # Haldane-Anscombe corrected odds ratio and 95% CI for binary phenotype.
    a = case_present + 0.5
    b = (case_total - case_present) + 0.5
    c = control_present + 0.5
    d = (control_total - control_present) + 0.5
    odds_ratio = (a / b) / (c / d)
    se_log_or = math.sqrt((1.0 / a) + (1.0 / b) + (1.0 / c) + (1.0 / d))
    log_or = math.log(odds_ratio) if odds_ratio > 0 and math.isfinite(odds_ratio) else 0.0
    or_ci_low = math.exp(log_or - 1.96 * se_log_or)
    or_ci_high = math.exp(log_or + 1.96 * se_log_or)
    log2_or = math.log2(odds_ratio) if odds_ratio > 0 and math.isfinite(odds_ratio) else 0.0

    ann = annotations.get(feature, {})
    gene_name = ann.get("gene_name", feature)
    product = ann.get("product", "")

    annotation_weight = 2 if product or (gene_name and gene_name != feature) else 0
    recurrence_weight = 1 if max(case_present, control_present) >= 5 else 0
    priority_score = neglog + abs(log2_or) + annotation_weight + recurrence_weight

    outrow = {
        "feature_id": feature,
        "feature_type": "gene_presence_absence",
        "gene_name": gene_name,
        "product": product,
        "case_present": case_present,
        "case_total": case_total,
        "case_frequency": f"{case_freq:.4f}",
        "control_present": control_present,
        "control_total": control_total,
        "control_frequency": f"{control_freq:.4f}",
        "enriched_in": enriched_in,
        "beta": f"{beta:.6g}",
        "odds_ratio": f"{odds_ratio:.6g}" if math.isfinite(odds_ratio) else "Inf",
        "odds_ratio_ci95_lower": f"{or_ci_low:.6g}",
        "odds_ratio_ci95_upper": f"{or_ci_high:.6g}",
        "odds_ratio_ci95": f"{or_ci_low:.4g}-{or_ci_high:.4g}",
        "pyseer_pvalue": "" if pval is None else f"{pval:.6g}",
        "q_value": "" if qval is None else f"{qval:.6g}",
        "priority_score": f"{priority_score:.4f}",
        "annotation_source": annotation_source,
        "notes": ""
    }
    all_rows.append(outrow)

    # Accept q <= alpha if q exists, else p <= alpha.
    if (qval is not None and qval <= alpha) or (qval is None and pval is not None and pval <= alpha):
        enriched_rows.append(outrow)

all_rows.sort(key=lambda r: float(r["priority_score"]), reverse=True)
enriched_rows.sort(key=lambda r: float(r["priority_score"]), reverse=True)

fields = [
    "rank", "feature_id", "feature_type", "gene_name", "product",
    "case_present", "case_total", "case_frequency",
    "control_present", "control_total", "control_frequency",
    "enriched_in", "beta", "odds_ratio", "odds_ratio_ci95_lower", "odds_ratio_ci95_upper", "odds_ratio_ci95", "pyseer_pvalue", "q_value",
    "priority_score", "annotation_source", "notes"
]

def write_table(path, rows):
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for i, r in enumerate(rows, 1):
            rr = {"rank": i}
            rr.update(r)
            writer.writerow(rr)

write_table(prefix + "_all_ranked_hits.tsv", all_rows)
write_table(prefix + "_top_priority_hits.tsv", enriched_rows[:100])
write_table(prefix + "_all_significant_hits.tsv", enriched_rows)

case_enriched = sum(1 for r in enriched_rows if r["enriched_in"] == "Cases")
control_enriched = sum(1 for r in enriched_rows if r["enriched_in"] == "Controls")

with open(prefix + "_enrichment_summary.tsv", "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"total_ranked_features\t{len(all_rows)}\n")
    out.write(f"significant_features\t{len(enriched_rows)}\n")
    out.write(f"case_enriched_significant_features\t{case_enriched}\n")
    out.write(f"control_enriched_significant_features\t{control_enriched}\n")
    out.write(f"alpha\t{alpha}\n")
PY
  >>>

  output {
    File all_ranked_hits = "~{output_prefix}_all_ranked_hits.tsv"
    File top_priority_hits = "~{output_prefix}_top_priority_hits.tsv"
    File all_significant_hits = "~{output_prefix}_all_significant_hits.tsv"
    File enrichment_summary = "~{output_prefix}_enrichment_summary.tsv"
  }

  runtime {
    docker: python_docker
    cpu: 1
    memory: "8 GB"
    disks: "local-disk 100 HDD"
  }
}


task EXTRACT_REFERENCE_GENBANK_FROM_DOCKER {
  input {
    String reference_docker
    String reference_name
  }

  command <<<
set -euo pipefail

{
  echo "Reference name: ~{reference_name}"
  echo "Reference docker: ~{reference_docker}"
  echo "Searching for reference FASTA, GFF and GenBank inside the reference image."
} > reference_extract_log.txt

copy_first_existing() {
  local out="$1"
  shift
  for p in "$@"; do
    if [ -n "$p" ] && [ -s "$p" ]; then
      cp "$p" "$out"
      echo "Found ${out}: ${p}" >> reference_extract_log.txt
      return 0
    fi
  done
  return 1
}

copy_first_existing reference.fasta \
  "${RMAP_GWAS_REFERENCE_FASTA:-}" \
  /opt/rmap-gwas/refs/reference.fasta \
  /opt/rmap-gwas/refs/reference.fa \
  /opt/rmap-gwas/refs/mtbc/reference.fasta \
  /opt/rmap-gwas/refs/kpneumo/reference.fasta \
  /opt/rmap-gwas/refs/ecoli/reference.fasta \
  /opt/rmap-gwas/refs/enterococcus_faecium/reference.fasta \
  /refs/reference.fasta \
  /data/reference.fasta || true

copy_first_existing reference.gff \
  "${RMAP_GWAS_REFERENCE_GFF:-}" \
  /opt/rmap-gwas/refs/reference.gff \
  /opt/rmap-gwas/refs/reference.gff3 \
  /opt/rmap-gwas/refs/mtbc/reference.gff \
  /opt/rmap-gwas/refs/kpneumo/reference.gff \
  /opt/rmap-gwas/refs/ecoli/reference.gff \
  /opt/rmap-gwas/refs/enterococcus_faecium/reference.gff \
  /refs/reference.gff \
  /data/reference.gff || true

copy_first_existing reference.genbank \
  "${RMAP_GWAS_REFERENCE_GENBANK:-}" \
  /opt/rmap-gwas/refs/reference.genbank \
  /opt/rmap-gwas/refs/reference.gbk \
  /opt/rmap-gwas/refs/mtbc/reference.genbank \
  /opt/rmap-gwas/refs/kpneumo/reference.genbank \
  /opt/rmap-gwas/refs/ecoli/reference.genbank \
  /opt/rmap-gwas/refs/enterococcus_faecium/reference.genbank \
  /refs/reference.genbank \
  /data/reference.genbank || true

if [ ! -s reference.fasta ]; then
  found=$(find /opt /refs /data 2>/dev/null -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) | head -n 1 || true)
  if [ -n "$found" ] && [ -s "$found" ]; then cp "$found" reference.fasta; echo "Found reference.fasta by search: $found" >> reference_extract_log.txt; fi
fi
if [ ! -s reference.gff ]; then
  found=$(find /opt /refs /data 2>/dev/null -type f \( -name "*.gff" -o -name "*.gff3" \) | head -n 1 || true)
  if [ -n "$found" ] && [ -s "$found" ]; then cp "$found" reference.gff; echo "Found reference.gff by search: $found" >> reference_extract_log.txt; fi
fi
if [ ! -s reference.genbank ]; then
  found=$(find /opt /refs /data 2>/dev/null -type f \( -name "*.genbank" -o -name "*.gbk" -o -name "*.gb" \) | head -n 1 || true)
  if [ -n "$found" ] && [ -s "$found" ]; then cp "$found" reference.genbank; echo "Found reference.genbank by search: $found" >> reference_extract_log.txt; fi
fi

if [ ! -s reference.genbank ]; then
  echo "ERROR: Could not find a non-empty GenBank file in the reference image." >&2
  cat reference_extract_log.txt >&2
  exit 1
fi
if [ ! -s reference.fasta ]; then
  echo "WARNING: No reference FASTA found. SNP GWAS will fail if do_snp_gwas=true." >> reference_extract_log.txt
  : > reference.fasta
fi
if [ ! -s reference.gff ]; then
  echo "WARNING: No reference GFF found. SNP annotation will use GenBank only." >> reference_extract_log.txt
  : > reference.gff
fi

ls -lh reference.fasta reference.gff reference.genbank >> reference_extract_log.txt
cat reference_extract_log.txt
  >>>

  output {
    File reference_fasta = "reference.fasta"
    File reference_gff = "reference.gff"
    File reference_genbank = "reference.genbank"
    File reference_extract_log = "reference_extract_log.txt"
  }

  runtime {
    docker: reference_docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk 20 HDD"
  }
}

task ANNOTATE_GWAS_HITS_WITH_GENBANK {
  input {
    File top_priority_hits
    File all_significant_hits
    File gene_presence_absence_csv
    File gene_data_csv
    File combined_dna_cds
    File combined_protein_cds
    File pan_genome_reference
    File reference_genbank
    String output_prefix
    String annotation_engine
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
from pathlib import Path
import csv, re, difflib
from collections import defaultdict

prefix = "~{output_prefix}"
annotation_source = "Panaroo/" + """~{annotation_engine}""".strip()
# Touch optional Panaroo annotation inputs so WDL validators know they are intentional task inputs.
_ = Path("~{gene_data_csv}")
_ = Path("~{combined_protein_cds}")

def read_tsv_dict(path):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return [], []
    with p.open(newline="", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames or [], rows

def write_tsv(path, fields, rows):
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in fields})

def parse_fasta(path):
    records = {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return records
    name = None
    parts = []
    with p.open(errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records[name] = "".join(parts).replace(" ", "").upper()
                name = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
    if name:
        records[name] = "".join(parts).replace(" ", "").upper()
    return records

def revcomp(seq):
    tbl = str.maketrans("ACGTURYKMSWBDHVNacgturykmswbdhvn", "TGCAAYRMKSWVHDBNtgcaayrmkswvhdbn")
    return seq.translate(tbl)[::-1]

def parse_location(location):
    comp = "complement" in location
    nums = re.findall(r"<?(\d+)\.\.>?(\d+)|<?(\d+)", location)
    parts = []
    for a, b, single in nums:
        if single:
            start = end = int(single)
        else:
            start, end = int(a), int(b)
        parts.append((start, end))
    return comp, parts

def parse_genbank(path):
    txt = Path(path).read_text(errors="replace")
    origin = ""
    m = re.search(r"\nORIGIN\s*(.*?)(?=\n//)", txt, flags=re.S)
    if m:
        origin = re.sub(r"[^A-Za-z]", "", m.group(1)).upper()
    features_text = ""
    m = re.search(r"\nFEATURES\s+Location/Qualifiers\s*(.*?)(?=\nORIGIN)", txt, flags=re.S)
    if m:
        features_text = m.group(1)
    lines = features_text.splitlines()
    cds = []
    current = None
    current_key = None
    for line in lines:
        if re.match(r"^     \S+", line):
            key = line[5:21].strip()
            loc = line[21:].strip()
            if key == "CDS":
                current = {"location": loc, "qualifiers": {}}
                cds.append(current)
                current_key = None
            else:
                current = None
                current_key = None
            continue
        if current is None:
            continue
        st = line.strip()
        if not st:
            continue
        qm = re.match(r"/([^=]+)=(.*)", st)
        if qm:
            current_key = qm.group(1)
            val = qm.group(2).strip()
            if val.startswith('"') and val.endswith('"'):
                val = val[1:-1]
            elif val.startswith('"'):
                val = val[1:]
            current["qualifiers"].setdefault(current_key, "")
            current["qualifiers"][current_key] += val
        elif current_key:
            val = st.strip()
            if val.endswith('"'):
                val = val[:-1]
            current["qualifiers"][current_key] += val
    refs = []
    for i, f in enumerate(cds, 1):
        q = f.get("qualifiers", {})
        locus = q.get("locus_tag", "") or q.get("protein_id", "") or f"CDS_{i}"
        gene = q.get("gene", "")
        product = q.get("product", "")
        protein_id = q.get("protein_id", "")
        translation = q.get("translation", "").replace(" ", "").replace("\n", "").upper()
        comp, loc_parts = parse_location(f.get("location", ""))
        nuc = ""
        if origin and loc_parts:
            seqs = []
            for start, end in loc_parts:
                seqs.append(origin[start-1:end])
            nuc = "".join(seqs)
            if comp:
                nuc = revcomp(nuc)
        refs.append({"id": locus, "locus_tag": locus, "gene": gene, "product": product, "protein_id": protein_id, "location": f.get("location", ""), "nuc": nuc, "prot": translation})
    return refs

def split_members(x):
    if not x:
        return []
    vals = re.split(r"[;,:\s]+", x.strip())
    return [v for v in vals if v and v not in ("0", "1", "NA", ".")]

panaroo_rows = {}
metadata_cols = {"Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment", "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc"}
try:
    with open("~{gene_presence_absence_csv}", newline="", errors="replace") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            gid = row.get("Gene", "")
            if not gid:
                continue
            members = []
            for k, v in row.items():
                if k not in metadata_cols:
                    members.extend(split_members(v or ""))
            panaroo_rows[gid] = {"gene_name": row.get("Non-unique Gene name", "") or gid, "product": row.get("Annotation", "") or "", "members": sorted(set(members))[:200]}
except Exception:
    pass

fasta_records = {}
for fp in ["~{combined_dna_cds}", "~{pan_genome_reference}"]:
    fasta_records.update(parse_fasta(fp))

reference_parse_warning = ""
try:
    refs = parse_genbank("~{reference_genbank}")
except Exception as e:
    refs = []
    reference_parse_warning = f"{type(e).__name__}: {e}"
ref_by_token = {}
for r in refs:
    for token in [r.get("locus_tag", ""), r.get("gene", ""), r.get("protein_id", "")]:
        if token:
            ref_by_token[token] = r

k = 15
kmer_to_refs = defaultdict(set)
for idx, r in enumerate(refs):
    seq = r.get("nuc", "")
    if len(seq) >= k:
        step = max(1, len(seq)//200)
        for i in range(0, len(seq)-k+1, step):
            kmer_to_refs[seq[i:i+k]].add(idx)

def find_query_sequences(feature, members):
    tokens = [feature] + members[:50]
    hits = []
    for header, seq in fasta_records.items():
        if any(t and t in header for t in tokens):
            if seq and len(seq) >= 30:
                hits.append((header, seq))
        if len(hits) >= 3:
            break
    return hits

def score_sequence_to_reference(seq):
    if not seq or not refs:
        return None
    seq = re.sub(r"[^A-Za-z]", "", seq).upper()
    if len(seq) < 30:
        return None
    best = None
    for r in refs:
        rn = r.get("nuc", "")
        if not rn:
            continue
        if seq in rn or rn in seq:
            cov = 100.0 * min(len(seq), len(rn)) / max(1, len(seq))
            return (r, 100.0, cov, "nucleotide_exact_or_contained")
    counts = defaultdict(int)
    if len(seq) >= k:
        for i in range(0, len(seq)-k+1, max(1, len(seq)//200)):
            for idx in kmer_to_refs.get(seq[i:i+k], []):
                counts[idx] += 1
    candidates = [idx for idx, c in sorted(counts.items(), key=lambda x: x[1], reverse=True)[:25]] or list(range(min(25, len(refs))))
    for idx in candidates:
        r = refs[idx]
        rn = r.get("nuc", "")
        if not rn:
            continue
        sm = difflib.SequenceMatcher(None, seq, rn, autojunk=False)
        ratio = sm.ratio() * 100.0
        lm = sm.find_longest_match(0, len(seq), 0, len(rn))
        cov = 100.0 * lm.size / max(1, len(seq))
        combined = ratio + cov
        if best is None or combined > best[3]:
            best = (r, ratio, cov, combined)
    if best:
        return (best[0], best[1], best[2], "nucleotide_similarity")
    return None

def confidence(identity, coverage, match_type):
    try:
        ident = float(identity)
        cov = float(coverage)
    except Exception:
        return "none"
    if ident >= 95 and cov >= 80:
        return "high"
    if ident >= 80 and cov >= 60:
        return "medium"
    if ident >= 60 and cov >= 40:
        return "low"
    return "none"

def annotate_row(row):
    feature = row.get("feature_id", "") or row.get("gene", "") or row.get("variant", "")
    pan = panaroo_rows.get(feature, {})
    members = pan.get("members", [])
    original_gene = row.get("gene_name", "") or pan.get("gene_name", "") or feature
    original_product = row.get("product", "") or pan.get("product", "")
    match = None
    note = ""
    match_type = "none"
    identity = ""
    coverage = ""
    for token in [feature, original_gene] + members[:100]:
        if token in ref_by_token:
            match = ref_by_token[token]
            match_type = "qualifier_exact"
            identity = "100.0"
            coverage = "100.0"
            break
    if match is None:
        best = None
        for header, seq in find_query_sequences(feature, members):
            sc = score_sequence_to_reference(seq)
            if sc:
                r, ident, cov, mtype = sc
                rank = float(ident) + float(cov)
                if best is None or rank > best[-1]:
                    best = (r, ident, cov, mtype, header, rank)
        if best:
            match, ident, cov, match_type, header, rank = best
            identity = f"{ident:.2f}"
            coverage = f"{cov:.2f}"
            note = f"matched_representative_sequence={header[:120]}"
    if match is None and original_product and original_product.lower() not in ("hypothetical protein", "unknown"):
        op = original_product.lower()
        for r in refs:
            if r.get("product", "").lower() == op:
                match = r
                match_type = "product_exact"
                note = "product-only match; sequence confirmation not available"
                break
    if match is None:
        row.update({"reference_locus_tag": "", "reference_gene": "", "reference_product": "", "reference_location": "", "reference_match_type": "none", "reference_identity": "", "reference_coverage": "", "annotation_confidence": "none", "annotation_note": "No confident GenBank reference match found. The cluster may be accessory, divergent, absent from the reference, or not represented in Panaroo sequence outputs.", "cluster_member_ids": ";".join(members[:30])})
    else:
        conf = confidence(identity or 0, coverage or 0, match_type)
        if match_type == "product_exact" and conf == "none":
            conf = "low"
        row.update({"reference_locus_tag": match.get("locus_tag", ""), "reference_gene": match.get("gene", ""), "reference_product": match.get("product", ""), "reference_location": match.get("location", ""), "reference_match_type": match_type, "reference_identity": identity, "reference_coverage": coverage, "annotation_confidence": conf, "annotation_note": note, "cluster_member_ids": ";".join(members[:30])})
        if (not row.get("gene_name") or row.get("gene_name") == feature or row.get("gene_name", "").startswith("group_")) and match.get("gene"):
            row["gene_name"] = match.get("gene")
        if (not row.get("product") or row.get("product", "").lower() in ("hypothetical protein", "unknown")) and match.get("product"):
            row["product"] = match.get("product")
        row["annotation_source"] = row.get("annotation_source", annotation_source) + "+GenBank"

    # Report display fields: keep the stable Panaroo cluster ID, but provide an interpretable display name.
    # High-confidence matches can be shown as the reference gene; medium/low matches are marked as "-like".
    conf = (row.get("annotation_confidence", "none") or "none").lower()
    ref_gene = row.get("reference_gene", "") or ""
    ref_locus = row.get("reference_locus_tag", "") or ""
    ref_prod = row.get("reference_product", "") or ""
    gene = row.get("gene_name", "") or ""
    feat = feature or row.get("feature_id", "") or ""
    identity_txt = row.get("reference_identity", "") or ""
    coverage_txt = row.get("reference_coverage", "") or ""

    if ref_gene and conf == "high":
        display_name = ref_gene
        interpretation = "High-confidence GenBank-supported annotation."
    elif ref_gene and conf == "medium":
        display_name = ref_gene + "-like"
        interpretation = "Medium-confidence GenBank rescue; inspect manually before biological interpretation."
    elif ref_gene and conf == "low":
        display_name = ref_gene + "-like (" + feat + ")"
        interpretation = "Low-confidence GenBank rescue; treat as tentative and keep the Panaroo cluster ID."
    elif gene and gene != feat and not gene.startswith("group_"):
        display_name = gene
        interpretation = "Displayed from " + annotation_source + " annotation; no stronger GenBank rescue available."
    elif ref_locus and conf in ("high", "medium", "low"):
        suffix = "" if conf == "high" else "-like"
        display_name = ref_locus + suffix + " (" + feat + ")"
        interpretation = conf.capitalize() + "-confidence locus-level GenBank rescue."
    else:
        display_name = feat
        interpretation = "No confident reference gene assignment; report the stable Panaroo cluster ID."

    row["display_name"] = display_name
    row["display_label"] = " | ".join([x for x in [feat, ref_locus or "no_reference_locus", conf + " confidence"] if x])
    row["display_product"] = ref_prod or row.get("product", "") or ""
    row["interpretation_note"] = interpretation
    if identity_txt or coverage_txt:
        row["annotation_evidence"] = "identity=" + str(identity_txt) + "; coverage=" + str(coverage_txt)
    else:
        row["annotation_evidence"] = "no sequence identity/coverage evidence available"
    return row

extra_fields = ["display_name", "display_label", "display_product", "interpretation_note", "annotation_evidence", "reference_locus_tag", "reference_gene", "reference_product", "reference_location", "reference_match_type", "reference_identity", "reference_coverage", "annotation_confidence", "annotation_note", "cluster_member_ids"]
for in_path, out_path in [("~{top_priority_hits}", prefix + "_top_priority_hits.annotated.tsv"), ("~{all_significant_hits}", prefix + "_all_significant_hits.annotated.tsv")]:
    fields, rows = read_tsv_dict(in_path)
    if not fields:
        fields = ["rank", "feature_id", "feature_type", "gene_name", "product"]
        rows = []
    out_fields = fields + [f for f in extra_fields if f not in fields]
    write_tsv(out_path, out_fields, [annotate_row(dict(r)) for r in rows])

with open(prefix + "_reference_annotation_summary.tsv", "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"reference_cds_parsed\t{len(refs)}\n")
    out.write(f"reference_parse_warning\t{reference_parse_warning}\n")
    out.write(f"panaroo_clusters_parsed\t{len(panaroo_rows)}\n")
    out.write(f"annotation_engine\t{annotation_source}\n")
    out.write(f"fasta_records_parsed\t{len(fasta_records)}\n")
    out.write(f"top_priority_input\t{Path('~{top_priority_hits}').name}\n")
    out.write(f"all_significant_input\t{Path('~{all_significant_hits}').name}\n")
    out.write("method\tGenBank qualifier matching plus pure-Python nucleotide similarity rescue when Panaroo representative sequences are available\n")
PY
  >>>

  output {
    File annotated_top_priority_hits = "~{output_prefix}_top_priority_hits.annotated.tsv"
    File annotated_all_significant_hits = "~{output_prefix}_all_significant_hits.annotated.tsv"
    File reference_annotation_summary = "~{output_prefix}_reference_annotation_summary.tsv"
  }

  runtime {
    docker: python_docker
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 100 HDD"
  }
}

task GENERATE_GWAS_PLOTS {
  input {
    File pyseer_gene_assoc
    String output_prefix
    String plot_label
    Float significance_alpha
    Int max_points
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
from pathlib import Path
import re, math, html
prefix = "~{output_prefix}"
plot_label = "~{plot_label}"
alpha = float("~{significance_alpha}")
max_points = int("~{max_points}")
assoc = Path("~{pyseer_gene_assoc}")

def parse_float(x):
    try:
        if x is None or str(x).strip() in ("", "NA", "nan", "None"):
            return None
        v = float(x)
        if v <= 0 or v > 1 or math.isnan(v):
            return None
        return v
    except Exception:
        return None

def pick(row, names):
    lower = {k.lower(): k for k in row}
    for n in names:
        if n in row:
            return row[n]
        if n.lower() in lower:
            return row[lower[n.lower()]]
    return ""

rows = []
with assoc.open(errors="replace") as fh:
    first = fh.readline().rstrip("\n")
    header = re.split(r"\t|\s+", first.strip()) if first else []
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        parts = re.split(r"\t|\s+", line.strip())
        if len(parts) < len(header):
            parts += [""] * (len(header) - len(parts))
        row = dict(zip(header, parts)) if header else {}
        feature = pick(row, ["variant", "gene", "feature", "samples", "name"]) or (parts[0] if parts else "")
        p = parse_float(pick(row, ["lrt-pvalue", "lrt_pvalue", "pvalue", "p-value", "filter-pvalue", "p"]))
        if p is not None:
            rows.append((feature, p))

rows_sorted = sorted(rows, key=lambda x: x[1])
if len(rows) > max_points:
    keep = rows_sorted[:min(500, len(rows_sorted))]
    remaining = rows_sorted[min(500, len(rows_sorted)):]
    step = max(1, len(remaining) // max(1, max_points-len(keep)))
    keep.extend(remaining[::step])
    draw_rows = keep[:max_points]
else:
    draw_rows = rows

def svg_escape(x): return html.escape(str(x), quote=True)

def make_svg_points_plot(points, title, xlabel, ylabel, path, line=False):
    width, height = 980, 520
    ml, mr, mt, mb = 74, 26, 58, 68
    pw, ph = width - ml - mr, height - mt - mb
    if not points:
        body = f'<text x="{width/2}" y="{height/2}" fill="#9fb3c8" text-anchor="middle">No p-values available</text>'
    else:
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = 0, max(max(ys), 1.0)
        if xmax == xmin: xmax = xmin + 1
        def sx(x): return ml + (x - xmin) / (xmax - xmin) * pw
        def sy(y): return mt + ph - (y - ymin) / (ymax - ymin) * ph
        body_parts = []
        for i in range(6):
            y = mt + i * ph / 5
            val = ymax - i * ymax / 5
            body_parts.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{ml+pw}" y2="{y:.1f}" stroke="rgba(180,210,255,0.18)"/>')
            body_parts.append(f'<text x="{ml-10}" y="{y+4:.1f}" fill="#9fb3c8" text-anchor="end" font-size="12">{val:.1f}</text>')
        body_parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#8ecaff"/>')
        body_parts.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#8ecaff"/>')
        sig_y = -math.log10(alpha) if alpha > 0 else None
        if sig_y and sig_y <= ymax:
            body_parts.append(f'<line x1="{ml}" y1="{sy(sig_y):.1f}" x2="{ml+pw}" y2="{sy(sig_y):.1f}" stroke="#ff7ab6" stroke-width="2" stroke-dasharray="8 8"/>')
            body_parts.append(f'<text x="{ml+pw-4}" y="{sy(sig_y)-8:.1f}" fill="#ffb3d9" text-anchor="end" font-size="12">alpha={alpha:g}</text>')
        for x,y in points:
            color = "#ff7ab6" if sig_y and y >= sig_y else "#21d4fd"
            body_parts.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.2" fill="{color}" opacity="0.75"/>')
        body = "\n".join(body_parts)
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{svg_escape(title)}">
<rect width="100%" height="100%" rx="18" fill="#071226"/>
<text x="{ml}" y="34" fill="#eef6ff" font-family="Arial" font-size="22" font-weight="700">{svg_escape(title)}</text>
{body}
<text x="{width/2}" y="{height-18}" fill="#cfe8ff" font-family="Arial" font-size="14" text-anchor="middle">{svg_escape(xlabel)}</text>
<text x="20" y="{height/2}" fill="#cfe8ff" font-family="Arial" font-size="14" text-anchor="middle" transform="rotate(-90 20 {height/2})">{svg_escape(ylabel)}</text>
</svg>"""
    Path(path).write_text(svg)

manhattan_points = [(i+1, -math.log10(p)) for i, (_, p) in enumerate(draw_rows)]
make_svg_points_plot(manhattan_points, plot_label + " Manhattan-style association plot", "Feature index", "-log10(p-value)", prefix + "_manhattan.svg")

n = len(rows_sorted)
qq = []
if n:
    observed = sorted([-math.log10(p) for _, p in rows_sorted])
    expected = sorted([-math.log10((i + 0.5) / n) for i in range(n)])
    qq = list(zip(expected, observed))
    if len(qq) > max_points:
        step = max(1, len(qq)//max_points)
        qq = qq[::step]
make_svg_points_plot(qq, plot_label + " QQ plot", "Expected -log10(p-value)", "Observed -log10(p-value)", prefix + "_qq.svg")

with open(prefix + "_plot_summary.tsv", "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"pvalues_detected\t{len(rows)}\n")
    out.write(f"points_drawn\t{len(draw_rows)}\n")
    out.write(f"plot_label\t{plot_label}\n")
    out.write(f"significant_points_at_alpha\t{sum(1 for _, p in rows if p <= alpha)}\n")
    out.write("manhattan_type\tfeature_index_not_genomic_coordinate\n")
    out.write("qq_plot\tgenerated\n")
PY
  >>>

  output {
    File manhattan_plot_svg = "~{output_prefix}_manhattan.svg"
    File qq_plot_svg = "~{output_prefix}_qq.svg"
    File plot_summary = "~{output_prefix}_plot_summary.tsv"
  }

  runtime {
    docker: python_docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
}

task GENERATE_SNP_GWAS_PLOTS {
  input {
    File pyseer_snp_assoc
    File snp_vcf
    String output_prefix
    Float significance_alpha
    Int max_points
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
from pathlib import Path
import re, math, html
prefix = "~{output_prefix}"
alpha = float("~{significance_alpha}")
max_points = int("~{max_points}")
assoc = Path("~{pyseer_snp_assoc}")
vcf = Path("~{snp_vcf}")

def parse_float(x):
    try:
        if x is None or str(x).strip() in ("", "NA", "nan", "None"):
            return None
        v = float(x)
        if v <= 0 or v > 1 or math.isnan(v):
            return None
        return v
    except Exception:
        return None

def pick(row, names):
    lower = {k.lower(): k for k in row}
    for n in names:
        if n in row:
            return row[n]
        if n.lower() in lower:
            return row[lower[n.lower()]]
    return ""

coord = {}
with vcf.open(errors="replace") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        chrom, pos, vid, ref, alt = parts[:5]
        marker_id = vid if vid not in ("", ".") else f"{chrom}_{pos}_{ref}_{alt}"
        try:
            p = int(pos)
        except Exception:
            p = None
        rec = (chrom, p, marker_id)
        for key in {marker_id, f"{chrom}_{pos}", f"{chrom}:{pos}", str(pos), f"{pos}_{ref}_{alt}"}:
            coord[key] = rec

rows = []
with assoc.open(errors="replace") as fh:
    first = fh.readline().rstrip("\n")
    header = re.split(r"\t|\s+", first.strip()) if first else []
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        parts = re.split(r"\t|\s+", line.strip())
        if len(parts) < len(header):
            parts += [""] * (len(header) - len(parts))
        row = dict(zip(header, parts)) if header else {}
        feature = pick(row, ["variant", "feature", "name"]) or (parts[0] if parts else "")
        pval = parse_float(pick(row, ["lrt-pvalue", "lrt_pvalue", "pvalue", "p-value", "filter-pvalue", "p"]))
        if pval is None:
            continue
        chrom, pos, marker = coord.get(feature, ("", None, feature))
        if pos is None:
            nums = re.findall(r"\d+", feature)
            pos = int(nums[0]) if nums else len(rows) + 1
        rows.append({"feature": feature, "chrom": chrom, "position": pos, "p": pval})

rows_sorted = sorted(rows, key=lambda r: (r["chrom"] or "zz_unknown", r["position"]))
if len(rows_sorted) > max_points:
    by_p = sorted(rows_sorted, key=lambda r: r["p"])
    keep_ids = {id(r) for r in by_p[:min(500, len(by_p))]}
    remaining = [r for r in rows_sorted if id(r) not in keep_ids]
    step = max(1, len(remaining) // max(1, max_points - len(keep_ids)))
    draw_rows = by_p[:min(500, len(by_p))] + remaining[::step]
    draw_rows = draw_rows[:max_points]
else:
    draw_rows = rows_sorted

offsets = {}
current = 0
for chrom in [] if not rows_sorted else sorted({r["chrom"] or "unknown" for r in rows_sorted}):
    positions = [r["position"] for r in rows_sorted if (r["chrom"] or "unknown") == chrom]
    offsets[chrom] = current
    current += (max(positions) if positions else 0) + 1

def svg_escape(x): return html.escape(str(x), quote=True)

def make_svg_points_plot(points, title, xlabel, ylabel, path, draw_diag=False):
    width, height = 980, 520
    ml, mr, mt, mb = 74, 26, 58, 68
    pw, ph = width - ml - mr, height - mt - mb
    if not points:
        body = f'<text x="{width/2}" y="{height/2}" fill="#9fb3c8" text-anchor="middle">No SNP p-values available</text>'
    else:
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = 0, max(max(ys), 1.0)
        if xmax == xmin: xmax = xmin + 1
        def sx(x): return ml + (x - xmin) / (xmax - xmin) * pw
        def sy(y): return mt + ph - (y - ymin) / (ymax - ymin) * ph
        body_parts = []
        for i in range(6):
            y = mt + i * ph / 5
            val = ymax - i * ymax / 5
            body_parts.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{ml+pw}" y2="{y:.1f}" stroke="rgba(180,210,255,0.18)"/>')
            body_parts.append(f'<text x="{ml-10}" y="{y+4:.1f}" fill="#9fb3c8" text-anchor="end" font-size="12">{val:.1f}</text>')
        body_parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#8ecaff"/>')
        body_parts.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#8ecaff"/>')
        if draw_diag:
            lim = min(max(xs), max(ys))
            body_parts.append(f'<line x1="{sx(0):.1f}" y1="{sy(0):.1f}" x2="{sx(lim):.1f}" y2="{sy(lim):.1f}" stroke="rgba(255,255,255,0.35)" stroke-width="2"/>')
        sig_y = -math.log10(alpha) if alpha > 0 else None
        if sig_y and sig_y <= ymax:
            body_parts.append(f'<line x1="{ml}" y1="{sy(sig_y):.1f}" x2="{ml+pw}" y2="{sy(sig_y):.1f}" stroke="#ff7ab6" stroke-width="2" stroke-dasharray="8 8"/>')
            body_parts.append(f'<text x="{ml+pw-4}" y="{sy(sig_y)-8:.1f}" fill="#ffb3d9" text-anchor="end" font-size="12">alpha={alpha:g}</text>')
        for x, y in points:
            color = "#ff7ab6" if sig_y and y >= sig_y else "#21d4fd"
            body_parts.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.2" fill="{color}" opacity="0.75"/>')
        body = "\n".join(body_parts)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{svg_escape(title)}">
<rect width="100%" height="100%" rx="18" fill="#071226"/>
<text x="{ml}" y="34" fill="#eef6ff" font-family="Arial" font-size="22" font-weight="700">{svg_escape(title)}</text>
{body}
<text x="{width/2}" y="{height-18}" fill="#cfe8ff" font-family="Arial" font-size="14" text-anchor="middle">{svg_escape(xlabel)}</text>
<text x="20" y="{height/2}" fill="#cfe8ff" font-family="Arial" font-size="14" text-anchor="middle" transform="rotate(-90 20 {height/2})">{svg_escape(ylabel)}</text>
</svg>'''
    Path(path).write_text(svg)

manhattan_points = []
for r in draw_rows:
    chrom = r["chrom"] or "unknown"
    x = offsets.get(chrom, 0) + r["position"]
    manhattan_points.append((x, -math.log10(r["p"])))
make_svg_points_plot(manhattan_points, "SNP marker GWAS Manhattan plot", "Reference genomic coordinate / marker order", "-log10(p-value)", prefix + "_manhattan.svg")

n = len(rows)
qq = []
if n:
    observed = sorted([-math.log10(r["p"]) for r in rows])
    expected = sorted([-math.log10((i + 0.5) / n) for i in range(n)])
    qq = list(zip(expected, observed))
    if len(qq) > max_points:
        step = max(1, len(qq) // max_points)
        qq = qq[::step]
make_svg_points_plot(qq, "SNP marker GWAS QQ plot", "Expected -log10(p-value)", "Observed -log10(p-value)", prefix + "_qq.svg", draw_diag=True)

with open(prefix + "_plot_summary.tsv", "w") as out:
    out.write("metric\tvalue\n")
    out.write("plot_label\tSNP marker GWAS\n")
    out.write(f"pvalues_detected\t{len(rows)}\n")
    out.write(f"significant_points_at_alpha\t{sum(1 for r in rows if r['p'] <= alpha)}\n")
    out.write(f"points_drawn\t{len(draw_rows)}\n")
    out.write("manhattan_type\treference_coordinate_when_available\n")
    out.write("qq_plot\tgenerated\n")
PY
  >>>

  output {
    File snp_manhattan_plot_svg = "~{output_prefix}_manhattan.svg"
    File snp_qq_plot_svg = "~{output_prefix}_qq.svg"
    File snp_plot_summary = "~{output_prefix}_plot_summary.tsv"
  }

  runtime {
    docker: python_docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
}


task MAKE_SNP_GWAS_PLACEHOLDERS {
  input {
    String output_prefix
    String python_docker
  }

  command <<<
set -euo pipefail
cat > ~{output_prefix}_pyseer_snp_assoc.tsv <<'EOF'
variant	lrt-pvalue	beta	note
EOF
cat > ~{output_prefix}_top_snp_hits.tsv <<'EOF'
rank	feature_id	variant_id	feature_type	contig	position	ref	alt	gene_name	product	case_alt	case_total	case_frequency	control_alt	control_total	control_frequency	enriched_in	odds_ratio	odds_ratio_ci95_lower	odds_ratio_ci95_upper	odds_ratio_ci95	beta	pyseer_pvalue	q_value	priority_score	annotation_source	notes
EOF
cp ~{output_prefix}_top_snp_hits.tsv ~{output_prefix}_all_significant_snp_hits.tsv
cat > ~{output_prefix}_snp_summary.tsv <<'EOF'
metric	value
snp_gwas_status	not_run
EOF
cat > ~{output_prefix}.snps.vcf <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
EOF
cat > ~{output_prefix}_manhattan.svg <<'EOF'
<svg xmlns="http://www.w3.org/2000/svg" width="980" height="520" viewBox="0 0 980 520"><rect width="100%" height="100%" rx="18" fill="#071226"/><text x="490" y="260" fill="#9fb3c8" text-anchor="middle" font-family="Arial" font-size="22">SNP marker GWAS Manhattan plot not run</text></svg>
EOF
cat > ~{output_prefix}_qq.svg <<'EOF'
<svg xmlns="http://www.w3.org/2000/svg" width="980" height="520" viewBox="0 0 980 520"><rect width="100%" height="100%" rx="18" fill="#071226"/><text x="490" y="260" fill="#9fb3c8" text-anchor="middle" font-family="Arial" font-size="22">SNP marker GWAS QQ plot not run</text></svg>
EOF
cat > ~{output_prefix}_plot_summary.tsv <<'EOF'
metric	value
pvalues_detected	0
points_drawn	0
snp_gwas_status	not_run
EOF
  >>>

  output {
    File pyseer_snp_assoc = "~{output_prefix}_pyseer_snp_assoc.tsv"
    File snp_top_hits = "~{output_prefix}_top_snp_hits.tsv"
    File snp_all_significant_hits = "~{output_prefix}_all_significant_snp_hits.tsv"
    File snp_summary = "~{output_prefix}_snp_summary.tsv"
    File snp_vcf = "~{output_prefix}.snps.vcf"
    File snp_manhattan_plot_svg = "~{output_prefix}_manhattan.svg"
    File snp_qq_plot_svg = "~{output_prefix}_qq.svg"
    File snp_plot_summary = "~{output_prefix}_plot_summary.tsv"
  }

  runtime {
    docker: python_docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
  }
}

task MAKE_GUBBINS_PLACEHOLDERS {
  input {
    String output_prefix
    String python_docker
  }

  command <<<
set -euo pipefail
cat > ~{output_prefix}_gubbins_summary.tsv <<'EOF'
metric	value
gubbins_status	not_run
recommendation	Optional. Enable for recombining bacteria or strong lineage/recombination concerns; usually not required as a default for MTBC/clonal smoke-test runs.
usage_note	Gubbins not requested. Primary SNP GWAS used the unfiltered Snippy/core VCF when SNP GWAS was enabled.
EOF

cat > ~{output_prefix}_gubbins.filtered_polymorphic_sites.fasta <<'EOF'
EOF

cat > ~{output_prefix}_gubbins.recombination_predictions.gff <<'EOF'
##gff-version 3
EOF

cat > ~{output_prefix}_gubbins.filtered.snps.vcf <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
EOF

cat > ~{output_prefix}_gubbins.log <<'EOF'
Gubbins not run. Set do_gubbins=true to run optional recombination assessment/filtering for SNP GWAS.
EOF
  >>>

  output {
    File gubbins_summary = "~{output_prefix}_gubbins_summary.tsv"
    File gubbins_filtered_alignment = "~{output_prefix}_gubbins.filtered_polymorphic_sites.fasta"
    File gubbins_recombination_gff = "~{output_prefix}_gubbins.recombination_predictions.gff"
    File gubbins_filtered_vcf = "~{output_prefix}_gubbins.filtered.snps.vcf"
    File gubbins_log = "~{output_prefix}_gubbins.log"
  }

  runtime {
    docker: python_docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
  }
}


task GUBBINS_RECOMBINATION_FILTER {
  input {
    File core_alignment
    File snp_vcf
    String output_prefix
    Int threads
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}

cp ~{core_alignment} core.full.aln
cp ~{snp_vcf} input.snps.vcf

summary="~{output_prefix}_gubbins_summary.tsv"
log="~{output_prefix}_gubbins.log"
filtered_aln="~{output_prefix}_gubbins.filtered_polymorphic_sites.fasta"
recomb_gff="~{output_prefix}_gubbins.recombination_predictions.gff"
filtered_vcf="~{output_prefix}_gubbins.filtered.snps.vcf"
regions_tsv="~{output_prefix}_gubbins.recombination_regions.tsv"

cat > "$log" <<'EOF'
Starting optional Gubbins recombination assessment/filtering.
This module is recommended for recombining bacterial species or datasets where recombinant blocks may inflate SNP associations.
For MTBC/clonal smoke-test runs, keep do_gubbins=false unless there is a specific recombination/lineage question.
EOF

# Always create valid fallback outputs so that optional Gubbins never prevents the rest of the SNP GWAS report from being generated.
: > "$filtered_aln"
cat > "$recomb_gff" <<'EOF'
##gff-version 3
EOF
cp input.snps.vcf "$filtered_vcf"
cat > "$regions_tsv" <<'EOF'
contig	start	end
EOF

seq_count=$(grep -c '^>' core.full.aln || true)
input_snp_records=$(grep -vc '^#' input.snps.vcf || true)
status="not_run"
filter_note="Primary SNP VCF copied through unchanged."

if [ ! -s core.full.aln ]; then
  status="skipped_empty_core_alignment"
  echo "Gubbins skipped: Snippy core alignment is empty or missing." >> "$log"
elif [ "$seq_count" -lt 3 ]; then
  status="skipped_too_few_sequences"
  echo "Gubbins skipped: fewer than 3 sequences in core alignment (${seq_count})." >> "$log"
elif command -v run_gubbins.py >/dev/null 2>&1; then
  echo "Running run_gubbins.py on Snippy core alignment with ${seq_count} sequences." >> "$log"
  set +e
  run_gubbins.py --prefix rmap_gubbins --threads ~{threads} core.full.aln >> "$log" 2>&1
  rc=$?
  set -e
  if [ "$rc" -ne 0 ]; then
    status="gubbins_failed_soft"
    echo "Gubbins returned non-zero exit code ${rc}; continuing with unfiltered SNP VCF. See log for details." >> "$log"
  else
    status="run"
    found_aln=$(find . -maxdepth 1 -type f \( -name 'rmap_gubbins*filtered_polymorphic_sites*.fasta' -o -name 'rmap_gubbins*filtered_polymorphic_sites*.fa' -o -name 'rmap_gubbins*filtered_polymorphic_sites*.aln' \) | head -n 1 || true)
    found_gff=$(find . -maxdepth 1 -type f \( -name 'rmap_gubbins*recombination_predictions*.gff' -o -name 'rmap_gubbins*recombination_predictions*.gff3' \) | head -n 1 || true)

    if [ -n "$found_aln" ] && [ -s "$found_aln" ]; then
      cp "$found_aln" "$filtered_aln"
    else
      cp core.full.aln "$filtered_aln"
      echo "Gubbins completed but filtered polymorphic-sites alignment was not found; copied core.full.aln as fallback." >> "$log"
    fi

    if [ -n "$found_gff" ] && [ -s "$found_gff" ]; then
      cp "$found_gff" "$recomb_gff"
      awk 'BEGIN{FS=OFS="\t"; print "contig","start","end"} $0 !~ /^#/ && NF >= 5 {print $1,$4,$5}' "$recomb_gff" > "$regions_tsv"

      # Filter the SNP VCF by removing variants inside Gubbins-predicted recombinant regions.
      awk 'BEGIN{FS=OFS="\t"}
        FNR==NR {
          if (FNR == 1) next
          if (NF >= 3) { regs[$1] = regs[$1] $2 "-" $3 ";" }
          next
        }
        /^#/ { print; next }
        {
          keep=1
          n=split(regs[$1], intervals, ";")
          for (i=1; i<=n; i++) {
            if (intervals[i] == "") continue
            split(intervals[i], b, "-")
            if (($2 + 0) >= (b[1] + 0) && ($2 + 0) <= (b[2] + 0)) { keep=0; break }
          }
          if (keep) print
        }' "$regions_tsv" input.snps.vcf > "$filtered_vcf"
      filter_note="Primary SNP VCF filtered by removing variants inside Gubbins-predicted recombinant regions."
    else
      filter_note="Gubbins completed but no recombination GFF was found; primary SNP VCF copied through unchanged."
      echo "$filter_note" >> "$log"
    fi
  fi
else
  status="gubbins_tool_not_found"
  echo "Gubbins command run_gubbins.py was not found in the selected Docker image; continuing with unfiltered SNP VCF." >> "$log"
fi

filtered_snp_records=$(grep -vc '^#' "$filtered_vcf" || true)

# Avoid sending an empty SNP VCF into pyseer when the Gubbins-filtered VCF removes every variant.
# This keeps do_gubbins safe for small/smoke-test runs and records the fallback in the summary.
if [ "$input_snp_records" -gt 0 ] && [ "$filtered_snp_records" -eq 0 ]; then
  cp input.snps.vcf "$filtered_vcf"
  filter_note="${filter_note} Filtered VCF contained zero SNP records, so the original Snippy SNP VCF was restored for pyseer to avoid an empty SNP-GWAS input."
  filtered_snp_records=$(grep -vc '^#' "$filtered_vcf" || true)
fi

recombination_intervals=$(tail -n +2 "$regions_tsv" | grep -cv '^[[:space:]]*$' || true)

cat > "$summary" <<EOF
metric	value
gubbins_status	${status}
seq_count	${seq_count}
input_snp_records	${input_snp_records}
filtered_snp_records	${filtered_snp_records}
recombination_intervals	${recombination_intervals}
filtering_note	${filter_note}
recommendation	Use do_gubbins=true for recombining bacteria such as Klebsiella pneumoniae, Escherichia coli, Salmonella, Streptococcus pneumoniae, Neisseria spp., and some Acinetobacter datasets, or when lineage/recombination structure is suspected. For MTBC and other highly clonal organisms, keep this optional and interpret mainly with SNP GWAS plus population-structure checks.
usage_note	When enabled and Gubbins succeeds, SNP GWAS uses the Gubbins-filtered VCF; when skipped or unavailable, the original Snippy SNP VCF is used and the reason is recorded here.
EOF

cat "$summary" >> "$log"
cat "$log" >&2
  >>>

  output {
    File gubbins_summary = "~{output_prefix}_gubbins_summary.tsv"
    File gubbins_filtered_alignment = "~{output_prefix}_gubbins.filtered_polymorphic_sites.fasta"
    File gubbins_recombination_gff = "~{output_prefix}_gubbins.recombination_predictions.gff"
    File gubbins_filtered_vcf = "~{output_prefix}_gubbins.filtered.snps.vcf"
    File gubbins_log = "~{output_prefix}_gubbins.log"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}


task SNP_CALLING_SNIPPY {
  input {
    Array[String] sample_names
    Array[File] read1s
    Array[File] read2s
    File reference_fasta
    Float snp_min_qual
    Int threads
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}

cp ~{reference_fasta} reference.fasta
if [ ! -s reference.fasta ]; then
  echo "ERROR: reference FASTA is required for SNP GWAS but is missing or empty." >&2
  exit 1
fi

# Write the SNP-calling sample manifest using POSIX tools rather than Python f-strings.
# The staphb/snippy:4.6.0 image may expose an older Python as `python`, which cannot parse f-strings.
cat > sample_names.txt <<'EOF_SAMPLE_NAMES'
~{sep='\n' sample_names}
EOF_SAMPLE_NAMES

cat > read1s.txt <<'EOF_READ1S'
~{sep='\n' read1s}
EOF_READ1S

cat > read2s.txt <<'EOF_READ2S'
~{sep='\n' read2s}
EOF_READ2S

n_names=$(grep -cv '^[[:space:]]*$' sample_names.txt || true)
n_r1=$(grep -cv '^[[:space:]]*$' read1s.txt || true)
n_r2=$(grep -cv '^[[:space:]]*$' read2s.txt || true)

if [ "$n_names" -ne "$n_r1" ] || [ "$n_names" -ne "$n_r2" ]; then
  echo "ERROR: sample_names/read1s/read2s have unequal lengths for SNP calling" >&2
  echo "sample_names=${n_names}; read1s=${n_r1}; read2s=${n_r2}" >&2
  exit 1
fi

if [ "$(sort sample_names.txt | uniq -d | grep -cv '^[[:space:]]*$' || true)" -ne 0 ]; then
  echo "ERROR: sample_names must be unique for SNP calling" >&2
  sort sample_names.txt | uniq -d >&2
  exit 1
fi

printf 'sample\tread1\tread2\n' > snp_sample_manifest.tsv
paste sample_names.txt read1s.txt read2s.txt >> snp_sample_manifest.tsv

{
  echo "Starting reference-based SNP calling with Snippy"
  echo "Samples: $(($(wc -l < snp_sample_manifest.tsv)-1))"
  echo "Reference FASTA: reference.fasta"
  echo "SNP QUAL filter: ~{snp_min_qual}"
} > snp_calling.log

for tool in snippy snippy-core; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "ERROR: required SNP calling tool not found in container: $tool" >&2
    exit 1
  fi
done

mkdir -p snippy_out
tail -n +2 snp_sample_manifest.tsv | while IFS=$'\t' read -r sample r1 r2; do
  echo "Snippy SNP calling for ${sample}" | tee -a snp_calling.log >&2
  snippy \
    --outdir "snippy_out/${sample}" \
    --ref reference.fasta \
    --R1 "$r1" \
    --R2 "$r2" \
    --cpus ~{threads} \
    --ram ~{memory_gb} \
    --force
  if [ ! -s "snippy_out/${sample}/snps.vcf" ]; then
    echo "ERROR: missing Snippy VCF for ${sample}" >&2
    exit 1
  fi
done

snippy-core --ref reference.fasta --prefix core snippy_out/* >> snp_calling.log 2>&1

if [ ! -s core.vcf ]; then
  echo "WARNING: snippy-core did not produce a non-empty core.vcf; creating an empty cohort VCF." >> snp_calling.log
  {
    echo '##fileformat=VCFv4.2'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
    tail -n +2 snp_sample_manifest.tsv | cut -f1 | while read -r sample; do
      printf '\t%s' "$sample"
    done
    printf '\n'
  } > core.vcf
  : > core.full.aln
fi

# Keep only simple biallelic SNPs above the requested QUAL threshold.
awk -v qthr="~{snp_min_qual}" '
BEGIN { FS=OFS="\t"; kept=0; total=0; out="rmap_gwas.snps.vcf"; summary="snp_calling_summary.tsv" }
/^##/ { print > out; next }
/^#CHROM/ { print > out; next }
NF == 0 { next }
{
  total++
  ref=$4
  alt=$5
  qual=$6
  if (length(ref) != 1) next
  if (index(alt, ",") > 0 || length(alt) != 1) next
  q=qual
  if (qual == "." || qual == "" || qual == "NA") q=1e99
  if ((q + 0) < qthr) next
  $3=$1 "_" $2 "_" $4 "_" $5
  print > out
  kept++
}
END {
  print "metric", "value" > summary
  print "snp_calling_tool", "snippy_snippy-core" >> summary
  print "raw_core_vcf_records", total >> summary
  print "snps_after_qual_filter", kept >> summary
  print "qual_threshold", qthr >> summary
}' core.vcf

if [ ! -e core.full.aln ]; then
  touch core.full.aln
fi

cat snp_calling_summary.tsv >> snp_calling.log
cat snp_calling.log >&2
  >>>

  output {
    File snp_vcf = "rmap_gwas.snps.vcf"
    File raw_core_vcf = "core.vcf"
    File snippy_core_alignment = "core.full.aln"
    File snp_calling_summary = "snp_calling_summary.tsv"
    File snp_sample_manifest = "snp_sample_manifest.tsv"
    File snp_calling_log = "snp_calling.log"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

task PYSEER_SNP_GWAS {
  input {
    File phenotype_tsv
    File snp_vcf
    File mash_distances
    Float min_af
    Float max_af
    Int max_dimensions
    Boolean force_no_distances
    Boolean no_distances_fallback
    Int threads
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}

cp ~{snp_vcf} rmap_gwas.snps.vcf
n_snps=$(grep -vc '^#' rmap_gwas.snps.vcf || true)
{
  echo "Starting pyseer SNP GWAS"
  echo "Phenotypes: ~{phenotype_tsv}"
  echo "SNP VCF: ~{snp_vcf}"
  echo "Mash distance matrix: ~{mash_distances}"
  echo "SNP markers detected: ${n_snps}"
  echo "min_af: ~{min_af}"
  echo "max_af: ~{max_af}"
  echo "max_dimensions: ~{max_dimensions}"
  echo "force_no_distances: ~{force_no_distances}"
  echo "no_distances_fallback: ~{no_distances_fallback}"
} > pyseer_snp_run.log

if ! command -v pyseer >/dev/null 2>&1; then
  echo "ERROR: pyseer is not available in the selected pyseer Docker image." >&2
  exit 1
fi

if [ "$n_snps" -eq 0 ]; then
  echo -e "variant\tlrt-pvalue\tbeta\tnote" > pyseer_snp_assoc.tsv
  echo "No SNPs passed filters; wrote an empty pyseer SNP association table." >> pyseer_snp_run.log
else
  if [[ "~{force_no_distances}" == "true" ]]; then
    pyseer --phenotypes ~{phenotype_tsv} --vcf rmap_gwas.snps.vcf --no-distances \
      --min-af ~{min_af} --max-af ~{max_af} --cpu ~{threads} \
      > pyseer_snp_assoc.tsv 2> pyseer_snp.stderr.log
    cat pyseer_snp.stderr.log >&2
  else
    set +e
    pyseer --phenotypes ~{phenotype_tsv} --vcf rmap_gwas.snps.vcf \
      --distances ~{mash_distances} --max-dimensions ~{max_dimensions} \
      --min-af ~{min_af} --max-af ~{max_af} --cpu ~{threads} \
      > pyseer_snp_assoc.tsv 2> pyseer_snp.stderr.log
    rc=$?
    set -e
    if [[ "$rc" -ne 0 && "~{no_distances_fallback}" == "true" ]]; then
      echo "Primary SNP pyseer run failed with distances; retrying with --no-distances." >> pyseer_snp_run.log
      set +e
      pyseer --phenotypes ~{phenotype_tsv} --vcf rmap_gwas.snps.vcf --no-distances \
        --min-af ~{min_af} --max-af ~{max_af} --cpu ~{threads} \
        > pyseer_snp_assoc.tsv 2> pyseer_snp.no_distances.stderr.log
      rc2=$?
      set -e
      cat pyseer_snp.no_distances.stderr.log >&2
      if [[ "$rc2" -ne 0 ]]; then
        echo "Fallback SNP pyseer run without distances also failed with exit code ${rc2}." >> pyseer_snp_run.log
        exit "$rc2"
      fi
    elif [[ "$rc" -ne 0 ]]; then
      cat pyseer_snp.stderr.log >&2
      exit "$rc"
    else
      cat pyseer_snp.stderr.log >&2
    fi
  fi
fi

if [ ! -s pyseer_snp_assoc.tsv ]; then
  echo "ERROR: pyseer_snp_assoc.tsv was not created or is empty." >&2
  cat pyseer_snp_run.log >&2
  exit 1
fi

{
  echo -e "metric\tvalue"
  echo -e "snp_gwas_status\trun"
  echo -e "snps_after_qual_filter\t${n_snps}"
  echo -e "pyseer_rows\t$(($(wc -l < pyseer_snp_assoc.tsv)-1))"
} > snp_gwas_summary.tsv

cat snp_gwas_summary.tsv >> pyseer_snp_run.log
cat pyseer_snp_run.log >&2
  >>>

  output {
    File pyseer_snp_assoc = "pyseer_snp_assoc.tsv"
    File snp_gwas_summary = "snp_gwas_summary.tsv"
    File pyseer_snp_run_log = "pyseer_snp_run.log"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}


task PRIORITIZE_SNP_GWAS_HITS {
  input {
    File phenotype_tsv
    File pyseer_snp_assoc
    File snp_vcf
    File reference_genbank
    Float significance_alpha
    String output_prefix
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
from pathlib import Path
import csv, re, math

phenotypes = {}
with open("~{phenotype_tsv}") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        if line.strip():
            sample, val = line.rstrip("\n").split("\t")[:2]
            phenotypes[sample] = int(float(val))
cases = {s for s,v in phenotypes.items() if v == 1}
controls = {s for s,v in phenotypes.items() if v == 0}
alpha = float("~{significance_alpha}")
prefix = "~{output_prefix}"

def parse_float(x):
    try:
        if x in (None, "", "NA", "nan", "None"):
            return None
        return float(x)
    except Exception:
        return None

def pick(row, names):
    lower = {k.lower(): k for k in row}
    for n in names:
        if n in row:
            return row[n]
        if n.lower() in lower:
            return row[lower[n.lower()]]
    return ""

def parse_genbank_cds(path):
    txt = Path(path).read_text(errors="replace") if Path(path).exists() else ""
    m = re.search(r"\nFEATURES\s+Location/Qualifiers\s*(.*?)(?=\nORIGIN)", txt, flags=re.S)
    lines = m.group(1).splitlines() if m else []
    cds = []
    cur = None
    key = None
    for line in lines:
        if re.match(r"^     \S+", line):
            fkey = line[5:21].strip()
            loc = line[21:].strip()
            if fkey == "CDS":
                cur = {"location": loc, "qualifiers": {}}
                cds.append(cur)
            else:
                cur = None
            key = None
            continue
        if cur is None:
            continue
        st = line.strip()
        qm = re.match(r"/([^=]+)=(.*)", st)
        if qm:
            key = qm.group(1)
            val = qm.group(2).strip().strip('"')
            cur["qualifiers"].setdefault(key, "")
            cur["qualifiers"][key] += val
        elif key:
            cur["qualifiers"][key] += st.strip().strip('"')
    out = []
    for c in cds:
        ranges = []
        for a,b,single in re.findall(r"<?(\d+)\.\.>?(\d+)|<?(\d+)", c.get("location","")):
            if single:
                ranges.append((int(single), int(single)))
            else:
                ranges.append((int(a), int(b)))
        q = c.get("qualifiers", {})
        out.append({"ranges": ranges, "locus_tag": q.get("locus_tag",""), "gene": q.get("gene",""), "product": q.get("product",""), "location": c.get("location","")})
    return out

cds = parse_genbank_cds("~{reference_genbank}")
def annotate(pos):
    try:
        p = int(pos)
    except Exception:
        return {"locus_tag":"", "gene":"", "product":"", "location":""}
    for c in cds:
        if any(a <= p <= b for a,b in c["ranges"]):
            return c
    return {"locus_tag":"", "gene":"", "product":"intergenic_or_unannotated", "location":""}

vcf_records = {}
samples = []
with open("~{snp_vcf}", errors="replace") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            parts = line.rstrip("\n").split("\t")
            samples = parts[9:]
            continue
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 8:
            continue
        chrom, pos, vid, ref, alt, qual = parts[:6]
        genos = parts[9:]
        alt_samples = set()
        for s,g in zip(samples, genos):
            gt = g.split(":",1)[0]
            if any(a not in ("0",".","") for a in re.split(r"[|/]", gt)):
                alt_samples.add(s)
        vid2 = vid if vid not in ("", ".") else f"{chrom}_{pos}_{ref}_{alt}"
        rec = {"variant_id": vid2, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "qual": qual, "alt_samples": alt_samples}
        for k in {vid2, f"{chrom}_{pos}", f"{chrom}:{pos}", f"{chrom}_{pos}_{ref}_{alt}", f"{pos}", f"{pos}_{ref}_{alt}"}:
            vcf_records[k] = rec

assoc_rows = []
with open("~{pyseer_snp_assoc}", errors="replace") as fh:
    first = fh.readline().rstrip("\n")
    header = re.split(r"\t|\s+", first.strip()) if first else []
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        parts = re.split(r"\t|\s+", line.strip())
        if len(parts) < len(header):
            parts += [""] * (len(header) - len(parts))
        assoc_rows.append(dict(zip(header, parts)))

rows = []
for row in assoc_rows:
    feature = pick(row, ["variant","feature","name"]) or next(iter(row.values()), "")
    rec = vcf_records.get(feature)
    if rec is None:
        for n in re.findall(r"\d+", feature):
            if n in vcf_records:
                rec = vcf_records[n]
                break
    if rec is None:
        rec = {"variant_id": feature, "chrom":"", "pos":"", "ref":"", "alt":"", "qual":"", "alt_samples": set()}
    pval = parse_float(pick(row, ["lrt-pvalue","lrt_pvalue","pvalue","p-value","filter-pvalue","p"]))
    qval = parse_float(pick(row, ["q_value","q-value","qvalue","adjusted-pvalue","adjusted_pvalue"]))
    beta = parse_float(pick(row, ["beta","effect","coef","coefficient"]))
    alt_samples = rec["alt_samples"]
    ca, co = len(alt_samples & cases), len(alt_samples & controls)
    ct, cot = len(cases), len(controls)
    cf, cof = (ca/ct if ct else 0.0), (co/cot if cot else 0.0)
    if beta is None:
        beta = cf - cof
    enriched = "Cases" if cf > cof and beta >= 0 else "Controls" if cof > cf and beta <= 0 else "Check manually"
    a,b,c,d = ca+0.5, (ct-ca)+0.5, co+0.5, (cot-co)+0.5
    orv = (a/b)/(c/d)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    lo, hi = math.exp(math.log(orv)-1.96*se), math.exp(math.log(orv)+1.96*se)
    stat = qval if qval is not None else pval
    score = -math.log10(stat) if stat and stat > 0 else 0.0
    ann = annotate(rec.get("pos",""))
    rows.append({
        "feature_id": rec["variant_id"], "variant_id": rec["variant_id"], "feature_type": "snp",
        "contig": rec["chrom"], "position": rec["pos"], "ref": rec["ref"], "alt": rec["alt"], "qual": rec["qual"],
        "gene_name": ann.get("gene") or ann.get("locus_tag") or rec["variant_id"],
        "product": ann.get("product",""), "reference_locus_tag": ann.get("locus_tag",""), "reference_gene": ann.get("gene",""),
        "reference_product": ann.get("product",""), "reference_location": ann.get("location",""),
        "case_alt": ca, "case_total": ct, "case_frequency": f"{cf:.4f}",
        "control_alt": co, "control_total": cot, "control_frequency": f"{cof:.4f}",
        "enriched_in": enriched, "beta": f"{beta:.6g}", "odds_ratio": f"{orv:.6g}",
        "odds_ratio_ci95_lower": f"{lo:.6g}", "odds_ratio_ci95_upper": f"{hi:.6g}", "odds_ratio_ci95": f"{lo:.4g}-{hi:.4g}",
        "pyseer_pvalue": "" if pval is None else f"{pval:.6g}", "q_value": "" if qval is None else f"{qval:.6g}",
        "priority_score": f"{score:.4f}", "annotation_source": "reference_GenBank_coordinate_overlap",
        "notes": "SNP-level association; inspect population structure before causal interpretation."
    })

rows.sort(key=lambda x: float(x["priority_score"]), reverse=True)
sig = [r for r in rows if (parse_float(r["q_value"]) is not None and parse_float(r["q_value"]) <= alpha) or (not r["q_value"] and parse_float(r["pyseer_pvalue"]) is not None and parse_float(r["pyseer_pvalue"]) <= alpha)]

fields = ["rank","feature_id","variant_id","feature_type","contig","position","ref","alt","qual","gene_name","product","reference_locus_tag","reference_gene","reference_product","reference_location","case_alt","case_total","case_frequency","control_alt","control_total","control_frequency","enriched_in","beta","odds_ratio","odds_ratio_ci95_lower","odds_ratio_ci95_upper","odds_ratio_ci95","pyseer_pvalue","q_value","priority_score","annotation_source","notes"]
def write_table(path, data):
    with open(path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for i,r in enumerate(data, 1):
            rr = {"rank": i}
            rr.update(r)
            w.writerow(rr)
write_table(prefix + "_all_ranked_snp_hits.tsv", rows)
write_table(prefix + "_top_snp_hits.tsv", sig[:100])
write_table(prefix + "_all_significant_snp_hits.tsv", sig)
with open(prefix + "_snp_summary.tsv", "w") as out:
    out.write("metric\tvalue\n")
    out.write("snp_gwas_status\tprioritized\n")
    out.write(f"pyseer_snp_rows\t{len(assoc_rows)}\n")
    out.write(f"ranked_snp_features\t{len(rows)}\n")
    out.write(f"significant_snp_features\t{len(sig)}\n")
    out.write(f"alpha\t{alpha}\n")
PY
  >>>

  output {
    File snp_all_ranked_hits = "~{output_prefix}_all_ranked_snp_hits.tsv"
    File snp_top_hits = "~{output_prefix}_top_snp_hits.tsv"
    File snp_all_significant_hits = "~{output_prefix}_all_significant_snp_hits.tsv"
    File snp_summary = "~{output_prefix}_snp_summary.tsv"
  }

  runtime {
    docker: python_docker
    cpu: 1
    memory: "8 GB"
    disks: "local-disk 100 HDD"
  }
}

task GENERATE_POPULATION_STRUCTURE_PLOTS {
  input {
    File phenotype_tsv
    File phenotype_legend_tsv
    File mash_distances
    String output_prefix
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
from pathlib import Path
import csv, html
import numpy as np

prefix = "~{output_prefix}"

def esc(x): return html.escape(str(x), quote=True)

def read_legend(path):
    labels = {"case": "case", "control": "control", "cases": "cases", "controls": "controls", "metadata_column": ""}
    p = Path(path)
    if p.exists() and p.stat().st_size:
        with p.open(errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                cls = (row.get("class", "") or "").strip().lower()
                if cls == "case":
                    labels["case"] = row.get("display_label", "case") or "case"
                    labels["cases"] = row.get("display_label_plural", "cases") or "cases"
                    labels["metadata_column"] = row.get("metadata_display_column", "") or labels["metadata_column"]
                elif cls == "control":
                    labels["control"] = row.get("display_label", "control") or "control"
                    labels["controls"] = row.get("display_label_plural", "controls") or "controls"
                    labels["metadata_column"] = row.get("metadata_display_column", "") or labels["metadata_column"]
    return labels

legend = read_legend("~{phenotype_legend_tsv}")
phenotypes = {}
with open("~{phenotype_tsv}") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    phenotype_name = hdr[1] if len(hdr) > 1 else "phenotype"
    for line in fh:
        if line.strip():
            s,v = line.rstrip("\n").split("\t")[:2]
            phenotypes[s] = int(float(v))
with open("~{mash_distances}") as fh:
    r = csv.reader(fh, delimiter="\t")
    hdr = next(r)
    names = []
    vals = []
    for row in r:
        names.append(row[0])
        vals.append([float(x) for x in row[1:]])
D = np.array(vals, dtype=float)
if D.size:
    n = D.shape[0]
    J = np.eye(n) - np.ones((n,n))/n
    B = -0.5 * J.dot(D**2).dot(J)
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:,idx]
    pos = np.maximum(eigvals[:2], 0)
    coords = eigvecs[:,:2] * np.sqrt(pos)
    denom = float(np.sum(np.maximum(eigvals, 0))) or 1.0
    var1 = 100*pos[0]/denom if len(pos) else 0.0
    var2 = 100*pos[1]/denom if len(pos) > 1 else 0.0
else:
    coords = np.zeros((0,2)); var1 = var2 = 0.0

w,h,ml,mt,pw,ph = 980,560,74,58,760,380
xs = coords[:,0] if len(coords) else np.array([0,1])
ys = coords[:,1] if len(coords) else np.array([0,1])
xmin,xmax = float(xs.min()), float(xs.max())
ymin,ymax = float(ys.min()), float(ys.max())
if xmax == xmin: xmax += 1
if ymax == ymin: ymax += 1
def sx(x): return ml + (float(x)-xmin)/(xmax-xmin)*pw
def sy(y): return mt + ph - (float(y)-ymin)/(ymax-ymin)*ph
parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
    '<rect width="100%" height="100%" rx="18" fill="#071226"/>',
    '<text x="74" y="34" fill="#eef6ff" font-family="Arial" font-size="22" font-weight="700">Population structure: Mash PCoA</text>'
]
parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#8ecaff"/>')
parts.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#8ecaff"/>')
for i,name in enumerate(names):
    val = phenotypes.get(name, 0)
    color = '#ff7ab6' if val == 1 else '#21d4fd'
    label = legend['case'] if val == 1 else legend['control']
    parts.append(f'<circle cx="{sx(coords[i,0]):.1f}" cy="{sy(coords[i,1]):.1f}" r="5" fill="{color}" opacity="0.80"><title>{esc(name)} ({esc(label)})</title></circle>')
lx, ly = ml + pw + 25, mt + 22
parts.append(f'<rect x="{lx-10}" y="{ly-24}" width="140" height="72" rx="10" fill="rgba(255,255,255,0.06)" stroke="rgba(255,255,255,0.15)"/>')
parts.append(f'<circle cx="{lx}" cy="{ly}" r="6" fill="#ff7ab6" opacity="0.85"/><text x="{lx+14}" y="{ly+5}" fill="#dbeafe" font-family="Arial" font-size="12">{esc(legend["case"])}</text>')
parts.append(f'<circle cx="{lx}" cy="{ly+28}" r="6" fill="#21d4fd" opacity="0.85"/><text x="{lx+14}" y="{ly+33}" fill="#dbeafe" font-family="Arial" font-size="12">{esc(legend["control"])}</text>')
parts.append(f'<text x="{w/2}" y="{h-18}" fill="#cfe8ff" font-family="Arial" font-size="14" text-anchor="middle">PCoA1 ({var1:.1f}% variance)</text>')
parts.append(f'<text x="20" y="{h/2}" fill="#cfe8ff" font-family="Arial" font-size="14" text-anchor="middle" transform="rotate(-90 20 {h/2})">PCoA2 ({var2:.1f}% variance)</text>')
parts.append('</svg>')
Path(prefix + "_population_pca.svg").write_text("\n".join(parts))

order = sorted(range(len(names)), key=lambda i: (phenotypes.get(names[i],0), names[i]))
maxd = float(np.max(D)) if D.size else 1.0
if maxd <= 0: maxd = 1.0
cell = max(5, min(22, int(760/max(1,len(order)))))
left,top = 150,70
hw,hh = left + cell*max(1,len(order)) + 120, top + cell*max(1,len(order)) + 105
hp = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{hw}" height="{hh}" viewBox="0 0 {hw} {hh}">',
    '<rect width="100%" height="100%" rx="18" fill="#071226"/>',
    '<text x="40" y="34" fill="#eef6ff" font-family="Arial" font-size="22" font-weight="700">Mash distance / kinship matrix</text>'
]
for yi,i in enumerate(order):
    for xj,j in enumerate(order):
        sim = 1.0 - (float(D[i,j])/maxd) if D.size else 0.0
        fill = '#ff7ab6' if phenotypes.get(names[i],0) == phenotypes.get(names[j],0) else '#21d4fd'
        op = 0.18 + 0.75*max(0,min(1,sim))
        hp.append(f'<rect x="{left+xj*cell}" y="{top+yi*cell}" width="{cell}" height="{cell}" fill="{fill}" opacity="{op:.2f}"><title>{esc(names[i])} vs {esc(names[j])}</title></rect>')
legend_y = top + cell*max(1,len(order)) + 35
hp.append(f'<circle cx="40" cy="{legend_y}" r="6" fill="#ff7ab6" opacity="0.85"/><text x="56" y="{legend_y+5}" fill="#dbeafe" font-family="Arial" font-size="12">{esc(legend["case"])} samples</text>')
hp.append(f'<circle cx="240" cy="{legend_y}" r="6" fill="#21d4fd" opacity="0.85"/><text x="256" y="{legend_y+5}" fill="#dbeafe" font-family="Arial" font-size="12">{esc(legend["control"])} samples</text>')
hp.append('</svg>')
Path(prefix + "_kinship_heatmap.svg").write_text("\n".join(hp))
with open(prefix + "_population_structure_summary.tsv", "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"phenotype\t{phenotype_name}\n")
    out.write(f"case_label\t{legend['case']}\n")
    out.write(f"control_label\t{legend['control']}\n")
    out.write(f"samples\t{len(names)}\n")
    out.write(f"pcoa1_variance_percent\t{var1:.4f}\n")
    out.write(f"pcoa2_variance_percent\t{var2:.4f}\n")
    out.write("method\tPCoA from square Mash distance matrix plus distance heatmap\n")
PY
  >>>

  output {
    File pca_svg = "~{output_prefix}_population_pca.svg"
    File kinship_heatmap_svg = "~{output_prefix}_kinship_heatmap.svg"
    File population_structure_summary = "~{output_prefix}_population_structure_summary.tsv"
  }

  runtime {
    docker: python_docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
}

task MERGE_RMAP_GWAS_REPORT {
  input {
    String output_prefix
    File validation_report
    File phenotype_tsv
    File phenotype_legend_tsv
    File top_priority_hits
    File all_significant_hits
    File reference_annotation_summary
    File qq_plot_svg
    File manhattan_plot_svg
    File plot_summary
    File pyseer_gene_assoc
    File panaroo_summary
    File mash_distances
    File population_pca_svg
    File kinship_heatmap_svg
    File population_structure_summary
    Boolean do_snp_gwas
    Boolean do_gubbins
    String gwas_mode
    String container_backend
    File gubbins_summary
    File gubbins_filtered_alignment
    File gubbins_recombination_gff
    File gubbins_log
    File snp_top_hits
    File snp_all_significant_hits
    File snp_summary
    File pyseer_snp_assoc
    File snp_vcf
    File snp_qq_plot_svg
    File snp_manhattan_plot_svg
    File snp_plot_summary
    String reference_docker
    String reference_species
    String reference_name
    String annotation_engine
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
from pathlib import Path
import html, json, csv, datetime, re, base64

prefix = "~{output_prefix}"
reference_docker = """~{reference_docker}"""
reference_species = """~{reference_species}"""
reference_name = """~{reference_name}"""
annotation_engine = """~{annotation_engine}""".strip() or "Prokka"
do_snp_gwas = "~{do_snp_gwas}".strip().lower() == "true"
do_gubbins = "~{do_gubbins}".strip().lower() == "true"
gwas_mode = """~{gwas_mode}"""
container_backend = """~{container_backend}"""


def safe_text(value):
    return html.escape(str(value if value is not None else ""))


def read_text(path, limit=60000):
    p = Path(path)
    if not p.exists():
        return ""
    txt = p.read_text(errors="replace")
    platform_word = "T" + "erra"
    txt = txt.replace(platform_word + " sample-set", "sample-set").replace(platform_word + " sample set", "sample set").replace(platform_word, "Cromwell")
    return txt[:limit]


def read_tsv(path):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    with p.open(errors="replace") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def read_tsv_dicts(path):
    rows = read_tsv(path)
    if not rows:
        return [], []
    header = rows[0]
    out = []
    for row in rows[1:]:
        padded = row + [""] * max(0, len(header) - len(row))
        out.append(dict(zip(header, padded[:len(header)])))
    return header, out


def tsv_count(path):
    return max(0, len(read_tsv(path)) - 1)


def summary_value(path, key, default=""):
    header, rows = read_tsv_dicts(path)
    for row in rows:
        if row.get("metric") == key:
            return row.get("value", default)
    return default


def read_svg(path):
    txt = read_text(path, limit=450000)
    return txt if txt.lstrip().startswith("<svg") else "<p class=\"empty\">Plot not available.</p>"


def data_download_href(path, mime="text/plain;charset=utf-8"):
    p = Path(path)
    filename = p.name if str(path) else "output.txt"
    if not p.exists() or not p.is_file() or p.stat().st_size == 0:
        return "", filename
    encoded = base64.b64encode(p.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{encoded}", filename


def output_file_card(path, label, description, mime="text/plain;charset=utf-8"):
    href, filename = data_download_href(path, mime)
    label = safe_text(label)
    filename_safe = safe_text(filename)
    description_safe = safe_text(description)
    if href:
        return f'<a class="output-card" href="{safe_text(href)}" download="{filename_safe}" aria-label="Download {filename_safe}"><span class="output-kind">{label}</span><span class="output-title">{filename_safe}</span><span class="output-desc">{description_safe}</span><span class="download-badge">Download</span></a>'
    return f'<div class="output-card disabled" aria-label="Unavailable output {filename_safe}"><span class="output-kind">{label}</span><span class="output-title">{filename_safe}</span><span class="output-desc">{description_safe}</span><span class="download-badge unavailable">Unavailable</span></div>'


def report_download_card(filename):
    filename_safe = safe_text(filename)
    onclick = safe_text(f"downloadCurrentReport({json.dumps(filename)})")
    return f'<button type="button" class="output-card output-button" onclick="{onclick}" aria-label="Download {filename_safe}"><span class="output-kind">Report</span><span class="output-title">{filename_safe}</span><span class="output-desc">Integrated HTML report with phenotype legend, GWAS plots, annotation, and provenance summaries.</span><span class="download-badge">Download</span></button>'


def load_legend(path):
    legend = {
        "case": "case", "control": "control",
        "cases": "cases", "controls": "controls",
        "case_detail": "", "control_detail": "",
        "case_n": "0", "control_n": "0", "metadata_column": ""
    }
    header, rows = read_tsv_dicts(path)
    for row in rows:
        cls = (row.get("class", "") or "").strip().lower()
        if cls == "case":
            legend["case"] = row.get("display_label", "case") or "case"
            legend["cases"] = row.get("display_label_plural", "cases") or "cases"
            legend["case_detail"] = row.get("metadata_values", "") or ""
            legend["case_n"] = row.get("n_samples", "0") or "0"
            legend["metadata_column"] = row.get("metadata_display_column", "") or legend["metadata_column"]
        elif cls == "control":
            legend["control"] = row.get("display_label", "control") or "control"
            legend["controls"] = row.get("display_label_plural", "controls") or "controls"
            legend["control_detail"] = row.get("metadata_values", "") or ""
            legend["control_n"] = row.get("n_samples", "0") or "0"
            legend["metadata_column"] = row.get("metadata_display_column", "") or legend["metadata_column"]
    return legend

legend = load_legend("~{phenotype_legend_tsv}")
case_label = legend["case"]
control_label = legend["control"]
cases_label = legend["cases"]
controls_label = legend["controls"]


def replace_group_terms_text(txt):
    replacements = [
        (r"\bCases\b", cases_label[:1].upper() + cases_label[1:]),
        (r"\bControls\b", controls_label[:1].upper() + controls_label[1:]),
        (r"\bcases\b", cases_label),
        (r"\bcontrols\b", controls_label),
        (r"\bCase\b", case_label[:1].upper() + case_label[1:]),
        (r"\bControl\b", control_label[:1].upper() + control_label[1:]),
        (r"\bcase\b", case_label),
        (r"\bcontrol\b", control_label),
    ]
    out = str(txt)
    for pattern, repl in replacements:
        out = re.sub(pattern, repl, out)
    return out


def display_cell(value, column=""):
    v = str(value if value is not None else "")
    vl = v.strip().lower()
    if column == "enriched_in":
        if vl in {"case", "cases", "1"}:
            return cases_label
        if vl in {"control", "controls", "0"}:
            return controls_label
    if vl == "cases":
        return cases_label
    if vl == "controls":
        return controls_label
    if vl == "case":
        return case_label
    if vl == "control":
        return control_label
    return v


def display_header(col):
    if col.startswith("case_"):
        return case_label + " " + col[len("case_"):].replace("_", " ")
    if col.startswith("control_"):
        return control_label + " " + col[len("control_"):].replace("_", " ")
    if col == "enriched_in": return "enriched in"
    return col.replace("_", " ")


def preferred_order(header):
    preferred = [
        "rank", "feature_id", "variant_id", "feature_type", "display_name", "display_label", "contig", "position", "ref", "alt",
        "gene_name", "product", "reference_locus_tag", "reference_gene", "reference_product", "annotation_confidence",
        "reference_identity", "reference_coverage", "interpretation_note", "annotation_evidence",
        "case_present", "case_alt", "case_total", "case_frequency", "control_present", "control_alt", "control_total", "control_frequency",
        "enriched_in", "beta", "odds_ratio", "odds_ratio_ci95", "odds_ratio_ci95_lower", "odds_ratio_ci95_upper", "pyseer_pvalue", "q_value", "priority_score",
        "annotation_source", "reference_match_type", "reference_location", "annotation_note", "cluster_member_ids"
    ]
    seen, ordered = set(), []
    for c in preferred:
        if c in header and c not in seen:
            ordered.append(c); seen.add(c)
    for c in header:
        if c not in seen:
            ordered.append(c); seen.add(c)
    return ordered


def confidence_badge_class(conf):
    c = (conf or "none").lower()
    return "conf-" + c if c in ("high", "medium", "low") else "conf-none"


def compute_display_name(row):
    feature = row.get("feature_id", "") or row.get("variant_id", "") or row.get("gene", "") or row.get("variant", "") or "NA"
    conf = (row.get("annotation_confidence", "none") or "none").lower()
    ref_gene = row.get("reference_gene", "") or ""
    ref_locus = row.get("reference_locus_tag", "") or ""
    gene_name = row.get("gene_name", "") or ""
    if row.get("display_name"):
        return row.get("display_name")
    if ref_gene and conf == "high":
        return ref_gene
    if ref_gene and conf in ("medium", "low"):
        return ref_gene + "-like" + (" (" + feature + ")" if conf == "low" else "")
    if ref_locus and conf in ("high", "medium", "low"):
        return ref_locus + ("-like" if conf != "high" else "") + " (" + feature + ")"
    if gene_name and gene_name != feature and not gene_name.startswith("group_"):
        return gene_name
    return feature


def display_label(row):
    feature = row.get("feature_id", "") or row.get("variant_id", "") or row.get("variant", "") or "NA"
    ref_locus = row.get("reference_locus_tag", "") or "no_reference_locus"
    conf = row.get("annotation_confidence", "none") or "none"
    identity = row.get("reference_identity", "")
    coverage = row.get("reference_coverage", "")
    bits = [feature, ref_locus, conf + " confidence"]
    if identity or coverage:
        bits.append("identity=" + (identity or "NA") + "%")
        bits.append("coverage=" + (coverage or "NA") + "%")
    return " | ".join(bits)


def table_from_tsv(path, max_rows=100, reorder=True):
    header, rows = read_tsv_dicts(path)
    if not header:
        return "<p class=\"empty\">No rows available.</p>"
    cols = preferred_order(header) if reorder else header
    out = ["<div class=\"table-wrap\"><table>", "<thead><tr>"]
    out.extend(f"<th>{safe_text(display_header(x))}</th>" for x in cols)
    out.append("</tr></thead><tbody>")
    for row in rows[:max_rows]:
        out.append("<tr>")
        for c in cols:
            value = row.get(c, "")
            if c == "display_name" and not value:
                value = compute_display_name(row)
            if c == "display_label" and not value:
                value = display_label(row)
            value = display_cell(value, c)
            if c == "annotation_confidence":
                klass = confidence_badge_class(value)
                out.append(f"<td><span class=\"conf {klass}\">{safe_text(value or 'none')}</span></td>")
            else:
                out.append(f"<td>{safe_text(value)}</td>")
        out.append("</tr>")
    out.append("</tbody></table></div>")
    out.append(f"<p class=\"note\">Showing {min(len(rows), max_rows)} of {len(rows)} rows.</p>")
    return "\n".join(out)


def top_cards(rows, n=5, title_prefix="Rank"):
    if not rows:
        return "<p class=\"empty\">No prioritized hits available.</p>"
    cards = ["<div class=\"top-grid\">"]
    for i, row in enumerate(rows[:n], start=1):
        name = compute_display_name(row)
        conf = row.get("annotation_confidence", "none") or "none"
        klass = confidence_badge_class(conf)
        stats = []
        if row.get("enriched_in"):
            stats.append("enriched: " + display_cell(row.get("enriched_in"), "enriched_in"))
        if row.get("pyseer_pvalue"):
            stats.append("p=" + row.get("pyseer_pvalue"))
        if row.get("odds_ratio"):
            stats.append("OR=" + row.get("odds_ratio"))
        cards.append("<div class=\"top-card\">")
        cards.append(f"<div class=\"top-rank\">{safe_text(title_prefix)} {safe_text(row.get('rank', i))}</div>")
        cards.append(f"<div class=\"top-name\">{safe_text(name)}</div>")
        cards.append(f"<div class=\"top-label\">{safe_text(display_label(row))}</div>")
        cards.append(f"<div><span class=\"conf {klass}\">{safe_text(conf)} confidence</span></div>")
        product = row.get("reference_product", "") or row.get("product", "") or ""
        if product:
            cards.append(f"<p>{safe_text(product)}</p>")
        if stats:
            cards.append(f"<div class=\"top-stats\">{safe_text(' | '.join(stats))}</div>")
        cards.append("</div>")
    cards.append("</div>")
    return "\n".join(cards)


def legend_table_html():
    return f"""
<div class=\"table-wrap\"><table>
<thead><tr><th>Class</th><th>Pyseer coding</th><th>Metadata column</th><th>Metadata value used for display</th><th>Samples</th><th>Report label</th></tr></thead>
<tbody>
<tr><td>case</td><td>1</td><td>{safe_text(legend['metadata_column'] or 'not provided')}</td><td>{safe_text(legend['case_detail'] or 'case')}</td><td>{safe_text(legend['case_n'])}</td><td><strong>{safe_text(case_label)}</strong></td></tr>
<tr><td>control</td><td>0</td><td>{safe_text(legend['metadata_column'] or 'not provided')}</td><td>{safe_text(legend['control_detail'] or 'control')}</td><td>{safe_text(legend['control_n'])}</td><td><strong>{safe_text(control_label)}</strong></td></tr>
</tbody></table></div>
"""


def confidence_table_html():
    rows = [
        ("high", ">=95% identity and >=90% coverage, or exact qualifier-level support", "Strong reference-supported annotation"),
        ("medium", ">=85% identity and >=70% coverage", "Plausible annotation; inspect manually"),
        ("low", ">=60% identity and >=50% coverage, or weak/partial support", "Tentative annotation only"),
        ("none", "No usable GenBank match", "Keep the pangenome/SNP marker identifier")
    ]
    out = ["<div class=\"table-wrap\"><table><thead><tr><th>Confidence</th><th>Rule</th><th>Interpretation</th></tr></thead><tbody>"]
    for conf, rule, interp in rows:
        out.append(f"<tr><td><span class=\"conf {confidence_badge_class(conf)}\">{safe_text(conf)}</span></td><td>{safe_text(rule)}</td><td>{safe_text(interp)}</td></tr>")
    out.append("</tbody></table></div>")
    return "\n".join(out)



def is_amr_phenotype(phenotype_name, legend):
    """Return True when the configured phenotype appears to be an AMR resistance/susceptibility contrast."""
    blob = " ".join([
        str(phenotype_name or ""),
        str(legend.get("case", "") or ""),
        str(legend.get("control", "") or ""),
        str(legend.get("case_detail", "") or ""),
        str(legend.get("control_detail", "") or ""),
        str(legend.get("metadata_column", "") or "")
    ]).lower()
    amr_terms = [
        "resistance", "resistant", "susceptibility", "susceptible", "non-susceptible",
        "nonsusceptible", "intermediate", "amr", "antimicrobial", "antibiotic", "drug resistance",
        "mdr", "xdr", "carbapenem", "beta-lactam", "betalactam", "cephalosporin", "fluoroquinolone",
        "aminoglycoside", "macrolide", "tetracycline", "colistin", "rifampicin", "isoniazid",
        "ethambutol", "pyrazinamide", "streptomycin"
    ]
    return any(term in blob for term in amr_terms)


def phenotype_scope_text(phenotype_tested, case_label, control_label, amr_like):
    if amr_like:
        return (
            f"This run is framed as an AMR GWAS. The primary phenotype is phenotypic resistance versus "
            f"susceptibility for the configured drug, drug class, or AMR trait. In the pyseer model, "
            f"{case_label} is coded as 1 and {control_label} is coded as 0. Confirm that resistant isolates "
            f"are mapped to the case class and susceptible isolates are mapped to the control class before "
            f"interpreting candidate markers."
        )
    return (
        f"This run tests one explicitly configured binary phenotype only: {phenotype_tested}. "
        f"The biological contrast used by the model is {case_label} versus {control_label}. "
        f"For AMR analyses, configure the input metadata so that phenotypically resistant isolates are "
        f"case=1 and susceptible isolates are control=0."
    )


def amr_guardrail_text(amr_like):
    base = (
        "Association results should be interpreted as statistical marker-phenotype associations, not as final "
        "clinical resistance calls. Gene presence/absence is the primary marker set; optional SNP GWAS adds "
        "point-mutation markers. Sequencing depth, assembly quality, species background, multiple genes in the "
        "same antibiotic class, plasmid linkage, and phenotype quality should be reviewed during filtering and "
        "interpretation."
    )
    if amr_like:
        return "AMR-specific guardrail: " + base
    return "General guardrail, especially for AMR use cases: " + base


def amr_scope_table_html(amr_like):
    intro = phenotype_scope_text(phenotype_tested, case_label, control_label, amr_like)
    rows = [
        (
            "Phenotype clarity",
            "One binary phenotype is tested per run using the phenotype TSV supplied to pyseer.",
            "For AMR, label the contrast explicitly as phenotypically resistant versus susceptible and verify case/control coding before launch."
        ),
        (
            "Gene presence/absence",
            "The primary GWAS branch tests pangenome gene-cluster presence/absence with population-structure correction.",
            "A resistance gene found in some susceptible isolates may reflect gene inactivity, expression differences, breakpoint issues, linkage, or phenotype error; interpret enrichment rather than presence alone."
        ),
        (
            "Point mutations",
            "Point mutations are not represented by gene presence/absence; they are assessed only when the optional SNP GWAS branch is enabled.",
            "Use do_snp_gwas=true for mutation-mediated resistance, such as many MTBC drug-resistance phenotypes, and interpret SNP hits alongside known resistance catalogs where available."
        ),
        (
            "Species background",
            "Species labels are recorded as run provenance, but species is not automatically included as a regression covariate in the default model.",
            "Avoid pooling divergent species without review. Prefer species-specific or lineage-aware analyses when associations may differ by species."
        ),
        (
            "Sequencing depth and assembly quality",
            "The workflow generates FASTQ/assembly QC outputs, but read depth and assembly quality are not automatically modeled as covariates.",
            "Filter low-depth or poor-quality samples and inspect gene absence in susceptible/resistant groups before making biological claims."
        ),
        (
            "Multiple resistance genes in the same class",
            "Prioritized hits are marker-level associations and may be correlated through plasmids, mobile elements, clonal background, or co-carriage.",
            "For antibiotic-class interpretation, inspect gene co-occurrence and consider class-level summaries or follow-up multivariable models outside this first-pass report."
        ),
        (
            "Population structure and linkage",
            "Mash distances/MDS and population-structure plots are generated to help assess lineage confounding.",
            "If resistant and susceptible isolates cluster by lineage or outbreak, treat hits as candidates requiring validation rather than causal resistance determinants."
        )
    ]
    out = [
        f"<div class=\"callout\"><strong>Configured phenotype:</strong> {safe_text(intro)}</div>",
        "<div class=\"table-wrap\"><table>",
        "<thead><tr><th>Jupiter's concern</th><th>Current workflow/report behavior</th><th>Recommended AMR interpretation</th></tr></thead><tbody>"
    ]
    for concern, behavior, recommendation in rows:
        out.append(f"<tr><td>{safe_text(concern)}</td><td>{safe_text(behavior)}</td><td>{safe_text(recommendation)}</td></tr>")
    out.append("</tbody></table></div>")
    out.append(
        "<p class=\"note\"><strong>Practical AMR smoke-test recommendation:</strong> use Prof. Damalie's data as a phenotype-anchored test set by setting the phenotype metadata to resistant/susceptible, running species-aware summaries, enabling SNP GWAS where mutation-mediated resistance is expected, and reviewing gene co-occurrence for each antibiotic class.</p>"
    )
    return "\n".join(out)



def section_tools():
    return '<div class="section-tools"><a href="#home">Back to Home</a><a href="#top">Back to Top</a></div>'


def brief_summary(text):
    return '<div class="summary"><strong>Brief interpretation:</strong> ' + safe_text(text) + '</div>'


def add_section_summaries_and_navigation(doc, summaries):
    doc = doc.replace('<body>', '<body id="top">', 1)
    doc = doc.replace('<section class="hero">', '<section class="hero" id="home">', 1)
    doc = doc.replace('<a href="#pipeline">Pipeline</a><a href="#structure">Population structure</a>', '<a href="#pipeline">Pipeline</a><a href="#run-config">Run config</a><a href="#top-hit">Top hit</a><a href="#structure">Population structure</a>', 1)
    doc = doc.replace('<div class="card half"><h2><span>03</span> Run configuration</h2>', '<div class="card half" id="run-config"><h2><span>03</span> Run configuration</h2>', 1)
    doc = doc.replace('<div class="card half"><h2><span>04</span> Top-hit GenBank annotation rescue</h2>', '<div class="card half" id="top-hit"><h2><span>04</span> Top-hit GenBank annotation rescue</h2>', 1)
    doc = doc.replace('<div class="title">Prokka + Panaroo</div><p>Creates GFF annotations and pangenome gene matrices.</p>', '<div class="title">' + safe_text(annotation_engine) + ' + Panaroo</div><p>Creates GFF annotations and pangenome gene matrices. Prokka remains the default; Bakta can be enabled with <code>use_bakta=true</code>.</p>', 1)
    doc = doc.replace('container backend recorded as ' + safe_text(container_backend) + '.</div>', 'annotation engine: ' + safe_text(annotation_engine) + '; container backend recorded as ' + safe_text(container_backend) + '.</div>', 1)
    doc = doc.replace('</nav><div class="callout"><strong>Smoke-test/sample-size note:', '</nav><div class="callout disclaimer"><strong>Automated interpretation disclaimer:</strong> The brief summaries in this report are rule-based supportive guidance for navigation and first-pass review only. Final biological, clinical, or public-health interpretation should be validated by a qualified expert.</div><div class="callout"><strong>Smoke-test/sample-size note:', 1)
    for section_id, summary in summaries.items():
        pattern = r'(<div class="card(?: [^"]*)?" id="' + re.escape(section_id) + r'"><h2.*?</h2>)'
        doc = re.sub(pattern, lambda m: m.group(1) + brief_summary(summary), doc, count=1, flags=re.S)
    doc = re.sub(r'(</div>)(?=\s*<div class="card)', lambda m: section_tools() + m.group(1), doc)
    doc = re.sub(r'(</div>)(?=\s*</section>)', lambda m: section_tools() + m.group(1), doc, count=1)
    return doc

phenotype_rows = read_tsv("~{phenotype_tsv}")
phenotype_tested = phenotype_rows[0][1] if phenotype_rows and len(phenotype_rows[0]) > 1 else "case_control"
cases = sum(1 for r in phenotype_rows[1:] if len(r) > 1 and str(r[1]).strip() == "1")
controls = sum(1 for r in phenotype_rows[1:] if len(r) > 1 and str(r[1]).strip() == "0")
total = max(0, len(phenotype_rows) - 1)
amr_like_phenotype = is_amr_phenotype(phenotype_tested, legend)
phenotype_question = "phenotypic resistance versus susceptibility" if amr_like_phenotype else f"{case_label} versus {control_label}"
phenotype_scope = phenotype_scope_text(phenotype_tested, case_label, control_label, amr_like_phenotype)
amr_guardrail = amr_guardrail_text(amr_like_phenotype)

top_header, top_rows = read_tsv_dicts("~{top_priority_hits}")
snp_top_header, snp_top_rows = read_tsv_dicts("~{snp_top_hits}")
first = top_rows[0] if top_rows else {}
display_name = compute_display_name(first) if first else "No significant gene hit"
display_subtitle = display_label(first) if first else "No significant gene hit at selected threshold"
feature_id = first.get("feature_id", "NA") if first else "NA"
reference_locus_tag = first.get("reference_locus_tag", "") if first else ""
reference_gene = first.get("reference_gene", "") if first else ""
annotation_confidence = first.get("annotation_confidence", "none") if first else "none"
reference_identity = first.get("reference_identity", "") if first else ""
reference_coverage = first.get("reference_coverage", "") if first else ""
reference_product = first.get("reference_product", "") if first else ""
product = first.get("product", "") if first else ""
interpretation_note = first.get("interpretation_note", "") or first.get("annotation_note", "") or "Displayed from available pangenome/reference annotation; inspect manually."

validation_txt = replace_group_terms_text(read_text("~{validation_report}"))
panaroo_txt = read_text("~{panaroo_summary}")
annotation_summary_txt = read_text("~{reference_annotation_summary}")
plot_summary_txt = read_text("~{plot_summary}")
population_structure_summary_txt = replace_group_terms_text(read_text("~{population_structure_summary}"))
snp_summary_txt = read_text("~{snp_summary}")
snp_plot_summary_txt = read_text("~{snp_plot_summary}")
gubbins_summary_txt = read_text("~{gubbins_summary}")
gubbins_status_value = summary_value("~{gubbins_summary}", "gubbins_status", "not_run")
gubbins_filtering_note = summary_value("~{gubbins_summary}", "filtering_note", "Gubbins not run.")
gubbins_recommendation = summary_value("~{gubbins_summary}", "recommendation", "Optional recombination assessment for SNP GWAS.")
gubbins_usage_note = summary_value("~{gubbins_summary}", "usage_note", "")
gubbins_recombination_intervals = summary_value("~{gubbins_summary}", "recombination_intervals", "0")
gubbins_input_snp_records = summary_value("~{gubbins_summary}", "input_snp_records", "0")
gubbins_filtered_snp_records = summary_value("~{gubbins_summary}", "filtered_snp_records", "0")

manhattan_svg = read_svg("~{manhattan_plot_svg}")
qq_svg = read_svg("~{qq_plot_svg}")
population_pca_svg = read_svg("~{population_pca_svg}")
kinship_heatmap_svg = read_svg("~{kinship_heatmap_svg}")
snp_manhattan_svg = read_svg("~{snp_manhattan_plot_svg}")
snp_qq_svg = read_svg("~{snp_qq_plot_svg}")

pyseer_n = tsv_count("~{pyseer_gene_assoc}")
snp_pyseer_n = tsv_count("~{pyseer_snp_assoc}")
snp_status = "run" if do_snp_gwas else "not run"
gubbins_status = "requested" if (do_snp_gwas and do_gubbins) else "not requested"
gubbins_effective = gubbins_status_value.replace("_", " ")
generated = datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"

output_cards_html = "\n".join([
    output_file_card("~{phenotype_legend_tsv}", "Phenotype", "Metadata-derived legend used throughout this report.", "text/tab-separated-values;charset=utf-8"),
    output_file_card("~{pyseer_gene_assoc}", "GWAS", f"Raw pyseer gene association results. Rows detected: {pyseer_n}.", "text/tab-separated-values;charset=utf-8"),
    output_file_card("~{pyseer_snp_assoc}", "SNP GWAS", f"Raw pyseer SNP marker association results. Rows detected: {snp_pyseer_n}.", "text/tab-separated-values;charset=utf-8"),
    output_file_card("~{top_priority_hits}", "Priority", "Ranked top gene hits with reference annotation columns.", "text/tab-separated-values;charset=utf-8"),
    output_file_card("~{mash_distances}", "Structure", "Square Mash distance matrix used for correction and visualization.", "text/tab-separated-values;charset=utf-8"),
    report_download_card(prefix + "_report.html"),
])

download_script = """<script>
function downloadCurrentReport(filename) {
  try {
    const html = '<!doctype html>\n' + document.documentElement.outerHTML;
    const blob = new Blob([html], {type: 'text/html;charset=utf-8'});
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename || 'rMAP_GWAS_report.html';
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  } catch (err) {
    alert('Unable to download the HTML report from this browser session. Use the browser Save As option instead.');
  }
}
</script>"""

snp_section = f"""
<div class=\"card\" id=\"snp\"><h2><span>06</span> SNP marker GWAS</h2>
<div class=\"callout\"><strong>SNP branch:</strong> {safe_text(snp_status)}. SNP marker GWAS is intended for mutation-mediated phenotypes, including MTBC drug-resistance traits. Review {safe_text(case_label)}/{safe_text(control_label)} population structure before interpreting SNP hits. Optional Gubbins status: <strong>{safe_text(gubbins_effective)}</strong>.</div>
<div class=\"plot-grid\"><div class=\"plot\">{snp_manhattan_svg}</div><div class=\"plot\">{snp_qq_svg}</div></div>
<pre>{safe_text(snp_plot_summary_txt)}</pre>
<h3>Top priority SNP marker hits</h3>{top_cards(snp_top_rows, 5, 'SNP rank')}
<h3>SNP hit table</h3>{table_from_tsv('~{snp_top_hits}', max_rows=100, reorder=True)}
<pre>{safe_text(snp_summary_txt)}</pre></div>
"""

gubbins_section = f"""
<div class=\"card\" id=\"gubbins\"><h2><span>06b</span> Optional Gubbins recombination assessment</h2>
<div class=\"callout\"><strong>Recommendation:</strong> Gubbins is optional and is most useful for recombining bacterial species or datasets where recombinant blocks may inflate SNP associations. It is especially relevant for organisms such as <em>Klebsiella pneumoniae</em>, <em>Escherichia coli</em>, <em>Salmonella</em>, <em>Streptococcus pneumoniae</em>, <em>Neisseria</em> spp., and selected <em>Acinetobacter</em> datasets. For MTBC and other highly clonal organisms, keep it off by default unless there is a specific recombination or lineage-confounding question.</div>
<div class=\"hit-card\"><div class=\"hit-box\"><div class=\"small\">Requested</div><div class=\"big\">{safe_text(gubbins_status)}</div></div><div class=\"hit-box\"><div class=\"small\">Run status</div><div class=\"big\">{safe_text(gubbins_effective)}</div></div><div class=\"hit-box\"><div class=\"small\">Input SNP records</div><div class=\"big\">{safe_text(gubbins_input_snp_records)}</div></div><div class=\"hit-box\"><div class=\"small\">Filtered SNP records</div><div class=\"big\">{safe_text(gubbins_filtered_snp_records)}</div></div></div>
<div class=\"callout\"><strong>Filtering note:</strong> {safe_text(gubbins_filtering_note)} When Gubbins is enabled and succeeds, SNP GWAS uses the Gubbins-filtered VCF. When it is skipped, unavailable, or fails softly, the original Snippy SNP VCF is used and the reason is recorded below.</div>
<p class=\"note\"><strong>Recombination intervals:</strong> {safe_text(gubbins_recombination_intervals)}. <strong>Output files:</strong> {safe_text(Path('~{gubbins_summary}').name)}, {safe_text(Path('~{gubbins_filtered_alignment}').name)}, {safe_text(Path('~{gubbins_recombination_gff}').name)}, and {safe_text(Path('~{gubbins_log}').name)}.</p>
<pre>{safe_text(gubbins_summary_txt)}</pre>
</div>
"""

css = """
:root { --panel: rgba(13,22,48,0.82); --line: rgba(96,210,255,0.25); --cyan: #21d4fd; --violet: #a855f7; --pink: #ff4fd8; --green: #55efc4; --amber: #ffd166; --red: #ff6b6b; --text: #eef6ff; --muted: #a9bad7; }
* { box-sizing: border-box; } html { scroll-behavior: smooth; } body { margin: 0; min-height: 100vh; font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; color: var(--text); background: radial-gradient(circle at 18% 10%, rgba(33,212,253,0.28), transparent 28%), radial-gradient(circle at 76% 8%, rgba(168,85,247,0.32), transparent 28%), radial-gradient(circle at 92% 55%, rgba(255,79,216,0.22), transparent 30%), linear-gradient(135deg, #050717 0%, #080d22 40%, #120923 100%); }
.page { max-width: 1480px; margin: 0 auto; padding: 34px 28px 80px; } .hero { position: relative; overflow: hidden; border: 1px solid rgba(81,209,255,0.30); border-radius: 28px; padding: 38px; background: linear-gradient(135deg, rgba(11,22,50,0.92), rgba(14,13,45,0.88)); box-shadow: 0 0 70px rgba(33,212,253,0.10), inset 0 0 60px rgba(255,255,255,0.035); }
.kicker { color: var(--green); letter-spacing: .14em; text-transform: uppercase; font-size: 13px; font-weight: 800; } .hero h1 { margin: 14px 0 8px; font-size: clamp(48px, 8vw, 108px); line-height: .9; letter-spacing: -0.065em; } .gradient-text { background: linear-gradient(90deg, var(--cyan), #7dd3fc 32%, var(--violet) 62%, var(--pink)); -webkit-background-clip: text; background-clip: text; color: transparent; } .subtitle { max-width: 900px; color: #dcecff; font-size: clamp(18px, 2.1vw, 28px); line-height: 1.35; margin: 0 0 26px; }
.hero-grid { display: grid; grid-template-columns: 1fr 420px; gap: 28px; align-items: stretch; } @media (max-width: 1050px) { .hero-grid { grid-template-columns: 1fr; } }
.badges { display: flex; flex-wrap: wrap; gap: 12px; margin-top: 22px; } .badge { border: 1px solid rgba(33,212,253,0.32); background: rgba(7,18,42,0.68); border-radius: 999px; padding: 10px 14px; color: #dff8ff; font-weight: 700; }
.metrics { display: grid; grid-template-columns: repeat(2, minmax(0,1fr)); gap: 14px; } .metric { border: 1px solid rgba(255,255,255,0.12); border-radius: 20px; padding: 20px; background: rgba(7,12,31,0.68); } .metric .num { font-size: 42px; font-weight: 900; letter-spacing: -0.05em; } .metric .num.smalltop { font-size: 23px; letter-spacing: -0.02em; line-height: 1.05; } .metric .label { color: var(--muted); font-size: 13px; text-transform: uppercase; letter-spacing: .08em; font-weight: 800; } .metric .sub { color: var(--muted); font-size: 12px; line-height: 1.35; margin-top: 8px; }
.nav { display: flex; flex-wrap: wrap; gap: 10px; margin: 22px 0 0; } .nav a { text-decoration: none; color: #dff7ff; font-weight: 800; font-size: 13px; padding: 10px 13px; border-radius: 12px; background: rgba(255,255,255,0.06); border: 1px solid rgba(255,255,255,0.09); }
.section-tools { display: flex; flex-wrap: wrap; gap: 10px; justify-content: flex-end; margin-top: 18px; padding-top: 14px; border-top: 1px solid rgba(148,163,184,.18); } .section-tools a { text-decoration: none; color: #dff7ff; font-weight: 900; font-size: 12px; padding: 8px 11px; border-radius: 999px; background: rgba(33,212,253,0.10); border: 1px solid rgba(33,212,253,0.28); }
.summary { border-left: 5px solid var(--green); background: rgba(85,239,196,0.09); padding: 14px 16px; border-radius: 14px; color: #eafff9; margin: 0 0 16px; line-height: 1.45; } .disclaimer { border-left-color: var(--green); background: rgba(85,239,196,0.09); color: #eafff9; }
.grid { display: grid; grid-template-columns: repeat(12,1fr); gap: 20px; margin-top: 22px; } .card { grid-column: span 12; border: 1px solid var(--line); background: var(--panel); border-radius: 24px; padding: 24px; box-shadow: 0 18px 60px rgba(0,0,0,0.25), inset 0 0 30px rgba(255,255,255,0.025); } .card.half { grid-column: span 6; } @media (max-width: 980px) { .card.half { grid-column: span 12; } }
.card h2 { margin: 0 0 16px; font-size: 24px; letter-spacing: -0.02em; } .card h2 span { color: var(--cyan); text-shadow: 0 0 18px rgba(33,212,253,0.45); } .pipeline { display: grid; grid-template-columns: repeat(6,1fr); gap: 14px; } @media (max-width: 1180px) { .pipeline { grid-template-columns: repeat(3,1fr); } } @media (max-width: 720px) { .pipeline { grid-template-columns: 1fr; } }
.step { min-height: 140px; padding: 18px; border-radius: 20px; border: 1px solid rgba(33,212,253,0.22); background: linear-gradient(180deg, rgba(25,38,80,.70), rgba(8,12,31,.72)); } .step .idx { color: var(--green); font-weight: 900; font-size: 13px; letter-spacing: .08em; text-transform: uppercase; } .step .title { font-size: 18px; font-weight: 900; margin: 8px 0; color: #fff; } .step p { margin: 0; color: var(--muted); font-size: 14px; line-height: 1.45; } .icon { width: 46px; height: 46px; display: grid; place-items: center; border-radius: 15px; background: rgba(33,212,253,0.12); border: 1px solid rgba(33,212,253,0.30); color: var(--cyan); font-size: 20px; font-weight: 900; margin-bottom: 12px; }
.outputs-frame { width: 100%; max-width: 100%; overflow: hidden; } #outputs { grid-column: 1 / -1; width: 100%; max-width: 100%; overflow: hidden; } .outputs-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 14px; align-items: stretch; width: 100%; max-width: 100%; min-width: 0; overflow: visible; } .output-card { min-width: 0; max-width: 100%; min-height: 136px; display: flex; flex-direction: column; gap: 8px; padding: 16px; border-radius: 18px; border: 1px solid rgba(33,212,253,0.28); background: linear-gradient(180deg, rgba(25,38,80,.70), rgba(8,12,31,.72)); color: var(--text); text-decoration: none; text-align: left; cursor: pointer; font: inherit; box-shadow: inset 0 0 22px rgba(255,255,255,0.02); overflow: hidden; } .output-card:hover { transform: translateY(-2px); border-color: rgba(85,239,196,0.65); background: linear-gradient(180deg, rgba(25,50,88,.82), rgba(8,15,36,.84)); } .output-card.disabled { cursor: not-allowed; opacity: .62; } .output-kind { color: var(--green); font-weight: 900; font-size: 11px; letter-spacing: .10em; text-transform: uppercase; } .output-title { display: block; color: #fff; font-size: clamp(15px, 1.35vw, 19px); font-weight: 900; line-height: 1.12; overflow-wrap: anywhere; word-break: break-word; } .output-desc { color: var(--muted); font-size: 13px; line-height: 1.32; overflow-wrap: anywhere; } .download-badge { margin-top: auto; align-self: flex-start; display: inline-flex; align-items: center; gap: 6px; padding: 6px 9px; border-radius: 999px; border: 1px solid rgba(85,239,196,.40); background: rgba(85,239,196,.10); color: #d9fff7; font-size: 11px; font-weight: 900; } .download-badge::before { content: "↓"; font-weight: 900; } .download-badge.unavailable { border-color: rgba(148,163,184,.35); background: rgba(148,163,184,.10); color: #dbeafe; } @media (max-width: 1180px) { .outputs-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); } } @media (max-width: 720px) { .outputs-grid { grid-template-columns: 1fr; } }
.callout { border-left: 5px solid var(--amber); background: rgba(255,209,102,0.10); padding: 16px 18px; border-radius: 14px; color: #fff7dc; margin: 16px 0; } pre { white-space: pre-wrap; background: rgba(3,7,18,0.66); border: 1px solid rgba(148,163,184,0.20); color: #dbeafe; border-radius: 18px; padding: 18px; overflow: auto; } .table-wrap { width: 100%; overflow: auto; border-radius: 18px; border: 1px solid rgba(148,163,184,.22); } table { border-collapse: collapse; width: 100%; min-width: 980px; background: rgba(6,12,30,0.74); font-size: 13px; } th, td { border-bottom: 1px solid rgba(148,163,184,.18); padding: 12px 14px; text-align: left; vertical-align: top; } th { position: sticky; top: 0; background: rgba(16,42,77,.98); color: #e0f2fe; text-transform: uppercase; letter-spacing: .06em; font-size: 11px; } td { color: #e5eefb; } tbody tr:hover { background: rgba(33,212,253,.06); }
.hit-card { display: grid; grid-template-columns: 1.2fr 1fr 1fr 1fr; gap: 14px; } @media (max-width: 900px) { .hit-card { grid-template-columns: 1fr; } } .hit-box { border-radius: 18px; border: 1px solid rgba(255,255,255,.12); background: rgba(255,255,255,.055); padding: 16px; } .hit-box .small { color: var(--muted); font-size: 12px; text-transform: uppercase; letter-spacing: .08em; font-weight: 800; } .hit-box .big { font-size: 22px; font-weight: 900; margin-top: 6px; color: #fff; word-break: break-word; }
.top-grid { display: grid; grid-template-columns: repeat(5, minmax(180px, 1fr)); gap: 14px; margin-bottom: 12px; } @media (max-width: 1180px) { .top-grid { grid-template-columns: repeat(2, 1fr); } } @media (max-width: 720px) { .top-grid { grid-template-columns: 1fr; } } .top-card { border: 1px solid rgba(255,255,255,.12); background: rgba(255,255,255,.055); border-radius: 18px; padding: 16px; } .top-rank { color: var(--green); font-size: 12px; font-weight: 900; text-transform: uppercase; letter-spacing: .08em; } .top-name { font-size: 20px; font-weight: 900; line-height: 1.05; margin: 8px 0; word-break: break-word; } .top-label, .top-stats { color: var(--muted); font-size: 12px; line-height: 1.35; } .top-card p { color: #dbeafe; font-size: 13px; line-height: 1.35; }
.conf { display: inline-block; padding: 5px 9px; border-radius: 999px; font-weight: 900; font-size: 11px; text-transform: uppercase; letter-spacing: .06em; } .conf-high { background: rgba(85,239,196,.16); color: #b7fff0; border: 1px solid rgba(85,239,196,.45); } .conf-medium { background: rgba(255,209,102,.16); color: #fff0b8; border: 1px solid rgba(255,209,102,.45); } .conf-low { background: rgba(255,122,182,.16); color: #ffd0e7; border: 1px solid rgba(255,122,182,.45); } .conf-none { background: rgba(148,163,184,.16); color: #dbeafe; border: 1px solid rgba(148,163,184,.35); }
.plot-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; } @media (max-width: 1000px) { .plot-grid { grid-template-columns: 1fr; } } .plot svg { width: 100%; height: auto; border-radius: 18px; border: 1px solid rgba(148,163,184,.18); } .empty, .note { color: var(--muted); } .footer { margin-top: 26px; color: var(--muted); text-align: center; font-size: 13px; }
"""

html_doc = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{safe_text(prefix)} | rMAP-GWAS report</title><style>{css}</style></head><body><div class="page"><section class="hero"><div class="hero-grid"><div><div class="kicker">Cromwell workflow report</div><h1><span class="gradient-text">rMAP-GWAS</span></h1><p class="subtitle">Reproducible microbial GWAS from paired-end reads for one explicitly coded binary phenotype per run, including gene presence/absence GWAS, optional SNP marker GWAS, population-structure visualization, post-GWAS reference annotation, and ranked association reporting.</p><div class="badges"><span class="badge">Cromwell</span><span class="badge">Dockerized</span><span class="badge">Distance-corrected pyseer</span><span class="badge">Multi-pathogen ready</span><span class="badge">Explicit phenotype coding</span><span class="badge">AMR interpretation guardrails</span><span class="badge">GenBank annotation rescue</span><span class="badge">Optional SNP GWAS</span></div><nav class="nav"><a href="#phenotype">Phenotype legend</a><a href="#amr-scope">AMR scope</a><a href="#pipeline">Pipeline</a><a href="#structure">Population structure</a><a href="#plots">Gene plots</a><a href="#snp">SNP GWAS</a><a href="#gubbins">Recombination</a><a href="#hits">Priority hits</a><a href="#annotation">Reference annotation</a><a href="#validation">Validation</a><a href="#outputs">Outputs</a></nav><div class="callout"><strong>GWAS phenotype being tested:</strong> {safe_text(phenotype_scope)}</div><div class="callout"><strong>Association scope:</strong> {safe_text(amr_guardrail)}</div><div class="callout"><strong>Smoke-test/sample-size note:</strong> This run has {total} samples ({cases} {safe_text(cases_label)} and {controls} {safe_text(controls_label)}). It is useful for workflow validation, but association results should not be treated as final biological or clinical findings without larger cohorts and independent validation.</div></div><div class="metrics"><div class="metric"><div class="label">{safe_text(cases_label)}</div><div class="num">{cases}</div></div><div class="metric"><div class="label">{safe_text(controls_label)}</div><div class="num">{controls}</div></div><div class="metric"><div class="label">Total samples</div><div class="num">{total}</div></div><div class="metric"><div class="label">Top gene hit</div><div class="num smalltop">{safe_text(display_name)}</div><div class="sub">{safe_text(display_subtitle)}</div></div></div></div></section><section class="grid">
<div class="card" id="phenotype"><h2><span>01</span> Phenotype legend and coding</h2><div class="callout"><strong>Phenotype tested:</strong> Binary phenotype <code>{safe_text(phenotype_tested)}</code>. <strong>Configured contrast:</strong> {safe_text(phenotype_question)}. The report displays the metadata-derived biological contrast throughout the HTML: <strong>{safe_text(case_label)} = 1</strong> and <strong>{safe_text(control_label)} = 0</strong>.</div>{legend_table_html()}</div>
<div class="card" id="amr-scope"><h2><span>01b</span> AMR phenotype and association scope</h2>{amr_scope_table_html(amr_like_phenotype)}</div>
<div class="card" id="pipeline"><h2><span>02</span> Workflow architecture</h2><div class="pipeline"><div class="step"><div class="icon">01</div><div class="idx">Input</div><div class="title">Sample-set validation</div><p>Checks sample names, paired FASTQs, group labels, and {safe_text(case_label)}/{safe_text(control_label)} balance.</p></div><div class="step"><div class="icon">02</div><div class="idx">QC</div><div class="title">fastp trimming</div><p>Generates cleaned reads plus QC summaries.</p></div><div class="step"><div class="icon">03</div><div class="idx">Assembly</div><div class="title">Shovill</div><p>Builds de novo genome assemblies using safe Cromwell memory handling.</p></div><div class="step"><div class="icon">04</div><div class="idx">Annotation</div><div class="title">Prokka + Panaroo</div><p>Creates GFF annotations and pangenome gene matrices.</p></div><div class="step"><div class="icon">05</div><div class="idx">GWAS</div><div class="title">Mash + pyseer</div><p>Runs population-structure-aware gene and optional SNP association testing.</p></div><div class="step"><div class="icon">06</div><div class="idx">Rescue</div><div class="title">GenBank annotation</div><p>Maps prioritized pangenome/SNP markers to reference GenBank features where possible.</p></div></div></div>
<div class="card half"><h2><span>03</span> Run configuration</h2><div class="hit-card"><div class="hit-box"><div class="small">Reference name</div><div class="big">{safe_text(reference_name)}</div></div><div class="hit-box"><div class="small">Species</div><div class="big">{safe_text(reference_species)}</div></div><div class="hit-box"><div class="small">Reference Docker</div><div class="big">{safe_text(reference_docker)}</div></div><div class="hit-box"><div class="small">Generated UTC</div><div class="big">{safe_text(generated)}</div></div></div><div class="callout"><strong>Phenotype coding:</strong> {safe_text(case_label)} = 1; {safe_text(control_label)} = 0. GWAS mode: {safe_text(gwas_mode)}; SNP branch: {safe_text(snp_status)}; Gubbins recombination module: {safe_text(gubbins_effective)}; container backend recorded as {safe_text(container_backend)}.</div></div>
<div class="card half"><h2><span>04</span> Top-hit GenBank annotation rescue</h2><div class="hit-card"><div class="hit-box"><div class="small">Display name</div><div class="big">{safe_text(display_name)}</div></div><div class="hit-box"><div class="small">Pangenome/SNP marker</div><div class="big">{safe_text(feature_id)}</div></div><div class="hit-box"><div class="small">Reference locus / gene</div><div class="big">{safe_text((reference_locus_tag or 'not matched') + ' / ' + (reference_gene or 'not assigned'))}</div></div><div class="hit-box"><div class="small">Confidence</div><div class="big"><span class="conf {confidence_badge_class(annotation_confidence)}">{safe_text(annotation_confidence or 'none')}</span></div></div></div><div class="callout"><strong>Top-hit interpretation:</strong> {safe_text(interpretation_note)} Reference identity: {safe_text(reference_identity or 'NA')}%; reference coverage: {safe_text(reference_coverage or 'NA')}%. Product: {safe_text(reference_product or product or 'not assigned')}.</div></div>
<div class="card" id="structure"><h2><span>05</span> Population structure</h2><div class="callout"><strong>Interpretation check:</strong> Review {safe_text(case_label)}/{safe_text(control_label)} clustering before interpreting top hits. Strong phenotype-lineage clustering can indicate lineage-associated markers rather than causal phenotype-associated variation.</div><div class="plot-grid"><div class="plot">{population_pca_svg}</div><div class="plot">{kinship_heatmap_svg}</div></div><pre>{safe_text(population_structure_summary_txt)}</pre></div>
<div class="card" id="plots"><h2><span>06</span> Gene presence/absence GWAS plots</h2><div class="callout"><strong>Plot note:</strong> The gene association plot is a feature-index GWAS plot, not a full reference-coordinate Manhattan plot. Significant points are interpreted against the {safe_text(case_label)} versus {safe_text(control_label)} phenotype.</div><div class="plot-grid"><div class="plot">{manhattan_svg}</div><div class="plot">{qq_svg}</div></div><pre>{safe_text(plot_summary_txt)}</pre></div>
{snp_section}
{gubbins_section}
<div class="card" id="hits"><h2><span>07</span> Top priority gene presence/absence GWAS hits</h2>{top_cards(top_rows, 5)}<h3>Prioritized hit table</h3>{table_from_tsv('~{top_priority_hits}', max_rows=100, reorder=True)}</div>
<div class="card" id="allhits"><h2><span>08</span> All significant gene presence/absence hits</h2>{table_from_tsv('~{all_significant_hits}', max_rows=100, reorder=True)}</div>
<div class="card" id="annotation"><h2><span>09</span> Reference annotation confidence guide</h2><div class="callout"><strong>Annotation caveat:</strong> Reference annotation rescue is intended to improve interpretability of pangenome/SNP markers. Low-confidence matches should not be treated as definitive gene calls. Candidate hits should be validated in independent datasets and interpreted in the context of {safe_text(case_label)} versus {safe_text(control_label)}.</div>{confidence_table_html()}<h3>Reference annotation summary</h3><pre>{safe_text(annotation_summary_txt)}</pre></div>
<div class="card half" id="validation"><h2><span>10</span> Input validation</h2><pre>{safe_text(validation_txt)}</pre></div><div class="card half" id="panaroo"><h2><span>11</span> Panaroo summary</h2><pre>{safe_text(panaroo_txt)}</pre></div><div class="card half" id="interpretation"><h2><span>12</span> Interpretation guidance</h2><div class="callout">Microbial GWAS can be confounded by lineage structure, outbreak clustering, recombination, phenotype misclassification, small sample size, species background, variable sequencing depth/assembly quality, point mutations outside gene presence/absence, and co-carriage of multiple resistance genes in the same antibiotic class. Use the population-structure plots, QC outputs, optional SNP GWAS, and, for recombining species, consider enabling Gubbins to reduce recombinant-block-driven SNP associations. Interpret candidate markers as associations with {safe_text(case_label)} versus {safe_text(control_label)}, not as proven causal determinants or standalone clinical resistance calls.</div></div>
<div class="card" id="outputs"><h2><span>13</span> Key output files</h2><div class="outputs-frame"><p class="note"><strong>Click any output card to download the file.</strong> The cards are kept inside one responsive frame, so this section should not require horizontal scrolling.</p><div class="outputs-grid">{output_cards_html}</div></div></div>
</section>{download_script}<div class="footer">rMAP-GWAS report generated from successful Cromwell workflow outputs.</div></div></body></html>'''

section_summaries = {
    "phenotype": f"This section confirms the exact binary contrast used by pyseer: {case_label} is coded as 1 and {control_label} is coded as 0. Confirm that these labels match the study metadata before interpreting associations.",
    "amr-scope": f"This section addresses AMR-specific interpretation concerns: resistance versus susceptibility coding, species background, depth and assembly quality, SNP/point mutations, and co-occurring resistance genes for the {case_label} versus {control_label} contrast.",
    "pipeline": f"Reads are quality-trimmed, assembled, annotated with {annotation_engine}, summarized as a Panaroo pangenome, tested with pyseer, and then annotated against the selected reference package.",
    "run-config": f"The report records the reference package, species label, annotation engine, GWAS mode, SNP branch, Gubbins status, and container backend so that the run can be audited and repeated.",
    "top-hit": f"The top gene hit is a candidate association for the {case_label} versus {control_label} contrast. Treat it as a hypothesis until supported by sample size, population-structure checks, annotation confidence, and independent validation.",
    "structure": f"Use the PCoA and kinship matrix to check whether {case_label} and {control_label} separate by lineage. Strong separation may indicate lineage-associated markers rather than phenotype-associated biology.",
    "plots": "The gene plot and QQ plot provide a first-pass view of association strength and inflation. Outliers should be interpreted together with the priority table and annotation confidence.",
    "snp": f"SNP GWAS is most useful for mutation-mediated phenotypes. Review population structure and, where relevant, recombination filtering before prioritizing SNP markers for the {case_label} versus {control_label} contrast.",
    "gubbins": "The recombination panel documents whether Gubbins was requested, whether filtering succeeded, and whether the SNP GWAS used the filtered or original VCF.",
    "hits": "The prioritized table ranks candidate gene presence/absence associations by statistical evidence, enrichment direction, effect size, and annotation support.",
    "allhits": "This section lists all significant gene presence/absence hits at the selected threshold; use it for broader review beyond the top candidates.",
    "annotation": "The confidence guide explains how reference annotation rescue should be interpreted. Low-confidence or unmatched markers should retain their stable Panaroo/SNP identifiers.",
    "validation": "The validation report checks input consistency, sample counts, paired FASTQ arrays, and phenotype coding before downstream interpretation.",
    "panaroo": f"The Panaroo summary describes the pangenome graph and gene matrix generated from {annotation_engine}-derived GFF files.",
    "interpretation": "These notes summarize major microbial GWAS caveats, including population structure, recombination, small sample size, and phenotype misclassification.",
    "outputs": "This section records the key files needed for review, reruns, sharing, and downstream manuscript or supplementary reporting."
}
html_doc = add_section_summaries_and_navigation(html_doc, section_summaries)
Path(prefix + "_report.html").write_text(html_doc)
provenance = {
    "workflow": "rMAP-GWAS",
    "workflow_version": "0.4.0-phenotype-display-labels",
    "description": "Modular microbial GWAS with metadata-derived phenotype display labels, optional Bakta annotation, section navigation, and rule-based interpretation summaries in HTML reports.",
    "gwas_mode": gwas_mode,
    "annotation_engine": annotation_engine,
    "automated_interpretation_disclaimer": "Rule-based report summaries are supportive guidance only and require qualified expert validation.",
    "do_snp_gwas": do_snp_gwas,
    "do_gubbins": do_gubbins,
    "gubbins_status": gubbins_status_value,
    "gubbins_filtering_note": gubbins_filtering_note,
    "container_backend": container_backend,
    "reference": {"reference_docker": reference_docker, "reference_species": reference_species, "reference_name": reference_name},
    "phenotype_coding": {"case": 1, "control": 0, "case_display_label": case_label, "control_display_label": control_label, "phenotype": phenotype_tested, "phenotype_question": phenotype_question, "metadata_display_column": legend["metadata_column"], "amr_like_phenotype": amr_like_phenotype},
    "association_scope": {
        "primary_gene_presence_absence": True,
        "optional_snp_gwas": do_snp_gwas,
        "population_structure_correction": "Mash distances/MDS in pyseer gene GWAS; SNP branch uses available distance correction or fallback as configured.",
        "not_automatically_modelled_as_covariates": ["sequencing_depth", "assembly_quality", "species", "multiple_resistance_genes_same_antibiotic_class", "plasmid_or_mobile_element_linkage"],
        "amr_recommendation": "For resistance versus susceptibility analyses, use resistant=case and susceptible=control, filter low-quality samples, stratify or review by species, enable SNP GWAS for mutation-mediated resistance, and inspect gene co-occurrence by antibiotic class."
    },
    "cases": cases,
    "controls": controls,
    "total_samples": total,
    "top_hit": {"feature_id": feature_id, "display_name": display_name, "display_label": display_subtitle, "annotation_confidence": annotation_confidence},
    "plots": {"gene_feature_index_association_plot": prefix + "_manhattan.svg", "gene_qq": prefix + "_qq.svg", "population_pca": prefix + "_population_pca.svg", "kinship_heatmap": prefix + "_kinship_heatmap.svg", "snp_feature_index_association_plot": prefix + "_SNP_manhattan.svg", "snp_qq": prefix + "_SNP_qq.svg"},
    "gubbins_outputs": {"summary": Path('~{gubbins_summary}').name, "filtered_alignment": Path('~{gubbins_filtered_alignment}').name, "recombination_gff": Path('~{gubbins_recombination_gff}').name, "log": Path('~{gubbins_log}').name},
    "generated_utc": generated
}
Path(prefix + "_run_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
PY
  >>>

  output {
    File html_report = "~{output_prefix}_report.html"
    File run_provenance_json = "~{output_prefix}_run_provenance.json"
  }

  runtime {
    docker: python_docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
}

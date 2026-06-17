# rMAP-GWAS

**rMAP-GWAS** is a Dockerized WDL/Cromwell workflow for reproducible **microbial genome-wide association studies (GWAS)** from paired-end bacterial isolate reads.

The current workflow accepts one sample set with explicit binary grouping, runs read trimming, genome assembly, annotation, pangenome construction, population-structure estimation, gene presence/absence GWAS, optional reference-based SNP GWAS, optional Gubbins recombination assessment for SNP analyses, post-GWAS GenBank annotation rescue & an integrated portable HTML report.

---

## Current implementation status

The current `rMAP_GWAS.wdl` implements:

- **Sample-set input validation** using `sample_names`, `read1s`, `read2s`, & `groups`.
- **Binary phenotype coding** from `groups` where cases are coded as `1` & controls as `0`.
- Optional biological display labels using `phenotype_name` & `phenotype_display_values` so the report can show contrasts such as `case (blood)` versus `control (sputum)`.
- **fastp** read trimming.
- **Shovill** bacterial isolate assembly.
- **QUAST** assembly QC.
- **Prokka** annotation by default, with optional **Bakta** annotation using `use_bakta=true`.
- **Panaroo** pangenome construction & gene presence/absence matrix generation.
- **Mash** distance matrix generation for population-structure correction & visualization.
- **pyseer** gene presence/absence GWAS.
- Prioritized gene association tables with case/control enrichment summaries & odds ratios.
- Reference GenBank annotation rescue for prioritized pangenome markers.
- Optional **Snippy + pyseer** reference-based SNP GWAS using `do_snp_gwas=true`.
- Optional **Gubbins** recombination assessment/filtering for SNP GWAS using `do_gubbins=true`.
- SVG Manhattan-style & QQ plots for gene & SNP GWAS results.
- A downloadable, offline HTML report with interpretation guardrails & key output files.

The current WDL **does not** run unitig/k-mer GWAS, Scoary, FastQC, MultiQC, or snpEff/bcftools csq. Those should not be described as active workflow steps unless they are added to the WDL.

## Visual summary of the rMAP-GWAS workflow
<p align="center">
  <img src="docs/assets/workflow/rMAP_GWAS.png"
       alt="rMAP-GWAS workflow"
       width="100%">
</p>

## Detailed workflow logic & sample routing

```text
Sample set / paired-end FASTQ inputs
        |
        |-- sample_names
        |-- read1s
        |-- read2s
        |-- groups
        |       case / control
        |
        v
Metadata-only input validation & cohort checks
        |
        |-- Confirm equal array lengths:
        |       sample_names
        |       read1_count
        |       read2_count
        |       groups
        |
        |-- Confirm unique sample IDs
        |
        |-- Confirm valid group labels:
        |       case / control
        |       1 / 0
        |       true / false
        |       yes / no
        |
        |-- Confirm at least one case and one control
        |
        |-- Generate validation report
        |
        |-- Important:
        |       This validation task does not localize FASTQ files.
        |       It checks metadata & array counts only.
        |
        `-- Downstream tasks start only after validation succeeds.
                         ║
                         ║
                         ▼
Reference package extraction
        |
        |-- Input:
        |       reference_docker
        |       reference_name
        |       reference_species
        |
        |-- Extract species-specific reference files from Docker image:
        |       reference.fasta
        |       reference.gff
        |       reference.genbank
        |
        |-- Used for:
        |       run provenance
        |       optional SNP GWAS
        |       post-GWAS GenBank annotation rescue
                         ║
                         ║
                         ▼
Phenotype table generation
        |
        |-- Case samples coded as:      1
        |-- Control samples coded as:   0
        |
        |-- Outputs:
        |       <output_prefix>_phenotypes.tsv
        |       <output_prefix>_sample_groups.tsv
                         ║
                         ║
                         ▼
All selected samples routed together
        |
        |-- No case/control separation for read processing
        |-- No species-exclusion branch at this stage
        |-- All samples proceed through the same QC, assembly, annotation,
        |   pangenome &  distance-matrix workflow
        |
        `-- Case/control labels are retained for phenotype coding,
            enrichment summaries, plots & GWAS interpretation.
                         ║
                         ║
                         ▼
Per-sample processing scatter
        |
        |-- For each selected sample:
        |
        |       Paired-end FASTQ files
        |              |
        |              v
        |       Read trimming & QC with fastp
        |              |
        |              |-- Trimmed R1 FASTQ
        |              |-- Trimmed R2 FASTQ
        |              |-- fastp HTML report
        |              `-- fastp JSON report
        |              |
        |              v
        |       *De novo* genome assembly with Shovill
        |              |
        |              |-- Contigs FASTA
        |              `-- Shovill log
        |              |
        |              v
        |       Assembly quality control with QUAST
        |              |
        |              `-- QUAST report TSV
        |              |
        |              v
        |       Genome annotation with Prokka
        |              |
        |              |-- GFF annotation
        |              |-- GenBank file
        |              |-- Protein FASTA
        |              `-- Nucleotide feature FASTA
        |
        `-- End per-sample scatter
                         ║
                         ║
                         ▼
Cohort-level pangenome construction
        |
        |-- Collect all Prokka GFF files
        |-- Run Panaroo pangenome analysis
        |-- Clean pangenome graph
        |-- Remove invalid genes where supported
        |-- Generate gene presence/absence matrix
        |
        |-- Outputs:
        |       gene_presence_absence.csv
        |       gene_presence_absence.Rtab
        |       gene_data.csv
        |       combined_DNA_CDS.fasta
        |       combined_protein_CDS.fasta
        |       pan_genome_reference.fa
        |       panaroo_summary.txt
                         ║
                         ║
                         ▼
Cohort-level genome distance estimation
        |
        |-- Collect all Shovill assemblies
        |-- Rename assemblies by sample ID
        |-- Create Mash sketches
        |-- Compute pairwise Mash distances
        |-- Convert long Mash output to a square pyseer-compatible matrix
        |
        |-- Output:
        |       mash_distances.tsv
                         ║
                         ║
                         ▼
Population-structure visualization
        |
        |-- Input:
        |       phenotype table
        |       Mash distance matrix
        |
        |-- Generate:
        |       Mash distance / kinship heatmap
        |       PCoA population-structure plot
        |
        `-- Used to assess whether phenotype labels cluster by lineage
            or genetic background.
                         ║
                         ║
                         ▼
Gene presence/absence GWAS
        |
        |-- Inputs:
        |       phenotype table
        |       Panaroo gene presence/absence matrix
        |       Mash distance matrix
        |
        |-- Run pyseer gene presence/absence GWAS
        |
        |-- Apply allele-frequency filters:
        |       min_af
        |       max_af
        |
        |-- Correct for population structure using Mash distances
        |       pyseer_max_dimensions
        |       pyseer_force_no_distances
        |       pyseer_no_distances_fallback
        |
        |-- Output:
        |       pyseer_gene_assoc.tsv
                         ║
                         ║
                         ▼
Gene GWAS plotting
        |
        |-- Generate:
        |       <output_prefix>_gene_gwas_qq.svg
        |       <output_prefix>_gene_gwas_manhattan.svg
        |       <output_prefix>_gene_gwas_plot_summary.tsv
        |
        |-- Note:
        |       The gene GWAS association plot is feature-index based.
        |       It is not a true reference-coordinate Manhattan plot unless
        |       gene clusters are confidently mapped to reference coordinates.
                         ║
                         ║
                         ▼
GWAS hit prioritization
        |
        |-- Parse pyseer association results
        |-- Link association hits to Panaroo/Prokka annotations
        |-- Calculate case and control gene frequencies
        |-- Estimate enrichment direction:
        |       case-enriched
        |       control-enriched
        |       mixed / check manually
        |
        |-- Rank features using:
        |       p-value or q-value
        |       effect direction
        |       odds ratio
        |       annotation availability
        |       recurrence across samples
        |
        |-- Outputs:
        |       <output_prefix>_all_ranked_hits.tsv
        |       <output_prefix>_top_priority_hits.tsv
        |       <output_prefix>_all_significant_hits.tsv
        |       <output_prefix>_enrichment_summary.tsv
                         ║
                         ║
                         ▼
Post-GWAS GenBank annotation rescue
        |
        |-- Inputs:
        |       prioritized GWAS hit tables
        |       Panaroo gene_presence_absence.csv
        |       Panaroo gene_data.csv
        |       Panaroo combined DNA/protein CDS files
        |       Panaroo pan-genome reference
        |       species-specific reference GenBank
        |
        |-- For each prioritized Panaroo cluster:
        |       identify representative sequence
        |       compare against reference GenBank CDS features
        |       rescue gene/locus/product annotation where possible
        |
        |-- Add annotation columns:
        |       reference_locus_tag
        |       reference_gene
        |       reference_product
        |       reference_location
        |       reference_match_type
        |       reference_identity
        |       reference_coverage
        |       annotation_confidence
        |       annotation_note
        |
        |-- Confidence interpretation:
        |       high    = strong reference-supported annotation
        |       medium  = plausible annotation; inspect manually
        |       low     = tentative annotation only
        |       none    = no usable GenBank match
        |
        |-- Outputs:
        |       <output_prefix>_top_priority_hits.annotated.tsv
        |       <output_prefix>_all_significant_hits.annotated.tsv
        |       <output_prefix>_reference_annotation_summary.tsv
                         ║
                         ║
                         ▼
Optional SNP GWAS branch
        |
        |-- Controlled by:
        |       do_snp_gwas = true / false
        |       do_gubbins = true / false
        |
        |-- If do_snp_gwas = false:
        |       create SNP placeholder files
        |       create Gubbins placeholder files
        |       report that SNP GWAS & Gubbins were not run
        |
        |-- If do_snp_gwas = true:
        |       use reference FASTA from reference_docker
        |       call SNPs with Snippy
        |       generate Snippy SNP VCF & core alignment
        |
        |       Optional Gubbins recombination assessment / filtering
        |              |
        |              |-- Controlled by:
        |              |       do_gubbins = true / false
        |              |
        |              |-- If do_gubbins = false:
        |              |       keep original Snippy SNP VCF for pyseer
        |              |       create Gubbins placeholder summary/log files
        |              |
        |              |-- If do_gubbins = true:
        |              |       run Gubbins on the Snippy core alignment
        |              |       identify predicted recombinant regions
        |              |       remove SNPs inside predicted recombinant blocks
        |              |       create a Gubbins-filtered SNP VCF
        |              |
        |              |-- Safe fallback behavior:
        |              |       if the core alignment is empty
        |              |       if fewer than 3 sequences are available
        |              |       if run_gubbins.py is unavailable
        |              |       if Gubbins fails softly
        |              |       if filtering removes all SNP records
        |              |
        |              `-- In fallback cases, continue with the original Snippy SNP VCF
        |                  and record the reason in the Gubbins summary/log.
        |
        |       run pyseer SNP GWAS using:
        |              Gubbins-filtered VCF when available
        |              otherwise original Snippy SNP VCF
        |
        |       prioritize SNP associations
        |       generate SNP QQ & Manhattan-style plots
        |
        |-- Outputs when enabled:
        |       <output_prefix>_SNP.vcf
        |       <output_prefix>_SNP_pyseer_assoc.tsv
        |       <output_prefix>_SNP_top_hits.tsv
        |       <output_prefix>_SNP_all_significant_hits.tsv
        |       <output_prefix>_SNP_summary.tsv
        |       <output_prefix>_SNP_qq.svg
        |       <output_prefix>_SNP_manhattan.svg
        |       <output_prefix>_SNP_gubbins_summary.tsv
        |       <output_prefix>_SNP_gubbins.filtered_polymorphic_sites.fasta
        |       <output_prefix>_SNP_gubbins.recombination_predictions.gff
        |       <output_prefix>_SNP_gubbins.filtered.snps.vcf
        |       <output_prefix>_SNP_gubbins.log
                         ║
                         ║
                         ▼
Integrated reporting & provenance
        |
        |-- Merge:
        |       validation report
        |       phenotype table
        |       sample group table
        |       Panaroo summary
        |       Mash population-structure plots
        |       gene GWAS QQ plot
        |       gene GWAS feature-index association plot
        |       prioritized annotated gene hits
        |       optional SNP GWAS outputs
        |       optional Gubbins recombination summary/log/filtering outputs
        |       reference annotation summary
        |
        |-- Display:
        |       case/control composition
        |       workflow architecture
        |       top-hit GenBank annotation rescue
        |       top 5 prioritized GWAS hits
        |       all significant hits
        |       annotation confidence guide
        |       optional Gubbins status and recommendation
        |       Gubbins input-vs-filtered SNP record counts when available
        |       small-sample-size caution
        |       Panaroo cluster-ID caveat
        |       PE/PPE & repetitive-region caution where relevant
        |
        |-- Record provenance:
        |       reference_docker
        |       reference_species
        |       reference_name
        |       container backend
        |       GWAS mode
        |       SNP GWAS status
        |       do_gubbins setting
        |       Gubbins run status and filtering note
        |       Gubbins output filenames
        |       pyseer distance-correction settings
        |
        |-- Final outputs:
        |       <output_prefix>_report.html
        |       <output_prefix>_run_provenance.json
        |       <output_prefix>_SNP_gubbins_summary.tsv
        |       <output_prefix>_SNP_gubbins.log
```


## Input model

The workflow uses a **single sample-set design**, not separate case/control arrays.

### Required inputs

| Input | Type | Description |
|---|---:|---|
| `sample_names` | `Array[String]+` | Unique sample IDs. Sample names must not contain whitespace. |
| `read1s` | `Array[File]+` | Forward paired-end FASTQ files. |
| `read2s` | `Array[File]+` | Reverse paired-end FASTQ files. |
| `groups` | `Array[String]+` | Binary group labels. Accepted case labels include `case`, `cases`, `1`, `true`, `yes`. Accepted control labels include `control`, `controls`, `0`, `false`, `no`. |

The four required arrays must have equal length & be in the same sample order.

### Optional phenotype display inputs

| Input | Type | Default | Description |
|---|---:|---:|---|
| `phenotype_name` | `String` | `case_control` | Name shown in phenotype & report tables. Use the biological contrast column name when available, for example `specimen_source`, `drug_resistance`, or `carriage_status`. |
| `phenotype_display_values` | `Array[String]` | `[]` | Optional metadata values used only to make the HTML report biologically interpretable. If provided, this array must have the same length & order as `sample_names`. |

Example: if `groups` is `case/control` but `phenotype_display_values` is `blood/sputum`, the report can display the contrast as **case (blood)** versus **control (sputum)** while preserving binary GWAS coding.

---

## Minimal input JSON

```json
{
  "rMAP_GWAS.sample_names": ["sample_001", "sample_002", "sample_003", "sample_004"],
  "rMAP_GWAS.read1s": [
    "~/sample_001_R1.fastq.gz",
    "~/sample_002_R1.fastq.gz",
    "~/sample_003_R1.fastq.gz",
    "~/sample_004_R1.fastq.gz"
  ],
  "rMAP_GWAS.read2s": [
    "~/sample_001_R2.fastq.gz",
    "~/sample_002_R2.fastq.gz",
    "~/sample_003_R2.fastq.gz",
    "~/sample_004_R2.fastq.gz"
  ],
  "rMAP_GWAS.groups": ["case", "case", "control", "control"]
}
```

### Example with biological display labels

```json
{
  "rMAP_GWAS.sample_names": ["sample_001", "sample_002", "sample_003", "sample_004"],
  "rMAP_GWAS.read1s": [
    "~/sample_001_R1.fastq.gz",
    "~/sample_002_R1.fastq.gz",
    "~/sample_003_R1.fastq.gz",
    "~/sample_004_R1.fastq.gz"
  ],
  "rMAP_GWAS.read2s": [
    "~/sample_001_R2.fastq.gz",
    "~/sample_002_R2.fastq.gz",
    "~/sample_003_R2.fastq.gz",
    "~/sample_004_R2.fastq.gz"
  ],
  "rMAP_GWAS.groups": ["case", "case", "control", "control"],
  "rMAP_GWAS.phenotype_name": "specimen_source",
  "rMAP_GWAS.phenotype_display_values": ["blood", "blood", "sputum", "sputum"]
}
```

### Terra sample-set mapping

When running on a Terra Cloud platform sample set, map the workflow inputs to the entity-set member attributes. For example, if the entity set is `gwasmtb_set` & the member entity type is `gwasmtb`:

```text
sample_names              = this.gwasmtbs.gwasmtb_id
read1s                    = this.gwasmtbs.read1
read2s                    = this.gwasmtbs.read2
groups                    = this.gwasmtbs.group
phenotype_display_values  = this.gwasmtbs.specimen_source
phenotype_name            = "specimen_source"
```

The `groups` column should contain the binary GWAS labels. The `phenotype_display_values` column should contain the biological label that makes the report readable.

---

## Common analysis controls

| Input | Type | Default | Description |
|---|---:|---:|---|
| `min_af` | `Float` | `0.01` | Minimum allele/feature frequency passed to pyseer. |
| `max_af` | `Float` | `0.99` | Maximum allele/feature frequency passed to pyseer. |
| `significance_alpha` | `Float` | `0.05` | Significance threshold used for prioritized hit tables & plot summaries. |
| `pyseer_max_dimensions` | `Int` | `2` | Maximum MDS dimensions used with Mash distances. Small tests may need fewer dimensions. |
| `pyseer_force_no_distances` | `Boolean` | `false` | If `true`, pyseer runs without Mash distance correction. |
| `pyseer_no_distances_fallback` | `Boolean` | `true` | If the distance-corrected pyseer model fails, retry without distances. |
| `output_prefix` | `String` | `rMAP_GWAS` | Prefix used for report & output file names. |
| `plot_max_points` | `Int` | `5000` | Maximum number of points drawn in SVG GWAS plots. |

`gwas_mode` is recorded in the report/provenance. In the current WDL, the gene presence/absence GWAS branch always runs; the SNP branch is controlled by `do_snp_gwas`.

---

## Optional SNP GWAS & Gubbins controls

| Input | Type | Default | Description |
|---|---:|---:|---|
| `do_snp_gwas` | `Boolean` | `false` | Enables reference-based SNP calling with Snippy & SNP association testing with pyseer. |
| `do_gubbins` | `Boolean` | `false` | Enables optional Gubbins recombination assessment/filtering for the SNP branch. This is only used when `do_snp_gwas=true`. |
| `snp_min_qual` | `Float` | `20.0` | Minimum SNP quality retained from the Snippy/core VCF. |
| `container_backend` | `String` | `docker` | Backend label recorded in the report/provenance. |

Use SNP GWAS for mutation-mediated phenotypes, such as many MTBC drug-resistance traits. For recombining bacteria, enable Gubbins when recombinant blocks may inflate SNP associations. For highly clonal organisms or test runs, keeping `do_gubbins=false` is usually simpler.

---

## Annotation engine

The default annotation engine is **Prokka**. To use Bakta, set:

```json
{
  "rMAP_GWAS.use_bakta": true,
  "rMAP_GWAS.bakta_docker": "gmboowa/rmap-gwas-bakta-db:light-0.1",
  "rMAP_GWAS.bakta_db": "/opt/bakta/db-light"
}
```

Bakta mode requires a compatible Bakta database available inside the selected container or mounted at the path given by `bakta_db`.

---

## Reference package inputs

The workflow always extracts reference files from a species-specific reference Docker image. A non-empty GenBank file is required for reference annotation rescue. A reference FASTA is additionally required when `do_snp_gwas=true`.

| Input | Type | Default | Description |
|---|---:|---:|---|
| `reference_docker` | `String` | `gmboowa/rmap-gwas-mtbc-refs:2026.06` | Docker image containing the species reference package. |
| `reference_species` | `String` | `*Mycobacterium tuberculosis*` | Species label recorded in the report. |
| `reference_name` | `String` | `MTBC_2026_06` | Reference package name recorded in the report. |

The WDL searches common paths such as:

```text
/opt/rmap-gwas/refs/reference.fasta
/opt/rmap-gwas/refs/reference.gff
/opt/rmap-gwas/refs/reference.genbank
/opt/rmap-gwas/refs/mtbc/reference.fasta
/opt/rmap-gwas/refs/mtbc/reference.gff
/opt/rmap-gwas/refs/mtbc/reference.genbank
/opt/rmap-gwas/refs/kpneumo/reference.fasta
/opt/rmap-gwas/refs/kpneumo/reference.gff
/opt/rmap-gwas/refs/kpneumo/reference.genbank
/refs/reference.fasta
/refs/reference.gff
/refs/reference.genbank
/data/reference.fasta
/data/reference.gff
/data/reference.genbank
```

The reference image may also define these environment variables:

```text
RMAP_GWAS_REFERENCE_FASTA
RMAP_GWAS_REFERENCE_GFF
RMAP_GWAS_REFERENCE_GENBANK
```

---

## Runtime controls

The WDL exposes task-level CPU, memory & disk settings. Defaults are test friendly & may need to be increased for larger cohorts.

| Resource group | Default threads | Default memory | Default disk |
|---|---:|---:|---:|
| fastp | `4` | `8 GB` | `50 GB` |
| assembly | `4` | `32 GB` | `200 GB` |
| annotation | `4` | `16 GB` | `100 GB` |
| pangenome | `4` | `32 GB` | `300 GB` |
| GWAS | `4` | `32 GB` | `300 GB` |
| SNP GWAS | uses `gwas_threads` | uses `gwas_memory_gb` | `300 GB` |
| Gubbins | `4` | `16 GB` | `300 GB` |

For full microbial GWAS cohorts, Panaroo, Mash, pyseer, SNP calling & Gubbins are the most likely steps to require additional memory, disk, or runtime.

---

## Default Docker images

| Step | Default image |
|---|---|
| fastp | `quay.io/biocontainers/fastp:0.23.4--hadf994f_2` |
| Shovill | `quay.io/biocontainers/shovill:1.1.0--hdfd78af_1` |
| QUAST | `staphb/quast:5.2.0` |
| Prokka | `staphb/prokka:1.14.6` |
| Bakta | `gmboowa/rmap-gwas-bakta-db:light-0.1` |
| Panaroo | `quay.io/biocontainers/panaroo:1.5.2--pyhdfd78af_0` |
| Mash / pyseer / Python utilities | `gmboowa/rmap-gwas-pyseer-annotate:0.2` |
| Snippy | `staphb/snippy:4.6.0` |
| Gubbins | `staphb/gubbins:latest` |
| Reference package | `gmboowa/rmap-gwas-mtbc-refs:2026.06` |

---

## Running locally with Cromwell

Validate the WDL if you have `womtool` available:

```bash
java -jar womtool.jar validate rMAP_GWAS.wdl
```

Run with Cromwell:

```bash
java -jar cromwell.jar run rMAP_GWAS.wdl -i inputs.json
```

For local Cromwell runs, make sure Docker is running & that your input FASTQ paths are accessible to the execution environment.

---

## Main workflow outputs

| Output | Description |
|---|---|
| `phenotype_tsv` | Binary phenotype table used by pyseer. |
| `phenotype_legend_tsv` | Case/control display legend used by the report. |
| `sample_groups_tsv` | Sample-level group and display-label table. |
| `trimmed_read1s`, `trimmed_read2s` | fastp-trimmed paired-end reads. |
| `assemblies` | Shovill contig FASTA files. |
| `quast_reports` | Per-sample QUAST report TSV files. |
| `gffs`, `gbks`, `faas`, `ffns` | Per-sample annotation outputs from Prokka or Bakta. |
| `gene_presence_absence_csv` | Panaroo pangenome gene presence/absence CSV. |
| `gene_presence_absence_rtab` | pyseer-compatible gene presence/absence matrix. |
| `mash_distances` | Square Mash distance matrix used by pyseer. |
| `population_pca_svg` | Mash-based PCoA plot. |
| `kinship_heatmap_svg` | Mash distance/kinship heatmap. |
| `population_structure_summary` | Population-structure summary table. |
| `pyseer_gene_assoc` | Raw pyseer gene presence/absence association results. |
| `raw_top_priority_hits` | Unannotated prioritized gene hits. |
| `raw_all_significant_hits` | Unannotated significant gene hits. |
| `top_priority_hits` | GenBank-annotated top-priority gene hits. |
| `all_significant_hits` | GenBank-annotated significant gene hits. |
| `reference_annotation_summary` | Summary of GenBank annotation rescue. |
| `enrichment_summary` | Case/control enrichment summary for gene hits. |
| `qq_plot_svg` | Gene GWAS QQ plot. |
| `manhattan_plot_svg` | Gene presence/absence Manhattan-style plot. |
| `plot_summary` | Gene GWAS plot summary. |
| `reference_fasta`, `reference_gff`, `reference_genbank` | Reference files extracted from the reference Docker image. |
| `snp_vcf` | SNP VCF from the SNP branch, or a placeholder if SNP GWAS was not run. |
| `pyseer_snp_assoc` | Raw pyseer SNP association results, or a placeholder if SNP GWAS was not run. |
| `snp_top_hits` | Prioritized SNP hits, or a placeholder if SNP GWAS was not run. |
| `snp_all_significant_hits` | All significant SNP hits, or a placeholder if SNP GWAS was not run. |
| `snp_summary` | SNP GWAS summary. |
| `snp_manhattan_plot_svg` | SNP GWAS Manhattan plot, or a not-run placeholder. |
| `snp_qq_plot_svg` | SNP GWAS QQ plot, or a not-run placeholder. |
| `snp_plot_summary` | SNP plot summary. |
| `gubbins_summary` | Gubbins status & filtering summary. |
| `gubbins_filtered_alignment` | Gubbins filtered alignment, or a placeholder if not run. |
| `gubbins_recombination_gff` | Gubbins recombination GFF, or a placeholder if not run. |
| `gubbins_log` | Gubbins log, or a not-run log. |
| `gubbins_filtered_vcf` | Gubbins-filtered SNP VCF, or a placeholder if not run. |
| `html_report` | Integrated offline HTML report. |
| `run_provenance` | JSON provenance file recording workflow configuration & selected reference/container settings. |

---

## HTML report contents

The generated report is designed for collaborators who need interpretable results rather than raw workflow logs. Current report sections include:

1. **Phenotype legend & coding**
2. **AMR phenotype and association scope**
3. **Workflow architecture**
4. **Run configuration**
5. **Top-hit GenBank annotation rescue**
6. **Population structure**
7. **Gene presence/absence GWAS plots**
8. **SNP marker GWAS**
9. **Optional Gubbins recombination assessment**
10. **Top priority gene presence/absence GWAS hits**
11. **All significant gene presence/absence hits**
12. **Reference annotation confidence guide**
13. **Input validation**
14. **Panaroo summary**
15. **Interpretation guidance**
16. **Key output files**

The key output file cards in the HTML report are downloadable from within the report.

---

## Interpretation guidance

rMAP-GWAS reports **statistical associations**, not proven causal mechanisms or clinical diagnostic calls.

Important interpretation checks include:

- Confirm that `groups` encodes the intended binary phenotype correctly.
- Confirm that `phenotype_display_values` matches the biological contrast shown in the report.
- Review case/control balance before interpreting associations.
- Review Mash PCoA & kinship plots for lineage or outbreak clustering.
- For AMR phenotypes, distinguish phenotypic resistance/susceptibility from gene carriage.
- For mutation-mediated resistance, consider enabling SNP GWAS.
- For recombining bacteria, consider enabling Gubbins for SNP analyses.
- Treat low-confidence GenBank rescue annotations as tentative.
- Validate candidate markers in independent datasets before biological or public-health interpretation.

The workflow warns when case or control counts are small. Very small test runs are useful for testing execution & report generation but should not be treated as final GWAS evidence.

---

## Recommended cohort metadata

A practical Terra metadata table should include at least:

| Column | Example | Purpose |
|---|---|---|
| sample ID | `sample_001` | Entity ID / sample name. |
| read1 | `~/.../sample_001_R1.fastq.gz` | Forward read. |
| read2 | `~/.../sample_001_R2.fastq.gz` | Reverse read. |
| group | `case` or `control` | Binary GWAS class. |
| phenotype display column | `blood`, `sputum`, `resistant`, `susceptible` | Biological label for the report. |
| species or lineage | `MTBC`, `*K. pneumoniae*`, lineage/ST | Helps interpret population structure. |
| source country/site/date | optional | Helps detect sampling imbalance & outbreak clusters. |

Only the first four fields are required by the WDL. Additional metadata should be preserved for interpretation & reporting.

---

## Troubleshooting

### Input validation fails

Check that:

- `sample_names`, `read1s`, `read2s` &  `groups` have the same length.
- Sample names are unique.
- Sample names do not contain whitespace.
- `groups` contains recognizable binary labels.
- `phenotype_display_values`, if provided, has the same length as `sample_names`.

### pyseer fails with distances

For small or highly imbalanced test cohorts, the pyseer null model can fail with distance correction. The default `pyseer_no_distances_fallback=true` retries with `--no-distances` so execution can complete. For final analyses, inspect population structure carefully & consider larger, better-balanced cohorts.

### SNP GWAS fails

SNP GWAS requires a non-empty reference FASTA extracted from `reference_docker`. Confirm that the reference image contains a valid FASTA & GenBank file.

### Bakta fails

Confirm that `use_bakta=true` is intentional & that `bakta_db` points to a database path available inside the selected Bakta container.

### Gubbins is not run

Gubbins only runs when both of these are true:

```json
{
  "rMAP_GWAS.do_snp_gwas": true,
  "rMAP_GWAS.do_gubbins": true
}
```

If SNP GWAS is disabled, the workflow produces Gubbins placeholder outputs for report completeness.

---

## Repository layout suggestion

```text
rMAP-GWAS/
├── rMAP_GWAS.wdl
├── README.md
├── inputs/
│   └── example.inputs.json
├── docs/
│   └── index.html
└── test-data/
    └── README.md
```

---

## Citation

If you use rMAP-GWAS, cite this repository & the underlying tools used in the workflow, including fastp, Shovill, QUAST, Prokka or Bakta, Panaroo, Mash, pyseer, Snippy & Gubbins where applicable.


---

## License

Add the project license here, for example MIT, Apache-2.0, or another license selected by the repository owner.

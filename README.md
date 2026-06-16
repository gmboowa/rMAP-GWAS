# rMAP-GWAS

**rMAP-GWAS** — **rapid Microbial Analysis Pipeline for Genome-Wide Association Studies** — is a portable WDL/Cromwell workflow for microbial case-control genome-wide association studies from paired-end Illumina reads. It is designed to generate interpretable, reproducible association results with annotated top-priority loci, case/control enrichment, statistical evidence, plots & a self-contained HTML report.

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
        |       It checks metadata and array counts only.
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
        |
        |-- If do_snp_gwas = false:
        |       create SNP placeholder files
        |       report that SNP GWAS was not run
        |
        |-- If do_snp_gwas = true:
        |       use reference FASTA from reference_docker
        |       call SNPs with Snippy
        |       run pyseer SNP GWAS
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
        |       reference annotation summary
        |
        |-- Display:
        |       case/control composition
        |       workflow architecture
        |       top-hit GenBank annotation rescue
        |       top 5 prioritized GWAS hits
        |       all significant hits
        |       annotation confidence guide
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
        |       pyseer distance-correction settings
        |
        |-- Final outputs:
        |       <output_prefix>_report.html
        |       <output_prefix>_run_provenance.json
```

### Sample routing logic

rMAP-GWAS uses a cohort-level case/control design. All selected samples are routed together through read QC, assembly, annotation, pangenome construction & distance estimation. The `groups` input is used to create the GWAS phenotype table, not to split samples into separate processing branches.

Each sample must have:

```text
sample ID
read1 FASTQ
read2 FASTQ
group label
```

The default group labels are:

```text
case     -> phenotype value 1
control  -> phenotype value 0
```

The workflow also accepts common equivalents such as `1/0`, `true/false`, & `yes/no`, depending on the configured `case_label` & `control_label`.

### Validation-ordering design

The current workflow validates metadata before expensive tasks are allowed to proceed. The validation task receives:

```text
sample_names
read1_count
read2_count
groups
case_label
control_label
```

It does not receive the actual FASTQ files. This prevents unnecessary FASTQ localization during validation & avoids wasting compute if the input table has incorrect sample IDs, group labels, or array lengths.

### Reference configuration

The same WDL can be used across species by changing:

```text
reference_docker
reference_name
reference_species
```

Each species-specific reference Docker image should provide:

```text
reference.fasta
reference.gff
reference.genbank
```

The GenBank file is used for post-GWAS annotation rescue of prioritized Panaroo gene clusters.

Example MTBC configuration:

```text
reference_docker  = "gmboowa/rmap-gwas-mtbc-refs:2026.06"
reference_name    = "MTBC_2026_06"
reference_species = "Mycobacterium tuberculosis complex"
```

Example *Klebsiella pneumoniae* configuration:

```text
reference_docker  = "gmboowa/rmap-gwas-kpneumo-refs:2026.06"
reference_name    = "KPNEUMO_2026_06"
reference_species = "Klebsiella pneumoniae"

```

The workflow also accepts common equivalents such as `1/0`, `true/false`, & `yes/no`.


### Case–control interpretation

The phenotype table is generated internally from the selected sample set. Samples labelled as cases are coded as `1`, while controls are coded as `0`.

```text
sample          case_control
SRRxxxxxx       1
SRRyyyyyy       0
```

This phenotype table is passed to pyseer together with the Panaroo gene presence/absence matrix & Mash distance matrix.

### Interpretation notes

Panaroo cluster IDs such as `group_2270` are pangenome feature identifiers, not stable biological gene names. These IDs may change across runs depending on the input cohort & pangenome clustering.

The GenBank annotation rescue step improves interpretability by mapping prioritized clusters back to a species-specific reference where possible. Low-confidence matches should be treated as tentative & should retain the original Panaroo cluster ID in reports.

For small smoke-test runs, the workflow can validate execution, reporting & end-to-end integration, but association results should not be interpreted as final biological or clinical findings without larger cohorts & independent validation.


### Population structure handling

Microbial GWAS can be confounded by clonal population structure, relatedness, outbreak clusters, lineage effects & uneven case/control sampling. rMAP-GWAS therefore computes pairwise Mash distances from assembled genomes & supplies this distance matrix to pyseer during gene-based association testing.

```text
Assemblies
    ↓
Mash sketch
    ↓
Mash distance matrix
    ↓
pyseer population-structure-aware association testing
```

### Reference & provenance tracking

For MTBC runs, the workflow records the intended reference package:

```text
`reference_docker`  = `gmboowa/rmap-gwas-mtbc-refs:2026.06`  
`reference_species` = *Mycobacterium tuberculosis* complex  
`reference_name`    = `MTBC_2026_06`
```


These values are written into the final report and provenance JSON. They document the species/reference configuration used or intended for the analysis.

### Current annotation mode

This workflow version uses Prokka annotation before Panaroo pangenome construction.

```text
Shovill contigs
    ↓
Prokka annotation
    ↓
GFF files
    ↓
Panaroo pangenome
```

A Bakta-based version can be implemented separately by replacing the annotation task & ensuring downstream GFF compatibility with Panaroo.

### Final workflow outputs

The major outputs are:

```text
Phenotype table
Trimmed FASTQ files
Genome assemblies
QUAST assembly QC reports
Prokka annotation files
Panaroo pangenome matrices
Mash distance matrix
pyseer gene association results
Prioritized GWAS hit tables
Enrichment summary
Interactive HTML report
Run provenance JSON
```

### Interpretation note

rMAP-GWAS identifies statistical associations between gene presence/absence patterns & case–control phenotype labels. These associations should be interpreted carefully alongside sample size, population structure, lineage distribution, outbreak clustering, phenotype definition & biological plausibility. Candidate hits should be validated in independent datasets or by targeted experimental/epidemiological follow-up where possible.


Microbial genome-wide association studies (mGWAS) can identify bacterial genetic variants, unitigs, genes & accessory-genome features associated with phenotypes such as antimicrobial resistance, virulence, host source, outbreak status, colonization, infection, or clinical outcome. However, microbial GWAS requires careful handling of population structure, case-control imbalance, phenotype quality, feature annotation & reproducible reporting.

`rMAP-GWAS` aims to provide an end-to-end cloud-ready workflow that starts from clearly designated **case** & **control** paired-end reads & produces:


- read QC & trimming summaries
- de novo assemblies
- assembly quality-control metrics
- genome annotation
- pangenome/gene presence-absence matrix
- unitig/k-mer feature matrix
- population-structure correction inputs
- microbial GWAS results using `pyseer`
- annotated significant loci & genes
- case/control enrichment summaries
- Manhattan & QQ plots
- top-priority GWAS hits table
- portable offline HTML report

---

## Intended use cases

`rMAP-GWAS` is designed for microbial isolate datasets where samples can be assigned into two phenotype groups:

- resistant vs susceptible isolates
- invasive vs colonizing isolates
- outbreak vs non-outbreak isolates
- hypervirulent vs non-hypervirulent isolates
- carbapenemase-positive vs carbapenemase-negative isolates
- convergent MDR-hvKp vs non-convergent isolates
- case vs control definitions supplied by the user

The workflow is organism-agnostic in principle, but users should interpret results in the context of species biology, sampling structure, recombination, clonal expansion & phenotype quality.

---

## Key features

- **Case-control aware input design**  
  Users provide separate arrays of case & control paired-end reads in a Cromwell/Terra JSON file.

- **Microbial GWAS engine**  
  Primary association testing is performed using `pyseer`.

- **Feature types supported**
  - unitigs/k-mers
  - pangenome gene presence/absence
  - optional reference-based SNPs when a reference genome & annotation are provided

- **Population-structure correction**  
  The workflow generates distance/covariate inputs to reduce false-positive associations caused by clonal population structure.

- **Annotated top hits**  
  Significant hits are annotated with gene names, products, genomic context, case/control frequencies, enrichment direction, *p*-values, q-values & priority scores.

- **Portable reporting**  
  A final self-contained HTML report summarizes the full analysis & can be shared without requiring external JavaScript, notebooks, or web services.

---

## Workflow structure

```text
rMAP_GWAS
├── VALIDATE_CASE_CONTROL_INPUTS
├── PREPARE_PHENOTYPE_TABLE
├── FASTP_TRIMMING
├── FASTQC
├── MULTIQC
├── ASSEMBLE_GENOMES
├── ASSEMBLY_QC
├── ANNOTATE_GENOMES
├── PANAROO_PANGENOME
├── BUILD_GENE_MATRIX
├── UNITIG_CALLER
├── MASH_DISTANCE_MATRIX
├── PYSEER_UNITIG_GWAS
├── PYSEER_GENE_GWAS
├── optional PYSEER_SNP_GWAS
├── ANNOTATE_SIGNIFICANT_HITS
├── PRIORITIZE_GWAS_HITS
├── MAKE_GWAS_PLOTS
└── MERGE_RMAP_GWAS_REPORT
```

---

## Primary tools

| Stage | Recommended tool |
|---|---|
| Read QC | FastQC |
| Read trimming | fastp |
| QC aggregation | MultiQC |
| Assembly | Shovill or SPAdes |
| Assembly QC | QUAST |
| Genome annotation | Bakta or Prokka |
| Pangenome | Panaroo |
| Unitig generation | unitig-caller |
| Population distance | Mash |
| GWAS engine | pyseer |
| Optional gene-only scan | Scoary |
| Optional SNP workflow | Snippy, bcftools, snpEff/bcftools csq |
| Reporting | Python, pandas, matplotlib, Jinja2 |

---

## Docker/container strategy

The workflow uses public BioContainers, Staph-B, or maintained project images where possible & custom `gmboowa/rmap-gwas-*` images for integration/reporting tasks.

### Public images to use where possible

| Task | Example image |
|---|---|
| fastp | `quay.io/biocontainers/fastp:<tag>` |
| FastQC | `staphb/fastqc:<tag>` |
| MultiQC | `multiqc/multiqc:<tag>` |
| Shovill | `quay.io/biocontainers/shovill:<tag>` |
| SPAdes | `quay.io/biocontainers/spades:<tag>` |
| QUAST | `staphb/quast:<tag>` |
| Panaroo | `quay.io/biocontainers/panaroo:<tag>` |
| unitig-caller | `quay.io/biocontainers/unitig-caller:<tag>` |
| pyseer | `quay.io/biocontainers/pyseer:<tag>` |
| Mash | `quay.io/biocontainers/mash:<tag>` |
| Scoary | `quay.io/biocontainers/scoary:<tag>` |
| Snippy | `staphb/snippy:<tag>` |

### Custom images 

| Image | Purpose |
|---|---|
| `gmboowa/rmap-gwas-pyseer-annotate` | pyseer execution, hit annotation, enrichment summaries, and prioritization |
| `gmboowa/rmap-gwas-report` | final offline HTML report generation |
| `gmboowa/rmap-gwas-bakta-db` | optional Bakta database image or reference bundle |
| `gmboowa/rmap-gwas-reference-bundle:<species-tag>` | optional species-specific reference FASTA/GFF/GenBank resources |

---

## Input JSON example

```json
{
  "rMAP_GWAS.case_sample_names": ["case_001", "case_002"],
  "rMAP_GWAS.case_read1s": ["~/case_001_R1.fastq.gz", "~/case_002_R1.fastq.gz"],
  "rMAP_GWAS.case_read2s": ["~/case_001_R2.fastq.gz", "~/case_002_R2.fastq.gz"],

  "rMAP_GWAS.control_sample_names": ["control_001", "control_002"],
  "rMAP_GWAS.control_read1s": ["~/control_001_R1.fastq.gz", "~/control_002_R1.fastq.gz"],
  "rMAP_GWAS.control_read2s": ["~/control_001_R2.fastq.gz", "~/control_002_R2.fastq.gz"],

  "rMAP_GWAS.phenotype_name": "case_control",
  "rMAP_GWAS.case_label": "case",
  "rMAP_GWAS.control_label": "control",

  "rMAP_GWAS.do_trimming": true,
  "rMAP_GWAS.do_assembly": true,
  "rMAP_GWAS.do_annotation": true,
  "rMAP_GWAS.do_unitig_gwas": true,
  "rMAP_GWAS.do_gene_gwas": true,
  "rMAP_GWAS.do_snp_gwas": false,

  "rMAP_GWAS.reference_fasta": "~/reference.fasta",
  "rMAP_GWAS.reference_gff": "~/reference.gff",
  "rMAP_GWAS.covariates_tsv": "~/covariates.tsv",

  "rMAP_GWAS.min_af": 0.01,
  "rMAP_GWAS.max_af": 0.99,
  "rMAP_GWAS.min_cases": 10,
  "rMAP_GWAS.min_controls": 10,
  "rMAP_GWAS.significance_alpha": 0.05,
  "rMAP_GWAS.max_cpus": 16,
  "rMAP_GWAS.max_memory_gb": 64
}
```

---

## Phenotype encoding

The workflow internally creates a phenotype table:

```text
sample      case_control
case_001    1
case_002    1
control_001 0
control_002 0
```

Interpretation:

```text
positive beta = enriched in cases
negative beta = enriched in controls
```

The final report should always display the phenotype coding used for the analysis.

---

## Main outputs

```text
rMAP_GWAS_report.html
rMAP_GWAS_top_priority_hits.tsv
rMAP_GWAS_all_significant_hits.tsv
rMAP_GWAS_pyseer_unitig_assoc.tsv.gz
rMAP_GWAS_pyseer_gene_assoc.tsv.gz
rMAP_GWAS_gene_presence_absence.Rtab
rMAP_GWAS_gene_presence_absence.csv
rMAP_GWAS_phenotypes.tsv
rMAP_GWAS_population_structure_distances.tsv
rMAP_GWAS_manhattan.svg
rMAP_GWAS_qqplot.svg
rMAP_GWAS_top_hits_barplot.svg
rMAP_GWAS_enrichment_summary.tsv
rMAP_GWAS_run_provenance.json
```

---

## Top-priority hits table

The key interpreted output is:

```text
rMAP_GWAS_top_priority_hits.tsv
```

Recommended columns:

```text
rank
feature_id
feature_type
gene_name
product
contig
position
nearest_gene
case_present
case_total
case_frequency
control_present
control_total
control_frequency
enriched_in
beta
odds_ratio
pyseer_pvalue
q_value
bonferroni_threshold
annotation_source
samples_with_feature
notes
priority_score
```

---

## Enrichment interpretation

A feature is reported as **case-enriched** when:

```text
case_frequency > control_frequency
beta > 0
```

A feature is reported as **control-enriched** when:

```text
control_frequency > case_frequency
beta < 0
```

If the direction from frequencies & model beta disagree, the feature is flagged as:

```text
Check manually
```

---

## Priority scoring

The final report should not rank hits by p-value alone. A suggested priority score is:

```text
priority_score =
  -log10(q_value)
  + abs(log2_odds_ratio)
  + annotation_weight
  + recurrence_weight
```

Suggested weights:

```text
annotation_weight = 2 if feature is inside a named CDS or known AMR/virulence/MGE gene
annotation_weight = 1 if feature is near a gene
recurrence_weight = 1 if feature is present in at least 5 cases or 5 controls
```

Features with very low frequency, poor annotation, or inconsistent enrichment direction should be deprioritized.

---

## HTML report sections

The final `rMAP_GWAS_report.html` should include:

1. Run overview
2. Input cohort summary
3. Case/control balance
4. Read QC & trimming summary
5. Assembly QC summary
6. Annotation & pangenome summary
7. Population structure summary
8. GWAS model summary
9. Manhattan plot
10. QQ plot
11. Top-priority hits
12. Case-enriched hits
13. Control-enriched hits
14. Annotated gene-level associations
15. Annotated unitig/SNP-level associations
16. Warnings & limitations
17. Downloads & provenance

---

## Minimum viable version

The first stable version should implement:

```text
FASTP
FASTQC
MultiQC
Shovill
QUAST
Bakta or Prokka
Panaroo
Mash
unitig-caller
pyseer unitig GWAS
pyseer gene GWAS
hit annotation
priority scoring
portable HTML reporting
```

The reference-based SNP GWAS branch can be added later as an optional module.

---

## Important limitations

Microbial GWAS results can be confounded by:

- clonal population structure
- recombination
- case/control imbalance
- phenotype misclassification
- outbreak overrepresentation
- low sample size
- low-frequency features
- poor genome assemblies
- incomplete gene annotation
- plasmid fragmentation in short-read assemblies

The workflow should warn users when:

```text
number of cases < 10
number of controls < 10
```

The workflow should also issue a caution when:

```text
number of cases < 50 or number of controls < 50
```

because underpowered microbial GWAS can produce unstable associations.

---

## Repository structure

```text
rMAP-GWAS/
├── README.md
├── LICENSE
├── rMAP_GWAS.wdl
├── inputs/
│   ├── example_inputs.json
│   └── example_covariates.tsv
├── docker/
│   ├── rmap-gwas-pyseer-annotate/
│   │   └── Dockerfile
│   └── rmap-gwas-report/
│       └── Dockerfile
├── scripts/
│   ├── validate_inputs.py
│   ├── prepare_phenotypes.py
│   ├── annotate_hits.py
│   ├── prioritize_hits.py
│   ├── make_gwas_plots.py
│   └── build_html_report.py
├── templates/
│   └── rmap_gwas_report_template.html
├── docs/
│   ├── workflow_overview.md
│   ├── input_specification.md
│   └── output_interpretation.md
├── test_data/
│   └── README.md
└── examples/
    └── submission_notes.md
```

---

## Installation & execution

This repository is designed for execution with WDL-compatible engines such as:

- Cromwell
- wdl
- Docker

Example Cromwell command:

```bash
java -jar cromwell.jar run rMAP_GWAS.wdl -i inputs/example_inputs.json
```
---
## Docker images used by rMAP-GWAS

rMAP-GWAS is fully containerized. Each major workflow step is controlled by a Docker image input, so users can run the same WDL across local Cromwell, cloud Cromwell, or other WDL-compatible execution environments.

### Core workflow Docker images

| WDL input variable      | Default Docker image                                | Main tools / role                     | Workflow stage                                                                    |
| ----------------------- | --------------------------------------------------- | ------------------------------------- | --------------------------------------------------------------------------------- |
| `fastp_docker`          | `quay.io/biocontainers/fastp:0.23.4--hadf994f_2`    | `fastp`                               | Read trimming & FASTQ quality control                                           |
| `shovill_docker`        | `quay.io/biocontainers/shovill:1.1.0--hdfd78af_1`   | `shovill`, SPAdes backend             | De novo genome assembly                                                           |
| `quast_docker`          | `staphb/quast:5.2.0`                                | `QUAST`                               | Assembly quality assessment                                                       |
| `prokka_docker`         | `staphb/prokka:1.14.6`                              | `Prokka`                              | Genome annotation & GFF generation                                              |
| `panaroo_docker`        | `quay.io/biocontainers/panaroo:1.5.2--pyhdfd78af_0` | `Panaroo`                             | Pangenome construction & gene presence/absence matrix generation                |
| `mash_docker`           | `gmboowa/rmap-gwas-pyseer-annotate:0.2`             | `Mash`, Python utilities              | Pairwise genome distance matrix generation                                        |
| `pyseer_docker`         | `gmboowa/rmap-gwas-pyseer-annotate:0.2`             | `pyseer`                              | Gene presence/absence GWAS                                                        |
| `python_docker`         | `gmboowa/rmap-gwas-pyseer-annotate:0.2`             | Python reporting and table utilities  | Validation, phenotype generation, prioritization, plots, and HTML report creation |
| `hit_annotation_docker` | `gmboowa/rmap-gwas-pyseer-annotate:0.2`             | Python post-GWAS annotation utilities | GenBank-based annotation rescue of top GWAS hits                                  |

### Species-specific reference Docker images

Species-specific reference images provide curated reference files used for provenance tracking & post-GWAS annotation rescue. Each reference image is expected to expose a reference GenBank file, preferably through:

```text
RMAP_GWAS_REFERENCE_GENBANK=/opt/rmap-gwas/refs/reference.genbank
```

The same rMAP-GWAS WDL can therefore be used across bacterial species by changing the `reference_docker`, `reference_name`, & `reference_species` inputs. The value in the `Species / complex` column can be used as `reference_species`.

| Pathogen group | Species / complex                    | Suggested `reference_docker`                  | Suggested `reference_name`                 |
| -------------- | ------------------------------------ | --------------------------------------------- | ------------------------------------------ |
| MTBC           | *Mycobacterium tuberculosis* complex | `gmboowa/rmap-gwas-mtbc-refs:2026.06`         | `MTBC_2026_06` (`GCF_000195955.2`)         |
| ESKAPEE        | *Enterococcus faecium*               | `gmboowa/rmap-gwas-efaecium-refs:2026.06`     | `EFAECIUM_2026_06` (`GCF_000174395.2`)     |
| ESKAPEE        | *Staphylococcus aureus*              | `gmboowa/rmap-gwas-saureus-refs:2026.06`      | `SAUREUS_2026_06` (`GCF_000013425.1`)      |
| ESKAPEE        | *Klebsiella pneumoniae*              | `gmboowa/rmap-gwas-kpneumo-refs:2026.06`      | `KPNEUMO_2026_06` (`GCF_000240185.1`)      |
| ESKAPEE        | *Acinetobacter baumannii*            | `gmboowa/rmap-gwas-abaumannii-refs:2026.06`   | `ABAUMANNII_2026_06` (`GCF_000015425.1`)   |
| ESKAPEE        | *Pseudomonas aeruginosa*             | `gmboowa/rmap-gwas-paeruginosa-refs:2026.06`  | `PAERUGINOSA_2026_06` (`GCF_000006765.1`)  |
| ESKAPEE        | *Enterobacter cloacae*               | `gmboowa/rmap-gwas-enterobacter-refs:2026.06` | `ENTEROBACTER_2026_06` (`GCF_000025565.1`) |
| ESKAPEE        | *Escherichia coli*                   | `gmboowa/rmap-gwas-ecoli-refs:2026.06`        | `ECOLI_2026_06` (`GCF_000005845.2`)        |

### Reference image contents

Each species-specific reference image should contain at minimum:

| File type                    | Recommended path                              | Purpose                                                                  |
| ---------------------------- | --------------------------------------------- | ------------------------------------------------------------------------ |
| Reference FASTA              | `/opt/rmap-gwas/refs/reference.fasta`         | Reference genome sequence                                                |
| Reference GFF                | `/opt/rmap-gwas/refs/reference.gff`           | Reference feature coordinates                                            |
| Reference GenBank            | `/opt/rmap-gwas/refs/reference.genbank`       | Gene, locus tag & product annotation for post-GWAS hit interpretation |
| Optional trusted annotations | `/opt/rmap-gwas/refs/trusted_annotations.tsv` | Curated gene/product names or priority AMR/virulence annotations         |

Recommended environment variables inside each reference image:

```text
RMAP_GWAS_REFERENCE_FASTA=/opt/rmap-gwas/refs/reference.fasta
RMAP_GWAS_REFERENCE_GFF=/opt/rmap-gwas/refs/reference.gff
RMAP_GWAS_REFERENCE_GENBANK=/opt/rmap-gwas/refs/reference.genbank
RMAP_GWAS_TRUSTED_ANNOTATIONS=/opt/rmap-gwas/refs/trusted_annotations.tsv
```

### Example reference configuration

For an MTBC run:

```text
reference_docker  = "gmboowa/rmap-gwas-mtbc-refs:2026.06"
reference_name    = "MTBC_2026_06"
reference_species = "Mycobacterium tuberculosis complex"
```

For a *Klebsiella pneumoniae* run:

```text
reference_docker  = "gmboowa/rmap-gwas-kpneumo-refs:2026.06"
reference_name    = "KPNEUMO_2026_06"
reference_species = "Klebsiella pneumoniae"
```

### Notes

* The workflow uses Prokka-generated GFF files for Panaroo pangenome construction.
* Panaroo gene clusters may appear as IDs such as `group_2270`.
* The post-GWAS annotation rescue step maps prioritized Panaroo gene clusters back to the supplied species-specific GenBank reference, where possible.
* The final HTML report includes GWAS hit tables, reference annotation rescue results, QQ plot, Manhattan-style feature plot & run provenance.
* Docker image inputs can be updated as and when needed.


---

## Citation

If you use `rMAP-GWAS`, please cite this repository & the core tools used in your analysis. A formal citation will be added once the workflow is released.



---

## License

This project is released under the MIT License.

---

## Disclaimer

`rMAP-GWAS` is intended for research & surveillance support. GWAS associations require careful interpretation & should be validated using independent datasets, biological evidence & where appropriate, experimental or epidemiological follow-up.

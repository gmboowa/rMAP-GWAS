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
Terra sample set / paired-end FASTQ inputs
        │
        ├── sample_names = this.gwasmtbs.gwasmtb_id
        ├── read1s       = this.gwasmtbs.read1
        ├── read2s       = this.gwasmtbs.read2
        └── groups       = this.gwasmtbs.group
                         │
                         ▼
Input validation and cohort checks
        │
        ├── Confirm equal array lengths:
        │       sample_names, read1s, read2s, groups
        │
        ├── Confirm unique sample IDs
        │
        ├── Confirm valid group labels:
        │       case / control
        │       1 / 0
        │       true / false
        │       yes / no
        │
        ├── Confirm at least one case and one control
        │
        └── Generate validation report
                         │
                         ▼
Phenotype table generation
        │
        ├── Case samples coded as:      1
        ├── Control samples coded as:   0
        │
        ├── Output:
        │       <output_prefix>_phenotypes.tsv
        │
        └── Output:
                <output_prefix>_sample_groups.tsv
                         │
                         ▼
All selected samples routed together
        │
        ├── No species exclusion branch at this stage
        ├── No case/control separation for assembly or annotation
        └── Case/control labels are retained only for GWAS phenotype coding
                         │
                         ▼
Per-sample processing scatter
        │
        ├── For each selected sample:
        │
        │       Paired-end FASTQ files
        │              │
        │              ▼
        │       Read trimming and QC with fastp
        │              │
        │              ├── Trimmed R1 FASTQ
        │              ├── Trimmed R2 FASTQ
        │              ├── fastp HTML report
        │              └── fastp JSON report
        │              │
        │              ▼
        │       De novo genome assembly with Shovill
        │              │
        │              ├── Contigs FASTA
        │              └── Shovill log
        │              │
        │              ▼
        │       Assembly quality control with QUAST
        │              │
        │              └── QUAST report TSV
        │              │
        │              ▼
        │       Genome annotation with Prokka
        │              │
        │              ├── GFF annotation
        │              ├── GenBank file
        │              ├── Protein FASTA
        │              └── Nucleotide feature FASTA
        │
        └── End per-sample scatter
                         │
                         ▼
Cohort-level pangenome construction
        │
        ├── Collect all Prokka GFF files
        │
        ├── Run Panaroo pangenome analysis
        │
        ├── Clean pangenome graph using strict mode
        │
        ├── Remove invalid genes
        │
        └── Generate gene presence/absence matrix
                         │
                         ├── gene_presence_absence.csv
                         ├── gene_presence_absence.Rtab
                         └── panaroo_summary.txt
                         │
                         ▼
Cohort-level genome distance estimation
        │
        ├── Collect all Shovill assemblies
        ├── Rename assemblies by sample ID
        ├── Create Mash sketches
        └── Compute pairwise Mash distance matrix
                         │
                         └── mash_distances.tsv
                         │
                         ▼
Gene-based microbial GWAS
        │
        ├── Inputs:
        │       phenotype table
        │       Panaroo gene presence/absence matrix
        │       Mash distance matrix
        │
        ├── Run pyseer gene presence/absence GWAS
        │
        ├── Apply allele-frequency filters:
        │       min_af
        │       max_af
        │
        └── Correct for population structure using Mash distances
                         │
                         └── pyseer_gene_assoc.tsv
                         │
                         ▼
GWAS hit prioritization
        │
        ├── Parse pyseer association results
        ├── Link association hits to Panaroo/Prokka annotations
        ├── Calculate case and control gene frequencies
        ├── Estimate enrichment direction:
        │       case-enriched
        │       control-enriched
        │       check manually
        │
        ├── Rank features using:
        │       p-value or q-value
        │       effect direction
        │       odds ratio
        │       annotation availability
        │       recurrence across samples
        │
        └── Generate prioritized GWAS tables
                         │
                         ├── <output_prefix>_all_ranked_hits.tsv
                         ├── <output_prefix>_top_priority_hits.tsv
                         ├── <output_prefix>_all_significant_hits.tsv
                         └── <output_prefix>_enrichment_summary.tsv
                         │
                         ▼
Integrated reporting and provenance
        │
        ├── Merge validation report
        ├── Summarize case/control composition
        ├── Summarize Panaroo outputs
        ├── Display prioritized association hits
        ├── Record reference/provenance settings
        │       reference_docker
        │       reference_species
        │       reference_name
        │
        └── Generate final report files
                         │
                         ├── <output_prefix>_report.html
                         └── <output_prefix>_run_provenance.json
```

### Sample routing logic

rMAP-GWAS uses a cohort-level case–control design. Unlike rMAP-TB, which routes samples into MTBC and non-MTBC branches, rMAP-GWAS routes all selected samples through the same processing path. The `group` column is not used to split samples into separate computational branches during read processing, assembly, or annotation. Instead, it is used to construct the phenotype table required for microbial GWAS.

Each selected sample must have:

```text
sample ID
read1 FASTQ
read2 FASTQ
group label
```

The group label is interpreted as:

```text
case      → phenotype value 1
control   → phenotype value 0
```

The workflow also accepts common equivalents such as `1/0`, `true/false`, and `yes/no`.

### Terra sample-set routing

When running in Terra using a `gwasmtb_set`, the workflow expects the following mappings:

```text
sample_names = this.gwasmtbs.gwasmtb_id
read1s       = this.gwasmtbs.read1
read2s       = this.gwasmtbs.read2
groups       = this.gwasmtbs.group
```

This allows a Terra sample table to be used directly without manually constructing separate case and control file arrays.

### Case–control interpretation

The phenotype table is generated internally from the selected Terra sample set. Samples labelled as cases are coded as `1`, while controls are coded as `0`.

```text
sample          case_control
SRRxxxxxx       1
SRRyyyyyy       0
```

This phenotype table is passed to pyseer together with the Panaroo gene presence/absence matrix and Mash distance matrix.

### Population structure handling

Microbial GWAS can be confounded by clonal population structure, relatedness, outbreak clusters, lineage effects, and uneven case/control sampling. rMAP-GWAS therefore computes pairwise Mash distances from assembled genomes and supplies this distance matrix to pyseer during gene-based association testing.

```text
Assemblies
    ↓
Mash sketch
    ↓
Mash distance matrix
    ↓
pyseer population-structure-aware association testing
```

### Reference and provenance tracking

For MTBC runs, the workflow records the intended reference package:

```text
reference_docker  = gmboowa/rmap-gwas-mtbc-refs:2026.06
reference_species = Mycobacterium tuberculosis complex
reference_name    = MTBC_2026_06
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

A Bakta-based version can be implemented separately by replacing the annotation task and ensuring downstream GFF compatibility with Panaroo.

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

rMAP-GWAS identifies statistical associations between gene presence/absence patterns and case–control phenotype labels. These associations should be interpreted carefully alongside sample size, population structure, lineage distribution, outbreak clustering, phenotype definition, and biological plausibility. Candidate hits should be validated in independent datasets or by targeted experimental/epidemiological follow-up where possible.

## Overview

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
4. Read QC and trimming summary
5. Assembly QC summary
6. Annotation and pangenome summary
7. Population structure summary
8. GWAS model summary
9. Manhattan plot
10. QQ plot
11. Top-priority hits
12. Case-enriched hits
13. Control-enriched hits
14. Annotated gene-level associations
15. Annotated unitig/SNP-level associations
16. Warnings and limitations
17. Downloads and provenance

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

## Development status

This repository is under active development.

Planned milestones:

- [ ] Create WDL skeleton
- [ ] Add case/control input validation
- [ ] Add phenotype table generation
- [ ] Add read QC & trimming
- [ ] Add assembly & assembly QC
- [ ] Add Bakta/Prokka annotation
- [ ] Add Panaroo gene matrix generation
- [ ] Add unitig-caller feature generation
- [ ] Add Mash distance matrix generation
- [ ] Add pyseer unitig GWAS
- [ ] Add pyseer gene GWAS
- [ ] Add hit annotation & prioritization
- [ ] Add portable HTML report
- [ ] Add test dataset
- [ ] Add Example JSON

---

## Citation

If you use `rMAP-GWAS`, please cite this repository & the core tools used in your analysis. A formal citation will be added once the workflow is released.

Suggested citation format for now:

```text
Mboowa, G. rMAP-GWAS: rapid Microbial Analysis Pipeline for Genome-Wide Association Studies. GitHub: https://github.com/gmboowa/rMAP-GWAS
```

---

## License

This project is released under the MIT License.

---

## Disclaimer

`rMAP-GWAS` is intended for research & surveillance support. GWAS associations require careful interpretation & should be validated using independent datasets, biological evidence & where appropriate, experimental or epidemiological follow-up.

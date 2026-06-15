version 1.0

workflow rMAP_GWAS {
  input {
    # Terra sample-set inputs
    # Suggested Terra mappings when running on a gwasmtb_set:
    #   sample_names = this.gwasmtbs.gwasmtb_id
    #   read1s       = this.gwasmtbs.read1
    #   read2s       = this.gwasmtbs.read2
    #   groups       = this.gwasmtbs.group
    Array[String]+ sample_names
    Array[File]+ read1s
    Array[File]+ read2s
    Array[String]+ groups

    # Phenotype coding
    String phenotype_name = "case_control"
    String case_label = "case"
    String control_label = "control"

    # Analysis controls
    Float min_af = 0.01
    Float max_af = 0.99
    # Pyseer population-structure controls.
    # For small smoke tests, using too many MDS dimensions can make the null model singular.
    Int pyseer_max_dimensions = 2
    Boolean pyseer_force_no_distances = false
    Boolean pyseer_no_distances_fallback = true
    Float significance_alpha = 0.05
    String output_prefix = "rMAP_GWAS"

    # Runtime controls
    # These defaults are Terra-smoke-test friendly. For full cohorts, increase pangenome/gwas resources as needed.
    Int fastp_threads = 4
    Int assembly_threads = 4
    Int annotation_threads = 4
    Int pangenome_threads = 4
    Int gwas_threads = 4

    Int fastp_memory_gb = 8
    Int assembly_memory_gb = 32
    Int annotation_memory_gb = 16
    Int pangenome_memory_gb = 32
    Int gwas_memory_gb = 32

    Int fastp_disk_gb = 50
    Int assembly_disk_gb = 200
    Int pangenome_disk_gb = 300
    Int gwas_disk_gb = 300

    # Docker images
    String fastp_docker = "quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    String shovill_docker = "quay.io/biocontainers/shovill:1.1.0--hdfd78af_1"
    String quast_docker = "staphb/quast:5.2.0"
    String prokka_docker = "staphb/prokka:1.14.6"
    String panaroo_docker = "quay.io/biocontainers/panaroo:1.5.2--pyhdfd78af_0"
    # Combined linux/amd64 image for local Colima testing and Terra/Cromwell execution.
    # Contains pyseer, mash, Python, pandas, numpy, scipy, statsmodels, scikit-learn and tqdm.
    String mash_docker = "gmboowa/rmap-gwas-pyseer-annotate:0.2"
    String pyseer_docker = "gmboowa/rmap-gwas-pyseer-annotate:0.2"
    String python_docker = "gmboowa/rmap-gwas-pyseer-annotate:0.2"

    # Reference/provenance settings. This records the species-specific reference package used or intended for the run.
    String reference_docker = "gmboowa/rmap-gwas-mtbc-refs:2026.06"
    String reference_species = "Mycobacterium tuberculosis complex"
    String reference_name = "MTBC_2026_06"
  }

  # Use Terra sample-set arrays directly.
  Array[String] all_sample_names = sample_names
  Array[File] all_read1s = read1s
  Array[File] all_read2s = read2s

  # Shovill checks usable RAM inside the VM and fails if --ram is >= available RAM.
  # Terra VMs expose slightly less usable RAM than the WDL runtime request, so keep
  # the Shovill command-level RAM below the task runtime memory.
  Int assembly_shovill_ram_gb = if assembly_memory_gb > 16 then assembly_memory_gb - 8 else if assembly_memory_gb > 8 then assembly_memory_gb - 4 else assembly_memory_gb

  call VALIDATE_SAMPLE_SET_INPUTS {
    input:
      sample_names = sample_names,
      read1s = read1s,
      read2s = read2s,
      groups = groups,
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
      output_prefix = output_prefix,
      python_docker = python_docker
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
        docker = fastp_docker
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

    call PROKKA_ANNOTATE {
      input:
        sample_name = all_sample_names[i],
        assembly = SHOVILL_ASSEMBLE.contigs_fasta,
        threads = annotation_threads,
        memory_gb = annotation_memory_gb,
        docker = prokka_docker
    }
  }

  call PANAROO_PANGENOME {
    input:
      gffs = PROKKA_ANNOTATE.gff,
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
      python_docker = python_docker
  }

  call MERGE_RMAP_GWAS_REPORT {
    input:
      output_prefix = output_prefix,
      validation_report = VALIDATE_SAMPLE_SET_INPUTS.validation_report,
      phenotype_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_tsv,
      top_priority_hits = PRIORITIZE_GWAS_HITS.top_priority_hits,
      all_significant_hits = PRIORITIZE_GWAS_HITS.all_significant_hits,
      pyseer_gene_assoc = PYSEER_GENE_GWAS.pyseer_gene_assoc,
      panaroo_summary = PANAROO_PANGENOME.panaroo_summary,
      mash_distances = MASH_DISTANCE_MATRIX.mash_distances,
      reference_docker = reference_docker,
      reference_species = reference_species,
      reference_name = reference_name,
      python_docker = python_docker
  }

  output {
    File phenotype_tsv = PREPARE_PHENOTYPE_TABLE.phenotype_tsv
    Array[File] trimmed_read1s = FASTP_TRIM.trimmed_read1
    Array[File] trimmed_read2s = FASTP_TRIM.trimmed_read2
    Array[File] assemblies = SHOVILL_ASSEMBLE.contigs_fasta
    Array[File] quast_reports = QUAST_ASSEMBLY_QC.quast_report_tsv
    Array[File] gffs = PROKKA_ANNOTATE.gff
    File gene_presence_absence_csv = PANAROO_PANGENOME.gene_presence_absence_csv
    File gene_presence_absence_rtab = PANAROO_PANGENOME.gene_presence_absence_rtab
    File mash_distances = MASH_DISTANCE_MATRIX.mash_distances
    File pyseer_gene_assoc = PYSEER_GENE_GWAS.pyseer_gene_assoc
    File top_priority_hits = PRIORITIZE_GWAS_HITS.top_priority_hits
    File all_significant_hits = PRIORITIZE_GWAS_HITS.all_significant_hits
    File enrichment_summary = PRIORITIZE_GWAS_HITS.enrichment_summary
    File html_report = MERGE_RMAP_GWAS_REPORT.html_report
    File run_provenance = MERGE_RMAP_GWAS_REPORT.run_provenance_json
  }
}

task VALIDATE_SAMPLE_SET_INPUTS {
  input {
    Array[String]+ sample_names
    Array[File]+ read1s
    Array[File]+ read2s
    Array[String]+ groups
    String case_label
    String control_label
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
sample_names = """~{sep='\n' sample_names}""".strip().splitlines()
read1s = """~{sep='\n' read1s}""".strip().splitlines()
read2s = """~{sep='\n' read2s}""".strip().splitlines()
groups = """~{sep='\n' groups}""".strip().splitlines()
case_label = "~{case_label}"
control_label = "~{control_label}"

def norm(x):
    return str(x).strip().lower()

errors = []
warnings = []

if not (len(sample_names) == len(read1s) == len(read2s) == len(groups)):
    errors.append(
        f"sample_names/read1s/read2s/groups have unequal lengths: "
        f"sample_names={len(sample_names)}, read1s={len(read1s)}, read2s={len(read2s)}, groups={len(groups)}"
    )

if len(sample_names) != len(set(sample_names)):
    errors.append("Sample names must be unique within the selected Terra sample set.")

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
    out.write("rMAP-GWAS Terra sample-set input validation report\n")
    out.write("================================================\n")
    out.write(f"Cases: {len(case_names)}\n")
    out.write(f"Controls: {len(control_names)}\n")
    out.write(f"Total samples: {len(sample_names)}\n")
    out.write(f"Case label: {case_label}\n")
    out.write(f"Control label: {control_label}\n\n")
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
    String output_prefix
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
sample_names = """~{sep='\n' sample_names}""".strip().splitlines()
groups = """~{sep='\n' groups}""".strip().splitlines()
phenotype = "~{phenotype_name}"
case_label = "~{case_label}"
control_label = "~{control_label}"

def norm(x):
    return str(x).strip().lower()

valid_case = {norm(case_label), "case", "cases", "1", "true", "yes"}
valid_control = {norm(control_label), "control", "controls", "0", "false", "no"}

rows = []
errors = []
for sample, group in zip(sample_names, groups):
    g = norm(group)
    if g in valid_case:
        rows.append((sample, case_label, 1))
    elif g in valid_control:
        rows.append((sample, control_label, 0))
    else:
        errors.append(f"{sample}:{group}")

if len(sample_names) != len(groups):
    raise SystemExit(f"sample_names and groups have unequal lengths: {len(sample_names)} vs {len(groups)}")
if errors:
    raise SystemExit("Unrecognized group labels: " + ", ".join(errors[:10]))

with open("~{output_prefix}_phenotypes.tsv", "w") as out:
    out.write("sample\t" + phenotype + "\n")
    for sample, group_label, value in rows:
        out.write(f"{sample}\t{value}\n")

with open("~{output_prefix}_sample_groups.tsv", "w") as out:
    out.write("sample\tgroup\tphenotype_value\n")
    for sample, group_label, value in rows:
        out.write(f"{sample}\t{group_label}\t{value}\n")
PY
  >>>

  output {
    File phenotype_tsv = "~{output_prefix}_phenotypes.tsv"
    File sample_groups_tsv = "~{output_prefix}_sample_groups.tsv"
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

task PANAROO_PANGENOME {
  input {
    Array[File]+ gffs
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
  >>>

  output {
    File gene_presence_absence_csv = "panaroo_out/gene_presence_absence.csv"
    File gene_presence_absence_rtab = "panaroo_out/gene_presence_absence.Rtab"
    File panaroo_summary = "panaroo_summary.txt"
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
    Array[String]+ sample_names
    Array[File]+ assemblies
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

    odds_num = (case_present + 0.5) / (case_total - case_present + 0.5)
    odds_den = (control_present + 0.5) / (control_total - control_present + 0.5)
    odds_ratio = odds_num / odds_den if odds_den else float("inf")
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
        "pyseer_pvalue": "" if pval is None else f"{pval:.6g}",
        "q_value": "" if qval is None else f"{qval:.6g}",
        "priority_score": f"{priority_score:.4f}",
        "annotation_source": "Panaroo/Prokka",
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
    "enriched_in", "beta", "odds_ratio", "pyseer_pvalue", "q_value",
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

task MERGE_RMAP_GWAS_REPORT {
  input {
    String output_prefix
    File validation_report
    File phenotype_tsv
    File top_priority_hits
    File all_significant_hits
    File pyseer_gene_assoc
    File panaroo_summary
    File mash_distances
    String reference_docker
    String reference_species
    String reference_name
    String python_docker
  }

  command <<<
set -euo pipefail
export PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${PATH:-}
python <<'PY'
from pathlib import Path
import html, json, csv, datetime

prefix = "~{output_prefix}"
reference_docker = "~{reference_docker}"
reference_species = "~{reference_species}"
reference_name = "~{reference_name}"

def read_text(path, limit=20000):
    p = Path(path)
    if not p.exists():
        return ""
    txt = p.read_text(errors="replace")
    return txt[:limit]

def table_from_tsv(path, max_rows=50):
    p = Path(path)
    if not p.exists():
        return "<p>File not found.</p>"
    with p.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        return "<p>No rows.</p>"
    header = rows[0]
    body = rows[1:max_rows+1]
    out = ["<table>", "<thead><tr>"]
    out.extend(f"<th>{html.escape(x)}</th>" for x in header)
    out.append("</tr></thead><tbody>")
    for row in body:
        out.append("<tr>")
        out.extend(f"<td>{html.escape(x)}</td>" for x in row)
        out.append("</tr>")
    out.append("</tbody></table>")
    if len(rows) - 1 > max_rows:
        out.append(f"<p><em>Showing first {max_rows} of {len(rows)-1} rows.</em></p>")
    return "\n".join(out)

phenotype_rows = []
with open("~{phenotype_tsv}") as fh:
    header = fh.readline().strip().split("\t")
    for line in fh:
        if line.strip():
            phenotype_rows.append(line.strip().split("\t"))

cases = sum(1 for r in phenotype_rows if len(r) > 1 and r[1] == "1")
controls = sum(1 for r in phenotype_rows if len(r) > 1 and r[1] == "0")

html_doc = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{html.escape(prefix)} rMAP-GWAS report</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.45; color: #1f2933; }}
h1, h2, h3 {{ color: #102a43; }}
.section {{ margin-top: 32px; padding-top: 12px; border-top: 1px solid #d9e2ec; }}
.badge {{ display:inline-block; padding: 4px 8px; border-radius: 6px; background:#e0f2fe; margin-right:8px; }}
.warning {{ background:#fff7ed; border-left: 5px solid #fb923c; padding: 12px; }}
table {{ border-collapse: collapse; width: 100%; font-size: 12px; }}
th, td {{ border: 1px solid #bcccdc; padding: 6px; text-align: left; vertical-align: top; }}
th {{ background: #dbeafe; }}
pre {{ background:#f8fafc; border:1px solid #e2e8f0; padding:12px; overflow-x:auto; }}
</style>
</head>
<body>
<h1>rMAP-GWAS report</h1>
<p><strong>Run prefix:</strong> {html.escape(prefix)}</p>
<p><strong>Generated:</strong> {datetime.datetime.utcnow().isoformat()} UTC</p>
<p><strong>Reference package:</strong> {html.escape(reference_name)} | {html.escape(reference_species)} | {html.escape(reference_docker)}</p>

<div class="section">
<h2>1. Cohort summary</h2>
<span class="badge">Cases: {cases}</span>
<span class="badge">Controls: {controls}</span>
<span class="badge">Total: {cases + controls}</span>
<p>Phenotype coding: <strong>cases = 1</strong>, <strong>controls = 0</strong>. Positive beta values are interpreted as case-enriched when supported by case/control frequencies.</p>
</div>

<div class="section">
<h2>2. Input validation</h2>
<pre>{html.escape(read_text("~{validation_report}"))}</pre>
</div>

<div class="section">
<h2>3. Top priority GWAS hits</h2>
{table_from_tsv("~{top_priority_hits}", max_rows=100)}
</div>

<div class="section">
<h2>4. All significant hits</h2>
{table_from_tsv("~{all_significant_hits}", max_rows=100)}
</div>

<div class="section">
<h2>5. Panaroo summary</h2>
<pre>{html.escape(read_text("~{panaroo_summary}"))}</pre>
</div>

<div class="section">
<h2>6. Interpretation notes</h2>
<div class="warning">
<p>Microbial GWAS can be confounded by clonal population structure, recombination, outbreak overrepresentation, phenotype misclassification, and small sample size. Significant associations should be interpreted alongside lineage structure and, where possible, validated in independent datasets.</p>
</div>
</div>

<div class="section">
<h2>7. Output files</h2>
<ul>
<li>{html.escape(Path("~{top_priority_hits}").name)}</li>
<li>{html.escape(Path("~{all_significant_hits}").name)}</li>
<li>{html.escape(Path("~{pyseer_gene_assoc}").name)}</li>
<li>{html.escape(Path("~{phenotype_tsv}").name)}</li>
<li>{html.escape(Path("~{mash_distances}").name)}</li>
</ul>
</div>

</body>
</html>
"""

report_path = Path(prefix + "_report.html")
report_path.write_text(html_doc)

provenance = {
    "workflow": "rMAP-GWAS",
    "workflow_version": "0.1.0",
    "description": "Gene presence/absence microbial GWAS from paired-end reads using fastp, shovill, prokka, panaroo, mash and pyseer.",
    "reference": {
        "reference_docker": reference_docker,
        "reference_species": reference_species,
        "reference_name": reference_name
    },
    "phenotype_coding": {"case": 1, "control": 0},
    "cases": cases,
    "controls": controls,
    "generated_utc": datetime.datetime.utcnow().isoformat()
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

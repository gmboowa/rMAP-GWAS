# rMAP-GWAS test dataset

This folder contains a small four-sample MRSA/MSSA test dataset for local Cromwell/Docker execution of rMAP-GWAS.

The test is intended as a workflow smoke test. It confirms that the workflow can run locally, that Docker images are accessible, that gene and SNP-GWAS branches can execute, that outputs are generated, and that the HTML report renders correctly.

This test is not powered for biological interpretation.

---

## Test design

| Sample | Group | Display label | FASTQ files |
|---|---|---|---|
| ERR11728856 | control | MSSA | `ERR11728856_1.fastq.gz`, `ERR11728856_2.fastq.gz` |
| ERR11728866 | control | MSSA | `ERR11728866_1.fastq.gz`, `ERR11728866_2.fastq.gz` |
| ERR11728974 | case | MRSA | `ERR11728974_1.fastq.gz`, `ERR11728974_2.fastq.gz` |
| ERR11728986 | case | MRSA | `ERR11728986_1.fastq.gz`, `ERR11728986_2.fastq.gz` |

The FASTQ files are downsampled paired-end reads and are provided only for testing workflow execution.

---

## Folder contents

```text
test/
├── README.md
├── cromwell_local_1job.conf
├── test_4_tiny_mrsa_mssa_local_lowram.json
└── fastqs/
    ├── ERR11728856_1.fastq.gz
    ├── ERR11728856_2.fastq.gz
    ├── ERR11728866_1.fastq.gz
    ├── ERR11728866_2.fastq.gz
    ├── ERR11728974_1.fastq.gz
    ├── ERR11728974_2.fastq.gz
    ├── ERR11728986_1.fastq.gz
    └── ERR11728986_2.fastq.gz
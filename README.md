# Classifiers Assemble: metagenomics Nextflow pipeline (SLURM)

This combines BLAST (outformat 7), KMA, Centrifuge, CLARK, MetaPhlAn4, and Kraken2 and Bracken into a single Nextflow DSL2 workflow. The different analyses are done in separate processes in parallel with SLURM.

The paired-end reads are first checked with FastQC, then trimmed with fastp and a summary report is generated with MultiQC. The trimmed reads are used for the classification analyses, BLAST takes a separate FASTA input.

FastQC, fastp, MultiQC, BLASTn, Centrifuge, KMA, Kraken2, Bracken, and MetaPhlAn will be installed into separate conda environments when the script is run

## Dependencies

* Python >=3.5
* Nextflow
* Conda or mamba with access to bioconda and conda-forge channels

## Packages (installed from conda or mamba)
* FastQC
* fastp
* MuliQC
* Kraken2
* Bracken
* Centrifuge
* CLARK
* MetaPhlAn4
* KMA
* BLAST


Before starting, make sure the cluster the samples being run on has Nextflow installed and the Nextflow and nodes can find the tools being used.

## Files

- `main.nf`: pipeline and process definitions
- `nextflow.config`: SLURM profile, resource labels, and default parameters
- `params.run.yml`: parameter overrides

## Requirements

- Nextflow (DSL2; Nextflow 23.10+)
- A working SLURM cluster
- Database files installed on storage accessible from compute nodes
- Paths in params.run.yml file must resolve from the directory where Nextflow is launched, or be absolute paths

## Run

Run from this directory:

```bash
nextflow run main.nf -profile slurm -params-file params.run.yml
```

Can change the output directory with `--outdir /path/to/results`.

## Enable analyses

Each analysis is disabled by default. Set the `run_*` parameter to `true` in the params.run.yml file.

### Quality Cntrol (FastQC, fastp, MultiQC)

- Enabled when `run_qc: true`, set to true by default. When `run_qc: false` the classification analysis is run on the raw reads



### BLAST

- Input: FASTA files matching `blast_input`
- Database prefix: `blast_db`
- Output: `results/blast/sample.tab`
- Preserves `-perc_identity 98`, `-qcov_hsp_perc 100`, and `-outfmt 7` this can be changed on the `params.run.yml` file.

### Centrifuge

Currently at the species level
Centrifuge requires python >=3.5 in order to be executable.
- Input: one centrifuge job per PE sample

### CLARK

- The CLARK helper script is copied from the conda environment to the working directory
- Dwnload the databases prior to running the script.

### KMA

- Input: paired trimmed reads, changeable in the `params.run.yml` file by setting `run_qc: false`.
- The output files are published in `results/kma/`.

### Kraken2 and Bracken

- Input: paired trimmed reads, changeable in the `params.run.yml` file by setting `run_qc: false`.
- Kraken2 output, report and unclassified folders and files are in `results/kraken2/`.
- Switch to `run_bracken: true`, so Bracken runs after Kraken2 it uses the Kraken2 report as input and emits results to `results/bracken/`.
- The Bracken taxonomy level defaults to `S` (species), this can be changed in the `params.run.yml` file.

### MetaPhlAn

- Input: paired trimmed reads, changeable in the `params.run.yml` file by setting `run_qc: false` files matching `metaphlan_input_dir/metaphlan_glob`, in the `params.run.yml` file willbe used instead.
- Output: profile and intermediate Bowtie2 are in `results/metaphlan/`.
- Confirm that the index name matches the database on the system.

## Downloading databases before running the pipeline

The pipeline copies the CLARK helper script automatically to the working directory, then runs `set_target.sh`, `updateTaxonomy.sh`, and optionally `classify_metagenome.sh`.

Download databases needed prior to running script

### BLAST Database
Use [Error_Model_Build](https://github.com/Tj-Idowu/Error_Model_Build) "build_reference_library.py" and the "BLAST Database building" in the README to build the database.

### Centrifuge Database

```bash
git clone https://github.com/DaehwanKimLab/centrifuge
cd centrifuge/indices
make p_compressed
```
### KMA Database
Select and download database file from [KMA and CCMetagen database] (https://figshare.unimelb.edu.au/projects/CCMetagen_databases/195755) choice options.

### Kraken2 and Bracken Database
Download database from [https://benlangmead.github.io/aws-indexes/k2] (https://benlangmead.github.io/aws-indexes/k2)

### MetaPhlAn Database
Download the database in a folder outside the conda environment.

```bash
metaphlan --install --bowtie2db database folder
```

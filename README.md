# Classifiers Assemble: metagenomics Nextflow pipeline (SLURM)

This combines BLAST (outformat 7), KMA, Centrifuge, CLARK, MetaPhlAn4, and Kraken2 and Bracken into a single Nextflow DSL2 workflow. The different analyses are done in separate processes in parallel with SLURM.

BLASTn, Centrifuge, KMA, Kraken2, Bracken, and MetaPhlAn will be installed into separate conda environments when the script is run

## Dependencies

* Python >=3.5
* Nextflow
* Bioconda
* Kraken2
* Bracken
* Centrifuge
* CLARK
* MetaPhlAn4
* KMA
* BLAST

Before starting, make sure the cluster the samples are being run on has Nextflow installed and that Nextflow and computer nodes can find the tools being used.

## Files

- `main.nf`: pipeline and process definitions
- `nextflow.config`: SLURM profile, resource labels, and default parameters
- `params.run.yml`: parameter overrides

## Requirements

- Nextflow (DSL2; Nextflow 23.10+)
- A working SLURM cluster
- Database files installed on storage accessible from compute nodes
- Paths in the params.run.yml file must resolve from the directory where Nextflow is launched, or be absolute paths

## Run

Run from this directory:

```bash
nextflow run main.nf -profile slurm -params-file params.run.yml
```

Can change the output directory with `--outdir /path/to/results`.

## Enable analyses

Each analysis is disabled by default. Set the `run_*` parameter to `true` in the params.run.yml file.

### BLAST

- Input: FASTA files matching `blast_input`
- Database prefix: `blast_db`
- Output: `results/blast/sample.tab`
- Preserves `-perc_identity 98`, `-qcov_hsp_perc 100`, and `-outfmt 7`.

### Centrifuge

Currently at the species level
Centrifuge requires Python >= 3.5 to be executable.
- Input: one existing centrifuge sample sheet
- Index prefix: `centrifuge_index`
- Check the centrifuge sample-sheet format and output naming. Check the output declarations in the .tsv file.

### CLARK

- `clark_home` points to the directory where CLARK is installed. It contains `clean.sh`, `updateTaxonomy.sh`, `set_targets.sh`, and `classify_metagenome.sh`. Set `clark_home` in the `params.run.yml` file.
- `clark_read1` and `clark_read2` inputs the paired reads. The default `clark_sample` is `zero`; you can change it in the `params.run.yml` file.
- Output located in `results/clark/sample_results.txt`.
- The database building has been removed, so do that before running the script.

### KMA

- Input: paired files `sample_R1.fastq` and `sample_R2.fastq`, set in the `params.run.yml` file.
- The output files are published in `results/kma/`.

### Kraken2 and Bracken

- Input: paired files named `sample_R1.fastq` and `sample_R2.fastq`, set in the `params.run.yml` file.
- Kraken2 output, report and unclassified folders and files are in `results/kraken2/`.
- Switch to `run_bracken: true`, so Bracken runs after Kraken2; it uses the Kraken2 report as input and emits results to `results/bracken/`.
- The Bracken taxonomy level defaults to `S` (species); this can be changed in the `params.run.yml` file.

### MetaPhlAn

- Input: files matching `metaphlan_input_dir/metaphlan_glob`, in the `params.run.yml` file.
- Output: profile and intermediate Bowtie2 are in `results/metaphlan/`.
- Confirm that the index name matches the database on the system.

## Downloading databases before running the pipeline

The pipeline builds the CLARK process automatically, running the three CLARK setup scripts. It copies CLARK into the task work directory before modifying/building the database.

Download databases needed before running the script

### BLAST Database
Use [Error_Model_Build](https://github.com/Tj-Idowu/Error_Model_Build) "build_reference_library.py" and the "BLAST Database building" in the README to build the database.

### Centrifuge Database

```bash
git clone https://github.com/DaehwanKimLab/centrifuge
cd centrifuge/indices
make p_compressed
```
### KMA Database
Select and download the database file from the [KMA and CCMetagen database] (https://figshare.unimelb.edu.au/projects/CCMetagen_databases/195755) choice options.

### Kraken2 and Bracken Database
Download database from [https://benlangmead.github.io/aws-indexes/k2] (https://benlangmead.github.io/aws-indexes/k2)

### MetaPhlAn Database
Download the database in a folder outside the conda environment.

```bash
metaphlan --install --bowtie2db database folder
```

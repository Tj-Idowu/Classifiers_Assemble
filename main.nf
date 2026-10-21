nextflow.enable.dsl=2

/*
 * Metagenomics workflow: do QC analysis (FatQC, fastp and MultiQC), then
 * metagenomic analysis with: Kraken2, Bracken, KMA, Centrifuge and BLAST (outformat 7) and CLARK.
 * Each process is run in their own Conda environment (env/*.yml)
 * Configure paths and enable analyses in nextflow.config or with -params-file (params.run.yml).
 * All input/database paths should be absolute paths visible from compute nodes.
 */

params.run_qc         = params.run_qc         ?: false
params.run_blast      = params.run_blast      ?: false
params.run_centrifuge = params.run_centrifuge ?: false
params.run_clark      = params.run_clark      ?: false
params.run_kma        = params.run_kma        ?: false
params.run_kraken2    = params.run_kraken2    ?: false
params.run_metaphlan  = params.run_metaphlan  ?: false

params.outdir = params.outdir ?: 'results'

/*
 * Channel for The Paired Sample Directory
 */

def readPairs(dir, r1_suffix, r2_suffix) {
    Channel.fromFilePairs("${dir}/*{${r1_suffix},${r2_suffix}}", flat: true, checkIfExists: true)
}

/*
 * Quality Control of Reads
 */

process FASTQC {
    conda "${projectDir}/envs/fastqc.yml"
    tag "${sample}"
    label 'fastqc'
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        "${read1}" "${read2}"
    """
}

process FASTP {
    conda "${projectDir}envs/fastp.yml"
    tag "${sample}"
    label 'fastp'
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}_trimmed_R1.fastq.gz"), path("{sample}_trimmed_R2.fastq.gz"), emit: reads
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.html", emit: html

    script:
    def other = params.fastp_other_args ?: ''
    """
    fastp --in1 "${reads}" --in2 "${red2}" --out1 "${sample}_trimmed_R1.fastq.gz" --out2 "${sample}_trimmed_R2.fastq.gz" \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred ${params.fastp_qualified_quality} \\
        --length_required ${params.fastp_length_required} \\
        --thread ${task.cpus} \\
        --json "${sample}.fastp.json" --html "${sample}.fastp.html"
        ${other}
    """
}

process MULTIQC{
    conda "${projectDir}/envs/multiqc.yml"
    tag 'multiqc'
    label 'multiqc'
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy', overwrite: true

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    """
    multiqc --force --filename multiqc_report.html .
    """
}

/*
 * Classifier Analysis of Reads
 */

process BLASTN {
    conda "${projectDir}/envs/blast.yml"
    tag "${sample}"
    label 'blast'
    publishDir "${params.outdir}/blast", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(fasta)
    val db_prefix

    output:
    tuple val(sample), path("${sample}.tab"), emit: results

    script:
    """
    blastn \\
      -query "${fasta}" \\
      -db "${db_prefix}" \\
      -num_threads ${task.cpus} \\
      -perc_identity ${params.blast_perc_identity} \\
      -qcov_hsp_perc ${params.blast_qcov_hsp_perc} \\
      -out "${sample}.tab" \\
      -outfmt 7
    """
}

process CENTRIFUGE {
    conda "${projectDir}/envs/centrifuge.yml"
    tag "${sample}"
    label 'centrifuge'
    publishDir "${params.outdir}/centrifuge", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)
    val index_prefix

    output:
    tuple val(sample), path("${sample}_centrifuge.tsv"), path("${sample}_centrifuge_report.tsv"), emit: results

    script:
    """
    centrifuge -x "${index_prefix}" -1 "${read1}" -2 "${read2}" -p ${task.cpus} -S "${sample}_centrifuge.tsv" \\
        --report-file "${sample}_centrifuge_report.tsv"
    """
}

process CLARK {
    conda "${projectDir}/envs/clark.yml"
    tag "${sample}"
    label 'clark'
    publishDir "${params.outdir}/clark", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path (read1), path (read2)
    val clark_db

    output:
    path "${sample}_clark.csv", emit: results
    path "${sample}_clark_run.log", emit: log

    script:
    def update_taxonomy = params.clark_update_taxonomy ? './updateTaxonomy.sh' : 'true'
    """
    set -euo pipefail

    # Decompress files so they can be read by clark
    r1="${read1}"; r2="${read2}"
    case "\$r1" in *.gz) zcat "\$r1" > "${sample}_R1.fastq"; r1="${sample}_R1.fastq";; esac
    case "\$r2" in *.gz) zcat "\$r2" > "${sample}_R2.fastq"; r2="${sample}_R2.fastq";; esac
    
    # Copy the CLARK helper script from CONDA_PREFIX/opt/clark to the working directory.
    cp -a "\$CONDA_PREFIX/opt/clark" clark_work
    ln -sfn "\$CONDA_PREFIX/bin" clark_work/exe
    cd clark_work

    {
    ./set_target.sh "${clark_db}" ${params.clark_targets} --${params.clark_rank}
    ${update_taxonomy}
    ./classify_metagenome.sh -P "../\$r1" "../\$r2" -R "../${sample}_clark" -n ${task.cpus}
    } 2>&1 | tee "../&{sample}_clark_run.log"
    """
}

process KMA {
    conda "${projectDir}/envs/kma.yml"
    tag "${sample}"
    label 'kma'
    publishDir "${params.outdir}/kma", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)
    val db_prefix

    output:
    tuple val(sample), path("${sample}.*"), emit: results

    script:
    """
    kma -ipe "${read1}" "${read2}" -o "${sample}" -t_db "${db_prefix}" \\
        -t ${task.cpus} -1t1 -mem_mode -and -apm f -ef -matrix
    """
}

process KRAKEN2 {
    conda "${projectDir}/envs/kraken2_bracken.yml"
    tag "${sample}"
    label 'kraken2'
    publishDir "${params.outdir}/kraken2", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)
    val db_dir

    output:
    tuple val(sample), path("${sample}_out.tsv"), path("${sample}_report.tsv"), path("${sample}.unclassified_#.fastq"), emit: results

    script:
    def compressed = read1.name.endsWith('.gz') ? '--gzip-compressed' : ''
    """
    kraken2 --db "${db_dir}" \\
      --memory-mapping \\
      ${compressed} \\
      --output "${sample}_out.tsv" \\
      --report "${sample}_report.tsv" \\
      --threads ${task.cpus} \\
      --paired "${read1}" "${read2}" \\
      --unclassified-out "${sample}.unclassified_#.fastq"
    """
}

process BRACKEN {
    conda "${projectDir}/envs/kraken2_bracken.yml"
    tag "${sample}"
    label 'bracken'
    publishDir "${params.outdir}/bracken", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(kraken_output), path(kraken_report), path(unclassified)
    val db_dir

    output:
    path "${sample}_bracken.tsv", emit: results

    script:
    """
    bracken -d "${db_dir}" -i "${kraken_report}" -o "${sample}_bracken.tsv" \\
      -r ${params.bracken_read_length} \\
      -l ${params.bracken_level} \\
      -t ${params.bracken_threshold}
    """
}

process METAPHLAN {
    conda "${projectDir}/envs/metaphlan.yml"
    tag "${sample}"
    label 'metaphlan'
    publishDir "${params.outdir}/metaphlan", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)
    val database_dir

    output:
    path "${sample}.txt", emit: profile
    path "${sample}.bowtie2.bz2", emit: bowtie2out

    script:
    def input_reads = [reads].flatten().collect { it.toString() }.join('.')
    """
    metaphlan "${iput_reads}" \\
      --bowtie2out "${sample}.bowtie2.bz2" --input_type fastq \\
      --nproc ${task.cpus} \\
      --unclassified_estimation \\
      --index "${params.metaphlan_index}" \\
      --bowtie2db "${database_dir}" \\
      -o "${sample}.txt"
    """
}

workflow {
    // Quality Control on Reads: trimmed reads go on to Classifier Analysis of Reads
    trimmed = null
    if (params.run_qc) {
        if (!params.qc_input_dir)
            error 'run_qc=true requires qc_input_dir'
        qc_r1_suffix = params.qc_r1_suffix ?: params.fastq_r1_suffix
        qc_r2_suffix = params.qc_r2_suffix ?: params.fastq_r2_suffix
        qc_pairs = readPairs(params.qc_input_dir, qc_r1_suffix, qc_r2_suffix)
        FASTQC(qc_pairs)
        FASTP(qc_pairs)
        MULTIQC(FASTQC.out.zip.mix(FASTP.out.json).collect())
        trimmed = FASTP.out.reads
    }

    // Run trimmed reads if QC ran and raw reads if QC was not run
    def pairedReads = { String tool, raw_dir ->
        if (trimmed != null) return trimmed
        if (!raw_dir) error "run_${tool}=true with run_qc=false requires ${tool}_input_dir"
        readPairs(raw_dir, params-fastq_r1_suffix, params.fastq_r2_suffix)
    }

    if (params.run_blast) {
        if (!params.blast_input || !params.blast_db) error 'run_blast=true requires blast_input and blast_db'
        blast_ch = Channel.fromPath(params.blast_input, checkIfExists: true)
            .map { f -> tuple(f.baseName, f) }
        BLASTN(blast_ch, file(params.blast_db).toString())
    }

    if (params.run_centrifuge) {
        if (!params.centrifuge_index) error 'run_centrifuge=true requires centrifuge_index'
        CENTRIFUGE(
            pairedReads('centrifuge', params.centrifuge_input_dir),
            file(params.centrifuge_index).toString()
        )
    }

    if (params.run_clark) {
        if (!params.clark_db) error 'run_clark=true requires clark_db'
        if (trimmed != null) {
            clark_reads = trimmed
        } else {
            if (!params.clark_read1 || !params.clark_read2)
                error 'run_clark=true with run_qc=false requires clark_read1 and clark_read2'
            clark_reads = Channel.of(tuple(
                params.clark_sample ?: 'clark_sample',
                file(params.clark_read1, checkIfExists: true),
                file(params.clark_read2, checkIfExists: true)
            ))
        }
        CLARK(clark_reads, file(params.clark_db).toString())
    }

    if (params.run_kma) {
        if (!params.kma_input_db) error 'run_kma=true requires kma_db'
        KMA(pairedReads('kma', params.kma_input_dir), file(params.kma_db).toString())
    }

    if (params.run_kraken2) {
        if (!params.kraken_db) error 'run_kraken2=true requires kraken_db'
        KRAKEN2(pairedReads('kraken', params.kraken_input_dir), file(params.kraken_db).toString())
        if (params.run_braken) {
            BRACKEN(KRAKEN2.out.results, file(params.kraken_db).toString())
        }
    }

    if (params.run_metaphlan) {
        if (!params.metaphlan_db) error 'run_metaphlan=true requires metaphlan_db'
        if (trimmed != null) {
            metaphlan_reads = trimmed.map {sample, r1, r2 -> tuple(sample, [r1, r2]) }
        } else {
            if (!params.metaphlan_input_dir)
                error 'run_metaphlan=true with run_qc=false requires metaphlan_input_dir'
        metaphlan_reads = Channel.fromPath("${params.metaphlan_input_dir}/${params.metaphlan_glob}", checkIfExists: true)
            .map { f ->
                def sample = f.name.replaceFirst(java.util.regex.Pattern.quote(params.metaphlan_suffix) + '$', '')
                tuple(sample, f)
            }
        }
        METAPHLAN(metaphlan_reads, file(params.metaphlan_db).toString())
    }
}

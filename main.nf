nextflow.enable.dsl=2

/*
 * Metagenomics workflow for running Kraken2, Bracken, KMA, Centrifuge, BLAST (outformat 7) and CLARK.
 * Configure paths and enable analyses in nextflow.config or with -params-file (params.run.yml).
 * All input/database paths should be absolute paths visible from compute nodes.
 */

params.run_blast      = params.run_blast      ?: false
params.run_centrifuge = params.run_centrifuge ?: false
params.run_clark      = params.run_clark      ?: false
params.run_kma        = params.run_kma        ?: false
params.run_kraken2    = params.run_kraken2    ?: false
params.run_metaphlan  = params.run_metaphlan  ?: false

params.outdir = params.outdir ?: 'results'

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
    conda "${projectDir}/envs/centrifuge.yml
    tag 'centrifuge'
    label 'centrifuge'
    publishDir "${params.outdir}/centrifuge", mode: 'copy', overwrite: true

    input:
    path sample_sheet
    val index_prefix

    output:
    path 'centrifuge*', optional: true, emit: results
    path '*.tsv', optional: true, emit: tables
    path '*.txt', optional: true, emit: text

    script:
    """
    centrifuge \
        -x "${index_prefix}" --sample-sheet "${sample_sheet}"
    """
}

process CLARK_ANALYSIS {
    label 'clark'
    publishDir "${params.outdir}/clark", mode: 'copy', overwrite: true

    input:
    path clark_home
    path read1
    path read2
    val sample

    output:
    path "${sample}_results.txt", emit: results
    path 'clark_run.log', emit: log

    script:
    """
    set -euo pipefail

    cp -a "${clark_home}" clark_work
    cd clark_work

    ./clean.sh
    ./updateTaxonomy.sh
    ./set_targets.sh CLARK_DB bacteria --species
    ./classify_metagenome.sh -P "../${read1}" "../${read2}" -R "../${sample}_results.txt"
    """
}

process KMA {
    conda "${projectDir}/envs/kma.yml
    tag "${sample}"
    label 'kma'
    publishDir "${params.outdir}/kma", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)
    val db_prefix
    path kma_exe

    output:
    tuple val(sample), path("${sample}.*"), emit: results

    script:
    """
    kma \\
        -ipe "${read1}" "${read2}" -o "${sample}" -t_db "${db_prefix}" \\
        -1t1 -mem_mode -and -apm f -ef -matrix
    """
}

process KRAKEN2 {
    conda "${projectDir}/envs/kraken2_bracken.yml
    tag "${sample}"
    label 'kraken2'
    publishDir "${params.outdir}/kraken2", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)
    val db_dir

    output:
    tuple val(sample), path("${sample}_out.tsv"), path("${sample}_report.tsv"), path("${sample}.unclassified_#.fastq"), emit: results

    script:
    """
    kraken2 --db "${db_dir}" \\
      --memory-mapping \\
      --output "${sample}_out.tsv" \\
      --report "${sample}_report.tsv" \\
      --threads ${task.cpus} \\
      --paired "${read1}" "${read2}" \\
      --unclassified-out "${sample}.unclassified_#.fastq"
    """
}

process BRACKEN {
    conda "${projectDir}/envs/kraken2_bracken.yml
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
    bracken -d "${db_dir}" \\
      -i "${kraken_report}" \\
      -o "${sample}_bracken.tsv" \\
      -r ${params.bracken_read_length} \\
      -l ${params.bracken_level} \\
      -t ${params.bracken_threshold}
    """
}

process METAPHLAN {
    conda "${projectDir}/envs/metaphlan.yml
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
    """
    metaphlan "${reads}" \\
      --bowtie2out "${sample}.bowtie2.bz2" \\
      --input_type fastq \\
      --unclassified_estimation \\
      --index "${params.metaphlan_index}" \\
      --bowtie2db "${database_dir}" \\
      -o "${sample}.txt"
    """
}

workflow {
    if (params.run_blast) {
        if (!params.blast_input || !params.blast_db) error 'run_blast=true requires blast_input and blast_db'
        blast_ch = Channel.fromPath(params.blast_input, checkIfExists: true)
            .map { f -> tuple(f.baseName, f) }
        BLASTN(blast_ch, Channel.value(new File(params.blast_db).absolutePath))
    }

    if (params.run_centrifuge) {
        if (!params.centrifuge_sample_sheet || !params.centrifuge_index || !params.centrifuge_exe)
            error 'run_centrifuge=true requires centrifuge_sample_sheet, centrifuge_index, and centrifuge_exe'
        CENTRIFUGE(
            Channel.value(file(params.centrifuge_sample_sheet)),
            Channel.value(params.centrifuge_index)
        )
    }

    if (params.run_clark) {
        if (!params.clark_home || !params.clark_read1 || !params.clark_read2)
            error 'run_clark=true requires clark_home, clark_read1, and clark_read2'
        clark_r1 = file(params.clark_read1)
        clark_r2 = file(params.clark_read2)
        clark_sample = params.clark_sample ?: 'clark_sample'
        CLARK_ANALYSIS(
            Channel.value(file(params.clark_home)),
            Channel.value(clark_r1),
            Channel.value(clark_r2),
            Channel.value(clark_sample)
        )
    }

    if (params.run_kma) {
        if (!params.kma_input_dir || !params.kma_db || !params.kma_exe)
            error 'run_kma=true requires kma_input_dir, kma_db, and kma_exe'
        kma_r1 = Channel.fromPath("${params.kma_input_dir}/*${params.fastq_r1_suffix}", checkIfExists: true)
            .map { r1 ->
                def sample = r1.name.replaceFirst(java.util.regex.Pattern.quote(params.fastq_r1_suffix) + '$', '')
                def r2path = r1.toString().replaceFirst(java.util.regex.Pattern.quote(params.fastq_r1_suffix) + '$', params.fastq_r2_suffix)
                def r2 = file(r2path)
                if (!r2.exists()) error "Missing paired read for ${r1}: expected ${r2path}"
                tuple(sample, r1, r2)
            }
        KMA(kma_r1, Channel.value(new File(params.kma_db).absolutePath), Channel.value(file(params.kma_exe)))
    }

    if (params.run_kraken2) {
        if (!params.kraken_input_dir || !params.kraken_db)
            error 'run_kraken2=true requires kraken_input_dir and kraken_db'
        kraken_pairs = Channel.fromPath("${params.kraken_input_dir}/*${params.fastq_r1_suffix}", checkIfExists: true)
            .map { r1 ->
                def sample = r1.name.replaceFirst(java.util.regex.Pattern.quote(params.fastq_r1_suffix) + '$', '')
                def r2path = r1.toString().replaceFirst(java.util.regex.Pattern.quote(params.fastq_r1_suffix) + '$', params.fastq_r2_suffix)
                def r2 = file(r2path)
                if (!r2.exists()) error "Missing paired read for ${r1}: expected ${r2path}"
                tuple(sample, r1, r2)
            }
        KRAKEN2(kraken_pairs, Channel.value(new File(params.kraken_db).absolutePath))
        if (params.run_bracken) {
            BRACKEN(KRAKEN2.out.results, Channel.value(new File(params.kraken_db).absolutePath))
        }
    }

    if (params.run_metaphlan) {
        if (!params.metaphlan_input_dir || !params.metaphlan_db)
            error 'run_metaphlan=true requires metaphlan_input_dir and metaphlan_db'
        metaphlan_reads = Channel.fromPath("${params.metaphlan_input_dir}/${params.metaphlan_glob}", checkIfExists: true)
            .map { f ->
                def sample = f.name.replaceFirst(java.util.regex.Pattern.quote(params.metaphlan_suffix) + '$', '')
                tuple(sample, f)
            }
        METAPHLAN(metaphlan_reads, Channel.value(new File(params.metaphlan_db).absolutePath))
    }
}

// adapted from nf-core/rnaseq
process SORTMERNA {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sortmerna:4.3.6--h9ee0642_0' :
        'biocontainers/sortmerna:4.3.6--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)
    path  fastas

    output:
    tuple val(meta), path("*non_rRNA.fastq.gz")    , emit: reads
    tuple val(meta), path("*aligned_rRNA.fastq.gz"), emit: aligned_reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        sortmerna \\
            ${'--ref '+fastas.join(' --ref ')} \\
            --reads $reads \\
            --threads $task.cpus \\
            --workdir . \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            $args

        mv non_rRNA_reads.f*q.gz ${prefix}_non_rRNA.fastq.gz
        mv rRNA_reads.f*q.gz ${prefix}_aligned_rRNA.fastq.gz
        mv rRNA_reads.log ${prefix}.sortmerna.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
        END_VERSIONS
        """
    } else {
        """
        sortmerna \\
            ${'--ref '+fastas.join(' --ref ')} \\
            --reads ${reads[0]} \\
            --reads ${reads[1]} \\
            --threads $task.cpus \\
            --workdir . \\
            --aligned rRNA_reads \\
            --fastx \\
            --other non_rRNA_reads \\
            --paired_in \\
            --out2 \\
            $args

        mv non_rRNA_reads_fwd.f*q.gz ${prefix}_1.non_rRNA.fastq.gz
        mv non_rRNA_reads_rev.f*q.gz ${prefix}_2.non_rRNA.fastq.gz
        mv rRNA_reads.log ${prefix}.sortmerna.log
        mv rRNA_reads_fwd.f*q.gz ${prefix}_1.aligned_rRNA.fastq.gz
        mv rRNA_reads_rev.f*q.gz ${prefix}_2.aligned_rRNA.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
        END_VERSIONS
        """
    }
}
process RIBODETECTOR {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribodetector:0.3.0--pyhdfd78af_0' :
        'biocontainers/ribodetector:0.3.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads)
    
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
    """
    ribodetector_cpu \\
        -i $reads \\
        -o ${prefix}_non_rRNA.fastq.gz \\
        -r ${prefix}_aligned_rRNA.fastq.gz \\
        --log ${prefix}.ribodetector.log \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribodetector_cpu: \$(ribodetector_cpu --version | sed -e "s/ribodetector_cpu v//g")
    END_VERSIONS
    """
}


process BEDTOOLS_GETFASTA {
    tag "$bed"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bed) // adjusted by me
    path fasta // genome.fasta

    output:
    tuple val(meta), path("*.fa"), emit: peak_fasta // adjusted by me
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${bed.baseName}"
    """
    bedtools \\
        getfasta \\
        $args \\
        -fi $fasta \\
        -bed $bed \\
        -fo ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
process PURECLIP {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pureclip:1.3.1--0' :
        'biocontainers/pureclip:1.3.1--0'}"

    
    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    path(fasta)
    
    output:
    tuple val(meta), path("*_sites.bed")    , emit: sites_bed
    tuple val(meta), path("*_regions.bed")  , emit: regions_bed
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pureclip \\
        -i $bam \\
        -bai $bai \\
        -g $fasta \\
        -ld \\
        -nt 8 \\
        -iv 'chr1;chr2;chr3;' \\
        $args \\
        -o PureCLIP.${prefix}.crosslink_sites.bed> \\
        -or PureCLIP.${prefix}.crosslink_regions.bed>

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pureclip: \$(pureclip --version | sed -e "s/pureclip v//g")
    END_VERSIONS
    """
}
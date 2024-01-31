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
    def control_bam = ''
    def control_bai = ''
    if (meta.control_bam != ' ' && meta.control_bai != ' '){
        control_bam = "-ibam ${meta.control_bam}"
        control_bai = "-ibai ${meta.control_bai}"
    } 

    """
    pureclip \\
        -i $bam \\
        -bai $bai \\
        -g $fasta \\
        -nt ${task.cpus} \\
        $control_bam \\
        $control_bai \\
        $args \\
        -o ${prefix}_pureclip_crosslink_sites_bed7.bed \\
        -or ${prefix}_pureclip_crosslink_regions.bed
    
    cat ${prefix}_pureclip_crosslink_sites_bed7.bed | cut -f 1,2,3,4,5,6 > ${prefix}_pureclip_crosslink_sites.bed
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pureclip: \$(pureclip --version | sed -e "s/pureclip v//g")
    END_VERSIONS
    """
}
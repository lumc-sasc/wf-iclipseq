process STREME {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meme:5.5.5--pl5321hda358d9_0' :
        'biocontainers/meme:5.5.5--0'}"

    input:
    tuple val(meta), path(fa)
    
    output:
    path "*.html"          , emit: html
    path "*.txt"           , emit: txt              , optional: true   
    path "*.xml"           , emit: xml              , optional: true        
    path "*_sequences.tsv" , emit: sequences_tsv    , optional: true   
    path "*_sites.tsv"     , emit: sites_tsv        , optional: true   
    path  "versions.yml"   , emit: versions         , optional: true   

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    streme \\
        --p $fa \\
        --oc . \\
        $args

    mv sequences.tsv ${prefix}_sequences.tsv
    mv sites.tsv ${prefix}_sites.tsv
    mv streme.html ${prefix}_streme.html
    mv streme.txt ${prefix}_streme.txt
    mv streme.xml ${prefix}_streme.xml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        streme: \$(streme --version | sed -e "s/streme v//g")
    END_VERSIONS
    """
}
process MULTIQC {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.17--pyhdfd78af_0' :
        'biocontainers/multiqc:1.17--pyhdfd78af_0' }"

    input:
    // path  multiqc_files, stageAs: "*/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)

    // paths to files

    // raw fastqc
    path ('fastqc/raw/*')

    // trimming
    path ('fastqc/trim/*')
    path ('trim_log/*')

    // filtering rRNA
    path ('sortmerna/*')
    path ('fastqc/sortmerna/*')

    // alignment and dedup
    path ('star/*')
    path ('star/genecounts/')
    path ('samtools/stats/*')
    path ('samtools/flagstat/*')
    path ('samtools/idxstats/*')
    path ('fastqc/dedup/*')

    // unmapped
    path ('fastqc/unmapped/*')

    // premapping
    path ('bowtie2/align/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    """
    multiqc \\
        --force \\
        $args \\
        $config \\
        $extra_config \\
        .    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_data
    touch multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
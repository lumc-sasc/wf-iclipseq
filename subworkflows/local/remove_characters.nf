// Decompress .gz, remove whitespaces and special characters and compress to .gz

include { GUNZIP as DECOMP                    } from '../../modules/nf-core/gunzip/main'
include { GZIP as COMP                        } from '../../modules/local/gzip'
include { LINUX_COMMAND as REMOVE_CHAR        } from '../../modules/goodwright/linux/command/main'

workflow REMOVE_CHARACTERS {
    take:
    reads                                   // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    // unzip to fastq format
    DECOMP ( reads )
    ch_versions = ch_versions.mix(DECOMP.out.versions)

    // remove special characters and whitespaces from read name according to Busch et al (2020)
    REMOVE_CHAR ( DECOMP.out.gunzip, [], false )
    ch_versions = ch_versions.mix(REMOVE_CHAR.out.versions)

    // convert back to fastq.gz format
    COMP ( REMOVE_CHAR.out.file )

    emit:
    reads          = COMP.out.gunzip       // channel: [ val(meta), [ reads ] ]
    versions       = ch_versions           // channel: [ versions.yml ]
}



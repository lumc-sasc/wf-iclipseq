// to increase the size of the binding sites and do motif analysis using STREME from the MEME suite

include { BEDTOOLS_SLOP               } from '../../modules/nf-core/bedtools/slop/main.nf'
include { BEDTOOLS_GETFASTA           } from '../../modules/nf-core/bedtools/getfasta/main.nf'
include { STREME                      } from '../../modules/local/streme.nf'

workflow MOTIF_ANALYSIS {
    take:
    peaksites_bed           // channel: [ val(meta), [ bed ] ]
    chrom_sizes             // channel: [ sizes ]
    fasta                   // channel: [ fasta ]

    main:
    ch_versions = Channel.empty()

    // resize binding site from 1 to chosen nt (check config file)
    BEDTOOLS_SLOP (
        peaksites_bed,
        chrom_sizes
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SLOP.out.versions)

    // convert bed to fasta format
    BEDTOOLS_GETFASTA (
        BEDTOOLS_SLOP.out.bed, 
        fasta
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)

   // perform motif analysis using STREME, generates a .html report (and optionally .tsv, .txt, .xml files)
    STREME (
        BEDTOOLS_GETFASTA.out.peak_fasta
    )

    emit:
    peak_fasta                 = BEDTOOLS_GETFASTA.out.peak_fasta     // channel: [ val(meta), [ fasta ] ]
    versions                   = ch_versions                          // channel: [ versions.yml ]
}
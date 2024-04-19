// based on rnaseq/subworkflows/local/prepare_genome.nf

// get modules
include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_STAR_INDEX         } from '../../modules/nf-core/untar/main'
include { CUSTOM_GETCHROMSIZES              } from '../../modules/nf-core/custom/getchromsizes/main'
include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/star/genomegenerate/main'

workflow PREPARE_GENOME {
    take:
    fasta           // file: path to genome.fasta
    gtf             // file: path to genome.gtf
    star_index      // directory: path to star index

    main:
    ch_versions = Channel.empty()

    // uncompress genome fasta file if required
    if (fasta) {
        if (fasta.endsWith('.gz')) {
            ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            ch_fasta = Channel.value(file(fasta))
        }
    }

    // uncompress GTF annotation file
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value(file(gtf))
        }
    }

    // create chromosome sizes file
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    
    // uncompress STAR index or generate from scratch if required
    ch_star_index = Channel.empty()
    if (star_index) { // if index already exists
        if (star_index.endsWith('.tar.gz')) {
            ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map { it[1] }
            ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
        } else {
            ch_star_index = Channel.value(file(star_index))
        }
    } else { // otherwise, make from scratch
            ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }



    emit:
    fasta               = ch_fasta          // channel: path(genome.fasta)
    gtf                 = ch_gtf            // channel: path(genome.gtf)
    chrom_sizes         = ch_chrom_sizes    // channel: path(genome.sizes)
    fai                 = ch_fai            // channel: path(genome.fai)
    star_index          = ch_star_index     // channel: path(star/index/)

    versions            = ch_versions.ifEmpty(null)

}
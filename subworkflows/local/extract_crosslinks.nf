// workflow proposed as in Busch et al 2020
// uses modules from nf-core and clipseq 
// partially adapted from clipseq

// nf-core modules
include { BEDTOOLS_BAMTOBED                            } from '../../modules/nf-core/bedtools/bamtobed/main.nf'
include { BEDTOOLS_SHIFT                               } from '../../modules/nf-core/bedtools/shift/main.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_POS } from '../../modules/nf-core/bedtools/genomecov/main.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_NEG } from '../../modules/nf-core/bedtools/genomecov/main.nf'

// from goodwright/clipseq
include { LINUX_COMMAND as SELECT_BED_POS              } from '../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as SELECT_BED_NEG              } from '../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as MERGE_AND_SORT              } from '../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as CROSSLINK_COVERAGE          } from '../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as CROSSLINK_NORMCOVERAGE      } from '../../modules/goodwright/linux/command/main.nf'

// doing it separately
include { LINUX_COMMAND as MERGE_AND_SORT_NEG          } from '../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as CROSSLINK_COVERAGE_NEG      } from '../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as CROSSLINK_NORMCOVERAGE_NEG  } from '../../modules/goodwright/linux/command/main.nf'

include { LINUX_COMMAND as MERGE_AND_SORT_POS          } from '../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as CROSSLINK_COVERAGE_POS      } from '../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as CROSSLINK_NORMCOVERAGE_POS  } from '../../modules/goodwright/linux/command/main.nf'



workflow EXTRACT_CROSSLINKS {
    take:
    bam         // channel: [ val(meta), [ bam ] ]
    chrom_sizes // channel: [ sizes ]

    main:
    ch_versions = Channel.empty()

    // MODULE from nf-core 
    // convert bam to bed
    BEDTOOLS_BAMTOBED ( bam )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)


    // MODULE from nf-core 
    // shift bed files by one base pair into 5' dir using bedtools shift
    BEDTOOLS_SHIFT (
        BEDTOOLS_BAMTOBED.out.bed,
        chrom_sizes.map { [ [:], it ] }
    )
    ch_versions = ch_versions.mix(BEDTOOLS_SHIFT.out.versions)


    // extract the 5' position of the new intervals and pile them up using bedtools genomecov
    // this needs to be run separately for each strand (positive, negative)

    // MODULE from nf-core 
    // depth at each position on the pos strand

    BEDTOOLS_GENOMECOV_POS (
        BEDTOOLS_SHIFT.out.bed.map{  [ it[0], it[1], 1 ] },
        chrom_sizes,
        'pos.bed'
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_POS.out.versions)

    // depth at each position on the neg strand
    BEDTOOLS_GENOMECOV_NEG (
        BEDTOOLS_SHIFT.out.bed.map{ [ it[0], it[1], 1 ] },
        chrom_sizes,
        'neg.bed'
    )

    // output consists of crosslink events as coverage tracks, # of crosslink events on each crosslinked nucleotide along the genome
    // bedgraph file format!
    // can normalize these by providing a scaling factor to bedtool genomecov

    // NOTE: from goodwright/clipseq

    /*
    * MODULE: Select columns in BED file using AWK
    */
    SELECT_BED_POS (
        BEDTOOLS_GENOMECOV_POS.out.genomecov,
        [],
        false
    )
    
    SELECT_BED_NEG (
        BEDTOOLS_GENOMECOV_NEG.out.genomecov,
        [],
        false
    )

    /*
    * CHANNEL: Join POS/NEG files into one channel so they can be merged in the next module
    */
    ch_merge_and_sort_input = SELECT_BED_POS.out.file
        .map{ [ it[0].id, it[0], it[1] ] }
        .join( SELECT_BED_NEG.out.file.map{ [ it[0].id, it[0], it[1] ] } )
        .map { [ it[1], [ it[2], it[4] ] ] }
    //EXAMPLE CHANNEL STRUCT: [ [id:test], [ BED(pos), BED(neg) ] ]
    //ch_merge_and_sort_input | view 

    /*
    * MODULE: Select columns in BED file using AWK
    */
    MERGE_AND_SORT (
        ch_merge_and_sort_input,
        [],
        false
    )

    /*
    * MODULE: Create coverage track using AWK
    */
    CROSSLINK_COVERAGE (
        MERGE_AND_SORT.out.file,
        [],
        false
    )

    /*
    * MODULE: Create normalised coverage track using AWK
    */
    CROSSLINK_NORMCOVERAGE (
        MERGE_AND_SORT.out.file,
        [],
        true
    )


    // doing it separately
    // SENSE / CODING STRAND
    MERGE_AND_SORT_POS ( SELECT_BED_POS.out.file, [], false )
    CROSSLINK_COVERAGE_POS ( MERGE_AND_SORT_POS.out.file, [], false)
    CROSSLINK_NORMCOVERAGE_POS ( MERGE_AND_SORT_POS.out.file, [], true )


    // ANTISENSE / TEMPLATE STRAND
    MERGE_AND_SORT_NEG ( SELECT_BED_NEG.out.file, [], false )
    CROSSLINK_COVERAGE_NEG ( MERGE_AND_SORT_NEG.out.file, [], false )
    CROSSLINK_NORMCOVERAGE_NEG ( MERGE_AND_SORT_NEG.out.file, [], true )


    emit:
        merged_bed                      = MERGE_AND_SORT.out.file           // channel: [ val(meta), [ bed ] ]
        merged_coverage                 = CROSSLINK_COVERAGE.out.file       // channel: [ val(meta), [ bedgraph ] ]
        merged_coverage_normalized      = CROSSLINK_NORMCOVERAGE.out.file   // channel: [ val(meta), [ bedgraph ] ]
        versions                        = ch_versions                       // channel: [ versions.yml ]

        pos_bed                  = MERGE_AND_SORT_POS.out.file           // channel: [ val(meta), [ bed ] ]
        pos_coverage             = CROSSLINK_COVERAGE_POS.out.file       // channel: [ val(meta), [ bedgraph ] ]
        pos_coverage_normalized  = CROSSLINK_NORMCOVERAGE_POS.out.file   // channel: [ val(meta), [ bedgraph ] ]

        neg_bed                  = MERGE_AND_SORT_NEG.out.file           // channel: [ val(meta), [ bed ] ]
        neg_coverage             = CROSSLINK_COVERAGE_NEG.out.file       // channel: [ val(meta), [ bedgraph ] ]
        neg_coverage_normalized  = CROSSLINK_NORMCOVERAGE_NEG.out.file   // channel: [ val(meta), [ bedgraph ] ]

}



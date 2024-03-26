#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
      LIST OF PARAMETERS
================================
            GENERAL
projectDir       : $projectDir
launchDir        : $launchDir
Results folder   : $params.outdir
================================
      INPUT & REFERENCES 
Reads            : $params.samplesheet
rRNA databases   : $params.ribo_databases
fasta            : $params.fasta
gtf              : $params.gtf
star_index dir   : $params.star_index
================================
      MAIN PROCESSES ENABLED?
REMOVE SPECIAL CHAR FROM READ NAME        : $params.remove_characters  
SKIP PREPARE GENOME                       : $params.skip_prepare_genome    
SKIP FASTQC                               : $params.skip_fastqc
SKIP TRIMGALORE                           : $params.skip_trimming
SKIP BARCODE EXTRACTION                   : $params.skip_umi_extract
WITH UMI DETECTION                        : $params.with_umi
READ TO DISCARD AFTER BARCODE EXTR.       : $params.umi_discard_read
MIN # TRIMMED READS SAMPLE REMOVAL        : $params.min_trimmed_reads
DO SORTMERNA                              : $params.remove_ribo_rna
SKIP RIBODETECTOR                         : $params.skip_ribodetector
SKIP FASTQC AFTER rRNA FILTERING          : $params.skip_fastqc_after_ribo_removal
SKIP ALIGNMENT (AND DEDUP)                : $params.skip_alignment
ANALYSE UNMAPPED                          : $params.analyse_unmapped 
SKIP FASTQC AFTER DEDUPLICATION           : $params.skip_fastqc_after_dedup
SKIP EXTRACT CROSSLINKS                   : $params.skip_extract_crosslinks
SKIP BEDGRAPH TO BIGWIG FILE CONVERSION   : $params.skip_bedgraphtobigwig
SKIP PURECLIP                             : $params.skip_pureclip
SKIP MACS2                                : $params.skip_macs2
SKIP MOTIF ANALYSIS                       : $params.skip_motif_detection
SKIP MULTIQC                              : $params.skip_multiqc
================================
      IMPORTANT PROCESS PARAMS          
umitools_bc_pattern                       : $params.umitools_bc_pattern
umitools_umi_separator                    : $params.umitools_umi_separator
================================
      SAVE INTERMEDIATE FILES?
save_edited_reads                         : $params.save_edited_reads
save_umi_intermeds (UMI-Tools)            : $params.save_umi_intermeds
save_trimmed (TrimGalore)                 : $params.save_trimmed
save ribo and non_ribo reads (SortMeRNA)  : $params.save_non_ribo_reads
save_reference (STAR)                     : $params.save_reference
save_unaligned (STAR)                     : $params.save_unaligned
save_align_intermeds (STAR)               : $params.save_align_intermeds
save_extract_crosslink_intermeds (BEDTools) : $params.save_extract_crosslink_intermeds
save_intermediate_peak_id (Bash)          : $params.save_intermediate_peak_id
save_binding_width_intermeds (BEDTools)   : $params.save_binding_width_intermeds

NOTE: In its current state the pipeline may NOT work entirely if some processes are not run. Still needs to be tested.
"""

// modules adapted from nf-core and already existing pipelines (rnaseq, clipseq etc)
// modules from existing pipelines are loose in subworkflows/nf-core 
include { INPUT_CHECK                           } from './subworkflows/local/input_check.nf'
include { FASTQC as FASTQC_AFTER_SORTMERNA      } from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_AFTER_DEDUP          } from './modules/nf-core/fastqc/main.nf'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE      } from './subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main.nf'
include { FASTQC as FASTQC_UNMAPPED             } from './modules/nf-core/fastqc/main.nf'
include { PREPARE_GENOME                        } from './subworkflows/local/prepare_genome.nf'
include { ALIGN_STAR                            } from './subworkflows/local/align_star.nf'
include { SORTMERNA                             } from './modules/local/sortmerna/main.nf'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_PURECLIP          } from './modules/local/homer/annotatepeaks/main.nf'
include { HOMER_ANNOTATEPEAKS as HOMER_ANNOTATEPEAKS_MACS2              } from './modules/local/homer/annotatepeaks/main.nf'
include { MACS2_CALLPEAK                        } from './modules/nf-core/macs2/callpeak/main.nf'
include { LINUX_COMMAND as UNIQUE_PEAK_NAME     } from './modules/goodwright/linux/command/main.nf'

// local modules
include { PURECLIP                              } from './modules/local/pureclip.nf'
include { MULTIQC                               } from './modules/local/multiqc/main.nf'
include { RIBODETECTOR                          } from './modules/local/ribodetector.nf' // not used at the moment
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_PURECLIP } from './modules/local/bedtools_intersect.nf' // adapted from nf-core
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_MACS2 } from './modules/local/bedtools_intersect.nf' // adapted from nf-core
include { STREME                                } from './modules/local/streme.nf'

// nf-core subworkflows only
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from './subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main.nf'

// local subworkflows
include { REMOVE_CHARACTERS                     } from './subworkflows/local/remove_characters.nf'
include { EXTRACT_CROSSLINKS                    } from './subworkflows/local/extract_crosslinks.nf'
include { RESIZE_SITES                          } from './subworkflows/local/resize_sites.nf'

// neg and pos separately
include { UCSC_BEDGRAPHTOBIGWIG as BEDGRAPHTOBIGWIG_RAW_POS }   from './modules/local/ucsc/bedgraphtobigwig/main.nf'
include { UCSC_BEDGRAPHTOBIGWIG as BEDGRAPHTOBIGWIG_NORM_POS }  from './modules/local/ucsc/bedgraphtobigwig/main.nf'
include { UCSC_BEDGRAPHTOBIGWIG as BEDGRAPHTOBIGWIG_RAW_NEG }   from './modules/local/ucsc/bedgraphtobigwig/main.nf'
include { UCSC_BEDGRAPHTOBIGWIG as BEDGRAPHTOBIGWIG_NORM_NEG }  from './modules/local/ucsc/bedgraphtobigwig/main.nf'



workflow {

      // check if mandatory input path parameters exist
      check_params = [
            samplesheet: params.samplesheet,    // sheet with sample information
            fasta: params.fasta,                // reference genome
            gtf: params.gtf                     // gene annotation
      ]
      for (param in check_params){
            if (!param.value) {
                  exit 1, "Missing mandatory parameter: ${param.key}"
            } 
            else {
                  file(param.value, checkIfExists: true)
            }
      }

      // empty channel for versions
      ch_versions            = Channel.empty()

      /*
            STEP 0: load data 
      */

      // SUBWORKFLOW INPUT_CHECK from nf-core/rnaseq
      // checks the samplesheet for correct format, extracts meta map and read path to channel
      INPUT_CHECK ( params.samplesheet )
      .reads // get emit: reads
      .set { ch_reads }
      ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

      /*
            STEP 1: preprocessing 
      */

      // SUBWORKFLOW (local)
      // replace whitespaces and special characters
      if (params.remove_characters) {
            REMOVE_CHARACTERS ( ch_reads )
            .reads
            .set { ch_reads }
            ch_versions = ch_versions.mix(REMOVE_CHARACTERS.out.versions)
      }

      // SUBWORKFLOW FASTQC + UMITOOLS + TRIMGALORE (CUTADAPT)
      // does fastqc + trimgalore (cutadapt) + umitools to remove adapters, relocate barcodes, low quality 3' ends
      ch_filtered_reads      = Channel.empty()
      ch_fastqc_raw_multiqc  = Channel.empty()
      ch_fastqc_trim_multiqc = Channel.empty()
      ch_trim_log_multiqc    = Channel.empty()
      ch_trim_read_count     = Channel.empty()

      FASTQ_FASTQC_UMITOOLS_TRIMGALORE ( 
            ch_reads, 
            params.skip_fastqc, 
            params.with_umi, 
            params.skip_umi_extract, 
            params.skip_trimming, 
            params.umi_discard_read, 
            params.min_trimmed_reads )
      .reads
      .set { ch_trimmed_reads } 

      // multiqc channels
      ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
      ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
      ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
      ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
      ch_versions            = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)


      // MODULE SORTMERNA from nf-core and nf-core/rnaseq
      // filter out rRNA
      ch_sortmerna_multiqc          = Channel.empty()
      ch_fastqc_filtered_multiqc = Channel.empty()
      ch_filtered_reads             = ch_trimmed_reads          // if rna filtering is skipped

      if (params.remove_ribo_rna) {
            // check for rRNA databases for sortmerna
            ch_ribo_db = file(params.ribo_databases, checkIfExists: true)
            if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
            // put databases in a channel
            ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()

            SORTMERNA ( ch_trimmed_reads, ch_sortmerna_fastas ) // run sortmerna
            .reads // get emit: reads
            .set {ch_filtered_reads} // rename read data (after rRNA filtering)
            
            ch_sortmerna_multiqc    = SORTMERNA.out.log
            ch_versions             = ch_versions.mix(SORTMERNA.out.versions.first())

            // MODULE FASTQC
            // do quality control on data after rRNA removal
            if (!params.skip_fastqc_after_ribo_removal) {
                  FASTQC_AFTER_SORTMERNA ( ch_filtered_reads )

                  ch_fastqc_filtered_multiqc    = FASTQC_AFTER_SORTMERNA.out.zip
                  ch_versions                   = ch_versions.mix(FASTQC_AFTER_SORTMERNA.out.versions.first())
            }
      }


      // SUBWORKFLOW STAR GENERATE GENOME
      // STAR genome indexing
      if (!params.skip_prepare_genome) {
            PREPARE_GENOME (
            params.fasta,
            params.gtf,
            params.star_index
            )
            ch_fasta       = PREPARE_GENOME.out.fasta             // channel: [ fasta ]
            ch_chrom_sizes = PREPARE_GENOME.out.chrom_sizes       // channel: [ chrom_sizes ]
            ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
      }

      /*
            STEP 2: post-processing 
      */

      // SUBWORKFLOW STAR ALIGN adapted from nf-core/rnaseq
      // STAR alignment

      // empty channels
      ch_genome_bam                 = Channel.empty()
      ch_genome_bam_index           = Channel.empty()
      ch_samtools_stats             = Channel.empty()
      ch_samtools_flagstat          = Channel.empty()
      ch_samtools_idxstats          = Channel.empty()
      ch_star_multiqc               = Channel.empty()
      ch_star_log_multiqc           = Channel.empty()
      ch_star_gene_multiqc          = Channel.empty()
      ch_unmapped                   = Channel.empty()
      ch_unmapped_multiqc           = Channel.empty()

      if (!params.skip_alignment) {
            ALIGN_STAR (
                  ch_filtered_reads,
                  PREPARE_GENOME.out.star_index.map { [ [:], it ] },
                  PREPARE_GENOME.out.gtf.map { [ [:], it ] },
                  params.star_ignore_sjdbgtf,
                  '',
                  params.seq_center ?: '',
                  PREPARE_GENOME.out.fasta.map { [ [:], it ] }
            )

            ch_genome_bam        = ALIGN_STAR.out.bam
            ch_genome_bam_index  = ALIGN_STAR.out.bai
            ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
            ch_samtools_stats    = ALIGN_STAR.out.stats
            ch_samtools_flagstat = ALIGN_STAR.out.flagstat
            ch_samtools_idxstats = ALIGN_STAR.out.idxstats
            ch_star_log_multiqc  = ALIGN_STAR.out.log_final
            ch_star_gene_multiqc = ALIGN_STAR.out.tab
            if (params.bam_csi_index) {
                  ch_genome_bam_index = ALIGN_STAR.out.csi
            }
            ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)


            // MODULE FASTQC from nf-core
            // analyse unmapped reads
            if (params.analyse_unmapped) {
                  // gather unmapped in channel
                  ch_unmapped         = ALIGN_STAR.out.fastq
                  FASTQC_UNMAPPED ( ch_unmapped )

                  ch_unmapped_multiqc = FASTQC_UNMAPPED.out.zip
                  ch_versions         = ch_versions.mix(FASTQC_UNMAPPED.out.versions.first())
            }


            // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
            // deduplication using UMI-Tools
            if (params.with_umi) {
                  // deduplicate genome BAM file before downstream analysis
                  BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME (
                        ch_genome_bam.join(ch_genome_bam_index, by: [0]),
                        params.umitools_dedup_stats
                  )
                  ch_genome_bam        = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam
                  ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bai
                  ch_samtools_stats    = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.stats
                  ch_samtools_flagstat = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.flagstat
                  ch_samtools_idxstats = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.idxstats

                  if (params.bam_csi_index) {
                        ch_genome_bam_index  = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.csi
                  }
                  ch_versions = ch_versions.mix(BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.versions)
            }
      }


      // MODULE FASTQC
      // do quality control on data after deduplication
      ch_fastqc_dedup_multiqc = Channel.empty()

      if (!params.skip_fastqc_after_dedup) {
            FASTQC_AFTER_DEDUP ( ch_genome_bam )

            ch_fastqc_dedup_multiqc       = FASTQC_AFTER_DEDUP.out.zip
            ch_versions                   = ch_versions.mix(FASTQC_AFTER_DEDUP.out.versions.first())
      }


      // SUBWORKFLOW from local, but based on/adapted from goodwright/clipseq, Busch et al (2020) and nf-core
      if (!params.skip_alignment && !params.skip_extract_crosslinks) {
            EXTRACT_CROSSLINKS (
                  ch_genome_bam,
                  ch_chrom_sizes
            )
            ch_merged_bed                  = EXTRACT_CROSSLINKS.out.merged_bed
            ch_merged_coverage             = EXTRACT_CROSSLINKS.out.merged_coverage
            ch_merged_coverage_normalized  = EXTRACT_CROSSLINKS.out.merged_coverage_normalized
            ch_versions                    = ch_versions.mix(EXTRACT_CROSSLINKS.out.versions)

            ch_pos_bed                  = EXTRACT_CROSSLINKS.out.pos_bed
            ch_pos_coverage             = EXTRACT_CROSSLINKS.out.pos_coverage
            ch_pos_coverage_normalized  = EXTRACT_CROSSLINKS.out.pos_coverage_normalized

            ch_neg_bed                  = EXTRACT_CROSSLINKS.out.neg_bed
            ch_neg_coverage             = EXTRACT_CROSSLINKS.out.neg_coverage
            ch_neg_coverage_normalized  = EXTRACT_CROSSLINKS.out.neg_coverage_normalized


            // MODULE from nf-core
            if (!params.skip_bedgraphtobigwig) {
                  BEDGRAPHTOBIGWIG_RAW_POS       ( ch_pos_coverage, ch_chrom_sizes )
                  BEDGRAPHTOBIGWIG_NORM_POS      ( ch_pos_coverage_normalized, ch_chrom_sizes )
                  BEDGRAPHTOBIGWIG_RAW_NEG       ( ch_neg_coverage, ch_chrom_sizes )
                  BEDGRAPHTOBIGWIG_NORM_NEG      ( ch_neg_coverage_normalized, ch_chrom_sizes )
                  ch_versions = ch_versions.mix(BEDGRAPHTOBIGWIG_RAW_POS.out.versions)
            }
      }

      /*
            STEP 3: analysis 
      */

      if (!params.skip_pureclip) {
            PURECLIP ( 
                  ch_genome_bam,
                  ch_genome_bam_index,
                  PREPARE_GENOME.out.fasta
            )
            pureclip_crosslink_sites      = PURECLIP.out.sites_bed
            pureclip_crosslink_regions    = PURECLIP.out.regions_bed
      }

      // MACS2, a broad peak caller
      if (!params.skip_macs2) {
            MACS2_CALLPEAK ( 
                  ch_genome_bam,  
                  'hs' // add as param?
            )
            ch_macs2_peaks    = MACS2_CALLPEAK.out.peak 
            ch_versions       = ch_versions.mix(MACS2_CALLPEAK.out.versions.first())
      }

      // resize binding sites for motif analysis and POSSIBLY gene annotation
      RESIZE_SITES (
            pureclip_crosslink_sites,
            ch_chrom_sizes,
            ch_fasta
      )

      // if true, uses resized PureCLIP binding sites for gene annotation
      if (params.resized_for_annotation) {
            pureclip_crosslink_sites = RESIZE_SITES.out.resized_bed
      }


      // bedtools intersection of bed file with gtf
      if (!params.skip_bedtools_annotation) {
            if (!params.skip_pureclip) {
                  // replace column 4 with unique peak names (macs2 already has these IDs)
                  UNIQUE_PEAK_NAME (
                        pureclip_crosslink_sites, [], false
                  )
                  pureclip_crosslink_sites = UNIQUE_PEAK_NAME.out.file

                  BEDTOOLS_INTERSECT_PURECLIP (
                        'pureclip',
                        pureclip_crosslink_sites,
                        PREPARE_GENOME.out.gtf.map { [ [:], it ] }
                  )
                  ch_versions       = ch_versions.mix(BEDTOOLS_INTERSECT_PURECLIP.out.versions.first())
            }

            if (!params.skip_macs2) {
                  BEDTOOLS_INTERSECT_MACS2 (
                        'macs2',
                        ch_macs2_peaks,
                        PREPARE_GENOME.out.gtf.map { [ [:], it ] }
                  )
                  ch_versions       = ch_versions.mix(BEDTOOLS_INTERSECT_MACS2.out.versions.first())
            }
      }

      // homer annotatepeaks
      if (!params.skip_homer_annotation) {
            if (!params.skip_pureclip){
                  HOMER_ANNOTATEPEAKS_PURECLIP (
                        'pureclip',
                        pureclip_crosslink_sites,
                        PREPARE_GENOME.out.fasta,
                        PREPARE_GENOME.out.gtf
                  )
                  ch_versions       = ch_versions.mix(HOMER_ANNOTATEPEAKS_PURECLIP.out.versions.first())
            }

            if (!params.skip_macs2) {
                  HOMER_ANNOTATEPEAKS_MACS2 (
                        'macs2',
                        ch_macs2_peaks,
                        PREPARE_GENOME.out.fasta,
                        PREPARE_GENOME.out.gtf
                  )
                  ch_versions       = ch_versions.mix(HOMER_ANNOTATEPEAKS_MACS2.out.versions.first())
            }
      }

      // perform motif analysis using STREME, generates a .html report (and optionally .tsv, .txt, .xml files)
      // uses pureclip results as pureclip is designed for iCLIP datasets specifically
      if (!params.skip_motif_detection){
            STREME (
            RESIZE_SITES.out.peak_fasta
            )
      }

      // MODULE MULTIQC
      // make multiqc report
      if (!params.skip_multiqc) {

            // config files
            ch_multiqc_config          = Channel.fromPath("$projectDir/multiqc_config.yml", checkIfExists: true )
            ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
            ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()

            MULTIQC (
                  // ch_multiqc_files.collect(),
                  ch_multiqc_config.toList(),
                  ch_multiqc_custom_config.toList(),
                  ch_multiqc_logo.toList(),

                  // paths
                  // raw fastqc
                  ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),

                  // trimming
                  ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
                  ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),

                  // filtering rRNA
                  ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),
                  ch_fastqc_filtered_multiqc.collect{it[1]}.ifEmpty([]),
                  
                  // alignment and dedup
                  ch_star_log_multiqc.collect{it[1]}.ifEmpty([]),
                  ch_star_gene_multiqc.collect{it[1]}.ifEmpty([]),
                  ch_samtools_stats.collect{it[1]}.ifEmpty([]),
                  ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
                  ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
                  ch_fastqc_dedup_multiqc.collect{it[1]}.ifEmpty([]),
                  
                  // analyse unmapped
                  ch_unmapped_multiqc.collect{it[1]}.ifEmpty([])

            ) // run multiqc
      }

}


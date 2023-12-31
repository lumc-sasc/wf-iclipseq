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
      PROCESSES ENABLED
REMOVE SPECIAL CHAR FROM READ NAME        : $params.remove_characters  
SKIP PREPARE GENOME                       : $params.skip_prepare_genome    
SKIP FASTQC                               : $params.skip_fastqc
SKIP TRIMGALORE                           : $params.skip_trimming
SKIP BARCODE EXTRACTION                   : $params.skip_umi_extract
WITH UMI DETECTION                        : $params.with_umi
READ TO DISCARD AFTER BARCODE EXTR.       : $params.umi_discard_read
MIN # TRIMMED READS SAMPLE REMOVAL        : $params.min_trimmed_reads
SKIP PRE-MAPPING                          : $params.skip_premapping
DO SORTMERNA                              : $params.remove_ribo_rna
SKIP FASTQC AFTER rRNA FILTERING          : $params.skip_fastqc_after_ribo_removal
SKIP ALIGNMENT (AND DEDUP)                : $params.skip_alignment
ANALYSE UNMAPPED                          : $params.analyse_unmapped 
SKIP FASTQC AFTER DEDUPLICATION           : $params.skip_fastqc_after_dedup
SKIP EXTRACT CROSSLINKS                   : $params.skip_extract_crosslinks
SKIP BEDGRAPH TO BIGWIG FILE CONVERSION   : $params.skip_bedgraphtobigwig
SKIP PURECLIP                             : $params.skip_pureclip
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
save_non_ribo_reads (SortMeRNA)           : $params.save_non_ribo_reads
save_reference (STAR)                     : $params.save_reference
save_unaligned (STAR)                     : $params.save_unaligned
save_align_intermeds (STAR)               : $params.save_align_intermeds
save_extract_crosslink_intermeds (BEDTools) : $params.save_extract_crosslink_intermeds

NOTE: In its current state the pipeline may NOT work entirely if some processes are not run. Still needs to be tested.
"""

// modules adapted from nf-core and already existing pipelines (rnaseq, etc)
// modules from existing pipelines are loose in subworkflows/nf-core 
include { INPUT_CHECK                           } from './subworkflows/nf-core/input_check.nf'
include { FASTQC as FASTQC_AFTER_SORTMERNA      } from './modules/nf-core/fastqc/main.nf'
include { FASTQC as FASTQC_AFTER_DEDUP          } from './modules/nf-core/fastqc/main.nf'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE      } from './subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main.nf'
include { FASTQC as FASTQC_UNMAPPED             } from './modules/nf-core/fastqc/main.nf'
include { PREPARE_GENOME                        } from './subworkflows/local/prepare_genome.nf'
include { BOWTIE2_BUILD                         } from './modules/nf-core/bowtie2/build/main.nf'
include { BOWTIE2_ALIGN                         } from './modules/nf-core/bowtie2/align/main.nf'
include { ALIGN_STAR                            } from './subworkflows/local/align_star.nf'
include { SORTMERNA                             } from './modules/local/sortmerna/main.nf'
include { HOMER_ANNOTATEPEAKS                   } from './modules/nf-core/homer/annotatepeaks/main.nf'

// local modules
include { PURECLIP                              } from './modules/local/pureclip.nf'
include { MULTIQC                               } from './modules/local/multiqc/main.nf'

// nf-core subworkflows only
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from './subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main.nf'

// local subworkflows
include { REMOVE_CHARACTERS                     } from './subworkflows/local/remove_characters.nf'
include { EXTRACT_CROSSLINKS                    } from './subworkflows/local/extract_crosslinks.nf'

// neg and pos separately
include { UCSC_BEDGRAPHTOBIGWIG as BEDGRAPHTOBIGWIG_RAW_POS }   from './modules/local/ucsc/bedgraphtobigwig/main.nf'
include { UCSC_BEDGRAPHTOBIGWIG as BEDGRAPHTOBIGWIG_NORM_POS }  from './modules/local/ucsc/bedgraphtobigwig/main.nf'
include { UCSC_BEDGRAPHTOBIGWIG as BEDGRAPHTOBIGWIG_RAW_NEG }   from './modules/local/ucsc/bedgraphtobigwig/main.nf'
include { UCSC_BEDGRAPHTOBIGWIG as BEDGRAPHTOBIGWIG_NORM_NEG }  from './modules/local/ucsc/bedgraphtobigwig/main.nf'


// workflow
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


      // SUBWORKFLOW INPUT_CHECK from nf-core/rnaseq
      // checks the samplesheet for correct format, extracts meta map and read path to channel
      INPUT_CHECK ( params.samplesheet )
      .reads // get emit: reads
      .set { ch_reads }
      ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


      // SUBWORKFLOW (local)
      // replace whitespaces and special characters
      if (params.remove_characters) {
            REMOVE_CHARACTERS ( ch_reads )
            .reads
            .set { ch_reads }
            ch_versions = ch_versions.mix(REMOVE_CHARACTERS.out.versions)
      }



      // SUBWORKFLOW FASTQC + UMITOOLS + TRIMGALORE (CUTADAPT)
      // does fastqc + trimgalore (cutadapt) + umitools to remove adapters, barcodes, low quality 3' ends
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


      // BOWTIE2 MODULES 
      // premapping
      ch_multiqc_bowtie2            = Channel.empty()
      ch_filtered_reads             = ch_trimmed_reads          // if premapping is skipped
      if (!params.skip_premapping) {
            ch_premap_fasta         = Channel.of(params.bw_fasta)

            // build index
            BOWTIE2_BUILD ( ch_premap_fasta.map { [ [:], it ] } )
            ch_version              = ch_versions.mix(BOWTIE2_BUILD.out.versions)

            // aligning
            BOWTIE2_ALIGN (
                  ch_trimmed_reads,
                  BOWTIE2_BUILD.out.index,
                  true, // save_unaligned
                  true // sort_bam
            )
            
            ch_filtered_reads       = BOWTIE2_ALIGN.out.fastq
            ch_multiqc_bowtie2      = BOWTIE2_ALIGN.out.log
            ch_aligned_bowtie2      = BOWTIE2_ALIGN.out.aligned
            ch_version              = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
      }


      // MODULE SORTMERNA from nf-core and nf-core/rnaseq
      // filter out rRNA
      ch_sortmerna_multiqc          = Channel.empty()

      if (params.remove_ribo_rna) {
            // check for rRNA databases for sortmerna
            ch_ribo_db = file(params.ribo_databases, checkIfExists: true)
            if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
            // put databases in a channel
            ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()

            SORTMERNA ( ch_filtered_reads, ch_sortmerna_fastas ) // run sortmerna
            .reads // get emit: reads
            .set {ch_filtered_reads} // rename read data (after rRNA filtering)
            
            ch_sortmerna_multiqc    = SORTMERNA.out.log
            ch_versions             = ch_versions.mix(SORTMERNA.out.versions.first())
      }


      // MODULE FASTQC
      // do quality control on data after rRNA removal
      ch_fastqc_filtered_multiqc = Channel.empty()

      if (!params.skip_fastqc_after_ribo_removal) {
            FASTQC_AFTER_SORTMERNA ( ch_filtered_reads )

            ch_fastqc_filtered_multiqc    = FASTQC_AFTER_SORTMERNA.out.zip
            ch_versions                   = ch_versions.mix(FASTQC_AFTER_SORTMERNA.out.versions.first())
      }


      // SUBWORKFLOW STAR GENERATE GENOME
      // STAR genome indexing
      if (!params.skip_prepare_genome) {
            PREPARE_GENOME (
            params.fasta,
            params.gtf,
            params.star_index
            )
            ch_chrom_sizes = PREPARE_GENOME.out.chrom_sizes
            ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
      }


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
      ch_aligner_pca_multiqc        = Channel.empty() // not used rn
      ch_aligner_clustering_multiqc = Channel.empty() // not used rn
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
            // deduplication using UMI-Tools (include mapped to multiple loci!)
            if (params.with_umi) {
                  // deduplicate genome BAM file before downstream analysis
                  BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME ( // TODO: does samtools manipulate the data?
                        ch_genome_bam.join(ch_genome_bam_index, by: [0]), // joining by which element?
                        params.umitools_dedup_stats // still cannot see these in the multiqc
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
      // TODO: extract crosslinks 
      if (!params.skip_alignment && !params.skip_extract_crosslinks) {
            EXTRACT_CROSSLINKS (
                  BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME.out.bam,
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
            // bigwig files can be used as input for IGV
            // https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/CNVVisualization%20and%20BigWigFiles_fDG_dtSW.htm
            if (!params.skip_bedgraphtobigwig) {
                  BEDGRAPHTOBIGWIG_RAW_POS       ( ch_pos_coverage, ch_chrom_sizes )
                  BEDGRAPHTOBIGWIG_NORM_POS      ( ch_pos_coverage_normalized, ch_chrom_sizes )
                  BEDGRAPHTOBIGWIG_RAW_NEG       ( ch_neg_coverage, ch_chrom_sizes )
                  BEDGRAPHTOBIGWIG_NORM_NEG      ( ch_neg_coverage_normalized, ch_chrom_sizes )
                  ch_versions = ch_versions.mix(BEDGRAPHTOBIGWIG_RAW_POS.out.versions)
            }
      }

      // TODO: library complexity

      // TODO?: check insertions and deletions (Busch et al) 

      // TODO: iCLIPro


      // TODO: peak calling
      // pureclip 
      if (!params.skip_pureclip) {
            PURECLIP ( 
                  ch_genome_bam,
                  ch_genome_bam_index,
                  PREPARE_GENOME.out.fasta
            )

            pureclip_crosslink_sites      = PURECLIP.out.sites_bed
            pureclip_crosslink_regions    = PURECLIP.out.regions_bed
      
      }

      // peakachu?
      if (!params.skip_peakachu) {
            placeholder = Channel.from(1,2)
      }

      // maybe paraclu


      


      // possibly clippy and PEKA



      // add dreme

      if (!params.skip_homer_annotation) {
      // homer annotatepeaks
            HOMER_ANNOTATEPEAKS (
                  pureclip_crosslink_regions, // peak regions
                  PREPARE_GENOME.out.fasta,
                  PREPARE_GENOME.out.gtf
            )
      }

      // TODO: further analysis





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
                  ch_unmapped_multiqc.collect{it[1]}.ifEmpty([]),

                  // premapping
                  ch_multiqc_bowtie2.collect{it[1]}.ifEmpty([])

            ) // run multiqc
      }

}


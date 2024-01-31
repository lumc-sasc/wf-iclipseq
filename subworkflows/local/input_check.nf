// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    meta.control_bam  = ' '
    meta.control_bai  = ' '

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }

    // add path(s) of the control file(s) to the meta map
    if (row.control_bam != '' && row.control_bai != ''){
        if (file(row.control_bam).exists() && file(row.control_bai).exists()) {
            meta.control_bam = file(row.control_bam)
            meta.control_bai = file(row.control_bai)
        } 
    }
        // fastq_meta[1].add(file(row.control_bam))
        // fastq_meta.merge(file(row.control_bam))
        // meta.control = 1

        // if (file(row.control_bai).exists()) {
        //     if (!file(row.control_bam).exists()) {
        //         exit 1, "ERROR: Could only find a .bai file. Make sure to add the bam file as well."    
        //     } else {
        //         meta.control_bai = file(row.control_bai)
        //         // fastq_meta[1].add(file(row.control_bai))\
        //         // fastq_meta.merge(file(row.control_bam))
        //         // meta.control = 2
        //     }
        // }
    //}

    return fastq_meta
}

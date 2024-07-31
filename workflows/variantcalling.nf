/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowVariantcalling.initialise(params, log)

// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.fasta, params.fai, params.interval, params.include_positions, params.exclude_positions ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
//if (params.input) { ch_input = Channel.fromPath(params.input) } else { exit 1, 'Input samplesheet not specified!' }
//if (params.fasta) { ch_fasta = Channel.fromPath(params.fasta) } else { exit 1, 'Reference fasta not specified!'   }

// Check optional parameters
//if (params.fai){
//    if( ( params.fasta.endsWith('.gz') && params.fai.endsWith('.fai') )
//       ||
//        ( !params.fasta.endsWith('.gz') && params.fai.endsWith('.gzi') )
//    ){
//      exit 1, 'Reference fasta and its index file format not matched!'
//    }
//    ch_fai = Channel.fromPath(params.fai)
//} else { 
//    ch_fai = Channel.empty()
//}

//if (params.interval){ ch_interval = Channel.fromPath(params.interval) } else { ch_interval = Channel.empty() }

//if (params.split_fasta_cutoff ) { split_fasta_cutoff = params.split_fasta_cutoff } else { split_fasta_cutoff = 100000 }

//if ( (params.include_positions) && (params.exclude_positions) ){
//    exit 1, 'Only one positions file can be given to include or exclude!'
//}else if (params.include_positions){ 
//    ch_positions = Channel.fromPath(params.include_positions) 
//} else if (params.exclude_positions){
//    ch_positions = Channel.fromPath(params.exclude_positions) 
//} else { 
//    ch_positions = [] 
//}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES:
//
include { BWA_INDEX          } from '../modules/local/bwaIndex'
include { BWA_MEM            } from '../modules/local/bwaMem'
include { SORT_BAM_SAMBAMBA   } from '../modules/local/sortBamSambamba'
include { MARK_DUPES_SAMBAMBA } from '../modules/local/markDupesSambamba'
include { INDEX_BAM_SAMBAMBA  } from '../modules/local/indexBamSambamba'
include { MOSDEPTH          } from '../modules/local/mosdepth'
include { INTERSECT_BEDS     } from '../modules/local/intersectBeds'
include { MERGE_SAMBAMBA     } from '../modules/local/mergeSambamba'
include { FREEBAYES         } from '../modules/local/freebayes'
include { BCFTOOLS_FILTER    } from '../modules/local/bcftoolsFilter'
include { GENERATE_FAIL_BED   } from '../modules/local/generateFailBed'
include { GENERATE_PASS_VCF   } from '../modules/local/generatePassVCF'
include { BEDTOOLS_SUBTRACT  } from '../modules/local/bedtoolsSubtract'
include { BCFTOOLS_SORT      } from '../modules/local/bcftoolsSort'
include { BCFTOOLS_INDEX     } from '../modules/local/bcftoolsIndex'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTCALLING {

     BWA_INDEX(params.genome)

    BWA_MEM(params.genome, 
        BWA_INDEX.out, 
        params.read_files
    )

    SORT_BAM_SAMBAMBA(BWA_MEM.out)

    MARK_DUPES_SAMBAMBA(SORT_BAM_SAMBAMBA.out)

    INDEX_BAM_SAMBAMBA(MARK_DUPES_SAMBAMBA.out.meta_bam)

    MOSDEPTH(MARK_DUPES_SAMBAMBA.out.meta_bam.join(INDEX_BAM_SAMBAMBA.out), 
        8
    )

    INTERSECT_BEDS(MOSDEPTH.out.collect(), 
        params.species
    )

    MERGE_SAMBAMBA(MARK_DUPES_SAMBAMBA.out.bam_only.collect(), 
        params.species
    )

    FREEBAYES(params.genome, 
        params.genome_index, 
        MERGE_SAMBAMBA.out, 
        INTERSECT_BEDS.out.overlap
    )

    BCFTOOLS_FILTER(params.genome, 
        FREEBAYES.out,
    )

    GENERATE_FAIL_BED(BCFTOOLS_FILTER.out)

    GENERATE_PASS_VCF(BCFTOOLS_FILTER.out)
    
    BEDTOOLS_SUBTRACT(INTERSECT_BEDS.out.overlap, 
        GENERATE_FAIL_BED.out
    )
    
    BCFTOOLS_SORT(GENERATE_PASS_VCF.out)

    BCFTOOLS_INDEX(BCFTOOLS_SORT.out)
}

workflow VARIANTCALLING_SANGER {

    ch_versions = Channel.empty()
    ch_fasta
     .map { fasta -> [ [ 'id': fasta.baseName -  ~/.fa\w*$/ ], fasta ] }
     .first()
     .set { ch_genome }

    //
    // check reference fasta index given or not
    //
    if( params.fai == null ){ 
   
       SAMTOOLS_FAIDX ( ch_genome,  [[], []] )
       ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )

       if( params.fasta.endsWith('.gz') ){
            ch_genome_index = SAMTOOLS_FAIDX.out.gzi
       }else{
            ch_genome_index = SAMTOOLS_FAIDX.out.fai
       }

    }else{
       ch_fai
        .map { fai -> [ [ 'id': fai.baseName ], fai ] }
        .first()
        .set { ch_genome_index }
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix( INPUT_CHECK.out.versions )


    //
    // SUBWORKFLOW: align reads if required
    //
    if( params.align ){

        if ( params.vector_db.endsWith( '.tar.gz' ) ) {

            UNTAR ( [ [:], params.vector_db ] ).untar
            | map { meta, file -> file }
            | set { ch_vector_db }
            ch_versions = ch_versions.mix ( UNTAR.out.versions )


        } else {

            Channel.fromPath ( params.vector_db )
            | set { ch_vector_db }

        }

        ALIGN_PACBIO (
            ch_genome,
            INPUT_CHECK.out.reads,
            ch_vector_db
        )
       ch_versions = ch_versions.mix( ALIGN_PACBIO.out.versions )

       ALIGN_PACBIO.out.cram
        .join( ALIGN_PACBIO.out.crai )
        .set{ ch_aligned_reads }

    } else {

        //
        // SUBWORKFLOW: merge the input reads by sample name
        //
        INPUT_MERGE (
            ch_genome,
            ch_genome_index,
            INPUT_CHECK.out.reads,
        )
        ch_versions = ch_versions.mix( INPUT_MERGE.out.versions )
        ch_aligned_reads = INPUT_MERGE.out.indexed_merged_reads

    }

    //
    // SUBWORKFLOW: split the input fasta file and filter input reads
    //
    INPUT_FILTER_SPLIT (
        ch_fasta,
        ch_aligned_reads,
        ch_interval,
        split_fasta_cutoff
    )
    ch_versions = ch_versions.mix( INPUT_FILTER_SPLIT.out.versions )


    //
    // SUBWORKFLOW: call deepvariant
    //
    DEEPVARIANT_CALLER (
        INPUT_FILTER_SPLIT.out.reads_fasta
    )
    ch_versions = ch_versions.mix( DEEPVARIANT_CALLER.out.versions )


    //
    // convert VCF channel meta id 
    // 
    DEEPVARIANT_CALLER.out.vcf 
     .map{ meta, vcf -> [ [ id: vcf.baseName ], vcf ] }
     .set{ vcf }

    //
    // process VCF output files
    //
    PROCESS_VCF( vcf, ch_positions )
    ch_versions = ch_versions.mix( PROCESS_VCF.out.versions )


    //
    // MODULE: Combine different version together
    // 
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//workflow.onComplete {
//    if (params.email || params.email_on_fail) {
//        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
//    }
//    NfcoreTemplate.summary(workflow, params, log)
//    if (params.hook_url) {
//        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
//    }
//}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_vcf = ''
params.output = params.output ?: "${workflow.launchDir}"
params.out_root = "${params.output}/PathogenicSplicingScreening" + new Date().format('yyyyMMdd-HHmmss')

include { SPLICEAI; VEP } from './modules/processes.nf'

workflow {
    Channel.fromPath(params.input_vcf)
        | map { it -> 
                tuple(it, "${it}.tbi", "${params.reference}", "${params.annotation_gtf}") }
        | view
        | SPLICEAI
        | map { it -> tuple(it, "${params.reference}") }
        | view
        // | VEP
}

workflow.onComplete {
      println ""
      println "~*~*~*~*~~*~*~*~*~~*~*~*~*~~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~"
      println "Pipeline completed at: $workflow.complete"
      println "Execution time       : $workflow.duration"
      println "Execution status     : ${ workflow.success ? 'OK' : 'failed' }"
      println "~*~*~*~*~~*~*~*~*~~*~*~*~*~~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~"
}

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_vcf = ''

process SPLICEAI {
    container 'betelgeuse:5000/library/utsu/spliceai:1.3.1'
    
    input:
        tuple path(input_vcf), path(input_tbi), path(reference_fasta), path(annotation_gtf)

    output:
        path 'splai.vcf'

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh && \\
    conda activate spliceai && \\
    spliceai \\
      -I ${input_vcf} \\
      -O splai.vcf \\
      -R ${reference_fasta} \\
      -A ${annotation_gtf} \\
      -D 4999 \\
      -M 0
    """
}

process VEP {
    container 'betelgeuse:5000/library/utsu/vep:112.2'
    containerOptions "-u 0 -v ${params.vep_data}:/data -v ${params.vep_plugin_resources}:/plugin_resources"

    input: 
        tuple path(input_vcf), path(reference_fasta)
    
    output:
        path 'splai.vep.vcf' 

    script:
    """
    /opt/vep/src/ensembl-vep/vep \\
      --dir_cache /data \\
      --cache \\
      --offline \\
      --no_stats \\
      --gencode_basic \\
      --variant_class \\
      --canonical \\
      --symbol \\
      --numbers \\
      --vcf \\
      --pick_allele \\
      --force_overwrite \\
      --use_given_ref \\
      --assembly ${params.assembly} \\
      --fasta ${reference_fasta} \\
      --plugin MaxEntScan,/plugin_resources/maxentscan/fordownload \\
      --plugin LoF,loftee_path:/vep/loftee,human_ancestor_fa:/plugin_resources/loftee/human_ancestor.fa.gz,conservation_file:/plugin_resources/loftee/phylocsf_gerp.sql \\
      -i ${input_vcf} \\
      -o splai.vep.vcf \\
      
    """
}

// process CALCULATE_PS {
//     container 'utsuno/prioritize:latest'

//     input:
//         tuple path(input_vcf), path(reference_fasta), path(annotation_gtf)

//     output:
//         path 'prioritize.tsv'

//     script:
//     """
//     conda activate prioritize; \\
//     prioritize \\
//     -I ${input_vcf} \\
//     -R ${reference_fasta} \\
//     -A ${annotation_gtf} \\
//     -O prioritize.tsv
//     """
// }

// workflow.onComplete {
//       println ""
//       println "~*~*~*~*~~*~*~*~*~~*~*~*~*~~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~"
//       println "Pipeline completed at: $workflow.complete"
//       println "Execution time       : $workflow.duration"
//       println "Execution status     : ${ workflow.success ? 'OK' : 'failed' }"
//       println "~*~*~*~*~~*~*~*~*~~*~*~*~*~~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~"
// }

workflow {
    Channel.fromPath(params.input)
      | map { it -> tuple(it, "${it}.tbi", "${params.reference_fasta}", "${params.annotation_gtf}") }
      | SPLICEAI
      | map { it -> tuple(it, "${params.reference_fasta}") }
      | VEP
      | view()
}

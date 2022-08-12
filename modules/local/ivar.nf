process IVAR_CONSENSUS {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*.consensus.fa"
    tag "${sample}"
    container params.container_ivar
    memory {2.GB * task.attempt}
    errorStrategy {'retry'}
    maxRetries 2

    input:
        tuple val(sample), file(bamfile)

    output:
        tuple val(sample), path("${task.process}/${sample}.consensus.fa"), emit: consensus
        path("${task.process}/${sample}.consensus.fa"), emit: consensus_collect
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        samtools --version >> $log_file
        ivar --version >> $log_file

        samtools mpileup -aa -A -Q 0 !{bamfile} 2>> $err_file | \
        ivar consensus -p !{task.process}/!{sample}.consensus \
            -t 0.6 \
            -m !{params.mincov} \
            -i !{sample} \
            2>> $err_file >> $log_file
    '''
}

process IVAR_VARIANTS {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*.variants.tsv"
    tag "${sample}"
    container params.container_ivar
    memory {2.GB * task.attempt}
    errorStrategy {'retry'}
    maxRetries 2

    input:
        tuple val(sample), file(bamfile)
        file(ref_genome)
        file(ref_gff)

    output:
        path("${task.process}/${sample}.variants.tsv"), emit: ivar_tsvs
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        samtools --version >> $log_file

        samtools mpileup -aa -A -B -Q 0 --reference !{ref_genome} !{bamfile} 2>> $err_file | \
            ivar variants -t 0.05 \
            -m !{params.mincov} \
            -r !{ref_genome} \
            -g !{ref_gff} \
            -p !{task.process}/!{sample}.variants \
            2>> $err_file >> $log_file
    '''
}

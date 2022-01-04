params.pangolin_options = ''
process PANGOLIN {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*/lineage_report.csv"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_pangolin

    input:
        tuple val(sample), file(fasta)

    output:
        path("${task.process}/${sample}/lineage_report.csv"), emit: pangolin
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        pangolin --version >> $log_file
        pangolin --pangoLEARN-version >> $log_file
        pangolin --pango-designation-version >> $log_file

        pangolin !{params.pangolin_options} \
            --outdir !{task.process}/!{sample} \
            --threads !{task.cpus} \
            !{fasta} \
            2>> $err_file >> $log_file
    '''
}

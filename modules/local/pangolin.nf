params.pangolin_options = '--max-ambig 0.5 --min-length 10000'
process PANGOLIN {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/lineage_report.csv"
    echo false
    cpus params.medcpus
    container params.container_pangolin

    input:
        file(fasta)

    output:
        path("${task.process}/lineage_report.csv"), emit: pangolin
        path("logs/${task.process}/${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        pangolin --all-versions >> $log_file

        pangolin !{params.pangolin_options} \
            --outdir !{task.process}/ \
            --threads !{task.cpus} \
            !{fasta} \
            2>> $err_file >> $log_file
    '''
}

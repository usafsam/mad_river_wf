process VADR {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*"
    cpus params.medcpus
    container params.container_vadr

    input:
        file(fasta)
    
    output:
        path("${task.process}/*"), optional: true, emit: vadr
        path("${task.process}/${file(params.outdir).getSimpleName()}/${file(params.outdir).getSimpleName()}.vadr.sqa"), optional: true, emit: vadr_file
        path("logs/${task.process}/${workflow.sessionId}.{log,err}")
    
    shell:
    '''
        mkdir -p logs/!{task.process}
        log_file=logs/!{task.process}/!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        v-annotate.pl -h | head -n 2 | tail -n 1 >> $log_file

        mkdir '!{task.process}'

        fasta-trim-terminal-ambigs.pl \
            !{params.vadr_trim_options} \
            !{fasta} > '!{task.process}'/trimmed_!{fasta} 2>> $err_file

        cd '!{task.process}'

        if [[ -s trimmed_!{fasta} ]]
        then
            v-annotate.pl !{params.vadr_options} \
                --cpu !{task.cpus} \
                --mkey !{params.vadr_reference} \
                --mdir !{params.vadr_mdir} \
                trimmed_!{fasta} \
                "$(basename '!{params.outdir}')" \
                2>> ../$err_file >> ../$log_file
        fi
    '''
}

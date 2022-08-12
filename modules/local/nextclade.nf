process GRAB_NEXTCLADE_DATA {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    container params.container_nextclade

    output:
        path("data/${params.nextclade_reference_name}/reference.fasta"), emit: reference_genome
        path("data/${params.nextclade_reference_name}/"), emit: reference_nextclade

    shell:
    '''
        mkdir -p logs/
        log_file=logs/!{task.process}.log
        err_file=logs/!{task.process}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        nextclade --version >> $log_file

        nextclade dataset get \
            --name='!{params.nextclade_reference_name}' \
            --output-dir='data/!{params.nextclade_reference_name}'
    '''
}

process NEXTCLADE {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/nextclade.*"
    cpus params.medcpus
    container params.container_nextclade

    input:
        file(fasta)
        file(ref_nextclade)

    output:
        path("${task.process}/nextclade.csv"), emit: nextclade_csv
        path("${task.process}/nextclade.json"), emit: nextclade_json
        path("${task.process}/nextclade.auspice.json"), emit: nextclade_auspice_json
        path("${task.process}/nextclade.*"), emit: nextclade_out_files
        path("logs/${task.process}/${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        nextclade --version >> $log_file

        nextclade run \
            --in-order \
            --jobs !{task.cpus} \
            --input-dataset !{ref_nextclade} \
            --output-all !{task.process}/ \
            --output-basename nextclade \
            !{params.nextclade_options} \
            !{fasta} \
            2>> $err_file >> $log_file
    '''
}

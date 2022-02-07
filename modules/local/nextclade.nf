process GRAB_NEXTCLADE_DATA {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    echo false
    container params.container_nextclade

    output:
        path("data/sars-cov-2_MN908947/reference.fasta"), emit: reference_genome
        path("data/sars-cov-2_MN908947/"), emit: reference_nextclade

    shell:
    '''
        mkdir -p logs/
        log_file=logs/!{task.process}.log
        err_file=logs/!{task.process}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        nextclade --version >> $log_file

        nextclade dataset get \
            --name='sars-cov-2' \
            --reference='MN908947' \
            --output-dir='data/sars-cov-2_MN908947'
    '''
}

params.nextclade_options = ''
process NEXTCLADE {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/nextclade.*"
    echo false
    cpus params.medcpus
    container params.container_nextclade

    input:
        file(fasta)
        file(ref_nextclade)
        file(primers_csv)

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
            --jobs !{task.cpus} \
            -i !{fasta} \
            --input-dataset !{ref_nextclade} \
            --input-pcr-primers !{primers_csv} \
            --output-dir !{task.process}/ \
            --output-basename nextclade \
            --output-csv !{task.process}/nextclade.csv \
            --output-json !{task.process}/nextclade.json \
            --output-tree !{task.process}/nextclade.auspice.json \
            !{params.nextclade_options} \
            2>> $err_file >> $log_file
    '''
}

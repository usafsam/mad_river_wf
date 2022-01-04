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
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_nextclade

    input:
        tuple val(sample), file(fasta), file(ref_nextclade), file(primers_csv)

    output:
        path("${task.process}/${sample}_nextclade.csv"), emit: nextclade
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        nextclade --version >> $log_file

        nextclade run \
            --jobs !{task.cpus} \
            -i !{fasta} \
            --input-dataset !{ref_nextclade} \
            --input-pcr-primers !{primers_csv} \
            --output-dir !{task.process}/ \
            --output-basename !{sample}_nextclade \
            --output-csv !{task.process}/!{sample}_nextclade.csv \
            --output-json !{task.process}/!{sample}_nextclade.json \
            --output-tree !{task.process}/!{sample}_nextclade.auspice.json \
            !{params.nextclade_options} \
            2>> $err_file >> $log_file
    '''
}

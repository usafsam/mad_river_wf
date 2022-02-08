params.vadr_options = '--split --glsearch -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn'
params.vadr_reference = 'sarscov2'
params.vadr_mdir = '/opt/vadr/vadr-models'
params.vadr_trim_options = '--minlen 50 --maxlen 30000'
process VADR {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*"
    echo false
    cpus params.medcpus
    container 'staphb/vadr:latest'

    input:
        file(fasta)
    
    output:
        path("${task.process}/*"), optional: true, emit: vadr
        path("${task.process}/vadr.vadr.sqa"), optional: true, emit: vadr_file
        path("logs/${task.process}/${workflow.sessionId}.{log,err}")
    
    shell:
    '''
        mkdir -p logs/!{task.process}
        log_file=logs/!{task.process}/!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        v-annotate.pl -h | head -n 2 | tail -n 1 >> $log_file

        fasta-trim-terminal-ambigs.pl \
            !{params.vadr_trim_options} \
            !{fasta} > trimmed_!{fasta} 2>> $err_file
        
        if [[ -s trimmed_!{fasta} ]]
        then
            v-annotate.pl !{params.vadr_options} \
                --cpu !{task.cpus} \
                --mkey !{params.vadr_reference} \
                --mdir !{params.vadr_mdir} \
                trimmed_!{fasta} \
                !{task.process} \
                2>> $err_file >> $log_file
        fi
    '''
}
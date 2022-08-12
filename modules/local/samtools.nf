process SAMTOOLS_SORT_INDEX {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_trimclip.sorted.bam"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_trimclip.sorted.bam.bai"
    tag "${sample}"
    cpus params.medcpus
    container params.container_samtools

    input:
        tuple val(sample), file(samfile)

    output:
        tuple val(sample), path("${task.process}/${sample}_trimclip.sorted.bam"), emit: bamfile
        path("${task.process}/${sample}_trimclip.sorted.bam.bai"), emit: bamfile_index
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        samtools --version >> $log_file

        samtools view -bShu !{samfile} | \
            samtools sort \
                -o !{task.process}/!{sample}_trimclip.sorted.bam \
                -@ !{task.cpus} \
                2>> $err_file >> $log_file
        samtools index \
            !{task.process}/!{sample}_trimclip.sorted.bam \
            !{task.process}/!{sample}_trimclip.sorted.bam.bai \
            2>> $err_file >> $log_file
    '''
}

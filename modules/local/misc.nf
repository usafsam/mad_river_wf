process DATA_PREP {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    container params.container_python

    input:
        file(gff_file_gzip)
        file(primer_tsv)

    output:
        path("*.gff"), emit: gff_file_ivar
        path("*.fasta"), emit: primer_fasta
        path("*.primer.bed"), emit: primer_bed
        path("*.amplicon.bed"), emit: amplicon_bed

    shell:
    '''
        mkdir -p logs/
        log_file=logs/!{task.process}.log
        err_file=logs/!{task.process}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo $(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*$//' >> $log_file

        gunzip -f !{gff_file_gzip} \
            2>> $err_file >> $log_file

        echo "" >> $log_file
        python3 --version >> $log_file

        process_primer_reference.py !{primer_tsv} \
            --schema_name SARSCoV_1200 \
            --reference $(cat $(basename !{gff_file_gzip} .gz) | sed '/^#/d' | cut -f1 | head -n1)
    '''
}

process LINEAGE_EXCEL {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_lineage.xlsx"
    conda "${workflow.projectDir}/env/performance_lineage_excel.yml"

    input:
        file(pangolin_csv)
        file(nextclade_csv)
        file(vadr_annotations)

    output:
        path("${file(params.outdir).getSimpleName()}_lineage.xlsx")
        path("logs/${task.process}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/
        log_file=logs/!{task.process}.log
        err_file=logs/!{task.process}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        python3 --version >> $log_file

        performance_lineage_excel.py $(basename !{params.outdir}) \
            lineage \
            --pangolin_summary !{pangolin_csv} \
            --nextclade_summary !{nextclade_csv} \
            --vadr_annotations !{vadr_annotations} \
            2>> $err_file >> $log_file
    '''
}

process PERFORMANCE_EXCEL {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_performance.xlsx"
    conda "${workflow.projectDir}/env/performance_lineage_excel.yml"

    input:
        file(stats_json)
        file(coverage_summary)
        file(run_info)

    output:
        path("${file(params.outdir).getSimpleName()}_performance.xlsx")
        path("logs/${task.process}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/
        log_file=logs/!{task.process}.log
        err_file=logs/!{task.process}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        python3 --version >> $log_file

        performance_lineage_excel.py $(basename !{params.outdir}) \
            performance \
            --stats_json !{stats_json} \
            --coverage_summary !{coverage_summary} \
            --run_info !{run_info} \
            --protocol !{params.protocol} \
            2>> $err_file >> $log_file
    '''
}

params.spike_gene_coverage_sample_column = 'seqName'
params.spike_gene_coverage_file_out = "xlsx"
process SPIKE_GENE_COVERAGE {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_21563-25384.${params.spike_gene_coverage_file_out}"
    conda "${workflow.projectDir}/env/nextalign_gene_coverage.yml"

    input:
        file(nextclade_csv)

    output:
        path("${file(params.outdir).getSimpleName()}_nextclade_21563-25384.${params.spike_gene_coverage_file_out}")
        path("logs/${task.process}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/
        log_file=logs/!{task.process}.log
        err_file=logs/!{task.process}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        Rscript --version 2>> $log_file

        mv !{nextclade_csv} $(basename !{params.outdir})_nextclade.csv
        nextalign_gene_coverage.R $(basename !{params.outdir})_nextclade.csv \
            21563 25384 \
            --sample_id !{params.spike_gene_coverage_sample_column} \
            --out_type !{params.spike_gene_coverage_file_out} \
            2>> $err_file >> $log_file
    '''
}

params.traffic_light_plot_file_out = "pdf"
process TRAFFIC_LIGHT_PLOT {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.summary.${params.traffic_light_plot_file_out}"
    conda "${workflow.projectDir}/env/traffic_light_plot.yml"

    input:
        file(ivar_tsv)
        file(coverage_summary)

    output:
        path("${file(params.outdir).getSimpleName()}.summary.${params.traffic_light_plot_file_out}")
        path("logs/${task.process}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/
        log_file=logs/!{task.process}.log
        err_file=logs/!{task.process}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        Rscript --version 2>> $log_file

        traffic_light_plot.R $(basename !{params.outdir}) \
            --ivar_tsv !{ivar_tsv} \
            --cov !{coverage_summary} \
            --out_type !{params.traffic_light_plot_file_out} \
            2>> $err_file >> $log_file
    '''
}

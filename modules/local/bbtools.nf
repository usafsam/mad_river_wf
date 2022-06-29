process BBMAP_ALIGN {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.maxcpus
    container params.container_bbtools
    memory {4.GB * task.attempt}
    errorStrategy {'retry'}
    maxRetries 2

    input:
        tuple val(sample), file(specimen)
        file(ref_genome)

    output:
        tuple val(sample), path("${task.process}/${sample}_mapped.sam.gz"), emit: aligned_sams
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        case $(echo "!{task.memory}" | cut -d' ' -f2) in
            [gG]B*) UNIT=1073741824;;
            [mM]B*) UNIT=1048576;;
            [kK]B*) UNIT=1024;;
            B*) UNIT=1;;
        esac
        MEMORY=$(( $(echo "!{task.memory}" | cut -d '.' -f1 | cut -d ' ' -f1) * $UNIT ))
        if [[ -z $MEMORY ]]; then
            XMX_FLAG="-Xmx$MEMORY"
        fi

        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools bbmap.sh: $(bbmap.sh -h | grep 'Last modified')" >> $log_file

        bbmap.sh \
            in=!{specimen} \
            outm=!{task.process}/!{sample}_mapped.sam.gz \
            ref=!{ref_genome} \
            nodisk \
            local \
            maxindel=500 \
            -Xms$MEMORY \
            threads=!{task.cpus} \
            ow \
            k=12 \
            2>> $err_file >> $log_file
    '''
}

process BBMERGE {
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_merge_ihist.txt"
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 1.GB

    input:
        tuple val(sample), file(reads)

    output:
        tuple val(sample), path("${task.process}/${sample}_merged.fq.gz"), emit: bbmerged_specimens
        path("${task.process}/${sample}_merge_ihist.txt")
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools bbmerge.sh: $(bbmerge.sh -h | grep 'Last modified')" >> $log_file

        bbmerge.sh \
            in1=!{reads[0]} \
            in2=!{reads[1]} \
            out=!{task.process}/!{sample}_merged.fq.gz \
            adapter=default \
            strict \
            ihist=!{task.process}/!{sample}_merge_ihist.txt \
            2>> $err_file >> $log_file
    '''
}

process PILEUP {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_basecov_border5.txt"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_covstats.txt"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_hist.txt"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 2.GB

    input:
        tuple val(sample), file(samfile)

    output:
        path("${task.process}/${sample}_basecov_border5.txt"), emit: basecov
        path("${task.process}/${sample}_covstats.txt"), emit: covstats
        path("${task.process}/${sample}_hist.txt"), emit: hist
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        case $(echo "!{task.memory}" | cut -d' ' -f2) in
            [gG]B*) UNIT=1073741824;;
            [mM]B*) UNIT=1048576;;
            [kK]B*) UNIT=1024;;
            B*) UNIT=1;;
        esac
        MEMORY=$(( $(echo "!{task.memory}" | cut -d '.' -f1 | cut -d ' ' -f1) * $UNIT ))
        if [[ -z $MEMORY ]]; then
            XMX_FLAG="-Xmx$MEMORY"
        fi

        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools pileup.sh: $(pileup.sh -h | grep 'Last modified')" >> $log_file

        pileup.sh \
            in=!{samfile} \
            basecov=!{task.process}/!{sample}_basecov_border5.txt \
            ow \
            border=5 \
            out=!{task.process}/!{sample}_covstats.txt \
            hist=!{task.process}/!{sample}_hist.txt \
            -Xms$MEMORY \
            threads=!{task.cpus} \
            2>> $err_file >> $log_file
    '''
}

process REMOVE_JUNK_DELS {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 8.GB
    errorStrategy 'retry'
    maxRetries 2

    input:
        tuple val(sample), file(samfile)
        file(ref_genome)

    output:
        tuple val(sample), path("${task.process}/${sample}_filtered.sam.gz"), optional: true, emit: filtered_junk_sams
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        case $(echo "!{task.memory}" | cut -d' ' -f2) in
            [gG]B*) UNIT=1073741824;;
            [mM]B*) UNIT=1048576;;
            [kK]B*) UNIT=1024;;
            B*) UNIT=1;;
        esac
        MEMORY=$(( $(echo "!{task.memory}" | cut -d '.' -f1 | cut -d ' ' -f1) * $UNIT ))
        if [[ -z $MEMORY ]]; then
            XMX_FLAG="-Xmx$MEMORY"
        fi

        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools filtersam.sh: $(filtersam.sh -h | grep 'Last modified')" >> $log_file

        if [[ $(stat -L -c "%s" !{samfile}) -lt 500 ]]; then
            echo "!{samfile} has a size of $(stat -L -c "%s" !{samfile}) bytes, which is less than the allowable threshold. Bailing out." >> $log_file;
            exit 0;
        fi
        filtersam.sh \
            ref=!{ref_genome} \
            ow \
            in=!{samfile} \
            out=!{task.process}/!{sample}_filtered.sam.gz \
            mbad=2 \
            del \
            sub=f \
            border=3 \
            mbv=0 \
            -Xms$MEMORY \
            threads=!{task.cpus} \
            2>> $err_file >> $log_file
    '''
}

process REMOVE_SINGLETONS {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 8.GB
    errorStrategy 'retry'
    maxRetries 2

    input:
        tuple val(sample), file(samfile)
        file(ref_genome)

    output:
        tuple val(sample), path("${task.process}/${sample}_filtered2.sam.gz"), optional: true, emit: filtered_singletons_sams
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        case $(echo "!{task.memory}" | cut -d' ' -f2) in
            [gG]B*) UNIT=1073741824;;
            [mM]B*) UNIT=1048576;;
            [kK]B*) UNIT=1024;;
            B*) UNIT=1;;
        esac
        MEMORY=$(( $(echo "!{task.memory}" | cut -d '.' -f1 | cut -d ' ' -f1) * $UNIT ))
        if [[ -z $MEMORY ]]; then
            XMX_FLAG="-Xmx$MEMORY"
        fi

        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools filtersam.sh: $(filtersam.sh -h | grep 'Last modified')" >> $log_file

        if [[ $(stat -L -c "%s" !{samfile}) -lt 500 ]]; then
            echo "!{samfile} has a size of $(stat -L -c "%s" !{samfile}) bytes, which is less than the allowable threshold. Bailing out." >> $log_file;
            exit 0;
        fi
        filtersam.sh \
            ref=!{ref_genome} \
            ow \
            in=!{samfile} \
            out=!{task.process}/!{sample}_filtered2.sam.gz \
            mbad=1 \
            sub \
            mbv=2 \
            -Xms$MEMORY \
            threads=!{task.cpus} \
            2>> $err_file >> $log_file
    '''
}

process SOFT_TRIM {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 2.GB

    input:
        tuple val(sample), file(samfile)

    output:
        tuple val(sample), path("${task.process}/${sample}_trimclip.sam.gz"), emit: soft_trimmed_sams
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        case $(echo "!{task.memory}" | cut -d' ' -f2) in
            [gG]B*) UNIT=1073741824;;
            [mM]B*) UNIT=1048576;;
            [kK]B*) UNIT=1024;;
            B*) UNIT=1;;
        esac
        MEMORY=$(( $(echo "!{task.memory}" | cut -d '.' -f1 | cut -d ' ' -f1) * $UNIT ))
        if [[ -z $MEMORY ]]; then
            XMX_FLAG="-Xmx$MEMORY"
        fi

        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools bbduk.sh: $(bbduk.sh -h | grep 'Last modified')" >> $log_file

        bbduk.sh \
            in=!{samfile} \
            out=!{task.process}/!{sample}_trimclip.sam.gz \
            -Xms$MEMORY \
            threads=!{task.cpus} \
            ow \
            2>> $err_file >> $log_file
    '''
}

process SUMMARIZE_COVERAGE {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "coverageSummary.txt"
    echo false
    cpus params.medcpus
    container params.container_bbtools

    input:
        file(basecov_files)

    output:
        path("coverageSummary.txt"), emit: coverage_summary
        path("logs/${task.process}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/
        log_file=logs/!{task.process}.log
        err_file=logs/!{task.process}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools summarizecoverage.sh: $(summarizecoverage.sh -h | grep 'Last modified')" >> $log_file

        summarizecoverage.sh !{basecov_files} out=coverageSummary.txt 2>> $err_file >> $log_file
    '''
}

process TAKE_VIRAL {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 2.GB

    input:
        tuple val(sample), file(specimen)
        file(ref_genome)

    output:
        tuple val(sample), path("${task.process}/${sample}_viral.fq.gz"), emit: viral_specimens
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        case $(echo "!{task.memory}" | cut -d' ' -f2) in
            [gG]B*) UNIT=1073741824;;
            [mM]B*) UNIT=1048576;;
            [kK]B*) UNIT=1024;;
            B*) UNIT=1;;
        esac
        MEMORY=$(( $(echo "!{task.memory}" | cut -d '.' -f1 | cut -d ' ' -f1) * $UNIT ))
        if [[ -z $MEMORY ]]; then
            XMX_FLAG="-Xmx$MEMORY"
        fi

        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools bbduk.sh: $(bbduk.sh -h | grep 'Last modified')" >> $log_file

        bbduk.sh \
            ow \
            -Xms$MEMORY \
            threads=!{task.cpus} \
            ref=!{ref_genome} \
            in=!{specimen} \
            outm=!{task.process}/!{sample}_viral.fq.gz \
            k=25 \
            2>> $err_file >> $log_file
    '''
}

process TRIM_ADAPTERS {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 2.GB

    input:
        tuple val(sample), file(specimen), file(adapters)

    output:
        tuple val(sample), path("${task.process}/${sample}_trimmed.fq.gz"), emit: adapter_trimmed_specimens
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        case $(echo "!{task.memory}" | cut -d' ' -f2) in
            [gG]B*) UNIT=1073741824;;
            [mM]B*) UNIT=1048576;;
            [kK]B*) UNIT=1024;;
            B*) UNIT=1;;
        esac
        MEMORY=$(( $(echo "!{task.memory}" | cut -d '.' -f1 | cut -d ' ' -f1) * $UNIT ))
        if [[ -z $MEMORY ]]; then
            XMX_FLAG="-Xmx$MEMORY"
        fi

        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools bbduk.sh: $(bbduk.sh -h | grep 'Last modified')" >> $log_file

        bbduk.sh \
            in=!{specimen} \
            out=!{task.process}/!{sample}_trimmed.fq.gz \
            minlen=60 \
            ktrim=r \
            k=19 \
            mink=9 \
            hdist=2 \
            hdist2=1 \
            ref=!{adapters} \
            maq=14 \
            qtrim=r \
            trimq=10 \
            maxns=0 \
            tbo \
            tpe \
            -Xms$MEMORY \
            threads=!{task.cpus} \
            ftm=5 \
            ow \
            2>> $err_file >> $log_file
    '''
}

process TRIM_PRIMERS {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 2.GB

    input:
        tuple val(sample), file(specimen), file(primers)

    output:
        tuple val(sample), path("${task.process}/${sample}_trimmed2.fq.gz"), emit: primer_trimmed_specimens
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        case $(echo "!{task.memory}" | cut -d' ' -f2) in
            [gG]B*) UNIT=1073741824;;
            [mM]B*) UNIT=1048576;;
            [kK]B*) UNIT=1024;;
            B*) UNIT=1;;
        esac
        MEMORY=$(( $(echo "!{task.memory}" | cut -d '.' -f1 | cut -d ' ' -f1) * $UNIT ))
        if [[ -z $MEMORY ]]; then
            XMX_FLAG="-Xmx$MEMORY"
        fi

        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools bbduk.sh: $(bbduk.sh -h | grep 'Last modified')" >> $log_file

        bbduk.sh \
            in=!{specimen} \
            out=!{task.process}/!{sample}_trimmed2.fq.gz \
            ref=!{primers} \
            ktrim=l \
            restrictleft=30 \
            hdist=3 \
            qhdist=1 \
            rcomp=f \
            mm=f \
            -Xms$MEMORY \
            threads=!{task.cpus} \
            2>> $err_file >> $log_file
    '''
}

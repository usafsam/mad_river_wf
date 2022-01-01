#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.skip_performance_lineage_excel = false

params.reads = workflow.launchDir + '/reads'
params.outdir = workflow.launchDir + '/results'
if (!params.skip_performance_lineage_excel) {
    params.run_info = workflow.launchDir + '/RunInfo.xml'
    params.sample_sheet = workflow.launchDir + '/SampleSheet.csv'
    params.stats_json = workflow.launchDir + '/Stats/Stats.json'
}

params.ivar_gff_gzip_file = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz"
params.primer_fasta = workflow.projectDir + "/reference/primers.fasta"
params.adapter_fasta = workflow.projectDir + "/reference/adapters.fa.gz"
params.nextclade_primers = workflow.projectDir + "/reference/midnight_primers_nextclade.csv"

params.container_biocontainers = 'biocontainers/biocontainers:v1.2.0_cv1'
params.container_bbtools = 'staphb/bbtools:latest'
params.container_ivar = 'staphb/ivar:latest'
params.container_nextclade = 'nextstrain/nextclade:latest'
params.container_pangolin = 'staphb/pangolin:latest'
params.container_samtools = 'staphb/samtools:latest'

params.mincov = 15
if (!params.skip_performance_lineage_excel) {
    params.protocol = '1200Midnight_XT'
}

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUs used in this workflow is ${params.maxcpus}")
if (params.maxcpus < 8) {
    params.medcpus = params.maxcpus
} else {
    params.medcpus = 8
}

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

process GRAB_IVAR_GFF {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    echo false
    container params.container_biocontainers

    input:
        file(gff_file_gzip)

    output:
        path("*.gff"), emit: gff_file_ivar

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
    '''
}

process INTERLEAVE {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 300.MB

    input:
        tuple val(sample), file(reads)

    output:
        tuple val(sample), path("${task.process}/${sample}_interleaved.fq.gz"), optional: true, emit: interleaved_specimens
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        echo "BBTools reformat.sh: $(reformat.sh -h | grep 'Last modified')" >> $log_file

        reformat.sh \
            in1=!{reads[0]} \
            in2=!{reads[1]} \
            out=!{task.process}/!{sample}_interleaved.fq.gz \
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
        tuple val(sample), file(specimen)

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
            in=!{specimen} \
            out=!{task.process}/!{sample}_merged.fq.gz \
            adapter=default \
            strict \
            ihist=!{task.process}/!{sample}_merge_ihist.txt \
            2>> $err_file >> $log_file
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

process SAMTOOLS_SORT_INDEX {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_trimclip.sorted.bam.gz"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_trimclip.sorted.bam.gz.bai"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_samtools

    input:
        tuple val(sample), file(samfile)

    output:
        tuple val(sample), path("${task.process}/${sample}_trimclip.sorted.bam.gz"), emit: bamfile
        path("${task.process}/${sample}_trimclip.sorted.bam.gz.bai"), emit: bamfile_index
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
                -o !{task.process}/!{sample}_trimclip.sorted.bam.gz \
                -@ !{task.cpus} \
                2>> $err_file >> $log_file
        samtools index \
            !{task.process}/!{sample}_trimclip.sorted.bam.gz \
            !{task.process}/!{sample}_trimclip.sorted.bam.gz.bai \
            2>> $err_file >> $log_file
    '''
}

process IVAR_CONSENSUS {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*.consensus.fa"
    tag "${sample}"
    echo false
    container params.container_ivar
    memory {2.GB * task.attempt}
    errorStrategy {'retry'}
    maxRetries 2

    input:
        tuple val(sample), file(bamfile)

    output:
        tuple val(sample), path("${task.process}/${sample}.consensus.fa"), emit: consensus
        path("${task.process}/${sample}.consensus.fa"), emit: consensus_collect
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        samtools --version >> $log_file
        ivar --version >> $log_file

        samtools mpileup -aa -A -Q 0 !{bamfile} 2>> $err_file | \
        ivar consensus -p !{task.process}/!{sample}.consensus \
            -t 0.6 \
            -m !{params.mincov} \
            2>> $err_file >> $log_file
    '''
}

process IVAR_VARIANTS {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*.variants.tsv"
    tag "${sample}"
    echo false
    container params.container_ivar
    memory {2.GB * task.attempt}
    errorStrategy {'retry'}
    maxRetries 2

    input:
        tuple val(sample), file(bamfile), file(ref_genome), file(ref_gff)

    output:
        path("${task.process}/${sample}.variants.tsv"), emit: ivar_tsvs
        path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        samtools --version >> $log_file

        samtools mpileup -aa -A -B -Q 0 --reference !{ref_genome} !{bamfile} 2>> $err_file | \
            ivar variants -t 0.05 \
            -m !{params.mincov} \
            -r !{ref_genome} \
            -g !{ref_gff} \
            -p !{task.process}/!{sample}.variants \
            2>> $err_file >> $log_file
    '''
}

params.traffic_light_plot_file_out = "pdf"
process TRAFFIC_LIGHT_PLOT {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.summary.${params.traffic_light_plot_file_out}"
    echo false
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

params.spike_gene_coverage_sample_column = 'seqName'
params.spike_gene_coverage_file_out = "xlsx"
process SPIKE_GENE_COVERAGE {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_21563-25384.${params.spike_gene_coverage_file_out}"
    echo false
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

process PERFORMANCE_LINEAGE_EXCEL {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_{lineage,performance}.xlsx"
    echo false
    conda "${workflow.projectDir}/env/performance_lineage_excel.yml"

    input:
        file(stats_json)
        file(coverage_summary)
        file(pangolin_csv)
        file(nextclade_csv)
        file(run_info)
        file(sample_sheet)

    output:
        path("${file(params.outdir).getSimpleName()}_{lineage,performance}.xlsx")
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

        performance_lineage_excel.py $(basename !{params.outdir}) \
            lineage \
            --sample_sheet !{sample_sheet} \
            --pangolin_summary !{pangolin_csv} \
            --nextclade_summary !{nextclade_csv} \
            2>> $err_file >> $log_file
    '''
}

workflow {
    println("Currently using the Mad River workflow.")

    // This is where the results will be
    println("The files and directory for the results is " + params.outdir)

    Channel
        .fromPath(params.ivar_gff_gzip_file, type:'file')
        .filter { fh ->
            fh.exists()
        }
        .ifEmpty {
            exit 1, "No GFF file was selected!\nDid you forget to set 'params.ivar_gff_gzip_file'?"
        }
        .view { "iVar GFF file for Reference Genome: $it" }
        .set { ivar_gff_gzip_file_checked }

    Channel
        .fromPath(params.primer_fasta, type:'file')
        .filter { fh ->
            fh.exists()
        }
        .ifEmpty {
            exit 1, "A FASTA file for primers is required!\nDid you forget to set 'params.primer_fasta'?"
        }
        .view { "Primer FASTA: $it" }
        .set { primer_fasta }

    Channel
        .fromPath(params.adapter_fasta, type:'file')
        .filter { fh ->
            fh.exists()
        }
        .ifEmpty {
            exit 1, "A FASTA.gz file for adapters is required!\nDid you forget to set 'params.adapter_fasta'?"
        }
        .view { "Adapter FASTA.gz: $it" }
        .set { adapter_fasta }

    Channel
        .fromPath(params.nextclade_primers, type:'file')
        .filter { fh ->
            fh.exists()
        }
        .ifEmpty {
            exit 1, "Nextclade requires a CSV file for primers!\nDid you forget to set 'params.nextclade_primers'?"
        }
        .view { "Nextclade primers CSV: $it"} 
        .set { nextclade_primers }

    Channel
        .fromFilePairs(["${params.reads}/*_R{1,2}*.fastq.gz", "${params.reads}/*_{1,2}.fastq*"], size: 2)
        .map { reads ->
            tuple(reads[0].replaceAll(~/_S[0-9]+(_L[0-9]+)?/,""), reads[1])
        }
        .filter { reads ->
            reads[1][0].exists() && reads[1][1].exists()
        }
        .ifEmpty {
            exit 1, "No fastq or fastq.gz files were found at ${params.reads}!\nDid you forget to set 'params.reads'?"
        }
        .set { paired_reads }

    if (!params.skip_performance_lineage_excel) {
        Channel
            .fromPath(params.stats_json, type:'file')
            .filter { fh ->
                fh.exists()
            }
            .ifEmpty {
                exit 1, "The performance excel file needs data from a Stats.json file!\nIt should have been generated when you invoked bcl2fastq.\nDid you forget to set 'params.stats_json'?"
            }
            .view { "bcl2fastq Stats.json file: $it" }
            .set { stats_json }

        Channel
            .fromPath(params.run_info, type:'file')
            .filter { fh ->
                fh.exists()
            }
            .ifEmpty {
                exit 1, "The performance excel file needs data from a RunInfo.xml file!\nIt should have been generated on the Illumina machine.\nDid you forget to set 'params.run_info'?"
            }
            .view { "Illumina RunInfo.xml file: $it" }
            .set { run_info_xml }

        Channel
            .fromPath(params.sample_sheet, type:'file')
            .filter { fh ->
                fh.exists()
            }
            .ifEmpty {
                exit 1, "The lineage excel file needs data from a SampleSheet.csv file!\nIt should have been generated for the Illumina machine.\nDid you forget to set 'params.sample_sheet'?"
            }
            .view { "SampleSheet.csv file: $it" }
            .set { sample_sheet_csv }
    }

    println("")

    // GRAB DATA
    // =========
    GRAB_NEXTCLADE_DATA() // out: reference_genome, reference_nextclade
    GRAB_IVAR_GFF(ivar_gff_gzip_file_checked) // out: gff_file_ivar

    // PRE-ALIGNMENT
    // =============
    INTERLEAVE(paired_reads) // out: interleaved_specimens, _
    BBMERGE(INTERLEAVE.out.interleaved_specimens) // out: bbmerged_specimens, _
    TAKE_VIRAL(BBMERGE.out.bbmerged_specimens, GRAB_NEXTCLADE_DATA.out.reference_genome) // out: viral_specimens, _
    ch_trim_adapters = TAKE_VIRAL.out.viral_specimens | combine(adapter_fasta)
    TRIM_ADAPTERS(ch_trim_adapters) // out: adapter_trimmed_specimens, _
    ch_trim_primers = TRIM_ADAPTERS.out.adapter_trimmed_specimens | combine(primer_fasta)
    TRIM_PRIMERS(ch_trim_primers) // out: primer_trimmed_specimens, _

    BBMAP_ALIGN(TRIM_PRIMERS.out.primer_trimmed_specimens, GRAB_NEXTCLADE_DATA.out.reference_genome) // out: aligned_sams, _

    // POST-ALIGNMENT
    // ==============
    REMOVE_JUNK_DELS(BBMAP_ALIGN.out.aligned_sams, GRAB_NEXTCLADE_DATA.out.reference_genome) // out: filtered_junk_sams, _
    REMOVE_SINGLETONS(REMOVE_JUNK_DELS.out.filtered_junk_sams, GRAB_NEXTCLADE_DATA.out.reference_genome) // out: filtered_singletons_sams, _
    SOFT_TRIM(REMOVE_SINGLETONS.out.filtered_singletons_sams) // out: soft_trimmed_sams, _
    PILEUP(SOFT_TRIM.out.soft_trimmed_sams) // out: basecov, covstats, hist, _
    SUMMARIZE_COVERAGE(PILEUP.out.basecov.toSortedList()) // out: coverage_summary, _

    // SAMTOOLS AND IVAR
    // =================
    SAMTOOLS_SORT_INDEX(SOFT_TRIM.out.soft_trimmed_sams) // out: bamfile, bamfile_index, _
    IVAR_CONSENSUS(SAMTOOLS_SORT_INDEX.out.bamfile) // out: consensus, consensus_collect, _
    ch_ivar_variants = SAMTOOLS_SORT_INDEX.out.bamfile |
        combine(GRAB_NEXTCLADE_DATA.out.reference_genome) |
        combine(GRAB_IVAR_GFF.out.gff_file_ivar)
    IVAR_VARIANTS(ch_ivar_variants) // out: ivar_tsvs, _
    TRAFFIC_LIGHT_PLOT(IVAR_VARIANTS.out.ivar_tsvs.toSortedList(), SUMMARIZE_COVERAGE.out.coverage_summary) // out: _, _

    // PANGOLIN AND NEXTCLADE
    // ======================
    PANGOLIN(IVAR_CONSENSUS.out.consensus) // out: pangolin, _
    ch_combined_pangolin_report = PANGOLIN.out.pangolin |
        collectFile(
            name: "combined_lineage_report.csv",
            keepHeader: true,
            newLine: true,
            sort: true,
            storeDir: "${params.outdir}"
        )
    ch_nextclade = IVAR_CONSENSUS.out.consensus |
        combine(GRAB_NEXTCLADE_DATA.out.reference_nextclade) |
        combine(nextclade_primers)
    NEXTCLADE(ch_nextclade) // out: nextclade, _
    ch_combined_nextclade_report = NEXTCLADE.out.nextclade |
        collectFile(
            name: "combined_nextclade_report.csv",
            keepHeader: true,
            newLine: true,
            sort: true,
            storeDir: "${params.outdir}"
        )
    SPIKE_GENE_COVERAGE(ch_combined_nextclade_report)
    if (!params.skip_performance_lineage_excel) {
        PERFORMANCE_LINEAGE_EXCEL(
            stats_json,
            SUMMARIZE_COVERAGE.out.coverage_summary,
            ch_combined_pangolin_report,
            ch_combined_nextclade_report,
            run_info_xml,
            sample_sheet_csv
        )
    }
}

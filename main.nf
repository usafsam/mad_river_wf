#!/usr/bin/env nextflow

println("Currently using the Mad River workflow.")
println("")

params.reads = workflow.launchDir + '/reads'
params.outdir = workflow.launchDir + '/results'
params.run_info = workflow.launchDir + '/RunInfo.xml'
params.sample_sheet = workflow.launchDir + '/SampleSheet.csv'
params.stats_json = workflow.launchDir + '/Stats/Stats.json'

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
params.protocol = '1200Midnight_XT'

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUs used in this workflow is ${params.maxcpus}")
if (params.maxcpus < 5) {
    params.medcpus = params.maxcpus
} else {
    params.medcpus = 8
}

// This is where the results will be
println("The files and directory for the results is " + params.outdir)

Channel
    .fromPath(params.ivar_gff_gzip_file, type:'file')
    .ifEmpty{
        println("No GFF file was selected!")
        println("Did you forget to set 'params.ivar_gff_gzip_file'?")
        exit 1
    }
    .view { "iVar GFF file for Reference Genome : $it" }
    .set {ivar_gff_gzip_file_checked}

Channel
    .fromPath(params.primer_fasta, type:'file')
    .filter{ fh ->
        fh.exists()
    }
    .ifEmpty{
        println("A FASTA file for primers is required!")
        println("Did you forget to set 'params.primer_fasta'?")
        exit 1
    }
    .view { "Primer FASTA : $it" }
    .set { primer_fasta }

Channel
    .fromPath(params.adapter_fasta, type:'file')
    .filter{ fh ->
        fh.exists()
    }
    .ifEmpty{
        println("A FASTA.gz file for adapters is required!")
        println("Did you forget to set 'params.adapter_fasta'?")
        exit 1
    }
    .view { "Adapter FASTA.gz : $it"}
    .set { adapter_fasta }

Channel
    .fromPath(params.nextclade_primers, type:'file')
    .filter{ fh ->
        fh.exists()
    }
    .ifEmpty{
        println("Nextclade requires a CSV file for primers!")
        println("Did you forget to set 'params.nextclade_primers'?")
        exit 1
    }
    .view { "Nextclade primers CSV: $it"}
    .set { nextclade_primers }

Channel
    .fromFilePairs(["${params.reads}/*_R{1,2}*.fastq.gz", "${params.reads}/*_{1,2}.fastq*"], size: 2)
    .map{reads -> tuple(reads[0].replaceAll(~/_S[0-9]+(_L[0-9]+)?/,""), reads[1])}
    .filter{ reads ->
        reads[1][0].exists() && reads[1][1].exists()
    }
    .ifEmpty{
        println("No fastq or fastq.gz files were found at ${params.reads}!")
        println("Did you forget to set 'params.reads'?")
        exit 1
    }
    .set {paired_reads}


Channel
    .fromPath(params.stats_json, type:'file')
    .filter{ fh ->
        fh.exists()
    }
    .ifEmpty {
        println("The performance excel file needs data from a Stats.json file!")
        println("It should have been generated when you invoked bcl2fastq.")
        println("Did you forget to set 'params.stats_json'?")
        exit 1
    }
    .view {"bcl2fastq Stats.json file: $it"}
    .set {stats_json}

Channel
    .fromPath(params.run_info, type:'file')
    .filter{ fh ->
        fh.exists()
    }
    .ifEmpty {
        println("The performance excel file needs data from a RunInfo.xml file!")
        println("It should have been generated on the Illumina machine.")
        println("Did you forget to set 'params.run_info'?")
        exit 1
    }
    .view {"Illumina RunInfo.xml file: $it"}
    .set {run_info_xml}

Channel
    .fromPath(params.sample_sheet, type:'file')
    .filter{ fh ->
        fh.exists()
    }
    .ifEmpty {
        println("The lineage excel file needs data from a SampleSheet.csv file!")
        println("It should have been generated for the Illumina machine.")
        println("Did you forget to set 'params.sample_sheet'?")
        exit 1
    }
    .view {"SampleSheet.csv file: $it"}
    .set {sample_sheet_csv}

println("")

process grab_nextclade_data {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    echo false
    container params.container_nextclade

    output:
    file("data/sars-cov-2_MN908947/reference.fasta") into reference_genome
    file("data/sars-cov-2_MN908947/") into reference_nextclade

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

reference_genome
    .into { reference_genome_take_viral; reference_genome_bbmap_align; reference_genome_remove_junk_dels; reference_genome_remove_singletons; reference_genome_ivar_variants }

process grab_ivar_gff {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    echo false
    container params.container_biocontainers

    input:
    file(gff_file_gzip) from ivar_gff_gzip_file_checked

    output:
    file("*.gff") into gff_file_ivar

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

process interleave {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 300.MB

    input:
    set val(sample), file(reads) from paired_reads

    output:
    tuple sample, file("${task.process}/${sample}_interleaved.fq.gz") optional true into interleaved_specimens
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

process bbmerge {
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_merge_ihist.txt"
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 1.GB

    input:
    set val(sample), file(specimen) from interleaved_specimens

    output:
    tuple sample, file("${task.process}/${sample}_merged.fq.gz") into bbmerged_specimens
    file("${task.process}/${sample}_merge_ihist.txt")
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

bbmerged_specimens
    .combine(reference_genome_take_viral)
    .set { take_viral_channel }

process take_viral {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 2.GB

    input:
    set val(sample), file(specimen), file(ref_genome) from take_viral_channel

    output:
    tuple sample, file("${task.process}/${sample}_viral.fq.gz") into viral_specimens
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

viral_specimens
    .combine(adapter_fasta)
    .set { trim_adapters_channel }

process trim_adapters {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 2.GB

    input:
    set val(sample), file(specimen), file(adapters) from trim_adapters_channel

    output:
    tuple sample, file("${task.process}/${sample}_trimmed.fq.gz") into adapter_trimmed_specimens
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

adapter_trimmed_specimens
    .combine(primer_fasta)
    .set { trim_primers_channel }

process trim_primers {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 2.GB

    input:
    set val(sample), file(specimen), file(primers) from trim_primers_channel

    output:
    tuple sample, file("${task.process}/${sample}_trimmed2.fq.gz") into primer_trimmed_specimens
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

primer_trimmed_specimens
    .combine(reference_genome_bbmap_align)
    .set { bbmap_align_channel }

process bbmap_align {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.maxcpus
    container params.container_bbtools
    memory {4.GB * task.attempt}
    errorStrategy {'retry'}
    maxRetries 2

    input:
    set val(sample), file(specimen), file(ref_genome) from bbmap_align_channel

    output:
    tuple sample, file("${task.process}/${sample}_mapped.sam.gz") into aligned_sams
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

aligned_sams
    .combine(reference_genome_remove_junk_dels)
    .set { remove_junk_dels_channel }

process remove_junk_dels {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 8.GB
    errorStrategy 'retry'
    maxRetries 2

    input:
    set val(sample), file(samfile), file(ref_genome) from remove_junk_dels_channel

    output:
    tuple sample, file("${task.process}/${sample}_filtered.sam.gz") optional true into filtered_junk_sams
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

filtered_junk_sams
    .combine(reference_genome_remove_singletons)
    .set { remove_singletons_channel }

process remove_singletons {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 8.GB
    errorStrategy 'retry'
    maxRetries 2

    input:
    set val(sample), file(samfile), file(ref_genome) from remove_singletons_channel

    output:
    tuple sample, file("${task.process}/${sample}_filtered2.sam.gz") optional true into filtered_singletons_sams
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

process soft_trim {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_bbtools
    memory 2.GB

    input:
    set val(sample), file(samfile) from filtered_singletons_sams

    output:
    tuple sample, file("${task.process}/${sample}_trimclip.sam.gz") into soft_trimmed_sams
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

soft_trimmed_sams
    .into { soft_trimmed_sams_pileup; soft_trimmed_sams_samtools_sort_index }

process pileup {
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
    set val(sample), file(samfile) from soft_trimmed_sams_pileup

    output:
    file("${task.process}/${sample}_basecov_border5.txt") into basecov
    file("${task.process}/${sample}_covstats.txt")
    file("${task.process}/${sample}_hist.txt")
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

process summarize_coverage {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "coverageSummary.txt"
    echo false
    cpus params.medcpus
    container params.container_bbtools

    input:
    file(basecov_files) from basecov.toSortedList()

    output:
    file("coverageSummary.txt") into coverageSummary
    file("logs/${task.process}.{log,err}")

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

coverageSummary
    .into {coverageSummary_traffic_light_plot; coverageSummary_performance_lineage_excel}

process samtools_sort_index {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_trimclip.sorted.bam.gz"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*_trimclip.sorted.bam.gz.bai"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_samtools

    input:
    set val(sample), file(samfile) from soft_trimmed_sams_samtools_sort_index

    output:
    tuple sample, file("${task.process}/${sample}_trimclip.sorted.bam.gz") into bamfile
    file("${task.process}/${sample}_trimclip.sorted.bam.gz")
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

bamfile
    .into { bamfile_consensus ; bamfile_variants }

process ivar_consensus {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*.consensus.fa"
    tag "${sample}"
    echo false
    container params.container_ivar
    memory {2.GB * task.attempt}
    errorStrategy {'retry'}
    maxRetries 2

    input:
    set val(sample), file(bamfile) from bamfile_consensus

    output:
    tuple val(sample), file("${task.process}/${sample}.consensus.fa") into consensus_nextclade, consensus_pangolin
    file("${task.process}/${sample}.consensus.fa") into consensus_collect
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

consensus_collect
    .collectFile(name: "${file(params.outdir).getSimpleName()}_consensus.fasta", sort: true, storeDir: "${params.outdir}")

bamfile_variants
    .combine(reference_genome_ivar_variants)
    .combine(gff_file_ivar)
    .set { ivar_variants_channel }

process ivar_variants {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*.variants.tsv"
    tag "${sample}"
    echo false
    container params.container_ivar
    memory {2.GB * task.attempt}
    errorStrategy {'retry'}
    maxRetries 2

    input:
    set val(sample), file(bamfile), file(ref_genome), file(ref_gff) from ivar_variants_channel

    output:
    file("${task.process}/${sample}.variants.tsv") into ivar_tsvs
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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
process traffic_light_plot {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.summary.${params.traffic_light_plot_file_out}"
    echo false
    conda "${workflow.projectDir}/env/traffic_light_plot.yml"

    input:
    file(ivar_tsv) from ivar_tsvs.toSortedList()
    file(coverage_summary) from coverageSummary_traffic_light_plot

    output:
    file("${file(params.outdir).getSimpleName()}.summary.${params.traffic_light_plot_file_out}")
    file("logs/${task.process}.{log,err}")

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
process pangolin {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*/lineage_report.csv"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_pangolin

    input:
    set val(sample), file(fasta) from consensus_pangolin

    output:
    file("${task.process}/${sample}/lineage_report.csv") into pangolin_files
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
        mkdir -p !{task.process} logs/!{task.process}
        log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
        err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

        # time stamp + capturing tool versions
        date | tee -a $log_file $err_file > /dev/null
        pangolin --version >> $log_file
        pangolin --pangoLEARN-version >> $log_file

        pangolin !{params.pangolin_options} \
            --outdir !{task.process}/!{sample} \
            !{fasta} \
            2>> $err_file >> $log_file
    '''
}

pangolin_files
    .collectFile(name: "combined_lineage_report.csv", keepHeader: true, newLine: true, sort: true, storeDir: "${params.outdir}")
    .set { combined_pangolin_report_performance_lineage_excel }

consensus_nextclade
    .combine(reference_nextclade)
    .combine(nextclade_primers)
    .set { nextclade_channel }

params.nextclade_options = ''
process nextclade {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/*"
    tag "${sample}"
    echo false
    cpus params.medcpus
    container params.container_nextclade

    input:
    set val(sample), file(fasta), file(ref_nextclade), file(primers_csv) from nextclade_channel

    output:
    file("${task.process}/${sample}_nextclade.csv") into nextclade_files
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

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

nextclade_files
    .collectFile(name: "combined_nextclade_report.csv", keepHeader: true, newLine: true, sort: true, storeDir: "${params.outdir}")
    .into{combined_nextclade_report_spike_gene_coverage; combined_nextclade_report_performance_lineage_excel}

params.spike_gene_coverage_sample_column = 'seqName'
params.spike_gene_coverage_file_out = "xlsx"
process spike_gene_coverage {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_21563-25384.${params.spike_gene_coverage_file_out}"
    echo false
    conda "${workflow.projectDir}/env/nextalign_gene_coverage.yml"

    input:
    file(nextclade_csv) from combined_nextclade_report_spike_gene_coverage

    output:
    file("${file(params.outdir).getSimpleName()}_nextclade_21563-25384.${params.spike_gene_coverage_file_out}")
    file("logs/${task.process}.{log,err}")

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

process performance_lineage_excel {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}.{log,err}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*_{lineage,performance}.xlsx"
    echo false
    conda "${workflow.projectDir}/env/performance_lineage_excel.yml"

    input:
    file(stats_json) from stats_json
    file(coverage_summary) from coverageSummary_performance_lineage_excel
    file(pangolin_csv) from combined_pangolin_report_performance_lineage_excel
    file(nextclade_csv) from combined_nextclade_report_performance_lineage_excel
    file(run_info) from run_info_xml
    file(sample_sheet) from sample_sheet_csv

    output:
    file("${file(params.outdir).getSimpleName()}_{lineage,performance}.xlsx")
    file("logs/${task.process}.{log,err}")

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

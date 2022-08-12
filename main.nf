#!/usr/bin/env nextflow

println("Currently using the Mad River workflow.")
println("Author: Padraic Fanning")
println("")

nextflow.enable.dsl=2

// # Params

// ## Input Params

params.skip_performance_excel = false
params.run_as_single_reads = false

params.reads = workflow.launchDir + '/reads'
params.outdir = workflow.launchDir + '/results'

println("The files and directory for the results is " + params.outdir)

if (!params.skip_performance_excel) {
    params.run_info = workflow.launchDir + '/RunInfo.xml'
    params.stats_json = workflow.launchDir + '/Stats/Stats.json'
}

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUs used in this workflow is ${params.maxcpus}")

// if (params.maxcpus < 2) {
//     params.halfmaxcpus = 1
// } else {
//     params.halfmaxcpus = Math.ceil(params.maxcpus / 2)
// }

if (params.maxcpus < 8) {
    params.medcpus = params.maxcpus
} else {
    params.medcpus = 8
}

if (params.maxcpus < 4) {
    params.halfmedcpus = params.maxcpus
} else {
    params.halfmedcpus = Math.floor(params.medcpus / 2)
}

// ## Species Params
params.species = 'sarscov2'
switch (params.species) {
    case 'sarscov2':
        println("Using the subworkflow and parameters for SARS-CoV-2")
        // Nextclade params for SARS-CoV-2
        params.do_nextclade = true
        params.nextclade_reference_name = 'sars-cov-2'
        // Primers
        params.primer_tsv = workflow.projectDir + "/reference/IDT_Midnight_Primers_v2.0.tsv"
        // iVar GFF
        params.ivar_gff_gzip_file = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz"
        // VADR params for SARS-CoV-2
        params.do_vadr = true
        params.vadr_reference = 'sarscov2'
        params.vadr_options = '--split --glsearch -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn'
        params.vadr_trim_options = '--minlen 50 --maxlen 30000'
        params.vadr_mdir = '/opt/vadr/vadr-models'
        // Pangolin params
        params.do_pangolin = true
        params.pangolin_options = '--max-ambig 0.5 --min-length 10000'
        // Other params
        params.do_spike_gene_coverage = true
        if (!params.skip_performance_excel) {
            params.protocol = '1200Midnight_XT'
        }
        break
    case 'hmpxv':
        println("Using the subworkflow and parameters for Monkeypox Virus")
        // Nextclade params for Monkeypox
        params.do_nextclade = true
        params.nextclade_reference_name = 'hMPXV'
        // Primers
        params.primer_tsv = workflow.projectDir + "/reference/hMPXV_YSPH_Primers_v1.0.tsv"
        // iVar GFF
        params.ivar_gff_gzip_file = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/621/545/GCF_014621545.1_ASM1462154v1/GCF_014621545.1_ASM1462154v1_genomic.gff.gz"
        // VADR params for Monkeypox
        params.do_vadr = true
        params.vadr_reference = 'mpxv'
        params.vadr_options = '--split --glsearch -s -r --nomisc --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin'
        params.vadr_trim_options = '--minlen 50 --maxlen 210000'
        params.vadr_mdir = '/opt/vadr/vadr-models'
        // Pangolin params
        params.do_pangolin = false
        params.pangolin_options = ''
        // Other params
        params.do_spike_gene_coverage = false
        if (!params.skip_performance_excel) {
            params.protocol = '2000YSPH_XT'
        }
        break
    default:
        // Nextclade params
        params.do_nextclade = false
        params.nextclade_reference_name = ''
        // Primers
        params.primer_tsv = ''
        // iVar GFF (directly from NCBI)
        params.ivar_gff_gzip_file = ''
        // VADR params
        params.do_vadr = false
        params.vadr_reference = ''
        params.vadr_options = ''
        params.vadr_trim_options = ''
        params.vadr_mdir = ''
        // Pangolin params
        params.do_pangolin = false
        params.pangolin_options = ''
        // Other params
        params.do_spike_gene_coverage = false
        break
}

include { analysis_prep; } from './subworkflows/analysis_prep.nf' addParams(
    outdir: params.outdir,
    container_nextclade: params.container_nextclade,
    container_python: params.container_python
)
include { analysis_core; } from './subworkflows/analysis_core.nf' addParams(
    outdir: params.outdir,
    maxcpus: params.maxcpus,
    // halfmaxcpus: params.halfmaxcpus,
    medcpus: params.medcpus,
    halfmedcpus: params.halfmedcpus,
    container_bbtools: params.container_bbtools,
    container_ivar: params.container_ivar,
    container_samtools: params.container_samtools,
    mincov: params.mincov
)
if (!params.skip_performance_excel) {
    include { performance_excel; } from './subworkflows/performance_excel.nf' addParams(
        outdir: params.outdir,
        container_python: params.container_python,
        protocol: params.protocol
    )
}
include { sarscov2; } from './subworkflows/sars-cov-2.nf' addParams(
    outdir: params.outdir,
    medcpus: params.medcpus,
    container_nextclade: params.container_nextclade,
    container_pangolin: params.container_pangolin,
    container_vadr: params.container_vadr,
    nextclade_options: params.nextclade_options,
    pangolin_options: params.pangolin_options,
    vadr_options: params.vadr_options,
    vadr_reference: params.vadr_reference,
    vadr_trim_options: params.vadr_trim_options,
    vadr_mdir: params.vadr_mdir
)
include { hmpxv; } from './subworkflows/hMPXV.nf' addParams(
    outdir: params.outdir,
    medcpus: params.medcpus,
    container_nextclade: params.container_nextclade,
    container_vadr: params.container_vadr,
    nextclade_options: params.nextclade_options,
    vadr_options: params.vadr_options,
    vadr_reference: params.vadr_reference,
    vadr_trim_options: params.vadr_trim_options,
    vadr_mdir: params.vadr_mdir
)

Channel
    .fromPath(params.ivar_gff_gzip_file, type:'file')
    .filter { fh ->
        fh.exists()
    }
    .ifEmpty {
        exit 1, """No GFF file was selected!
        Did you forget to set 'params.ivar_gff_gzip_file'?""".stripIndent()
    }
    .set { ivar_gff_gzip_file_checked }

Channel
    .fromPath(params.primer_tsv, type:'file')
    .filter { fh ->
        fh.exists()
    }
    .ifEmpty {
        exit 1, """A TSV file for primers is required!
        Did you forget to set 'params.primer_tsv'?""".stripIndent()
    }
    .set { primer_tsv }

Channel
    .fromPath(params.adapter_fasta, type:'file')
    .filter { fh ->
        fh.exists()
    }
    .ifEmpty {
        exit 1, """A FASTA.gz file for adapters is required!
        Did you forget to set 'params.adapter_fasta'?""".stripIndent()
    }
    .first()
    .set { adapter_fasta }

if (!params.run_as_single_reads) {
    Channel
        .fromFilePairs([
                "${params.reads}/*_R{1,2}*.fastq.gz",
                "${params.reads}/*_{1,2}.fastq*"
            ],
            size: 2
        )
        .map { reads ->
            tuple(reads[0].replaceAll(~/_S[0-9]+(_L[0-9]+)?/,""), reads[1])
        }
        .filter { reads ->
            reads[1][0].exists() && reads[1][1].exists()
        }
        .ifEmpty {
            exit 1, """No fastq or fastq.gz files were found at ${params.reads}!
            Did you forget to set 'params.reads'?""".stripIndent()
        }
        .set { paired_reads }
    Channel.empty().set { single_reads }
}
else {
    Channel.empty().set { paired_reads }
    Channel
        .fromPath("${params.reads}/*.fastq.gz")
        .map { reads -> tuple(reads.simpleName, reads) }
        .filter { reads -> reads[1].exists() }
        .ifEmpty {
            exit 1, """No fastq or fastq.gz files were found at ${params.reads}!
            Did you forget to set 'params.reads'?""".stripIndent()
        }
        .set { single_reads }
}

if (!params.skip_performance_excel) {
    Channel
        .fromPath(params.stats_json, type:'file')
        .filter { fh ->
            fh.exists()
        }
        .ifEmpty {
            exit 1, """The performance excel file needs data from a Stats.json file!
            It should have been generated when you invoked bcl2fastq.
            Did you forget to set 'params.stats_json'?""".stripIndent()
        }
        .first()
        .set { stats_json }

    Channel
        .fromPath(params.run_info, type:'file')
        .filter { fh ->
            fh.exists()
        }
        .ifEmpty {
            exit 1, """The performance excel file needs data from a RunInfo.xml file!
            It should have been generated on the Illumina machine.
            Did you forget to set 'params.run_info'?""".stripIndent()
        }
        .first()
        .set { run_info_xml }
}

workflow {
    ivar_gff_gzip_file_checked.view {
        "iVar GFF file for Reference Genome: $it"
    }
    primer_tsv.view { "Primer TSV: $it" }
    adapter_fasta.view { "Adapter FASTA.gz: $it" }
    if (!params.skip_performance_excel) {
        stats_json.view { "bcl2fastq Stats.json file: $it" }
        run_info_xml.view { "Illumina RunInfo.xml file: $it" }
    }

    analysis_prep(ivar_gff_gzip_file_checked, primer_tsv)
    analysis_core(
        paired_reads,
        single_reads,
        analysis_prep.out.reference_genome,
        adapter_fasta,
        analysis_prep.out.primer_fasta.first(),
        analysis_prep.out.gff_file_ivar.first()
    )
    if (!params.skip_performance_excel) {
        performance_excel(
            stats_json,
            analysis_core.out.coverage_summary,
            run_info_xml
        )
    }
    switch (params.species) {
        case 'sarscov2':
            sarscov2(
                analysis_core.out.combined_consensus_fasta,
                analysis_prep.out.reference_nextclade
            )
            break
        case 'hmpxv':
            hmpxv(
                analysis_core.out.combined_consensus_fasta,
                analysis_prep.out.reference_nextclade
            )
            break
        default:
            break
    }
}

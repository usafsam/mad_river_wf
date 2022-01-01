#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = workflow.launchDir + '/reads'
params.outdir = workflow.launchDir + '/results'
// params.run_info = workflow.launchDir + '/RunInfo.xml'
// params.sample_sheet = workflow.launchDir + '/SampleSheet.csv'
// params.stats_json = workflow.launchDir + '/Stats/Stats.json'

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

    // Channel
    //     .fromPath(params.stats_json, type:'file')
    //     .filter { fh ->
    //         fh.exists()
    //     }
    //     .ifEmpty {
    //         exit 1, "The performance excel file needs data from a Stats.json file!\nIt should have been generated when you invoked bcl2fastq.\nDid you forget to set 'params.stats_json'?"
    //     }
    //     .view { "bcl2fastq Stats.json file: $it" }
    //     .set { stats_json }

    // Channel
    //     .fromPath(params.run_info, type:'file')
    //     .filter { fh ->
    //         fh.exists()
    //     }
    //     .ifEmpty {
    //         exit 1, "The performance excel file needs data from a RunInfo.xml file!\nIt should have been generated on the Illumina machine.\nDid you forget to set 'params.run_info'?"
    //     }
    //     .view { "Illumina RunInfo.xml file: $it" }
    //     .set { run_info_xml }

    // Channel
    //     .fromPath(params.sample_sheet, type:'file')
    //     .filter { fh ->
    //         fh.exists()
    //     }
    //     .ifEmpty {
    //         exit 1, "The lineage excel file needs data from a SampleSheet.csv file!\nIt should have been generated for the Illumina machine.\nDid you forget to set 'params.sample_sheet'?"
    //     }
    //     .view { "SampleSheet.csv file: $it" }
    //     .set { sample_sheet_csv }

    println("")

}

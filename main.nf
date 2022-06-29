#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.skip_performance_excel = false
params.run_as_single_reads = false

params.reads = workflow.launchDir + '/reads'
params.outdir = workflow.launchDir + '/results'
if (!params.skip_performance_excel) {
    params.run_info = workflow.launchDir + '/RunInfo.xml'
    params.stats_json = workflow.launchDir + '/Stats/Stats.json'
}

params.ivar_gff_gzip_file = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz"
params.primer_tsv = workflow.projectDir + "/reference/IDT_Midnight_Primers_v2.0.tsv"
params.adapter_fasta = workflow.projectDir + "/reference/adapters.fa.gz"

params.container_bbtools = 'staphb/bbtools:latest'
params.container_ivar = 'staphb/ivar:latest'
params.container_nextclade = 'nextstrain/nextclade:latest'
params.container_pangolin = 'staphb/pangolin:latest'
params.container_python = 'python:latest'
params.container_samtools = 'staphb/samtools:latest'
params.container_vadr = 'staphb/vadr:latest'

params.mincov = 15
if (!params.skip_performance_excel) {
    params.protocol = '1200Midnight_XT'
}

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUs used in this workflow is ${params.maxcpus}")
if (params.maxcpus < 8) {
    params.medcpus = params.maxcpus
} else {
    params.medcpus = 8
}

include {
    BBMAP_ALIGN;
    BBMERGE;
    PILEUP;
    REMOVE_JUNK_DELS;
    REMOVE_SINGLETONS;
    SOFT_TRIM;
    SUMMARIZE_COVERAGE;
    TAKE_VIRAL;
    TRIM_ADAPTERS;
    TRIM_PRIMERS
} from './modules/local/bbtools'
include { IVAR_CONSENSUS; IVAR_VARIANTS } from './modules/local/ivar'
include {
    DATA_PREP;
    LINEAGE_EXCEL;
    PERFORMANCE_EXCEL;
    SPIKE_GENE_COVERAGE;
    TRAFFIC_LIGHT_PLOT
} from './modules/local/misc'
include { GRAB_NEXTCLADE_DATA; NEXTCLADE } from './modules/local/nextclade'
include { PANGOLIN } from './modules/local/pangolin'
include { SAMTOOLS_SORT_INDEX } from './modules/local/samtools'
include { VADR } from './modules/local/vadr'

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
        .fromPath(params.primer_tsv, type:'file')
        .filter { fh ->
            fh.exists()
        }
        .ifEmpty {
            exit 1, "A TSV file for primers is required!\nDid you forget to set 'params.primer_tsv'?"
        }
        .view { "Primer TSV: $it" }
        .set { primer_tsv }

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

    if (!params.run_as_single_reads) {
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
        Channel.empty().set { single_reads }
    }
    else {
        Channel.empty().set { paired_reads }
        Channel
            .fromPath("${params.reads}/*.fastq.gz")
            .map { reads -> tuple(reads.simpleName, reads) }
            .filter { reads -> reads[1].exists() }
            .ifEmpty {
                exit 1, "No fastq or fastq.gz files were found at ${params.reads}!\nDid you forget to set 'params.reads'?"
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
    }

    println("")

    // GRAB DATA
    // =========
    GRAB_NEXTCLADE_DATA() // out: reference_genome, reference_nextclade
    DATA_PREP(ivar_gff_gzip_file_checked, primer_tsv) // out: gff_file_ivar, primer_fasta

    // PRE-ALIGNMENT
    // =============
    BBMERGE(paired_reads) // out: bbmerged_specimens, _
    ch_take_viral = BBMERGE.out.bbmerged_specimens | concat(single_reads)
    TAKE_VIRAL(ch_take_viral, GRAB_NEXTCLADE_DATA.out.reference_genome) // out: viral_specimens, _
    ch_trim_adapters = TAKE_VIRAL.out.viral_specimens | combine(adapter_fasta)
    TRIM_ADAPTERS(ch_trim_adapters) // out: adapter_trimmed_specimens, _
    ch_trim_primers = TRIM_ADAPTERS.out.adapter_trimmed_specimens | combine(DATA_PREP.out.primer_fasta)
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
    ch_combined_consensus_fasta = IVAR_CONSENSUS.out.consensus_collect |
        collectFile(
            name: "${file(params.outdir).getSimpleName()}_consensus.fasta",
            sort: { it.getSimpleName() },
            storeDir: "${params.outdir}"
        )
    ch_ivar_variants = SAMTOOLS_SORT_INDEX.out.bamfile |
        combine(GRAB_NEXTCLADE_DATA.out.reference_genome) |
        combine(DATA_PREP.out.gff_file_ivar)
    IVAR_VARIANTS(ch_ivar_variants) // out: ivar_tsvs, _
    TRAFFIC_LIGHT_PLOT(IVAR_VARIANTS.out.ivar_tsvs.toSortedList(), SUMMARIZE_COVERAGE.out.coverage_summary) // out: _, _

    // PANGOLIN, NEXTCLADE, AND VADR
    // =============================
    PANGOLIN(ch_combined_consensus_fasta) // out: pangolin, _
    NEXTCLADE(ch_combined_consensus_fasta, GRAB_NEXTCLADE_DATA.out.reference_nextclade) // out: nextclade_csv, nextclade_json, nextclade_auspice_json, nextclade_out_files, _
    VADR(ch_combined_consensus_fasta) // out: vadr, vadr_file, _
    SPIKE_GENE_COVERAGE(NEXTCLADE.out.nextclade_csv)
    LINEAGE_EXCEL(
        PANGOLIN.out.pangolin,
        NEXTCLADE.out.nextclade_csv,
        VADR.out.vadr_file
    )
    if (!params.skip_performance_excel) {
        PERFORMANCE_EXCEL(
            stats_json,
            SUMMARIZE_COVERAGE.out.coverage_summary,
            run_info_xml
        )
    }
}

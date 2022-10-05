include {
    BBMAP_ALIGN;
    BBMERGE;
    INTERLEAVE;
    PILEUP;
    REMOVE_JUNK_DELS;
    REMOVE_SINGLETONS;
    SOFT_TRIM;
    SUMMARIZE_COVERAGE;
    TAKE_VIRAL;
    TRIM_ADAPTERS;
    TRIM_PRIMERS
} from '../modules/local/bbtools' addParams(
    outdir: params.outdir,
    maxcpus: params.maxcpus,
    // halfmaxcpus: params.halfmaxcpus,
    medcpus: params.medcpus,
    halfmedcpus: params.halfmedcpus,
    container_bbtools: params.container_bbtools
)
include { 
    IVAR_CONSENSUS;
    IVAR_VARIANTS
} from '../modules/local/ivar' addParams(
    outdir: params.outdir,
    container_ivar: params.container_ivar,
    mincov: params.mincov
)
include { TRAFFIC_LIGHT_PLOT; } from '../modules/local/misc' addParams(
    outdir: params.outdir
)
include { SAMTOOLS_SORT_INDEX; } from '../modules/local/samtools' addParams(
    outdir: params.outdir,
    medcpus: params.medcpus,
    container_samtools: params.container_samtools
)

workflow analysis_core {
    take:
        paired_reads
        single_reads
        reference
        adapter_fasta
        primer_fasta
        gff_file_ivar

    main:
        // PRE-ALIGNMENT
        // =============
        INTERLEAVE(paired_reads)
        ch_take_viral = INTERLEAVE.out.interleaved_specimens | concat(single_reads)
        TAKE_VIRAL(ch_take_viral, reference)
        TRIM_ADAPTERS(TAKE_VIRAL.out.viral_specimens, adapter_fasta)
        TRIM_PRIMERS(TRIM_ADAPTERS.out.adapter_trimmed_specimens, primer_fasta)
        BBMERGE(TRIM_PRIMERS.out.primer_trimmed_specimens)

        BBMAP_ALIGN(BBMERGE.out.bbmerged_specimens, reference)

        // POST-ALIGNMENT
        // ==============
        REMOVE_JUNK_DELS(BBMAP_ALIGN.out.aligned_sams, reference)
        REMOVE_SINGLETONS(REMOVE_JUNK_DELS.out.filtered_junk_sams, reference)
        SOFT_TRIM(REMOVE_SINGLETONS.out.filtered_singletons_sams)
        PILEUP(SOFT_TRIM.out.soft_trimmed_sams)
        SUMMARIZE_COVERAGE(PILEUP.out.basecov.toSortedList())

        // SAMTOOLS AND IVAR
        // =================
        SAMTOOLS_SORT_INDEX(SOFT_TRIM.out.soft_trimmed_sams)
        IVAR_CONSENSUS(SAMTOOLS_SORT_INDEX.out.bamfile)
        IVAR_VARIANTS(SAMTOOLS_SORT_INDEX.out.bamfile, reference, gff_file_ivar)
        TRAFFIC_LIGHT_PLOT(IVAR_VARIANTS.out.ivar_tsvs.toSortedList(), SUMMARIZE_COVERAGE.out.coverage_summary)

    emit:
        combined_consensus_fasta = IVAR_CONSENSUS.out.consensus_collect |
            collectFile(
                name: "${file(params.outdir).getSimpleName()}_consensus.fasta",
                sort: { it.getSimpleName() },
                storeDir: "${params.outdir}"
            )
        coverage_summary = SUMMARIZE_COVERAGE.out.coverage_summary
}

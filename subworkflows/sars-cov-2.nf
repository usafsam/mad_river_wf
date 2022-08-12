include {
    LINEAGE_EXCEL;
    SPIKE_GENE_COVERAGE
} from '../modules/local/misc' addParams(
    outdir: params.outdir
)
include { NEXTCLADE; } from '../modules/local/nextclade' addParams(
    outdir: params.outdir,
    medcpus: params.medcpus,
    container_nextclade: params.container_nextclade,
    nextclade_options: params.nextclade_options
)
include { PANGOLIN; } from '../modules/local/pangolin' addParams(
    outdir: params.outdir,
    medcpus: params.medcpus,
    container_pangolin: params.container_pangolin,
    pangolin_options: params.pangolin_options
)
include { VADR; } from '../modules/local/vadr' addParams(
    outdir: params.outdir,
    medcpus: params.medcpus,
    container_vadr: params.container_vadr,
    vadr_options: params.vadr_options,
    vadr_reference: params.vadr_reference,
    vadr_trim_options: params.vadr_trim_options,
    vadr_mdir: params.vadr_mdir
)

workflow sarscov2 {
    take:
        combined_consensus_fasta
        reference_nextclade

    main:
        PANGOLIN(combined_consensus_fasta)
        NEXTCLADE(combined_consensus_fasta, reference_nextclade)
        VADR(combined_consensus_fasta)
        SPIKE_GENE_COVERAGE(NEXTCLADE.out.nextclade_csv)
        LINEAGE_EXCEL(
            PANGOLIN.out.pangolin,
            NEXTCLADE.out.nextclade_csv,
            VADR.out.vadr_file
        )
}

include { DATA_PREP; } from '../modules/local/misc' addParams(
    outdir: params.outdir,
    container_python: params.container_python
)
include { GRAB_NEXTCLADE_DATA; } from '../modules/local/nextclade' addParams(
    outdir: params.outdir,
    container_nextclade: params.container_nextclade,
    nextclade_reference_name: params.nextclade_reference_name
)

workflow analysis_prep {
    take:
        ivar_gff_gzip_file_checked
        primer_tsv

    main:
        GRAB_NEXTCLADE_DATA()
        DATA_PREP(ivar_gff_gzip_file_checked, primer_tsv)

    emit:
        reference_genome = GRAB_NEXTCLADE_DATA.out.reference_genome
        reference_nextclade = GRAB_NEXTCLADE_DATA.out.reference_nextclade
        gff_file_ivar = DATA_PREP.out.gff_file_ivar
        primer_fasta = DATA_PREP.out.primer_fasta
}

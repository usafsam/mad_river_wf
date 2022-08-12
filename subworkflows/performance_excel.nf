include { PERFORMANCE_EXCEL; } from '../modules/local/misc' addParams(
    outdir: params.outdir,
    container_python: params.container_python,
    protocol: params.protocol
)

workflow performance_excel {
    take:
        stats_json
        coverage_summary
        run_info_xml

    main:
        PERFORMANCE_EXCEL(stats_json, coverage_summary, run_info_xml)
}

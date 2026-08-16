#' TMScal Package Load Message
#'
#' @keywords internal
.onAttach <- function(libname, pkgname) {
    msg <- c(
        "\n╔═══════════════════════════════════════════════════════════════════╗",
        "║                                                                   ║",
        "║                    TMScal loaded successfully!                    ║",
        "║                                                                   ║",
        "║  Trinucleotide Mutation Spectrum for Prognostic Stratification    ║",
        "║                                                                   ║",
        "╚═══════════════════════════════════════════════════════════════════╝",
        "",
        "  Main Functions:",
        "    run_all_pipeline()       - Complete analysis pipeline",
        "    train()              - Train TMS model",
        "    predict()            - Single patient prediction",
        "",
        "  Optimization:",
        "    optimize_min_mutations() - Auto-optimize min_mutations",
        "    view_optimization()      - View optimization results",
        "",
        "  Data I/O:",
        "    convert_tcga_maf()     - convert MAF file to TMScal format",
        "    prepare_data()- convert MAF and clinical files to TMScal format",
        "",
        "  Documentation:",
        "    vignette('TMScal_intro') - Introduction vignette",
        "    ?run_tms_standard        - Function documentation",
        "    https://github.com/cxue/TMScal - GitHub repository",
        ""
    )
    packageStartupMessage(paste(msg, collapse = "\n"))
}
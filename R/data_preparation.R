#' Prepare data for TMScal analysis
#' 
#' @param maf_file Path to MAF file
#' @param clin_file Path to clinical file
#' @param sample_file Optional sample info file
#' @param output_dir Output directory
#' @return List with prepared data
#' @export
prepare_data <- function(maf_file, clin_file, sample_file = NULL,
                         output_dir = "./tmscal_prepared") {
    
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Load MAF
    log_message("Loading MAF file...")
    maf <- fread(maf_file, data.table = FALSE)
    
    # Detect format
    if (all(c("chromosome", "start_pos", "ref_allele", "alt_allele",
              "sample_id", "variant_class", "context") %in% colnames(maf))) {
        log_message("TMScal format detected")
    } else if (all(c("Hugo_Symbol", "Chromosome", "Start_Position",
                     "Reference_Allele", "Tumor_Seq_Allele2",
                     "Tumor_Sample_Barcode", "Variant_Classification",
                     "CONTEXT") %in% colnames(maf))) {
        log_message("TCGA format detected, converting...")
        maf <- convert_tcga_maf(maf)
    } else {
        stop("Unrecognized MAF format")
    }
    
    # Process clinical
    log_message("Processing clinical data...")
    clin <- process_clinical_data(clin_file)
    
    # Filter MSS if sample_file provided
    if (!is.null(sample_file)) {
        log_message("Filtering MSS samples...")
        sample_info <- fread(sample_file, data.table = FALSE, skip = 4)
        mss_patients <- sample_info$PATIENT_ID[
            sample_info$MSI_SCORE_MANTIS < 0.4 & 
            !is.na(sample_info$MSI_SCORE_MANTIS)
        ]
        clin <- clin[clin$sample_id %in% mss_patients, ]
    }
    
    # Ensure compatible output format (sample_id, os_time, os_event)
    clin_out <- data.frame(
        sample_id = clin$sample_id,
        os_time = clin$dfs_time,
        os_event = ifelse(clin$dfs_status == 1, 
                          "1:Recurred/Progressed", "0:DiseaseFree"),
        stringsAsFactors = FALSE
    )
    
    # Save
    fwrite(maf, file.path(output_dir, "maf_processed.txt"), sep = "\t")
    fwrite(clin_out, file.path(output_dir, "clin_processed.txt"), sep = "\t")
    
    log_message(sprintf("Prepared data saved to: %s", output_dir))
    log_message(sprintf("  MAF: %d rows", nrow(maf)))
    log_message(sprintf("  Clinical: %d samples", nrow(clin_out)))
    
    return(list(maf = maf, clin = clin_out))
}

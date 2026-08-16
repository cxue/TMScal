#' Check data quality for TMScal
#' 
#' @param maf MAF data
#' @param clin Clinical data
#' @return Quality metrics
#' @export
check_data_quality <- function(maf, clin) {
    
    n_samples <- nrow(clin)
    n_events <- sum(clin$dfs_status)
    event_rate <- n_events / n_samples
    
    muts_per_sample <- table(maf$sample_id)
    median_muts <- median(muts_per_sample)
    
    quality <- list(
        n_samples = n_samples,
        n_events = n_events,
        event_rate = event_rate,
        median_mutations = median_muts,
        suitable = n_samples >= 300 & n_events >= 50 & median_muts >= 20
    )
    
    cat("\n═══ Data Quality Check ═══\n")
    cat(sprintf("Samples: %d\n", n_samples))
    cat(sprintf("Events: %d (%.1f%%)\n", n_events, event_rate*100))
    cat(sprintf("Median mutations/sample: %.0f\n", median_muts))
    cat(sprintf("Suitable for TMScal: %s\n", ifelse(quality$suitable, "YES", "NO")))
    
    if (!quality$suitable) {
        cat("\nWarnings:\n")
        if (n_samples < 300) cat("  - Sample size < 300\n")
        if (n_events < 50) cat("  - Events < 50\n")
        if (median_muts < 20) cat("  - Low mutation burden\n")
    }
    
    return(quality)
}

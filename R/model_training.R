#' Train TMS model
#' 
#' @param maf_file Processed MAF file
#' @param clin_file Processed clinical file
#' @param output_dir Output directory
#' @return TMS model object
#' @export
train <- function(maf_file, clin_file, output_dir = "./tmscal_model") {
    
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Build pathway matrix
    log_message("Building pathway matrix...")
    maf <- fread(maf_file, data.table = FALSE)
    clin <- process_clinical_data(clin_file)
    clin$dfs_time <- as.numeric(clin$dfs_time)
    clin$dfs_status <- as.numeric(clin$dfs_status)
    
    pr <- build_pathway_matrix(maf)
    X <- pr$proportions
    common <- intersect(rownames(X), clin$sample_id)
    X <- X[common, , drop = FALSE]
    clin <- clin[match(common, clin$sample_id), ]
    
    # Split
    log_message("Splitting data...")
    split <- stratified_split(clin, 0.6, 42)
    X_train <- X[split$train, , drop = FALSE]
    clin_train <- clin[split$train, ]
    
    # Feature selection
    log_message("Selecting features...")
    stability <- compute_stability(X_train, clin_train)
    cox_res <- univariate_cox(X_train, clin_train)
    selection <- adaptive_select(stability, cox_res$p_value, cox_res$coef,
                                 X_train, clin_train)
    
    # Build model
    model <- list(
        pathways = colnames(X_train)[selection$idx],
        weights = cox_res$coef[selection$idx],
        tertiles = selection$train_tertiles,
        stability_cutoff = 0.07,
        p_cutoff = selection$p_cut,
        train_eval = selection$train_eval
    )
    names(model$weights) <- model$pathways
    class(model) <- "tmscal"
    
    saveRDS(model, file.path(output_dir, "model.rds"))
    
    return(model)
}
#' Predict TMS scores for new samples
#' 
#' @param maf_file MAF file for new samples
#' @param model Trained TMScal model
#' @return Data frame with predictions
#' @export
predict <- function(maf_file, model) {
    
    maf <- fread(maf_file, data.table = FALSE)
    pr <- build_pathway_matrix(maf)
    X <- pr$proportions
    
    scores <- calc_tms(X, model$pathways, model$weights)
    groups <- cut(scores,
                  breaks = c(-Inf, model$tertiles[1], model$tertiles[2], Inf),
                  labels = c("Low", "Medium", "High"))
    
    result <- data.frame(
        sample_id = rownames(X),
        TMS_score = scores,
        TMS_group = groups
    )
    
    return(result)
}
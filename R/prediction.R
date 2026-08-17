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
    scores[!is.finite(scores)] <- 0

    # 检查tertiles是否有效
    if (is.null(model$tertiles) || any(is.na(model$tertiles)) ||
        length(model$tertiles) != 2) {
        warning("Invalid tertiles in model, using internal tertiles")
        tertiles <- quantile(scores, probs=c(1/3, 2/3), na.rm=TRUE)
    } else {
        tertiles <- model$tertiles
    }

    groups <- cut(scores,
                  breaks = c(-Inf, tertiles[1], tertiles[2], Inf),
                  labels = c("Low", "Medium", "High"),
                  include.lowest = TRUE)

    result <- data.frame(
        sample_id = rownames(X),
        TMS_score = scores,
        TMS_group = groups
    )

    return(result)
}

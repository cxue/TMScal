#' Summary method for TMScal model
#'
#' @param object TMScal model object
#' @param ... Additional arguments
#' @return Invisible NULL
#' @export
summary.tmscal <- function(object, ...) {
    cat("\n")
    cat("═══ TMScal Model Summary ═══\n")
    cat("─────────────────────────────\n")
    
    # 兼容新旧字段名
    if (!is.null(object$selected_features)) {
        n_pathways <- length(object$selected_features)
        pathway_names <- object$selected_features
        pathway_weights <- object$feature_weights
    } else if (!is.null(object$pathways)) {
        n_pathways <- length(object$pathways)
        pathway_names <- object$pathways
        pathway_weights <- object$weights
    } else {
        n_pathways <- 0
        pathway_names <- character(0)
        pathway_weights <- numeric(0)
    }
    
    cat(sprintf("Pathways: %d\n", n_pathways))
    
    if (!is.null(object$freq_threshold)) {
        cat(sprintf("Frequency threshold: %.2f\n", object$freq_threshold))
    }
    if (!is.null(object$alpha)) {
        cat(sprintf("Alpha: %.2f\n", object$alpha))
    }
    if (!is.null(object$lambda)) {
        cat(sprintf("Lambda: %.6f\n", object$lambda))
    }
    
    cat("─────────────────────────────\n")
    
    if (!is.null(object$train_cindex)) {
        cat(sprintf("Train C-index: %.3f\n", object$train_cindex))
    }
    
    if (!is.null(object$train_tertiles)) {
        cat(sprintf("Train tertiles: %.4f, %.4f\n", 
                    object$train_tertiles[1], object$train_tertiles[2]))
    }
    
    cat("─────────────────────────────\n")
    
    if (n_pathways > 0) {
        cat("\nSelected Pathways:\n")
        for (i in 1:n_pathways) {
            cat(sprintf("  %2d. %s (weight=%.4f)\n", 
                        i, pathway_names[i], pathway_weights[i]))
        }
    }
    
    cat("\n")
    invisible(object)
}

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
    cat(sprintf("Pathways: %d\n", length(object$pathways)))
    cat(sprintf("Stability cutoff: %.2f\n", object$stability_cutoff))
    cat(sprintf("P cutoff: %.2f\n", object$p_cutoff))
    cat("─────────────────────────────\n")
    
    if (!is.null(object$train_eval)) {
        cat(sprintf("Train HR: %.2f (P=%.4f)\n", 
                    object$train_eval$hr, object$train_eval$hr_p))
        cat(sprintf("Train C-index: %.3f\n", object$train_eval$cindex))
        cat(sprintf("Train Logrank P: %.4f\n", object$train_eval$logrank_p))
        cat("─────────────────────────────\n")
    }
    
    if (!is.null(object$test_eval)) {
        cat(sprintf("Test HR: %.2f (P=%.4f)\n", 
                    object$test_eval$hr, object$test_eval$hr_p))
        cat(sprintf("Test C-index: %.3f\n", object$test_eval$cindex))
        cat(sprintf("Test Logrank P: %.4f\n", object$test_eval$logrank_p))
        cat("─────────────────────────────\n")
    }
    
    cat("\nSelected Pathways:\n")
    for (i in seq_along(object$pathways)) {
        cat(sprintf("  %2d. %s (weight=%.3f)\n", 
                    i, object$pathways[i], object$weights[i]))
    }
    cat("\n")
    
    invisible(object)
}

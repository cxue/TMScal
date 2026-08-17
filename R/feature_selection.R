# Adaptive feature selection

#' Adaptive feature selection based on stability and P-value
#' @param stability Vector of stability scores
#' @param p_value Vector of p-values
#' @param coefs Vector of Cox coefficients
#' @param X_train Training matrix
#' @param clin_train Training clinical data
#' @return Selection results
#' @keywords internal
adaptive_select <- function(stability, p_value, coefs, X_train, clin_train,
                            stab_cut=0.07, p_min=0.05, p_max=0.50,
                            min_pathways=5, max_pathways=40) {
    
    best_result <- NULL
    
    for (p_cut in seq(p_min, p_max, 0.05)) {
        selected_idx <- which(stability < stab_cut & p_value < p_cut)
        n_selected <- length(selected_idx)
        
        if (n_selected < min_pathways) next
        if (n_selected > max_pathways) next
        
        selected_pathways <- colnames(X_train)[selected_idx]
        selected_coefs <- coefs[selected_idx]
        names(selected_coefs) <- selected_pathways
        
        train_scores <- calc_tms(X_train, selected_pathways, selected_coefs)
        train_tertiles <- quantile(train_scores, probs=c(1/3, 2/3), na.rm=TRUE)
        train_eval <- evaluate_model(clin_train, train_scores, train_tertiles)
        
        score <- -log10(train_eval$logrank_p + 1e-10) + log(max(train_eval$hr, 0.1))
        
        if (is.null(best_result) || score > best_result$score) {
            best_result <- list(
                idx = selected_idx,
                n = n_selected,
                p_cut = p_cut,
                score = score,
                train_eval = train_eval,
                train_scores = train_scores,
                train_tertiles = train_tertiles
            )
        }
    }
    
    if (is.null(best_result)) {
        selected_idx <- which(stability < stab_cut & p_value < 0.50)
        selected_pathways <- colnames(X_train)[selected_idx]
        selected_coefs <- coefs[selected_idx]
        names(selected_coefs) <- selected_pathways
        train_scores <- calc_tms(X_train, selected_pathways, selected_coefs)
        train_tertiles <- quantile(train_scores, probs=c(1/3, 2/3), na.rm=TRUE)
        
        best_result <- list(
            idx = selected_idx,
            n = length(selected_idx),
            p_cut = 0.50,
            score = NA,
            train_eval = evaluate_model(clin_train, train_scores, train_tertiles),
            train_scores = train_scores,
            train_tertiles = train_tertiles
        )
    }
    
    return(best_result)
}

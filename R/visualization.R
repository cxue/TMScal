#' Plot KM curve for TMS groups
#' 
#' @param model TMScal model
#' @param clin_data Clinical data with TMS scores
#' @param output_file Output file path
#' @return ggplot object
#' @export
plot_km_curve <- function(model, clin_data, output_file = NULL) {
    
    suppressPackageStartupMessages({
        library(survival)
        library(survminer)
        library(ggplot2)
    })
    
    # Calculate scores if not present
    if (!"TMS_score" %in% colnames(clin_data)) {
        scores <- predict_tscores(clin_data, model)
        clin_data$TMS_score <- scores
    }
    
    # Group
    clin_data$TMS_group <- cut(clin_data$TMS_score,
                                breaks = c(-Inf, model$tertiles[1], 
                                          model$tertiles[2], Inf),
                                labels = c("Low", "Medium", "High"))
    
    # Fit KM
    fit <- survfit(Surv(dfs_time, dfs_status) ~ TMS_group, data = clin_data)
    
    # Colors
    TMS_COLORS <- c("Low" = "#2166AC", "Medium" = "#FDDBC7", "High" = "#B2182B")
    
    # Plot
    p <- ggsurvplot(fit, data = clin_data,
                    palette = TMS_COLORS,
                    pval = TRUE,
                    xlab = "Time (months)",
                    ylab = "Survival Probability",
                    legend.title = "TMS Group",
                    legend.labs = c("Low", "Medium", "High"),
                    ggtheme = theme_minimal(base_size = 12),
                    risk.table = TRUE)
    
    if (!is.null(output_file)) {
        ggsave(output_file, p$plot, width = 8, height = 6, dpi = 300)
    }
    
    return(p)
}

#' Plot pathway stability distribution
#' 
#' @param stability Vector of stability scores
#' @param cutoff Stability cutoff
#' @param output_file Output file
#' @export
plot_stability <- function(stability, cutoff = 0.07, output_file = NULL) {
    
    library(ggplot2)
    
    df <- data.frame(stability = stability)
    
    p <- ggplot(df, aes(x = stability)) +
        geom_histogram(bins = 50, fill = "lightblue", color = "black") +
        geom_vline(xintercept = cutoff, color = "red", linetype = "dashed") +
        labs(title = "Pathway Stability Distribution",
             x = "Stability Score", y = "Count") +
        theme_minimal()
    
    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 8, height = 5, dpi = 300)
    }
    
    return(p)
}

#' Plot feature importance
#' 
#' @param model TMScal model
#' @param output_file Output file
#' @export
plot_feature_importance <- function(model, output_file = NULL) {
    
    library(ggplot2)
    
    df <- data.frame(
        pathway = model$pathways,
        weight = model$weights
    )
    df$pathway <- factor(df$pathway, levels = rev(df$pathway))
    
    p <- ggplot(df, aes(x = pathway, y = weight, fill = weight > 0)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
        labs(title = "Feature Contributions to TMS",
             x = "Pathway", y = "Weight") +
        theme_minimal() +
        theme(legend.position = "none")
    
    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 7, height = 6, dpi = 300)
    }
    
    return(p)
}

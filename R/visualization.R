#' Plot KM curve for TMS groups
#' 
#' @param prediction_file Path to prediction CSV file (sample_id, TMS_score, TMS_group)
#' @param clin_file Path to clinical data file
#' @param output_file Output file path
#' @return ggsurvplot object
#' @export
plot_km_curve <- function(prediction_file, clin_file, output_file = NULL) {
    
    suppressPackageStartupMessages({
        library(survival)
        library(survminer)
        library(ggplot2)
        library(data.table)
    })
    
    # Load prediction
    predictions <- fread(prediction_file, data.table = FALSE)
    
    # Load clinical data
    clin_data <- fread(clin_file, data.table = FALSE)
    colnames(clin_data) <- tolower(colnames(clin_data))
    clin_data$dfs_time <- as.numeric(clin_data$os_time)
    clin_data$dfs_status <- ifelse(grepl("1:|Recurred|Progressed", clin_data$os_event), 1, 0)
    clin_data <- clin_data[!is.na(clin_data$dfs_time) & clin_data$dfs_time > 0, ]
    
    # Merge
    clin_data <- merge(clin_data, predictions[, c("sample_id", "TMS_score", "TMS_group")],
                       by = "sample_id", all.x = TRUE)
    
    # Factor
    clin_data$TMS_group <- factor(clin_data$TMS_group, 
                                   levels = c("Low", "Medium", "High"))
    
    # Fit KM
    fit <- survfit(Surv(dfs_time, dfs_status) ~ TMS_group, data = clin_data)
    
    # Colors
    TMS_COLORS <- c("Low" = "#2166AC", "Medium" = "#FDDBC7", "High" = "#B2182B")
    
    # Calculate HR
    hl <- clin_data$TMS_group %in% c("Low", "High")
    clin_hl <- clin_data[hl, ]
    groups_hl <- factor(clin_hl$TMS_group, levels = c("Low", "High"))
    fit_hl <- coxph(Surv(dfs_time, dfs_status) ~ groups_hl, data = clin_hl)
    hr <- exp(coef(fit_hl))
    hr_ci <- exp(confint(fit_hl))
    hr_p <- summary(fit_hl)$coefficients[1, "Pr(>|z|)"]
    
    # Plot
    p <- ggsurvplot(fit, data = clin_data,
                    palette = TMS_COLORS,
                    pval = TRUE,
                    pval.size = 4,
                    xlab = "Time (months)",
                    ylab = "DFS Probability",
                    legend.title = "TMS Group",
                    legend.labs = c("Low", "Medium", "High"),
                    ggtheme = theme_minimal(base_size = 12),
                    risk.table = TRUE,
                    risk.table.col = "strata",
                    risk.table.height = 0.25,
                    conf.int = FALSE)
    
    # Add HR annotation
    hr_text <- sprintf("High vs Low: HR = %.2f (%.2f-%.2f)\nP = %.4f",
                       hr, hr_ci[1], hr_ci[2], hr_p)
    
    p$plot <- p$plot + 
        ggplot2::annotate("text", x = 15, y = 0.2, label = hr_text,
                          size = 4, hjust = 0)
    
    if (!is.null(output_file)) {
        # 根据扩展名保存
        if (grepl("\\.pdf$", output_file)) {
            pdf(output_file, width = 8, height = 6)
            print(p)
            dev.off()
        } else {
            png(output_file, width = 8, height = 6, units = "in", res = 300)
            print(p)
            dev.off()
        }
        cat(sprintf("  KM curve saved to: %s\n", output_file))
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
    
    # 兼容新旧字段名
    if (!is.null(model$selected_features)) {
        pathways <- model$selected_features
        weights <- model$feature_weights
    } else if (!is.null(model$pathways)) {
        pathways <- model$pathways
        weights <- model$weights
    } else {
        stop("Model does not contain features")
    }
    
    df <- data.frame(pathway = pathways, weight = weights)
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

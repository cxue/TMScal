#' Predict TMS scores for new samples
#'
#' @param maf_file MAF file for new samples
#' @param model Trained TMScal model
#' @param output_file Optional output CSV file path
#' @return Data frame with predictions
#' @export
predict <- function(maf_file, model, output_file = NULL) {

    maf <- fread(maf_file, data.table = FALSE)
    mat <- build_raw_matrix(maf)
    
    # TMB
    tmb <- rowSums(mat)
    tmb_log <- log10(tmb + 1)
    tmb_scaled <- (tmb_log - model$tmb_mean) / model$tmb_sd
    
    # 使用模型的高频pathway
    available_pathways <- intersect(model$high_freq_pathways, colnames(mat))
    X <- log10(mat[, available_pathways, drop = FALSE] + 1)
    X <- cbind(X, TMB = tmb_scaled)
    
    # 选择模型特征
    selected <- model$selected_features
    available_selected <- intersect(selected, colnames(X))
    
    if (length(available_selected) == 0) {
        stop("No selected features available in test data")
    }
    
    # 获取权重
    weights <- model$feature_weights[match(available_selected, selected)]
    
    # 计算分数
    scores <- X[, available_selected, drop = FALSE] %*% weights
    
    # 标准化（使用训练集统计量）
    scores_scaled <- (scores - model$train_mean) / model$train_sd
    
    # 分组（使用训练集tertiles）
    groups <- cut(scores_scaled,
                  breaks = c(-Inf, model$train_tertiles[1], model$train_tertiles[2], Inf),
                  labels = c("Low", "Medium", "High"),
                  include.lowest = TRUE)
    
    result <- data.frame(
        sample_id = rownames(mat),
        TMS_score = as.numeric(scores_scaled),
        TMS_group = groups,
        stringsAsFactors = FALSE
    )
    
    # 保存结果
    if (!is.null(output_file)) {
        fwrite(result, output_file)
        log_message(sprintf("Predictions saved to: %s", output_file))
    }
    
    return(result)
}

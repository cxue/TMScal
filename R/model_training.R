#' Train TMS model with internal CV tuning
#'
#' @param maf_file Processed MAF file (training set)
#' @param clin_file Processed clinical file (training set)
#' @param output_dir Output directory
#' @param seed Random seed for reproducibility
#' @return TMS model object
#' @export
train <- function(maf_file, clin_file, output_dir = "./tmscal_model", seed = 42) {

    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # Load data
    log_message("Loading training data...")
    maf <- fread(maf_file, data.table = FALSE)
    clin <- fread(clin_file, data.table = FALSE)
    colnames(clin) <- tolower(colnames(clin))
    clin$dfs_time <- as.numeric(clin$os_time)
    clin$dfs_status <- ifelse(grepl("1:|Recurred|Progressed", clin$os_event), 1, 0)
    clin <- clin[!is.na(clin$dfs_time) & clin$dfs_time > 0, ]

    # Build pathway matrix
    log_message("Building pathway matrix...")
    mat <- build_raw_matrix(maf)
    common <- intersect(rownames(mat), clin$sample_id)
    mat <- mat[common, , drop = FALSE]
    clin <- clin[match(common, clin$sample_id), ]

    y_train <- Surv(clin$dfs_time, clin$dfs_status)

    # TMB
    tmb_train <- rowSums(mat)
    tmb_log <- log10(tmb_train + 1)

    log_message(sprintf("Training: %d samples × %d pathways", nrow(mat), ncol(mat)))

    # ============================================
    # Internal CV: frequency threshold
    # ============================================
    log_message("Internal CV: frequency threshold...")
    
    freq_thresholds <- c(0.05, 0.08, 0.10, 0.12, 0.15)
    threshold_results <- data.frame(
        threshold = freq_thresholds,
        cv_cindex = NA,
        n_features = NA
    )
    
    set.seed(seed)  # 使用外部seed
    folds <- sample(rep(1:5, length.out = nrow(mat)))
    
    for (i in 1:length(freq_thresholds)) {
        thresh <- freq_thresholds[i]
        pathway_freq <- apply(mat > 0, 2, mean)
        high_freq <- pathway_freq > thresh
        X_freq <- log10(mat[, high_freq, drop = FALSE] + 1)
        
        fold_cindices <- numeric(5)
        
        for (fold in 1:5) {
            tr_idx <- which(folds != fold)
            val_idx <- which(folds == fold)
            
            X_tr <- X_freq[tr_idx, , drop = FALSE]
            X_val <- X_freq[val_idx, , drop = FALSE]
            
            tmb_mean_f <- mean(tmb_log[tr_idx])
            tmb_sd_f <- sd(tmb_log[tr_idx])
            X_tr <- cbind(X_tr, TMB = (tmb_log[tr_idx] - tmb_mean_f) / tmb_sd_f)
            X_val <- cbind(X_val, TMB = (tmb_log[val_idx] - tmb_mean_f) / tmb_sd_f)
            
            set.seed(seed + 100)  # 变体seed
            cv_fold <- tryCatch({
                cv.glmnet(X_tr, y_train[tr_idx], family = "cox", 
                          alpha = 0.1, nfolds = 5, cox.ties = "breslow")
            }, error = function(e) NULL)
            
            if (!is.null(cv_fold)) {
                coef_fold <- as.numeric(coef(cv_fold, s = "lambda.min"))
                sel_fold <- which(coef_fold != 0)
                if (length(sel_fold) > 0) {
                    scores_val <- X_val[, sel_fold, drop = FALSE] %*% coef_fold[sel_fold]
                    fold_cindices[fold] <- calculate_cindex(y_train[val_idx], scores_val)
                }
            }
        }
        
        threshold_results$cv_cindex[i] <- mean(fold_cindices, na.rm = TRUE)
        threshold_results$n_features[i] <- sum(high_freq)
        log_message(sprintf("  threshold %.2f: CV=%.3f, features=%d",
                            thresh, threshold_results$cv_cindex[i],
                            threshold_results$n_features[i]))
    }
    
    best_threshold_idx <- which.max(threshold_results$cv_cindex)
    best_threshold <- threshold_results$threshold[best_threshold_idx]
    log_message(sprintf("Best threshold: %.2f", best_threshold))

    # ============================================
    # Internal CV: alpha
    # ============================================
    log_message("Internal CV: alpha...")
    
    alpha_values <- c(0.05, 0.10, 0.20, 0.30, 0.50)
    alpha_results <- data.frame(alpha = alpha_values, cv_cindex = NA)
    
    pathway_freq <- apply(mat > 0, 2, mean)
    high_freq <- pathway_freq > best_threshold
    X_base <- log10(mat[, high_freq, drop = FALSE] + 1)
    
    for (i in 1:length(alpha_values)) {
        alpha_val <- alpha_values[i]
        fold_cindices <- numeric(5)
        
        for (fold in 1:5) {
            tr_idx <- which(folds != fold)
            val_idx <- which(folds == fold)
            
            X_tr <- X_base[tr_idx, , drop = FALSE]
            X_val <- X_base[val_idx, , drop = FALSE]
            
            tmb_mean_f <- mean(tmb_log[tr_idx])
            tmb_sd_f <- sd(tmb_log[tr_idx])
            X_tr <- cbind(X_tr, TMB = (tmb_log[tr_idx] - tmb_mean_f) / tmb_sd_f)
            X_val <- cbind(X_val, TMB = (tmb_log[val_idx] - tmb_mean_f) / tmb_sd_f)
            
            set.seed(seed + 200)  # 变体seed
            cv_fold <- tryCatch({
                cv.glmnet(X_tr, y_train[tr_idx], family = "cox", 
                          alpha = alpha_val, nfolds = 5, cox.ties = "breslow")
            }, error = function(e) NULL)
            
            if (!is.null(cv_fold)) {
                coef_fold <- as.numeric(coef(cv_fold, s = "lambda.min"))
                sel_fold <- which(coef_fold != 0)
                if (length(sel_fold) > 0) {
                    scores_val <- X_val[, sel_fold, drop = FALSE] %*% coef_fold[sel_fold]
                    fold_cindices[fold] <- calculate_cindex(y_train[val_idx], scores_val)
                }
            }
        }
        
        alpha_results$cv_cindex[i] <- mean(fold_cindices, na.rm = TRUE)
        log_message(sprintf("  alpha %.2f: CV=%.3f", alpha_val, alpha_results$cv_cindex[i]))
    }
    
    best_alpha_idx <- which.max(alpha_results$cv_cindex)
    best_alpha <- alpha_results$alpha[best_alpha_idx]
    log_message(sprintf("Best alpha: %.2f", best_alpha))

    # ============================================
    # Final lambda search
    # ============================================
    log_message("Determining lambda...")
    
    X_final <- cbind(X_base, TMB = (tmb_log - mean(tmb_log)) / sd(tmb_log))
    
    set.seed(seed + 300)  # 变体seed
    cv_final <- cv.glmnet(X_final, y_train, family = "cox", 
                          alpha = best_alpha, nfolds = 5, cox.ties = "breslow")
    
    lambda_min <- cv_final$lambda.min
    lambda_1se <- cv_final$lambda.1se
    lambda_range <- exp(seq(log(lambda_min), log(lambda_1se), length.out = 10))
    
    lambda_results <- data.frame()
    
    for (lam in lambda_range) {
        coef_val <- as.numeric(coef(cv_final, s = lam))
        n_feat <- sum(coef_val != 0)
        
        if (n_feat >= 10 && n_feat <= 60) {
            closest_idx <- which.min(abs(cv_final$lambda - lam))
            cv_cidx <- cv_final$cvm[closest_idx]
            lambda_results <- rbind(lambda_results, data.frame(
                lambda = lam, n_features = n_feat, cv_cindex = cv_cidx
            ))
        }
    }
    
    valid_lambda <- lambda_results[
        lambda_results$n_features >= 15 & lambda_results$n_features <= 40, ]
    
    if (nrow(valid_lambda) > 0) {
        best_lambda <- valid_lambda$lambda[which.max(valid_lambda$cv_cindex)]
    } else {
        best_lambda <- lambda_min
    }
    
    log_message(sprintf("Best lambda: %.6f", best_lambda))

    # ============================================
    # Build final model
    # ============================================
    log_message("Building final model...")
    
    coef_final <- as.numeric(coef(cv_final, s = best_lambda))
    names(coef_final) <- rownames(coef(cv_final, s = best_lambda))
    selected_final <- which(coef_final != 0)
    
    feature_names <- names(coef_final)[selected_final]
    feature_weights <- coef_final[selected_final]
    
    log_message(sprintf("Selected %d features", length(selected_final)))
    
    train_scores <- X_final[, selected_final, drop = FALSE] %*% feature_weights
    train_mean <- mean(train_scores)
    train_sd <- sd(train_scores)
    train_scores_scaled <- (train_scores - train_mean) / train_sd
    train_tertiles <- quantile(train_scores_scaled, probs = c(1/3, 2/3))
    train_cindex <- calculate_cindex(y_train, train_scores_scaled)
    
    log_message(sprintf("Train C-index: %.3f", train_cindex))
    
    model <- list(
        freq_threshold = best_threshold,
        alpha = best_alpha,
        lambda = best_lambda,
        selected_features = feature_names,
        feature_weights = feature_weights,
        train_mean = train_mean,
        train_sd = train_sd,
        train_tertiles = train_tertiles,
        tmb_mean = mean(tmb_log),
        tmb_sd = sd(tmb_log),
        high_freq_pathways = names(pathway_freq)[high_freq],
        train_cindex = train_cindex,
        seed = seed,
        created_date = Sys.Date()
    )
    
    class(model) <- "tmscal"
    
    saveRDS(model, file.path(output_dir, "model.rds"))
    log_message(sprintf("Model saved to: %s/model.rds", output_dir))
    
    return(model)
}

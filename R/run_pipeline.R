#' Run complete TMScal pipeline
#'
#' @param maf_file MAF file
#' @param clin_file Clinical file
#' @param output_dir Output directory
#' @param split_ratio Train ratio (default 0.7)
#' @param seed Random seed
#' @return List with model and predictions
#' @export
run_all_pipeline <- function(maf_file, clin_file, output_dir = "./tmscal_output",
                             split_ratio = 0.7, seed = 764) {

    log_message("Starting TMScal pipeline...")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # Step 1: Split data
    log_message("Step 1/4: Splitting data...")
    
    maf_all <- fread(maf_file, data.table = FALSE)
    clin_all <- fread(clin_file, data.table = FALSE)
    colnames(clin_all) <- tolower(colnames(clin_all))
    clin_all$dfs_time <- as.numeric(clin_all$os_time)
    clin_all$dfs_status <- ifelse(grepl("1:|Recurred|Progressed", clin_all$os_event), 1, 0)
    clin_all <- clin_all[!is.na(clin_all$dfs_time) & clin_all$dfs_time > 0, ]

    set.seed(seed)
    train_samples <- c(); test_samples <- c()
    for (s in unique(clin_all$dfs_status)) {
        idx <- which(clin_all$dfs_status == s)
        n <- length(idx)
        n_train <- round(n * split_ratio)
        shuffled <- sample(clin_all$sample_id[idx])
        train_samples <- c(train_samples, shuffled[1:n_train])
        test_samples <- c(test_samples, shuffled[(n_train+1):n])
    }
    
    # Save split files
    train_dir <- file.path(output_dir, "train")
    test_dir <- file.path(output_dir, "test")
    dir.create(train_dir, showWarnings = FALSE)
    dir.create(test_dir, showWarnings = FALSE)
    
    fwrite(maf_all[maf_all$sample_id %in% train_samples, ], 
           file.path(train_dir, "maf_train.txt"), sep = "\t")
    fwrite(clin_all[clin_all$sample_id %in% train_samples, c("sample_id", "os_time", "os_event")], 
           file.path(train_dir, "clin_train.txt"), sep = "\t")
    fwrite(maf_all[maf_all$sample_id %in% test_samples, ], 
           file.path(test_dir, "maf_test.txt"), sep = "\t")
    fwrite(clin_all[clin_all$sample_id %in% test_samples, c("sample_id", "os_time", "os_event")], 
           file.path(test_dir, "clin_test.txt"), sep = "\t")
    
    log_message(sprintf("  Train: %d, Test: %d", length(train_samples), length(test_samples)))

    # Step 2: Train model
    log_message("Step 2/4: Training model...")
    model <- train(
        maf_file = file.path(train_dir, "maf_train.txt"),
        clin_file = file.path(train_dir, "clin_train.txt"),
        output_dir = file.path(output_dir, "model"),
        seed = seed
    )

    # Step 3: Predict test set
    log_message("Step 3/4: Predicting test set...")
    predictions <- predict(
        maf_file = file.path(test_dir, "maf_test.txt"),
        model = model
    )

    # Step 4: Evaluate
    log_message("Step 4/4: Evaluating...")
    
    clin_test <- fread(file.path(test_dir, "clin_test.txt"), data.table = FALSE)
    colnames(clin_test) <- tolower(colnames(clin_test))
    clin_test$dfs_time <- as.numeric(clin_test$os_time)
    clin_test$dfs_status <- ifelse(grepl("1:|Recurred|Progressed", clin_test$os_event), 1, 0)
    clin_test <- clin_test[!is.na(clin_test$dfs_time) & clin_test$dfs_time > 0, ]
    
    clin_merged <- merge(clin_test, predictions[, c("sample_id", "TMS_score", "TMS_group")],
                         by = "sample_id", all.x = TRUE)
    clin_merged$TMS_group <- factor(clin_merged$TMS_group, 
                                     levels = c("Low", "Medium", "High"))
    
    # HR
    hl <- clin_merged$TMS_group %in% c("Low", "High")
    clin_hl <- clin_merged[hl, ]
    groups_hl <- factor(clin_hl$TMS_group, levels = c("Low", "High"))
    fit_hl <- coxph(Surv(dfs_time, dfs_status) ~ groups_hl, data = clin_hl)
    hr <- exp(coef(fit_hl))
    hr_ci <- exp(confint(fit_hl))
    hr_p <- summary(fit_hl)$coefficients[1, "Pr(>|z|)"]
    
    # C-index
    fit_cont <- coxph(Surv(dfs_time, dfs_status) ~ TMS_score, data = clin_merged)
    cindex <- summary(fit_cont)$concordance[1]
    
    # Log-rank
    logrank <- survdiff(Surv(dfs_time, dfs_status) ~ TMS_group, data = clin_merged)
    logrank_p <- 1 - pchisq(logrank$chisq, 2)
    
    log_message(sprintf("  Test: HR=%.2f (P=%.4f), C-index=%.3f", hr, hr_p, cindex))
    
    # Save predictions
    fwrite(clin_merged, file.path(output_dir, "test_predictions.csv"))
    
    return(list(
        model = model,
        predictions = predictions,
        clin_test = clin_merged,
        evaluation = list(hr = hr, hr_ci = hr_ci, hr_p = hr_p, 
                          cindex = cindex, logrank_p = logrank_p)
    ))
}

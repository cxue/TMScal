# Pathway analysis core functions

#' Build pathway mutation matrix
#' @param maf MAF data frame
#' @param min_mutations Minimum mutations per sample
#' @return List with proportions and raw counts
#' @keywords internal
build_pathway_matrix <- function(maf, min_mutations = 0) {
    maf_missense <- maf[grepl("Missense", maf$variant_class, ignore.case=TRUE), ]
    if (nrow(maf_missense) == 0) stop("No missense mutations found")
    
    maf_missense$pathway <- mapply(get_pathway,
                                   maf_missense$ref_allele,
                                   maf_missense$alt_allele,
                                   maf_missense$context)
    maf_missense <- maf_missense[!is.na(maf_missense$pathway), ]
    
    pathway_names <- generate_pathway_names()
    counts <- as.data.table(table(maf_missense$sample_id, maf_missense$pathway))
    setnames(counts, c("sample_id", "pathway", "n"))
    
    patients <- unique(counts$sample_id)
    mat <- matrix(0, nrow=length(patients), ncol=192)
    rownames(mat) <- patients
    colnames(mat) <- pathway_names
    
    for (i in 1:nrow(counts)) {
        if (counts$pathway[i] %in% pathway_names) {
            mat[counts$sample_id[i], counts$pathway[i]] <- counts$n[i]
        }
    }
    
    total_muts <- rowSums(mat)
    valid <- total_muts >= min_mutations
    mat <- mat[valid, , drop=FALSE]
    
    if (nrow(mat) == 0) stop("No samples pass filter")
    
    epsilon <- 0.5
    mat_eps <- mat + epsilon
    props <- mat_eps / rowSums(mat_eps)
    props[props <= 0] <- .Machine$double.eps
    
    clr <- t(apply(props, 1, function(row) {
        gm <- exp(mean(log(row)))
        log(row / gm)
    }))
    clr[!is.finite(clr)] <- 0
    
    return(list(proportions=clr, patients=rownames(mat),
                raw_counts=mat, mutation_counts=total_muts[valid]))
}

#' Process clinical data
#' @param clin_file Path to clinical data file
#' @return Processed clinical data frame
#' @keywords internal
process_clinical_data <- function(clin_file) {
    # 检测文件格式（是否TCGA格式，有#注释行）
    first_line <- readLines(clin_file, n = 1)
    
    if (grepl("^#", first_line)) {
        # TCGA格式：前4行是注释
        clin <- fread(clin_file, data.table = FALSE, skip = 4)
    } else {
        # 标准格式
        clin <- fread(clin_file, data.table = FALSE)
    }
    
    colnames(clin) <- tolower(colnames(clin))
    
    # 识别样本ID列
    id_col <- NULL
    for (candidate in c("sample_id", "patient_id")) {
        if (candidate %in% colnames(clin)) {
            id_col <- candidate
            break
        }
    }
    if (is.null(id_col)) stop("Cannot identify sample ID column")
    clin$sample_id <- clin[[id_col]]
    
    # 识别生存时间列
    time_col <- NULL
    for (candidate in c("dfs_time", "dfs_months", "os_time", "os_months")) {
        if (candidate %in% colnames(clin)) {
            time_col <- candidate
            break
        }
    }
    if (is.null(time_col)) stop("Cannot identify survival time column")
    clin$dfs_time <- as.numeric(clin[[time_col]])
    
    # 识别生存状态列
    status_col <- NULL
    for (candidate in c("dfs_status", "dfs_event", "os_status", "os_event")) {
        if (candidate %in% colnames(clin)) {
            status_col <- candidate
            break
        }
    }
    if (is.null(status_col)) stop("Cannot identify survival status column")
    clin$dfs_status <- ifelse(
        grepl("^1|DECEASED|DEAD|RECURR|PROGRESS", 
              toupper(clin[[status_col]])),
        1, 0
    )
    
    clin <- clin[!is.na(clin$dfs_time) & clin$dfs_time > 0, ]
    clin <- clin[clin$sample_id != "", ]
    
    return(clin)
}

#' Compute pathway stability via bootstrap
#' @param X Pathway proportions matrix
#' @param clin Clinical data
#' @param n_bootstrap Number of bootstrap iterations
#' @param seed Random seed
#' @return Vector of stability scores
#' @keywords internal
compute_stability <- function(X, clin, n_bootstrap=200, seed=42) {
    set.seed(seed)
    bootstrap_means <- matrix(NA, nrow=n_bootstrap, ncol=ncol(X))
    colnames(bootstrap_means) <- colnames(X)
    
    for (b in 1:n_bootstrap) {
        boot_idx <- c()
        for (s in unique(clin$dfs_status)) {
            idx <- which(clin$dfs_status == s)
            boot_idx <- c(boot_idx, sample(idx, length(idx), replace=TRUE))
        }
        bootstrap_means[b, ] <- colMeans(X[boot_idx, , drop=FALSE])
    }
    
    stability <- apply(bootstrap_means, 2, sd) / (abs(colMeans(X)) + 0.1)
    stability[!is.finite(stability)] <- 999
    return(stability)
}

#' Univariate Cox regression
#' @param X Pathway proportions matrix
#' @param clin Clinical data
#' @return List with coefficients and p-values
#' @keywords internal
univariate_cox <- function(X, clin) {
    y <- Surv(clin$dfs_time, clin$dfs_status)
    coefs <- numeric(ncol(X))
    p_values <- numeric(ncol(X))
    
    for (i in 1:ncol(X)) {
        x <- X[, i]
        if (sd(x, na.rm=TRUE) < 1e-10) {
            coefs[i] <- 0
            p_values[i] <- 1
            next
        }
        fit <- tryCatch(coxph(y ~ x), error=function(e) NULL)
        if (!is.null(fit)) {
            coefs[i] <- coef(fit)
            p_values[i] <- summary(fit)$coefficients[1, "Pr(>|z|)"]
        } else {
            coefs[i] <- 0
            p_values[i] <- 1
        }
    }
    
    return(list(coef=coefs, p_value=p_values))
}

#' Stratified split
#' @param clin Clinical data
#' @param ratio Train ratio
#' @param seed Random seed
#' @return List with train and test indices
#' @keywords internal
stratified_split <- function(clin, ratio, seed) {
    set.seed(seed)
    train_idx <- c(); test_idx <- c()
    for (s in unique(clin$dfs_status)) {
        idx <- which(clin$dfs_status == s)
        n_train <- round(length(idx) * ratio)
        train_s <- sample(idx, n_train)
        test_s <- setdiff(idx, train_s)
        train_idx <- c(train_idx, train_s)
        test_idx <- c(test_idx, test_s)
    }
    return(list(train=train_idx, test=test_idx))
}

#' Calculate TMS score
#' @param X Pathway proportions matrix
#' @param pathways Selected pathway names
#' @param weights Pathway weights
#' @return Vector of TMS scores
#' @keywords internal
calc_tms <- function(X, pathways, weights) {
    names(weights) <- pathways
    scores <- rep(0, nrow(X))
    for (p in pathways) {
        if (p %in% colnames(X)) {
            scores <- scores + X[, p] * weights[p]
        }
    }
    scores[!is.finite(scores)] <- 0
    return(scores)
}

#' Evaluate model performance
#' @param clin Clinical data
#' @param scores TMS scores
#' @param tertiles Cutoffs
#' @return Evaluation metrics
#' @keywords internal
evaluate_model <- function(clin, scores, tertiles) {
    groups <- cut(scores,
                  breaks=c(-Inf, tertiles[1], tertiles[2], Inf),
                  labels=c("Low", "Medium", "High"))
    clin$group <- factor(groups, levels=c("Low", "Medium", "High"))
    
    hl <- clin$group %in% c("Low", "High")
    clin_hl <- clin[hl, ]
    groups_hl <- factor(clin_hl$group, levels=c("Low", "High"))
    fit_hl <- coxph(Surv(dfs_time, dfs_status) ~ groups_hl, data=clin_hl)
    hr <- exp(coef(fit_hl))
    hr_ci <- exp(confint(fit_hl))
    hr_p <- summary(fit_hl)$coefficients[1, "Pr(>|z|)"]
    
    fit_cont <- coxph(Surv(dfs_time, dfs_status) ~ scores, data=clin)
    cindex <- summary(fit_cont)$concordance[1]
    
    logrank <- survdiff(Surv(dfs_time, dfs_status) ~ group, data=clin)
    logrank_p <- 1 - pchisq(logrank$chisq, 2)
    
    trend_fit <- coxph(Surv(dfs_time, dfs_status) ~ as.numeric(group), data=clin)
    trend_p <- summary(trend_fit)$coefficients[,"Pr(>|z|)"]
    
    event_rates <- sapply(c("Low", "Medium", "High"), function(g) {
        sub <- clin[clin$group == g, ]
        if (nrow(sub) > 0) sum(sub$dfs_status)/nrow(sub) else NA
    })
    
    return(list(hr=hr, hr_ci=hr_ci, hr_p=hr_p, cindex=cindex,
                logrank_p=logrank_p, trend_p=trend_p, 
                event_rates=event_rates, group_sizes=table(clin$group)))
}

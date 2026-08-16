#' Convert TCGA MAF to TMScal format
#' 
#' @param maf TCGA format MAF
#' @return TMScal format MAF
#' @export
convert_tcga_maf <- function(maf) {
    
    # Extract missense
    missense <- maf[grepl("Missense", maf$Variant_Classification, ignore.case=TRUE), ]
    
    # Extract 11-nt context
    extract_11nt <- function(context, ref) {
        if (is.na(context) || is.na(ref)) return(NA)
        if (nchar(ref) != 1 || !ref %in% c("A","C","G","T")) return(NA)
        if (nchar(context) < 11 || nchar(context) > 100) return(NA)
        
        ref_pos <- regexpr(ref, context, fixed=TRUE)[1]
        if (ref_pos < 6) return(NA)
        
        start <- ref_pos - 5
        end <- ref_pos + 5
        if (start < 1 || end > nchar(context)) return(NA)
        
        return(substr(context, start, end))
    }
    
    missense$context_len <- nchar(missense$CONTEXT)
    valid <- missense$context_len >= 11 & missense$context_len <= 100
    
    missense$context_11 <- NA
    missense$context_11[valid] <- mapply(
        extract_11nt,
        missense$CONTEXT[valid],
        missense$Reference_Allele[valid]
    )
    
    # Build TMScal format
    result <- data.frame(
        chromosome = missense$Chromosome,
        start_pos = missense$Start_Position,
        ref_allele = missense$Reference_Allele,
        alt_allele = missense$Tumor_Seq_Allele2,
        sample_id = missense$Tumor_Sample_Barcode,
        variant_class = "Missense_Mutation",
        context = missense$context_11,
        stringsAsFactors = FALSE
    )
    
    # Convert sample ID
    result$sample_id <- sub("(TCGA-\\w{2}-\\w{4})-\\d+", "\\1", result$sample_id)
    
    # Clean
    result <- result[result$ref_allele %in% c("A","C","G","T"), ]
    result <- result[result$alt_allele %in% c("A","C","G","T"), ]
    result <- result[!is.na(result$context), ]
    
    return(result)
}


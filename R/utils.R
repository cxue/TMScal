# TMScal utility functions

#' Log message with timestamp
#' @param msg Message to log
#' @param level Log level
#' @export
log_message <- function(msg, level = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
}

#' Generate 192 pathway names
#' @return Character vector of 192 pathway names
#' @export
generate_pathway_names <- function() {
    bases <- c("A", "C", "G", "T")
    pathways <- character(192)
    idx <- 1
    for (left2 in bases) {
        for (left1 in bases) {
            for (ref in bases) {
                possible_alts <- setdiff(bases, ref)
                for (alt in possible_alts) {
                    pathways[idx] <- paste0(left2, left1, ref, "_", ref, "_", alt)
                    idx <- idx + 1
                }
            }
        }
    }
    return(pathways)
}

#' Extract pathway from mutation context
#' @param ref Reference allele
#' @param alt Alternate allele
#' @param context Trinucleotide context
#' @return Pathway string
#' @export
get_pathway <- function(ref, alt, context) {
    if (is.null(context) || is.na(context)) return(NA)
    context <- gsub("\\s+", "", context)
    
    if (nchar(context) == 3) {
        tri <- context
    } else if (nchar(context) == 11) {
        tri <- substr(context, 4, 6)
    } else if (nchar(context) >= 6) {
        tri <- substr(context, 4, 6)
    } else {
        return(NA)
    }
    
    if (nchar(ref) != 1 || !ref %in% c("A","C","G","T")) return(NA)
    if (nchar(alt) != 1 || !alt %in% c("A","C","G","T")) return(NA)
    
    return(paste0(tri, "_", ref, "_", alt))
}
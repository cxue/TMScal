#' Run complete TMScal pipeline
#' 
#' @param maf_file MAF file
#' @param clin_file Clinical file
#' @param sample_file Optional sample info
#' @param output_dir Output directory
#' @return List with all results
#' @export
run_all_pipeline <- function(maf_file, clin_file, sample_file = NULL,
                             output_dir = "./tmscal_output") {
    
    log_message("Starting TMScal pipeline...")
    
    # Prepare
    prepared <- prepare_data(maf_file, clin_file, sample_file,
                             file.path(output_dir, "prepared"))
    
    # Train
    model <- train(file.path(output_dir, "prepared", "maf_processed.txt"),
                   file.path(output_dir, "prepared", "clin_processed.txt"),
                   file.path(output_dir, "model"))
    
    # Predict
    predictions <- predict(file.path(output_dir, "prepared", "maf_processed.txt"),
                          model)
    
    return(list(prepared = prepared, model = model, predictions = predictions))
}
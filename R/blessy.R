#' Main blessy Function
#'
#' This function serves as the main entry point for the blessy package. It takes file paths for the transcript count data and domain mapping data,
#' processes them using a series of functions, and generates the final outputs: a DoCo count matrix and a final dictionary data frame.
#'
#' @param tx_count_file_path A string specifying the path to the transcript count file.
#' @param mapping_file_path A string specifying the path to the domain mapping file.
#'
#' @return A list containing two elements: 
#' \describe{
#'   \item{doco_count}{The aggregated DoCo count matrix.}
#'   \item{final_dict}{The final dictionary data frame with gene, transcript, and DoCo information.}
#' }
#' @examples
#' # Example usage:
#' results <- blessy("path/to/transcript_counts.txt", "path/to/domain_mapping.txt")
#' head(results$doco_count)
#' head(results$final_dict)
#' @export
blessy <- function(tx_count_file_path, mapping_file_path) {
    # Step 1: Load domain mapping data
    gr <- blessy.load_transcript_domains(mapping_file_path)
    
    # Step 2: Annotate domain data with DoCo information
    doco_se <- blessy.annotate_tx_with_doco(gr)
    
    # Step 3: Load transcript count data
    tx_count_se <- blessy.load_transcript_counts(tx_count_file_path)
    
    # Step 4: Annotate transcript counts with DoCo and gene information
    final_se <- blessy.annotate_doco_with_gene(tx_count_se, doco_se)
    
    # Step 5: Generate DoCo count matrix
    doco_count <- blessy.generate_doco_count(final_se)
    
    # Step 6: Generate final dictionary data frame
    final_dict <- blessy.generate_final_dict(final_se)
    
    # Return the results as a list
    return(list(doco_count = doco_count, final_dict = final_dict))
}

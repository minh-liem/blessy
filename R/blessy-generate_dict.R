#' Generate a Final Dictionary 
#'
#' This function creates a dictionary data frame containing gene, transcript, and DoCo information from a SummarizedExperiment object.
#'
#' @param final_se A SummarizedExperiment object that contains transcript counts and DoCo information.
#'
#' @return A data frame with three columns: Gene, Transcript, and DoCo.
#' @export
#'
#' @examples
#' # Assuming `final_se` is your SummarizedExperiment object
#' #final_dict <- blessy.generate_final_dict(final_se)
#' head(final_dict)
#' @export
blessy.generate_final_dict <- function(final_se) {
    # Extract the relevant data from the SummarizedExperiment object
    gene_ids <- rowData(final_se)$GENEID
    transcript_ids <- rownames(final_se)
    doco_vector <- rowData(final_se)$FinalDoCo
    
    # Create the final dictionary data frame
    final_dict <- data.frame(Gene = gene_ids, Transcript = transcript_ids, DoCo = doco_vector)
    
    # Return the final dictionary data frame
    return(final_dict)
}

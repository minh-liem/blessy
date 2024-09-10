#' Generate the DoCo Count Matrix
#'
#' This function aggregates the read counts of transcripts belonging to the same DoCo from a SummarizedExperiment object.
#'
#' @param final_se A SummarizedExperiment object that contains transcript counts and DoCo information.
#'
#' @return A data frame representing the DoCo count matrix, where each row corresponds to a unique DoCo with its aggregated read counts.
#' @export
#'
#' @examples
#' # Assuming `final_se` is your SummarizedExperiment object
#' #doco_count_matrix <- blessy.generate_doco_count(final_se)
#' head(doco_count_matrix)
#' @export
blessy.generate_doco_count <- function(final_se) {
    # Extract the count matrix and rowData from the SummarizedExperiment object
    count_matrix <- assays(final_se)$counts
    doco_vector <- rowData(final_se)$FinalDoCo
    
    # Aggregate the counts by DoCo
    doco_count <- aggregate(count_matrix, by = list(doco = doco_vector), FUN = sum)
    
    # Return the aggregated DoCo count matrix
    return(doco_count)
}

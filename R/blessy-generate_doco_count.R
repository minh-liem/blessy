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
    doco_vector <- rowData(final_se)$DoCo
    
    # Ensure that the length of doco_vector matches the number of rows in count_matrix
    if (nrow(count_matrix) != length(doco_vector)) {
        stop("Mismatch between the number of rows in the count matrix and the length of the DoCo vector.")
    }
    
    # Convert count_matrix to a data frame for easier manipulation if it's not already a data frame
    count_matrix_df <- as.data.frame(count_matrix)
    
    # Create a data frame combining the DoCo vector and the count matrix
    combined_df <- cbind(doco = doco_vector, count_matrix_df)
    
    # Aggregate the counts by DoCo
    doco_count <- aggregate(. ~ doco, data = combined_df, FUN = sum)
    
    # Return the aggregated DoCo count matrix
    return(doco_count)
}


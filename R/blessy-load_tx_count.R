#' Load Transcript Counts into a SummarizedExperiment Object
#'
#' This function reads a transcript count file and loads the data into a `SummarizedExperiment` object.
#'
#' @param file_path A string specifying the path to the transcript count file (in `.txt` format).
#'
#' @return A `SummarizedExperiment` object containing the transcript counts with associated gene IDs.
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' # Example usage:
#' # count_se <- blessy.load_transcript_counts("./data/tx_count.txt")
#' @export
blessy.load_transcript_counts <- function(file_path) {
    # Read the .txt file
    count_data <- read.table(file_path, header = TRUE, sep = "\t", row.names = 2)
    
    # Separate counts from gene metadata
    gene_ids <- count_data$GENEID  # Store GENEID separately
    count_matrix <- count_data[, -c(1, 2)]  # Remove the GENEID column from the counts
    
    # Create a SummarizedExperiment object
    se <- SummarizedExperiment(
        assays = list(counts = as.matrix(count_matrix)),
        rowData = DataFrame(GENEID = gene_ids)
    )
    
    return(se)
}
#' Load Transcript Counts into a SummarizedExperiment Object
#'
#' This function reads a transcript count file and loads the data into a `SummarizedExperiment` object.
#'
#' @param file_path A string specifying the path to the transcript count file (in `.txt` format).
#'
#' @return A `SummarizedExperiment` object containing the transcript counts with associated transcript IDs and gene IDs.
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @import S4Vectors
#' @examples
#' # Example usage:
#' count_se <- blessy.load_transcript_counts("./data/tx_count.txt")
#' @export
blessy.load_transcript_counts <- function(file_path) {
    # Read the .txt file
    count_data <- read.table(file_path, header = TRUE, sep = "\t")
    
    # Extract TXNAME (transcript ID) and GENEID
    tx_names <- count_data$TXNAME  # Store TXNAME separately
    gene_ids <- count_data$GENEID  # Store GENEID separately
    
    # Remove the TXNAME and GENEID columns to get the count matrix
    count_matrix <- count_data[, -c(1, 2)]  # Assuming TXNAME is column 1 and GENEID is column 2
    
    # Create a SummarizedExperiment object with both TXNAME and GENEID
    se <- SummarizedExperiment(
        assays = list(counts = as.matrix(count_matrix)),
        rowData = DataFrame(TXNAME = tx_names, GENEID = gene_ids)
    )
    
    return(se)
}

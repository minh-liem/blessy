#' Annotate transcripts with domain combinations (DoCo) using the mapping GRanges object
#'
#' @param gr A GRanges object containing transcript IDs and domains with their genomic coordinates.
#' @return A SummarizedExperiment object containing the DoCo for each transcript.
#' @import dplyr
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @examples
#' doco_se <- blessy.annotate_tx_with_doco(gr)
#' @export
blessy.annotate_tx_with_doco <- function(gr) {
    # Convert GRanges object to a data frame for easier manipulation
    gr_df <- as.data.frame(gr)
    
    # Order domains based on strand and coordinates
    ordered_df <- gr_df %>%
        group_by(mcols.TranscriptID, strand) %>%
        arrange(
            mcols.TranscriptID, 
            # Fix to handle strand sorting correctly
            if_else(strand == "-", -end, start)
        ) %>%
        summarise(
            doco = paste(
                paste(mcols.Domain, "::", seqnames, ":", start, "-", end, "(", strand, ")", sep = ""),
                collapse = ", "
            ),
            .groups = 'drop'
        )
    
    # Create a new SummarizedExperiment object to hold the DoCo information
    se <- SummarizedExperiment(
        assays = list(doco = matrix(ordered_df$doco, ncol = 1)),
        rowData = DataFrame(TranscriptID = ordered_df$mcols.TranscriptID)
    )
    
    return(se)
}


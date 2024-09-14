#' Annotate Transcript Counts with Domain Combinations (DoCo) and Gene Information
#'
#' This function takes transcript count data and domain combination (DoCo) data,
#' and adds the domain combination information to each transcript, along with the associated gene.
#'
#' @param count_se A `SummarizedExperiment` object containing transcript counts with associated transcript IDs and gene IDs.
#' @param doco_se A `SummarizedExperiment` object containing domain combination information (DoCo) for each transcript.
#'
#' @return A `SummarizedExperiment` object with annotated DoCo information.
#' @importFrom dplyr left_join
#' @examples
#' # Example usage:
#' final_se <- blessy.annotate_doco_with_gene(count_se, doco_se)
#' @export
blessy.annotate_doco_with_gene <- function(count_se, doco_se) {
    
    # Extract transcript IDs without version numbers
    transcript_ids <- sub("\\.\\d+$", "", rowData(count_se)$TXNAME)  # Extract TXNAME without versions
    
    # Create a data frame to map transcript IDs without version numbers to gene IDs
    count_df <- data.frame(
        TranscriptID = transcript_ids,
        GeneID = rowData(count_se)$GENEID,
        stringsAsFactors = FALSE
    )
    
    # Extract the DoCo information from doco_se
    doco_df <- data.frame(
        TranscriptID = rowData(doco_se)$TranscriptID,  # Extract TranscriptID from doco_se
        DoCo = assays(doco_se)$doco,                   # Extract DoCo from assays
        stringsAsFactors = FALSE
    )
    
    # Remove version numbers from TranscriptID in doco_df
    doco_df$TranscriptID <- sub("\\.\\d+$", "", doco_df$TranscriptID)
    
    # Merge DoCo with the count data based on the TranscriptID
    annotated_df <- dplyr::left_join(count_df, doco_df, by = "TranscriptID")
    
    # Handle missing DoCo information: for transcripts with no DoCo, set the DoCo to ";;GeneID"
    annotated_df$DoCo <- ifelse(
        is.na(annotated_df$DoCo),
        paste(";;", annotated_df$GeneID, sep = ""),
        paste(annotated_df$DoCo, annotated_df$GeneID, sep = ";;")
    )
    
    # Add the DoCo information to the count_se object
    rowData(count_se)$DoCo <- annotated_df$DoCo
    
    return(count_se)
}

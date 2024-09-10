#' Annotate the DoCo groups with gene ID associated with each transcript
#'
#' @param count_se A SummarizedExperiment object containing transcript counts with gene IDs.
#' @param doco_se A SummarizedExperiment object containing DoCo information based on transcript IDs.
#' 
#' @return A SummarizedExperiment object with transcript count with DoCo fully annotated (with gene IDs)
#' @import dplyr
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' # Example usage:
#' # final_se <- annotate_tx_with_doco(count_se, doco_se)
#' @export
blessy.annotate_doco_with_gene <- function(count_se, doco_se) {

    # Extract transcript IDs without version numbers
    transcript_ids <- sub("\\.\\d+$", "", rownames(count_se))
    
    # Create a data frame to map transcript IDs without version numbers to gene IDs
    count_df <- data.frame(
        TranscriptID = transcript_ids,
        GeneID = rowData(count_se)$GENEID,
        stringsAsFactors = FALSE
    )
    
    # Extract the DoCo information from doco_se
    doco_df <- data.frame(
        TranscriptID = rowData(doco_se)$TranscriptID,
        DoCo = assays(doco_se)$doco,
        stringsAsFactors = FALSE
    )
    
    # Remove version numbers from TranscriptID in doco_df
    doco_df$TranscriptID <- sub("\\.\\d+$", "", doco_df$TranscriptID)
    
    # Merge DoCo with the count data based on the TranscriptID
    annotated_df <- left_join(count_df, doco_df, by = "TranscriptID")
    
    # Handle missing DoCo information: for transcripts with no DoCo, set the FinalDoCo to ";;GeneID"
    annotated_df$FinalDoCo <- ifelse(
        is.na(annotated_df$DoCo),
        paste(";;", annotated_df$GeneID, sep = ""),
        paste(annotated_df$DoCo, annotated_df$GeneID, sep = ";;")
    )
    
    # Add the FinalDoCo information to the count_se object
    rowData(count_se)$FinalDoCo <- annotated_df$FinalDoCo
    
    return(count_se)
}




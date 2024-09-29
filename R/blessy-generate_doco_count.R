#' Generate the DoCo Count Matrix using dict_se, count_se, and GTF file for missing genes
#'
#' This function aggregates the read counts of transcripts belonging to the same DoCo, using
#' the DoCo information from `dict_se`, the transcript counts from `count_se`, and a GTF file for missing gene annotations.
#'
#' @param dict_se A SummarizedExperiment object containing the transcript-to-DoCo mapping.
#' @param count_se A SummarizedExperiment object containing transcript counts.
#' @param gtf_file_path A string specifying the path to the gene annotation GTF file.
#'
#' @return A list containing two elements:
#' \describe{
#'   \item{doco_count}{A data frame representing the DoCo count matrix, where each row corresponds to a unique DoCo with its aggregated read counts.}
#'   \item{final_dict_df}{A data frame representing the final dictionary with `TranscriptID`, `DoCo`, and `GeneID` columns.}
#' }
#' @importFrom SummarizedExperiment assays rowData
#' @importFrom dplyr left_join
#' @importFrom rtracklayer import
#' @importFrom utils write.csv
#' @examples
#' # Assuming `dict_se` and `count_se` are your SummarizedExperiment objects
#' results <- blessy.generate_doco_count(dict_se, count_se, "../../data/Homo_sapiens.GRCh38.111.gtf")
#' head(results$doco_count)
#' head(results$final_dict_df)
#' @export
blessy.generate_doco_count <- function(dict_se, count_se, gtf_file_path) {
    # Step 1: Extract DoCo and Gene information from dict_se, making sure GeneID is included
    doco_df <- data.frame(
        TranscriptID = rowData(dict_se)$TranscriptID,  # Extract TranscriptID from dict_se
        DoCo = assays(dict_se)$DoCo,                   # Extract DoCo from dict_se
        GeneID = rowData(dict_se)$GeneID               # Extract GeneID from rowData
    )
    
    # Step 2: Extract transcript counts from count_se using rowData to get the TranscriptID
    transcript_counts_df <- data.frame(
        TranscriptID = rowData(count_se)$TXNAME,        # Use TXNAME from rowData of count_se
        as.data.frame(assays(count_se)$counts)          # Extract count matrix from count_se
    )
    
    # Step 3: Perform a left join on TranscriptID to combine the DoCo information with transcript counts
    combined_df <- dplyr::left_join(transcript_counts_df, doco_df, by = "TranscriptID")
    
    # Step 4: Identify transcripts that are missing DoCo information
    missing_doco <- is.na(combined_df$DoCo)
    
    # Step 5: Import the GTF file to map genes for the missing transcripts
    gtf_data <- rtracklayer::import(gtf_file_path, format = "gtf")  # Use the rtracklayer namespace directly
    gtf_transcripts <- gtf_data[gtf_data$type == "transcript"]
    
    # Step 6: Create a data frame for gene-to-transcript mapping from the GTF file
    gtf_df <- data.frame(
        TranscriptID = sub("\\.\\d+$", "", mcols(gtf_transcripts)$transcript_id),  # Clean transcript ID
        GeneID = mcols(gtf_transcripts)$gene_id
    )
    
    # Step 7: Map the GeneID from the GTF file to the missing transcripts
    missing_transcripts <- combined_df$TranscriptID[missing_doco]
    gene_map <- gtf_df[gtf_df$TranscriptID %in% missing_transcripts, ]
    
    # Step 8: Assign the missing transcripts to DoCo as ";;GeneID", only if a matching gene is found, and update combined_df
    for (i in seq_along(missing_transcripts)) {
        tx_id <- missing_transcripts[i]
        gene_id <- gene_map$GeneID[gene_map$TranscriptID == tx_id]
        
        if (length(gene_id) > 0) {
            combined_df$DoCo[combined_df$TranscriptID == tx_id] <- paste0(";;", gene_id)
            combined_df$GeneID[combined_df$TranscriptID == tx_id] <- gene_id  # Update GeneID in combined_df
        } else {
            combined_df <- combined_df[combined_df$TranscriptID != tx_id, ]  # Remove transcripts with no matching gene
        }
    }
    
    # Step 9: Create the final dictionary data frame with TranscriptID, DoCo, and GeneID columns
    final_dict_df <- combined_df[, c("TranscriptID", "DoCo", "GeneID")]
    
    # Step 10: Save the final dictionary to "dict.csv"
    write.csv(final_dict_df, file = "./out/dict.csv", row.names = FALSE)
    print("Final dictionary saved to ./out/dict.csv.")
    
    # Step 11: Exclude non-numeric columns before aggregation
    combined_df <- combined_df[, !colnames(combined_df) %in% c("TranscriptID", "GeneID")]
    
    # Step 12: Aggregate the counts by DoCo (only numeric columns remain)
    doco_count <- aggregate(. ~ DoCo, data = combined_df, FUN = sum)
    
    # Step 13: Return the aggregated DoCo count matrix and final dictionary data frame
    return(doco_count)
}

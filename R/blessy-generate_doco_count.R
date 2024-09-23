#' Generate the DoCo Count Matrix using dict_se, count_se, and GTF file for missing genes
#'
#' This function aggregates the read counts of transcripts belonging to the same DoCo, using
#' the DoCo information from `dict_se`, the transcript counts from `count_se`, and a GTF file for missing gene annotations.
#'
#' @param dict_se A SummarizedExperiment object containing the transcript-to-DoCo mapping.
#' @param count_se A SummarizedExperiment object containing transcript counts.
#' @param gtf_file_path A string specifying the path to the gene annotation GTF file.
#'
#' @return A data frame representing the DoCo count matrix, where each row corresponds to a unique DoCo with its aggregated read counts.
#' @importFrom SummarizedExperiment assays rowData
#' @importFrom rtracklayer import
#' @importFrom dplyr left_join
#' @importFrom utils write.csv
#' @examples
#' # Assuming `dict_se` and `count_se` are your SummarizedExperiment objects
#' doco_count_matrix <- blessy.generate_doco_count(dict_se, count_se, "../../data/Homo_sapiens.GRCh38.111.gtf")
#' head(doco_count_matrix)
#' @export
blessy.generate_doco_count <- function(dict_se, count_se, gtf_file_path) {
    # Step 1: Extract DoCo information from dict_se
    doco_df <- data.frame(
        TranscriptID = rowData(dict_se)$TranscriptID,  # Extract TranscriptID from dict_se
        DoCo = assays(dict_se)$DoCo                    # Extract DoCo from dict_se
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
    gtf_data <- import(gtf_file_path, format = "gtf")
    gtf_transcripts <- gtf_data[gtf_data$type == "transcript"]
    
    # Step 6: Create a data frame for gene-to-transcript mapping from the GTF file
    gtf_df <- data.frame(
        TranscriptID = sub("\\.\\d+$", "", mcols(gtf_transcripts)$transcript_id),  # Clean transcript ID
        GeneID = mcols(gtf_transcripts)$gene_id
    )
    
    # Step 7: Map the GeneID from the GTF file to the missing transcripts
    missing_transcripts <- combined_df$TranscriptID[missing_doco]
    gene_map <- gtf_df[gtf_df$TranscriptID %in% missing_transcripts, ]
    
    # Step 8: Assign the missing transcripts to DoCo as ";;GeneID", only if a matching gene is found
    for (i in seq_along(missing_transcripts)) {
        tx_id <- missing_transcripts[i]
        gene_id <- gene_map$GeneID[gene_map$TranscriptID == tx_id]
        
        if (length(gene_id) > 0) {
            combined_df$DoCo[combined_df$TranscriptID == tx_id] <- paste0(";;", gene_id)
        } else {
            combined_df <- combined_df[combined_df$TranscriptID != tx_id, ]  # Remove transcripts with no matching gene
        }
    }
    
    # Step 9: Aggregate the counts by DoCo
    doco_count <- aggregate(. ~ DoCo, data = combined_df[, -1], FUN = sum)  # Remove TranscriptID column
    
    # Step 10: Return the aggregated DoCo count matrix
    return(doco_count)
}


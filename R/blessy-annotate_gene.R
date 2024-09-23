#' Annotate transcripts with Gene information using the GTF file
#'
#' This function maps transcripts in the DoCo SummarizedExperiment to genes using a GTF annotation file,
#' adds GeneID information to the domain combination (DoCo), and saves the result as a CSV file.
#'
#' @param doco_se A SummarizedExperiment object containing DoCo (domain combinations) for each transcript.
#' @param gtf_file_path A string specifying the path to the gene annotation GTF file.
#'
#' @return A SummarizedExperiment object containing the transcript-to-gene mapping with updated DoCo information. 
#' A CSV file "dict.csv" is saved in the `./out` directory.
#' 
#' @importFrom rtracklayer import
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr left_join
#' @importFrom utils write.csv
#' @examples
#' dict_se <- blessy.annotate_gene(doco_se, "../../data/Homo_sapiens.GRCh38.111.gtf")
#'
#' @export
blessy.annotate_gene <- function(doco_se, gtf_file_path) {
    # Step 1: Remove version numbers from TranscriptIDs in doco_se
    transcript_ids <- rowData(doco_se)$TranscriptID
    transcript_ids_clean <- sub("\\.\\d+$", "", transcript_ids)
    
    # Step 2: Import the GTF file and extract GeneID and TranscriptID
    gtf_data <- import(gtf_file_path, format = "gtf")
    
    # Filter GTF data for transcripts and genes only
    gtf_transcripts <- gtf_data[gtf_data$type == "transcript"]
    
    # Extract transcript and gene information from GTF
    gtf_df <- data.frame(
        TranscriptID = sub("\\.\\d+$", "", mcols(gtf_transcripts)$transcript_id),  # Remove transcript version from GTF
        GeneID = mcols(gtf_transcripts)$gene_id
    )
    
    # Step 3: Map GeneID to each transcript in doco_se by TranscriptID
    doco_df <- data.frame(
        TranscriptID = transcript_ids_clean,
        DoCo = assays(doco_se)$doco
    )
    
    # Perform a left join to map GeneIDs to transcripts
    annotated_df <- left_join(doco_df, gtf_df, by = "TranscriptID")
    
    # Step 4: For each transcript, append ";;GeneID" to the DoCo string
    annotated_df$DoCo <- paste(annotated_df$DoCo, annotated_df$GeneID, sep = ";;")
    
    # Step 5: Create a new SummarizedExperiment object with the updated DoCo information
    row_data <- DataFrame(TranscriptID = annotated_df$TranscriptID, GeneID = annotated_df$GeneID)
    assays <- list(DoCo = matrix(annotated_df$DoCo, ncol = 1))
    dict_se <- SummarizedExperiment(assays = assays, rowData = row_data)
    
    # Step 6: Write the result to dict.csv in the ./out directory
    if (!dir.exists("./out")) {
        dir.create("./out")
    }
    
    output_file <- "./out/dict.csv"
    write.csv(annotated_df, file = output_file, row.names = FALSE)
    
    # Print confirmation message
    print(paste("DoCo dictionary saved to:", output_file))
    
    # Return the dict_se object
    return(dict_se)
}

#' Main blessy Function
#'
#' This function serves as the main entry point for the blessy package. It takes file paths for the transcript count data, 
#' generates domain and transcript mappings, processes them using a series of functions, and generates the final outputs: 
#' a DoCo count matrix and a final dictionary data frame.
#'
#' @param tx_count_file_path A string specifying the path to the transcript count file.
#' @param domain_track A string specifying the domain track name (e.g., "unipDomain").
#' @param transcript_track A string specifying the transcript track name (e.g., "wgEncodeGencodeBasicV44").
#' @param genome A string specifying the genome version (e.g., "hg38").
#'
#' @return A list containing two elements: 
#' \describe{
#'   \item{doco_count}{The aggregated DoCo count matrix.}
#'   \item{final_dict}{The final dictionary data frame with gene, transcript, and DoCo information.}
#' }
#' @examples
#' # Example usage:
#' results <- blessy("path/to/transcript_counts.txt", "unipDomain", "wgEncodeGencodeBasicV44", "hg38")
#' head(results$doco_count)
#' head(results$final_dict)
#' @export
blessy <- function(tx_count_file_path, domain_track, transcript_track, genome) {
    # Step 1: Get domain track and save as a BED file
    blessy.get_domain_track(genome, domain_track)
    
    # Step 2: Get transcript track and save as a BED file
    blessy.get_transcript_track(genome, transcript_track)
    
    # Step 3: Map transcripts to domains and get the SummarizedExperiment object
    se <- blessy.transcript_domain_mapping(
        paste0("./out/", domain_track, "_", genome, "_domain.bed"),
        paste0("./out/", transcript_track, "_", genome, "_transcript.bed")
    )
    
    # Step 4: Load domain mapping data using the SummarizedExperiment object
    gr <- blessy.load_transcript_domain_mapping(se)
    
    # Step 5: Annotate domain data with DoCo information
    doco_se <- blessy.annotate_tx_with_doco(gr)
    
    # Step 6: Load transcript count data
    tx_count_se <- blessy.load_transcript_counts(tx_count_file_path)
    
    # Step 7: Annotate transcript counts with DoCo and gene information
    final_se <- blessy.annotate_doco_with_gene(tx_count_se, doco_se)
    
    # Step 8: Generate DoCo count matrix
    doco_count <- blessy.generate_doco_count(final_se)
    
    # Step 9: Generate final dictionary data frame
    final_dict <- blessy.generate_final_dict(final_se)
    
    # Return the results as a list
    return(list(doco_count = doco_count, final_dict = final_dict))
}
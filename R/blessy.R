#' Main blessy Function
#'
#' This function serves as the main entry point for the blessy package. It takes file paths for the transcript count data, 
#' generates domain and transcript mappings, processes them using a series of functions, and generates the final output: 
#' a DoCo count matrix.
#'
#' @param tx_count_file_path A string specifying the path to the transcript count file.
#' @param domain_track A string specifying the domain track name (e.g., "unipDomain").
#' @param transcript_track A string specifying the transcript track name (e.g., "wgEncodeGencodeBasicV44").
#' @param genome A string specifying the genome version (e.g., "hg38").
#' @param gtf_file_path A string specifying the path to the gene annotation GTF file.
#'
#' @return A data frame representing the DoCo count matrix.
#' @examples
#' # Example usage:
#' results <- blessy("path/to/transcript_counts.txt", "unipDomain", "wgEncodeGencodeBasicV44", "hg38", "path/to/gtf_file.gtf")
#' head(results)
#' @export
blessy <- function(tx_count_file_path, domain_track, transcript_track, genome, gtf_file_path) {
    
    # Step 1: Get domain track and store in domain_df
    domain_df <- blessy.get_domain_track(genome, domain_track)
    
    # Step 2: Get transcript track and store in transcript_df
    transcript_df <- blessy.get_transcript_track(genome, transcript_track)
    
    # Step 3: Map transcripts to domains and get the SummarizedExperiment object (using data frames, not file paths)
    se <- blessy.transcript_domain_mapping(domain_df, transcript_df)
    
    # Step 4: Load domain mapping data using the SummarizedExperiment object
    gr <- blessy.load_transcript_domain_mapping(se)
    
    # Step 5: Annotate domain data with DoCo information
    doco_se <- blessy.annotate_tx_with_doco(gr)
    
    # Step 6: Annotate transcripts with gene information using the GTF file
    dict_se <- blessy.annotate_gene(doco_se, gtf_file_path)
    
    # Step 7: Load transcript count data
    count_se <- blessy.load_transcript_counts(tx_count_file_path)
    
    # Step 8: Generate DoCo count matrix using dict_se and count_se
    doco_count <- blessy.generate_doco_count(dict_se, count_se, gtf_file_path)
    
    # Return the DoCo count matrix
    return(doco_count)
}

#' Map Transcripts to Domains and Save Results
#'
#' This function maps transcripts to domains by finding overlaps between transcript and domain BED files,
#' and saves the result as a CSV file. It also returns a `SummarizedExperiment` object.
#'
#' @param domain_bed_file A string specifying the path to the domain BED file.
#' @param transcript_bed_file A string specifying the path to the transcript BED file.
#'
#' @return A `SummarizedExperiment` object containing the mapping of transcripts to domains. A CSV file
#' with the transcript-domain mapping is saved in the `./out` directory as `transcript_domain_mapping.csv`.
#'
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom utils write.csv
#'
#' @examples
#' # Example usage:
#' blessy.transcript_domain_mapping("./out/unipDomain_hg38_domain.bed", "./out/wgEncodeGencodeBasicV44_hg38_transcript.bed")
#'
#' @export
blessy.transcript_domain_mapping <- function(domain_bed_file, transcript_bed_file) {
    # Step 1: Import the BED files as GRanges objects
    domain_bed <- import(domain_bed_file, format = "BED")
    transcript_bed <- import(transcript_bed_file, format = "BED")
    
    # Step 2: Find common sequence levels between the two BED files
    common_seqlevels <- intersect(seqlevels(transcript_bed), seqlevels(domain_bed))
    
    # Step 3: Restrict both GRanges objects to only the common sequence levels
    transcript_bed <- keepSeqlevels(transcript_bed, common_seqlevels, pruning.mode = "coarse")
    domain_bed <- keepSeqlevels(domain_bed, common_seqlevels, pruning.mode = "coarse")
    
    # Step 4: Find overlaps between transcripts and domains
    overlaps <- findOverlaps(transcript_bed, domain_bed)
    
    # Step 5: Extract the overlapping transcripts and domains
    transcript_hits <- as.data.frame(transcript_bed[queryHits(overlaps)])
    domain_hits <- as.data.frame(domain_bed[subjectHits(overlaps)])
    
    # Step 6: Format the domain information as required: domain_id:chr:start-end(strand)
    domain_info <- paste0(
        domain_hits$name, ":", 
        domain_hits$seqnames, ":", 
        domain_hits$start, "-", 
        domain_hits$end, "(", 
        domain_hits$strand, ")"
    )
    
    # Step 7: Create the final output data frame with transcript_id and domain_info
    output_df <- data.frame(
        transcript_id = transcript_hits$name,
        domain_info = domain_info
    )
    
    # Step 8: Create output directory if it doesn't exist
    if (!dir.exists("./out")) {
        dir.create("./out")
    }
    
    # Step 9: Save the result as "transcript_domain_mapping.csv" in the ./out directory
    output_file <- "./out/transcript_domain_mapping.csv"
    write.csv(output_df, file = output_file, row.names = FALSE)
    
    # Print confirmation message
    print(paste("Transcript-domain mapping saved to:", output_file))
    
    # Step 10: Convert the result into a SummarizedExperiment object
    rowData <- DataFrame(transcript_id = output_df$transcript_id)
    assays <- list(domain_info = as.matrix(output_df$domain_info))
    se <- SummarizedExperiment(assays = assays, rowData = rowData)
    
    # Return the SummarizedExperiment object
    return(se)
}

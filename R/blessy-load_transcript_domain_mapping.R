#' Load Mapping from SummarizedExperiment into a GRanges Object
#'
#' This function reads a `SummarizedExperiment` object containing transcript IDs and protein domain information
#' with genomic coordinates, and loads the data into a `GRanges` object from the `GenomicRanges` package using regex.
#'
#' @param se A `SummarizedExperiment` object containing transcript and domain information.
#'
#' @return A GRanges object containing the genomic ranges of the protein domains with transcript IDs as metadata.
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame
#' @examples
#' # Example usage:
#' # Load and test the GRanges object from a SummarizedExperiment
#' gr <- blessy.load_transcript_domain_mapping(se)
#' @export
blessy.load_transcript_domain_mapping <- function(se) {
    # Extract domain and transcript information from the SummarizedExperiment
    transcript_ids <- rowData(se)$transcript_id
    domain_info <- as.vector(assays(se)$domain_info)
    
    # Parse the domain info using the corrected regex
    matches <- regmatches(domain_info, regexec("^(.*):chr([^:]+):(\\d+)-(\\d+)\\((\\+|-)\\)$", domain_info))
    
    # Extract relevant components
    domain_names <- sapply(matches, function(x) x[2])  # Domain name
    chromosomes <- sapply(matches, function(x) x[3])   # Chromosome
    starts <- sapply(matches, function(x) x[4])        # Start position
    ends <- sapply(matches, function(x) x[5])          # End position
    strands <- sapply(matches, function(x) x[6])       # Strand
    
    # Create a GRanges object
    gr <- GRanges(
        seqnames = chromosomes,
        ranges = IRanges(start = as.numeric(starts), end = as.numeric(ends)),
        strand = strands,
        mcols = DataFrame(TranscriptID = transcript_ids, Domain = domain_names)
    )
    
    # Return the GRanges object
    return(gr)
}




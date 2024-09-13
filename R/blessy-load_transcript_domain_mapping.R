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
#' gr <- blessy.load_transcript_domain(se)
#' @export
blessy.load_transcript_domain_mapping <- function(se) {
    # Extract domain and transcript information from the SummarizedExperiment
    transcript_ids <- rowData(se)$transcript_id
    domain_info <- as.vector(assays(se)$domain_info)
    
    # Parse the domain info using a more flexible regex
    matches <- regmatches(domain_info, regexec("(.*)::(chr[^:]+):(\\d+)-(\\d+)\\((\\+|-)\\)", domain_info))
    
    # Extract relevant components
    domain_names <- sapply(matches, function(x) if (length(x) == 6) x[2] else NA)  # Domain name
    chromosomes <- sapply(matches, function(x) if (length(x) == 6) x[3] else NA)   # Chromosome
    starts <- sapply(matches, function(x) if (length(x) == 6) x[4] else NA)        # Start position
    ends <- sapply(matches, function(x) if (length(x) == 6) x[5] else NA)          # End position
    strands <- sapply(matches, function(x) if (length(x) == 6) x[6] else NA)       # Strand
    
    # Identify rows with NAs (indicating unmatched data)
    na_rows <- is.na(starts) | is.na(ends)
    
    # Print or handle rows with unmatched data
    if (any(na_rows)) {
        warning("Some rows have unmatched data and will be excluded:")
        print(data.frame(transcript_ids[na_rows], domain_info[na_rows]))
    }
    
    # Remove rows with NAs
    valid_rows <- !na_rows
    
    # Create a GRanges object
    gr <- GRanges(
        seqnames = chromosomes[valid_rows],
        ranges = IRanges(start = as.numeric(starts[valid_rows]), end = as.numeric(ends[valid_rows])),
        strand = strands[valid_rows],
        mcols = DataFrame(TranscriptID = transcript_ids[valid_rows], Domain = domain_names[valid_rows])
    )
    
    return(gr)
}

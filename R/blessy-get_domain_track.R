#' Get Domain Track Data and Save as BED File
#'
#' This function fetches domain track data from the UCSC Genome Browser API and saves it as a BED file.
#'
#' @param genome A string specifying the genome version (e.g., "hg38").
#' @param track A string specifying the track name (e.g., "unipDomain").
#'
#' @return A BED file saved in the `./out` directory with the domain track data. The function prints a message confirming the save.
#'
#' @import httr
#' @import jsonlite
#' @importFrom dplyr select
#'
#' @examples
#' # Example usage:
#' domain_df <- blessy.get_domain_track("hg38", "unipDomain")
#'
#' @export
blessy.get_domain_track <- function(genome, track) {
    # Step 1: Define the UCSC Genome Browser API URL using the input parameters
    url <- paste0("https://api.genome.ucsc.edu/getData/track?genome=", genome, ";track=", track, ";maxItemsOutput=-1")
    
    # Step 2: Make the API request to fetch the track data
    response <- GET(url)
    
    # Step 3: Parse the JSON response
    content <- content(response, as = "text")
    data <- fromJSON(content)
    
    # Step 4: Extract the track information (e.g., unipDomain)
    track_data <- data[[track]]
    
    # Step 5: Convert the track data to a data frame
    track_df <- as.data.frame(track_data)
    
    # Step 6: Select the necessary columns for the BED file: chrom, chromStart, chromEnd, name, score, strand
    bed_df <- track_df %>%
        select(chrom, chromStart, chromEnd, name, score, strand)
    
    # Step 7: Create output directory if it doesn't exist
    if (!dir.exists("./out")) {
        dir.create("./out")
    }
    
    # Step 8: Save the data as a BED file in the ./out directory
    output_file <- paste0("./out/", track, "_", genome, "_domain.bed")
    write.table(bed_df, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Optional: Print confirmation message
    print(paste("Domain track data saved as:", output_file))
    
    return(bed_df)
}

#' Get Transcript Track Data and Save as BED File
#'
#' This function fetches transcript track data from the UCSC Genome Browser API and saves it as a BED file.
#' It also returns the transcript data frame (bed_df).
#'
#' @param genome A string specifying the genome version (e.g., "hg38").
#' @param track A string specifying the track name (e.g., "wgEncodeGencodeBasicV44").
#'
#' @return A BED file saved in the `./out` directory with the transcript track data, and the transcript data frame.
#'
#' @import httr
#' @import jsonlite
#' @importFrom dplyr select bind_rows
#'
#' @examples
#' tx_df <- blessy.get_transcript_track("hg38", "wgEncodeGencodeBasicV44")
#'
#' @export
blessy.get_transcript_track <- function(genome, track) {
    # Step 1: Define UCSC Genome Browser API URL using the input parameters
    url <- paste0("https://api.genome.ucsc.edu/getData/track?genome=", genome, ";track=", track, ";maxItemsOutput=-1")
    
    # Step 2: Make the API request to fetch the transcript data
    response <- GET(url)
    
    # Step 3: Parse the JSON response
    content <- content(response, as = "text", encoding = "UTF-8")
    data <- fromJSON(content)
    
    # Step 4: Initialize an empty list to store data frames for each chromosome
    all_transcripts <- list()
    
    # Step 5: Loop over each chromosome in the track data
    for (chrom in names(data[[track]])) {
        # Extract data for the current chromosome
        chrom_data <- data[[track]][[chrom]]
        
        # Check if chrom_data is a data frame and is not empty
        if (is.data.frame(chrom_data) && nrow(chrom_data) > 0) {
            # Convert it to a data frame (already a data frame, kept for consistency)
            chrom_df <- as.data.frame(chrom_data)
            
            # Add the chromosome name as a column
            chrom_df$chromosome <- chrom
            
            # Append to the list
            all_transcripts[[chrom]] <- chrom_df
        }
    }
    
    # Step 6: Combine all individual chromosome data frames into one
    transcripts_df <- bind_rows(all_transcripts)
    
    # Step 7: Select necessary columns for BED file (chromosome, chromStart, chromEnd, name, score, strand)
    bed_df <- transcripts_df %>%
        select(chromosome = chrom, chromStart = txStart, chromEnd = txEnd, name = name, score = score, strand)
    
    # Step 8: Ensure the score field is numeric and set to 0 if missing
    bed_df$score[is.na(bed_df$score)] <- 0
    
    # Step 9: Create output directory if it doesn't exist
    if (!dir.exists("./out")) {
        dir.create("./out")
    }
    
    # Step 10: Save the data as a BED file in the ./out directory
    output_file <- paste0("./out/", track, "_", genome, "_transcript.bed")
    write.table(bed_df, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Print confirmation message
    print(paste("Transcript track data saved as:", output_file))
    
    # Step 11: Return the transcript data frame
    return(bed_df)
}

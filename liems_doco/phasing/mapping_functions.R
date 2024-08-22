#' Extract genomic coordinates from domain string
#'
#' @param domain The domain mapped to each transcript in the mapping.txt file
#' @return A list containing chromosome, start, end, and strand information
extract_coordinates <- function(domain) {
  library(dplyr)
  library(tidyr)
  match <- regmatches(domain, regexec(".*::chr(\\d+|X|Y):(\\d+)-(\\d+)\\((.)\\)", domain))
  if (length(match[[1]]) == 5) {
    chr <- match[[1]][2]
    start <- as.numeric(match[[1]][3])
    end <- as.numeric(match[[1]][4])
    strand <- match[[1]][5]
    return(list(chr = chr, start = start, end = end, strand = strand))
  }
  return(list(chr = NA, start = NA, end = NA, strand = NA))
}




#' Process mapping data from a file
#'
#' @param filepath Path to the mapping data file
#' @import dplyr
#' @import tidyr
#' @return A data frame of processed mapping data
process_mapping_data <- function(filepath) {
  mapping <- read.table(filepath)
  mapping <- mapping %>%
    rowwise() %>%
    mutate(coordinates = list(extract_coordinates(V2))) %>%
    unnest_wider(coordinates)
  
  mapping_processed <- mapping %>%
    group_by(V1, strand) %>%
    arrange(V1, if_else(strand == "-", desc(end), start)) %>%
    summarise(V2 = paste(V2, collapse = ", "), .groups = 'drop')
  
  colnames(mapping_processed) <- c("transcript_id", "strand", "doco")
  return(mapping_processed)
}


# Example usage
mapping_processed <- process_mapping_data("./data/mapping.txt")
write.csv(mapping_processed, "./data/partial_dict.csv", row.names = FALSE)


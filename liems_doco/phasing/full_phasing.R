#' Load and Format Partial Dictionary
#'
#' @param filepath Path to the CSV file containing the partial dictionary
#' @return A data frame with formatted partial dictionary
load_and_format_dict <- function(filepath) {
  dict <- read.csv(filepath)
  dict$transcript_id <- sub("\\..*", "", dict$transcript_id) # Compatibility adjustment
  return(dict)
}


#' Load and Format Transcript Count Data
#'
#' @param filepath Path to the text file containing transcript counts
#' @return A data frame with loaded transcript counts
load_transcript_count <- function(filepath) {
  count <- read.table(filepath, header = TRUE)
  return(count)
}


#' Annotate Transcript Counts with DOCO Information
#'
#' @param count Data frame of counts
#' @param dict Data frame of dictionary entries
#' @return Annotated count data frame
annotate_doco <- function(count, dict) {
  count$doco <- ""
  
  for (transcript in count$TXNAME) {
    matches <- dict[dict$transcript_id == transcript, ]
    gene <- count$GENEID[count$TXNAME == transcript]
    count$doco[count$TXNAME == transcript] <- paste(matches$doco, gene, sep = ';;')
  }
  
  return(count)
}


#' Aggregate DOCO Information
#'
#' @param count Data frame with DOCO annotations
#' @return Aggregated DOCO count
aggregate_doco_count <- function(count) {
  doco_count <- aggregate(count[3:10], by = list(count$doco), FUN = sum)
  names(doco_count)[1] <- "doco"
  return(doco_count)
}


#' Generate Final Dictionary Containing Gene, DoCo and Tx
#'
#' @param phasing_count Data frame of phasing counts
#' @return Data frame with gene IDs added and reordered
generate_final_dict <- function(count) {
  final_dict <- data.frame(count$GENEID, count$TXNAME, count$doco)
  colnames(final_dict) <- c("Gene", "Transcript", "DoCo")
  return(final_dict)
}


# Example usage of the functions
dict <- load_and_format_dict("./data/partial_dict.csv")
count <- load_transcript_count("./data/tx_count.txt")
count <- annotate_doco(count, dict)
final_dict <- generate_final_dict(count)
doco_count <- aggregate_doco_count(count)


# Writing outputs
write.csv(count, "./out/all_feat_count.csv", row.names = FALSE)
write.csv(final_dict, "./out/final_dict.csv", row.names = FALSE)
write.csv(doco_count, "./out/doco_count.csv", row.names = FALSE)



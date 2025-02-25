% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blessy-get_domain_track.R
\name{blessy.get_domain_track}
\alias{blessy.get_domain_track}
\title{Get Domain Track Data and Save as BED File}
\usage{
blessy.get_domain_track(genome, track)
}
\arguments{
\item{genome}{A string specifying the genome version (e.g., "hg38").}

\item{track}{A string specifying the track name (e.g., "unipDomain").}
}
\value{
A BED file saved in the \code{./out} directory with the domain track data. The function prints a message confirming the save.
}
\description{
This function fetches domain track data from the UCSC Genome Browser API and saves it as a BED file.
}
\examples{
# Example usage:
domain_df <- blessy.get_domain_track("hg38", "unipDomain")

}

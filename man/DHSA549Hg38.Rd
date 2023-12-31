% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{DHSA549Hg38}
\alias{DHSA549Hg38}
\title{DNase hypersensitivity (DHS) peaks in A549 cell example data, lifted to hg38}
\description{
An example dataset containing narrowPeak file from ENCODE. Retrieve
record with \code{object[["AH22505"]]} on Annotation Hub. Construction
script is in 'inst/script/DHSA549Hg38.R'. Function returns a GRanges
object with metadata score, signal value, p/q value and peak.
}
\examples{

dhs <- DHSA549Hg38()
dhs

}

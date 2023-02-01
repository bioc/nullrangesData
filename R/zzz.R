#' ExperimentHub datasets for the nullranges package
#'
#' DNase hypersensitivity sites (DHS), CTCF binding sites, and
#' CTCF genomic interactions for demonstration of functions in the
#' nullranges package.
#'
#' @examples
#'
#' suppressPackageStartupMessages(library(GenomicRanges))
#' dhs <- DHSA549Hg38()
#' dhs
#'
#' @importFrom utils read.csv
#' @importFrom ExperimentHub createHubAccessors
#' @import GenomicRanges InteractionSet
#' @docType package
#' @name nullrangesData
NULL

.onLoad <- function(libname, pkgname) {
    fl <- system.file("extdata", "metadata.csv", package=pkgname)
    titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
    ExperimentHub::createHubAccessors(pkgname, titles)
}

#' DNase hypersensitivity (DHS) peaks in A549 cell example data, lifted to hg38
#'
#' An example dataset containing narrowPeak file from ENCODE. Retrieve
#' record with \code{object[["AH22505"]]} on Annotation Hub. Construction
#' script is in 'inst/script/DHSA549Hg38.R'. Function returns a GRanges
#' object with metadata score, signal value, p/q value and peak.
#'
#' @name DHSA549Hg38
#'
#' @examples
#'
#' dhs <- DHSA549Hg38()
#' dhs
#'
NULL

#' 10Kb bins from hg19 with GM12878 metadata annotation features
#'
#' 10Kb bins were tiled across hg19 and annotated with CTCF and DNase
#' site features from GM12878. Feature annotations for each bin
#' include 1) the number of CTCF sites, 2) the CTCF signal strength
#' (from peak calls), 3) the number of DNase sites, 4) the DNase
#' signal strength (from signal tracks), and finally 5) the
#' presence/absence of a loop to any other bin. Function returns
#' a GRanges object with covariate metadata
#'
#' @examples
#'
#' bins <- hg19_10kb_bins()
#' bins
#'
#' @name hg19_10kb_bins
NULL

#' CTCF-bound 10Kb paired genomic interactions
#'
#' 10Kb bins were tiled across hg19 then subset by those which
#' contained CTCF sites. All pairs of CTCF-bound 10Kb bins were
#' generated and annotated with feature overlaps from GM12878. Feature
#' annotations include 1) presence/absence of a loop between
#' bin-pairs, 2) the total CTCF signal from both bin-pairs, 3) the
#' number of CTCF sites from both bin-pairs, 4) the distance between
#' bin-pairs, and finally 4) whether a convergent set of CTCF sites
#' exists between bin-pairs. Function returns
#' a GInteractions object with covariate metadata
#'
#' @examples
#'
#' binPairs <- hg19_10kb_ctcfBoundBinPairs()
#' binPairs
#'
#' @name hg19_10kb_ctcfBoundBinPairs
NULL

#' Chromium Single Cell Multiome ATAC + Gene Expression assay
#'
#' An example dataset of Single cell Gene and ATAC-seq data for
#' Construction script is in 'inst/script/make-multiome-data.R'. Function returns a gene expression GRanges
#' object with pseudo bulk counts stored in metadata column "counts1" as NumericList.
#'
#' @examples
#'
#' data("sc_rna")
#'
#' @name sc_rna
NULL

#' Chromium Single Cell Multiome ATAC + Gene Expression assay
#'
#' An example dataset of Single cell Gene and ATAC-seq data for
#' Construction script is in 'inst/script/make-multiome-data.R'. Function returns a promoter GRanges
#' object with pseudo bulk counts stored in metadata column "counts2" as NumericList.
#'
#' @examples
#'
#' data("sc_promoter")
#'
#' @name sc_promoter
NULL

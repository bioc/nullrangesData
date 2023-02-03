library(excluderanges)
suppressMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(GenomicRanges))
library(nullranges)
ah <- AnnotationHub()

## Derive exclude regions from ENCODE
query_data <- query(ah, c("excluderanges", "hg38", "Exclusion regions"))
query_data
exclude_hg38_Kundaje <- query_data[["AH95917"]]

## Derive telomere,centromere from UCSC and rCGH
query_data2 <- query(ah, c("excluderanges", "UCSC", "Homo Sapiens", "hg38"))
query_data2
telomere <- query_data2[["AH95938"]]
telomere <- trim(telomere)

suppressPackageStartupMessages(library(rCGH))
# hg38 # data.frame
# Adjust chromosome names
hg38$chrom[hg38$chrom == 23] <- "X"
hg38$chrom[hg38$chrom == 24] <- "Y"
hg38$chrom <- paste0("chr", hg38$chrom)
# Make GRanges object
hg38.UCSC.centromere <- makeGRangesFromDataFrame(
  hg38,
  seqnames.field = "chrom",
  start.field = "centromerStart",
  end.field = "centromerEnd")
# Assign seqinfo data
seqlengths(hg38.UCSC.centromere) <- hg38$length
genome(hg38.UCSC.centromere) <- "hg38"
# Resulting object
hg38.UCSC.centromere

## Combining ENCODE with centromere and telomere
exclude_hg38_all <- reduce(c(exclude_hg38_Kundaje, telomere, hg38.UCSC.centromere))
exclude_hg38_all
summary(width(exclude_hg38_all))
## Remove small pieces for bootstrap speed up
exclude <- exclude_hg38_all %>%
  plyranges::filter(width(exclude_hg38_all) >= 500)  

### Based on gene density
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86
filt <- AnnotationFilterList(GeneIdFilter("ENSG", "startsWith"))
g <- genes(edb, filter = filt)

library(GenomeInfoDb)
g <- keepStandardChromosomes(g, pruning.mode = "coarse")
seqlevels(g, pruning.mode="coarse") <- setdiff(seqlevels(g), c("MT"))
# normally we would assign a new style, but for recent host issues
## seqlevelsStyle(g) <- "UCSC" 
seqlevels(g) <- paste0("chr", seqlevels(g))
genome(g) <- "hg38"

g <- sortSeqlevels(g)
g <- sort(g)

## Segmentation with tiling block length 2e6
L_s <- 2e6
seg_cbs <- segmentDensity(g, n = 3, L_s = L_s,
                          exclude = exclude, type = "cbs")

seg_hmm <- segmentDensity(g, n = 3, L_s = L_s,
                          exclude = exclude, type = "hmm")

save(seg_cbs, file = "data/seg_cbs.rda", compress = "xz")
save(seg_hmm, file = "data/seg_hmm.rda", compress = "xz")
save(exclude_hg38_all, file = "data/exclude_hg38_all.rda", compress = "xz")

plot <- lapply(c("ranges", "barplot", "boxplot"), function(x) plotSegment(seg_cbs,type=x,exclude=exclude))
plot[[1]]
plot[[2]]
plot[[3]]

meta <- data.frame(
  Title = c("DHSA549Hg38","hg19_10kb_bins","hg19_10kb_ctcfBoundBinPairs"),
  Description =
    c("DNase hypersensitive peaks in A549 cell example data, lifted to hg38",
      "10Kb bins from hg19 with GM12878 metadata annotation features",
      "CTCF-bound 10Kb paired genomic interactions"),
  BiocVersion = c("3.14","3.14","3.14"),
  Genome = c("hg38","hg19","hg19"),
  SourceType = c("BED","TXT","TXT"),
  SourceUrl =
    c("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/wgEncodeAwgDnaseUwdukeA549UniPk.narrowPeak.gz",
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525",
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525"),
  SourceVersion = c("v1","v1","v1"),
  Species = rep("Homo sapiens", 3),
  TaxonomyId = rep(9606, 3),
  Coordinate_1_based = rep(TRUE, 3),
  DataProvider = c("UCSC","Aiden Lab","Aiden Lab"),
  Maintainer = rep("Michael Love <michaelisaiahlove@gmail.com>", 3),
  RDataClass = c("GenomicRanges","GenomicRanges","InteractionSet"),
  DispatchClass = rep("Rda", 3),
  RDataPath = file.path("nullrangesData","v1",c("DHSA549Hg38.rda","hg19_10kb_bins.rda","hg19_10kb_ctcfBoundBinPairs.rda")),
  Tags = c("wgEncode:DnaseSeq:A549 cell","DnaseSeq:CTCF:GM12878 cell","HiC:CTCF:GM12878 cell"),
  Notes = c("","",""))

write.csv(meta, file="metadata.csv", row.names=FALSE)

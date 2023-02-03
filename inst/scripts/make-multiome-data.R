suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(plyranges))
# Load multi-modal Seurat object
seurat <- readRDS(url("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/seurat.rds"))
passIndex <- seurat@meta.data$pass_accQC==TRUE & seurat@meta.data$pass_rnaQC==TRUE
## use cells pass both assay
seurat <- seurat %>%
  .[,seurat@meta.data$pass_accQC==TRUE & seurat@meta.data$pass_rnaQC==TRUE]

# Loda feature metadata
feature_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/filtered_feature_bc_matrix/features.tsv.gz") %>%
  setnames(c("ens_id","gene","view","chr","start","end"))
feature_metadata.rna <- feature_metadata[view=="Gene Expression"]

feature_metadata.atac <- feature_metadata[view=="Peaks"] %>%
  .[,ens_id:=NULL] %>% setnames("gene","peak")
# head(feature_metadata.atac,n=3)

## classify peaks into distal and promoter, select promoter
foo <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/atac_peak_annotation.tsv") %>%
  .[,c("peak","peak_type")] %>%
  .[peak_type%in%c("distal", "promoter")]
feature_metadata.atac <- feature_metadata.atac %>%
  merge(foo,by="peak",all.x=TRUE)

## construct RNA Granges with counts as CompressedNumericList
feature_metadata.rna2 <-feature_metadata.rna[feature_metadata.rna$gene %in% (seurat@assays$RNA@counts %>% rownames())]
rna_name <- feature_metadata.rna2$gene[feature_metadata.rna2$chr%in% paste0("chr",1:22)]
rna_Granges<-feature_metadata.rna2[which(chr %in% paste0("chr",1:22))] %>% .[,c("chr","start","end","gene")] %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)

## construct promoter Granges with counts as compressedNumericList
promoter_name<-feature_metadata.atac[which(chr %in% paste0("chr",1:22) & peak_type=="promoter")]$peak
promoter_Granges <-feature_metadata.atac[which(chr %in% paste0("chr",1:22)& peak_type=="promoter")] %>% .[,c("chr","start","end","peak")] %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)

groups <- seurat@meta.data$celltype
levels(groups)<-setdiff(unique(groups),NA)
# Aggregate across cluster-sample groups
getPseudobulk <- function(mat.sparse, celltype) {
  mat.summary <- do.call(cbind, lapply(levels(celltype), function(ct) {
    cells <- which(celltype==ct)
    pseudobulk <- rowSums(mat.sparse[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}

## call function
rna.summary <- getPseudobulk(GetAssayData(seurat)[rna_name,], groups)
promoter.summary <- getPseudobulk(seurat@assays$ATAC[promoter_name,], groups)

rna.sd <-apply(rna.summary, 1, sd)
### remove rna that no changs across cell types
rna.summary <- rna.summary[-which(rna.sd==0),]
promoter.sd <-apply(promoter.summary, 1, sd)
### remove promoter that no changs across cell types
promoter.summary <- promoter.summary[-which(promoter.sd==0),]

library(edgeR)
rna.scaled <-edgeR::cpm(rna.summary,log = T)
promoter.scaled <-edgeR::cpm(promoter.summary,log = T)

## split sparse count matrix into NumericList
sc_rna <- rna_Granges[-which(rna.sd==0)] %>%
  mutate(counts1 = NumericList(asplit(rna.scaled, 1))) %>% sort()
sc_promoter <- promoter_Granges[-which(promoter.sd==0)] %>%
  mutate(counts2 = NumericList(asplit(promoter.scaled, 1))) %>% sort()

## use join_overlap_inner to separate read counts all=0 from non-overlaping
sc_rna$counts1 <- NumericList(lapply(sc_rna$counts1, function(z) (z-mean(z))/sd(z)))
# rna$counts1 <- NumericList(lapply(rna$counts1, scale))
sc_promoter$counts2 <- NumericList(lapply(sc_promoter$counts2, function(z) (z-mean(z))/sd(z)))

save(sc_rna, file = "data/sc_rna.rda", compress = "xz")
save(sc_promoter, file = "data/sc_promoter.rda", compress = "xz")

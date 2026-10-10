# Supply an existing Seurat object and exact cell-ID annotation file.
attach_pgam5_v6 <- function(object, annotation_csv) {
  metadata <- read.csv(gzfile(annotation_csv), stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(!anyDuplicated(metadata$cell_id), !anyDuplicated(colnames(object)))
  idx <- match(colnames(object), metadata$cell_id)
  if (anyNA(idx)) stop("Unmatched cell IDs: confirm HCC scope and sample/barcode prefixes")
  fields <- c("cell_type", "PGAM5_RNA_status", "PGAM5_macrophage_annotation",
              "PGAM5_related_state_annotation", "dominant_program", "candidate_program_score",
              "candidate_state_member", "usable_for_validated_TCGA_PGAM5_cell_abundance")
  added <- metadata[idx, fields, drop = FALSE]
  rownames(added) <- colnames(object)
  names(added) <- paste0("pgam5_v6_", names(added))
  SeuratObject::AddMetaData(object, metadata = added)
}

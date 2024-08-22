create_loupe_file <- function(h5ad_filename, output_dir, output_name, projection_names, categorical_names) {
    library(Seurat)
    library(SeuratData)
    library(SeuratDisk)
    library(loupeR)
    library(rhdf5)
    
    Convert(paste(h5ad_filename, ".h5ad", sep=""), dest = paste(h5ad_filename, ".h5seurat", sep=""), overwrite = TRUE)
    sobj <- LoadH5Seurat(paste(h5ad_filename, ".h5seurat", sep=""), assays="RNA", meta.data = FALSE, misc = FALSE)

    clusters <- list()
    projections <- list()

    for (projection in projection_names) {
        projections[[projection]] <- Embeddings(sobj, reduction = projection)
    }

    for (group in categorical_names) {
        codes <- h5read(paste(h5ad_filename,".h5ad", sep=""), paste("/obs", group, "codes", sep="/"))
        categories = h5read(paste(h5ad_filename, ".h5ad", sep=""), paste("/obs", group, "categories", sep="/"))

        unique_codes = unique(codes)

        if (length(unique_codes) == length(categories) + 1) {
            if ("-1" %in% unique_codes) {
                categories <- c("NA", categories)
            } else {
                stop("Wrong number of categories (non-NA).")
            }
            levels = categories
        } else if (length(unique_codes) != length(categories)) {
            stop("Wrong number of categories.")
        }

        clusters[[group]] <- factor(
            codes, labels=categories,
        )
    }

    create_loupe(
        count_mat=sobj@assays$RNA$counts,
        clusters=clusters,
        projections=projections,
        output_dir=output_dir,
        output_name=output_name,
        feature_ids=rownames(sobj@assays$RNA),
        force=TRUE
    )
}


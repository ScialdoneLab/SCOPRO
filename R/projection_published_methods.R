#' findCellTypesSeurat
#'
#' @param queryObj Seurat object of the in vitro dataset
#' @param referenceObj Seurat object of the in vivo dataset
#' @param k.anchor k.anchor parameter in the function \emph{FindTransferAnchors}
#' @param k.filter k.filter parameter in the function \emph{FindTransferAnchors}
#' @param k.weight k.weight parameter in the function \emph{TransferData}
#' @param namedLabels refdata parameter in the function \emph{FindTransferAnchors}
#' @return  Seurat object
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#' @description Projection done by functions \emph{FindTransferAnchors} and \emph{TransferData} from package \emph{Seurat}
#' @export findCellTypesSeurat
#' @seealso \url{https://CRAN.R-project.org/package=Seurat}



findCellTypesSeurat <- function(queryObj, referenceObj, k.anchor = 5, k.filter = 200, namedLabels, k.weight = 50){

  if (!(requireNamespace("Seurat", quietly = TRUE))) {
    stop("Package Seurat needed for this function to work. Please install it: install.packages('Seurat')")
  }

  anchors <- Seurat::FindTransferAnchors(reference = referenceObj, query = queryObj,k.anchor = k.anchor,k.filter = k.filter)
  predictions <- Seurat::TransferData(anchorset = anchors, refdata = namedLabels, k.weight = k.weight)
  queryObj <- Seurat::AddMetaData(object = queryObj, metadata = predictions)

  return(queryObj)
}



#' findCellTypes_scmap
#'
#' @inheritParams findCellTypesSeurat
#' @param n_features selectFeatures parameter in the function \emph{selectFeatures}
#' @param threshold threshold parameter in the function \emph{scmapCluster}
#' @return  Seurat object
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#' @description Projection done by function \emph{scmapCluster} from package \emph{scmap}
#' @export findCellTypes_scmap
#' @seealso \url{https://bioconductor.org/packages/release/bioc/html/scmap.html}


findCellTypes_scmap <- function(queryObj, referenceObj, namedLabels, n_features = 500, threshold = 0.7) {


  if (sum(!unlist(lapply(list("scmap","SingleCellExperiment","SummarizedExperiment","S4Vectors"),requireNamespace, quietly = TRUE)))>0){
    stop("Package scmap and SingleCellExperiment needed for this function to work. Please install them: BiocManager::install('scmap'), BiocManager::install('SingleCellExperiment'), BiocManager::install('SummarizedExperiment') and BiocManager::install('S4Vectors')")
  }
  norm_vivo <- as.matrix(Seurat::GetAssayData(referenceObj, slot = "data",assay="RNA"))
  raw_vivo <- as.matrix(Seurat::GetAssayData(referenceObj, slot = "counts",assay="RNA"))

  norm_vitro <- as.matrix(Seurat::GetAssayData(queryObj, slot = "data",assay="RNA"))

  raw_vitro <- as.matrix(Seurat::GetAssayData(queryObj, slot = "counts",assay="RNA"))

  sce_pro <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = as.matrix(norm_vitro)))

  SingleCellExperiment::logcounts(sce_pro) <- (SingleCellExperiment::normcounts(sce_pro))
  SummarizedExperiment::rowData(sce_pro)$feature_symbol <- rownames(sce_pro)
  SingleCellExperiment::counts(sce_pro) <- as.matrix(raw_vitro)



  cell_type1 <- as.data.frame(namedLabels)
  row.names(cell_type1) <- colnames(norm_vivo)
  colnames(cell_type1) <- "cells_type1"

  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts =  as.matrix(norm_vivo)), colData = cell_type1)
  SingleCellExperiment::logcounts(sce) <- (SingleCellExperiment::normcounts(sce))
  SingleCellExperiment::counts(sce) <- (raw_vivo)


  SummarizedExperiment::rowData(sce)$feature_symbol <- rownames(sce)
  sce <- sce[!duplicated(rownames(sce)), ]
  sce
  sce <- scmap::selectFeatures(sce, suppress_plot = TRUE,n_features = n_features)



  sce <- scmap::indexCluster(sce,cluster_col = "cells_type1")



  scmapCluster_results <- scmap::scmapCluster(
    projection = sce_pro,
    index_list = list(
      yan = S4Vectors::metadata(sce)$scmap_cluster_index
    ),threshold = threshold




  )

  result_label <- as.vector(scmapCluster_results$scmap_cluster_labs)

  queryObj <- Seurat::AddMetaData(object = queryObj, metadata = result_label, col.name = "predictions")

  return(queryObj)
}


#' findCellTypesScibet
#'
#' @inheritParams findCellTypesSeurat
#' @return  Seurat object
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#' @description Projection done by function \emph{SciBet_R} from package \emph{scibetR}
#' @export findCellTypesScibet
#' @seealso \url{https://github.com/zwj-tina/scibetR}


findCellTypesScibet <- function(queryObj, referenceObj, namedLabels){
  if (!(requireNamespace(c("scibetR"), quietly = TRUE))) {
    stop("Package scibetR needed for this function to work. Please install it: devtools::install_github('PaulingLiu/scibet')")
  }
  norm_vitro <- as.matrix(Seurat::GetAssayData(queryObj, slot = "data",assay="RNA"))
  norm_vitro <- t(norm_vitro)
  norm_vitro <- apply(norm_vitro, 2, as.numeric)
  norm_vitro <- as.data.frame(norm_vitro)
  norm_vivo <- as.data.frame(Seurat::GetAssayData(referenceObj, slot = "data",assay="RNA"))
  norm_vivo <- t(norm_vivo)
  norm_vivo <- apply(norm_vivo, 2, as.numeric)
  norm_vivo <- as.data.frame(norm_vivo)
  norm_vivo_full <- cbind(norm_vivo,as.vector(namedLabels))
  colnames(norm_vivo_full)=c(colnames(norm_vivo),"label")
  norm_vivo_full <- as.data.frame(norm_vivo_full)
  row.names(norm_vivo_full)=NULL

  predictions <- scibetR::SciBet_R(norm_vivo_full,norm_vitro)
  queryObj <- Seurat::AddMetaData(object = queryObj, metadata = predictions, col.name = "predictions")
  return(queryObj)
}




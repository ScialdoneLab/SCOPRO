
#' plot_score_genes
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @param markers_to_plot Character vector with the names of the genes to plot.
#' @param label_1 Character value. Label for the in vitro dataset
#' @param label_2 Character value. Label for the in vivo dataset
#' @param final_name Character vector with the names of the genes to show in the plot.
#' @param max_size Numeric value, specifying the size of the dot.
#' @param text_size Numeric value, specifying the size of the text in the plot.
#' @param title_name Character value.

#' @return ggplot2::ggplot2 object showing balloon plot.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export plot_score_genes
#'


plot_score_genes <- function(markers_to_plot ,label_1, label_2, norm_vitro,norm_vivo, cluster_vitro,cluster_vivo, final_name, max_size = 9, text_size= 9.5, title_name)
{

  markers_to_plot <- markers_to_plot[markers_to_plot %in% row.names(norm_vitro)]



  all_markers_common_name <- final_name

  label_mouse <- c(rep(label_1,length(colnames(norm_vitro))),rep(label_2,length(colnames(norm_vivo))))
  mouse_norm_common=cbind(norm_vitro[markers_to_plot,],norm_vivo[markers_to_plot,])




  cluster_vitro_name <- paste(label_1, factor(cluster_vitro), sep="-")
  cluster_vivo_name <- paste(label_2, factor(cluster_vivo), sep="-")
  mouse_norm_common <- cbind(norm_vitro[markers_to_plot, ], norm_vivo[markers_to_plot, ])
  label_mouse <- c(cluster_vitro_name, cluster_vivo_name)
  out_2 <- plot_balons(mouse_norm_common, label_mouse,markers_to_plot, length(markers_to_plot), length(colnames(mouse_norm_common)), max_size = max_size,text_size = text_size, all_markers_common_name, keep_order = FALSE)
  out_2 <- out_2 + ggplot2::ggtitle(title_name)

  return(out_2)
}


#' plot_score
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams  select_common_genes
#' @inheritParams plot_score_genes
#' @param y_name Character value
#' @param fill_name Character value.

#' @return ggplot2::ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export plot_score
#'

plot_score <- function(SCOPRO_output, marker_stages, marker_stages_filter, selected_stages, name_vivo, y_name, fill_name, title_name){
  final_score <- c(as.vector(unlist(SCOPRO_output[[5]])))
  label <- c(names(unlist(SCOPRO_output[[5]])))
  df <- data.frame(final_score,label)

  index_specific <- which(selected_stages %in% name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = label, y = final_score, fill = label)) +
    ggplot2::geom_bar(stat = "identity") + ggplot2::theme_minimal()
  p + ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
             axis.title.x = ggplot2::element_blank()) + ggplot2::ylab(y_name) + ggplot2::labs(fill=fill_name) +  ggplot2::ylim(c(0,1)) + ggplot2::ggtitle(paste0(title_name," (",sum(marker_stages_filter %in% marker_specific)," genes)",sep=" "))
}


#' plot_score_sc
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams  select_common_genes
#' @inheritParams plot_score_genes
#' @inheritParams plot_score
#' @param SCOPRO_output output given by function \emph{SCOPRO_single_cell}

#' @return ggplot2::ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export plot_score_sc
#'

plot_score_sc <- function (SCOPRO_output, cluster_vitro, marker_stages, marker_stages_filter, selected_stages, name_vivo, y_name, fill_name, title_name)
{
  final_score <- c(as.vector(unlist(SCOPRO_output[[5]])))
  label <- cluster_vitro
  df <- data.frame(final_score, label)
  index_specific <- which(selected_stages %in% name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = label, y = final_score,
                                        fill = label)) + ggplot2::geom_boxplot()+ggplot2::geom_jitter(position=ggplot2::position_jitter(0.2)) +
    ggplot2::theme_minimal()
  p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         vjust = 1, hjust = 1), axis.title.x = ggplot2::element_blank()) +
    ggplot2::ylab(y_name) + ggplot2::labs(fill = fill_name) +
    ggplot2::ylim(c(-0.1, 1)) + ggplot2::ggtitle(paste0(title_name,
                                                        " (", sum(marker_stages_filter %in% marker_specific),
                                                        " genes)", sep = " "))
}




#' plot_balons
#' @noRd
#'
plot_balons <- function(norm_counts, final_cluster, genes_to_plot, max_number, total_cells, max_size=5, text_size=7,label_final, keep_order=TRUE){

  fattore <- factor(final_cluster,levels=unique(final_cluster))
  if (keep_order == FALSE){
    fattore <- factor(final_cluster)}
  livelli <- levels(factor(fattore))



  all_markers <- genes_to_plot

  all_markers <- lapply(all_markers, function(x) {x[x!=0]})


  all_markers_values <- rep(list(0),length(livelli))
  for ( i in 1:length(livelli)){
    gene_values_0 <- apply(norm_counts[genes_to_plot,final_cluster==livelli[i]], 1, mean)
    all_markers_values[[i]] <- gene_values_0
  }


  values <- rep(list(0),length(livelli))
  for ( i in 1:length(livelli)){
    valore=apply(norm_counts[genes_to_plot,final_cluster==livelli[i]],1,function(x){sum(x>0.1)/length(x)})
    values[[i]]=valore
  }

  cluster_label_all=rep(list(0),length(livelli))
  for ( i in 1:length(livelli)){
    cluster_label_all[[i]] <- rep(livelli[i],length(genes_to_plot))
  }


  all_markers_plot <- genes_to_plot
  cluster_label_all <- unlist(cluster_label_all)
  values <- unlist(values)
  all_markers_values <- unlist(all_markers_values)
  balloon_melted <- data.frame(all_markers_plot,cluster_label_all,values,label_final)
  colnames(balloon_melted) <- c("genes_plot","cluster_label_all","values","label_all_final")

  p <- ggplot2::ggplot(balloon_melted, ggplot2::aes(x =factor(cluster_label_all,levels=unique(cluster_label_all)), y = factor(balloon_melted[,4],levels=(unique(label_final)))))
  p + ggplot2::geom_point( ggplot2::aes(size=values,colour=as.numeric(all_markers_values))) + ggplot2::theme(panel.background=ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "blue", fill=NA, size=3))+ggplot2::scale_size_area(max_size=max_size)+ ggplot2::scale_colour_gradient(low = "grey", high = "brown", na.value = NA)+ggplot2::labs(colour="Mean",size="Value")+ggplot2::xlab("Cluster")+ggplot2::ylab("Markers")+ggplot2::theme(text = ggplot2::element_text(size=text_size),axis.text.x = ggplot2::element_text(angle=45,vjust=1,hjust=1),axis.title.x=ggplot2::element_blank())
}


#' cellTypesPerClusterBalloonPlot
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams  select_common_genes
#' @inheritParams plot_score_genes
#' @description All the stages in \emph{cluster_vivo_factor} are showns in the balloon plot
#' @param obj Seurat object given as output from function \emph{findCellTypesSeurat}
#' @param cluster_vivo_factor Factor vector specifyng the cluster partition of the in vivo datataset
#' @param order_label_vivo Character vector specifyng the order of the columns in the balloon plot. The names must be present in \emph{cluster_vivo_factor}
#' @param label_size Numeric value specifyng the label.size parameter in the function \emph{balloonplot} from package \emph{gplots}.
#' @param thresold Numeric value. Cells with a predicted score below or equal to \emph{thresold} are labelled as unassigned. Only used when \emph{method} is Seurat
#' @param method Character name. Must be one of Seurat, scibetR or scmap

#' @return  balloon plot given by function \emph{balloonplot} from package \emph{gplots}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export cellTypesPerClusterBalloonPlot
#' @seealso \url{https://CRAN.R-project.org/package=gplots}

cellTypesPerClusterBalloonPlot <- function (obj, cluster_vivo_factor, order_label_vivo, method = "Seurat", title_name, text_size = 0.4, label_size = 0.4, thresold = 0.5)
{if (method == "Seurat"){
  if (!(requireNamespace("gplots", quietly = TRUE))) {
    stop("Package gplots needed for this function to work. Please install it: install.packages('gplots')")
  }
  obj@meta.data$predicted.id[obj@meta.data$prediction.score.max <=
                               thresold] <- "Unassigned"
  levels_vitro <- c(levels(cluster_vivo_factor))
  levels_vitro_factor <- factor(levels_vitro, levels = order_label_vivo)
  if (sum(obj@meta.data$predicted.id == "Unassigned") > 0) {
    levels_vitro <- c(levels(cluster_vivo_factor), "Unassigned")
    levels_vitro_factor <- factor(levels_vitro, levels = c(order_label_vivo,
                                                           "Unassigned"))
  }
  cellTypes <- rep(levels(levels_vitro_factor), each = length(unique(as.vector(obj$seurat_clusters))))
  clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(levels(levels_vitro_factor)))
  nElements <- numeric()
  for (cellType in levels(levels_vitro_factor)) {
    nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predicted.id) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
  }
  data <- data.frame(cellTypes, clusters, nElements)
  for (cluster in unique(data$clusters)) {
    data$nElements[which(data$clusters == cluster)] <- 100 *
      data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters == cluster)])
  }
  return(gplots::balloonplot(data$cellTypes, data$clusters, data$nElements, label.size = label_size, text.size = text_size,
                             ylab = "", xlab = "", main = title_name))
}
  if (method == "scibetR"){
    if (!(requireNamespace("gplots", quietly = TRUE))) {
      stop("Package gplots needed for this function to work. Please install it: install.packages('gplots')")
    }
    levels_vitro <- c(levels(cluster_vivo_factor))
    levels_vitro_factor <- factor(levels_vitro, levels = order_label_vivo)
    cellTypes <- rep(levels(levels_vitro_factor), each = length(unique(as.vector(obj$seurat_clusters))))
    clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(levels(levels_vitro_factor)))
    nElements <- numeric()
    for (cellType in levels(levels_vitro_factor)) {
      nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predictions) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
    }
    data <- data.frame(cellTypes, clusters, nElements)
    for (cluster in unique(data$clusters)) {
      data$nElements[which(data$clusters == cluster)] <- 100 *
        data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters ==  cluster)])}
    return(gplots::balloonplot(data$cellTypes, data$clusters,
                               data$nElements, label.size = label_size, text.size = text_size,
                               ylab = "", xlab = "", main = title_name))
  }
  if (method == "scmap"){
    if (!(requireNamespace("gplots", quietly = TRUE))) {
      stop("Package gplots needed for this function to work. Please install it: install.packages('gplots')")
    }
    levels_vitro <- c(levels(cluster_vivo_factor), "unassigned")
    levels_vitro_factor <- factor(levels_vitro, levels = c(order_label_vivo,
                                                           "unassigned"))
    cellTypes <- rep(levels(levels_vitro_factor), each = length(unique(as.vector(obj$seurat_clusters))))
    clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(levels(levels_vitro_factor)))
    nElements <- numeric()
    for (cellType in levels(levels_vitro_factor)) {
      nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predictions) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
    }
    data <- data.frame(cellTypes, clusters, nElements)
    for (cluster in unique(data$clusters)) {
      data$nElements[which(data$clusters == cluster)] <- 100 *
        data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters == cluster)])
    }
    return(gplots::balloonplot(data$cellTypes, data$clusters,
                               data$nElements, label.size = label_size, text.size = text_size,
                               ylab = "", xlab = "", main = title_name))

  }
  if (method != "Seurat" & method != "scibetR" & method != "scmap"){
    stop("method must be equal to Seurat, ScibetR or scmap.")
  }

}





#' cellTypesPerClusterBalloonPlot_small
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams  select_common_genes
#' @inheritParams plot_score_genes
#' @inheritParams cellTypesPerClusterBalloonPlot
#' @return  balloon plot given by function \emph{balloonplot} from package \emph{gplots}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#' @description Only the stages in \emph{cluster_vivo_factor} with an assignment greater than zero are showns in the balloon plot
#' @export cellTypesPerClusterBalloonPlot_small
#' @seealso \url{https://CRAN.R-project.org/package=gplots}




cellTypesPerClusterBalloonPlot_small <- function (obj, cluster_vivo_factor, order_label_vivo, title_name, method = "Seurat",text_size = 0.4,
                                                  label_size = 0.7, thresold = 0.5)
{
  if (method == "Seurat"){
    if (!(requireNamespace("gplots", quietly = TRUE))) {
      stop("Package gplots needed for this function to work. Please install it: install.packages('gplots')")
    }
    obj@meta.data$predicted.id[obj@meta.data$prediction.score.max <= thresold] <- "Unassigned"
    levels_vitro <- c(levels(cluster_vivo_factor))
    levels_vitro_factor <- factor(levels_vitro, levels = c(order_label_vivo))
    if (sum(obj@meta.data$predicted.id == "Unassigned") > 0) {
      levels_vitro <- c(levels(cluster_vivo_factor), "Unassigned")
      levels_vitro_factor <- factor(levels_vitro, levels = c(order_label_vivo,"Unassigned"))
    }
    cellTypes <- rep(levels_vitro[levels_vitro %in% unique(as.vector(obj@meta.data$predicted.id))], each = length(unique(as.vector(obj$seurat_clusters))))
    clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(unique(as.vector(obj@meta.data$predicted.id))))
    nElements <- numeric()
    for (cellType in levels_vitro[levels_vitro %in%
                                  unique(as.vector(obj@meta.data$predicted.id))]) {
      nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predicted.id) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
    }
    data <- data.frame(cellTypes, clusters, nElements)
    for (cluster in unique(data$clusters)) {
      data$nElements[which(data$clusters == cluster)] <- 100 * data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters == cluster)])
    }
    return(gplots::balloonplot(data$cellTypes, data$clusters, data$nElements, label.size = label_size, text.size = text_size,
                               ylab = "", xlab = "", main = title_name))
  }
  if (method == "scibetR"){
    if (!(requireNamespace("gplots", quietly = TRUE))) {
      stop("Package gplots needed for this function to work. Please install it: install.packages('gplots')")
    }
    cellTypes <- rep(levels(cluster_vivo_factor)[levels(cluster_vivo_factor) %in% unique(as.vector(obj@meta.data$predictions))], each = length(unique(as.vector(obj$seurat_clusters))))
    clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(unique(as.vector(obj@meta.data$predictions))))
    nElements <- numeric()
    for (cellType in levels(cluster_vivo_factor)[levels(cluster_vivo_factor) %in% unique(as.vector(obj@meta.data$predictions))]) {
      nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predictions) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
    }
    data <- data.frame(cellTypes, clusters, nElements)
    for (cluster in unique(data$clusters)) {
      data$nElements[which(data$clusters == cluster)] <- 100 * data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters == cluster)])
    }
    print(data)
    return(gplots::balloonplot(data$cellTypes, data$clusters,
                               data$nElements, label.size = label_size, text.size = text_size,
                               ylab = "", xlab = "", main = title_name))
  }
  if (method == "scmap"){
    levels_vitro <- c(levels(cluster_vivo_factor), "unassigned")
    levels_vitro_factor <- factor(levels_vitro, levels =  c(order_label_vivo,"unassigned"))
    cellTypes <- rep(levels_vitro[levels_vitro %in%
                                    unique(as.vector(obj@meta.data$predictions))], each = length(unique(as.vector(obj$seurat_clusters))))
    clusters <- rep(unique(as.vector(obj$seurat_clusters)), length(unique(as.vector(obj@meta.data$predictions))))
    nElements <- numeric()
    for (cellType in levels_vitro[levels_vitro %in% unique(as.vector(obj@meta.data$predictions))]) {
      nElements <- c(nElements, table(factor(as.vector(obj$seurat_clusters)[which(as.vector(obj@meta.data$predictions) == cellType)], levels = unique(as.vector(obj$seurat_clusters)))))
    }
    data <- data.frame(cellTypes, clusters, nElements)
    for (cluster in unique(data$clusters)) {
      data$nElements[which(data$clusters == cluster)] <- 100 *
        data$nElements[which(data$clusters == cluster)]/sum(data$nElements[which(data$clusters ==
                                                                                   cluster)])
    }
    return(gplots::balloonplot(data$cellTypes, data$clusters,
                               data$nElements, label.size = label_size, text.size = text_size,
                               ylab = "", xlab = "", main = title_name))
  }
  if (method != "Seurat" & method != "scibetR" & method != "scmap"){
    stop("method must be equal to Seurat, scibetR or scmap")
  }
}



#' plot_boxplot
#'
#' @inheritParams SCOPRO
#' @inheritParams select_top_markers
#' @inheritParams select_common_genes
#' @inheritParams plot_score_genes
#' @inheritParams plot_score
#' @inheritParams cellTypesPerClusterBalloonPlot
#' @param norm_vivo Norm count matrix (n_genes X n_cells)
#' @param cluster_vivo cluster partition
#' @param gene_name Character name of the gene of interest
#' @param x_name Character name of x axis.
#' @return  ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#' @description Only the stages in \emph{cluster_vivo_factor} with an assignment greater than zero are showns in the balloon plot
#' @export plot_boxplot
#' @seealso \url{https://CRAN.R-project.org/package=gplots}


plot_boxplot <- function(norm_vivo, gene_name, cluster_vivo, order_label_vivo, title_name, fill_name, x_name, y_name){
  gene_expr <- norm_vivo[gene_name,]
  cluster_vivo <- factor(cluster_vivo,levels = order_label_vivo)
  data_common <- data.frame(cluster_vivo,gene_expr)
  ggplot2::ggplot(data_common, ggplot2::aes(x=cluster_vivo, y=gene_expr,fill= cluster_vivo)) +
    ggplot2::geom_boxplot() + ggplot2::geom_jitter(position=ggplot2::position_jitter(0.2)) + ggplot2::xlab(x_name) + ggplot2::ylab(y_name) + ggplot2::ggtitle(title_name) + ggplot2::labs(fill = fill_name)
}







#' plot_in_vivo_markers
#'
#' @inheritParams SCOPRO
#' @inheritParams plot_score_genes
#' @param coordinate_umap data frame with umap coordinates (UMAP 1 and UMAP2).
#' @param threshold Numeric value.
#' @description Cells are coloured according to the number of markers of the in vivo stage \emph{name_vivo} that are epxressed above \emph{threshold}
#' @return ggplot2::ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#'
#'
#'
#' @export plot_in_vivo_markers
#'
plot_in_vivo_markers <- function (coordinate_umap, norm_vitro, marker_stages, selected_stages, name_vivo, threshold = 0, title_name){
  if (sum(!unlist(lapply(list(c("grDevices","gplots")),requireNamespace, quietly = TRUE)))>0) {
    stop("Package gplots and grDevices needed for this function to work. Please install them: install.packages('gplots') and install.packages('grDevices')")
  }
  rank_intersect = detect_expressed_genes(norm_vitro, marker_stages, selected_stages, name_vivo, 0)
  ramp <- grDevices::colorRamp(c("white", "darkorange3"))
  length_ramp <- length(unique(rank_intersect))
  ramp.list <- grDevices::rgb( ramp(seq(0, 1, length = length_ramp)), max = 255)
  index_color=round(length(ramp.list)/2,0)
  breaks <- seq(0,1,length.out=1000)
  gradient1 <- gplots::colorpanel( sum( breaks[-1]<= as.numeric(stats::quantile(breaks,0.50))), "#FFFFFF",ramp.list[index_color])
  gradient2 <- gplots::colorpanel( sum( breaks[-1] > as.numeric(stats::quantile(breaks,0.50)) ), ramp.list[index_color], "#CD6600" )
  hm.colors_esc <- c(gradient1,gradient2)
  index_sort <- order(rank_intersect)
  row.names(coordinate_umap) <- colnames(norm_vitro)
  coordinate_umap <- coordinate_umap[colnames(norm_vitro)[index_sort], ]
  rank_intersect <- sort(rank_intersect)
  umap_plot <- ggplot2::ggplot(coordinate_umap, ggplot2::aes(coordinate_umap[,
                                                                             1], coordinate_umap[, 2])) + ggplot2::geom_point(ggplot2::aes(colour = rank_intersect), size = 1) + ggplot2::xlab("UMAP 1") + ggplot2::ylab("UMAP 2") + ggplot2::ggtitle(title_name) + ggplot2::labs(col= paste0("Number of markers")) + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), text = ggplot2::element_text(size = 14, family = "Helvetica")) + ggplot2::scale_colour_gradientn(colors = hm.colors_esc)
  return((umap_plot))
}


#' plot_in_vivo_markers_interactive
#'
#' @inheritParams SCOPRO
#' @inheritParams plot_score_genes
#' @inheritParams plot_in_vivo_markers
#' @param max_number Integer. Maximum number of genes name to show.
#' @param min_x Set the min limit on the x axis.
#' @param max_x Set the max limit on the x axis.
#' @param min_y Set the min limit on the y axis.
#' @param max_y Set the min limit on the y axis.
#' @description Cells are coloured according to the number of markers of the in vivo stage \emph{name_vivo} that are epxressed above \emph{threshold}.
#' It is based on plotly library.
#' @return plotly object given by \emph{plot_ly function} (from library \emph{plotly}).
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://plotly.com/r/}
#'
#'
#'
#' @export plot_in_vivo_markers_interactive
#'



plot_in_vivo_markers_interactive <- function (coordinate_umap, norm_vitro, marker_stages, selected_stages, name_vivo, threshold = 0, max_number = 5, min_x = NULL, max_x = NULL, min_y = NULL, max_y = NULL){
  if (sum(!unlist(lapply(list(c("grDevices","plotly")),requireNamespace, quietly = TRUE)))>0) {
    stop("Package plotly and grDevices needed for this function to work. Please install them: install.packages('plotly') and install.packages('grDevices')")
  }
  rank_intersect <- detect_expressed_genes(norm_vitro, marker_stages, selected_stages, name_vivo, threshold)
  ramp <- grDevices::colorRamp(c("white", "darkorange3"))
  length_ramp <- length(unique(rank_intersect))
  ramp.list <- grDevices::rgb( ramp(seq(0, 1, length = length_ramp)), max = 255)
  index_sort <- order(rank_intersect)
  list_intersect <- detect_names_genes(norm_vitro, marker_stages, selected_stages, name_vivo, threshold)
  index_specific <- which(selected_stages %in% name_vivo)
  marker_specific <- marker_stages[[index_specific]]
  marker_specific <- marker_specific[marker_specific %in% row.names(norm_vitro)]
  text <- names_marker_genes(list_intersect, marker_specific, max_number)
  text <- text[index_sort]
  row.names(coordinate_umap) <- colnames(norm_vitro)
  coordinate_umap <- coordinate_umap[colnames(norm_vitro)[index_sort], ]
  rank_intersect <- sort(rank_intersect)
  colnames(coordinate_umap) <- c("UMAP_1", "UMAP_2")
  fig <- plotly::plot_ly(data = coordinate_umap, x = ~UMAP_1,
                         y = ~UMAP_2, color = ~rank_intersect, type = "scatter",
                         mode = "markers", marker = list(size = 5, width = 2,
                                                         line = list(color = "black", width = 0.5)), text = ~text,
                         hoverinfo = "text", colors = ramp.list)
  if (!is.null(min_x)) {
    fig <- fig %>% plotly::layout(xaxis = list(range = c(min_x,
                                                         max_x)), yaxis = list(range = c(min_y, max_y))) %>%
      plotly::layout(plot_bgcolor = "rgb(254, 247, 234)") %>%
      plotly::layout(paper_bgcolor = "rgb(254, 247, 234)")
  }
  return(fig)
}




#' heatmap_markers_vivo
#'
#' @inheritParams SCOPRO
#' @inheritParams plot_score_genes
#' @inheritParams plot_in_vivo_markers
#' @param max_number Integer. Maximum number of genes to show for each \emph{selected_stages}.
#' @param color_vector Character vector with colour assignment.
#' @param median_profile Logical value. If TRUE, then the median expression profile for the given markers in each stage is shown
#' @return Heatmap class object
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
#'
#'
#'
#' @export heatmap_markers_vivo
#'

heatmap_markers_vivo <- function(norm_vivo, marker_stages, selected_stages, max_number = 10, cluster_vivo, median_profile = FALSE, color_vector = NULL){
  if (median_profile == FALSE){
    if (sum(!unlist(lapply(list(c("ComplexHeatmap","circlize")),requireNamespace, quietly = TRUE)))>0) {
      stop("Package circlize and ComplexHeatmap needed for this function to work. Please install them: install.packages('circlize') and BiocManager::install('ComplexHeatmap')")
    }
    for (i in selected_stages){
      index_specific <- which(selected_stages %in% i)
      marker_specific <- marker_stages[[index_specific]]
      if (length(marker_specific) == 0) {
        stop(paste0("There are no markers for the stage: ", i))
      }
    }
    marker_stages_top <- rep(list(0),length(selected_stages))
    for (i in 1:length(selected_stages)){
      if (length(marker_stages[[i]]) >= max_number) {
        marker_stages_top[[i]] <- marker_stages[[i]][1: max_number]
      }
      else {marker_stages_top[[i]] <- marker_stages[[i]]}
    }
    marker_all <- unlist(marker_stages_top)
    cluster_vivo <- factor(cluster_vivo, levels = selected_stages)
    index_list <- rep(list(0), length(levels(cluster_vivo)))
    for ( i in 1:length(levels(cluster_vivo))){
      index_list[[i]] <- which(cluster_vivo == levels(cluster_vivo)[i])
    }
    index_all <- unlist(index_list)
    heatdata <- norm_vivo[marker_all, index_all]
    cluster_vivo <- cluster_vivo[index_all]
    cluster_unique <- unique(cluster_vivo)
    row.names(heatdata) <- marker_all
    if(is.null(color_vector)){
      color_cluster=rep(0,length(unique(cluster_vivo)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- gg_color_hue(length(color_cluster))[i]
      }
      names(color_cluster) <- as.character(cluster_unique)}
    if(!is.null(color_vector)){
      color_cluster <- rep(0,length(unique(cluster_vivo)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- color_vector[i]
      }
      names(color_cluster) <- as.character(cluster_unique)
    }
    color_condition <- "green"
    batch <- rep("In vivo dataset", length(cluster_vivo))
    data_heatmap=data.frame(Cluster = cluster_vivo, Condition = batch)
    haCol1 <- ComplexHeatmap::HeatmapAnnotation(df = data_heatmap,col=list(Cluster=color_cluster,Condition = c("In vivo dataset" = "green 1")), show_legend = T)
    ht21 <- ComplexHeatmap::Heatmap(as.matrix(heatdata)
                                    , col= circlize::colorRamp2(c(0, round(max(heatdata),4)), c("grey", "brown"))
                                    , name = "log norm counts"
                                    #, column_title = "Absolute values"
                                    #, cluster_columns = as.dendrogram(my.tree)
                                    ,cluster_columns = F
                                    , cluster_rows =F
                                    , top_annotation = haCol1
                                    , row_names_gp = grid::gpar(fontsize = 6)
                                    , show_column_names = F
                                    , show_row_names = T
    )
    ComplexHeatmap::draw(ht21)
  }
  if (median_profile == TRUE){
    if (sum(!unlist(lapply(list(c("ComplexHeatmap","circlize")),requireNamespace, quietly = TRUE)))>0) {
      stop("Package circlize and ComplexHeatmap needed for this function to work. Please install them: install.packages('circlize') and BiocManager::install('ComplexHeatmap')")
    }
    for (i in selected_stages){
      index_specific <- which(selected_stages %in% i)
      marker_specific <- marker_stages[[index_specific]]
      if (length(marker_specific) == 0) {
        stop(paste0("There are no markers for the stage: ", i))
      }
    }
    marker_stages_top <- rep(list(0),length(selected_stages))
    for (i in 1:length(selected_stages)){
      if (length(marker_stages[[i]]) >= max_number) {
        marker_stages_top[[i]] = marker_stages[[i]][1: max_number]
      }
      else {marker_stages_top[[i]] = marker_stages[[i]]}
    }
    marker_all <- unlist(marker_stages_top)
    median_es <- find_mean_cluster(cluster_vivo, selected_stages, norm_vivo, marker_all, method = "median")
    median_es_matrix <- matrix(unlist(median_es), ncol = length(median_es), byrow = FALSE)
    colnames(median_es_matrix) <- selected_stages
    row.names(median_es_matrix) <- marker_all
    heatdata <- median_es_matrix
    cluster_vivo <- selected_stages
    cluster_unique=unique(cluster_vivo)
    row.names(heatdata) <- marker_all
    if(is.null(color_vector)){
      color_cluster <- rep(0,length(unique(cluster_vivo)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- gg_color_hue(length(color_cluster))[i]
      }
      names(color_cluster) <- as.character(cluster_unique)}
    if(!is.null(color_vector)){
      color_cluster=rep(0,length(unique(cluster_vivo)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- color_vector[i]
      }
      names(color_cluster) <- as.character(cluster_unique)
    }
    color_condition <- "green"
    batch <- rep("In vivo dataset", length(cluster_vivo))
    data_heatmap=data.frame(Cluster = cluster_vivo, Condition = batch)
    haCol1 <- ComplexHeatmap::HeatmapAnnotation(df = data_heatmap,col=list(Cluster=color_cluster,Condition = c("In vivo dataset" = "green 1")), show_legend = T)
    ht21 <- ComplexHeatmap::Heatmap(as.matrix(heatdata)
                                    , col= circlize::colorRamp2(c(0, round(max(heatdata),4)), c("grey", "brown"))
                                    , name = "log norm counts"
                                    #, column_title = "Absolute values"
                                    #, cluster_columns = as.dendrogram(my.tree)
                                    ,cluster_columns = F
                                    , cluster_rows =F
                                    , top_annotation = haCol1
                                    , row_names_gp = grid::gpar(fontsize = 6)
                                    , show_column_names = F
                                    , show_row_names = T
    )
    ComplexHeatmap::draw(ht21)
  }
  if (!is.logical(median_profile)) {
    stop("Parameter median_profile must be logical")
  }
}








#' heatmap_markers_vitro
#'
#' @inheritParams SCOPRO
#' @inheritParams plot_score_genes
#' @inheritParams plot_in_vivo_markers
#' @inheritParams heatmap_markers_vivo
#' @param marker_plot Character vector with the names of the genes to plot.
#' @param mean_profile Logical value. If TRUE, then the mean expression profile for\emph{marker_plot} in each cluster is shown
#' @return Heatmap class object
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
#'
#'
#'
#' @export heatmap_markers_vitro
#'

heatmap_markers_vitro <- function(norm_vitro, marker_plot, cluster_vitro, mean_profile = FALSE, color_vector = NULL){
  if (mean_profile == FALSE){
    if (sum(!unlist(lapply(list(c("ComplexHeatmap","circlize")),requireNamespace, quietly = TRUE)))>0) {
      stop("Package circlize and ComplexHeatmap needed for this function to work. Please install them: install.packages('circlize') and BiocManager::install('ComplexHeatmap')")
    }
    if (sum(marker_plot %in% row.names(norm_vitro)) == 0){
      stop("No provided genes are in the dataset norm_vitro. Please change marker_plot")
    }
    marker_plot <- marker_plot[marker_plot %in% row.names(norm_vitro)]
    cluster_vitro <- factor(cluster_vitro)
    index_list = rep(list(0), length(levels(cluster_vitro)))
    for ( i in 1:length(levels(cluster_vitro))){
      index_list[[i]] = which(cluster_vitro == levels(cluster_vitro)[i])
    }
    index_all <- unlist(index_list)
    heatdata <- norm_vitro[marker_plot, index_all]
    cluster_vitro <- cluster_vitro[index_all]
    cluster_unique <- unique(cluster_vitro)
    row.names(heatdata) <- marker_plot
    if(is.null(color_vector)){
      color_cluster=rep(0,length(unique(cluster_vitro)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- gg_color_hue(length(color_cluster))[i]
      }
      names(color_cluster) <- as.character(cluster_unique)}
    if(!is.null(color_vector)){
      color_cluster <- rep(0,length(unique(cluster_vitro)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- color_vector[i]
      }
      names(color_cluster) <- as.character(cluster_unique)
    }
    color_condition <- "green"
    batch <- rep("In vitro dataset", length(cluster_vitro))
    data_heatmap=data.frame(Cluster = cluster_vitro, Condition = batch)
    haCol1 <- ComplexHeatmap::HeatmapAnnotation(df = data_heatmap,col=list(Cluster=color_cluster,Condition = c("In vitro dataset" = "yellow 1")), show_legend = T)
    ht21 <- ComplexHeatmap::Heatmap(as.matrix(heatdata)
                                    , col= circlize::colorRamp2(c(0, round(max(heatdata),4)), c("grey", "brown"))
                                    , name = "log norm counts"
                                    #, column_title = "Absolute values"
                                    #, cluster_columns = as.dendrogram(my.tree)
                                    ,cluster_columns = F
                                    , cluster_rows =F
                                    , top_annotation = haCol1
                                    , row_names_gp = grid::gpar(fontsize = 6)
                                    , show_column_names = F
                                    , show_row_names = T
    )
    ComplexHeatmap::draw(ht21)
  }
  if (mean_profile == TRUE){
    if (sum(!unlist(lapply(list(c("ComplexHeatmap","circlize")),requireNamespace, quietly = TRUE)))>0) {
      stop("Package circlize and ComplexHeatmap needed for this function to work. Please install them: install.packages('circlize') and BiocManager::install('ComplexHeatmap')")
    }
    if (sum(marker_plot %in% row.names(norm_vitro)) == 0){
      stop("No provided genes are in the dataset norm_vitro. Please change marker_plot")
    }
    marker_plot <- marker_plot[marker_plot %in% row.names(norm_vitro)]
    selected_stages <- levels(factor(cluster_vitro))
    median_es <- find_mean_cluster(cluster_vitro, selected_stages, norm_vitro, marker_plot, method = "mean")
    median_es_matrix <- matrix(unlist(median_es), ncol = length(median_es), byrow = FALSE)
    colnames(median_es_matrix) <- selected_stages
    row.names(median_es_matrix) <- marker_plot

    heatdata <- median_es_matrix

    cluster_vitro <- selected_stages
    cluster_unique <- unique(cluster_vitro)
    if(is.null(color_vector)){
      color_cluster <- rep(0,length(unique(cluster_vitro)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- gg_color_hue(length(color_cluster))[i]
      }
      names(color_cluster)=as.character(cluster_unique)}
    if(!is.null(color_vector)){
      color_cluster=rep(0,length(unique(cluster_vitro)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- color_vector[i]
      }
      names(color_cluster) <- as.character(cluster_unique)
    }
    color_condition <- "green"
    batch <- rep("In vitro dataset", length(cluster_vitro))
    data_heatmap=data.frame(Cluster = cluster_vitro, Condition = batch)
    haCol1 <- ComplexHeatmap::HeatmapAnnotation(df= data_heatmap,col=list(Cluster=color_cluster,Condition = c("In vitro dataset" = "yellow 1")), show_legend = T)
    ht21 <- ComplexHeatmap::Heatmap(as.matrix(heatdata)
                                    , col= circlize::colorRamp2(c(0, round(max(heatdata),4)), c("grey", "brown"))
                                    , name = "log norm counts"
                                    #, column_title = "Absolute values"
                                    #, cluster_columns = as.dendrogram(my.tree)
                                    ,cluster_columns = F
                                    , cluster_rows =F
                                    , top_annotation = haCol1
                                    , row_names_gp = grid::gpar(fontsize = 6)
                                    , show_column_names = F
                                    , show_row_names = T
    )
    ComplexHeatmap::draw(ht21)
  }
  if (!is.logical(mean_profile)) {
    stop("Parameter mean_profile must be logical")
  }
}



#' plot_score_genes_heatmap
#'
#' @inheritParams SCOPRO
#' @inheritParams plot_score_genes
#' @inheritParams plot_in_vivo_markers
#' @inheritParams heatmap_markers_vivo
#' @inheritParams heatmap_markers_vitro
#' @param marker_plot_name Character vector with the names of the rows to show in the heatmap
#' @return Heatmap class object
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
#'
#'
#'
#' @export plot_score_genes_heatmap
#'


plot_score_genes_heatmap <- function(norm_vivo, cluster_vivo, norm_vitro, cluster_vitro, marker_plot, marker_plot_name, mean_profile = FALSE, color_vector = NULL){
  if (mean_profile == FALSE){
    if (sum(!unlist(lapply(list(c("ComplexHeatmap","circlize")),requireNamespace, quietly = TRUE)))>0) {
      stop("Package circlize and ComplexHeatmap needed for this function to work. Please install them: install.packages('circlize') and BiocManager::install('ComplexHeatmap')")
    }
    if (sum(marker_plot %in% row.names(norm_vitro)) == 0){
      stop("No provided genes are in the dataset norm_vitro. Please change marker_plot")
    }
    if (sum(marker_plot %in% row.names(norm_vivo)) == 0){
      stop("No provided genes are in the dataset norm_vivo. Please change marker_plot")
    }
    marker_plot <- marker_plot[marker_plot %in% row.names(norm_vitro)]
    marker_plot <- marker_plot[marker_plot %in% row.names(norm_vivo)]
    cluster_vitro <- factor(cluster_vitro)
    index_list = rep(list(0), length(levels(cluster_vitro)))
    for ( i in 1:length(levels(cluster_vitro))){
      index_list[[i]] = which(cluster_vitro == levels(cluster_vitro)[i])
    }
    index_all <- unlist(index_list)
    norm_vitro <- norm_vitro[marker_plot, index_all]
    norm_vivo <- norm_vivo[marker_plot,]
    heatdata <- cbind(norm_vivo, norm_vitro)
    cluster_vitro <- cluster_vitro[index_all]
    cluster_all <- c(cluster_vivo, cluster_vitro)
    cluster_unique <- unique(cluster_all)
    row.names(heatdata) <- marker_plot_name
    if(is.null(color_vector)){
      color_cluster=rep(0,length(unique(cluster_all)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- gg_color_hue(length(color_cluster))[i]
      }
      names(color_cluster) <- as.character(cluster_unique)}
    if(!is.null(color_vector)){
      color_cluster <- rep(0,length(unique(cluster_vivo)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- color_vector[i]
      }
      names(color_cluster) <- as.character(cluster_unique)
    }
    batch_vivo <- rep("In vivo dataset", length(cluster_vivo))
    batch_vitro <- rep("In vitro dataset", length(cluster_vitro))
    batch_all <- c(batch_vivo, batch_vitro)
    data_heatmap=data.frame(Cluster = cluster_all, Condition = batch_all)
    haCol1 <- ComplexHeatmap::HeatmapAnnotation(df= data_heatmap,col=list(Cluster=color_cluster,Condition = c("In vivo dataset" = "green 1", "In vitro dataset" = "yellow 1" )), show_legend = T)
    ht21 <- ComplexHeatmap::Heatmap(as.matrix(heatdata)
                                    , col= circlize::colorRamp2(c(0, round(max(heatdata),4)), c("grey", "brown"))
                                    , name = "log norm counts"
                                    #, column_title = "Absolute values"
                                    #, cluster_columns = as.dendrogram(my.tree)
                                    ,cluster_columns = F
                                    , cluster_rows =F
                                    , top_annotation = haCol1
                                    , row_names_gp = grid::gpar(fontsize = 6)
                                    , show_column_names = F
                                    , show_row_names = T
    )
    ComplexHeatmap::draw(ht21)
  }
  if (mean_profile == TRUE){
    if (sum(!unlist(lapply(list(c("ComplexHeatmap","circlize")),requireNamespace, quietly = TRUE)))>0) {
      stop("Package circlize and ComplexHeatmap needed for this function to work. Please install them: install.packages('circlize') and BiocManager::install('ComplexHeatmap')")
    }
    if (sum(marker_plot %in% row.names(norm_vitro)) == 0){
      stop("No provided genes are in the dataset norm_vitro. Please change marker_plot")
    }
    if (sum(marker_plot %in% row.names(norm_vivo)) == 0){
      stop("No provided genes are in the dataset norm_vivo. Please change marker_plot")
    }
    marker_plot <- marker_plot[marker_plot %in% row.names(norm_vitro)]
    marker_plot <- marker_plot[marker_plot %in% row.names(norm_vivo)]
    selected_stages <- levels(factor(cluster_vitro))
    median_es <- find_mean_cluster(cluster_vitro, selected_stages, norm_vitro, marker_plot, method = "mean")
    median_es_matrix <- matrix(unlist(median_es), ncol = length(median_es), byrow = FALSE)
    colnames(median_es_matrix) <- selected_stages
    row.names(median_es_matrix) <- marker_plot

    selected_stages_vivo <- levels(factor(cluster_vivo))
    median_vivo <- find_mean_cluster(cluster_vivo, selected_stages_vivo, norm_vivo, marker_plot, method = "mean")
    median_vivo_matrix <- matrix(unlist(median_vivo), ncol = length(median_vivo), byrow = FALSE)
    colnames(median_vivo_matrix) <- selected_stages_vivo
    row.names(median_vivo_matrix) <- marker_plot
    heatdata <- cbind(median_vivo_matrix, median_es_matrix)
    row.names(heatdata) <- marker_plot_name
    cluster_all <- c(selected_stages_vivo, selected_stages)
    cluster_unique <- unique(cluster_all)
    if(is.null(color_vector)){
      color_cluster <- rep(0,length(unique(cluster_all)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- gg_color_hue(length(color_cluster))[i]
      }
      names(color_cluster)=as.character(cluster_unique)}
    if(!is.null(color_vector)){
      color_cluster=rep(0,length(unique(cluster_all)))
      for (i in 1:length(color_cluster) ){
        color_cluster[i] <- color_vector[i]
      }
      names(color_cluster) <- as.character(cluster_unique)
    }
    batch_vivo <- rep("In vivo dataset", dim(median_vivo_matrix)[2])
    batch_vitro <- rep("In vitro dataset", dim(median_es_matrix)[2])
    batch_all <- c(batch_vivo, batch_vitro)
    data_heatmap=data.frame(Cluster = cluster_all, Condition = batch_all)
    haCol1 <- ComplexHeatmap::HeatmapAnnotation(df= data_heatmap,col=list(Cluster=color_cluster,Condition = c("In vivo dataset" = "green 1", "In vitro dataset" = "yellow 1")), show_legend = T)
    ht21 <- ComplexHeatmap::Heatmap(as.matrix(heatdata)
                                    , col= circlize::colorRamp2(c(0, round(max(heatdata),4)), c("grey", "brown"))
                                    , name = "log norm counts"
                                    #, column_title = "Absolute values"
                                    #, cluster_columns = as.dendrogram(my.tree)
                                    ,cluster_columns = F
                                    , cluster_rows =F
                                    , top_annotation = haCol1
                                    , row_names_gp = grid::gpar(fontsize = 6)
                                    , show_column_names = F
                                    , show_row_names = T
    )
    ComplexHeatmap::draw(ht21)
  }
  if (!is.logical(mean_profile)) {
    stop("Parameter mean_profile must be logical")
  }
}


#' plot_umap
#' @param coordinate_umap Data frame with dimensionality reduction coordinates.
#' Number of rows must be equal to the number of cells
#' @param cluster Vector with cluster assignment. The length must be qual to the number of cells
#' @return ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://CRAN.R-project.org/package=ggplot2}
#' @export plot_umap



plot_umap <- function(coordinate_umap, cluster){
  umap_plot <- ggplot2::ggplot(coordinate_umap, ggplot2::aes(coordinate_umap[, 1], coordinate_umap[, 2] )) +
    ggplot2::geom_point(ggplot2::aes(colour = as.factor(cluster))) + ggplot2::xlab("UMAP 1") + ggplot2::ylab("UMAP 2") + ggplot2::labs(col = "Cluster")
  return(list(umap_plot))
}

#' plot_gene
#'
#' Cells are coloured according to the expression of \emph{gene_id} and plotted
#' according to \emph{coordinate_umap}.
#'
#' @param gene_id Character name of the gene.
#' @param title_name Character name.
#' @param norm_matrix Data frame with cells on the columns and genes on the rows
#' @inheritParams plot_umap
#' @return ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://CRAN.R-project.org/package=ggplot2}
#'
#' @export plot_gene


plot_gene <- function(norm_matrix, coordinate_umap, gene_id, title_name){


  gene_expr <- as.vector(norm_matrix[gene_id, ])
  index_sort <- order(gene_expr)
  row.names(coordinate_umap) <- colnames(norm_matrix)
  coordinate_umap <- coordinate_umap[colnames(norm_matrix)[index_sort], ]
  gene_expr <- sort(gene_expr)
  umap_plot <- ggplot2::ggplot(coordinate_umap, ggplot2::aes(coordinate_umap[, 1], coordinate_umap[, 2] )) +
    ggplot2::geom_point(ggplot2::aes(colour = (gene_expr)))  + ggplot2::xlab("UMAP 1") + ggplot2::ylab("UMAP 2") + ggplot2::labs(col = "Log norm counts")+
    ggplot2::scale_colour_gradient(low = "white", high = "blue") +
    ggplot2::ggtitle(title_name)
  return((umap_plot))
}



#' plot_qc_cluster
#'
#'
#' @param seurat_object Output of the function \emph{create_seurat_object}
#' @param mit_prefix Character name. Starting letters to identify mitochondrial genes
#' @param rib_prefix_1 Character name. Starting letters to identify ribosomial genes
#' @param rib_prefix_2 Character name. Starting letters to identify ribosomial genes
#' @inheritParams SCOPRO
#' @return ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://CRAN.R-project.org/package=ggplot2}
#'
#' @export plot_qc_cluster



plot_qc_cluster <- function(seurat_object,cluster_vitro,mit_prefix,rib_prefix_1,rib_prefix_2){

  raw_counts <- (Seurat::GetAssayData(seurat_object, slot = "counts",assay="RNA"))
  sum_umi <- apply(raw_counts,2,sum)
  sum_geni <- apply(raw_counts,2,function(x){
    x <- x[x!=0]
    return(length(x))
  })


  frac_mito <- apply(raw_counts,2,function(x){
    return(sum(x[grep(mit_prefix,row.names(raw_counts))])/sum(x))



  })


  frac_ribo <- apply(raw_counts,2,function(x){
    frac_ribo_s <- x[grep(rib_prefix_1, row.names(raw_counts))]
    frac_ribo_l <- x[grep(rib_prefix_2, row.names(raw_counts))]
    frac_ribo_all <- c(frac_ribo_s, frac_ribo_l)
    return(sum(frac_ribo_all) / sum(x))
  })


  cluster_vitro <- as.factor(cluster_vitro)

  data_plot <- data.frame(as.numeric(frac_mito), as.numeric(sum_geni), as.numeric(sum_umi), as.numeric(frac_ribo),cluster_vitro)


  plot_1 <- ggplot2::ggplot(data_plot, ggplot2::aes(x=cluster_vitro, y=frac_mito,fill=cluster_vitro)) + ggplot2::geom_boxplot() + ggplot2::xlab("clusters") + ggplot2::ylab("frac mito reads") + ggplot2::ggtitle("Frac mito reads")


  plot_2 <- ggplot2::ggplot(data_plot, ggplot2::aes(x=cluster_vitro, y=sum_umi,fill=cluster_vitro)) +
    ggplot2::geom_boxplot() + ggplot2::xlab("clusters") + ggplot2::ylab("UMI counts") + ggplot2::ggtitle("UMI counts")

  plot_3 <- ggplot2::ggplot(data_plot, ggplot2::aes(x=cluster_vitro, y=sum_geni,fill=cluster_vitro)) +
    ggplot2::geom_boxplot() + ggplot2::xlab("clusters") + ggplot2::ylab("Number of genes") + ggplot2::labs(color="Cluster") + ggplot2::ggtitle("Number of genes")

  plot_4 <- ggplot2::ggplot(data_plot, ggplot2::aes(x=cluster_vitro, y=frac_ribo, fill=cluster_vitro)) +
    ggplot2::geom_boxplot() + ggplot2::xlab("clusters") + ggplot2::ylab("Frac ribo reads") + ggplot2::labs(color="Cluster") + ggplot2::ggtitle("Frac ribo reads")
  return (list(plot_1, plot_2, plot_3, plot_4))
}


#' plot_barplot
#'
#'
#' @param fraction Logical. If TRUE, for each \emph{cluster_vitro}, the fraction of cells
#' mapped to the in vivo stages are shown. If FALSE, the number of cells mapped to the in vivo stages are shown
#' @inheritParams SCOPRO
#' @inheritParams findCellTypesSeurat
#' @inheritParams findCellTypesScibet
#' @inheritParams findCellTypes_scmap
#' @inheritParams cellTypesPerClusterBalloonPlot
#' @param obj Seurat object given as output from functions \emph{findCellTypesSeurat}, \emph{findCellTypesScibet} or
#' \emph{findCellTypes_scmap}
#' @return ggplot2 object.
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://CRAN.R-project.org/package=ggplot2}
#'
#' @export plot_barplot


plot_barplot <- function(obj, method, cluster_vitro, threshold, fraction = TRUE){
  if (method == "Seurat") {

    obj@meta.data$predicted.id[obj@meta.data$prediction.score.max <= threshold] <- "Unassigned"
    df=data.frame(obj@meta.data$predicted.id,cluster_vitro)
    colnames(df)=c("projection_result","cluster_vitro")
    p = create_barplot(df[,1],df[,2],paste0("Projection result ", method),fraction = fraction)
    return(p[[1]]+ggplot2::theme(text = ggplot2::element_text(size=14),axis.text.x = ggplot2::element_text(size=14,angle = 45, vjust = 1, hjust = 1)))
  }
  if (method == "scibetR") {
    df=data.frame(obj@meta.data$predictions,cluster_vitro)
    colnames(df)=c("projection_result","cluster_vitro")
    p = create_barplot(df[,1],df[,2],paste0("Projection result ", method),fraction = fraction)
    return(p[[1]]+ggplot2::theme(text = ggplot2::element_text(size=14),axis.text.x = ggplot2::element_text(size=14,angle = 45, vjust = 1, hjust = 1)))
  }
  if (method == "scmap") {
    df=data.frame(obj@meta.data$predictions,cluster_vitro)
    colnames(df)=c("projection_result","cluster_vitro")
    p = create_barplot(df[,1],df[,2],paste0("Projection result ", method),fraction = fraction)
    return(p[[1]]+ggplot2::theme(text = ggplot2::element_text(size=14),axis.text.x = ggplot2::element_text(size=14,angle = 45, vjust = 1, hjust = 1)))
  }
  if (method != "Seurat" & method != "scibetR" & method !=
      "scmap") {
    stop("method must be equal to Seurat, ScibetR or scmap.")
  }
}



#' difference_significance_vitro

#' @inheritParams SCOPRO
#' @inheritParams findCellTypesSeurat
#' @inheritParams findCellTypesScibet
#' @inheritParams findCellTypes_scmap
#' @inheritParams cellTypesPerClusterBalloonPlot
#' @param p_value Numeric value.
#' @return List with empirical p values
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @description For each cluster, a centroid is defined as the mean expression value of \emph{marker_stages_filter} across cells in the cluster.
#' Comparison of in vitro and in vivo is done starting from a distance matrix based on spearman correlation between cluster centroid.
#' For a cluster A in the in vitro dataset, the Spearman’s correlation based distance between the cells and its centroid is computed. These values define an empirical distribution. The Spearman’s correlation based distance between the centroid of cluster A in the in vitro dataset and centroid of cluster B in the in vivo dataset is also computed.  An empirical p value is assigned to cluster B, from the comparison with the empirical distribution.
#' The distance between a pair of clusters is significant if the empirical p value above defined is smaller or equal to \emph{p_value}

#'
#' @export difference_significance_vitro





difference_significance_vitro <- function(norm_vitro, p_value = 0.10, marker_stages, marker_stages_filter, cluster_vitro, cluster_vivo, selected_stages, norm_vivo){
  for (i in selected_stages){
    index_specific <- which(selected_stages %in% i)
    marker_specific <- marker_stages[[index_specific]]
    if (length(marker_specific) == 0) {
      stop(paste0("There are no markers for the stage: ", i))
    }
  }
  marker_stages_top <- rep(list(0),length(selected_stages))
  for (i in 1:length(selected_stages)){
    marker_stages_top[[i]] = marker_stages[[i]]
  }
  marker_all <- unlist(marker_stages_top)
  marker_all <- marker_all[marker_all%in%marker_stages_filter]
  genes_all <- marker_all

  mean_profile_vivo <- find_mean_cluster(cluster_vivo, selected_stages, norm_vivo, genes_all, method = "mean")

  selected_stages_vitro = levels(factor(cluster_vitro))
  mean_profile_vitro <- find_mean_cluster(cluster_vitro, selected_stages_vitro, norm_vitro, genes_all, method = "mean")


  name_vitro <- levels(factor(cluster_vitro))
  difference_vitro <- rep(list(0),length(name_vitro))
  for ( i in 1:length(name_vitro)){
    difference_vitro[[i]] <- difference_significance_vitro_single(norm_vitro,mean_profile_vitro, mean_profile_vivo, p_value = 0.10, name_vitro = name_vitro[i], genes_all, cluster_vitro, selected_stages)
    if (sum(difference_vitro[[i]] < p_value)){
      logic_index = which(difference_vitro[[i]] > p_value)
      cluster_vivo_different = selected_stages[logic_index]
      print(paste0("In vitro cluster ", name_vitro[i]," is different from in vivo clusters", cluster_vivo_different))
    }
    else{
      print(paste0("In vitro cluster ", name_vitro[i]," is not different from all in vivo clusters"))
    }
    names(difference_vitro[[i]]) = selected_stages
  }
  names(difference_vitro) = name_vitro
  return(difference_vitro)
}




#' plot_distance_vivo_vitro
#' @inheritParams SCOPRO
#' @inheritParams plot_score_genes
#' @inheritParams plot_in_vivo_markers
#' @inheritParams heatmap_markers_vivo
#' @inheritParams heatmap_markers_vitro
#' @inheritParams plot_score_genes_heatmap
#' @param text_size Integer number specifyng the size of the text.
#' @return Heatmap class object
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
#'
#'
#'
#' @export plot_distance_vivo_vitro






plot_distance_vivo_vitro <- function(norm_vivo, cluster_vivo, selected_stages, marker_stages, marker_stages_filter, norm_vitro, cluster_vitro, text_size = 10){
  if (sum(!unlist(lapply(list(c("ComplexHeatmap","circlize")),requireNamespace, quietly = TRUE)))>0) {
    stop("Package circlize and ComplexHeatmap needed for this function to work. Please install them: install.packages('circlize') and BiocManager::install('ComplexHeatmap')")
  }
  for (i in selected_stages){
    index_specific <- which(selected_stages %in% i)
    marker_specific <- marker_stages[[index_specific]]
    if (length(marker_specific) == 0) {
      stop(paste0("There are no markers for the stage: ", i))
    }
  }
  marker_stages_top <- rep(list(0),length(selected_stages))
  for (i in 1:length(selected_stages)){
    marker_stages_top[[i]] = marker_stages[[i]]
  }
  marker_all <- unlist(marker_stages_top)

  marker_all <- marker_all[marker_all%in%marker_stages_filter]
  genes <- marker_all
  mean_profile_vivo <- find_mean_cluster(cluster_vivo, selected_stages, norm_vivo, genes, method = "mean")

  selected_stages_vitro = levels(factor(cluster_vitro))
  mean_profile_vitro <- find_mean_cluster(cluster_vitro, selected_stages_vitro, norm_vitro, genes, method = "mean")


  median_vivo <- matrix(unlist(mean_profile_vivo), ncol = length(mean_profile_vivo), byrow = FALSE)
  median_vitro <- matrix(unlist(mean_profile_vitro), ncol = length(mean_profile_vitro), byrow = FALSE)

  median_all_matrix=cbind(median_vivo, median_vitro)
  names_vivo=rep(0,length(mean_profile_vivo))
  for ( i in 1:length(mean_profile_vivo)){
    j=selected_stages[i]
    names_vivo[i]=paste0("In vivo-",j)}

  names_vitro=rep(0,length(mean_profile_vitro))
  for ( i in 1:length(mean_profile_vitro)){
    j=selected_stages_vitro[i]
    names_vitro[i]=paste0("In vitro-",j)}

  colnames(median_all_matrix)=c(names_vivo,names_vitro)



  median_all_matrix<-as.matrix(median_all_matrix)
  dissimilarity <- sqrt((1-cor((median_all_matrix), method = "spearman"))/2)

  my.dist <- as.dist(dissimilarity)


  dist_mouse=as.matrix(my.dist)
  dist_mouse=dist_mouse[grepl("In vivo-",row.names(dist_mouse)),grepl("In vitro-",row.names(dist_mouse))]




  ht21 <- ComplexHeatmap::Heatmap(as.matrix(dist_mouse),
                                  cluster_rows = TRUE, col = circlize::colorRamp2(c(round(min(dist_mouse),2),
                                                                                    round(max(dist_mouse),2)), c("grey", "brown")),
                                  name = "spearman distance", cluster_columns = TRUE,
                                  row_names_gp = grid::gpar(fontsize = text_size), show_column_names = TRUE,
                                  show_row_names = TRUE)
  ComplexHeatmap::draw(ht21)


}



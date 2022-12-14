% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_score_genes_heatmap}
\alias{plot_score_genes_heatmap}
\title{plot_score_genes_heatmap}
\usage{
plot_score_genes_heatmap(
  norm_vivo,
  cluster_vivo,
  norm_vitro,
  cluster_vitro,
  marker_plot,
  marker_plot_name,
  mean_profile = FALSE,
  color_vector = NULL
)
}
\arguments{
\item{norm_vivo}{Norm count matrix (n_genes X n_cells) for in vivo dataset}

\item{cluster_vivo}{cluster for in vivo dataset}

\item{norm_vitro}{Norm count matrix (n_genes X n_cells) for in vitro dataset}

\item{cluster_vitro}{cluster for in vitro dataset}

\item{marker_plot}{Character vector with the names of the genes to plot.}

\item{marker_plot_name}{Character vector with the names of the rows to show in the heatmap}

\item{mean_profile}{Logical value. If TRUE, then the mean expression profile for\emph{marker_plot} in each cluster is shown}

\item{color_vector}{Character vector with colour assignment.}
}
\value{
Heatmap class object
}
\description{
plot_score_genes_heatmap
}
\seealso{
\url{https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap}
}
\author{
Gabriele Lubatti \email{gabriele.lubatti@helmholtz-muenchen.de}
}

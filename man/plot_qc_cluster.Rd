% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_qc_cluster}
\alias{plot_qc_cluster}
\title{plot_qc_cluster}
\usage{
plot_qc_cluster(
  seurat_object,
  cluster_vitro,
  mit_prefix,
  rib_prefix_1,
  rib_prefix_2
)
}
\arguments{
\item{seurat_object}{Output of the function \emph{create_seurat_object}}

\item{cluster_vitro}{cluster for in vitro dataset}

\item{mit_prefix}{Character name. Starting letters to identify mitochondrial genes}

\item{rib_prefix_1}{Character name. Starting letters to identify ribosomal genes}

\item{rib_prefix_2}{Character name. Starting letters to identify ribosomal genes}
}
\value{
ggplot2 object.
}
\description{
plot_qc_cluster
}
\seealso{
\url{https://CRAN.R-project.org/package=ggplot2}
}
\author{
Gabriele Lubatti \email{gabriele.lubatti@helmholtz-muenchen.de}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/projection_published_methods.R
\name{findCellTypes_scmap}
\alias{findCellTypes_scmap}
\title{findCellTypes_scmap}
\usage{
findCellTypes_scmap(
  queryObj,
  referenceObj,
  namedLabels,
  n_features = 500,
  threshold = 0.7
)
}
\arguments{
\item{queryObj}{Seurat object of the in vitro dataset}

\item{referenceObj}{Seurat object of the in vivo dataset}

\item{namedLabels}{refdata parameter in the function \emph{FindTransferAnchors}}

\item{n_features}{selectFeatures parameter in the function \emph{selectFeatures}}

\item{threshold}{threshold parameter in the function \emph{scmapCluster}}
}
\value{
Seurat object
}
\description{
Projection done by function \emph{scmapCluster} from package \emph{scmap}
}
\seealso{
\url{https://bioconductor.org/packages/release/bioc/html/scmap.html}
}
\author{
Gabriele Lubatti \email{gabriele.lubatti@helmholtz-muenchen.de}
}

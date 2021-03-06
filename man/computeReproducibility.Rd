% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{computeReproducibility}
\alias{computeReproducibility}
\title{Computes reproducibility}
\usage{
computeReproducibility(
  se,
  samples,
  features = NULL,
  log = TRUE,
  pseudocount = 1,
  pt.alpha = 0.5,
  pt.size = 0.9,
  contour = TRUE,
  density = FALSE,
  contour.alpha = 0.5,
  contour.size = 0.7,
  density.alpha = 0.5,
  spline.degree = 6
)
}
\arguments{
\item{se}{summarized experiment}

\item{samples}{samples to calculate reproducibility}

\item{features}{features to calculate metric for}

\item{log}{TRUE if data should be logged}

\item{pseudocount}{pseudocount for log}

\item{pt.alpha}{alpha for points}

\item{pt.size}{size for points}

\item{contour}{TRUE if contours should be displayed}

\item{density}{TRUE if density should be displayed}

\item{contour.alpha}{alpha for contour}

\item{contour.size}{size for contour}

\item{density.alpha}{alpha for density}

\item{spline.degree}{degree of spline fit}
}
\value{
a list object with plots and metrics
}
\description{
Computes reproducibility
}

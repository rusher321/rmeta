% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/figureMeta.R
\name{topTax}
\alias{topTax}
\title{topTax
plot the high abundance tax}
\usage{
topTax(
  metadata,
  K = 20,
  rmUnclass = F,
  sample_order = NULL,
  tax_order = NULL,
  tax_colour = NULL
)
}
\arguments{
\item{metadata}{, tax profile row is tax col is sample}

\item{K}{, top number}

\item{rmUnclass}{, if to remove the unclass tax}

\item{sample_order, }{vector of sample id}

\item{tax_order, }{vector of taxnomy}

\item{tax_colour, }{vector of colour}
}
\value{
figure
}
\description{
topTax
plot the high abundance tax
}

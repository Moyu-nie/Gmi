% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intertomain.R
\name{intertomain}
\alias{intertomain}
\title{Convert Interaction Term Names to Corresponding Main Effect Indices}
\usage{
intertomain(interNameList, p)
}
\arguments{
\item{interNameList}{A list of interaction term names.}

\item{p}{Integer. The total number of variables (dimensionality) in the design matrix.}
}
\value{
An integer vector giving the indices of main effects involved in the interaction terms.
}
\description{
convert an interaction name list into the corresponding main effect for enforcingstrong hierachy
}
\examples{
interNameList <- list("X1X3", "X2X4", "X1X2")
p <- 5
intertomain(interNameList, p)
}

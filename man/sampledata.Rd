\name{sampledata}
\alias{sampledata}
\docType{data}
\title{
Sample Data
}
\description{
An example for the data input.
}
\usage{data("sampledata")}
\format{
  A data frame with 74 observations on the following 10 variables.
  \describe{
    \item{\code{sig_type}}{a character vector}
    \item{\code{qtl_type}}{a character vector}
    \item{\code{cg}}{a character vector}
    \item{\code{rs}}{a character vector}
    \item{\code{rs.CHR}}{a character vector}
    \item{\code{rs.POS}}{a numeric vector}
    \item{\code{rs.gene}}{a character vector}
    \item{\code{cg.CHR}}{a numeric vector}
    \item{\code{cg.POS}}{a numeric vector}
    \item{\code{cg.gene}}{a character vector}
  }
}
\details{
The input data must contain cols that sampledata has
}
\source{
nothing
}
\references{
2016 tutorial
}
\examples{
data(sampledata)
## maybe str(sampledata) ; plot(sampledata) ...
}
\keyword{datasets}

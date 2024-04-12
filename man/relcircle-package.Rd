\name{relcircle-package}
\alias{relcircle-package}
\alias{relcircle}
\docType{package}
\title{
Draw Regulatory Relationships Between Genes
}
\description{
Connect relevant sites and label their relationships.
}
\details{

Relationship table, needs to include "rs","rs.CHR","rs.gene","rs.POS","cg","cg.CHR","cg.gene","cg.POS".
}
\author{
Ruilin Li

Maintainer:Ruilin Li<2837471734@qq.com>
}
\references{
2016 tutorial
}
\keyword{ package }
\seealso{
\code{\link{plot}}
}
\examples{
data(sampledata)
relacircle(sampledata)
}

\name{relacircle}
\alias{relacircle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw Regulatory Relationships Between Genes
}
\description{
Connect relevant sites and label their relationships.
}
\usage{
relacircle(table, type1 = FALSE, type2 = FALSE, line.col = FALSE,
pch = 1, pch.col = "blue", lty = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
Relationship table, needs to include "rs","rs.CHR","rs.gene","rs.POS","cg","cg.CHR","cg.gene","cg.POS".
}
  \item{type1}{
Classification of corresponding relationships, displayed by different colors, optional, default gray.
}
  \item{type2}{
The second classification is displayed according to different line types, which can be left blank and defaults to solid lines.
}
  \item{line.col}{
The color, vector, and length of the line connecting the points need to be consistent with the number of types in type1. If type1 is not set, enter one, and default to gray.
}
  \item{pch}{
The type of annotation site on the chromosome, default to a hollow circle.
}
  \item{pch.col}{
The color of the upper chromosomal site, default to blue.
}
  \item{lty}{
The type, array, and length of the line connecting the points need to be consistent with the number of types in type2. If type2 is not set, enter one. The default is a solid line.
}
}
\details{
Relationship table, needs to include "rs","rs.CHR","rs.gene","rs.POS","cg","cg.CHR","cg.gene","cg.POS".
}
\value{
return the graphic.
}
\references{
2016 tutorial
}
\author{
Ruilin Li
}
\note{
return the graphic.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot}}
}
\examples{
data(sampledata)
relacircle(sampledata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):

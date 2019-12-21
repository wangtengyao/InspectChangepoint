\name{plot.hdchangeseq}
\alias{plot.hdchangeseq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot function for 'hdchangeseq' class
}
\description{
Visualising the high-dimensional time series in an 'hdchangeseq' class object. The data matrix or its mean structure is visualised using a grid of coloured rectangles with colours corresponding to the value contained in corresponding coordinates. A heat-spectrum (red to white for values low to high) is used to convert values to colours. 
}
\usage{
\method{plot}{hdchangeseq}(x, noise = TRUE, shuffle = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of 'hdchangeseq' class
}
  \item{noise}{
If noise == TRUE, the data matrix is plotted, otherwise, only the mean structure is plotted.
}
  \item{shuffle}{
Whether to shuffle the rows of the plotted matrix.
}
  \item{\dots}{
Other graphical parameters are not used.}
}
\note{
Based on the \code{\link{image}} function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
n <- 2000; p <- 200; ks <- 40; zs <- c(500,1000,1500); varthetas <- c(0.1,0.15,0.2); overlap <- 0.5
obj <- multi.change(n, p, ks, zs, varthetas, overlap)
plot(obj, noise = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~aplot }% use one of  RShowDoc("KEYWORDS")

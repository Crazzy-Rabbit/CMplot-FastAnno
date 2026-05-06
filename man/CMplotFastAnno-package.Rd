\name{CMplotFastAnno-package}
\alias{CMplotFastAnno-package}
\alias{CMplotFastAnno}
\docType{package}
\title{CMplot-Compatible GWAS Visualization with Fast Annotation Support}
\description{
CMplotFastAnno provides a CMplot-compatible \code{CMplot()} function for
Manhattan plots, circular Manhattan plots, combined Manhattan and Q-Q plots,
Q-Q plots, and marker density plots.
}
\details{
The package keeps the original CMplot calling style and adds practical GWAS
visualization features, including external annotation tables, top and nearby
label placement, multi-lane annotation, already-scaled \code{-log10(P)} input,
display-only y-axis compression with \code{skip}, \code{cut}, and
\code{cutfactor}, and clearer significance-line controls.

Only base R, graphics, grDevices, stats, and utils are required.
}
\seealso{
\code{\link{CMplot}}
}
\keyword{package}

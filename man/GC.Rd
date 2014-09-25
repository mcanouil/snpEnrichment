\name{GC}
\alias{GC}

\title{
Full Garbage Collection
}

\description{
\code{\link{GC}} performs garbage collection until free memory indicators show no change.
}

\usage{
GC(verbose = getOption("verbose"), reset=FALSE)
}

\arguments{
    \item{verbose}{[logical]: if \code{TRUE}, the garbage collection prints statistics about cons cells and the space allocated for vectors.}
    \item{reset}{[logical]: if \code{TRUE} the values for maximum space used are reset to the current values.}
}

% \details{}

\value{
\code{\link{GC}} returns a matrix with rows \code{"Ncells"} (_cons cells_), usually
28 bytes each on 32-bit systems and 56 bytes on 64-bit systems,
and \code{"Vcells"} (_vector cells_, 8 bytes each), and columns
\code{"used"} and \code{"gc trigger"}, each also interpreted in megabytes
(rounded up to the next 0.1Mb).

If maxima have been set for either \code{"Ncells"} or \code{"Vcells"}, a
fifth column is printed giving the current limits in Mb (with \code{NA}
denoting no limit).

The final two columns show the maximum space used since the last
call to \code{GC(reset=TRUE)} (or since R started).
}

% \note{}

\seealso{
The \code{R Internals} manual.

\code{Memory} on R's memory management, and \code{gctorture} if you are an R
developer.

\code{reg.finalizer} for actions to happen at garbage collection.
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{GC}
\keyword{garbage}

\examples{
GC() # - do it now
x <- integer(100000); for(i in 1:18) x <- c(x,i)
GC(TRUE)
GC(reset=TRUE)
}

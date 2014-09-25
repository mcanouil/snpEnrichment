\name{mclapply2}
\alias{mclapply2}

\title{
Parallel Versions of \code{lapply} with cores and memory control
}

\description{
\code{\link{mclapply2}} is a mclapply modification from parallel package, to avoid a memory overload.
The maximum number of cores is computed depending on the amount of memory used by the parent process and the amount of free memory available on the machine.
Note: This number is under-estimated.
}

\usage{
mclapply2(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE,
          mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
          mc.cleanup = TRUE, mc.allow.recursive = TRUE)
}

\arguments{
    \item{X}{a vector (atomic or list) or an expressions vector.
        Other objects (including classed objects) will be coerced by \code{\link{as.list}}.}
    \item{FUN}{the function to be applied to (\code{mclapply}) each
        element of \code{X} or (\code{mcmapply}) in parallel to \code{\dots}.}
    \item{...}{For \code{mclapply}, optional arguments to \code{FUN}.}
    \item{mc.preschedule}{if set to \code{TRUE} then the computation is
        first divided to (at most) as many jobs are there are cores and then
        the jobs are started, each job possibly covering more than one
        value.  If set to \code{FALSE} then one job is forked for each value
        of \code{X}.  The former is better for short computations or large
        number of values in \code{X}, the latter is better for jobs that
        have high variance of completion time and not too many values of
        \code{X} compared to \code{mc.cores}.}
    \item{mc.set.seed}{See \code{mcparallel}.}
    \item{mc.silent}{if set to \code{TRUE} then all output on
        \file{stdout} will be suppressed for all parallel processes forked
        (\file{stderr} is not affected).}
    \item{mc.cores}{The number of cores to use, i.e. at most how many
        child processes will be run simultaneously.   The option is
        initialized from environment variable \env{MC_CORES} if set.  Must
        be at least one, and parallelization requires at least two cores.}
    \item{mc.cleanup}{if set to \code{TRUE} then all children that have
        been forked by this function will be killed (by sending
        \code{SIGTERM}) before this function returns.  Under normal
        circumstances \code{mclapply} waits for the children to deliver
        results, so this option usually has only effect when \code{mclapply}
        is interrupted. If set to \code{FALSE} then child processes are
        collected, but not forcefully terminated.  As a special case this
        argument can be set to the number of the signal that should be used
        to kill the children instead of \code{SIGTERM}.}
    \item{mc.allow.recursive}{Unless true, calling \code{mclapply} in a
        child process will use the child and not fork again.}
}


% \details{}

\value{
A list of the same length as \code{X} and named by \code{X}.
}

% \note{}

\seealso{
\code{mclapply} from package parallel.
}

\author{Mickael Canouil \email{mickael.canouil@good.ibl.fr}}

\keyword{mclapply2}
\keyword{parallel}
\keyword{core}

% \examples{}

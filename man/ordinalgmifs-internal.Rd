\name{ordinalgmifs-internal}
\alias{ordinalgmifs-internal}
\alias{du.adjcat}
\alias{du.bcr}
\alias{du.cum}
\alias{du.fcr}
\alias{du.stereo}
\alias{fn.acat}
\alias{fn.bcr}
\alias{fn.cum}
\alias{fn.fcr}
\alias{fn.stereo}
\alias{G}
\alias{gradient}

\title{
Functions Called by ordinalgmifs Functions, Not by the User
}
\description{
These functions are called my other \code{ordinalgmifs} functions and are not intended to be directly called by the user.
}
\details{
The \code{du.adjcat}, \code{du.bcr}, \code{du.cum}, \code{du.fcr}, and \code{du.stereo} functions calculate the derivatives at the current step for the adjacent category, backward CR, cumulative link, forward CR, and stereotype logit models, respectively, are used to identify which penalized parameter is updated. The \code{fn.acat}, \code{fn.bcr}, \code{fn.cum}, \code{fn.fcr}, and \code{fn.stereo} are the log-likelihood functions for the adjacent category, backward CR, cumulative link, forward CR, and stereotype logit models, respectively, are used to estimate the thresholds and non-penalized subset parameters (if included) at each step of the algorithm.
The \code{G} function returns the probability for the indicated link function. The \code{gradient} function returns the gradient of the log-likelihood for the cumulative link models and is used for the cumulative link constrained optimization.
}
\value{these functions are called for intermediate results used 
internally by user-invoked functions}
\author{
Kellie J. Archer, \email{archer.43@osu.edu}
}
\seealso{
See Also as \code{ \link{ordinalgmifs} }
}
\keyword{ internal }

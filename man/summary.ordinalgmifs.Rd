\name{summary.ordinalgmifs}
\alias{summary.ordinalgmifs}
\title{
Summarize an Ordinal GMIFS Object.
}
\description{
\code{summary} method for class \code{ordinalgmifs}.}
\usage{
\method{summary}{ordinalgmifs}(object, model.select = "AIC", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{ordinalgmifs} object.
}
  \item{model.select}{
when \code{x} is specified any model along the solution path can be selected. The default is \code{model.select="AIC"} which extracts the model having the lowest AIC. Other options are \code{model.select="BIC"} or any numeric value from the solution path.
}
  \item{\dots}{other arguments.}
}
\details{
Prints the following items extracted from the fitted \code{ordinalgmifs} object:
the probability model and link used and model parameter estimates. For models that include
\code{x}, the parameter estimates, AIC, BIC, and log-likelihood are printed for indicated \code{model.select} step or if \code{model.select} is not supplied the step at which the minimum AIC was observed.
}

\author{
Kellie J. Archer
}
\seealso{
See Also \code{\link{ordinalgmifs}}, \code{\link{coef.ordinalgmifs}}, \code{\link{plot.ordinalgmifs}}, \code{\link{predict.ordinalgmifs}}
}
%\examples{
%\donttest{
%data(hccframe)
%# To minimize processing time, MPO_E302_R is coerced into the model and only a subset of 
%# two CpG sites (DDIT3_P1313_R and HDAC9_P137_R) are included as penalized covariates
%# in this demonstration.
%hcc.fit <- ordinalgmifs(group ~ MPO_E302_R, x = c("DDIT3_P1313_R", "HDAC9_P137_R"), data = hccframe)
%summary(hcc.fit)
%}
%}
\keyword{ methods }

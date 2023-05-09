\name{plot.ordinalgmifs}
\alias{plot.ordinalgmifs}
\title{
Plot Solution Path for Ordinal GMIFS Fitted Model.
}
\description{
This function plots either the coefficient path, the AIC, or the log-likelihood for a fitted \code{ordinalgmifs} object.
}
\usage{
\method{plot}{ordinalgmifs}(x, type = "trace", xlab=NULL, ylab=NULL, main=NULL, ...)
}
\arguments{
  \item{x}{
an \code{ordinalgmifs} object.
}
  \item{type}{
default is \code{"trace"} which plots the coefficient path for the fitted object. Also available are \code{"AIC"}, \code{"BIC"}, and \code{"logLik"}.
}
  \item{xlab}{
 a default x-axis label will be used which can be changed by specifying a user-defined x-axis label.
}
  \item{ylab}{
 a default y-axis label will be used which can be changed by specifying a user-defined y-axis label.
}
  \item{main}{
 a default main title will be used which can be changed by specifying a user-defined main title.
}
  \item{\dots}{other arguments.}
}
\value{No return value, called for side effects}
\author{
Kellie J. Archer
}
\seealso{
See Also \code{\link{ordinalgmifs}}, \code{\link{coef.ordinalgmifs}}, \code{\link{summary.ordinalgmifs}}, \code{\link{predict.ordinalgmifs}}
}
%\examples{
%\donttest{
%data(hccframe)
%# To minimize processing time, MPO_E302_R is coerced into the model and only a subset of 
%# two CpG sites (DDIT3_P1313_R and HDAC9_P137_R) are included as penalized covariates
%# in this demonstration.
%hcc.fit <- ordinalgmifs(group ~ MPO_E302_R, x = c("DDIT3_P1313_R", "HDAC9_P137_R"), data = hccframe)
%plot(hcc.fit)
%}
%}
\keyword{ methods }

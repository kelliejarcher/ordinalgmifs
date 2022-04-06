\name{coef.ordinalgmifs}
\alias{coef.ordinalgmifs}
\title{
Extract Model Coefficients
}
\description{
\code{coef.ordinalgmifs} is a generic function which extracts the model coefficients from a fitted model object fit using \code{ordinalgmifs}
}
\usage{
\method{coef}{ordinalgmifs}(object, model.select = "AIC", ...)
}
\arguments{
  \item{object}{
an \code{ordinalgmifs} object.
}
  \item{model.select}{
when \code{x} is specified any model along the solution path can be selected. The default is \code{model.select="AIC"} which extracts the coefficients from the model having the lowest AIC. Other options are \code{model.select="BIC"} or any numeric value from the solution path.
}
  \item{\dots}{other arguments.}
}
\value{
Coefficients extracted from the model object.
}
\references{
Hastie T., Taylor J., Tibshirani R., and Walther G. (2007) Forward stagewise regression and the monotone lasso. \emph{Electronic Journal of Statistics}, 1, 1-29.
}
\author{
Kellie J. Archer
}
\seealso{
See Also \code{\link{ordinalgmifs}}, \code{\link{summary.ordinalgmifs}}, \code{\link{plot.ordinalgmifs}}, \code{\link{predict.ordinalgmifs}}
}
%\examples{
%\donttest{
%data(hccframe)
%# To minimize processing time, MPO_E302_R is coerced into the model and only a subset of 
%# two CpG sites (DDIT3_P1313_R and HDAC9_P137_R) are included as penalized covariates
%# in this demonstration.
%hcc.fit <- ordinalgmifs(group ~ MPO_E302_R, x = c("DDIT3_P1313_R", "HDAC9_P137_R"), data = hccframe)
%coef(hcc.fit)
%}
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

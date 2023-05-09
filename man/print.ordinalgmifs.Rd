\name{print.ordinalgmifs}
\alias{print.ordinalgmifs}
\title{
Print the Contents of an Ordinal GMIFS Fitted Object.
}
\description{
This function prints the names of the list objects from an \code{ordinalgmifs} fitted model.
}
\usage{
\method{print}{ordinalgmifs}(x, ...)
}
\arguments{
  \item{x}{an \code{ordinalgmifs} object.}
  \item{\dots}{other arguments.}
}
\value{returns the object names in the fitted ordinalgmifs object}
\author{
Kellie J. Archer
}
\note{
The contents of an \code{ordinalgmifs} fitted object differ depending upon whether \code{x} is specified in the \code{ordinalgmifs} model (i.e., penalized variables are included in
the model fit hence a solution path is returned) or only \code{terms} on the right hand side of the equation are included (unpenalized variables). In the
latter case, we recommend using the VGAM package.
}
\seealso{
See Also \code{\link{ordinalgmifs}}, \code{\link{coef.ordinalgmifs}}, \code{\link{summary.ordinalgmifs}}, \code{\link{plot.ordinalgmifs}}, \code{\link{predict.ordinalgmifs}}
}
%\examples{
%\donttest{
%data(hccframe)
%# To minimize processing time, MPO_E302_R is coerced into the model and only a subset of 
%# two CpG sites (DDIT3_P1313_R and HDAC9_P137_R) are included as penalized covariates
%# in this demonstration.
%hcc.fit <- ordinalgmifs(group ~ MPO_E302_R, x = c("DDIT3_P1313_R", "HDAC9_P137_R"), data = hccframe)
%print(hcc.fit)
%}
%}
\keyword{ methods }

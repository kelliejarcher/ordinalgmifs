\name{predict.ordinalgmifs}
\alias{predict.ordinalgmifs}
\alias{fitted.ordinalgmifs}
\title{
Predicted Probabilities and Class for Ordinal GMIFS Fit.
}
\description{
This function returns a list the includes the predicted probabilities as well as the predicted class for an \code{ordinalgmifs} fitted object.
}
\usage{
\method{predict}{ordinalgmifs}(object, neww = NULL, newdata, newx = NULL, model.select = "AIC", ...)
}
\arguments{
  \item{object}{an \code{ordinalgmifs} fitted object.}
  \item{neww}{an optional formula that includes the unpenalized variables to use for predicting the response. If omitted, the training data are used.}
  \item{newdata}{an optional data.frame that minimally includes the unpenalized variables to use for predicting the response. If omitted, the training data are used.}
  \item{newx}{an optional matrix of penalized variables to use for predicting the response. If omitted, the training data are used.}
  \item{model.select}{when \code{x} is specified any model along the solution path can be selected. The default is \code{model.select="AIC"} which calculates the predicted values using the coefficients from the model having the lowest AIC. Other options are \code{model.select="BIC"} or any numeric value from the solution path.}
  \item{\dots}{other arguments.}
}
\value{
  \item{predicted}{a matrix of predicted probabilities from the fitted model.}
  \item{class}{a vector containing the predicted class taken as that class having the largest predicted probability.}
  \item{\dots}{other arguments.}
}
\author{
Kellie J. Archer, Jiayi Hou, Qing Zhou, Kyle Ferber, John G. Layne, Amanda Gentry
}
\seealso{
See Also \code{\link{ordinalgmifs}}, \code{\link{coef.ordinalgmifs}}, \code{\link{summary.ordinalgmifs}}, \code{\link{plot.ordinalgmifs}} 
}
%\examples{
%\donttest{
%data(hccframe)
%# To minimize processing time, MPO_E302_R is coerced into the model and only a subset of 
%# two CpG sites (DDIT3_P1313_R and HDAC9_P137_R) are included as penalized covariates
%# in this demonstration.
%hcc.fit <- ordinalgmifs(group~MPO_E302_R, x = c("DDIT3_P1313_R", "HDAC9_P137_R"), data = hccframe)
%phat <- predict(hcc.fit)
%table(phat$class, hccframe$group)
%head(phat$predicted)
%# See the package vignette if an independent test set is used for prediction.
%}
%}
\keyword{ methods }

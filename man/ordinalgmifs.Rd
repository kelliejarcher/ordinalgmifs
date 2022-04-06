\name{ordinalgmifs}
\alias{ordinalgmifs}
\title{
Ordinal Generalized Monotone Incremental Forward Stagewise Regression
}
\description{
This function can fit a cumulative link, adjacent category, forward and backward continuation ratio, and stereotype ordinal response model when the number of parameters exceeds the sample size, using the the generalized monotone incremental forward stagewise method.
}
\usage{
ordinalgmifs(formula, data, x = NULL, subset, epsilon = 0.001, tol = 1e-05, 
	scale = TRUE, probability.model = "Cumulative", link = "logit", 
	verbose=FALSE, assumption=NULL, ...)
}

\arguments{
  \item{formula}{an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The left side of the formula is the ordinal outcome while the variables on the right side of the formula are the covariates that are not included in the penalization process. Note that if all variables in the model are to be penalized, an intercept only model formula should be specified.}
  \item{data}{an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. }
  \item{x}{an optional matrix of predictors that are to be penalized in the model fitting process.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{epsilon}{small incremental amount used to update a coefficient at a given step.}
  \item{tol}{the iterative process stops when the difference between successive log-likelihoods is less than this specified level of tolerance.}
  \item{scale}{logical, if TRUE the penalized predictors are centered and scaled.}
  \item{probability.model}{the type of ordinal response model to be fit. Can be \code{"Cumulative"}, \code{"AdjCategory"}, \code{"ForwardCR"}, \code{"BackwardCR"}, or \code{"Stereotype"}}
  \item{link}{the link function used. Allowable links for \code{"Cumulative"}, \code{"ForwardCR"}, and \code{"BackwardCR"} are \code{"logit"}, \code{"probit"}, and \code{"cloglog"}. For an \code{"AdjCategory"} model only a \code{"loge"} link is allowed; for a \code{"Stereotype"} model only a \code{"logit"} link is allowed.}
  \item{verbose}{logical, if TRUE the step number is printed to the console (default is FALSE).}
  \item{assumption}{integer, only use with \code{probability.model = "ForwardCR"} and \code{link = "cloglog"} to denote the assumption to use for discrete censored survival modeling. If \code{assumption = 1}, assume the observation was censored at the end of the discrete time interval in which the censoring occurred; if \code{assumption = 2}, assume the observation was censored at the beginning of the interval in which censoring occurred; if \code{assumption = 3}, assume constant hazard rate within the interval in which the censoring occurred; if no censoring occurs, do not specify a value for \code{assumption}.}
  \item{\dots}{additional arguments }
}
\details{
A model specified as \code{response~terms, x=penalized.terms} where \code{response} is the ordinal response vector and \code{terms} is the series of variables in the model that are not to be penalized and \code{x} is a matrix of variables that are to be penalized. For example, \code{terms} may include the variables age and gender while \code{x} includes hundreds to thousands of features from a high-throughput genomic experiment. In the event that no baseline demographic/clinical characteristics/subject level variables are available or needed in \code{terms} (all variables are to be penalized) then the model is specified as \code{response~1, x=penalized.terms}. 
}
\value{
\item{AIC}{a vector of AIC values for each step (if \code{x} is specified).}
\item{BIC}{a vector of BIC values for each step (if \code{x} is specified).}
\item{alpha}{the ordinal threshold estimates for the fitted model.}
\item{theta}{the coefficient estimates for the unpenalized variables (if \code{terms} are specified on the right hand side of the model formula).}
\item{beta}{the coefficient estimates for the penalized variables (if \code{x} is specified in the model).}
\item{phi}{the scaling coefficient estimates (if a \code{"Stereotype"} logit model is fit).}
\item{logLik}{a vector of log-likelihood values for each step(if \code{terms} are specified on the right hand side of the model formula).}
\item{link}{the link function used in the model fit.}
\item{model.select}{the step at which the minimum AIC was observed (if \code{terms} are specified on the right hand side of the model formula).}
\item{probability.model}{the model fit.}
\item{scale}{logical indicating whether penalized variables were centered and scaled.}
\item{w}{the unpenalized variables in the model (if any).}
\item{x}{the penalized variables in the model (if any).}
\item{y}{the ordinal response.}
}
\references{
Hastie T., Taylor J., Tibshirani R., and Walther G. (2007) Forward stagewise regression and the monotone lasso. \emph{Electronic Journal of Statistics}, 1, 1-29.
}
\author{
Kellie J. Archer, Jiayi Hou, Qing Zhou, Kyle Ferber, John G. Layne, Amanda Gentry
}
\seealso{
See Also \code{\link{coef.ordinalgmifs}}, \code{\link{summary.ordinalgmifs}}, \code{\link{plot.ordinalgmifs}}, \code{\link{predict.ordinalgmifs}}
}
\examples{
data(hccframe)
# To minimize processing time, MPO_E302_R is coerced into the model and only a subset of 
# two CpG sites (DDIT3_P1313_R and HDAC9_P137_R) are included as penalized covariates
# in this demonstration, and epsilon is set to 0.01
hcc.fit <- ordinalgmifs(group ~ MPO_E302_R, x = c("DDIT3_P1313_R", "HDAC9_P137_R"), 
	data = hccframe, epsilon = 0.01)
coef(hcc.fit)
summary(hcc.fit)
phat <- predict(hcc.fit)
head(phat$predicted)
table(phat$class, hccframe$group)
}
\keyword{ models }
\keyword{ regression }

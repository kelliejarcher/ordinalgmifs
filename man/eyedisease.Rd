\name{eyedisease}
\alias{eyedisease}
\docType{data}
\title{
Eye Disease Risk Factors
}
\description{
Eye Disease Risk Factors data from Section 9.1 of Agresti's Analysis of Ordinal Categorical Data. The primary data are from the Wisconsin Epidemiological Study of Diabetic Retinopathy. The primary outcome is severity of retinopathy which was measured in the left and right eye of every subject.  
}
\usage{data(eyedisease)}
\format{
  A data frame with 720 observations on the following 19 variables.
  \describe{
    \item{\code{rme}}{right eye macular oedema (absent = 0, present = 1)}
    \item{\code{lme}}{left eye macular oedema (absent = 0, present = 1)}
    \item{\code{rre}}{right eye refraction index}
    \item{\code{lre}}{left eye refraction index}
    \item{\code{riop}}{right eye intraocular eye pressure}
    \item{\code{liop}}{left eye intraocular eye pressure}
    \item{\code{age}}{age}
    \item{\code{diab}}{duration of diabetes (in years)}
    \item{\code{gh}}{glycosylated haemoglobin level}
    \item{\code{sbp}}{systolic blood pressure}
    \item{\code{dbp}}{diastolic blood pressure}
    \item{\code{bmi}}{body mass index}
    \item{\code{pr}}{pulse rate?}
    \item{\code{sex}}{gender (male=1, female=2)}
    \item{\code{prot}}{proteinuria (absent = 0, present = 1)}
    \item{\code{dose}}{a numeric vector}
    \item{\code{rerl}}{right eye severity of retinopathy, an ordered factor with levels \code{None} < \code{Mild} < \code{Moderate} < \code{Proliferative}}
    \item{\code{lerl}}{left eye severity of retinopathy, an ordered factor with levels \code{None} < \code{Mild} < \code{Moderate} < \code{Proliferative}}
    \item{\code{id}}{subject identifier}
  }
}
\references{
R. Klein and B.E.K. Klein and S.E. Moss and M.D. Davis and D.L. DeMets. (1984) The Wisconsin Epidemiologic Study of Diabetic Retinopathy II. Prevalence and risk of diabetic retinopathy when age at diagnosis is less than 30 years. \emph{Archives of Opthalmology} 101, 520-526.

J. Williamson and K. Kim. (1996) A global odds ratio regression model for bivariate ordered categorical data from opthalmologic studies. \emph{Statistics in Medicine} 15: 1507-1518.

A. Agresti. (2010) \emph{Analysis of Ordered Categorical Data}, Second Edition. Wiley. Hoboken, NJ. 
}
\seealso{
See Also as \code{ \link{ordinalgmifs} }
}
\examples{
data(eyedisease)
}
\keyword{datasets}

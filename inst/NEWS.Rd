\name{NEWS}
\title{News for Package ordinalgmifs}
\section{Changes in version 1.0.1}{
\itemize{
  \item Changed ordinal.gmifs to accommodate the rare instance when one predictor exhibits no variability. In the case, if scale=TRUE, an error occurs because the scale function attempts to divide by zero.
  \item Modified plot.ordinalgmifs so users can include their own xlab, ylab, or main title on plots rather than using the defaults.
}
}
\section{Changes in version 1.0.2}{
\itemize{
  \item Update predict.ordinalgmifs to enable N-fold CV.
  \item Added CITATION.
}
}
\section{Changes in version 1.0.3}{
\itemize{
  \item Added p>n to stopping criteria in ordinal.gmifs function.
  \item Added assumption=NULL parameter which can be changed for forward CR models that use Surv() as the response to 1, 2, or 3. Modified internal functions G, fn.fcr, du.fcr along with ordinalgmifs to fit the four different discrete censored survival models.
}
}
\section{Changes in version 1.0.4}{
\itemize{
  \item Change primary modeling function name from ordinal.gmifs to ordinalgmifs.
  \item Edited ordinalgmifs-package.Rd to ensure it will be in sync with the rest of the package.
  \item Remove unconditional use of suggests packages from BioConductor and remove example from vignette that depends on Biobase package.
  \item Housekeeping: Fix Rd file that had a duplicated alias.
  }
}

\section{Changes in version 1.0.5}{
\itemize{
  \item Included pgtol option in optim control argument to avoid infinite values.
  }
}

\section{Changes in version 1.0.6}{
\itemize{
  \item Fix potential PROTECT errors flagged by rchk.
  }
}

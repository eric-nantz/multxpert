\name{pvalci}
\alias{pvalci}

\title{Common \eqn{p}-value-based procedures: Simultaneous confidence intervals}
\description{
Computation of simultaneous confidence intervals for selected multiple testing
procedures based on univariate \eqn{p}-values (Bonferroni, Holm and 
fixed-sequence procedures).

}
\usage{
pvalci(rawp, est, stderror, weight, covprob, proc)
}
\arguments{
  \item{rawp}{Vector of raw \eqn{p}-values.}
  \item{est}{Vector of point estimates.}
  \item{stderror}{Vector of standard errors associated with the point estimates.}
  \item{weight}{Vector of hypothesis weights whose sum is equal to 1 (default is 
    a vector of equal weights).}
  \item{covprob}{Simultaneous coverage probability (default is 0.975).}
  \item{proc}{Vector of character strings containing the procedure name. 
    This vector should include any of the following: \code{"Bonferroni"},
    \code{"Holm"}, \code{"Fixed-sequence"}.}
}

\details{

This function computes one-sided simultaneous confidence limits for the Bonferroni, 
Holm (Holm, 1979) and fixed-sequence (Westfall and Krishen, 2001) procedures in 
in general one-sided hypothesis testing problems (equally or unequally weighted null hypotheses).   

The simultaneous confidence intervals are computed using the methods developed
in Hsu and Berger (1999), Strassburger and Bretz (2008) and Guilbaud (2008). 
For more information on the algorithms used in the function, see 
Dmitrienko et al. (2009, Section 2.6).

}
\value{
A data frame \code{result} with columns for the raw \eqn{p}-values,  point estimates,
standard errors, weights, adjusted \eqn{p}-values, and simultaneous confidence limits 
for each of the procedures. 
}

\references{

Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L., 
Tamhane, A.C., Hsu, J.C. (2009). Multiple testing methodology. 
\emph{Multiple Testing Problems in Pharmaceutical Statistics}. 
Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). Chapman and 
Hall/CRC Press, New York. \cr

Guilbaud, O. (2008). Simultaneous confidence regions corresponding to 
Holm's stepdown procedure and other closed-testing procedures. 
\emph{Biometrical Journal}. 5, 678--692. \cr

Holm, S. (1979). A simple sequentially rejective multiple test procedure. 
\emph{Scandinavian Journal of Statistics}. 6, 65--70. \cr

Hsu, J.C., Berger, R.L. (1999). Stepwise confidence intervals without 
multiplicity adjustment for dose-response and toxicity studies. 
\emph{Journal of the American Statistical Association}. 94, 468--482. \cr

Strassburger, K., Bretz, F. (2008). Compatible simultaneous lower confidence 
bounds for the Holm procedure and other Bonferroni based closed tests. 
\emph{Statistics in Medicine}. 27, 4914--4927. \cr

Westfall, P. H., Krishen, A. (2001). Optimally weighted, fixed
sequence, and gatekeeping multiple testing procedures. \emph{Journal of
Statistical Planning and Inference}. 99, 25--40. \cr

}

\source{
  \url{http://multxpert.com/wiki/MultXpert_package}
}

\seealso{\code{\link{parci}}}

\examples{

# Consider a clinical trial conducted to evaluate the effect of three
# doses of a treatment compared to a placebo with respect to a normally
# distributed endpoint 

# Three null hypotheses of no effect are tested in the trial:
# Null hypothesis H1: No difference between Dose 1 and Placebo
# Null hypothesis H2: No difference between Dose 2 and Placebo
# Null hypothesis H3: No difference between Dose 3 and Placebo

# Null hypotheses of no treatment effect are equally weighted
weight<-c(1/3,1/3,1/3)

# Treatment effect estimates (mean  dose-placebo differences)
est<-c(2.3,2.5,1.9)

# Pooled standard deviation
sd<-9.5

# Study design is balanced with 180 patients per treatment arm
n<-180

# Standard errors
stderror<-rep(sd*sqrt(2/n),3)

# T-statistics associated with the three dose-placebo tests
stat<-est/stderror

# Compute degrees of freedom
nu<-2*(n-1)

# Compute raw one-sided p-values
rawp<-1-pt(stat,nu)

# Compute lower one-sided simultaneous confidence limits 
# for the Bonferroni procedure 
pvalci(rawp,est,stderror,weight,covprob=0.975,proc="Bonferroni")

# Compute lower one-sided simultaneous confidence limits
# for the Holm and Fixed-sequence procedures
pvalci(rawp,est,stderror,weight,covprob=0.975,proc=c("Holm", "Fixed-sequence"))

}
\keyword{procedure}
\keyword{confidence limits}

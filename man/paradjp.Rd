\name{paradjp}
\alias{paradjp}

\title{Common parametric procedures: Adjusted \eqn{p}-values}
\description{

Computation of adjusted \eqn{p}-values for commonly used parametric
multiple testing procedures (single-step and step-down Dunnett procedures).

}
\usage{
paradjp(stat,n,proc)
}
\arguments{
\item{stat}{Vector of test statistics.}
\item{n}{Common sample size in each treatment group.}
\item{proc}{Vector of character strings containing the procedure name.
  This vector should include any of the following: \code{"Single-step Dunnett"}, \code{"Step-down Dunnett"}.}
}
\details{

This function computes adjusted \eqn{p}-values for the single-step Dunnett procedure
(Dunnett, 1955) and step-down Dunnett procedure (Naik, 1975; Marcus, Peritz
and Gabriel, 1976) in one-sided hypothesis testing problems with a balanced
one-way layout and equally weighted null hypotheses. For more information on the
algorithms used in the function, see Dmitrienko et al. (2009, Section 2.7).

}
\value{
A list with the following components:
\item{proc}{Name of procedure used.}
\item{result}{A data frame with columns for the test statistics,
  one-sided raw \eqn{p}-values, and one-sided adjusted \eqn{p}-values for the specified procedure.}
}
\references{

Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L.,
Tamhane, A.C., Hsu, J.C. (2009). Multiple testing methodology.
\emph{Multiple Testing Problems in Pharmaceutical Statistics}.
Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). Chapman and
Hall/CRC Press, New York. \cr

Dunnett, C.W. (1955). A multiple comparison procedure for
comparing several treatments with a control. \emph{Journal of the American
Statistical Association}. 50, 1096--1121. \cr

Marcus, R. Peritz, E., Gabriel, K.R. (1976). On closed testing
procedures with special reference to ordered analysis of variance.
\emph{Biometrika}. 63, 655--660. \cr

Naik, U.D. (1975). Some selection rules for comparing \eqn{p} processes
with a standard. \emph{Communications in Statistics. Series A}.
4, 519--535.
}

\source{
  \url{http://multxpert.com/wiki/MultXpert_package}
}

\seealso{\code{\link{pvaladjp}}}

\examples{

# Consider a clinical trial conducted to evaluate the effect of three
# doses of a treatment compared to a placebo with respect to a normally
# distributed endpoint

# Three null hypotheses of no effect are tested in the trial:
# Null hypothesis H1: No difference between Dose 1 and Placebo
# Null hypothesis H2: No difference between Dose 2 and Placebo
# Null hypothesis H3: No difference between Dose 3 and Placebo

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

# Compute one-sided adjusted p-values for the single-step Dunnett procedure
paradjp(stat, n, proc="Single-step Dunnett")

# Compute one-sided adjusted p-values for the single-step and 
# step-down Dunnett procedures
paradjp(stat, n, proc=c("Single-step Dunnett", "Step-down Dunnett"))

}
\keyword{procedure}
\keyword{\eqn{p}-value}

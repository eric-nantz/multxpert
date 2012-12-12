\name{pargateadjp}
\alias{pargateadjp}

\title{Multistage parallel gatekeeping procedures: Adjusted \eqn{p}-values}
\description{

Computation of adjusted \eqn{p}-values for multistage parallel
gatekeeping procedures.

} 

\usage{ pargateadjp(gateproc, independence, alpha, printDecisionRules)
} 

\arguments{ 

\item{gateproc}{List of gatekeeping procedure parameters
in each family of null hypotheses, including the family label, vector of
raw \eqn{p}-values, procedure name and procedure parameter
(\code{pargateadjp} function supports truncated and regular versions of
the Bonferroni, Holm, Hommel, Hochberg and fallback procedures).}
\item{independence}{Boolean indicator (TRUE, Independence condition is
imposed (i.e., inferences in earlier families are independent of inferences in later families); 
FALSE, Independence condition is not imposed).}
\item{alpha}{Global family-wise error rate (default is 0.05). Note that
this argument is not needed if the function is called to compute
adjusted \eqn{p}-values, i.e., if \code{printDecisionRules=FALSE}.}
\item{printDecisionRules}{Boolean indicator for printing the decision
rules for the gatekeeping procedure (default is FALSE).}  
} 

\details{

This function computes adjusted \eqn{p}-values and generates decision rules for multistage parallel
gatekeeping procedures in hypothesis testing problems with multiple families
of null hypotheses (null hypotheses are assumed to be equally weighted within
each family) based on the methodology presented in Dmitrienko, Tamhane
and Wiens (2008) and Dmitrienko, Kordzakhia and Tamhane (2011). For more
information on parallel gatekeeping procedures (computation of adjusted \eqn{p}-values, 
independence condition, etc), see Dmitrienko and Tamhane (2009, Section 5.4).

}
\value{
A data frame \code{result} with columns for the family labels, procedures, procedure
parameters (truncation parameters), raw \eqn{p}-values, and adjusted \eqn{p}-values.
}
\references{

Dmitrienko, A., Tamhane, A., Wiens, B. (2008). General multistage
gatekeeping procedures. \emph{Biometrical Journal}. 50, 667--677. \cr

Dmitrienko, A., Tamhane, A.C. (2009). Gatekeeping procedures in
clinical trials. \emph{Multiple Testing Problems in Pharmaceutical
Statistics}. Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors).
Chapman and Hall/CRC Press, New York. \cr

Dmitrienko, A., Kordzakhia, G., Tamhane, A.C. (2011). Multistage and mixture
parallel gatekeeping procedures in clinical trials. \emph{Journal of
Biopharmaceutical Statistics}. To appear.
}

\source{
  \url{http://multxpert.com/wiki/MultXpert_package}
}

\examples{

# Consider a clinical trial with two families of null hypotheses

# Family 1: Primary null hypotheses (one-sided p-values)
# H1 (Endpoint 1), p1=0.0082
# H2 (Endpoint 2), p2=0.0174

# Family 2: Secondary null hypotheses (one-sided p-values)
# H3 (Endpoint 3), p3=0.0042
# H4 (Endpoint 4), p4=0.0180

# Define family label and raw p-values in Family 1
label1<-"Primary endpoints"
rawp1<-c(0.0082,0.0174)

# Define family label and raw p-values in Family 2
label2<-"Secondary endpoints"
rawp2<-c(0.0042,0.0180)

# Independence condition is imposed (Families 1 and 2 are tested
# sequentually from first to last and thus adjusted p-values 
# in Family 1 do not depend on inferences in Family 2)
independence<-TRUE

# Define a two-stage parallel gatekeeping procedure which
# utilizes the truncated Holm procedure in Family 1 (truncation
# parameter=0.5) and regular Holm procedure in Family 2 (truncation
# parameter=1)

# Create a list of gatekeeping procedure parameters
family1<-list(label=label1, rawp=rawp1, proc="Holm", procpar=0.5)
family2<-list(label=label2, rawp=rawp2, proc="Holm", procpar=1)
gateproc<-list(family1,family2)

# Compute adjusted p-values
pargateadjp(gateproc, independence)

# Generate decision rules using a one-sided alpha=0.025
pargateadjp(gateproc, independence, alpha=0.025, printDecisionRules=TRUE)


}
\keyword{procedure}
\keyword{\eqn{p}-value}

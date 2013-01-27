#' Common parametric procedures: Adjusted \eqn{p}-values
#'  
#' Computation of adjusted \eqn{p}-values for commonly used parametric
#' multiple testing procedures (single-step and step-down Dunnett procedures).
#' 
#' @usage paradjp(stat,n,proc)
#' @param stat Vector of test statistics.
#' @param n Common sample size in each treatment group.
#' @param proc Vector of character strings containing the procedure name.  This vector should include any of the following: 
#'   \code{"Single-step Dunnett"}, \code{"Step-down Dunnett"}.
#' @return A list with the following components:
#'   \itemize{
#'     \item{proc}{Name of procedure used.}
#'     \item{result}{A data frame with columns for the test statistics, one-sided raw \eqn{p}-values, and one-sided 
#'        adjusted \eqn{p}-values for the specified procedure.}
#'   }
#' @references 
#'   Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L.,
#'   Tamhane, A.C., Hsu, J.C. (2009). Multiple testing methodology.
#'   \emph{Multiple Testing Problems in Pharmaceutical Statistics}.
#'   Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). Chapman and
#'   Hall/CRC Press, New York.
#'   
#'   Dunnett, C.W. (1955). A multiple comparison procedure for
#'   comparing several treatments with a control. \emph{Journal of the American
#'   Statistical Association}. 50, 1096--1121. \cr
#'
#'   Marcus, R. Peritz, E., Gabriel, K.R. (1976). On closed testing
#'   procedures with special reference to ordered analysis of variance.
#'   \emph{Biometrika}. 63, 655--660. \cr
#' 
#;   Naik, U.D. (1975). Some selection rules for comparing \eqn{p} processes
#'   with a standard. \emph{Communications in Statistics. Series A}.
#'   4, 519--535.
#' @seealso \code{\link{pvaladjp}}
#' @source \url{http://multxpert.com/wiki/MultXpert_package}
#' @export
#' @examples
#' 
#' # Consider a clinical trial conducted to evaluate the effect of three
#' # doses of a treatment compared to a placebo with respect to a normally
#' # distributed endpoint
#' 
#' # Three null hypotheses of no effect are tested in the trial:
#' # Null hypothesis H1: No difference between Dose 1 and Placebo
#' # Null hypothesis H2: No difference between Dose 2 and Placebo
#' # Null hypothesis H3: No difference between Dose 3 and Placebo
#' 
#' # Treatment effect estimates (mean  dose-placebo differences)
#' est<-c(2.3,2.5,1.9)
#' 
#' # Pooled standard deviation
#' sd<-9.5
#' 
#' # Study design is balanced with 180 patients per treatment arm
#' n<-180
#' 
#' # Standard errors
#' stderror<-rep(sd*sqrt(2/n),3)
#' 
#' # T-statistics associated with the three dose-placebo tests
#' stat<-est/stderror
#' 
#' # Compute one-sided adjusted p-values for the single-step Dunnett procedure
#' paradjp(stat, n, proc="Single-step Dunnett")
#' 
#' # Compute one-sided adjusted p-values for the single-step and 
#' # step-down Dunnett procedures
#' paradjp(stat, n, proc=c("Single-step Dunnett", "Step-down Dunnett"))
paradjp<-function(stat,n,proc=c("Single-step Dunnett", "Step-down Dunnett"))
# STAT, Vector of test statistics
# N, Common sample size in each treatment group
# PROC, Procedure name
{

	if (n<=0) stop("Sample size must be positive")

	# Number of null hypotheses
	m<-length(stat)

        if(m==0) stop("No test statistics are specified")

        if(length(n)==0) stop("No sample size specified")

	# Number of degrees of freedom
	nu<-(m+1)*(n-1)
	# Raw p-values
	rawp<-1-pt(stat,2*(n-1))
	# Adjusted p-values
	adjp<-rep(0,m)

        if(!all(proc %in% c("Single-step Dunnett", "Step-down Dunnett"))) stop("Procedure name is not recognized. ParAdjP function supports only the single-step Dunnett and step-down Dunnett procedures")

        # Number of procedures specified
        nproc <- length(proc)

        # Set up matrix to contain adjusted p-values
        adjp <- matrix(0,m,nproc)
        dimnames(adjp) <- list(NULL, paste(proc, ".adj.pvalue", sep=""))

	if (is.element("Single-step Dunnett", proc))
        {
        #for (i in 1:m) adjp[i]<-1-pdunnett(stat[i],nu,m)
        adjp[,"Single-step Dunnett.adj.pvalue"] <- sapply(stat, function(x) {1-pdunnett(x,nu,m)})
        }
	if (is.element("Step-down Dunnett", proc))
        {
        adjptmp <- rep(NA, length(stat))
		# Sort test statistics from largest to smallest
		or<-order(stat,decreasing=TRUE)
		sorted<-stat[or]
		for (i in 1:m)
		{
			if (i==1)
			{
				#adjp[1,"Step-down Dunnett.adj.pvalue"] <- 1-pdunnett(sorted[1],nu,m)
				adjptmp[1] <- 1-pdunnett(sorted[1],nu,m)
                                #print(1-pdunnett(sorted[1],nu,m))
				#maxp<-adjp[1, "Step-down Dunnett.adj.pvalue"]
				maxp<-adjptmp[1]
			}
			if (i>1 & i<m)
			{
				#adjp[i, "Step-down Dunnett.adj.pvalue"] <- max(maxp,1-pdunnett(sorted[i],nu,m-i+1))
				adjptmp[i] <- max(maxp,1-pdunnett(sorted[i],nu,m-i+1))

				#maxp<-max(maxp,adjp[i, "Step-down Dunnett.adj.pvalue"])
				maxp<-max(maxp,adjptmp[i])
			}
			if (i==m) {
                            #adjp[m, "Step-down Dunnett.adj.pvalue"] <- max(maxp,1-pt(sorted[m],nu))
                            adjptmp[m] <- max(maxp,1-pt(sorted[m],nu))
                        }
		}
		# Return to original ordering
		temp<-adjptmp
		adjptmp[or]<-temp
                adjp[,"Step-down Dunnett.adj.pvalue"] <- adjptmp
	}

	# Data frame returned by the function
	result<-data.frame(stat, round(rawp, 4), adjp)
	names(result)[1]<-"Test.statistic"
	names(result)[2]<-"Raw.pvalue"
	#names(result)[3]<-"Adj.pvalue"

	return(result=result)
}
# End of paradjp

#' Common parametric procedures: Simultaneous confidence intervals
#'  
#' Computation of simultaneous confidence intervals for commonly used parametric
#' multiple testing procedures (single-step and step-down Dunnett procedures).
#' 
#' @usage parci(stat, n, est, stderror, covprob, proc)
#' @param stat Vector of test statistics.
#' @param n Common sample size in each treatment group.
#' @param est Vector of point estimates
#' @param stderror Vector of standard errors associated with the point estimates.
#' @param covprob Simultaneous coverage probability (default is 0.975).
#' @param proc Vector of character strings containing the procedure name. This vector should include any of the following: 
#'   \code{"Single-step Dunnett"}, \code{"Step-down Dunnett"}.
#' @details 
#'   This function computes lower one-sided simultaneous confidence limits for the single-step Dunnett procedure
#'   (Dunnett, 1955) and step-down Dunnett procedure (Naik, 1975; Marcus, Peritz and Gabriel, 1976) in one-sided hypothesis testing 
#'   problems with a balanced one-way layout and equally weighted null hypotheses.
#'
#'   The simultaneous confidence intervals are computed using the methods developed
#'   in Bofinger (1987) and Stefansson, Kim and Hsu (1988). For more information on the
#'   algorithms used in the function, see Dmitrienko et al. (2009, Section 2.7).
#' @return A data frame \code{result} with columns for the test statistics, point estimates,
#'  standard errors, adjusted \eqn{p}-values, and lower simultaneous confidence limits
#'  for the specified procedure.
#' @references 
#'   Bofinger, E. (1987). Step-down procedures for comparison with a control.
#'   \emph{Australian Journal of Statistics}. 29, 348--364. \cr
#'   
#'   Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L.,
#'   Tamhane, A.C., Hsu, J.C. (2009). Multiple testing methodology.
#'   \emph{Multiple Testing Problems in Pharmaceutical Statistics}.
#'   Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). Chapman and
#'   Hall/CRC Press, New York. \cr
#'   
#'   Dunnett, C.W. (1955). A multiple comparison procedure for
#'   comparing several treatments with a control. \emph{Journal of the American
#'   Statistical Association}. 50, 1096--1121. \cr
#'   
#'   Marcus, R. Peritz, E., Gabriel, K.R. (1976). On closed testing
#'   procedures with special reference to ordered analysis of variance.
#'   \emph{Biometrika}. 63, 655--660. \cr
#'   
#'   Naik, U.D. (1975). Some selection rules for comparing \eqn{p} processes
#'   with a standard. \emph{Communications in Statistics. Series A}.
#'   4, 519--535. \cr
#'   
#'   Stefansson, G., Kim, W.-C., Hsu, J.C. (1988). On confidence sets in multiple
#'   comparisons. \emph{Statistical Decision Theory and Related Topics IV}. Gupta, S.S.,
#'   Berger, J.O. (editors). Academic Press, New York, 89--104.
#' @seealso \code{\link{parci}}
#' @source \url{http://multxpert.com/wiki/MultXpert_package}
#' @keywords procedure, confidence limits
#' @export
#' 
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
#' # Compute lower one-sided simultaneous confidence limits
#' # for the single-step Dunnett procedure
#' parci(stat,n,est,stderror,covprob=0.975,proc="Single-step Dunnett")
#' 
#' # Compute lower one-sided simultaneous confidence limits
#' # for the single-step and step-down Dunnett procedures
#' parci(stat,n,est,stderror,covprob=0.975,proc=c("Single-step Dunnett", "Step-down Dunnett"))
parci<-function(stat,n,est,stderror,covprob=0.975,proc) { 
# STAT, Vector of test statistics
# N, Common sample size in each treatment group
# EST, Vector of point estimates
# STDERROR, Vector of standard errors associated with the point estimates
# COVPROB, Simultaneous coverage probability
# PROC, Procedure name
    
    # Number of null hypotheses
    m<-length(stat)

    if (m==0) stop("No test statistics are specified")

    if (m!=length(est)) stop("STAT and EST vectors have different lengths")
    if (m!=length(stderror)) stop("STAT and STDERROR vectors have different lengths")

    if (covprob>=1) stop("Simultaneous coverage probability must be <1")
    if (covprob<=0) stop("Simultaneous coverage probability must be >0")

    if(!all(proc %in% c("Single-step Dunnett", "Step-down Dunnett"))) stop("Procedure name is not recognized. ParCI function supports only the single-step Dunnett and step-down Dunnett procedures")

    if (n<=0) stop("Sample size must be positive")

    # number of procedures specified
    nproc <- length(proc)

    # set up matrix to contain confidence limits
    cimat <- matrix(0,m,nproc)
    dimnames(cimat) <- list(NULL, paste(proc, ".conf.limit", sep=""))

    # set up matrix to contain adjusted p-values
    adjpmat <- matrix(0,m,nproc)
    dimnames(adjpmat) <- list(NULL, paste(proc, ".adj.pvalue", sep=""))

    # Degrees of freedon
    nu<-(m+1)*(n-1)

    # Compute adjusted p-values
    result <- paradjp(stat,n,proc)

    #adjpmat <- result[, grep(".adj.pvalue", names(result), value=TRUE)]

    # One-sided familywise error rate
    alpha<-1-covprob

    # Rejection/acceptance of null hypotheses
    #reject<-(adjp<=alpha)

    # Vectors of confidence limits
    ci<-rep(0,m)

    zero<-rep(0,m)

	if (is.element("Single-step Dunnett", proc)) {
            adjpmat[, "Single-step Dunnett.adj.pvalue"] <- round(result[, "Single.step.Dunnett.adj.pvalue"], 4)

            reject <- (result[, "Single.step.Dunnett.adj.pvalue"] <= alpha)
           # Critical value

           c<-qdunnett(1-alpha,nu,m)
           cimat[, "Single-step Dunnett.conf.limit"] <-round(est-c*stderror, 4)
        }

	if (is.element("Step-down Dunnett", proc)) {
            adjpmat[, "Step-down Dunnett.adj.pvalue"] <- round(result[, "Step.down.Dunnett.adj.pvalue"], 4)

            reject <- (result[, "Step.down.Dunnett.adj.pvalue"] <= alpha)

	    # All null hypotheses are rejected
  	    if (sum(reject)==m) {
               # Critical value
               c<-qt(1-alpha,nu)
               cimat[, "Step-down Dunnett.conf.limit"] <- round(pmax(zero,est-c*stderror), 4)
            }

            # Some null hypotheses are accepted
  	    if (sum(reject)<m) {
                for (i in 1:m) {
                   if (reject[i]==1) cimat[i, "Step-down Dunnett.conf.limit"]<-0
                   if (reject[i]==0) {
                      # Critical value
                      c<-qdunnett(1-alpha,nu,m-sum(reject))
                      cimat[i, "Step-down Dunnett.conf.limit"] <- round(est[i]-c*stderror[i],4)
                   }
                }
            }
        }


    # Data frame returned by the function
    result<-data.frame(stat, est, stderror, adjpmat, cimat)
    names(result)[1]<-"Test.statistic"
    names(result)[2]<-"Estimate"
    names(result)[3]<-"Std.error"
    #names(result)[4]<-"Adj.pvalue"
    #names(result)[5]<-"Conf.limit"

    return(result=result)

    }
# End of parci

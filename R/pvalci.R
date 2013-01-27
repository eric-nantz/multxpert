#' Common \eqn{p}-value-based procedures: Simultaneous confidence intervals
#' 
#' Computation of simultaneous confidence intervals for selected multiple testing
#' procedures based on univariate \eqn{p}-values (Bonferroni, Holm and fixed-sequence procedures).
#' 
#' @usage pvalci(rawp, est, stderror, weight, covprob, proc)
#' @param rawp Vector of raw \eqn{p}-values.
#' @param est Vector of point estimates
#' @param stderror Vector of standard errors associated with the point estimates.
#' @param weight Vector of hypothesis weights whose sum is equal to 1 (default is a vector of equal weights).
#' @param covprob Simultaneous coverage probability (default is 0.975).
#' @param proc Vector of character strings containing the procedure name. This vector should include any of the following: 
#'   \code{"Bonferroni"}, \code{"Holm"}, \code{"Fixed-sequence"}.
#'
#' @details
#'   This function computes one-sided simultaneous confidence limits for the Bonferroni, 
#'   Holm (Holm, 1979) and fixed-sequence (Westfall and Krishen, 2001) procedures in 
#'   in general one-sided hypothesis testing problems (equally or unequally weighted null hypotheses).   
#'   
#'   The simultaneous confidence intervals are computed using the methods developed
#'   in Hsu and Berger (1999), Strassburger and Bretz (2008) and Guilbaud (2008). 
#'   For more information on the algorithms used in the function, see 
#'   Dmitrienko et al. (2009, Section 2.6).
#'
#' @return A data frame \code{result} with columns for the raw \eqn{p}-values,  point estimates,
#'   standard errors, weights, adjusted \eqn{p}-values, and simultaneous confidence limits 
#'   for each of the procedures. 
#' 
#' @references
#'   Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L., 
#'   Tamhane, A.C., Hsu, J.C. (2009). Multiple testing methodology. 
#'   \emph{Multiple Testing Problems in Pharmaceutical Statistics}. 
#'   Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). Chapman and 
#'   Hall/CRC Press, New York. \cr
#'   
#'   Guilbaud, O. (2008). Simultaneous confidence regions corresponding to 
#'   Holm's stepdown procedure and other closed-testing procedures. 
#'   \emph{Biometrical Journal}. 5, 678--692. \cr
#'   
#'   Holm, S. (1979). A simple sequentially rejective multiple test procedure. 
#'   \emph{Scandinavian Journal of Statistics}. 6, 65--70. \cr
#'   
#'   Hsu, J.C., Berger, R.L. (1999). Stepwise confidence intervals without 
#'   multiplicity adjustment for dose-response and toxicity studies. 
#'   \emph{Journal of the American Statistical Association}. 94, 468--482. \cr
#'   
#'   Strassburger, K., Bretz, F. (2008). Compatible simultaneous lower confidence 
#'   bounds for the Holm procedure and other Bonferroni based closed tests. 
#'   \emph{Statistics in Medicine}. 27, 4914--4927. \cr
#'   
#'   Westfall, P. H., Krishen, A. (2001). Optimally weighted, fixed
#'   sequence, and gatekeeping multiple testing procedures. \emph{Journal of
#'   Statistical Planning and Inference}. 99, 25--40. \cr
#' 
#' @source \url{http://multxpert.com/wiki/MultXpert_package}
#' @seealso \code{\link{parci}}
#' @keywords procedure, confidence limits
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
#' # Null hypotheses of no treatment effect are equally weighted
#' weight<-c(1/3,1/3,1/3)
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
#' # Compute degrees of freedom
#' nu<-2*(n-1)
#' 
#' # Compute raw one-sided p-values
#' rawp<-1-pt(stat,nu)
#' 
#' # Compute lower one-sided simultaneous confidence limits 
#' # for the Bonferroni procedure 
#' pvalci(rawp,est,stderror,weight,covprob=0.975,proc="Bonferroni")
#' 
#' # Compute lower one-sided simultaneous confidence limits
#' # for the Holm and Fixed-sequence procedures
#' pvalci(rawp,est,stderror,weight,covprob=0.975,proc=c("Holm", "Fixed-sequence"))
pvalci<-function(rawp,est,stderror,weight=rep(1/length(rawp),length(rawp)),covprob=0.975,proc=c("Bonferroni", "Holm", "Fixed-sequence")) {
  # The PvalCI function computes one-sided multiplicity-adjusted confidence
  # intervals (simultaneous confidence intervals) for the Bonferroni, Holm and fixed-sequence
  # procedures in general hypothesis testing problems with equally or
  # unequally weighted null hypotheses
  
  # RAWP, Vector of raw p-values
  # EST, Vector of point estimates
  # STDERROR, Vector of standard errors associated with the point estimates
  # WEIGHT, Vector of hypothesis weights
  # COVPROB, Simultaneous coverage probability
  # PROC, Procedure name

    # Number of null hypotheses
    m<-length(rawp)

    if (m==0) stop("No p-values are specified")

    for (i in 1:m)
        {
        if (rawp[i]<0) stop("P-values must be positive")
        if (rawp[i]>1) stop("P-values must be less than 1")
        }

    if (m!=length(weight)) stop("RAWP and WEIGHT vectors have different lengths")
    if (m!=length(est)) stop("RAWP and EST vectors have different lengths")
    if (m!=length(stderror)) stop("RAWP and STDERROR vectors have different lengths")

    if (sum(weight)>1) stop("Sum of hypothesis weights must be <=1")

    for (i in 1:length(weight))
        {
        if (weight[i]<0) stop("Hypothesis weights must be >=0")
        }

    if (covprob>=1) stop("Simultaneous coverage probability must be <1")
    if (covprob<=0) stop("Simultaneous coverage probability must be >0")

    if (!all(proc %in% c("Bonferroni", "Holm", "Fixed-sequence")))
        stop("Procedure name is not recognized")
    #if (proc!="Bonferroni" & proc!="Holm" & proc!="Fixed-sequence") stop("Procedure name is not recognized")

    # number of procedures specified
    nproc <- length(proc)

        # set up matrix to contain confidence limits
    cimat <- matrix(0,m,nproc)
    dimnames(cimat) <- list(NULL, paste(proc, ".conf.limit", sep=""))

 	# One-sided familywise error rate
 	alpha<-1-covprob

    # Compute adjusted p-values
    result <- pvaladjp(rawp=rawp,weight=weight,alpha=alpha,proc=proc)
    adjpmat <- result[, grep(".adj.pvalue", names(result), value=TRUE)]
    #adjp<-out[c(-1,-2)]
    #print(out)

        # Rejection/acceptance of null hypotheses
	#reject<-(adjp<=alpha)

    # Vectors of confidence limits
    ci<-rep(0,m)

 	zero<-rep(0,m)


   # adjp<-out[,3]

    # Bonferroni procedure
	if (is.element("Bonferroni", proc)) {
           reject <- (result[,"Bonferroni.adj.pvalue"] <= alpha)
           cimat[,"Bonferroni.conf.limit"] <-est-(stderror*qnorm(1-(alpha*weight)))
        }

    # Holm procedure
	if (is.element("Holm", proc)) {
            reject <- (result[,"Holm.adj.pvalue"] <= alpha)
        bonfci<-est-(stderror*qnorm(1-(alpha*weight)))
        # All null hypotheses are rejected
 		if (sum(reject)==m) cimat[,"Holm.conf.limit"] <-pmax(zero,bonfci)
        # Some null hypotheses are accepted
        if (sum(reject)<m)
            {
            for(i in 1:m)
                {
                if (reject[i]==1) cimat[i, "Holm.conf.limit"]<-0
                if (reject[i]==0)
                    {
                    adjalpha<-(alpha*weight[i])/sum(weight[reject==0])
                    cimat[i,"Holm.conf.limit"]<-est[i]-(stderror[i]*qnorm(1-adjalpha))
                    }
                }
            }
        }

   	# Fixed-sequence procedure
	if (is.element("Fixed-sequence", proc)) {
            reject <- (result[,"Fixed.sequence.adj.pvalue"] <= alpha)
		# All null hypotheses are accepted
  		if (sum(reject)==0)
            {
            cimat[1,"Fixed-sequence.conf.limit"] <-est[1]-stderror[1]*qnorm(1-alpha)
            for(i in 2:m) cimat[i, "Fixed-sequence.conf.limit"]<-NA
            }
        # All null hypotheses are rejected
		if (sum(reject)==m)
            {
            temp1<-est-stderror*qnorm(1-alpha)
            cimat[,"Fixed-sequence.conf.limit"] <-min(temp1)
            }
        # Some null hypotheses are accepted and some are rejected
        if (sum(reject)>0 & sum(reject)<m)
            {
            cimat[1, "Fixed-sequence.conf.limit"]<-0
            for(i in 2:m)
                {
                if (reject[i]==1) cimat[i, "Fixed-sequence.conf.limit"]<-0
                if (reject[i]==0 & reject[i-1]==1) cimat[i, "Fixed-sequence.conf.limit"]<-est[i]-stderror[i]*qnorm(1-alpha)
                if (reject[i]==0 & reject[i-1]==0) cimat[i, "Fixed-sequence.conf.limit"]<-NA
                }
            }
        }

	# Data frame returned by the function
    #result<-data.frame(rawp,est,stderror,weight,adjp,ci)
    result<-data.frame(rawp,est,stderror,weight,adjpmat,cimat)
    names(result)[1]<-"Raw.pvalue"
    names(result)[2]<-"Estimate"
    names(result)[3]<-"Std.error"
    names(result)[4]<-"Weight"
    #names(result)[5]<-"Adj.pvalue"
    #names(result)[6]<-"Conf.limit"

    return(result=result)
}
# End of pvalci

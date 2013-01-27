#' Common \eqn{p}-value-based procedures: Adjusted \eqn{p}-values
#' 
#' Computation of adjusted \eqn{p}-values for commonly used multiple testing
#' procedures based on univariate \eqn{p}-values (Bonferroni, Holm, Hommel, Hochberg, 
#' fixed-sequence and fallback procedures).
#' 
#' @usage pvaladjp(rawp, weight, alpha, proc, printDecisionRules)
#' 
#' @param rawp Vector of raw \eqn{p}-values.
#' @param weight Vector of hypothesis weights whose sum is equal to 1 (default is a vector of equal weights).
#' @param alpha Familywise error rate (default is 0.05). Note that this argument is not needed if the function is called 
#'   to compute adjusted \eqn{p}-values, i.e., if \code{printDecisionRules=FALSE}.
#' @param proc Vector of character strings containing the procedure name. This vector should include any of the following: 
#'   \code{"Bonferroni"}, \code{"Holm"}, \code{"Hommel"}, \code{"Hochberg"}, \code{"Fixed-sequence"}, \code{"Fallback"}.
#' @param printDecisionRules Boolean indicator for printing the decision rules for each of the procedures specified 
#'   in \code{"proc"} (default is FALSE).
#' 
#' @details
#'   This function computes adjusted \eqn{p}-values and generates decision rules for the Bonferroni, 
#'   Holm (Holm, 1979), Hommel (Hommel, 1988), Hochberg (Hochberg, 1988), 
#'   fixed-sequence (Westfall and Krishen, 2001) and fallback (Wiens, 2003; Wiens and Dmitrienko, 2005) procedures. 
#'   
#'   The adjusted \eqn{p}-values are computed using the closure principle (Marcus, Peritz and Gabriel, 1976) in general 
#'   hypothesis testing problems (equally or unequally weighted null hypotheses).   The decision rules are generated only 
#'   in hypothesis testing problems with equally weighted null hypotheses.  For more information on the algorithms used 
#'   in the function, see Dmitrienko et al. (2009, Section 2.6).
#'
#' @return A data frame \code{result} with columns for the raw \eqn{p}-values, weights, and adjusted \eqn{p}-values for each of 
#' the procedures.
#' 
#' @references
#'   Dmitrienko, A., Bretz, F., Westfall, P.H., Troendle, J., Wiens, B.L., 
#'   Tamhane, A.C., Hsu, J.C. (2009). Multiple testing methodology. 
#'   \emph{Multiple Testing Problems in Pharmaceutical Statistics}. 
#'   Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors). Chapman and 
#'   Hall/CRC Press, New York. \cr
#'   
#'   Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple significance testing. 
#'   \emph{Biometrika}. 75, 800--802. \cr
#'   
#'   Holm, S. (1979). A simple sequentially rejective multiple test procedure. 
#'   \emph{Scandinavian Journal of Statistics}. 6, 65--70. \cr
#'   
#'   Hommel, G. (1988). A stagewise rejective multiple test procedure based on a 
#'   modified Bonferroni test. \emph{Biometrika}. 75, 383--386. \cr
#'   
#'   Marcus, R. Peritz, E., Gabriel, K.R. (1976). On closed testing
#'   procedures with special reference to ordered analysis of variance. 
#'   \emph{Biometrika}. 63, 655--660. \cr
#'   
#'   Westfall, P. H., Krishen, A. (2001). Optimally weighted, fixed
#'   sequence, and gatekeeping multiple testing procedures. \emph{Journal of
#'   Statistical Planning and Inference}. 99, 25--40. \cr
#'   
#'   Wiens, B. (2003). A fixed-sequence Bonferroni procedure for
#'   testing multiple endpoints. \emph{Pharmaceutical Statistics}. 2, 211--215. \cr
#'   
#'   Wiens, B., Dmitrienko, A. (2005). The fallback procedure for
#'   evaluating a single family of hypotheses. \emph{Journal of
#'   Biopharmaceutical Statistics}. 15, 929--942.
#' 
#' @source \url{http://multxpert.com/wiki/MultXpert_package}
#' @seealso \code{\link{paradjp}}
#' @keywords procedure, \eqn{p}-value
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
#' # Compute adjusted p-values for the Bonferroni procedure 
#' pvaladjp(rawp, weight, proc="Bonferroni")
#' 
#' # Compute adjusted p-values for the Hommel and Fallback procedures
#' pvaladjp(rawp, weight, proc=c("Hommel", "Fallback"))
#' 
#' # Generate decision rules for the Holm procedure 
#' # using a one-sided alpha=0.025
#' pvaladjp(rawp, weight, alpha=0.025, proc="Holm", printDecisionRules=TRUE)
pvaladjp<-function(rawp,weight=rep(1/length(rawp), length(rawp)),alpha=0.05,
		               proc=c("Bonferroni", "Holm", "Hommel", "Hochberg", "Fixed-sequence", "Fallback"),
				           printDecisionRules=FALSE) {
  
  # PValAdjP function computes adjusted p-values and generates decision rules
  # for the Bonferroni, Holm, Hommel, Hochberg, fixed-sequence and fallback procedures
  
  # RAWP, Vector of raw p-values
  # WEIGHT, Vector of hypothesis weights
  # ALPHA, Familywise error rate
  # PROC, Procedure name
  # PRINTDECISIONRULES: Boolean indicator which controls printing of decision rules

	# Number of null hypotheses
	m<-length(rawp)

	if (m==0) stop("No p-values are specified")

	for (i in 1:m)
	{
		if (rawp[i]<0) stop("P-values must be positive")
		if (rawp[i]>1) stop("P-values must be less than 1")
	}

	index <-order(rawp)

	if (alpha <= 0) stop("Alpha must be positive")
	if (alpha >= 1) stop("Alpha must be less than 1")

	# Number of weights
	nweis<-length(weight)

	if (m!=nweis) stop("RAWP and WEIGHT vectors have different lengths")

	if (sum(weight)>1) stop("Sum of hypothesis weights must be <=1")

	for (i in 1:nweis)
	{
		if (weight[i]<0) stop("Hypothesis weights must be >=0")
	}

	if(!all(proc %in% c("Bonferroni", "Holm", "Hommel", "Hochberg", "Fixed-sequence", "Fallback")))
		stop("Procedure name is not recognized. PvalAdjp function supports the Bonferroni, Holm, Hommel, Hochberg, Fixed-sequence, and Fallback procedures")

	# number of procedures specified
	nproc <- length(proc)

	# set up matrix to contain adjusted p-values
	adjp <- matrix(0,m,nproc)
	dimnames(adjp) <- list(NULL, paste(proc, ".adj.pvalue", sep=""))

	if (is.element("Bonferroni", proc)) {
		adjp[,"Bonferroni.adj.pvalue"]<-pvaltrunc(rawp,weight,"Holm",0)
	}
	if (is.element("Holm", proc)) {
		adjp[,"Holm.adj.pvalue"]<-pvaltrunc(rawp,weight,"Holm",1)
	}
	if (is.element("Hommel", proc)) {
		adjp[,"Hommel.adj.pvalue"]<-pvaltrunc(rawp,weight,"Hommel",1)
	}
	if (is.element("Hochberg", proc)) {
		adjp[,"Hochberg.adj.pvalue"]<-pvaltrunc(rawp,weight,"Hochberg",1)
	}
	if (is.element("Fixed-sequence", proc)) {
		adjp[,"Fixed-sequence.adj.pvalue"]<-pvaltrunc(rawp,weight,"Fixed-sequence",1)
	}
	if (is.element("Fallback", proc)) {
		adjp[,"Fallback.adj.pvalue"]<-pvaltrunc(rawp,weight,"Fallback",1)
	}

	# Data frame returned by the function
	result<-data.frame(round(rawp,4),round(weight,4), round(adjp,4))
	names(result)[1]<-"Raw.pvalue"
	names(result)[2]<-"Weight"

	if(printDecisionRules==TRUE) {
            if(length(unique(round(weight,3))) > 1) stop("Weights must be equal for decision rule calculations to be valid")
            pvalrule(rawp, alpha, proc)
        }

	return(result=result)
}
# End of pvaladjp

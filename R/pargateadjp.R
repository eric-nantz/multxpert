#' Multistage parallel gatekeeping procedures: Adjusted \eqn{p}-values
#' 
#' Computation of adjusted \eqn{p}-values for multistage parallel gatekeeping procedures.
#' 
#' @usage pargateadjp(gateproc, independence, alpha, printDecisionRules)
#' @param gateproc List of gatekeeping procedure parameters in each family of null hypotheses, including the family label, vector of
#'   raw \eqn{p}-values, procedure name and procedure parameter. (\code{pargateadjp} function supports truncated and regular versions of
#'   the Bonferroni, Holm, Hommel, Hochberg and fallback procedures).
#' @param independence Boolean indicator (TRUE, Independence condition is imposed (i.e., inferences in earlier families 
#'   are independent of inferences in later families); FALSE, Independence condition is not imposed).
#' @param alpha Global family-wise error rate (default is 0.05). Note that this argument is not needed if the function is called 
#'   to compute adjusted \eqn{p}-values, i.e., if \code{printDecisionRules=FALSE}.
#' @param printDecisionRules Boolean indicator for printing the decision rules for the gatekeeping procedure (default is FALSE).
#' 
#' @details
#'   This function computes adjusted \eqn{p}-values and generates decision rules for multistage parallel
#'   gatekeeping procedures in hypothesis testing problems with multiple families
#'   of null hypotheses (null hypotheses are assumed to be equally weighted within
#'   each family) based on the methodology presented in Dmitrienko, Tamhane
#'   and Wiens (2008) and Dmitrienko, Kordzakhia and Tamhane (2011). For more
#'   information on parallel gatekeeping procedures (computation of adjusted \eqn{p}-values, 
#'   independence condition, etc), see Dmitrienko and Tamhane (2009, Section 5.4).
#'   
#' @return A data frame \code{result} with columns for the family labels, procedures, procedure
#' parameters (truncation parameters), raw \eqn{p}-values, and adjusted \eqn{p}-values.
#' 
#' @references 
#'   Dmitrienko, A., Tamhane, A., Wiens, B. (2008). General multistage
#'   gatekeeping procedures. \emph{Biometrical Journal}. 50, 667--677. \cr
#'   
#'   Dmitrienko, A., Tamhane, A.C. (2009). Gatekeeping procedures in
#'   clinical trials. \emph{Multiple Testing Problems in Pharmaceutical
#'   Statistics}. Dmitrienko, A., Tamhane, A.C., Bretz, F. (editors).
#'   Chapman and Hall/CRC Press, New York. \cr
#'   
#'   Dmitrienko, A., Kordzakhia, G., Tamhane, A.C. (2011). Multistage and mixture
#'   parallel gatekeeping procedures in clinical trials. \emph{Journal of
#'   Biopharmaceutical Statistics}. To appear.
#'   
#' @source \url{http://multxpert.com/wiki/MultXpert_package}
#' @keywords procedure, \eqn{p}-value
#' 
#' @export
#' @examples
#' 
#' # Consider a clinical trial with two families of null hypotheses
#' 
#' # Family 1: Primary null hypotheses (one-sided p-values)
#' # H1 (Endpoint 1), p1=0.0082
#' # H2 (Endpoint 2), p2=0.0174
#' 
#' # Family 2: Secondary null hypotheses (one-sided p-values)
#' # H3 (Endpoint 3), p3=0.0042
#' # H4 (Endpoint 4), p4=0.0180
#' 
#' # Define family label and raw p-values in Family 1
#' label1<-"Primary endpoints"
#' rawp1<-c(0.0082,0.0174)
#' 
#' # Define family label and raw p-values in Family 2
#' label2<-"Secondary endpoints"
#' rawp2<-c(0.0042,0.0180)
#' 
#' # Independence condition is imposed (Families 1 and 2 are tested
#' # sequentually from first to last and thus adjusted p-values 
#' # in Family 1 do not depend on inferences in Family 2)
#' independence<-TRUE
#' 
#' # Define a two-stage parallel gatekeeping procedure which
#' # utilizes the truncated Holm procedure in Family 1 (truncation
#' # parameter=0.5) and regular Holm procedure in Family 2 (truncation
#' # parameter=1)
#' 
#' # Create a list of gatekeeping procedure parameters
#' family1<-list(label=label1, rawp=rawp1, proc="Holm", procpar=0.5)
#' family2<-list(label=label2, rawp=rawp2, proc="Holm", procpar=1)
#' gateproc<-list(family1,family2)
#' 
#' # Compute adjusted p-values
#' pargateadjp(gateproc, independence)
#' 
#' # Generate decision rules using a one-sided alpha=0.025
#' pargateadjp(gateproc, independence, alpha=0.025, printDecisionRules=TRUE)
pargateadjp<-function(gateproc, independence, alpha=0.05, printDecisionRules=FALSE) {
  # ParGateAdjP function computes adjusted p-values and generates decision rules
  # for multistage parallel gatekeeping procedures in hypothesis testing problems
  # with multiple families of null hypotheses (null hypotheses are assumed
  # to be equally weighted within each family)
  
  # GATEPROC, List of gatekeeping procedure parameters
  # INDEPENDENCE, Boolean indicator (TRUE, Independence condition is imposed; FALSE,
  # Independence condition is not imposed)
  # ALPHA: Global familywise error rate
  # PRINTDECISIONRULES: Boolean indicator which controls printing of decision rules

	# Number of families
	nfams<-length(gateproc)
	if (nfams<=1) stop("Function requires more than one family of null hypotheses")

	for (i in 1:nfams)
	{
		pr<-gateproc[[i]]$proc
		if (pr!="Bonferroni" & pr!="Holm" & pr!="Hommel" & pr!="Hochberg" & pr!="Fallback")
			stop("Procedure name is not recognized. ParGateAdjP function supports only the Bonferroni, Holm, Hommel, Hochberg and fallback procedures")
	}

	if (alpha <= 0) stop("Alpha must be positive")
	if (alpha >= 1) stop("Alpha must be less than 1")

	temp<-gateproc

	for (i in 1:nfams)
	{
		# Number of null hypotheses
		nhyps<-length(temp[[i]]$rawp)
		adjp<-rep(0,nhyps)
		# Placeholder for adjusted p-values
		gateproc[[i]][5]<-list(adjp=adjp)

		for (j in 1:nhyps)
		{

			# Find the lowest alpha level at which the current null hypothesis is rejected
			upper<-1
			lower<-0
			for (k in 1:20)
			{
				current<-(lower+upper)/2
				# Evaluate decision rules for the multistage parallel gatekeeping procedure
				res<-pargateeval(temp,current,independence)
				# Rejection decision for the current null hypothesis
				if (independence==TRUE | i==nfams) reject<-res[[i]][[7]][j]
				# Rejection decisions after retesting if the independence condition is not imposed
				if (independence==FALSE & i<nfams)
				{
					# If the current null hypothesis was retested
					if (res[[2*nfams-i]][[6]]>0) reject<-res[[2*nfams-i]][[7]][j] else reject<-res[[i]][[7]][j]
				}
				# Update the interval
				if (reject==TRUE) upper<-current
				if (reject==FALSE) lower<-current
			}

			# Global adjusted p-value
			gateproc[[i]][[5]][j]<-(lower+upper)/2
		}
	}

	# Build a data frame with the raw and global adjusted p-values
	count<-0
	for (i in 1:nfams) {
		count<-count+length(gateproc[[i]]$rawp)
	}
        result <- data.frame()
	k<-1
	for (i in 1:nfams)
	{
		# Number of null hypotheses
		nhyps<-length(gateproc[[i]]$rawp)
		for (j in 1:nhyps)
		{
			result[k,1]<-gateproc[[i]]$label
			result[k,2]<-gateproc[[i]]$proc
			result[k,3]<-gateproc[[i]]$procpar
			result[k,4]<-round(gateproc[[i]]$rawp[j], 4)
			result[k,5]<-round(gateproc[[i]][[5]][j], 4)
			k<-k+1
		}
	}
	names(result)[1]<-"Family"
	names(result)[2]<-"Procedure"
	names(result)[3]<-"Parameter"
	names(result)[4]<-"Raw.pvalue"
	names(result)[5]<-"Adj.pvalue"

	if(printDecisionRules==TRUE) { pargaterule(gateproc,alpha,independence)}

	return(result=result)
}
# End of pargateadjp

# Multxpert package implements commonly used p-value-based and parametric
# multiple testing procedures and parallel gatekeeping procedures

# For more information about the Multxpert package, visit the Multiplicity Expert web site
# http://multxpert.com/wiki/MultXpert_package

.onLoad <- function(lib, pkg) {
	if (interactive()) {
		packageStartupMessage("multxpert: Implementation of commonly used p-value based and parametric multiple testing procedures and parallel gatekeeping procedures. For more information visit http://multxpert.com/wiki/MultXpert_package")
	}
}
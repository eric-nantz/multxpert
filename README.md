# multxpert

The R package **multxpert** provides software implementation of commonly used p-value-based and parametric multiple testing procedures and gatekeeping procedures.  For more information, visit the [multxpert home site](http://multxpert.com/wiki/MultXpert_package).


## Installation

You can install the stable version on [CRAN](http://cran.r-project.org/package=multxpert):

```r
install.packages('multxpert', dependencies = TRUE)
```

Or install the development version of multxpert using any of the following methods:

 * Download the [zip ball](https://github.com/eric-nantz/multxpert/zipball/master) or [tar ball](https://github.com/eric-nantz/multxpert/tarball/master), decompress and run `R CMD INSTALL` on it.
 * Use the **devtools** package to install the absolutely latest version:

```r
## this package depends on R >= 2.1.0
## you may also need to update your packages: 
## options(repos = c(CRAN = 'http://cran.r-project.org'))
## update.packages()
library(devtools); install_github('multxpert', 'eric-nantz')
```

Note: Windows users have to first install [Rtools](http://cran.rstudio.com/bin/windows/Rtools/) in order to use the **devtools** package.


#-------------------------------------------------------
# script using devtools package to compile package
#-------------------------------------------------------

library(devtools)

# turn on dev_mode

dev_mode()

# load all function files
load_all("~/Dropbox/rpodcast_code/multxpert/")

# use this to generate manual files from R script files with Roxygen style 
document("~/Dropbox/rpodcast_code/multxpert/")

# check a package can be bult
check("~/Dropbox/rpodcast_code/multxpert/")

# turn off dev_mode

dev_mode()

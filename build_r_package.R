# 1 - just first time necessary
# building the r package
rm(list=ls())
code_files = c("package/bgev_functions.R","package/dist_check.R")
package_name = 'bgev'
path_for_package_files = 'bgev_cran'
package.skeleton(name=package_name, 
                 path=path_for_package_files, 
                 code_files = code_files)


# 2 
# build package with R: after buiding, or if you have already build and are editing just code for new version next step is to build package
# on terminal navigate to the folder where the package is and run e..g, R --vanilla CMD build bgev



# 3
# Check as CRAN:
# run R CMD check bgev_0.2.tar.gz --as-cran

# 4
# Upload the package to https://win-builder.r-project.org/upload.aspx and wait for report. Fix issues. 

# 5
# Last step: submission in R. 
# https://cran.r-project.org/submit.html

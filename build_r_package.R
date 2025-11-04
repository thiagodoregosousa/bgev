# building the r package
rm(list=ls())
code_files = c("package/bgev_functions.R","package/dist_check.R")
package_name = 'bgev'
path_for_package_files = 'bgev_cran'
package.skeleton(name=package_name, 
                 path=path_for_package_files, 
                 code_files = code_files)


# build package with R: after buiding, or if you have already build and are editing just code for new version next step is to build package
# on terminal navigate to the folder where the package is and run e..g, R --vanilla CMD build bgev


# Check as CRAN:
# run R CMD check bgev_0.2.tar.gz --as-cran
# bgev
bimodal gev distribution r code to compute density, distribution function, quantiles, random variables and estimation

# project organization

bgev_package_functions: functions that are going for the package (most recent version)

bgev_cran - files for the r package (the functions of bgev_cran/bgev/R are copied from bgev_package_functions). 

all other functions and files are tests and additional routines. Either for testing, 
or for future versions of this package. 



# Advices for building the package 

When using functions from other packages makesure to use them like DEoptim::DEoptim.control or stats::runif(n)
and put them  in imports like 

Imports:
    EnvStats,
    DEoptim,
    stats,
    MASS




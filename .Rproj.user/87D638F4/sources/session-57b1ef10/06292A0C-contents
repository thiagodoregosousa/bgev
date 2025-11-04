

################################
#Examples
# Density of a  bimodal generalized extreme value distribution with 
# mu=1, sigma=10, xi=0.3, delta=2){ 

dbgevd(.5,mu=1, sigma=10, xi=0.3, delta=2)
#[1] 0.02748573

# The 90'th percentile of a bimodal generalized extreme value distribution
qbgevd(.9,mu=1, sigma=10, xi=0.3, delta=2)
#[1] 3.538773

# Random sample of 4 observations from a bimodal generalized extreme value 
# distribution with mu=1, sigma=10, xi=0.3 and delta=2. 
# (Note: the call to set.seed simply allows you to reproduce this example.)

set.seed(20) 
rbgevd(4,mu=1, sigma=10, xi=0.3, delta=2) 
# [1]  3.4788506  3.2238194 -0.3633905  2.6166901



################################
#Examples

n=50
Z <- rdgevd(n, mu, sigma, xi, delta)
ebgevd(Z,R=10,mu=1,sigma=10,xi=0.3,delta=2,conf.level = 0.95)
#                parameter
#Estimative       mu                                     
#Ponctual       1.02279116583337                       
#Conf. Interval ( 1.01168448475556 ; 1.03389784691119 )
#parameter
#Estimative       sigma                                  
#Ponctual       16.9963296503858                       
#Conf. Interval ( 16.3867158696069 ; 17.6059434311648 )
#parameter
#Estimative       xi                                       
#Ponctual       0.139581984753769                        
#Conf. Interval ( 0.124364259639566 ; 0.154799709867972 )
#parameter
#Estimative       delta                                  
#Ponctual       2.49602823332784                       
#Conf. Interval ( 2.45331666498226 ; 2.53873980167342 )


library(dplyr)
library(mgcv)

#You need the multivariate normal RNG from the MASS package
mvrnorm = MASS::mvrnorm

# Functions to calculate the 1st and 2nd derivatives of a given time series with a given step size
# Uses two-point approximations for both 1st and 2nd derivs.
calc_1st_deriv = function(y,delta) (lead(y,1) - lag(y,1))/(2*delta)
calc_2nd_deriv = function(y,delta) (lead(y,1) + lag(y,1)-2*y)/delta^2

# Function to test if each the confidence intervals for the derivatives of a given 
# curve sastisfy the test for a point potentially being an extremum 
# The first deriv critera determines if the confidence intervals for the first 
# derivative at a given point overlaps zero. The second deriv criteria tests
# if the second deriv CI overlaps zero. A point is only considered to be a 
# candidate extremum if the first criteria is true and the second is false
find_candidate_peaks = function(deriv_1_bounds, deriv_2_bounds){
  deriv1_criteria = deriv_1_bounds[,1]<0&deriv_1_bounds[,2]>0
  deriv2_criteria  = !(deriv_2_bounds[,1]<0&deriv_2_bounds[,2]>0)
  is_candidate = deriv1_criteria&deriv2_criteria
  is_candidate
}


# function for generating underlying curve. Simple sum of a log-curve and a quadratic
# parameters a,b determine how large the log and quadratic components are, mid determines where the
# center of the quadratic is
fit_func = function(x, a,b,mid)  a*log(x)+b*(x-mid)^2 

n = 100         # number of data points
sigma = 1           # st. dev. of errors around mean value

# parameters for the true function
a   = 1      
b   = 1
mid = 2

#Parmaters determining how large the range of x is
low_lim = 0.1
high_lim =4


# data frame with x, values, random y values, and the true curve
training_data = data_frame(x = seq(low_lim,high_lim, length=n),y= rnorm(n,fit_func(x,a,b,mid),sigma),
                       true_val = fit_func(x,a,b,mid) )


# Extracts the actual minimum for the function (or the left boundary if the min
# is outside the data range)
true_min = optim(par = 1, fit_func, a=a,b=b,mid=mid,lower = low_lim,upper=high_lim,method = "L-BFGS-B")$par[1]


# The fitted model, usig a 20 basis function thin plate spline smoother with REML
# fitting criteria. m=3 specifies that the model should penalize squard third derivatives.
# This is important as if m=2 (the default) then prior simulations from the fit are too wiggly, and 
# end up with too wide a range of 2nd derivatives
mod = gam(y~s(x, bs= "tp", k=20,m = 3), data=training_data,method="REML")


# step size used for calculating derivatives
step_size = 0.01

# The test data, with one x per step unit across the range.
test_data = data_frame(x=seq(low_lim, high_lim,by= step_size))

#Simulate new functions from the posterior distribution of functions, using the
#test data and 500 simulations
n_sims = 500
mod_coef = coef(mod) # mean values for all basis functions
mod_vcov =vcov(mod)  # posterior variance-covariance matrix
mod_sims = mvrnorm(n_sims, mod_coef,mod_vcov) #random parameter draws
test_lp = predict.gam(mod,newdata = test_data,type = "lpmatrix") #the basis functions 
test_sims = test_lp %*% t(mod_sims) #random parameters times basis functions

#Calculates estimated first and second derivatives 
test_1st_deriv = apply(test_sims,MARGIN = 2,calc_1st_deriv, delta= step_size)
test_2nd_deriv = apply(test_sims,MARGIN = 2,calc_2nd_deriv, delta= step_size)

# 95% confidence intervals for the function, 1st, and 2nd derivatives
test_CI = t(apply(test_sims,
                  MARGIN = 1,
                  FUN = quantile,
                  probs=c(0.025,0.5,0.975),
                  na.rm=T))
test_1st_deriv_CI = t(apply(test_1st_deriv ,
                            MARGIN = 1,
                            FUN = quantile,
                            probs=c(0.025,0.5,0.975),
                            na.rm=T))
test_2nd_deriv_CI = t(apply(test_2nd_deriv ,
                            MARGIN = 1,
                            FUN = quantile,
                            probs=c(0.025,0.5, 0.975),
                            na.rm=T))

# Using the CIs for 1st and 2nd derivatives to test for peaks
candidate_peaks = as.vector(find_candidate_peaks(test_1st_deriv_CI[,c(1,3)], 
                                                 test_2nd_deriv_CI[,c(1,3)]))

candidate_peaks = ifelse(is.na(candidate_peaks), F, candidate_peaks)

#plotting curves ####


par(mfrow=c(3,1))

# Plot of raw data and model fit, with true function in blue and 
# estimated minima in red. Vertical blue line is the true minimum
plot(y~x, data= training_data)
points(true_val~x, data=training_data, col="blue",type="l",lwd=1)
matplot(test_data$x, test_CI,type="l",col="black",lty=c(2,1,2),add = T)
matplot(test_data$x[candidate_peaks], test_CI[candidate_peaks,],
        type="l",col="red",lty=c(2,1,2),add = T)
abline(v= true_min, col="blue",lty=2)

#plot of first derivatives plus CI
matplot(test_data$x,test_1st_deriv_CI,type="l",col="black",lty=c(2,1,2))
matplot(test_data$x[candidate_peaks], test_1st_deriv_CI[candidate_peaks,],
        type="l",col="red",lty=c(2,1,2),add = T)
abline(h=0,lty=3,col="red")

#Plot of estimated 2nd derivative plus CI
matplot(test_data$x,test_2nd_deriv_CI,type="l",col="black",lty=c(2,1,2))
matplot(test_data$x[candidate_peaks], test_2nd_deriv_CI[candidate_peaks,],
        type="l",col="red",lty=c(2,1,2),add = T)
abline(h=0, lty=3, col="red")


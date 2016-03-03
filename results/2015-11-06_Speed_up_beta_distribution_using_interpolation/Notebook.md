2015-11-06 Speed up beta distribution using interpolation
=============


Purpose: 
---------------------
Speed up the integration function for gamma distribution.

Conclusions
-----------------
Interpolation is not accurate enough, and sometimes it failed to optimize gamma or gamma is NAN because of "Method BFGS cannot handle constraints nor bounds".  
Sensitivity is very low because no variants are called. Specificity is very good.

	
Background
-----------------


Materials and Equipment
----
**1. Try cubic spline interpolation to approximate the integral.**  
[http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html](http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html)

	The maximal number of iterations maxit (set to 20 by the program)
	allowed for finding a smoothing spline with fp=s has been reached: s
	too small.
	There is an approximation returned but the corresponding weighted sum
	of squared residuals does not satisfy the condition abs(fp-s)/s < tol.
	  warnings.warn(message)

Fixed by:

	def EqlogGamma(gam, M): 
	    mu = np.arange(1e-3, 1-1e-3, 0.001)
	    tck = splrep(mu, kernel(mu, gam, M), s=0)
	    munew = [0, 1] 
	    temp = integ(munew, tck) 
	
	    integr = temp[-1]
	    integr = integr.astype(float)

	    logGamma = integr.astype(type('float', (float,), {}))
	    return logGamma

[https://github.com/scipy/scipy/issues/3691](https://github.com/scipy/scipy/issues/3691)

There could be an error " Method BFGS cannot handle constraints nor bounds."  
"WARNING:rvd3:Could not optimize gamma or gamma is NaN."


2.Faster integration using Ctypes:  
[http://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html](http://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html)



Results
-----------------
Sensitivity is not good because that no variants are called in some events.

Archived Samples
-------------------------


Archived Computer Data
------------------------------


Prepared by: _________ _Fan Zhang______ Date: ____________11/09/2014____________


Witnessed by: ________________________
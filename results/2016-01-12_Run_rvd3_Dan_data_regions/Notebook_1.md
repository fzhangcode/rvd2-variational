2016-01-12 Run rvd3 Dan data regions
=============


Purpose: 
---------------------
Run rvd3 (variational) on MTH1 from E1 on Dan's data.

Conclusions
-----------------
RVD3 gives good estimation of variant allele frequencies and cal variants.

Length of region and max iterations both effect results.

Background
-----------------
We focus on gene MTH1 (region: ChrIV:1014401-1015702, 1300 positions, ~2000 read depth) to see if there are new variants to be found.

Results of MCMC are in folder \fzhang\Research\rvd2\results\2014-06-17_call_mutations_in_MTH1_varying_case  
Results of calling concomitant mutations are in folder \fzhang\Research\rvd2\results\2014-06-26_call_concomitant_mutations


Materials and Equipment
----

Try different initial values and stop criterion for delta and gamma.  

###1. For gene MTH1 in E1:

Define optimization stopping criterion 

### Try different initial values: results (variant allele frequency) don't match Dan's results.
 
seed = 20160114

    MAXITER = 50
    ELBOTOLPCT = 0.001 *100
    MAXVARITER = 50     
    NORMTOL = 0.1

seed = 100  

    MAXITER = 100
    ELBOTOLPCT = 0.001 *100
    MAXVARITER = 100     
    NORMTOL = 0.1

### Different stop iterations for E-step and M-step give different ELBO: results are still not correct. But different stop MAXITER = 15 gives higher ELBO than MAXITER = 150, MAXITER = 20.

seed = 198605220525 -> c4-1_MTH1_198605220525_1.hdf5

    MAXITER = 150
    ELBOTOLPCT = 0.0001 *100
    MAXVARITER = 150     
    NORMTOL = 0.1


seed = 198605220525 -> c4-1_MTH1_198605220525_2.hdf5

    MAXITER = 20
    ELBOTOLPCT = 0.01 *100
    MAXVARITER = 20     
    NORMTOL = 0.1

seed = 198605220525 -> c4-1_MTH1_198605220525_3.hdf5

    MAXITER = 15
    ELBOTOLPCT = 0.01 *100
    MAXVARITER = 15     
    NORMTOL = 0.1


### Try different initial value intervals using [0.01, 100] instead of [0.0001, 10000]: results are still not corret. The results of these two scenario give us similar ELBO.

seed = 198605220525 -> c4-34_MTH1_198605220525_1.hdf5

    ## Define optimization stopping criterion
    MAXITER = 30
    ELBOTOLPCT = 0.01 *100
    MAXVARITER = 30    
    NORMTOL = 0.1

    delta = np.random.uniform(low = 0.01, high = 100, size = (N,J,2))
    gam = np.random.uniform(low=0.01, high=100, size = (J,2))

seed = 198605220525 -> c4-34_MTH1_198605220525_2.hdf5

    ## Define optimization stopping criterion
    MAXITER = 30
    ELBOTOLPCT = 0.001 *100
    MAXVARITER = 30    
    NORMTOL = 0.1

    delta = np.random.uniform(low = 0.01, high = 100, size = (N,J,2))
    gam = np.random.uniform(low=0.01, high=100, size = (J,2))


###Take first 400 positions to generate a new sequence: results are good.

   ## Define optimization stopping criterion

    MAXITER = 50
    ELBOTOLPCT = 0.001 *100
    MAXVARITER = 50
    NORMTOL = 0.1

    delta = np.random.uniform(low = 0.0001, high = 10000, size = (N,J,2))
    gam = np.random.uniform(low=0.0001, high=10000, size = (J,2))


seed = 1986 -> c4-1_MTH1_top400_1986.hdf5 -> looks correct from plotting the posterior distribution. The vaf looks correct of position 1014707 and 1014740. And also, rvd2 and rvd3 both call position 1014430 and 1014651 which are not identified by Dan. 

seed = 1986 -> c4-34_MTH1_top400_1986.hdf5  
seed = 20160201 -> c4-34_MTH1_top400_20160201.hdf5  
seed = 5 -> c4-34_MTH1_top400_5.hdf5

seed = 1986 -> c4-22_MTH1_top400_1986.hdf5
seed = 20160202 -> c4-22_MTH1_top400_20160202.hdf5 It works good that it detected many mutations that are the same as MCMC results.  

Try different initial values for c4-34_MTH1_top400, but they fail to estimate correct vaf for position 1014740. Dan reveals it is 0.983.


seed = 20160201 -> c4-42_MTH1_top400_20160201.hdf5 using less MAXITER to see if the ELBO is not over zero any more. It doesn't work.

    ## Define optimization stopping criterion
    MAXITER = 50
    ELBOTOLPCT = 0.001 *100
    MAXVARITER = 1    
    NORMTOL = 0.1

seed = 10 -> c4-15_MTH1_top400_10.hdf5 ?

    ## Define optimization stopping criterion
    MAXITER = 20
    ELBOTOLPCT = 0.001 *100
    MAXVARITER = 20    
    NORMTOL = 0.1
 
####Different seeds may cause different results.  
#### Pat: I don't think there's anything that says that ELBO has to be greater or less than zero. The ELBO is a lower bound on the log-likelihood. The likelihood has to be greater than zero so the log-likelihood is well-defined. But there isn't any other constraints and no upper bound as far as I know.

### For gene MTH1 and gene ADE16 in E1:
Depth chart files are from X:\fzhang\Research\rvd2\results\2014-06-26_call_concomitant_mutations


### Study position 1014740  
--
VAF > 90% at generation 322, 385, and 448.  
It might help to take a look at the parameters in detail at that position. In particular, what does the estimate of theta for each replicate look like?


seed = 20160201 -> c4-34_MTH1_top400_20160201.hdf5 

   ## Define optimization stopping criterion

    MAXITER = 50
    ELBOTOLPCT = 0.001 *100
    MAXVARITER = 50
    NORMTOL = 0.1

    delta = np.random.uniform(low = 0.0001, high = 10000, size = (N,J,2))
    gam = np.random.uniform(low=0.0001, high=10000, size = (J,2))

For the wrong estimated vaf for position 1014740, M0 is too big.
Two solutions: 
1. add a gamma prior on M0.
2. Divide M0 by 10 or 100 or 1000 or 10000 to see if the estimate of muj is getting correct.
The reason is that if we have more positions, M0 could get big, variance gets small.

It gives good estimation of VAF on position 1014740 when I tuned Mo to M0/100.

Results
-----------------


Archived Samples
-------------------------
 


Archived Computer Data
------------------------------


Prepared by: _________ _Fan Zhang______ Date: ____________01/25/2016____________


Witnessed by: ________________________


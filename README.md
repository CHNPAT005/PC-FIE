# Fourier instantaneous estimators and the Epps effect

## Authors:
- Patrick Chang

## Link to resources:

Link to paper: 

Link to the Dataset on ZivaHub: [Link](https://zivahub.uct.ac.za/articles/Malliavin-Mancino_estimators_implemented_with_the_non-uniform_fast_Fourier_transform_Dataset/11903442)

## Steps for Replication:
- Change directories for all the files under [/Scripts/Instantaneous](https://github.com/CHNPAT005/PC-FIE/tree/master/Scripts/Instantaneous). Currently the directories are set as: `cd("/Users/patrickchang1/PC-FIE")`. Change this to where you have stored the file `PC-FIE`. 

	- Run [/Scripts/Instantaneous/Synchronous](https://github.com/CHNPAT005/PC-FIE/blob/master/Scripts/Instantaneous/Synchronous) to reproduce Fig. 1.
	- Run [/Scripts/Instantaneous/Asynchronous](https://github.com/CHNPAT005/PC-FIE/blob/master/Scripts/Instantaneous/Asynchronous) to reproduce Fig. 2.
	- Run [/Scripts/Instantaneous/InstantaneousEpps](https://github.com/CHNPAT005/PC-FIE/blob/master/Scripts/Instantaneous/InstantaneousEpps) to reproduce Fig. 3 and 4.
	- Run [/Scripts/Instantaneous/OptimalTimeScale](https://github.com/CHNPAT005/PC-FIE/blob/master/Scripts/Instantaneous/OptimalTimeScale) to reproduce Fig. 5 and 6.

- To reproduce the Empirical analysis - download the processed dataset from ZivaHub and put the csv files into the folder `/Real Data`.
	- Run [/Scripts/Instantaneous/Empirical_Inst](https://github.com/CHNPAT005/PC-FIE/blob/master/Scripts/Instantaneous/Empirical_Inst) to reproduce Fig. 7-10. Note that the paper only reports surface and contour plots for 2019-06-25 and 2019-06-26, but the script file makes the plots for the week from 2019-06-24 to 2019-06-28.
	

## Using the functions for other purposes:
Some functions here include the 1 Dimensional Type 1 non-uniform fast Fourier transform and the Malliavin-Mancino integrated estimates. Instructions for using these functions can be found in previous work: https://github.com/CHNPAT005/PCEPTG-MM-NUFFT

### Cuchiero-Teichmann instantaneous estimator

The Cuchiero-Teichmann instantaneous estimator uses the specification of g(x) = cos(x). The estimator requires the data to be strictly synchronous, therefore asynchronous data needs to be synchronised beforehand using the previous tick interpolation.

The function requires 3 input variables:
- p: (n x 2) double float matrix of price observations.
- N: the cutting frequency (integer) used in the reconstruction of the spot estimates, the paper uses the notation M.
- outlength: the number of synchronous grid points to reconstruct the spot estimates.

#### Example

```julia

include("Functions/SDEs/Heston")
include("Functions/Instantaneous Estimators/MM-JR")

# Simulate some price observations from the Heston model.

nsim = 28800
P_Heston = Heston_CT(nsim, seed = 1, dt = nsim)
	# First variable in P_Heston is the price matrix, 
	# Second to fourth variable are the true volatility and co-volatility.

# Parameter settings
outlength = 1000	# length of output vector
M = 100	# Cutting freq.

# Output 
JR_Heston = MM_JR(P_Heston[1], M, outlength)
	# First variable is the volatility estimates of asset 1.
	# Second variable is the volatility estimates of asset 2.
	# Third variable is the co-volatility estimates of asset 1 and 2.

# Plot the results
tt = collect(0:1/outlength:1-(1/outlength))
t = t ./ nsim

p1 = plot(t, P_Heston[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)	# True volatility of asset 1.
plot!(p1, tt, JR_Heston[1], label = L"\textrm{Estimated CT } \hat{\Sigma}^{11}_{n,M}(t)", color = :red, line=(1, [:dash]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")
```

### Malliavin-Mancino instantaneous estimator

The Malliavin-Mancino instantaneous estimator uses the non-uniform fast Fourier transform to speed up the implementation. The estimator can deal with asynchronous data, therefore the observations require an associated time.

The function requires 2 input variables:
- p: (nx2) matrix of prices, with non-trade times represented as NaNs.
- t: (nx2) corresponding matrix of trade times, with non-trade times represented as NaNs.
and three optional input variables.
- N: (optional input) for the number of Fourier coefficients of the price process used in the convolution of the Malliavin-Mancino estimator (integer), controls the level of averaging and directly affects the time-scale investigated - defaults to the Nyquist frequency.
- M: (optional input) for the number of Fourier coefficients of the volatility process using in the reconstruction of the spot estimates - defaults to $`M = \frac{1}{8} \frac{1}{2\pi} \sqrt{n} \log n`$
- tol: tolerance requested - for the NUFFT implementations - defaults to 10^-12.



#### Example

```julia

include("Functions/SDEs/Heston")
include("Functions/Instantaneous Estimators/MM-Inst")

# Simulate some price observations from the Heston model.

nsim = 28800
P_Heston = Heston_CT(nsim, seed = 1, dt = nsim)
	# First variable in P_Heston is the price matrix, 
	# Second to fourth variable are the true volatility and co-volatility.
t = collect(1:1:nsim)

# Parameter settings
outlength = 1000	# length of output vector
M = 100	# Cutting freq.

# Output 
MM_Heston = MM_inst(P_Heston[1], [t t], outlength, M = M)
	# First variable is the volatility estimates of asset 1.
	# Second variable is the volatility estimates of asset 2.
	# Third variable is the co-volatility estimates of asset 1 and 2.

tt = collect(0:1/outlength:1-(1/outlength))
t = t ./ nsim

p1 = plot(t, P_Heston[2], label = L"\textrm{True } \Sigma^{11}(t)", color = :lightblue, line=(0.5, [:solid]), dpi = 300)
plot!(p1, tt, MM_Heston[1], label = L"\textrm{Estimated MM } \hat{\Sigma}^{11}_{n,N,M}(t)", color = :blue, line=(1, [:solid]))
xlabel!(p1, L"\textrm{Time}")
ylabel!(p1, L"\textrm{Instantaneous variance of } \Sigma^{11}(t)")

```
	

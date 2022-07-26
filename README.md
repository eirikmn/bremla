
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bremla

<!-- badges: start -->

[![R-CMD-check](https://github.com/eirikmn/bremla/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/eirikmn/bremla/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This Repository contains the Bremla package for uncertainty
quantification and synchronization for layer-counted proxy archives,
including the data and code to reproduce the results for

Myrvoll-Nilsen, E., Riechers, K., Rypdal, M. & Boers, N. (2022).
Comprehensive uncertainty estimation of the timing of Greenland warmings
of the Greenland Ice core records. Climate of the Past, 18, 1275-1294.
doi.org/10.5194/cp-18-1275-2022

and

Myrvoll-Nilsen, E., Riechers, K. & Boers, N. (202x). Tie-point paper
(title TBD)

## installation

The simplest way to install the package is to install the remotes
package and run

``` r
#install.packages("devtools")
devtools::install_github("eirikmn/bremla")
```

To get started, see the documentation ‘?bremla’ and the associated
example.

The package includes the NGRIP/GICC05 data set ‘NGRIP_d18O_and_dust_5cm’
downloaded from [Centre for Ice and Climate](iceandclimate.nbi.ku.dk) at
the Niels Bohr Institute of the University of Copenhagen, Denmark, and
the table of stadial-interstadial events presented by [Rasmussen et
al. (2014)](https://www.sciencedirect.com/science/article/pii/S0277379114003485).
Also includes the tie-point presented in [Adolphi et
al. (2018)](https://cp.copernicus.org/articles/14/1755/2018/) and
[Muscheler et
al. (2020)](https://www.cambridge.org/core/journals/radiocarbon/article/testing-and-improving-the-intcal20-calibration-curve-with-independent-records/D72D9214C47FE9441B5E730D33DCCE3D).
Other data and tie-points should work as well.

Warning: The provided real data example generates several thousand
simulated chronologies with individual lengths exceeding 18,000. More if
synchronized chronologies is also computed. This will require a
significant amount of memory. If needed, reduce the number of samples
generated or free up memory before use.

## Example

This is a real data example which produces 5000 simulations from a
synchronized time scale.

``` r
require(INLA)
library(bremla)
data("event_intervals")
data("events_rasmussen")
data("NGRIP_5cm")

age = NGRIP_5cm$age
depth = NGRIP_5cm$depth
proxy = NGRIP_5cm$d18O

data = data.frame(age=age,dy=c(NA,diff(age)),depth=depth,depth2=depth^2,proxy=proxy)
formula = dy~-1+depth2

events=list(locations = events_rasmussen$depth,
            locations_unit="depth",degree=1)

object = bremla(formula,data,reference.label="GICC05",
                nsims=5000,
                events=events,
                synchronization=list(method="adolphi"),
                control.fit=list(method="inla"),
                control.sim=list(synchronized=TRUE),
                print.progress=TRUE )
#> Initiating data formatting... completed!
#> Performing least squares fit... completed!
#> Performing INLA fit...
#> Computing remaining posteriors using Monte Carlo simulation...
#> INLA fit completed in 22.77553 seconds!
#> Simulating 5000 hyperparameters from INLA posterior... completed!
#> Sampling fixed coefficients... completed in 34.65476 seconds.
#> Simulating synchronized chronologies...
#> Synchronous age simulation 1000/5000. Elapsed time: 16.81051 seconds...
#> Synchronous age simulation 2000/5000. Elapsed time: 31.92003 seconds...
#> Synchronous age simulation 3000/5000. Elapsed time: 47.00688 seconds...
#> Synchronous age simulation 4000/5000. Elapsed time: 61.99271 seconds...
#> Synchronous age simulation 5000/5000. Elapsed time: 78.27454 seconds...
#> Computing posterior marginal mean and 95% credible intervals from chronology samples...
#>  completed in 2.852961seconds.
summary(object)
#> 
#> Call:
#> render("README.Rmd")
#> 
#> Time used:
#>   Model fitting Chron. sampling           Total 
#>         23.1660        112.9717        139.1153 
#> 
#> The fixed component is explained by linear predictor: 
#> dy ~ -1 + depth2 + psi_fill(degree=1, n_events=69)
#> 
#> The noise component is explained by an ar1 process.
#> 
#> The model is fitted using INLA, with following estimates for the hyperparameters:
#>                 mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#> sigma_epsilon 0.4348 0.0025     0.4298    0.4330   0.4347    0.4365     0.4397
#> phi           0.2159 0.0067     0.2027    0.2113   0.2159    0.2204     0.2289
#> 
#> Simulating 5000 chronologies, using GICC05 as reference.
#> 
#> 5000 synchronized chronologies sampled using 4 tie-point distributions (Adolphi).
#> Tie-points are fixed at GICC05 ages (yb2k):
#> 11050, 12050, 13050, 22050, 42050.
```

<img src="man/figures/README-plot-1.png" width="80%" />

## Attribution

This code is associated and written for the papers Myrvoll-Nilsen et
al. (2022) and Myrvoll-Nilsen et al. (202x) mentioned above. Feel free
to use the code, but please cite the accompanying papers.

## License

The code in this repository is made available under the terms of the GNU
(version \>=2) License. For details, see LICENSE.md file.

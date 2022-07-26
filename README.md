
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

<img src="man/figures/README-plot-1.png" width="80%" /> \## Attribution

This code is associated and written for the papers Myrvoll-Nilsen et
al. (2022) and Myrvoll-Nilsen et al. (202x) mentioned above. Feel free
to use the code, but please cite the accompanying papers.

## License

The code in this repository is made available under the terms of the GNU
(version \>=2) License. For details, see LICENSE.md file.



# Bremla
This Repository contains the Bremla package for uncertainty quantification and synchronization for layer-counted proxy archives, including the data and code to reproduce the results for

Myrvoll-Nilsen, E., Riechers, K., Rypdal, M. & Boers, N. (2022). Comprehensive uncertainty estimation of the timing of Greenland warmings of the Greenland Ice core records. Climate of the Past, x(x),x-x. doi:xxxxx

and 

Myrvoll-Nilsen, E., Riechers, K. & Boers, N. (202x). Tie-point paper (title TBD)

## installation

The simplest way to install the package is to install the remotes package and run

```r
#install.packages("remotes")
remotes::install_github("eirikmn/bremla")
```

To get started, see the documentation '?bremla' and the associated example.

The package includes the NGRIP/GICC05 data set 'NGRIP_d18O_and_dust_5cm' downloaded from [iceandclimate.nbi.ku.dk](Centre for Ice and Climate) at the Niels Bohr Institute of the University of Copenhagen, Denmark, and the table of stadial-interstadial events presented by [https://www.sciencedirect.com/science/article/pii/S0277379114003485](Rasmussen et al. (2014)). Also includes the tie-point presented in [https://cp.copernicus.org/articles/14/1755/2018/](Adolphi et al. (2018)) and [https://www.cambridge.org/core/journals/radiocarbon/article/testing-and-improving-the-intcal20-calibration-curve-with-independent-records/D72D9214C47FE9441B5E730D33DCCE3D](Muscheler et al. (2020)). Other data and tie-points should work as well.

Warning: The provided real data example generates several thousand simulated chronologies with individual lengths exceeding 18,000. More if synchronized chronologies is also computed. This will require a significant amount of memory. If needed, reduce the number of samples generated or free up memory before use.

## Example
This is a real data example which produces 5000 simulations from a synchronized time scale.

```r
require(INLA)
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
control.transition_dating=list(age.ref=age.reference,age.label="years (yb2k)")
```

```
## Error in eval(expr, envir, enclos): object 'age.reference' not found
```

```r
object = bremla(formula,data,reference.label="GICC05",
                  nsims=5000,
                  events=events,
                  synchronization=list(method="adolphi"),
                  control.fit=list(method="inla"),
                  control.sim=list(synchronized=TRUE),
                  print.progress=TRUE )
```

```
## Initiating data formatting... completed!
## Performing least squares fit... completed!
## Performing INLA fit...
## Computing remaining posteriors using Monte Carlo simulation...
## INLA fit completed in 39.79978 seconds!
```

```r
summary(object)
```

```
## 
## Call:
## knit(input = "readme.rmd", output = "readme.md")
## 
## Time used:
## Model fitting         Total 
##       40.0877       40.3451 
## 
## The fixed component is explained by linear predictor: 
## dy ~ -1 + depth2 + psi0_1 + psi1_1 + psi0_2 + psi1_2 + psi0_3 +     psi1_3 + psi0_4 + psi1_4 + psi0_5 + psi1_5 + psi0_6 + psi1_6 +     psi0_7 + psi1_7 + psi0_8 + psi1_8 + psi0_9 + psi1_9 + psi0_10 +     psi1_10 + psi0_11 + psi1_11 + psi0_12 + psi1_12 + psi0_13 +     psi1_13 + psi0_14 + psi1_14 + psi0_15 + psi1_15 + psi0_16 +     psi1_16 + psi0_17 + psi1_17 + psi0_18 + psi1_18 + psi0_19 +     psi1_19 + psi0_20 + psi1_20 + psi0_21 + psi1_21 + psi0_22 +     psi1_22 + psi0_23 + psi1_23 + psi0_24 + psi1_24 + psi0_25 +     psi1_25 + psi0_26 + psi1_26 + psi0_27 + psi1_27 + psi0_28 +     psi1_28 + psi0_29 + psi1_29 + psi0_30 + psi1_30 + psi0_31 +     psi1_31 + psi0_32 + psi1_32 + psi0_33 + psi1_33 + psi0_34 +     psi1_34 + psi0_35 + psi1_35 + psi0_36 + psi1_36 + psi0_37 +     psi1_37 + psi0_38 + psi1_38 + psi0_39 + psi1_39 + psi0_40 +     psi1_40 + psi0_41 + psi1_41 + psi0_42 + psi1_42 + psi0_43 +     psi1_43 + psi0_44 + psi1_44 + psi0_45 + psi1_45 + psi0_46 +     psi1_46 + psi0_47 + psi1_47 + psi0_48 + psi1_48 + psi0_49 +     psi1_49 + psi0_50 + psi1_50 + psi0_51 + psi1_51 + psi0_52 +     psi1_52 + psi0_53 + psi1_53 + psi0_54 + psi1_54 + psi0_55 +     psi1_55 + psi0_56 + psi1_56 + psi0_57 + psi1_57 + psi0_58 +     psi1_58 + psi0_59 + psi1_59 + psi0_60 + psi1_60 + psi0_61 +     psi1_61 + psi0_62 + psi1_62 + psi0_63 + psi1_63 + psi0_64 +     psi1_64 + psi0_65 + psi1_65 + psi0_66 + psi1_66 + psi0_67 +     psi1_67 + psi0_68 + psi1_68 + psi0_69 + psi1_69
## 
## The noise component is explained by an ar1 process.
## 
## The model is fitted using INLA, with following estimates for the hyperparameters:
##                  mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
## sigma_epsilon  4.4274 0.0645     4.2872    4.3857   4.4336    4.4744     4.5370
## phi           -0.9884 0.0003    -0.9889   -0.9886  -0.9884   -0.9882    -0.9876
```

```r
plot(object)
```

## Attribution

This code is associated and written for the papers Myrvoll-Nilsen et al. (2022) and Myrvoll-Nilsen et al. (202x) mentioned above. Feel free to use the code, but please cite the accompanying papers.

## License

The code in this repository is made available under the terms of the MIT License. For details, see LICENSE file.

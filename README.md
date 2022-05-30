This Repository contains the Bremla package for uncertainty quantification and synchronization for layer-counted proxy archives, including all the data/code to reproduce the results for

Myrvoll-Nilsen, E., Riechers, K., Rypdal, M. & Boers, N. (2022). Comprehensive uncertainty estimation of the timing of Greenland warmings of the Greenland Ice core records. Climate of the Past, x(x),x-x. doi:xxxxx

and 

Myrvoll-Nilsen, E., Riechers, K. & Boers, N. (202x). Tie-point paper (title TBD)

The simplest way to install the package is to install devtools and run 'devtools::install_github("eirikmn/bremla")'. To get started, see the documentation '?bremla' and the associated example.

The package includes the NGRIP/GICC05 data set 'NGRIP_d18O_and_dust_5cm' downloaded from iceandclimate.nbi.ku.dk, table of stadial-interstadial events (Table 2 from Rasmussen et al. (2014). Also includes the tie-point presented in Adolphi et al. (2018) and Muscheler et al. (2020). Other data and tie-points should work as well.

Warning: The default example generates several thousand simulated chronologies with individual lengths exceeding 18,000. More if synchronized chronologies is also computed. This will require a significant amount of memory. If needed, reduce the number of samples generated or free up memory before use.

Attribution

This code is associated and written for the papers Myrvoll-Nilsen et al. (2022) and Myrvoll-Nilsen et al. (202x) mentioned above. Feel free to use the code, but please cite the accompanying papers.

License

The code in this repository is made available under the terms of the MIT License. For details, see LICENSE file.
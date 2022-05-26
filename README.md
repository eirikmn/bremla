This Repository contains all the data/code to reproduce the results for

Myrvoll-Nilsen, E., Riechers, K., Rypdal, M. & Boers, N. (2022). Comprehensive uncertainty estimation of the timing of Greenland warmings of the Greenland Ice core records. Climate of the Past, x(x),x-x. doi:xxxxx

and 
Myrvoll-Nilsen, E., Riechers, K. & Boers, N. (202x). Tie-point paper (title TBD)

To use this code, install all files in the same directory. Use function main located in main.R as a starting point. All other functions are loaded from there.

The code should run a Bayesian regression analysis on the GICC05 time scale (NGRIP_d18O_and_dust_5cm.xls) from iceandclimate.nbi.ku.dk, but should be able to analyse other layer-counted proxy records as well. Updated to nclude framework for generating synchronized chronologies. Also includes tie-points from Adolphi et al. (2018) and Muscheler et al. (2020).

Warning: This code produces (by default) 30,000 simulated chronologies with individual lengths exceeding 18,000. Double if synchronized chronologies should also be computed. This will require a significant amount of memory. If needed, decrease the number of samples down or free up memory before use.

Attribution

This code is associated and written for the papers Myrvoll-Nilsen et al. (2022) and Myrvoll-Nilsen et al. (202x) mentioned above. Feel free to use the code, but please cite the accompanying papers.

License

The code in this repository is made available under the terms of the MIT License. For details, see LICENSE file.
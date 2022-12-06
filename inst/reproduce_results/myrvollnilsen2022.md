---
output: github_document
---

<!-- myrvollnilsen2022.md is generated from myrvollnilsen2022.Rmd. Please edit that file -->



# bremla

This file includes the data and code needed to reproduce the results for

Myrvoll-Nilsen, E., Riechers, K., Rypdal, M. & Boers, N. (2022). Comprehensive uncertainty estimation of the timing of Greenland warmings of the Greenland Ice core records. Climate of the Past, 18, 1275-1294. doi.org/10.5194/cp-18-1275-2022

# installation

The simplest way to install the package is to install the devtools package and run

```r
#install.packages("devtools")
devtools::install_github("eirikmn/bremla")
```


# Myrvoll-Nilsen et al. (2022)

The data used is the NGRIP/GICC05 data set 'NGRIP_d18O_and_dust_5cm', downloaded from [Centre for Ice and Climate](iceandclimate.nbi.ku.dk) at the Niels Bohr Institute of the University of Copenhagen, Denmark, and the table of stadial-interstadial events presented by [Rasmussen et al. (2014)](https://www.sciencedirect.com/science/article/pii/S0277379114003485). These are included in the package and can be imported as follows. 

```r
library(bremla)
library(dplyr)
library(ggplot2)
data("event_intervals")
data("events_rasmussen")
data("NGRIP_5cm")
age = NGRIP_5cm$age
depth = NGRIP_5cm$depth
depth2 = depth^2/depth[1]^2 #normalize for stability
proxy = NGRIP_5cm$d18O
n=length(diff(age))
data = data.frame(age=age,dy=c(NA,diff(age)),depth=depth,depth2=depth2,proxy=proxy)
events=list(locations = events_rasmussen$depth,
            locations_unit="depth",degree=1)


 plot(depth,proxy,xlab="Age (yb2k)",ylab=expression(paste(d^18,"O (permil)")),type="l",
      xlim=rev(range(depth)))+abline(v=events$locations,col="gray",lwd=0.7)
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-bremla-1.png" title="plot of chunk bremla" alt="plot of chunk bremla" width="100%" />

```
#> integer(0)
formula = dy~-1+depth2+proxy
nsims=3000
results = bremla(formula,data,reference.label="GICC05",
                nsims=nsims,
                events=events,
                control.fit=list(noise="ar1"),
                control.sim=list(synchronized=FALSE) 
                )

```

The results from the linear fit can be extracted by running:


```r
fit = results$fitting$LS$fit
{
  par(mfrow=c(2,1),mar=c(5,4,4,2)+0.1)
plot(results$data$depth,results$data$dy,type="l",xlab="Depth (m)",ylab="Layers per 5cm",main="(a) Least squares fit",xlim=rev(range(results$data$depth))) + lines(results$data$depth,fit$fitted.values,col="red")   
abline(v=events$locations,col="gray",lwd=0.8)
legend(1770,6.35,legend=c("Layer increments","Regression fit","Climate transitions"),
       col=c("black","red","gray"),lty=c(1,1,NA),cex=0.9,pch=c(NA,NA,"|"))
plot(results$data$depth,fit$residuals,type="l",xlab="Depth (m)",ylab="Residual errors per 5cm",main="(b) Least squares residuals",xlim=rev(range(results$data$depth)))
}
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-ls-1.png" title="plot of chunk ls" alt="plot of chunk ls" width="100%" />

```r
{
layout(mat=matrix(c(1,3,2,3),nrow=2))
  par(mar=c(5,3.5,4,2)+0.1)
hist(fit$residuals,xlab="Residual errors per 5cm",ylab="Density",main="(a) Histogram", col="orange")
{qqnorm(fit$residuals,main="(b) Q-Q Plot")
qqline(fit$residuals)}
acf(fit$residuals, lag.max = 20,main="(c) Autocorrelation function")
}
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-ls-2.png" title="plot of chunk ls" alt="plot of chunk ls" width="100%" />

```r
par(mfrow=c(1,1))
```

To compare the final age-depth uncertainties across the iid, AR(1) and AR(2) models run:

```r
resultsiid = bremla(formula,data,reference.label="GICC05",
                nsims=nsims,
                events=events,
                control.fit=list(noise="iid"),
                control.sim=list(synchronized=FALSE) 
                )
resultsar2 = bremla(formula,data,reference.label="GICC05",
                nsims=nsims,
                events=events,
                control.fit=list(noise="ar2"),
                control.sim=list(synchronized=FALSE) 
                )
gicc05 = results$original.chron$age[1+1:n]
{
plot(depth[1+1:n],resultsar2$simulation$summary$upper-gicc05,type="l",col="red", xlab="Depth (m)", ylab="Simulated time scale - GICC05 (years)",ylim=c(-260,260),xlim=rev(range(depth)), lwd=2)
lines(depth[1+1:n],results$simulation$summary$upper-gicc05,col="blue",lwd=2)
lines(depth[1+1:n],resultsiid$simulation$summary$upper-gicc05,col="black",lwd=2)
lines(depth[1+1:n],results$simulation$summary$mean-gicc05,col="gray",lwd=2)
abline(h=0,lty=3)
lines(depth[1+1:n],results$simulation$summary$lower-gicc05,col="blue",lwd=2)
lines(depth[1+1:n],resultsiid$simulation$summary$lower-gicc05,col="black",lwd=2)
lines(depth[1+1:n],resultsar2$simulation$summary$lower-gicc05,col="red",lwd=2)
legend(1640,275,legend=c("Mean","iid CI","AR(1) CI", "AR(2) CI"),
       col=c("gray","black","blue","red"),lty=c(1,1,1,1),cex=0.75)
}
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-iid_ar2-1.png" title="plot of chunk iid_ar2" alt="plot of chunk iid_ar2" width="100%" />

To generate the plot to compare different parameters for additional stochastic bias use the following code.


```r
results_biased = bremla(formula,data,reference.label="GICC05",
                nsims=nsims,
                events=events,
                control.fit=list(noise="ar1"),
                control.sim=list(synchronized=FALSE),
                control.bias=list(bias.model="uniform",nsims=nsims,
                                  biasparams = matrix(c(0.98,1.02,0.96,1.04),nrow=2)
                                  )
                )
{
  plot(depth[1+1:n],results_biased$biases$bias1$quant0.975-gicc05,type="l",col="blue",lty=2, xlab="depth", ylab="Simulated time scale - GICC05 (years)",xlim = rev(range(depth)), ylim = c(-2700,2700))
  lines(depth[1+1:n], results_biased$biases$bias2$quant0.975-gicc05,col="blue",lty=3)
  lines(depth[1+1:n], results$simulation$summary$upper-gicc05,col="blue",lty=1)
  abline(h=0,lty=3,col="gray")
  lines(depth[1+1:n], results$simulation$summary$lower-gicc05,col="blue",lty=1)
  lines(depth[1+1:n], results_biased$biases$bias1$quant0.025-gicc05,col="blue",lty=2)
  lines(depth[1+1:n], results_biased$biases$bias2$quant0.025-gicc05,col="blue",lty=3)
  lines(depth[1+1:n], NGRIP_5cm$MCE[1+1:n],col="black",lwd=2)
  legend(1660,2700,legend=c("MCE","Unbiased","2% bias rate", "4% bias rate"),
       col=c("black","blue","blue","blue"),lty=c(1,1,2,3),lwd=c(2,1,1,1),cex=0.65)
}
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-biased-1.png" title="plot of chunk biased" alt="plot of chunk biased" width="100%" />



## Abrupt warming transitions

We also import table D1 from Myrvoll-Nilsen et al. (2022) which gives suitable data windows for each abrupt warming transition that is analysed. This was found by varying the start and end point over a 2-dimensional grid and picking those where the noise in the linear ramp model yielded the lowest amplitude.





```r
eventnumber=13 #number between 1 and 29. specifies which transition to consider
#load data window and specifics to transition
lowerints = which.index(event_intervals$depth_int_lower.m, depth[2:length(depth)])
upperints = which.index(event_intervals$depth_int_upper.m, depth[2:length(depth)])
depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
age.reference = event_intervals$GICC_age.yb2k[eventnumber]
interval = lowerints[eventnumber]:upperints[eventnumber]
transitionlabel = str_sub(event_intervals$onsetlabel[eventnumber],
                      str_locate(event_intervals$onsetlabel[eventnumber],"GI")[1])
control.linramp=list(proxy=proxy,
                     interval=interval,
                     interval.unit="index",
                     depth.reference=event_intervals$NGRIP_depth_m[eventnumber]
                     )
results_event = linrampfitter(results,control.linramp)
control.transition_dating = list(label="Simulated transition")
results_event = events_depth_to_age(results_event,
                                    control.transition_dating=control.transition_dating)
#> Warning in events_depth_to_age(results_event, control.transition_dating =
#> control.transition_dating): Could not find synchronized chronology samples. Using
#> unsynchronized instead.
{
  layout(mat=matrix(c(1,3,2,3),nrow=2))
  
  
  plot(depth[interval],proxy[interval],type="l",col="gray",
     main=paste0("(a) ",transitionlabel," linear ramp fit"),
     xlab="Depth (m)", ylab=paste(expression(delta^18,"O (permil)")),
     xlim = rev(range(depth[interval]))
     )
lines(results_event$linramp$data$x,results_event$linramp$linrampfit$mean,lwd=1.5)
lines(results_event$linramp$data$x,results_event$linramp$linrampfit$q0.025,lwd=1.5,col="red")
lines(results_event$linramp$data$x,results_event$linramp$linrampfit$q0.975,lwd=1.5,col="red")
abline(v=depth.reference,lty=3,lwd=0.7,col="black")
ybottom = min(results_event$linramp$data$y)
ytop = min(results_event$linramp$linrampfit$q0.025)
margt0=results_event$linramp$param$t0$marg.t0
normt0.y = margt0[,2]/diff(range(margt0[,2]))*(ytop-ybottom)-min(margt0[,2])+ybottom
lines(x=margt0[,1],y=normt0.y,col="blue",lwd=2)
ybottom = min(results_event$linramp$data$y)
ytop = min(results_event$linramp$linrampfit$q0.025)
margt1 = results_event$linramp$param$t1$marginal
normt1.y = margt1[,2]/diff(range(margt1[,2]))*(ytop-ybottom)-min(margt1[,2])+ybottom
lines(x=margt1[,1],y=normt1.y,col="blue",lwd=2,lty=3)
plot(results_event$linramp$param$t0$marg.t0,type="l",xlab="Depth (m)",ylab="Density",
     main="(b) Marginal posterior of onset depth", 
     xlim=rev(range(depth.reference,results_event$linramp$param$t0$marg.t0[,1])))
abline(v=depth.reference,lty=3)
abline(v=results_event$linramp$param$t0$mean,lwd=1.5)
abline(v=results_event$linramp$param$t0$q0.025,col="gray")
abline(v=results_event$linramp$param$t0$q0.975,col="gray")
hist(results_event$event_dating$samples,main="(c) Histogram of onset age",
     xlab="Onset age (yr b2k)",col="orange",xlim=rev(range(results_event$event_dating$samples)), breaks=40)
abline(v=age.reference,lty=3)
abline(v=results_event$event_dating$mean,lwd=2)
par(mar=c(5,3.5,4,2)+0.1)
}
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_individual-1.png" title="plot of chunk events_individual" alt="plot of chunk events_individual" width="100%" />


We want to compare the dating onsets with those given by Rasmussen et al. (2014), Capron () and Buizert (). We start by importing those values


```r
library(readxl)
rasmussen_depths = event_intervals$NGRIP_depth_m
rasmussen_ages = event_intervals$GICC_age.yb2k #GISevents$`Age (a b2k)`
# buizert_onsets = read_excel("inst/Examples/Paper_results/data/Buizert_onsets.xlsx",
#                             col_types=c("text","numeric","numeric","numeric","numeric"))
buizert_depths = buizert_onsets$Depth; buizert_ages = buizert_onsets$Age
#buizert_depths[c(2:6,8:9,12:13,15:17,19:21)]
buizert_ages = buizert_ages[c(2:6,8:9,12:13,15:17,19:21)]
# capron_onsets_NGRIP_d18O = read_excel("inst/Examples/Paper_results/data/Capron_onsets.xls",
#                            sheet="d18O",n_max=25)
capron_ages_NGRIP_d18O = capron_onsets_NGRIP_d18O$`t1 (50%)`
# capron_onsets_NGRIP_dust = read_excel("inst/Examples/Paper_results/data/Capron_onsets.xls",
#                                 sheet="Ca2+",n_max=24)
capron_ages_NGRIP_dust = capron_onsets_NGRIP_dust$`t1 (50%)`
capind = numeric(29); buiind = numeric(29)
capind = c(NA,2,3, 4,5,6,NA,NA,7,8,NA,9,10,11,NA,NA,NA,NA,NA,12,13,14,NA,NA,15,NA,NA,16,17)
buiind = c(NA,2,NA,3,4,6,NA,NA,8,9,NA,11,12,13,NA,NA,NA,NA,NA,15,16,17,NA,NA,19,NA,NA,20,21)
```

We then generate a model fit using log-Ca as proxy


```r
dust = NGRIP_5cm$logdust
library(zoo)
dust_filled = na.approx(dust) #Impute missing data using linear interpolation
logdust = log(dust_filled)
data2=data
data2$proxy = logdust
results_dust = bremla(formula,data2,reference.label="GICC05",
                nsims=nsims,
                events=events,
                control.fit=list(noise="ar1"),
                control.sim=list(synchronized=FALSE) 
                )
```

We then estimate the dating uncertainty for all 29 events


```r
agesimmatrix_dust = matrix(NA,29,nsims)
agesimmatrix_d18O = matrix(NA,29,nsims)
depthlist_dust = c()
depthlist_d18O = c()
do.plot.rampfit = TRUE
do.plot.agehist = FALSE
depthstats = as.data.frame(matrix(NA,7,29))
colnames(depthstats) = c("true", "mean", "CIL","CIU","dustmean","dustCIL","dustCIU")
agestats = as.data.frame(matrix(NA,7,29))
colnames(agestats) = c("true", "mean", "CIL","CIU","dustmean","dustCIL","dustCIU")
# below is some fiddling with optimization parameters to improve the linear ramp fit for 
# some of the more challenging fits
steplengths = rep(0.005,29) 
steplengths[c(5,10)]=0.0001 #sometimes changing the steplength in the INLA optimization procedure can improve convergence (default is 0.005)
#compute.dust = rep(TRUE,29); compute.dust[c(8,27)]=FALSE #
#dtparams1 = rep(1,29); dtparams2 = rep(0.2,29) #prior-parameters for transition length
#t0params1 = rep(1,29); t0params2 = rep(0.2,29) #prior-parameters for transition onset
#dtparams1[10] = 20; dtparams2[10] = 0.1 #
#t0params1[10] = 120; t0params2[10] = 20
for(eventnumber in 1:29){
  if(eventnumber == 10){ #one of the transitions struggled to find good initial value 
    opt.params = c(40,NA,NA,NA) #for onset depth
    priorparams=NULL
  }else{
    opt.params = NULL
    priorparams=NULL
  }
    
    
  lowerints = which.index(event_intervals$depth_int_lower.m, results$data$depth)
  upperints = which.index(event_intervals$depth_int_upper.m, results$data$depth)
  interval_d18O = (lowerints[eventnumber]):upperints[eventnumber]+1
  
  lowerints = which.index(event_intervals$depth_int_lower.m, results_dust$data$depth)
  upperints = which.index(event_intervals$depth_int_upper.m, results_dust$data$depth)
  interval_dust = lowerints[eventnumber]:upperints[eventnumber]+1
  
  
  
  depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
  age.reference = event_intervals$GICC_age.yb2k[eventnumber]
  
  control.linramp_d18O = list(proxy=proxy,interval=interval_d18O,interval.unit="index",
                         depth.ref=depth.reference, #imp.fit=imp.fits[eventnumber],
                         h=steplengths[eventnumber],
                         opt.params=opt.params,
                         priorparams=priorparams,
                         rescale.y.factor=1,
                         imp.fit=TRUE,
                         silent=2L,
                         label=paste0(eventnumber," (d18O):", event_intervals[eventnumber,2])
                         )
  
  #rescale logdust proxy to avoid having to retune priors again
  rescale.y.factor = -20
  control.linramp_dust = list(proxy=logdust,interval=interval_dust,interval.unit="index",
                         depth.ref=depth.reference, #imp.fit=imp.fits[eventnumber],
                         h=steplengths[eventnumber],
                         opt.params=opt.params,
                         priorparams=priorparams,
                         rescale.y.factor=rescale.y.factor,
                         imp.fit=TRUE,
                         silent=2L,
                         label=paste0(eventnumber," (logCa):", event_intervals[eventnumber,2])
                         )
  
  temp_d18O = linrampfitter(results,control.linramp_d18O)
  if(do.plot.rampfit){
    {
      ybottom = min(temp_d18O$linramp$data$y)
      ytop = min(temp_d18O$linramp$linrampfit$q0.025)
      epsilon = (ytop-ybottom)*0.3
      
      margt0=temp_d18O$linramp$param$t0$marg.t0
      normt0.y = (margt0[,2]-min(margt0[,2]))/diff(range(margt0[,2]))*(ytop-ybottom)+ybottom-epsilon
      margt1 = temp_d18O$linramp$param$t1$marginal
      normt1.y = margt1[,2]/diff(range(margt1[,2]))*(ytop-ybottom)-min(margt1[,2])+ybottom-epsilon
      
      yrange = range(temp_d18O$linramp$data$y,temp_d18O$linramp$linrampfit$q0.025,
                     temp_d18O$linramp$linrampfit$q0.975,
                     normt0.y,normt1.y)
      
      plot(depth[interval_d18O],proxy[interval_d18O],type="l",col="gray",
          xlim=rev(range(depth[interval_d18O])), xlab="Depth (m)",ylab=expression(paste(delta^18,"O (permil)")),
          main=paste0(eventnumber," (d18O): ", event_intervals[eventnumber,2]),
          ylim=yrange)
      lines(temp_d18O$linramp$data$x,temp_d18O$linramp$linrampfit$mean)
      lines(temp_d18O$linramp$data$x,temp_d18O$linramp$linrampfit$q0.025,col="red")
      lines(temp_d18O$linramp$data$x,temp_d18O$linramp$linrampfit$q0.975,col="red")
      abline(v=rasmussen_depths[eventnumber],lty=3)
      
      lines(x=margt0[,1],y=normt0.y,col="blue",lwd=2)
      lines(x=margt1[,1],y=normt1.y,col="blue",lwd=2,lty=3)
      # ybottom = min(temp_d18O$linramp$data$y)
      # ytop = min(temp_d18O$linramp$linrampfit$q0.025)
      # margt0=temp_d18O$linramp$param$t0$marg.t0
      # normt0.y = margt0[,2]/diff(range(margt0[,2]))*(ytop-ybottom)-min(margt0[,2])+ybottom
      # lines(x=margt0[,1],y=normt0.y,col="blue",lwd=2)
      # ybottom = min(temp_d18O$linramp$data$y)
      # ytop = min(temp_d18O$linramp$linrampfit$q0.025)
      # margt1 = temp_d18O$linramp$param$t1$marginal
      # normt1.y = margt1[,2]/diff(range(margt1[,2]))*(ytop-ybottom)-min(margt1[,2])+ybottom
      # lines(x=margt1[,1],y=normt1.y,col="blue",lwd=2,lty=3)
    }
  }
  temp_dust = linrampfitter(results_dust,control.linramp_dust)
  if(do.plot.rampfit){
    {
      ybottom = min(temp_dust$linramp$data$y)
      ytop = min(temp_dust$linramp$linrampfit$q0.025)
      epsilon = (ytop-ybottom)*0.3
      
      margt0=temp_dust$linramp$param$t0$marg.t0
      normt0.y = (margt0[,2]-min(margt0[,2]))/diff(range(margt0[,2]))*(ytop-ybottom)+ybottom-epsilon
      margt1 = temp_dust$linramp$param$t1$marginal
      normt1.y = margt1[,2]/diff(range(margt1[,2]))*(ytop-ybottom)-min(margt1[,2])+ybottom-epsilon
      
      yrange = range(temp_dust$linramp$data$y,temp_dust$linramp$linrampfit$q0.025,
                     temp_dust$linramp$linrampfit$q0.975,
                     normt0.y,normt1.y)
      
      plot(temp_dust$linramp$data$x,temp_dust$linramp$data$y,type="l",col="gray",
          xlim=rev(range(depth[interval_dust])), xlab="Depth (m)",ylab=expression(paste("log(Ca)","")),
          main=paste0(eventnumber," (logdust): ", event_intervals[eventnumber,2]),
          ylim=yrange)
      lines(temp_d18O$linramp$data$x,temp_dust$linramp$linrampfit$mean)
      lines(temp_d18O$linramp$data$x,temp_dust$linramp$linrampfit$q0.025,col="red")
      lines(temp_d18O$linramp$data$x,temp_dust$linramp$linrampfit$q0.975,col="red")
      abline(v=rasmussen_depths[eventnumber],lty=3)
      
      lines(x=margt0[,1],y=normt0.y,col="blue",lwd=2)
      lines(x=margt1[,1],y=normt1.y,col="blue",lwd=2,lty=3)
       
      # plot(temp_dust$linramp$data$x,temp_dust$linramp$data$y,type="l",col="gray",xlim=rev(range(depth[interval_dust])), xlab="Depth (m)",ylab=expression(paste("log(Ca)","")),
      #     main=paste0(eventnumber," (logdust): ", event_intervals[eventnumber,2]))
      # # plot(depth[interval_d18O],logdust[interval_d18O],type="l",col="gray",xlim=rev(range(depth[interval_dust])), xlab="Depth (m)",ylab=expression(paste("log(Ca)","")),
      # #     main=paste0(eventnumber," (logdust): ", event_intervals[eventnumber,2]))
      # lines(temp_dust$linramp$data$x,temp_dust$linramp$linrampfit$mean)
      # lines(temp_dust$linramp$data$x,temp_dust$linramp$linrampfit$q0.025,col="red")
      # lines(temp_dust$linramp$data$x,temp_dust$linramp$linrampfit$q0.975,col="red")
      # abline(v=rasmussen_depths[eventnumber],lty=3)
      # 
      # ybottom = min(temp_dust$linramp$data$y)
      # ytop = min(temp_dust$linramp$linrampfit$q0.025)
      # margt0=temp_dust$linramp$param$t0$marg.t0
      # normt0.y = margt0[,2]/diff(range(margt0[,2]))*(ytop-ybottom)-min(margt0[,2])+ybottom
      # lines(x=margt0[,1],y=normt0.y,col="blue",lwd=2)
      # ybottom = min(temp_dust$linramp$data$y)
      # ytop = min(temp_dust$linramp$linrampfit$q0.025)
      # margt1 = temp_dust$linramp$param$t1$marginal
      # normt1.y = margt1[,2]/diff(range(margt1[,2]))*(ytop-ybottom)-min(margt1[,2])+ybottom
      # lines(x=margt1[,1],y=normt1.y,col="blue",lwd=2,lty=3)
    }
  }
  
  control.transition_dating=list(label="d18O",sync=FALSE)
  agetemp_d18O = events_depth_to_age(temp_d18O,control.transition_dating)
  
  
  agesimmatrix_d18O[eventnumber,] = agetemp_d18O$event_dating$samples 
  if(do.plot.agehist){
    hist(agesimmatrix_d18O[eventnumber,],freq=0,breaks=50,col="orange",main=paste0("d18O: ",eventnumber),xlab="Age (yb2k)")
    abline(v=age.reference,col="blue")
  }
  control.transition_dating_dust=list(label="dust",sync=FALSE)
  agetemp_dust = events_depth_to_age(temp_d18O,control.transition_dating_dust)
  agesimmatrix_dust[eventnumber,] = agetemp_dust$event_dating$samples
   if(do.plot.agehist){
    hist(agesimmatrix_dust[eventnumber,],freq=0,breaks=50,col="orange",main=paste0("dust: ",eventnumber),xlab="Age (yb2k)")
     abline(v=age.reference,col="blue")
  } 
    
  
  
  #compute summary statistics
  depthstats[1,eventnumber] = event_intervals$NGRIP_depth_m[eventnumber]
  depthstats[2,eventnumber] = agetemp_d18O$linramp$param$t0$mean
  depthstats[3,eventnumber] = agetemp_d18O$linramp$param$t0$q0.025
  depthstats[4,eventnumber] = agetemp_d18O$linramp$param$t0$q0.975
  depthstats[5,eventnumber] = agetemp_dust$linramp$param$t0$mean
  depthstats[6,eventnumber] = agetemp_dust$linramp$param$t0$q0.025
  depthstats[7,eventnumber] = agetemp_dust$linramp$param$t0$q0.975
  
  densage1 = density(agesimmatrix_d18O[eventnumber,]); densage1 = data.frame(x=densage1$x,y=densage1$y)
  densage2 = density(agesimmatrix_dust[eventnumber,]); densage2 = data.frame(x=densage2$x,y=densage2$y)
  zage1 = inla.zmarginal(densage1,silent=TRUE)
  zage2 = inla.zmarginal(densage2,silent=TRUE)
  
  agestats[1,eventnumber] = event_intervals$GICC_age.yb2k[eventnumber]
  agestats[2,eventnumber] = mean(agesimmatrix_d18O[eventnumber])
  agestats[3,eventnumber] = zage1$quant0.025
  agestats[4,eventnumber] = zage1$quant0.975
  agestats[5,eventnumber] = mean(agesimmatrix_dust[eventnumber])
  agestats[6,eventnumber] = zage1$quant0.025
  agestats[7,eventnumber] = zage1$quant0.975
  
}
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-1.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-2.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-3.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-4.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-5.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-6.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-7.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-8.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-9.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-10.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-11.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-12.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-13.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-14.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-15.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-16.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-17.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-18.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-19.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-20.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-21.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-22.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-23.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-24.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-25.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-26.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-27.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-28.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-29.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-30.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-31.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-32.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-33.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-34.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-35.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-36.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-37.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-38.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-39.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-40.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-41.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-42.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-43.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-44.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-45.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-46.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-47.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-48.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-49.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-50.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-51.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-52.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-53.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-54.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values

#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-55.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-56.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-57.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" /><img src="man/figures/reproduce/myrvollnilsen2022/results-events_all-58.png" title="plot of chunk events_all" alt="plot of chunk events_all" width="100%" />

## Results

To plot all posteriors with the comparison onsets run the code below.


```r
library(stringr)
par(mfrow=c(5,6),mar=c(2,2,1.5,1))
for(i in 1:29){
  dens1 = density(agesimmatrix_dust[i,]); dens1 = data.frame(x=dens1$x,y=dens1$y)#/max(dens1$y))
  dens2 = density(agesimmatrix_d18O[i,]); dens2 = data.frame(x=dens2$x,y=dens2$y)#/max(dens2$y))
  xlim = range(dens1$x,dens2$x);ylim = range(dens1$y,dens2$y)
  #par(mar=c(2,2,1.5,1))
  
  {
    plot(dens1,type="l",col=1,xlab="Onset year (b2k)",ylab="Density",
       main=str_sub(event_intervals[i,2],str_locate(event_intervals[i,2],"GI-")[1]),
       xlim=xlim,ylim=ylim)
  #Axis(side=1,labels=FALSE,at=dens1$x,tck=0.1)
  
  lines(dens2,type="l",col="gray")
  
  
  if(FALSE){#95% credible intervals
    abline(v=mean(agesimmatrix_dust[i,]),col=1)
    zdens1 = inla.zmarginal(dens1,silent = TRUE)
    abline(v=c(zdens1$quant0.025,zdens1$quant0.975),col=1,lty=3,lwd=0.8)
    abline(v=mean(agesimmatrix_d18O[i,]),col=2)
    zdens2 = inla.zmarginal(dens2,silent = TRUE)
    abline(v=c(zdens2$quant0.025,zdens2$quant0.975),col=2,lty=3,lwd=0.8)
  }
  
  abline(v=event_intervals$GICC_age.yb2k[i],col="blue")  
  if(!is.na(buiind[i])){
    abline(v=buizert_onsets$Age[buiind[i]],col="green")
    #abline(v=buizert_ages[buiind[i]],col="green")
  }
  
  
  # abline(v=capron_ages_NGRIP_d18O,col="orange")
  # abline(v=capron_ages_NGRIP_dust,col="orange",lty=3)
  if(!is.na(capind[i])){
    abline(v=capron_ages_NGRIP_d18O[capind[i]],col="red")
    abline(v=capron_ages_NGRIP_dust[capind[i]],col="pink",lty=1)
  }
  }
  
}
par(mar=c(0,0,0,0));plot(-1,axes=FALSE,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
legend(0,1,legend=c("d18O onset posterior","Ca2+ onset posterior",
                    "Rasmussen onset","Buizert onset",
                    "Capron d18O onset","Capron Ca2+ onset"),
       col=c("black","gray", "blue", "green", "red", "pink"),
       lty=c(1,1,1,1,1,1), cex=0.65,bty="n")
```

<img src="man/figures/reproduce/myrvollnilsen2022/results-events_plot-1.png" title="plot of chunk events_plot" alt="plot of chunk events_plot" width="100%" />

```r
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
rownames=c()
for(i in 1:29){
  rownames=c(rownames,str_sub(event_intervals[i,2],str_locate(event_intervals[i,2],"GI-")[1]))
}
```

To format the summary statistics of the onset depth and age into a data.frame object, run the code below.


```r
colnames = c("d18O_mean","d18O_q0.025","d18O_q0.975",
             "Ca_mean","Ca_q0.025","Ca_q0.975")
depthresults = t(depthstats)[,2:7]
colnames(depthresults)=colnames
rownames(depthresults)=rownames
ageresults = t(agestats[2:7,])
colnames(ageresults)=colnames
rownames(ageresults)=rownames

cat("Onset depth estimates:\n",sep="")
#> Onset depth estimates:
print(depthresults)
#>          d18O_mean d18O_q0.025 d18O_q0.975  Ca_mean Ca_q0.025 Ca_q0.975
#> GI-1d     1575.001    1574.870    1575.127 1575.001  1574.870  1575.127
#> GI-1e     1604.586    1604.505    1604.685 1604.586  1604.505  1604.685
#> GI-2.2    1794.002    1793.681    1794.357 1794.002  1793.681  1794.357
#> GI-3      1869.138    1868.966    1869.256 1869.138  1868.966  1869.256
#> GI-4      1891.683    1891.642    1891.718 1891.683  1891.642  1891.718
#> GI-5.2    1952.077    1952.001    1952.159 1952.077  1952.001  1952.159
#> GI-6      1974.444    1974.385    1974.498 1974.444  1974.385  1974.498
#> GI-7b     1997.376    1997.221    1997.531 1997.376  1997.221  1997.531
#> GI-7c     2009.804    2009.757    2009.856 2009.804  2009.757  2009.856
#> GI-8c     2070.039    2069.878    2070.139 2070.039  2069.878  2070.139
#> GI-9      2099.702    2099.671    2099.734 2099.702  2099.671  2099.734
#> GI-10     2124.460    2124.225    2124.703 2124.460  2124.225  2124.703
#> GI-11     2159.250    2159.209    2159.292 2159.250  2159.209  2159.292
#> GI-12c    2222.368    2222.187    2222.605 2222.368  2222.187  2222.605
#> GI-13b    2254.181    2254.032    2254.340 2254.181  2254.032  2254.340
#> GI-13c    2257.545    2257.160    2258.056 2257.545  2257.160  2258.056
#> GI-14b    2296.077    2295.910    2296.251 2296.077  2295.910  2296.251
#> GI-14c    2340.033    2339.888    2340.163 2340.033  2339.888  2340.163
#> GI-14d    2341.594    2341.552    2341.631 2341.594  2341.552  2341.631
#> GI-14e    2345.709    2345.580    2345.848 2345.709  2345.580  2345.848
#> GI-15.1   2355.398    2355.381    2355.412 2355.398  2355.381  2355.412
#> GI-15.2   2366.668    2366.347    2367.016 2366.668  2366.347  2367.016
#> GI-16.1b  2397.410    2397.033    2397.755 2397.410  2397.033  2397.755
#> GI-16.1c  2398.648    2398.423    2398.825 2398.648  2398.423  2398.825
#> GI-16.2   2402.352    2402.313    2402.394 2402.352  2402.313  2402.394
#> GI-17.1a  2409.624    2409.316    2410.038 2409.624  2409.316  2410.038
#> GI-17.1b  2411.301    2411.070    2411.485 2411.301  2411.070  2411.485
#> GI-17.1c  2414.894    2414.820    2414.988 2414.894  2414.820  2414.988
#> GI-17.2   2420.758    2420.655    2420.868 2420.758  2420.655  2420.868
cat("Onset age estimates:\n",sep="")
#> Onset age estimates:
print(ageresults)
#>          d18O_mean d18O_q0.025 d18O_q0.975  Ca_mean Ca_q0.025 Ca_q0.975
#> GI-1d     14098.85    14011.74    14142.80 14099.82  14011.74  14142.80
#> GI-1e     14710.88    14614.05    14768.25 14712.23  14614.05  14768.25
#> GI-2.2    23457.67    23261.45    23520.22 23449.70  23261.45  23520.22
#> GI-3      27846.77    27644.94    27927.34 27840.37  27644.94  27927.34
#> GI-4      28977.41    28770.15    29056.55 28976.63  28770.15  29056.55
#> GI-5.2    32604.51    32378.03    32685.70 32601.09  32378.03  32685.70
#> GI-6      33795.14    33582.09    33892.41 33795.37  33582.09  33892.41
#> GI-7b     35064.21    34875.03    35196.82 35064.75  34875.03  35196.82
#> GI-7c     35518.39    35343.81    35670.63 35518.20  35343.81  35670.63
#> GI-8c     38164.42    38054.56    38399.59 38155.90  38054.56  38399.59
#> GI-9      40084.84    39991.20    40347.91 40085.96  39991.20  40347.91
#> GI-10     41430.25    41312.53    41677.01 41416.63  41312.53  41677.01
#> GI-11     43399.22    43290.43    43662.77 43398.54  43290.43  43662.77
#> GI-12c    46810.81    46673.01    47062.67 46813.32  46673.01  47062.67
#> GI-13b    49083.04    48937.41    49338.70 49086.30  48937.41  49338.70
#> GI-13c    49270.28    49125.46    49527.45 49266.40  49125.46  49527.45
#> GI-14b    51603.31    51463.86    51874.30 51611.75  51463.86  51874.30
#> GI-14c    53864.79    53730.82    54151.67 53869.04  53730.82  54151.67
#> GI-14d    53955.31    53817.85    54237.97 53955.95  53817.85  54237.97
#> GI-14e    54157.16    54025.24    54448.85 54162.21  54025.24  54448.85
#> GI-15.1   54937.55    54791.71    55216.62 54937.21  54791.71  55216.62
#> GI-15.2   55767.84    55612.68    56046.33 55766.80  55612.68  56046.33
#> GI-16.1b  57910.54    57737.96    58180.21 57893.58  57737.96  58180.21
#> GI-16.1c  57980.53    57811.99    58250.45 57985.32  57811.99  58250.45
#> GI-16.2   58215.08    58048.70    58486.15 58216.62  58048.70  58486.15
#> GI-17.1a  58688.22    58549.18    58990.06 58712.15  58549.18  58990.06
#> GI-17.1b  58818.28    58654.16    59093.35 58816.17  58654.16  59093.35
#> GI-17.1c  59006.00    58851.00    59289.24 59007.44  58851.00  59289.24
#> GI-17.2   59408.67    59248.37    59689.21 59389.18  59248.37  59689.21
```


## Attribution

This code is associated and written for the papers Myrvoll-Nilsen et al. (2022) and Myrvoll-Nilsen et al. (202x) mentioned above. Feel free to use the code, but please cite the accompanying papers.

## License

The code in this repository is made available under the terms of the GNU (version >=2) License. For details, see LICENSE.md file.

---
output: github_document
---

<!-- results.md is generated from results.Rmd. Please edit that file -->



# bremla

This file includes the data and code needed to reproduce the results for

Myrvoll-Nilsen, E., Riechers, K. & Boers, N. (202x). Tie-point paper (title TBD)

# Myrvoll-Nilsen et al. (202x)

The data used is the NGRIP/GICC05 data set 'NGRIP_d18O_and_dust_5cm', downloaded from [Centre for Ice and Climate](iceandclimate.nbi.ku.dk) at the Niels Bohr Institute of the University of Copenhagen, Denmark, and the table of stadial-interstadial events presented by [Rasmussen et al. (2014)](https://www.sciencedirect.com/science/article/pii/S0277379114003485). These are included in the package and can be imported as follows. 

```r
#library(bremla)
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
# plot(depth,proxy,xlab="Age (yb2k)",ylab=expression(paste(d^18,"O (permil)")),type="l",
     # xlim=rev(range(depth)))+abline(v=events$locations,col="gray",lwd=0.7)
formula = dy~-1+depth2+proxy
```


We also use tie-points described by probability distributions obtained from Adolphi et al. (2018) and Muscheler et al. (2020). These can be imported as follows.


```r
library(ggplot2)
x.ref=c(11050,12050,13050,22050,42050)
adolphipdfs = adolphiloader(tieshifts=c(11050,12050,13050,22050,42050),plotdens = FALSE,
                            x.ref=x.ref)
lower = c(-77.53, -30.92, -12.37, 336.19, 101.51)+x.ref
upper = c(-41.09, 11.94, 21.99, 714.03, 509.08)+x.ref
means = c(-59.95, -7.04, 3.34, 521.44, 275.24)+x.ref
medians = c(-60.10, -6.42, 2.64, 532.36, 243.57)+x.ref
gg1 = ggplot(data=adolphipdfs$tie1,aes(x=x,y=y))+geom_line()+theme_bw()+xlab("GICC05 age (yb2k)")+
  ggtitle("Tie-point 1") + ylab("Density") + 
  geom_vline(aes(xintercept=x.ref[1]),color="black")+
  geom_vline(aes(xintercept=means[1]),color="blue")+
  geom_vline(aes(xintercept=lower[1]),color="red")+
  geom_vline(aes(xintercept=upper[1]),color="red")
gg2 = ggplot(data=adolphipdfs$tie2,aes(x=x,y=y))+geom_line()+theme_bw()+xlab("GICC05 age (yb2k)")+
  ggtitle("Tie-point 2") + ylab("Density") + 
  geom_vline(aes(xintercept=x.ref[2]),color="black")+
  geom_vline(aes(xintercept=means[2]),color="blue")+
  geom_vline(aes(xintercept=lower[2]),color="red")+
  geom_vline(aes(xintercept=upper[2]),color="red")
gg3 = ggplot(data=adolphipdfs$tie3,aes(x=x,y=y))+geom_line()+theme_bw()+xlab("GICC05 age (yb2k)")+
  ggtitle("Tie-point 3") + ylab("Density") + 
  geom_vline(aes(xintercept=x.ref[3]),color="black")+
  geom_vline(aes(xintercept=means[3]),color="blue")+
  geom_vline(aes(xintercept=lower[3]),color="red")+
  geom_vline(aes(xintercept=upper[3]),color="red")
gg4 = ggplot(data=adolphipdfs$tie4,aes(x=x,y=y))+geom_line()+theme_bw()+xlab("GICC05 age (yb2k)")+
  ggtitle("Tie-point 4") + ylab("Density") + 
  geom_vline(aes(xintercept=x.ref[4]),color="black")+
  geom_vline(aes(xintercept=means[4]),color="blue")+
  geom_vline(aes(xintercept=lower[4]),color="red")+
  geom_vline(aes(xintercept=upper[4]),color="red")
gg5 = ggplot(data=adolphipdfs$tie5,aes(x=x,y=y))+geom_line()+theme_bw()+xlab("GICC05 age (yb2k)")+
  ggtitle("Tie-point 5") + ylab("Density") + 
  geom_vline(aes(xintercept=x.ref[5]),color="black")+
  geom_vline(aes(xintercept=means[5]),color="blue")+
  geom_vline(aes(xintercept=lower[5]),color="red")+
  geom_vline(aes(xintercept=upper[5]),color="red")

library(ggpubr)
ggarrange(ggarrange(gg1,gg2,gg3,nrow=1),
          ggarrange(gg4,gg5,nrow=1),
          nrow=2)
```

<img src="man/figures/reproduce/results/results-adolphi-1.png" title="plot of chunk adolphi" alt="plot of chunk adolphi" width="100%" />

```r
# ggsave("adolphipdfs-4800x2000.eps", device=cairo_ps,width=4800,height=2400,units="px",dpi=500,limitsize=FALSE)
  
```

The model can be fitted using the \texttt{bremla} function, which estimates the posterior marginal distributions of the hyperparameters and also performs simulations of chronologies.



```r
nsims=10000
results = bremla(formula,data,reference.label="GICC05",
                nsims=nsims,
                events=events,
                control.fit=list(noise="ar1"),
                control.sim=list(synchronized=FALSE) ,
                print.progress=TRUE
                )
{
  
  par(mfrow=c(1,2))
  plot(results$fitting$inla$hyperparameters$posteriors$sigma_epsilon,type="l",xlab=expression(paste(sigma[epsilon])),ylab="Density")
  
  abline(v=results$fitting$inla$hyperparameters$results$sigma_epsilon$mean)
  abline(v=c(results$fitting$inla$hyperparameters$results$sigma_epsilon$quant0.025,
             results$fitting$inla$hyperparameters$results$sigma_epsilon$quant0.975),
             col="red")
         
  plot(results$fitting$inla$hyperparameters$posteriors$phi,type="l",xlab=expression(phi),ylab="Density")
  abline(v=results$fitting$inla$hyperparameters$results$phi$mean)
  abline(v=c(results$fitting$inla$hyperparameters$results$phi$quant0.025,
             results$fitting$inla$hyperparameters$results$phi$quant0.975),
             col="red")
 par(mfrow=c(1,1)) 
}


ggp = ggplot() + theme_bw() + xlim(rev(range(results$data$age)))+
  xlab("GICC05 (yb2k)")+ylab(expression(paste(delta^18,"O (permil)")))
nevents=results$.args$events$nevents
agestarts = numeric(nevents)
for(i in 1:nevents){
  indexes = which(results$data[[paste0("psi0_",i)]]>0)
  age_i = results$data$age[indexes]
  agestarts[i]=age_i[1]
  proxy_i = results$data$proxy[indexes]
  if(i %% 2 ==0 ){
    color="red"
  }else{
    color="blue"
  }
  ggp = ggp + geom_line(data=data.frame(age_i=age_i,proxy_i=proxy_i),
                        mapping=aes(x=age_i,y=proxy_i),col=color)
}
agestarts[nevents]=last(results$data$age)
ggp = ggp + geom_vline(xintercept = agestarts,col="black",size=0.15)
print(ggp)
# ggsave("proxyplot-4800x2000.eps", device=cairo_ps,width=4800,height=2000,units="px",dpi=500,limitsize=FALSE)

ggd = data.frame(age=data$age,depth=data$depth,proxy=data$proxy)
ggdr = data.frame(depths=events_rasmussen$depth[events_rasmussen$depth>1492.45 & events_rasmussen$depth<2426],ages=events_rasmussen$age[events_rasmussen$age>11703.07 & events_rasmussen$age < 59944.50])
library(ggplot2)
ggp2 = ggplot(data=ggd) + theme_bw()+
  geom_line(aes(x=age,y=proxy),col="black",size=0.5)+
  geom_vline(data=ggdr,mapping=aes(xintercept=ages),col="gray",size=0.5)+
  xlab("GICC05 age (yb2k)") + ylab(expression(paste(delta^18,"O (permil)")))
  
# ggsave("proxyplot-bw-4800x2000.eps", device=cairo_ps,width=4800,height=2000,units="px",dpi=500,limitsize=FALSE)

# ggplot(data=ggd)+geom_line(aes(x=depth,y=age))

```

The fitted model to the layer increments and the simulated chronologies (as compared to the GICC05 time scale) can be plotted as follows.


```r
{
  library(ggplot2)
  
  dsims = colDiffs(rbind(matrix(rep(results$initial$age,ncol(results$simulation$age)),nrow=1),results$simulation$age))
  dmeans = rowMeans(dsims)
  dsd = rowSds(dsims)
ggd = data.frame(y = results$data$dy,depth=results$data$depth,
                   mean=dmeans,
                   lower = dmeans-1.96*dsd,
                   upper = dmeans+1.96*dsd)
  ggp = ggplot(data=ggd,aes(x=depth)) + theme_bw()+ #xlim(rev(range(results$data$depth)))+
    geom_line(aes(y=y),col="gray")+
    geom_line(aes(y=mean),col="blue",size=0.3)+
    geom_ribbon(aes(ymin=lower,ymax=upper),color="red",alpha=0,size=0.05)+
    xlab("Depth (m)")+ylab("Layers per 5cm")
 print(ggp) 
 
 #ggsave("fitdiff-4800x2000.eps", device=cairo_ps,width=4800,height=2000,units="px",dpi=500,limitsize=FALSE)
}
```

<img src="man/figures/reproduce/results/results-unconstrained-1.png" title="plot of chunk unconstrained" alt="plot of chunk unconstrained" width="100%" />

```r
{
  gicc05=results$original.chron$age[1+1:n]
  ggd = data.frame(depth=results$data$depth,
                   mean=results$simulation$summary$mean-gicc05,
                   lower=results$simulation$summary$lower-gicc05,
                   upper=results$simulation$summary$upper-gicc05,
                   zeros = numeric(n))
  ggp = ggplot(data=ggd,aes(x=depth)) + theme_bw()+
    #geom_line(aes(y=zeros),color="gray",linetype="dashed")+
    geom_line(aes(y=mean),col="blue")+
    geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",col="red",alpha=0.3)+
    xlab("Depth (m)") + ylab("Simulated time scale - GICC05 (years)")
  print(ggp)
  
  #ggsave("unsync-4800x2000.eps", device=cairo_ps,width=4800,height=2000,units="px",dpi=500,limitsize=FALSE)
}
```

<img src="man/figures/reproduce/results/results-unconstrained-2.png" title="plot of chunk unconstrained" alt="plot of chunk unconstrained" width="100%" />


## Fixed tie-points


In the paper we demonstrate our model first using the Adolphi tie-points fixed to be equal to the mode of the probability distribution. We omit the first tie-point since it falls outside of our considered interval. The constrained chronologies can be computed as follows.


```r
fixedsims = cbind(#rep(inla.mmarginal(adolphipdfs$tie1),nsims),
                  rep(inla.mmarginal(adolphipdfs$tie2),nsims),
                  rep(inla.mmarginal(adolphipdfs$tie3),nsims),
                  rep(inla.mmarginal(adolphipdfs$tie4),nsims),
                  rep(inla.mmarginal(adolphipdfs$tie5),nsims)
                  )
synchronization = list(samples = fixedsims,nsims=nsims,
                       locations=c(#1492.50,
                                   1502.80,1532.70,1767.25,2132.40),
                       locations_unit="depth")
results_fixed = bremla(formula,data,reference.label="GICC05",
                nsims=nsims,
                events=events,
                synchronization=synchronization,
                control.fit=list(noise="ar1"),
                control.sim=list(synchronized=TRUE) ,
                print.progress=TRUE
                )
#> Initiating data formatting... completed!
#> Performing least squares fit... completed!
#> Performing INLA fit... completed.
#> Computing remaining posteriors using Monte Carlo simulation...
#> INLA fit completed in 9.303405 seconds!
#> Simulating 3000 hyperparameters from INLA posterior... completed!
#> Sampling fixed coefficients... completed in 20.93335 seconds.
#> Simulating synchronized chronologies...
#> Synchronous age simulation 1000/3000. Elapsed time: 13.7753 seconds...
#> Synchronous age simulation 2000/3000. Elapsed time: 26.8408 seconds...
#> Synchronous age simulation 3000/3000. Elapsed time: 40.04465 seconds...
#> Computing posterior marginal mean and 95% credible intervals from chronology samples...
#>  completed in 1.394116seconds.
{
  gicc05=results_fixed$original.chron$age[1+1:n]
  ggd = data.frame(depth=results_fixed$data$depth,
                   mean=results_fixed$simulation$summary_sync$mean-gicc05,
                   lower=results_fixed$simulation$summary_sync$lower-gicc05,
                   upper=results_fixed$simulation$summary_sync$upper-gicc05,
                   zeros = numeric(n))
  ggp = ggplot(data=ggd,aes(x=depth)) + theme_bw()+
    #geom_line(aes(y=zeros),color="gray",linetype="dashed")+
    geom_line(aes(y=mean),col="blue")+
    geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",col="red",alpha=0.3)+
    xlab("Depth (m)") + ylab("Simulated time scale - GICC05 (years)")
  print(ggp)
  
  #ggsave("fixed-4800x2000.eps", device=cairo_ps,width=4800,height=2000,units="px",dpi=500,limitsize=FALSE)
}
```

<img src="man/figures/reproduce/results/results-fixed-1.png" title="plot of chunk fixed" alt="plot of chunk fixed" width="100%" />

Synchronizing our chronologies with the Adolphi tie-points where the uncertainties is incorporated cam be done using Monte Carlo simulations. All of this can be done using the bremla function.


```r
library(ggplot2)
library(INLA)
synchronization = list(method="adolphi",nsims=nsims,
                       locations=c(#1492.50,
                                   1502.80,1532.70,1767.25,2132.40),
                       locations_unit="depth")
results_sync = bremla(formula,data,reference.label="GICC05",
                nsims=nsims,
                events=events,
                synchronization=synchronization,
                control.fit=list(noise="ar1"),
                control.sim=list(synchronized=2,summary=list(CI.type="hpd")) ,
                print.progress=TRUE
                )
#> Initiating data formatting... completed!
#> Performing least squares fit... completed!
#> Performing INLA fit... completed.
#> Computing remaining posteriors using Monte Carlo simulation...
#> INLA fit completed in 9.504719 seconds!
#> Simulating 3000 hyperparameters from INLA posterior... completed!
#> Simulating mean vector from fitted coefficients... completed in 21.61501 seconds!
#> Simulating chronologies...
#> Age simulation 1000/3000. Elapsed time: 5.246463 seconds...
#> Age simulation 2000/3000. Elapsed time: 10.44417 seconds...
#> Age simulation 3000/3000. Elapsed time: 15.6153 seconds...
#> Completed in 15.61997 seconds!
#> Computing posterior marginal mean and 95% credible intervals from chronology samples...
#>  completed in 76.86353seconds.
#> Sampling fixed coefficients... completed in 22.75464 seconds.
#> Simulating synchronized chronologies...
#> Synchronous age simulation 1000/3000. Elapsed time: 13.81857 seconds...
#> Synchronous age simulation 2000/3000. Elapsed time: 29.88376 seconds...
#> Synchronous age simulation 3000/3000. Elapsed time: 46.25854 seconds...
#> Computing posterior marginal mean and 95% credible intervals from chronology samples...
#>  completed in 102.344seconds.
tiepointranges = as.data.frame(matrix(NA,nrow=4,ncol=3))
gicc05=results_sync$original.chron$age[1+1:n]
gicc05tie = gicc05[results_sync$tie_points$locations_indexes]
colnames(tiepointranges)=c("lower","upper","depth")
for(i in 1:4){
  tiepointranges[i,1:2] = inla.hpdmarginal(0.95,adolphipdfs[[paste0("tie",i+1)]])-gicc05tie[i]
}
tiepointranges[,3] = depth[results_sync$tie_points$locations_indexes]
{
  
  ggd = data.frame(depth=results_sync$data$depth,
                   mean=results_sync$simulation$summary_sync$mean-gicc05,
                   lower=results_sync$simulation$summary_sync$lower-gicc05,
                   upper=results_sync$simulation$summary_sync$upper-gicc05,
                   zeros = numeric(n))
  ggp = ggplot(data=ggd,aes(x=depth)) + theme_bw()+
    #geom_line(aes(y=zeros),color="gray",linetype="dashed")+
    geom_line(aes(y=mean),col="blue")+
    geom_ribbon(aes(ymin=lower,ymax=upper),fill="red",col="red",alpha=0.3)+
    xlab("Depth (m)") + ylab("Simulated time scale - GICC05 (years)") +
    geom_linerange(data=tiepointranges,aes(x=depth,ymin=lower,ymax=upper),color="magenta")
  print(ggp)
  
  #ggsave("sync-4800x2000.eps", device=cairo_ps,width=4800,height=2000,units="px",dpi=500,limitsize=FALSE)
}
```

<img src="man/figures/reproduce/results/results-sync-1.png" title="plot of chunk sync" alt="plot of chunk sync" width="100%" />


# Applications




## Abrupt warmings

The onset of abrupt warmings is dated in two stages that are combined to express the complete dating uncertainty, pertaining to the uncertainty originating from both the age-depth relationship as well as estimating the depth in the record which represent the onset of the transition. The latter is estimated using a linear ramp model as described in Erhardt et al. (2019). The linear ramp fit can be performed and plotted as follows.




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
results_event = linrampfitter(results_sync,control.linramp)
control.transition_dating = list(label="Simulated transition")
results_event = events_depth_to_age(results_event,
                                    control.transition_dating=control.transition_dating)
{
  ggd1 = data.frame(depth=depth[interval],proxy=proxy[interval],
                    x=results_event$linramp$data$x,
                    mean=results_event$linramp$linrampfit$mean,
                    lower=results_event$linramp$linrampfit$q0.025,
                    upper=results_event$linramp$linrampfit$q0.975
                    )
  gg1 = ggplot(data=ggd1,aes(x=depth))+geom_line(aes(y=proxy),col="gray")+
    theme_bw()+ xlim(rev(range(depth[interval])))+xlab("Depth (m)")+
    ylab(expression(paste(delta^18,"O (permil)")))+ggtitle("(a) GI-11 Linear ramp fit")+
    geom_ribbon(aes(x=x,ymin=lower,ymax=upper),fill="red",col="red",alpha=0.3)+
    geom_line(aes(x=x,y=mean),col="blue",size=0.7)
    
  ybottom = min(results_event$linramp$data$y)
  ytop = min(results_event$linramp$linrampfit$q0.025)
  margt0=results_event$linramp$param$t0$marg.t0
  normt0.y = margt0[,2]/diff(range(margt0[,2]))*(ytop-ybottom)-min(margt0[,2])+ybottom
  gg1 = gg1 + geom_line(data=data.frame(x=margt0[,1],y=normt0.y),aes(x=x,y=y),col="blue")
  
  margt1 = results_event$linramp$param$t1$marginal
  normt1.y = margt1[,2]/diff(range(margt1[,2]))*(ytop-ybottom)-min(margt1[,2])+ybottom
  gg1 = gg1 + geom_line(data=data.frame(x=margt1[,1],y=normt1.y),aes(x=x,y=y),col="blue",
                        linetype="dashed")
  
  gg1 = gg1 + geom_vline(aes(xintercept=depth.reference),linetype="dotted",size=0.8)
  
  gd2 = as.data.frame(results_event$linramp$param$t0$marg.t0)
  gg2 = ggplot(data=gd2,mapping=aes(x=x,y=y)) + geom_line()+theme_bw()+
    xlab("Onset depth (m)") + ylab("Density") + ggtitle("(b) Marginal posterior of onset depth")+
    geom_vline(aes(xintercept = results_event$linramp$param$t0$mean),col="blue")+
    geom_vline(aes(xintercept = results_event$linramp$param$t0$q0.025),col="red")+
    geom_vline(aes(xintercept = results_event$linramp$param$t0$q0.975),col="red")
  
  
  gd3 = as.data.frame(results_event$linramp$param$t1$marginal)
  gg3 = ggplot(data=gd3,mapping=aes(x=x,y=y)) + geom_line()+theme_bw()+
    xlab("Onset depth (m)") + ylab("Density") + ggtitle("(c) Marginal posterior of end-point depth")+
    geom_vline(aes(xintercept = results_event$linramp$param$t1$mean),col="blue")+
    geom_vline(aes(xintercept = results_event$linramp$param$t1$q0.025),col="red")+
    geom_vline(aes(xintercept = results_event$linramp$param$t1$q0.975),col="red")
  
  
  ggarrange(gg1,ggarrange(gg2,gg3,nrow=1),nrow=2)
  #ggsave("GI-11-rampfit-4800x3000.eps", device=cairo_ps,width=4800,height=3000,units="px",dpi=500,limitsize=FALSE)
  
  gghist = ggplot(data = data.frame(ages=results_event$event_dating$samples)) + theme_bw()+
    geom_histogram(aes(x=ages,y=..density..),fill="orange",col="black")+xlab("Onset age (yb2k)")+
    ylab("Density") + ggtitle("GI-11 Onset age") +
    geom_vline(aes(xintercept = age.reference),lty=3,size=0.8)+
    geom_vline(aes(xintercept=results_event$event_dating$mean),col="blue",size=0.7)+
    geom_vline(aes(xintercept=results_event$event_dating$q0.025),col="red",size=0.7)+
    geom_vline(aes(xintercept=results_event$event_dating$q0.975),col="red",size=0.7)+
    xlim(rev(range(results_event$event_dating$samples)))
  
  #ggsave("GI-11-histogram-4800x3000.eps", device=cairo_ps,width=4800,height=3000,units="px",dpi=500,limitsize=FALSE)
  
}
```


We want to compare the dating onsets with those given by Rasmussen et al. (2014), Capron (2021) and Buizert et al. (2015). We start by importing those values


```r
##### IMPORT MYRVOLL-NILSEN 2022 DATA AND COMPARE POSTERIORS ####
library(readxl)
rasmussen_depths = event_intervals$NGRIP_depth_m
rasmussen_ages = event_intervals$GICC_age.yb2k #GISevents$`Age (a b2k)`
buizert_onsets = read_excel("data/Buizert_onsets.xlsx",
                            col_types=c("text","numeric","numeric","numeric","numeric"))
buizert_depths = buizert_onsets$Depth; buizert_ages = buizert_onsets$Age
#buizert_depths[c(2:6,8:9,12:13,15:17,19:21)]
buizert_ages = buizert_ages[c(2:6,8:9,12:13,15:17,19:21)]
capron_onsets_NGRIP_d18O = read_excel("data/Capron_onsets.xls",
                           sheet="d18O",n_max=25)
capron_ages_NGRIP_d18O = capron_onsets_NGRIP_d18O$`t1 (50%)`
capind = numeric(29); buiind = numeric(29)
capind = c(NA,2,3, 4,5,6,NA,NA,7,8,NA,9,10,11,NA,NA,NA,NA,NA,12,13,14,NA,NA,15,NA,NA,16,17)
buiind = c(NA,2,NA,3,4,6,NA,NA,8,9,NA,11,12,13,NA,NA,NA,NA,NA,15,16,17,NA,NA,19,NA,NA,20,21)
```

We then estimate the dating uncertainty for all 29 events

```r
agesimmatrix_d18O = matrix(NA,29,nsims)
agesimmatrix_d18O_unsync = matrix(NA,29,nsims)
depthlist_d18O = c()
do.plot.rampfit = TRUE
do.plot.agehist = FALSE
depthstats = as.data.frame(matrix(NA,4,29))
colnames(depthstats) = c("true", "mean", "CIL","CIU","dustmean","dustCIL","dustCIU")
agestats = as.data.frame(matrix(NA,4,29))
colnames(agestats) = c("true", "mean", "CIL","CIU","dustmean","dustCIL","dustCIU")
steplengths = rep(0.01,29)
steplengths[21]=0.001 #sometimes changing the steplength slightly might improve convergence (default is 0.01)
for(eventnumber in 1:29){
  transitionlabel = str_sub(event_intervals$onsetlabel[eventnumber],
                      str_locate(event_intervals$onsetlabel[eventnumber],"GI")[1])
  lowerints = which.index(event_intervals$depth_int_lower.m, results_sync$data$depth)
  upperints = which.index(event_intervals$depth_int_upper.m, results_sync$data$depth)
  interval = lowerints[eventnumber]:upperints[eventnumber]
  
  depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
  age.reference = event_intervals$GICC_age.yb2k[eventnumber]
  
  control.linramp = list(proxy=proxy,interval=interval,interval.unit="index",
                         depth.ref=depth.reference, 
                         label=paste0(eventnumber," (d18O):", event_intervals[eventnumber,2])
                         )
  
  temp_d18O = linrampfitter(results_sync,control.linramp)
  if(do.plot.rampfit){
    {plot(depth[interval],proxy[interval],type="l",col="gray",
          xlim=rev(range(depth[interval])), xlab="Depth (m)",ylab=expression(paste(delta^18,"O (permil)")),
          main=paste0(eventnumber,": ", transitionlabel))
      lines(temp_d18O$linramp$data$x,temp_d18O$linramp$linrampfit$mean)
      lines(temp_d18O$linramp$data$x,temp_d18O$linramp$linrampfit$q0.025,col="red")
      lines(temp_d18O$linramp$data$x,temp_d18O$linramp$linrampfit$q0.975,col="red")
      abline(v=rasmussen_depths[eventnumber],lty=3)
      
      ybottom = min(temp_d18O$linramp$data$y)
      ytop = min(temp_d18O$linramp$linrampfit$q0.025)
      margt0=temp_d18O$linramp$param$t0$marg.t0
      normt0.y = margt0[,2]/diff(range(margt0[,2]))*(ytop-ybottom)-min(margt0[,2])+ybottom
      lines(x=margt0[,1],y=normt0.y,col="blue",lwd=2)
      ybottom = min(temp_d18O$linramp$data$y)
      ytop = min(temp_d18O$linramp$linrampfit$q0.025)
      margt1 = temp_d18O$linramp$param$t1$marginal
      normt1.y = margt1[,2]/diff(range(margt1[,2]))*(ytop-ybottom)-min(margt1[,2])+ybottom
      lines(x=margt1[,1],y=normt1.y,col="blue",lwd=2,lty=3)
    }
  }
  
  
  control.transition_dating=list(label="d18O")
  agetemp_d18O = events_depth_to_age(temp_d18O,control.transition_dating)
  agetemp_d18O_unsync = events_depth_to_age(temp_d18O,
                                            control.transition_dating=list(label="d18O",
                                                                           sync=FALSE))
  
  
  agesimmatrix_d18O[eventnumber,] = agetemp_d18O$event_dating$samples 
  agesimmatrix_d18O_unsync[eventnumber,] = agetemp_d18O_unsync$event_dating$samples 
  if(do.plot.agehist){
    hist(agesimmatrix_d18O[eventnumber,],freq=0,breaks=50,col="orange",main=paste0("Onset age: ",transitionlabel),xlab="Age (yb2k)")
    abline(v=age.reference,col="blue")
  }
  
  #Store summary statistics
  depthstats[1,eventnumber] = event_intervals$NGRIP_depth_m[eventnumber]
  depthstats[2,eventnumber] = agetemp_d18O$linramp$param$t0$mean
  depthstats[3,eventnumber] = agetemp_d18O$linramp$param$t0$q0.025
  depthstats[4,eventnumber] = agetemp_d18O$linramp$param$t0$q0.975
  
  densage = density(agesimmatrix_d18O[eventnumber,]); densage = data.frame(x=densage$x,y=densage$y)
  
  zage = inla.zmarginal(densage,silent=TRUE)
  
  agestats[1,eventnumber] = event_intervals$GICC_age.yb2k[eventnumber]
  agestats[2,eventnumber] = mean(agesimmatrix_d18O[eventnumber])
  agestats[3,eventnumber] = zage$quant0.025
  agestats[4,eventnumber] = zage$quant0.975
  
}
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-1.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-2.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-3.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-4.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-5.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-6.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-7.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-8.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />

```
#> Warning in log(fit$par[2]): NaNs produced

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values

#> Warning in log(fit$par[2]): collapsing to unique 'x' values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-9.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-10.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-11.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-12.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-13.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-14.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-15.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-16.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-17.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-18.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-19.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-20.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-21.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-22.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-23.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-24.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-25.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />

```
#> Warning in regularize.values(x, y, ties, missing(ties)): collapsing to unique 'x'
#> values
```

<img src="man/figures/reproduce/results/results-compare_transitions_1-26.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-27.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-28.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" /><img src="man/figures/reproduce/results/results-compare_transitions_1-29.png" title="plot of chunk compare_transitions_1" alt="plot of chunk compare_transitions_1" width="100%" />



```r

  # library(cowplot)
  
  ggplots = c()

  for(i in 1:29){
    dens = density(agesimmatrix_d18O[i,])
    dens = data.frame(x=dens$x,y=dens$y)#/max(dens2$y))
    densunsync = density(agesimmatrix_d18O_unsync[i,])
    densunsync = data.frame(x=densunsync$x,y=densunsync$y)#/max(dens2$y))
    xlim = rev(range(dens$x,densunsync$x))
    
    ggd = data.frame(x=dens$x,y=dens$y,x0=densunsync$x,y0=densunsync$y)
    ggdrasmus = data.frame(rasmus=event_intervals$GICC_age.yb2k[i])
    title = str_sub(event_intervals[i,2],str_locate(event_intervals[i,2],"GI-")[1])
    
    gg = ggplot(data=ggd) + theme_bw()+xlim(xlim)+
      geom_line(aes(x=x0,y=y0),col="gray",size=0.7)+
      geom_line(aes(x=x,y=y),col="black",size=0.7)+
      xlab("Onset year (yb2k)")+ylab("Density") +
      ggtitle(title)+
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank())+
      # geom_vline(aes(xintercept = event_intervals$GICC_age.yb2k[i]),col="blue",size=0.7)
      geom_vline(data=ggdrasmus, aes(xintercept = rasmus),col="blue",size=0.7)
    
      
      

  if(!is.na(buiind[i])){
    ggdbuizert = data.frame(buizert=buizert_onsets$Age[buiind[i]])
    gg = gg + geom_vline(data=ggdbuizert,aes(xintercept = buizert),col="green",size=0.7)
    #abline(v=buizert_onsets$Age[buiind[i]],col="green")
    #abline(v=buizert_ages[buiind[i]],col="green")
  }

  if(!is.na(capind[i])){
    ggdcapron = data.frame(capron=capron_ages_NGRIP_d18O[capind[i]])
    gg = gg + geom_vline(data=ggdcapron,aes(xintercept = capron),col="red",size=0.7)
    # abline(v=capron_ages_NGRIP_d18O[capind[i]],col="red")
  }
    # ggplots = c(ggplots,list(gg))
    ggplots[[i]] = gg
  }
  
  ggl=ggplot(data.frame(x = 1,x0=1, y = 1,y0=1, rasmus=1,buizert=1,capron=1))+
    geom_line(aes(x=x,y=y,col="Synchronized onset posterior"),size=1)+
    geom_line(aes(x=x0,y=y0,col="Unsynchronized onset posterior"))+
    geom_line(aes(x=x,y=rasmus,col="Rasmussen onset"))+
    geom_line(aes(x=x,y=buizert,col="Buizert onset"))+
    geom_line(aes(x=x,y=capron,col="Capron onset"))+
                        # colour = 'Something2'), aes(x, y, fill = colour))+
  # geom_point(alpha=0, shape = 0)+ # completely transparent rectangular point 
  # scale_fill_manual(values='black', drop=FALSE) +
   # guides(color = guide_legend(override.aes = list(alpha=1, size = 1.3)))+ # showing the point in the legend
    scale_color_manual(values=c("Synchronized onset posterior"="black",
                                "Unsynchronized onset posterior"="gray",
                                "Rasmussen et al. (2014)"="blue",
                                "Buizert et al. (2015)"="green",
                                "Capron et al. (2021)"="red"),
                       labels=c("Synchronized onset posterior",
                                "Unsynchronized onset posterior",
                                "Rasmussen et al. (2014)",
                                "Buizert et al. (2015)",
                                "Capron et al. (2021)"))+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.5, 0.5), # move the legend to the center
        legend.title = element_blank(),
        # legend.text = element_text(size = 17),
        legend.key = element_rect(fill='NA'),
        # legend.key.size = unit(0.8,'cm'),
        panel.grid = element_blank(),
        # legend.spacing = unit(3,'cm'),
        # legend.spacing.y = unit(0.5,'cm'),
        panel.border = element_rect(colour = "white", fill='white', size=1)
  )
  ggplots[[30]] = ggl
  ggarrange(plotlist=ggplots,widths=c(6,5))
#> geom_path: Each group consists of only one observation. Do you need to adjust
#> the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to adjust
#> the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to adjust
#> the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to adjust
#> the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to adjust
#> the group aesthetic?
```

<img src="man/figures/reproduce/results/results-compare_transitions_2-1.png" title="plot of chunk compare_transitions_2" alt="plot of chunk compare_transitions_2" width="100%" />

```r
  
  # ggsave("allevents-4800x3000.eps",device=cairo_ps,width=4800,height=3000,units="px",limitsize=FALSE)
    

rownames=c()
for(i in 1:29){
  rownames=c(rownames,str_sub(event_intervals[i,2],str_locate(event_intervals[i,2],"GI-")[1]))
}
#colnames = colnames(depthstats)[1:4]
#colnames = c("name","d18O_mean","d18O_q0.025","d18O_q0.975")
depthresults = t(depthstats)[,2:4]
colnames(depthresults)=c("Z* mean", "Z* q0.025","Z* q0.975")
ageresults = t(agestats[2:4,])
colnames(ageresults)=c("Y* mean", "Y* q0.025","Y* q0.975")
allresults = cbind(depthresults,ageresults)
rownames(allresults)=rownames
latexstr = "Event & Z^* mean (m) & Z^* 95\\% CI (m) & Y^* mean (yb2k) & Y^* 95\\% CI (yb2k) \\\\ \n"
for(i in 1:29){
  latexstr = paste0(latexstr,rownames[i], "& ", format(round(allresults[i,1],2),nsmall=2)," & (", format(round(allresults[i,2],2),nsmall=2), ",
         ",format(round(allresults[i,3],2),nsmall=2),") & ", format(round(allresults[i,4],2),nsmall=2)," & (",format(round(allresults[i,5],2),nsmall=2),", ",format(round(allresults[i,6],2),nsmall=2),") \\\\ \n"
         )
}
cat(latexstr)
#> Event & Z^* mean (m) & Z^* 95\% CI (m) & Y^* mean (yb2k) & Y^* 95\% CI (yb2k) \\ 
#> GI-1d& 1575.04 & (1574.87,
#>          1575.29) & 14166.01 & (14115.33, 14231.56) \\ 
#> GI-1e& 1604.57 & (1604.50,
#>          1604.65) & 14863.67 & (14767.35, 14930.62) \\ 
#> GI-2.2& 1793.99 & (1793.87,
#>          1794.10) & 23887.78 & (23702.23, 24068.44) \\ 
#> GI-3& 1869.23 & (1869.00,
#>          1869.46) & 28217.23 & (28065.24, 28409.62) \\ 
#> GI-4& 1891.68 & (1891.37,
#>          1891.99) & 29346.89 & (29171.93, 29511.88) \\ 
#> GI-5.2& 1952.11 & (1951.99,
#>          1952.23) & 32884.78 & (32757.73, 33089.86) \\ 
#> GI-6& 1974.44 & (1974.38,
#>          1974.50) & 34084.88 & (33948.41, 34282.15) \\ 
#> GI-7b& 1997.67 & (1997.38,
#>          1997.99) & 35375.32 & (35241.80, 35583.23) \\ 
#> GI-7c& 2009.80 & (2009.75,
#>          2009.85) & 35820.55 & (35692.60, 36038.07) \\ 
#> GI-8c& 2069.84 & (2069.78,
#>          2069.93) & 38498.11 & (38355.32, 38723.83) \\ 
#> GI-9& 2099.70 & (2099.67,
#>          2099.75) & 40441.97 & (40284.32, 40670.93) \\ 
#> GI-10& 2124.46 & (2124.23,
#>          2124.69) & 41739.72 & (41588.32, 41991.35) \\ 
#> GI-11& 2159.25 & (2158.97,
#>          2159.52) & 43698.10 & (43553.58, 43974.95) \\ 
#> GI-12c& 2222.78 & (2222.50,
#>          2223.05) & 47104.87 & (46963.53, 47412.90) \\ 
#> GI-13b& 2254.16 & (2254.09,
#>          2254.24) & 49390.05 & (49196.67, 49648.56) \\ 
#> GI-13c& 2257.60 & (2257.26,
#>          2257.93) & 49585.62 & (49387.68, 49845.44) \\ 
#> GI-14b& 2296.21 & (2295.89,
#>          2296.64) & 51920.06 & (51730.32, 52199.74) \\ 
#> GI-14c& 2339.86 & (2339.80,
#>          2339.91) & 54227.68 & (53981.43, 54462.51) \\ 
#> GI-14d& 2341.52 & (2341.41,
#>          2341.67) & 54312.57 & (54073.11, 54554.72) \\ 
#> GI-14e& 2345.72 & (2345.58,
#>          2345.89) & 54526.32 & (54286.87, 54769.57) \\ 
#> GI-15.1& 2355.47 & (2355.37,
#>          2355.59) & 55298.39 & (55058.56, 55543.84) \\ 
#> GI-15.2& 2366.61 & (2366.36,
#>          2366.89) & 56136.49 & (55868.87, 56359.96) \\ 
#> GI-16.1b& 2397.40 & (2397.27,
#>          2397.53) & 58233.18 & (57998.84, 58497.48) \\ 
#> GI-16.1c& 2398.74 & (2398.59,
#>          2398.89) & 58304.78 & (58076.05, 58574.68) \\ 
#> GI-16.2& 2402.35 & (2402.32,
#>          2402.38) & 58539.94 & (58307.11, 58806.38) \\ 
#> GI-17.1a& 2409.50 & (2409.24,
#>          2409.77) & 59012.75 & (58799.03, 59301.86) \\ 
#> GI-17.1b& 2411.94 & (2411.75,
#>          2412.18) & 59164.78 & (58947.70, 59450.43) \\ 
#> GI-17.1c& 2414.98 & (2414.93,
#>          2415.03) & 59346.93 & (59113.55, 59617.88) \\ 
#> GI-17.2& 2420.76 & (2420.65,
#>          2420.87) & 59746.16 & (59506.63, 60010.72) \\
```



# Other stuff




To get started, see the documentation '?bremla' and the associated example.

The package includes the NGRIP/GICC05 data set 'NGRIP_d18O_and_dust_5cm' downloaded from [Centre for Ice and Climate](iceandclimate.nbi.ku.dk) at the Niels Bohr Institute of the University of Copenhagen, Denmark, and the table of stadial-interstadial events presented by [Rasmussen et al. (2014)](https://www.sciencedirect.com/science/article/pii/S0277379114003485). Also includes the tie-point presented in [Adolphi et al. (2018)](https://cp.copernicus.org/articles/14/1755/2018/) and [Muscheler et al. (2020)](https://www.cambridge.org/core/journals/radiocarbon/article/testing-and-improving-the-intcal20-calibration-curve-with-independent-records/D72D9214C47FE9441B5E730D33DCCE3D). Other data and tie-points should work as well.

Warning: The provided real data example generates several thousand simulated chronologies with individual lengths exceeding 18,000. More if synchronized chronologies is also computed. This will require a significant amount of memory. If needed, reduce the number of samples generated or free up memory before use.

## Illustration example

The data used to generate the illustrative examples in the paper is as follows.


```r
rm(list=ls())
require(stats)
library(ggplot2)
set.seed(123)
n=800
y0 = 0
z0 = 0
phi=0.8; sigma=1.2
depth = z0+seq(from=1,to=n,length.out=n)

covariate = arima.sim(model=list(ar=c(0.99)),n=n)
noise = sigma*arima.sim(model=list(ar=c(phi)),n=n,sd=sqrt(1-phi^2))

dy = 10 + depth*5/n+covariate+noise
y = y0+cumsum(dy)

data=data.frame(age=c(y0,y),depth=c(z0,depth),dy=c(0,as.numeric(dy)),covariate=c(0,covariate))
#data=rbind(c(y0,z0,0,0),df)
formula=dy~ 1 + depth + covariate 
nsims = 10000

# Unsynchronized
res = bremla(formula,data,reference.label="simulated time scale",nsims=nsims,
             control.fit=list(noise="ar1"),
             x.label="Depth",
             y.label="Ensemble age - observed age",
             control.sim = list(synchronized=FALSE,
                                nsims=nsims))

library(ggplot2)
orichron = res$data$age
ggd = data.frame(depth=res$data$depth,
                 mean = res$simulation$summary$mean-orichron,
                 lower = res$simulation$summary$lower-orichron,
                 upper = res$simulation$summary$upper-orichron)
ggp = ggplot(data=ggd,aes(depth))+theme_bw()+xlab("Depth")+ylab("Ensemble age - observed age")+
  geom_ribbon(aes(ymin=lower,ymax=upper),col="red",fill="red",alpha=0.3,size=0.7)+
  geom_line(aes(y=mean),col="blue",size=0.7)
# plot(res)
#ggsave("illustrative-unsync.eps",device=cairo_ps,width=4800,height=2000,units="px")

# Fixed tie-points
ntie = 3
tiesims = matrix(NA,nrow=nsims,ncol=ntie)
tiedepths = c(100,300,600)
tieage = c(450,3000,8500)

for(i in 1:ntie){
  tiesims[,i] = rep(tieage[i],nsims)
}

res2 = tiepointsimmer(res,synchronization = list(locations=tiedepths,locations_unit="depth",
                       samples = tiesims)
               )

res2 = bremla_synchronized_simulation(res2,control.sim=list(synchronized=TRUE))

ggd = data.frame(depth=res2$data$depth,
                 mean = res2$simulation$summary_sync$mean-orichron,
                 lower = res2$simulation$summary_sync$lower-orichron,
                 upper = res2$simulation$summary_sync$upper-orichron)

oritie = orichron[tiedepths]
ggd2 = data.frame(depth=res2$tie_points$locations,
                  mean=res2$simulation$summary_sync$mean[tiedepths]-oritie)

ggp = ggplot(data=ggd,aes(depth))+theme_bw()+xlab("Depth")+ylab("Ensemble age - observed age")+
  geom_ribbon(aes(ymin=lower,ymax=upper),col="red",fill="red",alpha=0.3,size=0.7)+
  geom_line(aes(y=mean),col="blue",size=0.7)+
  geom_point(data=ggd2,mapping=aes(x=depth,y=mean),col="magenta",size=1)
#plot(res)
#ggsave("illustrative-fixedsync.eps",device=cairo_ps,width=4800,height=2000,units="px")

tiemeans = c(450,3000,8500)
tiesds = c(10,20,30)
res3 = tiepointsimmer(res,synchronization = list(locations=tiedepths,
                                                locations_unit="depth",
                                                method="normal",
                                                params=list(mean=tiemeans,
                                                            sd=tiesds))
               )
res3 = bremla_synchronized_simulation(res3,control.sim=list(synchronized=TRUE))
ggd = data.frame(depth=res3$data$depth,
                 mean = res3$simulation$summary_sync$mean-orichron,
                 lower = res3$simulation$summary_sync$lower-orichron,
                 upper = res3$simulation$summary_sync$upper-orichron)

oritie = orichron[tiedepths]
ggd2 = data.frame(depth=res3$tie_points$locations,
                  lower=res3$simulation$summary_sync$lower[tiedepths]-oritie,
                  upper=res3$simulation$summary_sync$upper[tiedepths]-oritie)

ggp = ggplot(data=ggd,aes(depth))+theme_bw()+xlab("Depth")+ylab("Ensemble age - observed age")+
  geom_ribbon(aes(ymin=lower,ymax=upper),col="red",fill="red",alpha=0.3,size=0.7)+
  geom_line(aes(y=mean),col="blue",size=0.7)+
  geom_segment(data=ggd2,mapping=aes(x=depth,xend=depth,y=lower,yend=upper), 
               col="magenta",size=0.7)
  
#plot(res3)
#ggsave("illustrative-sync.eps",device=cairo_ps,width=4800,height=2000,units="px")
```

## Simulation example

To verify our model we can test it on simulated data where the true parameters are known.


```r
rm(list=ls())
require(stats)
library(ggplot2)
set.seed(123)
n=1000
y0 = 1000
z0 = 1200
phi=0.8; sigma=1.2
depth = z0+seq(from=0.5,500,length.out=1000)
depth2 = depth^2/z0 #normalize for stability

dnoise = sigma*arima.sim(n=n,model=list(ar=c(phi)),sd=sqrt(1-phi^2))
proxy = arima.sim(n=n,model=list(ar=0.9),sd=sqrt(1-0.9^2))
events=list(locations=c(1,250,750,1000))
a1 = c(-rep(1,249)*15,numeric(751))
a2 = c(numeric(249),rep(1,500)*50,numeric(251))
a3 = c(numeric(749),-rep(1,251)*110)
c1 = c(depth[1:249]/60,numeric(751))
c2 = c(numeric(249),-depth[1:500]/40,numeric(251))
c3 = c(numeric(749),depth[1:251]/10)
dy = dnoise + a1+a2+a3+c1+c2+c3 + proxy*0.5
age = y0+cumsum(dy)
df=data.frame(age=age,dy=as.numeric(dy),proxy=as.numeric(proxy),depth=depth,depth2=depth2); data=rbind(c(y0,z0,0,0,0),df)
formula=dy~-1+proxy#+depth2
nsims = 10000
synchronization = list(locations = c(200,400,800),locations_unit="index",method="gauss", nsims=nsims,params=list(mean=c(2200,5600,11500),sd=c(10,50,60)))
res = bremla(formula,data,reference.label="simulated time scale",nsims=nsims,
             events=list(locations=c(1,250,750,1000),locations_unit="index"),
             control.fit=list(noise="ar1"),
             synchronization=synchronization,
             control.sim = list(synchronized=TRUE,
                                nsims=nsims))
par(mfrow=c(1,2))
plot(res$fitting$inla$hyperparameters$posteriors$sigma_epsilon,type="l",xlab=expression(sigma),ylab="Density",main=expression(paste("(a) Posterior distribution of ",sigma)))
abline(v=res$fitting$inla$hyperparameters$results$sigma_epsilon$mean,lwd=1.5)
abline(v=c(res$fitting$inla$hyperparameters$results$sigma_epsilon$quant0.025,
           res$fitting$inla$hyperparameters$results$sigma_epsilon$quant0.975),lwd=1.5,col="red")
plot(res$fitting$inla$hyperparameters$posteriors$phi,type="l",xlab=expression(phi),ylab="Density",main=expression(paste("(b) Posterior distribution of ",phi)))
abline(v=res$fitting$inla$hyperparameters$results$phi$mean,lwd=1.5)
abline(v=c(res$fitting$inla$hyperparameters$results$phi$quant0.025,
           res$fitting$inla$hyperparameters$results$phi$quant0.975),lwd=1.5,col="red")
```

<img src="man/figures/reproduce/results/results-simex-1.png" title="plot of chunk simex" alt="plot of chunk simex" width="100%" />

```r
par(mfrow=c(1,1))
#plot results
  
  dsims = colDiffs(rbind(matrix(rep(res$initial$age,ncol(res$simulation$age)),nrow=1), res$simulation$age))
  dmeans = rowMeans(dsims)
  dsd = rowSds(dsims)
ggd = data.frame(y = res$data$dy,depth=res$data$depth,
                   mean=dmeans,
                   lower = dmeans-1.96*dsd,
                   upper = dmeans+1.96*dsd)
  ggp = ggplot(data=ggd,aes(x=depth)) + theme_bw()+ #xlim(rev(range(results$data$depth)))+
    geom_line(aes(y=y),col="gray")+
    geom_line(aes(y=mean),col="blue",size=0.8)+
    geom_ribbon(aes(ymin=lower,ymax=upper),color="red",alpha=0,size=0.8)+
    xlab("Depth (m)")+ylab("Layers per 5cm")+ggtitle("(a) Fitted layer increment model")
  
  
  ggd2 = data.frame(mean=res$simulation$summary_sync$mean-res$data$age,
                    lower=res$simulation$summary_sync$lower-res$data$age,
                    upper=res$simulation$summary_sync$upper-res$data$age,
                    depth=res$data$depth
                    )
  
 reference=res$original.chron$age[1+1:n]
 reftie = reference[res$tie_points$locations_indexes]
ggdtie = data.frame(mean=synchronization$params$mean-reftie,
                    upper=synchronization$params$mean+1.96*synchronization$params$sd-reftie,
                    lower=synchronization$params$mean-1.96*synchronization$params$sd-reftie,
                    locations = res$data$depth[synchronization$locations])
ggp2 = ggplot(data=ggd2,aes(x=depth)) + theme_bw()+
  geom_ribbon(aes(ymin=lower,ymax=upper),color="red",fill="red",alpha=0.3)+
  geom_line(aes(y=mean),color="blue") +
  geom_linerange(data=ggdtie,aes(x=locations,ymin=lower,ymax=upper),color="magenta")+
  xlab("Depth (m)") + ylab("Deviation from simulated reference (years)")+
  ggtitle("(b) Simulated chronologies")
library(gridExtra)
ggsave("fitdiffsim-5500x4000.eps", device=cairo_ps,width=5500,height=4000,units="px",dpi=500,limitsize=FALSE,arrangeGrob(ggp,ggp2))
str = ""
for(i in 1:length(res$fitting$inla$fit$summary.fixed$mean)){
  str=paste0(str," & & ",format(res$fitting$inla$fit$summary.fixed$mean[i],digits=7,scientific=FALSE), " & (",
             format(res$fitting$inla$fit$summary.fixed$`0.025quant`[i],digits=7,scientific=FALSE),", ",
             format(res$fitting$inla$fit$summary.fixed$`0.975quant`[i],digits=7,scientific=FALSE),") \\\\ \n")
}
cat(str)  
#>  & & 0.4588836 & (0.3520685, 0.5656528) \\ 
#>  & & 0.02584957 & (-1.934273, 1.985957) \\ 
#>  & & 0.003184008 & (0.001023886, 0.005337028) \\ 
#>  & & 0.1387808 & (-1.820403, 2.09766) \\ 
#>  & & 0.01248634 & (0.01064573, 0.01443456) \\ 
#>  & & -0.06266468 & (-2.02353, 1.898212) \\ 
#>  & & 0.01020538 & (0.008370453, 0.01208408) \\
```

## Real data example
This is a real data example which produces 5000 simulations from a synchronized time scale.

```r
require(INLA)
library(bremla)
data("event_intervals")
data("events_rasmussen")
data("NGRIP_5cm")
age = NGRIP_5cm$age
depth = NGRIP_5cm$depth
depth2 = depth^2/depth[1]^2 #normalize for stability
proxy = NGRIP_5cm$d18O
data = data.frame(age=age,dy=c(NA,diff(age)),depth=depth,depth2=depth2,proxy=proxy)
formula = dy~-1+depth2
events=list(locations = events_rasmussen$depth,
            locations_unit="depth",degree=1)
results = bremla(formula,data,reference.label="GICC05",
                nsims=nsims,
                events=events,
                synchronization=list(method="adolphi"),
                control.fit=list(method="inla"),
                control.sim=list(synchronized=TRUE) )
summary(results)
#> 
#> Call:
#> knit("inst/reproduce_results/results.Rmd", "inst/reproduce_results/results.md")
#> 
#> Time used:
#>   Model fitting Chron. sampling           Total 
#>          9.1723        184.1855        199.3541 
#> 
#> The fixed component is explained by linear predictor: 
#> dy ~ -1 + depth2 + psi_fill(degree=1, n_events=69)
#> 
#> The noise component is explained by an ar1 process.
#> 
#> The model is fitted using INLA, with following estimates for the hyperparameters:
#>                 mean     sd quant0.025 quant0.25 quant0.5 quant0.75 quant0.975
#> sigma_epsilon 0.4530 0.0025     0.4480    0.4513   0.4530    0.4547     0.4581
#> phi           0.2802 0.0071     0.2662    0.2754   0.2802    0.2850     0.2942
#> 
#> Simulating 10000 chronologies, using GICC05 as reference.
#> 
#> 10000 synchronized chronologies sampled using 4 tie-point distributions (Adolphi).
#> Tie-points are fixed at GICC05 ages (yb2k):
#> 12050.3, 13050.62, 22050, 42049.59.
```

<img src="man/figures/reproduce/results/results-plot-1.png" title="plot of chunk plot" alt="plot of chunk plot" width="100%" />

## Attribution

This code is associated and written for the papers Myrvoll-Nilsen et al. (2022) and Myrvoll-Nilsen et al. (202x) mentioned above. Feel free to use the code, but please cite the accompanying papers.

## License

The code in this repository is made available under the terms of the GNU (version >=2) License. For details, see LICENSE.md file.
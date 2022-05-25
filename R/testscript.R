if(FALSE){
  data("event_intervals")
  data("events_rasmussen")
  data("NGRIP_5cm")

  age = NGRIP_5cm$age
  depth = NGRIP_5cm$depth
  d18O = NGRIP_5cm$d18O
  proxy=d18O

  #load abrupt warming transitions
  # eventdata = read_ods("datasets_used/GISevents.ods")
  # GISevents = eventdata[(eventdata$`NGRIP depth (m)`>min(depth))&(eventdata$`NGRIP depth (m)`<max(depth)),]
  # event_intervals = read.table("datasets_used/event_intervals.txt")



  eventdepths = events_rasmussen$depth
  eventindexes = c(1,which.index(eventdepths, depth[2:length(depth)]) )
  eventindexes = unique(eventindexes[!is.na(eventindexes)])

  #plot d18O proxies as a function of depth and age (GICC05), respectively
  par(mfrow=c(1,1),mar=c(5,4.5,2,2)+0.1)
  plot(depth,d18O,type="l",xlab="Depth (m)",ylab=expression(paste(delta^18,"O (permil)")),xlim=rev(range(depth))); abline(v=eventdepths,lwd=0.7,col="gray")
  plot(age,d18O,type="l",xlab="Age (yb2k)",ylab=expression(paste(delta^18,"O (permil)")),xlim=rev(range(age))); abline(v=age[eventindexes],lwd=0.7,col="gray")


  nsims=10000

  #events used in regression model
  events=eventdepths #locations of transitions used in regression model. Pairs with 'eventmeasure' for finding the corresponding indices
  eventmeasure="depth"
  #reg.model specifies the structure for which the formula string should be created
  reg.model = list(
    const=FALSE,depth1=FALSE,depth2=TRUE,proxy=TRUE,psi0=TRUE,psi1=TRUE); method="inla";
  plots = list(posteriors=TRUE)


  eventnumber=13 #number between 1 and 29. specifies which transition to consider

  #load data window and specifics to transition
  lowerints = which.index(event_intervals$depth_int_lower.m, depth[2:length(depth)])
  upperints = which.index(event_intervals$depth_int_upper.m, depth[2:length(depth)])
  depth.reference = event_intervals$NGRIP_depth_m[eventnumber]
  age.reference = event_intervals$GICC_age.yb2k[eventnumber]
  interval = lowerints[eventnumber]:upperints[eventnumber]



  #list object specifying parameters for abrupt warming transition dating
  event.estimation = list(interval=interval,t1.sims=50000,rampsims=50000,label="GI-11",
                          depth.reference=event_intervals$NGRIP_depth_m[eventnumber],
                          age.reference=event_intervals$NGRIP_agee_m[eventnumber])

  reference.label="GICC05" #this string is used in plot_results

  method="inla" #should INLA be run in fitting procedure? alternatively "LS" for no (not quite as tested)
  CI.type="quantiles" #method used for credible intervals. Alternatives are 'quantiles' and 'hpd'. For Gaussian processes (or any unimodal and symmetric) these are equal


  print.progress=TRUE

  transform = "identity" #set equal to "log" for logarithmic transformation
  noise = "ar1" #iid, ar1 and ar2 supported

  proxy.type="d18O"
  bias = list(bias.model="uniform",biasparams=cbind( c(1,1),c(0.98,1.02),c(0.96,1.04) ),
              store.samples=FALSE)
  store.everything=FALSE #should fixed model components and stochastic both be stored (TRUE), or just the sum (FALSE)

  #run main function wrapper
  synchronization = list(locations=c(11050,12050,13050,22050,42050),locations.type="age",method="adolphi",samples=NULL)
  object = bremla(age,depth,proxy, events=eventdepths,nsims=nsims, eventmeasure = eventmeasure,proxy.type=proxy.type,
                  reference.label=reference.label,transform=transform,reg.model = reg.model,
                  noise=noise, method=method, CI.type=CI.type,
                  synchronization=synchronization,
                  event.estimation = event.estimation,store.everything=store.everything,
                  print.progress=print.progress,bias = bias
  )
  plot(object)
  summary(object)

}

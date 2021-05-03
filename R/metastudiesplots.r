
metastudies_funnel_plot<-function(X,sigma, critval=1.96){
  n=length(X)
  significant<-(abs(X/sigma)>critval)
  nooutlier= (sigma<30*mean(sigma))&(abs(X) < 30*abs(X));
  dat<-data.frame(X, sigma, as.factor(significant&nooutlier))
  names(dat)=c("xvar", "yvar","significant")
  rangeX=1.1*max(max(abs(X)), max(abs(sigma[nooutlier]))*critval)

  dat<-dat[order(dat$significant),]

  ggplot(dat, aes(x=xvar,y=yvar)) +
    xlab("X")+
    ylab(expression(sigma))+
    geom_abline(intercept = 0,slope=1/critval,color="grey")+
    geom_abline(intercept = 0,slope=-1/critval,color="grey")+
    geom_point(size = 4,aes(colour = significant,
                            fill = significant), alpha=min(.8,max(40/n,.3)))+
    #scale_fill_manual(values=c("grey", "blue")) +
    scale_colour_manual(values=c("grey50", "blue")) +
    scale_x_continuous(expand = c(0,0),limits = c(-rangeX,rangeX))+
    scale_y_continuous(expand = c(0,0), limits = c(0,rangeX/critval))+
    theme(legend.position="top",
          #aspect.ratio=1,
          #panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey95", colour = NA))

}

metastudies_abs_funnel_plot<-function(X,sigma, critval=1.96, log.scale=FALSE, remove.outlier = !log.scale){
  X = abs(X)
  sigma = abs(sigma)
  n=length(X)
  significant<-(abs(X/sigma)>critval)


  dat<-data.frame(X, sigma, as.factor(significant))

  if (remove.outlier) {
    outlier = (sigma>30*mean(sigma)) | (abs(X) > 30*mean(abs(X)));
    dat = dat[!outlier,]
  }

  names(dat)=c("X", "sigma","significant")
  rangeX=1.1*max(max(abs(dat$X)), max(abs(dat$sigma))*critval)

  dat<-dat[order(dat$significant),]

  gg = ggplot(dat, aes(x=X,y=sigma)) +
    xlab("X")+
    ylab(expression(sigma))+
    geom_point(size = 4,aes(colour = significant,
      fill = significant), alpha=min(.8,max(40/n,.3)))+
    #scale_fill_manual(values=c("grey", "blue")) +
    scale_colour_manual(values=c("grey50", "blue")) +
    theme(legend.position="top",
          #aspect.ratio=1,
          #panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey95", colour = NA))

  if (!log.scale) {
    gg = gg +
      geom_abline(intercept = 0,slope=1/critval,color="black")+
      scale_x_continuous(expand = c(0,0),limits = c(0,rangeX))+
      scale_y_continuous(expand = c(0,0), limits = c(0,rangeX/critval))
  } else {
    min.x = 0.000000001
    line.X = 10^seq(log(min.x,10), log(rangeX,10), length.out = 501)
    line.sigma = line.X/critval
    line.df = data.frame(X=line.X, sigma=line.sigma)

    gg = gg +
      geom_line(data = line.df, color="black") +
      scale_x_log10(limits = c(min.x,rangeX))+
      scale_y_log10(limits = c(min.x,rangeX/critval))
  }

  gg
}


z_histogram=function(X,sigma){
  Z=X/sigma
  n=length(Z)
  ll=floor(min(Z));
  uu=ceiling(max(Z));

    if (n>=30) {
      uu2<-ceiling(max((uu-.36)/.32,0))*.32+.36;
      ll2<-floor(min((ll+.36)/.32,0))*.32-.36;
      edges<-c(seq(from=ll2,
                   to=-0.36,
                   by=0.32), 0, seq(from=0.36,
                                    to=uu2,
                                    by=0.32));
    } else {
      uu2<-ceiling(max((uu-.68)/.64,0))*.64+.68;
      ll2<-floor(min((ll+.68)/.64,0))*.64-.68;
      edges<-c(seq(from=ll2,
                   to=-0.68,
                   by=0.64), 0, seq(from=0.68,
                                    to=uu2,
                                    by=0.64));
    }


  ggplot(data = as.data.frame(Z), aes(Z))+
    geom_histogram(aes(y = ..density..),
                   fill = 'blue',
                   breaks=edges)+
    geom_vline(xintercept =-1.96,color='grey')+
    geom_vline(xintercept =1.96, color='grey')+
    xlab('Z')+
    ylab('Density')+
    xlim(c(min(edges),max(edges)))+
    theme(#panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey95", colour = NA))

}

#' Plot
estimates_plot<-function(ms, X=ms$X, sigma=ms$sigma, cutoffs=ms$cutoffs, symmetric=ms$symmetric, model=ms$model){
  restore.point("estimates_plot")
  n=500
  psi_hat= ms$psi_hat
  rangeZ=3
  dens=data.frame(z=seq(-rangeZ,rangeZ,length.out =n))
  shift=as.integer(model=="t")

  Tpowers=Tpowers_fun(dens$z,cutoffs,symmetric)
  betap=as.vector(c(psi_hat[-(1:(2+shift))],  1))
  dens$p=Tpowers%*%betap

  if (model=="t") df=psi_hat[3]
    else df=Inf

  dens$f=dt(((dens$z - psi_hat[1])/ psi_hat[2]), df=df)/psi_hat[2]
  names(dens)[names(dens) == 'f'] <- 'density of true effect'
  names(dens)[names(dens) == 'p'] <- 'publication probability'

  dens=melt(dens, id="z")
  ggplot(dens, aes(x=z, y=value)) +
    xlab(paste("Z, ", intToUtf8(952)))+
    geom_line(size=2, color="blue") +
    facet_grid(variable ~ .,  scales = "free_y") +
    expand_limits(y = 0) +
    scale_x_continuous(breaks =-3:3) +
    theme(#panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey95", colour = NA))


}

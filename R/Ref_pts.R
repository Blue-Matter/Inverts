

tabTS = function(yrs, x, p, ref, refdir = "b", sf = 4){

  Mean = apply(x,2,mean,na.rm=T)
  qs = apply(x,2,quantile,p=p,na.rm=T)
  pnams = paste0("P", round(ref*100,1))
  Ps=NULL
  for(i in 1:length(ref))   Ps=cbind(Ps,apply(x<ref[i],2,mean,na.rm=T))
  if(refdir != "b")Ps = 1-Ps
  colnames(Ps) = paste0("P", round(ref*100,1))
  dat = round(data.frame(cbind(Mean,t(qs),Ps)),sf)
  rownames(dat) = yrs
  dat
}


RefPlot = function(x, nam="", col="blue", refcols ="red"){

  if(nam=="")nam= deparse(substitute(x))
  yrs = as.numeric(rownames(x))
  qs = x[,grepl("X", colnames(x))]
  nq = ncol(qs)
  medno = (nq+1)/2
  med = qs[,medno]
  qs = qs[,(1:nq)!=medno]
  Ps = x[,grepl("P",colnames(x)),drop=F]
  Pnams= names(Ps)
  Plevls = as.numeric(sapply(Pnams,function(x)strsplit(x,"P")[[1]][2]))/100

  maxy = max(qs)
  reflines = Plevls * maxy
  matplot(yrs,qs,col="white",xlab="",ylab="",ylim=c(0,maxy*1.025));grid(); abline(h=Plevls,col=refcols,lty=2)
  mtext(nam,2,line=2.8,col=col)
  coly = makeTransparent(col,50)
  nq = ncol(qs)
  for(qq in 1:nq) polygon(c(yrs,rev(yrs)),c(qs[,qq],rev(qs[,nq-qq+1])),col=coly,border=NA)
  lines(yrs,med,col=col)
  matplot(yrs,Ps*maxy,add=T,type="l",lty=1,lwd=2,col=makeTransparent(refcols,99))
  axis(4,at=seq(0,maxy,length.out=11),seq(0,100,10))


}


#' Summarize of biomass and F reference points
#'
#' Tables and figures for visualizing historical biomass and exploitation rate reference points
#'
#' @param Hist An object of class Hist created by runMSE(OM, Hist=T)
#' @param plot Logical, should the results be plotted?
#' @param p Vector of quantiles to plot (Default is 90%, 50% intervals plus median)
#' @param B0ref The SSB0 (unfished spawning biomass) reference levels (defaults to 0.2 and 0.4)
#' @param BMSYref The SSBMSY reference levels (defaults to 0.5 and 1.0)
#' @author T. Carruthers
#' @examples
#' Ref.Points(Hist)
#' @export
Ref.Points=function(Hist, plot = T, p = c(0.05,0.25,0.5,0.75,0.95),
                    B0ref = c(0.2, 0.4), BMSYref = c(0.4, 0.8)){

  yrs =  Hist@OM@CurrentYr - ((Hist@OM@nyears-1):0)
  SSB = apply(Hist@TSdata$SBiomass,1:2,sum)
  SSB0 = Hist@Ref$ReferencePoints$SSB0
  SSBMSY = Hist@Ref$ReferencePoints$SSBMSY

  FM = Hist@TSdata$Find
  FMSY = Hist@Ref$ReferencePoints$FMSY

  outlist = list(SSB_SSB0 = tabTS(yrs, SSB/SSB0, p=p, ref = B0ref),
                 SSB_SSBMSY = tabTS(yrs, SSB/SSBMSY, p=p, ref = BMSYref),
                 F_FMSY = tabTS(yrs, FM/FMSY, p=p, ref=1))

  if(plot){
    par(mfrow=c(3,1),mai=c(0.3,0.6, 0.05, 0.3),omi=c(0.4,0.01,0.01,0.3))
    for(i in 1:length(outlist))    RefPlot(x=outlist[[i]],nam=names(outlist)[i])
    mtext("P(X < level)",4,line=0.4,col="red",outer=T)
    mtext("Historical Year",1,line=0.5,outer=T)

  }

  return(outlist)
}


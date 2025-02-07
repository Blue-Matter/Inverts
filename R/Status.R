
# Status


getStat = function(Hist, digits=4){

  nyears = Hist@OM@nyears
  yrs =  Hist@OM@CurrentYr - ((Hist@OM@nyears-1):0)
  SSB = apply(Hist@TSdata$SBiomass,1:2,sum)
  SSB0 = Hist@Ref$ReferencePoints$SSB0
  SSBMSY = Hist@Ref$ReferencePoints$SSBMSY

  FM = Hist@TSdata$Find
  FMSY = Hist@Ref$ReferencePoints$FMSY

  F01 = Hist@Ref$ByYear$F01_YPR

  outlist = list(SSB_SSB0 = round(SSB[,nyears]/SSB0,digits),
                 SSB_SSBMSY = round(SSB[,nyears]/SSBMSY,digits),
                 F_F01 = round(FM[,nyears]/F01[,nyears],digits),
                 F_FMSY = round(FM[,nyears]/FMSY,digits))
}

StatusPlot = function(x,nam="",col){
  hist(x,xlab="",main="",ylab="",col=col,border='white',breaks=2+ceiling(length(x)^0.4));mtext(nam,1,line=2.2,cex=0.9)
}

#' Summarize current stock status
#'
#' Histograms of spawning biomass and fishing mortatlity rate relative to reference levels
#'
#' @param Hist An object of class Hist created by runMSE(OM, Hist=T)
#' @param plot Logical, should the results be plotted?
#' @param col Color of bars
#' @author T. Carruthers
#' @examples
#' Status(Hist)
#' @export
Status = function(Hist, plot=T, col="cadetblue", digits = 4){

  outlist = getStat(Hist,digits)
  nams = c("SSB relative to unfished", "SSB relative to SSBMSY",
           "Fishing mortality rate relative to F0.1 (apical)",
           "Fishing mortality rate relative FMSY")

  if(plot){
    par(mfrow=c(2,2),mai=c(0.6,0.4, 0.05, 0.05),omi=c(0.01,0.4,0.4,0.01))
    for(i in 1:length(outlist))    StatusPlot(outlist[[i]],nams[i],col=col)
    mtext("Rel. Freq.",2,line=0.4,outer=T)
    mtext(paste0("Status in ", Hist@OM@CurrentYr),line=0.4,outer=T,cex=1.1)
  }

  outlist
}


#' Summarize current stock status according to standard DFO reference regions
#'
#' Distribution of estimated spawning biomass and fishing mortality rate relative to reference levels
#'
#' @param Hist An object of class Hist created by runMSE(OM, Hist=T)
#' @param legpos Position of legend 1 (see ?legend)
#' @param legpos2 Position of legend 2 (see ?legend)
#' @author T. Carruthers
#' @examples
#' DFO_status(Hist)
#' @export
DFO_status=function(Hist, legpos='left',legpos2="bottomleft"){

  outlist = getStat(Hist)
  Yr = Hist@OM@CurrentYr
  par(mfrow = c(1, 1), mai = c(0.7, 0.7, 0.02,    0.02))
  DFO_stat_plot(Br = outlist$SSB_SSBMSY, Fr = outlist$F_FMSY, xlab = "", ylab = "",legpos=legpos,legpos2=legpos2)
  mtext(paste0("Status in ",Hist@OM@CurrentYr), 3, line = 0.8, outer = T, font = 2)
  mtext(paste0("SB/SBMSY (",Yr,")"), 1, line = 2.2, font = 2)
  mtext(paste0("F/FMSY (",Yr,")"), 2, line = 2.7,  font = 2)

}



DFO_stat_plot = function(Br, Fr, xlab = NA, ylab = NA,legpos = "left",legpos2 ="bottomleft"){

  if (is.na(xlab))  xlab = "B/BMSY"
  if (is.na(ylab))  ylab = "F/FMSY"
  Brange <- c(0, quantile(Br, 0.99))
  Frange <- c(0, quantile(Fr, 0.99))
  nsim <- length(Br)
  plot(Brange, Frange, axes = F, col = "white", xlab = "", ylab = "")
  textpos <- Frange[1] + 0.825 * (Frange[2] - Frange[1])
  MSEtool:::add_zones(textpos)
  xp <- pretty(seq(0, max(Brange), length.out = 10))
  yp <- pretty(seq(0, max(Frange), length.out = 10))
  axis(1, xp, xp)
  axis(2, yp, yp)
  mtext(xlab, 1, line = 2.5)
  mtext(ylab, 2, line = 2.5)
  pointcol <- makeTransparent("black", 80)
  linecol <- "black"
  points(Br, Fr, pch = 19, cex = 0.9, col = pointcol)
  MSEtool:::encircle(Br, Fr, col = linecol, xrange = Brange + c(0, 0.1), yrange = Frange + c(0, 0.1), perc = 0.1, lty = 2, lwd = 1.2)
  MSEtool:::encircle(Br, Fr, col = linecol, xrange = Brange + c(0, 0.1), yrange = Frange + c(0, 0.1), perc = 0.5, lwd = 1.2)
  points(quantile(Br, 0.5), quantile(Fr, 0.5), col = "black",pch = 3, cex = 1.2, lwd = 2)
  fracs <- c(sum(Br < 0.4 & Fr > 1), sum(0.4 < Br & Br < 0.8 &
                                           Fr > 1), sum(0.8 < Br & Fr > 1), sum(Br < 0.4 & Fr <
                                                                                  1), sum(0.4 < Br & Br < 0.8 & Fr < 1), sum(0.8 < Br &
                                                                                                                               Fr < 1))
  fracs <- round(fracs/nsim * 100, 1)
  fposH <- Frange[1] + 0.99 * (Frange[2] - Frange[1])
  fposL <- Frange[1] + 0.02 * (Frange[2] - Frange[1])
  text(0.2, fposH, paste(fracs[1], "%"), col = "white", cex = 0.9, font = 2)
  text(0.6, fposH, paste(fracs[2], "%"), col = "grey73", cex = 0.9,font = 2)
  text(1.1, fposH, paste(fracs[3], "%"), col = "grey73", cex = 0.9, font = 2)
  text(0.2, fposL, paste(fracs[4], "%"), col = "white", cex = 0.9, font = 2)
  text(0.6, fposL, paste(fracs[5], "%"), col = "grey73", cex = 0.9,font = 2)
  text(1.1, fposL, paste(fracs[6], "%"), col = "grey73", cex = 0.9,font = 2)
  legend(legpos, legend = c("A simulation", "Median"), pch = c(19,       3), col = pointcol, cex = 0.9)
  legend(legpos2, legend = c("50%", "90%"), lty = c(1, 2),      col = linecol, bty = "n", cex = 0.9)

}


#' Compare stock status among operating models relative to standard DFO reference regions
#'
#' Distribution of estimated spawning biomass and fishing mortality rate relative to reference levels
#'
#' @param Histlist A a list of objects of class Hist created by runMSE(OM, Hist=T)
#' @param p Vector of fractions (percentiles) for plotting bars
#' @param abbrev_no Positive integer controlling abbreviation (length) of y axis labels (OM names)
#' @param pad numerical value of additional (absolute y scale) padding for percentage probabilies on plot
#' @param legpad numerical value (y axis multiplier) of additional legend padding
#' @param boxwd positive numerical value of width of boxes (absolute y scale)
#' @author T. Carruthers
#' @examples
#' DFO_status_comp(list(Hist,Hist2))
#' @export
DFO_status_comp=function(Histlist,p=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),abbrev_no=6,pad=0.02,legpad=0.4, boxwd=0.08){

  Yr = Histlist[[1]]@OM@CurrentYr
  #Histlist = list(Hist, Hist2)
  par(mai=c(0.8,0.6,0.1,0.1))
  nams=abbreviate(sapply(Histlist,function(x)x@OM@Name),abbrev_no)

  Brels = lapply(Histlist,function(x)getStat(x)$SSB_SSBMSY)
  qs = t(sapply(Brels,quantile,p=p))
  PH = t(sapply(Brels,function(x){c(mean(x<0.4),mean(x>0.4 & x<0.6),mean(x>0.8))}))
  mus = sapply(Brels,mean)
  nc = nrow(qs)
  miny=2
  locs = miny+1:nc-0.5
  plot(c(0,max(qs)*1.025),c(miny,max(locs)*(1+legpad)),col='white',axes=F,xlab="",ylab="",yaxs ="i",xaxs ="i")
  textpos=max(locs)*(1+legpad*0.66)
  MSEtool:::add_zones(textpos);grid()
  nq = length(p)
  ni = nrow(qs)
  lwds=c(3,6,10)
  for(i in 1:ni){
    for(qq in 1:2)lines(c(qs[i,qq],qs[i,nq-qq+1]),rep(locs[i],2),lwd=lwds[qq])
    midp = ceiling(nq/2)
    lines(rep(qs[i,midp],2),locs[i]+c(-boxwd,boxwd),lwd=3,col='red')
    polygon(c(rep(qs[i,3],2),rep(qs[i,5],2)),locs[i]+c(-boxwd,boxwd,boxwd,-boxwd),lwd=2)
    lines(rep(mus[i],2),locs[i]+c(-boxwd,boxwd)*2.3,lwd=2,lty=3,col='red')
    text(c(0.2,0.6,1.1),locs[i]+boxwd*2+pad,paste0(round(PH[i,]*100,2),"%"),cex=0.7,font=2,col=c("white","darkgrey","darkgrey"))
  }
  axis(1,c(-1E10,1E10),rep("",2)); axis(1)
  axis(2,c(-1E10,1E10),rep("",2)); axis(2,locs,nams,cex=0.7)
  mtext(paste0("SSB/SSBMSY (",Yr,")"),1,line=2.2)
  col=c(rep('black',3),rep('red',2))
  legend('topright',lwd=c(lwds[1:2],1,3,2),lty=c(1,1,1,1,3),
         legend=c("95%","90%","50%","Med","Mean"),col=col,text.col=col,bty="n",cex=0.9)

  colnames(PH) = c("P(Crit.)","P(Caut.)","P(Healthy)")
  rownames(PH) = nams
  PH
}

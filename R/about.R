
#' Get a summary of OM properties
#'
#' A plot and optional list of properties
#'
#' @param obj An object of class In or OM or RCMfit
#' @param xlim Vector c(min, max) of longitudinal range to plot
#' @param ylim Vector c(min, max) of latitudinal range to plot
#' @param xfac Positive numerical value (probably above 1) that stretches the white space on the x axis around the plot
#' @param col Color of shaded area
#' @param plot.fig Logical Should the figure be plotted>
#' @param ret.dat Logical Should the summary data be returned?
#' @param p The quantiles to be calculated e.g. p=c(0.005, 0.995) is the 99 percent interval
#' @param plot.MAs Logical, should the management areas (statistical areas) be plotted (slower)?
#' @author T. Carruthers
#' @examples
#' about(OM.GD.7)
#' @export
about = function(obj, xlim = c(-135,-125.5), ylim = c(48, 55.5), xfac=1.7, col="red",
                 plot.fig=T, ret.dat=F, p = c(0.05, 0.95), plot.MAs=T){

  #old.par <- par(no.readonly = TRUE)
  catchy = Fd = NA

  if(class(obj)=="OM"){
    OM=obj
    if("Data"%in%names(OM@cpars))catchy=OM@cpars$Data@Cat[1,]
    if("Find"%in%names(OM@cpars))Fd = OM@cpars$Find
  }

  if(class(obj)=="In"){
    OM=obj[[2]]
    catchy = obj[[3]]@Chist
  }

  if(class(obj)=="RCMfit"){
    OM=obj@OM
    catchy=OM@cpars@Data@Cat[1,]
    Fd = OM@cpars$Find
  }
  yrs = OM@CurrentYr-(OM@nyears-1):0

  if(plot.fig){

    par(mai =rep(0.1,4))
    xfac = 1.7
    xd = xlim[2]-xlim[1]
    yd = ylim[2]-ylim[1]
    xmid = mean(xlim)
    allxlim = xmid+c(-1,1) *(xd/2*xfac)
    plot(1,col="white",axes=F,xlim=allxlim,ylim=ylim);grid()

    poly = OM@Misc$mapinfo$poly

    if(class(poly)[1]=="data.frame"){
      polygon(poly[,1],poly[,2],col=col,border=col)
    }else{ #  if(class(poly)[2]=="sfc")
      plot(poly, add=T, col=col)
    }

    if(plot.MAs)plot(MA_polys,add=T)

    map(xlim = xlim, ylim = ylim, add=T,fill=T,col=rep("lightgrey",100))
    #cnt = as.numeric(st_centroid(poly[[1]]))
  }
  if(!is.null(OM@cpars))cpars = OM@cpars
  if(is.null(OM@cpars))cpars = list()

  # OM properties
  cparnams = names(cpars)
  onam = c("M","K", "Linf","R0","Depln.")
  snam = c("M_ageArray","K","Linf","R0","D")
  unit = c(1,          1,    1, 1E6, 1)
  sigf = c(3,          2,    1, 2, 2)
  ns = length(onam)
  levs = (1:ns)[snam%in%cparnams]
  for(ss in levs)    if(snam[ss]%in%cparnams) assign(onam[ss], round(quantile(cpars[[snam[ss]]]/unit[ss], p),sigf[ss]))
  slabs = sapply(levs,function(X,onam,snam){paste(onam[X],"=",get(onam[X])[1],"-",get(onam[X])[2])},snam=snam,onam=onam)
  yrng = range(yrs)
  ylab = paste0(yrng[1]," - ", yrng[2])

  if(plot.fig){
    wtbg = function()polygon(c(-1E10,1E10,1E10,-1E10),c(-1E10,-1E10,1E10,1E10),col="white",border="white")
    doax1 = function(){axis(1,c(-1E10,1E10));axis(4,c(-1E10,1E10))}
    doax2 = function(){axis(1,190:300*10,rep("",111));axis(2);axis(1,c(-1E10,1E10));axis(2,c(-1E10,1E10));axis(3,c(-1E10,1E10));axis(4,c(-1E10,1E10))}
    legend('topright', legend = c(OM@Name, ylab, slabs), bty="n", text.col=c(col,rep('black',length(slabs)+1)), cex=1)
    # Data properties
    if(!is.na(catchy[1])){
      par(fig = c(0,0.36,0,0.25),mai = c(0.4,0.7,0.1,0.05),new=T)
      plot(yrs,catchy,pch=19,col="white",ylab="",xlab="",ylim=c(0,max(catchy)*1.025)); wtbg(); doax1()
      points(yrs,catchy,pch=19);grid(); lines(yrs,catchy,col=col)
      mtext("Catch",2,line=2.2)
    }

    if(!is.na(Fd[1])){
      par(fig = c(0,0.36,0.25,0.46),mai = c(0.1,0.7,0.1,0.05),new=T)
      Fq = apply(Fd,2,quantile,c(0.05,0.25,0.5,0.75,0.95))
      plot(range(yrs),c(0,max(Fq)*1.025),col="white",axes=F,xlab="",ylab="");wtbg();doax2();grid()
      matplot(yrs,t(Fq),col="NA",xlab="",ylab="",ylim=c(0,max(Fq)*1.025),add=T)
      coly = makeTransparent(col,50)
      polygon(c(yrs,rev(yrs)),c(Fq[1,],rev(Fq[5,])),col=coly,border=NA)
      polygon(c(yrs,rev(yrs)),c(Fq[2,],rev(Fq[4,])),col=coly,border=NA)
      lines(yrs,Fq[3,],col=col)
      mtext("Fish. Mort.",2,line=2.2)
    }

    par(fig=c(0,1,0,1))
  } # end of plot.fig

  if(ret.dat){
    outlist = list()
    outlist[[1]] = OM@Name
    outlist[[2]] = yrng
    names(outlist)[1:2] = c("Name","Years")
    j =2
    for(ss in levs){
      j = j +1
      outlist[[j]] = get(onam[ss])
      names(outlist)[j] = onam[ss]
    }
    return(outlist)
  }

}

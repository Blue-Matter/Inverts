
smoothy<-function(xx,plot=F,enp_mult,plotname="",xlab="x",ylab="y",x = NA){
  tofill<-!is.na(xx)
  xx[xx==0]<-1E3
  if(sum(tofill)>3){
    predout<-rep(NA,length(xx))
    if(is.na(x[1]))x = 1:length(xx)
    dat<-data.frame(x=x,y=log(xx))
    enp.target<-sum(tofill)*enp_mult
    out<-loess(y~x,dat=dat,enp.target=enp.target)
    predout[tofill]<-exp(predict(out))
    if(plot){
      plot(x,xx,type="p",xlab=xlab,ylab=ylab,main=plotname); grid()
      lines(x,predout,col="#ff000090",lwd=2)
      legend('topright',c("Observed","Smoothed"),text.col=c("black","red"),bty='n')
    }
  }else{
    predout = xx
  }
  predout
}

I_obs_freq = function(I_hist, I_freq, LHYr, CurYr, Year){ # filter out data to match futureu sampling frequency
  nI = nrow(I_hist)
  ny = ncol(I_hist)
  nkeep = sum(I_freq!=0)
  I_keep = I_smth = array(NA,c(nkeep,ny))
  Is = (1:nI)[I_freq>0]

  j = 0
  for(i in Is){
    j = j+1
    Ivec = I_hist[i,]
    pind = match((LHYr+1):CurYr,Year)
    if(!(is.na(pind[1]))){
      np = length(pind)
      makeNA = rep(c(rep(TRUE,I_freq[i]-1),FALSE),100)[1:np]
      Ivec[pind[makeNA]] = NA
      I_keep[j,] = Ivec
    }else{
      I_keep[j,] = Ivec
    }
  }
  I_keep
}


doRec = function(Rec, MPrec, mod, TACdelta, TACrng){

  # TAC change
  maxchng = TACdelta[2]
  minchng = TACdelta[1]

  if(mod > (1+maxchng)){
    mod = 1+maxchng    # above max TAC increase
  }else if(mod < (1-maxchng)){
    mod = 1-maxchng    # below max TAC reduction
  } else{
    mod = 1  # within min TAC change
  }

  # TAC constraint
  TrialTAC = MPrec * mod

  if(TrialTAC > TACrng[2]){
    Rec@TAC = TACrng[2] # max TAC
  }else if(TrialTAC < TACrng[1]){
    Rec@TAC = TACrng[1] # min TAC
  }else{
    Rec@TAC = TrialTAC
  }

  Rec
}

doHCR = function(trial_TAC, est, ref, CP = c(0,1), CPy = c(0,1)){
  lev = est / ref
  if(lev <= CP[1]) TAC = trial_TAC*CPy[1]
  if(lev >= CP[2]) TAC = trial_TAC*CPy[2]
  if(lev > CP[1] & lev < CP[2]) TAC = (trial_TAC*CPy[1])+trial_TAC * (CPy[2]-CPy[1])*(lev - CP[1]) / (CP[2] - CP[1])
  TAC
}


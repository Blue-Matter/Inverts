# OM_mod



#' Modify natural mortality rate and recruitment
#'
#' A wrapper function that modifies the cpars slot of the OM object to control the value and trajectory of natural mortality rate and recruitment strength for all years
#'
#' @param OM An object of class OM
#' @param M_val A single value of M, an lower/upper range of M (uniformly sampled) or a per-simulation specification of M (vector nsim long)
#' @param M_trend A single value that is the annual percentage change in M (a value of 1 is a 1 percent increase per year in all years) or a custom pattern in M for each year that is a factor (1 is no change)
#' @param R_trend A single value that is the annual percentage change in recruitment strength (a value of 1 is a 1 percent increase per year in all years) or a custom pattern in M for each year that is a factor (1 is no change)
#' @author T. Carruthers
#' @examples
#' OM_mod(OM.MC.E, M_val = c(0.2,0.3), M_trend = 1, Rtrend = -1)
#' @export
OM_mod = function(OM, M_val=NaN, M_trend = 0, R_trend = 0){

  nsim = OM@nsim
  py = OM@proyears
  ny = OM@nyears
  ay = py+ny
  na = OM@maxage
  Mdims = c(nsim ,na+1, ay)
  yind = 1:(ny+py)

  if(class(OM)!="OM")stop("The first argument must be an object of class OM")

  # M_val --------------------------------------------------

  if(!is.na(M_val[1])) {  # M_val is specified
    lM = length(M_val)
    if(!(lM%in%c(1,2,nsim,na+1)))stop("M_val must be either a single value, a range (vector length 2), a per sim vector (length OM@nsim),or an age vector (length OM@maxage + 1)")
       if(lM ==1){
         OM@cpars$M_ageArray = array(M_val ,Mdims)
         OM@cpars$M = NULL
         OM@M = rep(M_val,2)
         cat(paste0("A single value of age- and time-invariant M was specified: ",M_val," \n"))
       }else if(lM==2){
         Msamp = runif(nsim, M_val[1],M_val[2])
         OM@cpars$M_ageArray = array(Msamp, Mdims)
         OM@cpars$M = NULL
         OM@M = M_val
         cat(paste0("Age- and time-invariant M was specified in the range of ",M_val[1]," - ",M_val[2]," \n"))
       }else if(lM ==nsim){
         OM@cpars$M_ageArray = array(M_val, Mdims)
         OM@cpars$M = NULL
         OM@M = range(M_val)
         cat(paste0("Per simulation, age- and time-invariant M was specified in the range of ",round(min(M_val),4)," - ",round(max(M_val),4)," \n"))
       }else if(lM == na+1){
         OM@cpars$M_ageArray = array(rep(M_val,each=nsim), Mdims)
         OM@cpars$M = NULL
         OM@M = range(M_val)
         cat(paste0("Age-specific, time-invariant M was specified in the range of ",min(M_val)," - ",max(M_val)," \n"))
       }
  } # if M_val is specified

  if(is.null(OM@cpars$M_ageArray)){ # default level for M_ageArray
    if(is.null(OM@cpars$M)){
      Msamp = runif(nsim,OM@M[1],OM@M[2])
      OM@cpars$M_ageArray = array(Msamp, Mdims)
    }else{
      OM@cpars$M_ageArray = array(OM@cpars$M, Mdims)
    }
  }
  OM@cpars$M = NULL


  # Mtrend ------------------------------------

  if(!(length(M_trend)%in%c(1,ay))) stop("M_trend must either be a single value (% annual change in M) or must be a factor vector n historical years + n projection years long (length OM@nyears + OM@proyears)")
  if(length(M_trend)!=1 & M_trend[1]!=0){
    if(length(M_trend)==1){
      fac = (1+M_trend/100)^(1:ay)
      OM@cpars$M_ageArray[] = OM@cpars$M_ageArray[] * array(rep(fac,each=nsim*(na+1)), Mdims)
      cat(paste0("Natural mortality rate modified to change by ",M_trend, " per cent per year in all years \n"))
    }else{
      OM@cpars$M_ageArray[] = OM@cpars$M_ageArray[] * array(rep(M_trend,each=nsim*(na+1)),Mdims)
      cat(paste0("Natural mortality rate modified to change by a custom specified trend in all years \n"))
    }
  }

  # R_trend -----------------------------------

  if(!(length(R_trend)%in%c(1,py+na+ny))) stop("R_trend must either be a single value (% annual change in M) or must be a factor vector nprojection years long (length OM@proyears)")

  rind = 1:(py+na+ny)
  if(length(R_trend)!=1 & R_trend[1]!=0){
    if(length(R_trend)==1){
      fac = (1+R_trend/100)^rind
      OM@cpars$Perr_y[] = OM@cpars$Perr_y[] * array(rep(fac,each=nsim), c(nsim,py+na+ny))
      cat(paste0("Recruitment strength modified to change by an additional ",R_trend, " per cent per year in all years \n"))
    }else{
      OM@cpars$Perr_y[] = OM@cpars$Perr_y[] * array(rep(R_trend,each=nsim), c(nsim,py+na+ny))
      cat(paste0("Recruitment strength modified to change by an additional custom specified trend in all years \n"))
    }
  }

  return(OM)

}



#' Modify future natural mortality rate and recruitment in projection years
#'
#' A wrapper function that modifies the cpars slot of the OM object to control the value and trajectory of natural mortality rate and recruitment strength for projection years.
#'
#' @param OM An object of class OM
#' @param M_val A single value of M, an lower/upper range of M (uniformly sampled) or a per-simulation specification of M (vector nsim long)
#' @param M_trend A single value that is the annual percentage change in M (a value of 1 is a 1 percent increase per year in projection years) or a custom pattern in M for each projection year that is a factor (1 is no change)
#' @param R_trend A single value that is the annual percentage change in recruitment strength (a value of 1 is a 1 percent increase per year in projection years) or a custom pattern in M for each projection year that is a factor (1 is no change)
#' @author T. Carruthers
#' @examples
#' OM_proj_mod(OM.MC.E, M_val = c(0.2,0.3), M_trend = 1, Rtrend = -1)
#' @export
OM_proj_mod = function(OM, M_val=NaN, M_trend = 0, R_trend = 0){

  nsim = OM@nsim
  py = OM@proyears
  ny = OM@nyears
  ay = py+ny
  na = OM@maxage
  Mdims = c(nsim ,na+1, ay)
  pind = ny+1:py

  if(class(OM)!="OM")stop("The first argument must be an object of class OM")

  # M_val --------------------------------------------------

  if(!is.na(M_val[1])) {  # M_val is specified
     lM = length(M_val)
     if(!(lM%in%c(1,2,nsim,na+1))) stop("M_val must be either a single value, a range (vector length 2), a per sim vector (length OM@nsim),or an age vector (length OM@maxage + 1)")
     if(lM ==1){
       OM@cpars$M_ageArray[,,pind] = array(M_val ,c(nsim,na+1,py))
       OM@cpars$M = NULL
       OM@M = rep(M_val,2)
       cat(paste0("A single value of age- and time-invariant M was specified: ",M_val," \n"))
     }else if(lM==2){
       Msamp = runif(nsim, M_val[1],M_val[2])
       OM@cpars$M_ageArray[,,pind] = array(Msamp, c(nsim,na+1,py))
       OM@cpars$M = NULL
       OM@M = M_val
       cat(paste0("Age- and time- invariant M was specified in the range of ",M_val[1]," - ",M_val[2]," \n"))
     }else if(lM ==nsim){
       OM@cpars$M_ageArray[,,pind] = array(M_val, c(nsim,na+1,py))
       OM@cpars$M = NULL
       OM@M = range(M_val)
       cat(paste0("Per simulation, age- and time-invariant M was specified in the range of ",round(min(M_val),4)," - ",round(max(M_val),4)," \n"))
     }else if(lM == na+1){
       OM@cpars$M_ageArray[,,pind] = array(rep(M_val,each=nsim), c(nsim,na+1,py))
       OM@cpars$M = NULL
       OM@M = range(M_val)
       cat(paste0("Age-specific, time-invariant M was specified in the range of ",min(M_val)," - ",max(M_val)," \n"))
     }
  } # if M_val is specified

  if(is.null(OM@cpars$M_ageArray)){ # default level for M_ageArray
    if(is.null(OM@cpars$M)){
      Msamp = runif(nsim,OM@M[1],OM@M[2])
      OM@cpars$M_ageArray = array(Msamp, Mdims)
    }else{
      OM@cpars$M_ageArray = array(OM@cpars$M, Mdims)
    }
  }
  OM@cpars$M = NULL


  # Mtrend ------------------------------------

  if(!(length(M_trend)%in%c(1,py))) stop("M_trend must either be a single value (% annual change in M) or must be a factor vector nprojection years long (length OM@proyears)")

  if(length(M_trend)!=1 & M_trend[1]!=0){
    if(length(M_trend)==1){
      fac = (1+M_trend/100)^(1:py)
      OM@cpars$M_ageArray[,,pind] = OM@cpars$M_ageArray[,,pind] * array(rep(fac,each=nsim*(na+1)), c(nsim,na+1,py))
      cat(paste0("Natural mortality rate modified to change by ",M_trend, " per cent per year in the projection \n"))
    }else{
      OM@cpars$M_ageArray[,,pind] = OM@cpars$M_ageArray[,,pind] * array(rep(M_trend,each=nsim*(na+1)),c(nsim, na+1, py))
      cat(paste0("Natural mortality rate modified to change by a custom specified trend in the projection years \n"))
    }
  }

  # R_trend -----------------------------------

  if(!(length(R_trend)%in%c(1,py))) stop("R_trend must either be a single value (% annual change in M) or must be a factor vector nprojection years long (length OM@proyears)")

  if(length(R_trend)!=1 & R_trend[1]!=0){
    rind = na+ny+1:py
    if(length(R_trend)==1){
      fac = (1+R_trend/100)^(1:py)
      OM@cpars$Perr_y[,rind] = OM@cpars$Perr_y[,rind] * array(rep(fac,each=nsim), c(nsim,py))
      cat(paste0("Recruitment strength modified to change by an additional ",R_trend, " per cent per year in the projection \n"))
    }else{
      OM@cpars$Perr_y[,rind] = OM@cpars$Perr_y[,rind] * array(rep(R_trend,each=nsim), c(nsim,py))
      cat(paste0("Recruitment strength modified to change by an additional custom specified trend in the projection years \n"))
    }
  }

  OM

}

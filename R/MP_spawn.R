

#' Spawns multiple MPs from a single MP
#'
#' An 'illegal' R function that allows you to rapidly make a large array of MPs
#'
#' @param MPname The name of a function of class MP
#' @param spawnlist A list of argument names, each list are argument values nMPs long
#' @param spawnnames A character vector of MP names
#' @param defaults A list of default values for other arguments #' @author T. Carruthers
#' @export
spawnMP = function(MPname, spawnlist, spawnnames = NA, defaults=list(NA)){

  template = get(MPname)
  forms = names(formals(template))

  if(!all(names(defaults)%in%forms))stop(paste("At least one of the named defaults is not a formal of MP",MPname))

  if(!is.na(defaults[[1]][1])){ # only if there are defaults
    for(i in 1:length(defaults))formals(template)[names(defaults)[i]]=defaults[i] # assign new defaults
  }

  if(is.na(spawnnames[1])){
    combos = lapply(spawnlist)
    spawnnames = paste(MPname,names(spawnlist)[1],spawnlist[[1]],sep="_")
    cat(paste0("Inventing spawnnames (you might want to add an argument):", paste(spawnnames,collapse=", "),"/n"))
  }

  spawncheck = sapply(spawnlist,length)
  if(!all(spawncheck==mean(spawncheck)))stop(paste("At least one of the named spawnlist is not a formal of MP",MPname))

  ns = length(spawnlist[[1]])
  nt = length(spawnlist)

  for(i in 1:ns){ # levels of argument

    temp2 = template

    for(tt in 1:nt){ # types of argument

      formchng = names(spawnlist)[tt]
      if(!(formchng%in%forms))stop(paste("The spawnlist name",formchng, "does not match a formal of formal of MP",MPname))
      formals(temp2)[formchng] = spawnlist[[tt]][i]

      if(tt ==nt){ # final assignment
        class(temp2) = "MP"
        assign(spawnnames[i], temp2, envir=.GlobalEnv) # !!! not cran happy?
        #assign(spawnnames[i],temp2, envir = environment())
      }
    } # end of type
  } # end of spawn

  #for(i in 1:ns)try(do.call(spawnnames[i],args=list(x=1,Data=NULL,reps=1)),silent=T)


}

spawnSC = function(spawnlist,spawnnames = NA){

  snams = names(spawnlist)
  ns = length(snams)
  ni = length(spawnlist[[1]])
  for(i in 1:ni){
    BSD = 0.2; BAC = 0.5; Bbias = 1; B1 = 0; B2 = 1; F1 = 0; F2 = 1
    for(ss in 1:ns) assign(snams[ss],spawnlist[[ss]][i])
    assign(spawnnames[i],
           make_MP(Shortcut,HCR_segment,
              B_err = c(BSD, BAC, Bbias),
              OCP_type = "SSB_SSBMSY",
              Ftarget_type = "FMSY",
              OCP = c(B1, B2),
              relF = c(F1, F2)),
           envir=.GlobalEnv)

  }

}

spawnSP = function(spawnlist,spawnnames = NA){

  snams = names(spawnlist)
  ns = length(snams)
  ni = length(spawnlist[[1]])
  for(i in 1:ni){
   B1 = 0; B2 = 1; F1 = 0; F2 = 1
    for(ss in 1:ns) assign(snams[ss],spawnlist[[ss]][i])
    assign(spawnnames[i],
           make_MP(SP_SS,HCR_segment,
                   OCP_type = "SSB_SSBMSY",
                   Ftarget_type = "FMSY",
                   OCP = c(B1, B2),
                   relF = c(F1, F2)),
           envir=.GlobalEnv)

  }

}


spawnSCA = function(spawnlist,spawnnames = NA){

  snams = names(spawnlist)
  ns = length(snams)
  ni = length(spawnlist[[1]])
  for(i in 1:ni){
    B1 = 0; B2 = 1; F1 = 0; F2 = 1
    for(ss in 1:ns) assign(snams[ss],spawnlist[[ss]][i])
    assign(spawnnames[i],
           make_MP(SCA,HCR_segment,
                   OCP_type = "SSB_SSBMSY",
                   Ftarget_type = "FMSY",
                   OCP = c(B1, B2),
                   relF = c(F1, F2)),
           envir=.GlobalEnv)

  }

}

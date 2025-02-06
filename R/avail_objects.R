
#' What objects of this class are available
#'
#' Generic class finder Finds objects of the specified class in the global environment or the
#' Inverts packages.
#'
#' @param classy A class of object (character string, e.g. 'In', 'OM', 'MSE', 'MP')
#' @param package Optional. Names(s) of the package to search for object of class `classy`. String
#' Default is all `Inverts` packages (Inverts, Inverts.GD, Inverts.MC, Inverts.GSU, Inverts.RSC). Always searches the global environment as well.
#' @param msg Print messages?
#' @examples
#' objs("In", msg=FALSE)
#' @author T. Carruthers
#' @seealso \link{avail}
#' @examples
#' RCM_input_files <- objs("In")
#' OMs <- objs("OM")
#' @export
objs = function (classy, package = NULL, msg = TRUE) {
  temp <- try(class(classy), silent = TRUE)
  if (methods::is(temp, "try-error"))
    classy <- deparse(substitute(classy))
  if (temp == "function")
    classy <- deparse(substitute(classy))
  else {
    packages <- c("Inverts", "Inverts.GD", "Inverts.MC", "Inverts.GSU", "Inverts.RSC")
    if (is.null(package)) {
      package <- packages
      pkgs <- search()
      search_package <- paste0("package:", package)
      package <- package[search_package %in% pkgs]
    }
    global_funs <- ls(envir = .GlobalEnv)[vapply(ls(envir = .GlobalEnv),
                                                 getclass.inverts, logical(1), classy = classy)]
    temp <- global_funs

    if ('Inverts' %in% package) {
      MSEtool_funs <- getfuncs.inverts('Inverts', classy, msg)
      temp <- c(temp, MSEtool_funs)
    }
    if ('Inverts.GD' %in% package) {
      SAMtool_funs <- getfuncs.inverts('Inverts.GD', classy, msg)
      temp <- c(temp, SAMtool_funs)
    }
    if ('Inverts.MC' %in% package) {
      DLMtool_funs <- getfuncs.inverts('Inverts.MC', classy, msg)
      temp <- c(temp, DLMtool_funs)
    }
    if ('Inverts.GSU' %in% package) {
      MSEextra_funs <- getfuncs.inverts('Inverts.GSU', classy, msg)
      temp <- c(temp, MSEextra_funs)
    }
    if ('Inverts.RSC' %in% package) {
      MSEextra_funs <- getfuncs.inverts('Inverts.RSC', classy, msg)
      temp <- c(temp, MSEextra_funs)
    }

    packagex <- package[!package %in% packages]
    if (length(packagex) > 0) {
      other <- sapply(1:length(packagex), function(i) get_funcs(packagex[i],
                                                                classy, msg))
      other <- unlist(other)
      temp <- c(temp, other)
    }
    if (length(temp) < 1)
      stop("No objects of class '", classy, "' found",
           call. = FALSE)
    return(unique(temp))
  }
}

getclass.inverts = function (x, classy){
  return(any(class(get(x)) == classy))
}

getfuncs.inverts = function (package, classy, msg){
  pkgs <- search()
  search_package <- paste0("package:", package)
  funs <- NULL
  if (search_package %in% pkgs) {
    if (msg)
      message("Searching for objects of class ", classy,
              " in package: ", package)
    funs <- ls(search_package)[vapply(ls(search_package),
                                      getclass.inverts, logical(1), classy = classy)]
  }
  else {
    stop("Package ", package, " not loaded. Use `library(",
         package, ")`", call. = FALSE)
  }
  funs
}

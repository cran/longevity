
.rlang_s3_register_compat <- function(fn, try_rlang = TRUE) {
  # Compats that behave the same independently of rlang's presence
  out <- switch(
    fn,
    is_installed = return(function(pkg) requireNamespace(pkg, quietly = TRUE))
  )

  # Only use rlang if it is fully loaded (#1482)
  if (try_rlang &&
      requireNamespace("rlang", quietly = TRUE) &&
      environmentIsLocked(asNamespace("rlang"))) {
    switch(
      fn,
      is_interactive = return(rlang::is_interactive)
    )

    # Make sure rlang knows about "x" and "i" bullets
    if (utils::packageVersion("rlang") >= "0.4.2") {
      switch(
        fn,
        abort = return(rlang::abort),
        warn = return((rlang::warn)),
        inform = return(rlang::inform)
      )
    }
  }

  # Fall back to base compats

  is_interactive_compat <- function() {
    opt <- getOption("rlang_interactive")
    if (!is.null(opt)) {
      opt
    } else {
      interactive()
    }
  }

  format_msg <- function(x) paste(x, collapse = "\n")
  switch(
    fn,
    is_interactive = return(is_interactive_compat),
    abort = return(function(msg) stop(format_msg(msg), call. = FALSE)),
    warn = return(function(msg) warning(format_msg(msg), call. = FALSE)),
    inform = return(function(msg) message(format_msg(msg)))
  )

  stop(sprintf("Internal error in rlang shims: Unknown function `%s()`.", fn))
}


# Internal from the 'vctrs' package
.s3_register <- function (generic, class, method = NULL){
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)
  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]
  caller <- parent.frame()
  get_method_env <- function() {
    top <- topenv(caller)
    if (isNamespace(top)) {
      asNamespace(environmentName(top))
    }
    else {
      caller
    }
  }
  get_method <- function(method) {
    if (is.null(method)) {
      get(paste0(generic, ".", class), envir = get_method_env())
    }
    else {
      method
    }
  }
  register <- function(...) {
    envir <- asNamespace(package)
    method_fn <- get_method(method)
    stopifnot(is.function(method_fn))
    if (exists(generic, envir)) {
      registerS3method(generic, class, method_fn, envir = envir)
    }
    else if (identical(Sys.getenv("NOT_CRAN"), "true")) {
      warn <- .rlang_s3_register_compat("warn")
      warn(c(sprintf("Can't find generic `%s` in package %s to register S3 method.",
                     generic, package), i = "This message is only shown to developers using devtools.",
             i = sprintf("Do you need to update %s to the latest version?",
                         package)))
    }
  }
  setHook(packageEvent(package, "onLoad"), function(...) {
    register()
  })
  is_sealed <- function(pkg) {
    identical(pkg, "base") || environmentIsLocked(asNamespace(pkg))
  }
  if (isNamespaceLoaded(package) && is_sealed(package)) {
    register()
  }
  invisible()
}

.onLoad <- function(...) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
     .s3_register("ggplot2::autoplot", "elife_par")
     # .s3_register("ggplot2::autoplot", "elife_hazard")
     .s3_register("ggplot2::autoplot", "elife_northropcoleman")
     .s3_register("ggplot2::autoplot", "elife_tstab")
     .s3_register("ggplot2::autoplot", "elife_profile")
     # .s3_register("ggplot2::autoplot", "elife_ecdf")
     # .s3_register("ggplot2::autoplot", "elife_npar")
  }
}

#' @rdname plot.elife_par
#' @export
#' @param object an object of class \code{elife_par} containing the fitted parametric model
autoplot.elife_par <- function(object, ...){
  args <- list(...)
  args$plot.type <- "ggplot"
  args$x <- object
  do.call(plot.elife_par, args = args)
}

# #' @rdname plot.elife_hazard
# #' @export
# autoplot.elife_hazard <- function(object, ...){
#   args <- list(...)
#   args$plot.type <- "ggplot"
#   args$x <- object
#   do.call(plot.elife_hazard, args = args)
# }


#' @export
#' @rdname plot.elife_northropcoleman
#' @param object object of class \code{elife_northropcoleman}, with the fitted piecewise-constant generalized Pareto model
autoplot.elife_northropcoleman <- function(object, ...){
  args <- list(...)
  args$plot.type <- "ggplot"
  args$x <- object
  do.call(plot.elife_northropcoleman, args = args)
}

#' @export
#' @rdname plot.elife_tstab
#' @param object object of class \code{elife_tstab}, representing parameter estimates to draw threshold stability plots
autoplot.elife_tstab <- function(object, ...){
  args <- list(...)
  args$plot.type <- "ggplot"
  args$x <- object
  do.call(plot.elife_tstab, args = args)
}

#' @export
#' @rdname plot.elife_profile
#' @param object object of class \code{elife_profile}
autoplot.elife_profile <- function(object, ...){
  args <- list(...)
  args$x <- object
  do.call(plot.elife_profile, args = args)
}



#' Check default arguments
#'
#' Check arguments and override default values.
#' If a named list, \code{arguments}, is provided by the user,
#' it will override any default value.
#' If one of the argument is provided directly,
#' it will take precedence over the values in \code{arguments}, provided it is not a default value.
#'
#' @param func function whose parameters are to be superseded
#' @param call user call, obtained from \code{match.call(expand.dots = FALSE)}
#' @param arguments named list with arguments
#' @keywords internal
#' @export
#' @return a named list with all arguments
check_arguments <- function(func,
                            call,
                            arguments = NULL){
  if(!is.null(arguments)){
    stopifnot("\"arguments\" should be a named list" = is.list(arguments))
  }
  arguments.default <- formals(func)
  new.arguments <- arguments.default
  # Override any default argument with the values in arguments
  new.arguments[names(arguments)] <- arguments
  # Function from https://stackoverflow.com/a/51389399
  to_list <- function(x) {
    xex <- parse(text = x )[[1]]
    xex[[1]] <- quote(list)
    eval.parent(xex)
  }
  supplied.arguments <- to_list(deparse(call[names(call) != "arguments"]))
  new.arguments[names(supplied.arguments)] <- supplied.arguments
  stopifnot("Mandatory \"time\" vector not supplied." = !is.null(new.arguments$time))
  new.arguments$arguments <- NULL
  return(new.arguments)
}

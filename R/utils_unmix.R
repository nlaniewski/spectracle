## from AutoSpectral::unmix.ols.fast()
solve_spectra <- function(spectra){
  ## should probably add a class to spectra
  unmixing.mat <- as.matrix(spectra[, .SD, .SDcols = is.numeric])
  ##
  rownames(unmixing.mat) <- spectra[
    ,
    j = gsub("NA ","", paste(S, N))
  ]
  ##
  XtX <- tcrossprod(unmixing.mat)
  ##
  unmixing.mat <- solve.default(XtX, unmixing.mat)
}
##
parameters.from.spectra <- function(spectra){
  spectra[
    ,
    j = .(
      N = paste0(as.character(N), '-A'),
      S = as.character(S),
      DETECTOR = as.character(detector),
      TYPE = "Unmixed_Fluorescence",
      B = "32",
      E = "0,0",
      R = "4194304",
      DISPLAY = "LOG"
    )
  ]
}
##
parameters.unmixed <- function(flowstate, spectra){
  ##
  unmixing.mat <- solve_spectra(spectra)
  ##
  parameters.unmixed <- parameters.from.spectra(spectra)
  ##
  par.n <- flowstate$parameters[,.N]
  ##
  parameters.unmixed[
    ,
    j = par := paste0('$P', (seq(.N) + par.n))
  ]
  ##
  parameters.unmixed <- merge(
    parameters.unmixed,
    flowstate$parameters[
      ,
      j = .(DETECTOR = sub("-A", "", N), V)
    ],
    sort = F
  )
  ##
  parameters.unmixed <- rbind(
    flowstate$parameters,
    parameters.unmixed,
    fill = TRUE
  )
  ##
  attributes.og <- names(attributes(flowstate$parameters))
  attributes.add <- attributes.og[
    !attributes.og %in% names(attributes(parameters.unmixed))
  ]
  if(length(attributes.add) > 0){
    for(attr.n in attributes.add){
      data.table::setattr(
        x = parameters.unmixed,
        name = attr.n,
        value = attr(flowstate$parameters, which = attr.n)
      )
    }
  }
  ##
  invisible(parameters.unmixed)
}
##


get.vars <- function(flowstate) {
  cols.by <- flowstate$data[, names(.SD), .SDcols = is.factor]
  detectors.pn <- flowstate$parameters[TYPE == "Raw_Fluorescence", N]
  cols.detector <- names(flowstate$data)[flowstate$data[, sapply(.SD, attr, "N") %in% detectors.pn]]
  cols.scatter <- grep("[FS]SC", names(flowstate$data), value = T)
  cols.mdat <- (
    names(flowstate$keywords)
    [names(flowstate$keywords) %in%
        c("$CYT", "$CYTSN", "CREATOR", "$PROJ", "$DATE")
    ]
  )
  mdat <- flowstate$keywords[, unique(.SD), .SDcols = cols.mdat]
  ##
  list(
    cols.by = cols.by,
    detectors = cols.detector,
    scatter = cols.scatter,
    mdat = mdat
  )
}









remove.saturating.doublets <- function(flowstate){
  ## remove saturating events
  flowstate::select_nonsaturating(flowstate)
  flowstate <- subset(flowstate, select.nonsaturating)
  flowstate$data[, select.nonsaturating := NULL]
  ## remove doublet events
  flowstate::select_singlets(flowstate, quantiles = c(0.85, 0.975))#@params
  flowstate <- subset(flowstate, select.singlets)
  flowstate$data[, select.singlets := NULL]
  ##
  invisible(flowstate)
}

mtext.keywords <- c("$DATE", "$CYT", "$CYTSN", "CREATOR", "$PROJ")

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
  mdat <- flowstate$keywords[, .SD, .SDcols = c('sample.id', cols.mdat)]
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
  data.table::setindex(flowstate$parameters, NULL)
  ## remove doublet events
  flowstate::select_singlets(flowstate, quantiles = c(0.85, 0.975))#@params
  flowstate <- subset(flowstate, select.singlets)
  flowstate$data[, select.singlets := NULL]
  # data.table::setindex(flowstate$parameters, NULL)
  ##
  invisible(flowstate)
}

.af.signatures <- function(flowstate){
  vars <- get.vars(flowstate)
  ## autofluorescence -- generalized (for now);
  ## mean and median vectors (mean)
  flowstate$data[
    i = N == "AF",
    j = {
      .mean <- lapply(.SD, mean)
      .median <- lapply(.SD, stats::median)
      ##
      data.table::rbindlist(
        list(
          c(detector = names(which.max(.mean)), .mean, vector.type = "mean"),
          c(detector = names(which.max(.median)), .median, vector.type = "median")
        )
      )
    },
    .SDcols = vars$detectors,
    by = c(vars$cols.by)
  ][, detector := factor(detector, levels = vars$detectors)][]
}

## from pals::kelly(22)
colors.kelly <- c(
  #"#F2F3F4"
  "#222222", "#F3C300", "#875692", "#F38400", "#A1CAF1",
  "#BE0032", "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5",
  "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300", "#882D17",
  "#8DB600", "#654522", "#E25822", "#2B3D26"
)

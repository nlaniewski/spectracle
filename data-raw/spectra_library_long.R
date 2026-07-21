ref.library <- readRDS(list.files(
  system.file("extdata/spectra_library", package = "spectracle"),
  full.names = T,
  pattern = "Aurora.*wide"
))
data.table::setnames(ref.library, 'detector.peak', 'detector')
ref.library[, laser := factor(laser, levels = unique(laser))]
ref.library[
  ,
  j = detector := {
    res <- sub("\\).*$", "", sub("^.*\\(", "", unique(configuration)))
    res <- trimws(strsplit(res, ",")[[1]])
    detector.levels <- unlist(lapply(res, function(j){
      laser <- gsub("\\d", "", j)
      detector.n <- as.numeric(gsub("\\D", "", j))
      paste0(laser, seq(detector.n))
    }))
    factor(detector, levels = detector.levels)
  }
]
data.table::setorder(ref.library, laser, detector)

traces <- data.table::melt(
  data = ref.library,
  measure.vars = ref.library[, names(.SD), .SDcols = is.numeric],
  variable.name = "Detector"
)

#' @title Derive Spectra from Full Spectrum Cytometry Controls
#' @description
#' Manufacturer::Cytometer::Software tested:
#' \itemize{
#'   \item Cytek::Aurora(5L)::SpectroFlo 3.3.0
#' }
#'
#' As a fully data-driven/automated process -- for best performance -- `spectracle` relies on adherence to established naming conventions for appropriate identification of the unstained control, parsing of keyword-value pairs into factored metadata, and accurate annotation of derived spectra.
#'
#' @param raw.reference.controls a character vector of full directory/path names to raw single-color control .fcs files; one of [list.dirs] or [list.files].
#' @param top.expressing.override a named numeric vector -- default `NULL`; to override the automated selection of 'spectral.events', the name(s) of the supplied numeric vector should match to unique .fcs sample name(s). Example as follows: c("BUV615 (Cells)" = 20); for 'BUV615', only 20 top expressing events will be used to characterize/derive the final spectra (see example below).
#' @param print.timings logical -- default `FALSE`; if `TRUE`, a [tictoc][tictoc::tic.log] will print to the console, detailing execution times.
#' @param ... not defined; placeholder.
#'
#' @returns
#' A [data.table][data.table::data.table] of normalized `[0,1]` spectra; an additional *generalized* autofluoresence signature is appended to the table.
#' @export
#'
#' @examples
#' \dontrun{
#' raw.reference.controls <- fcs.files <- list.files(
#' system.file("extdata/expt1", package = "spectracle"),
#' full.names = TRUE
#' )
#'
#' spectra <- spectracle(raw.reference.controls)
#'
#' ## name fix due to no marker '$PnS' names in these source files
#' ## add empirically determined marker names
#' spectra[
#' i = N == '450',
#' j = c('N', 'S') := list('eFluor 450', 'TCRgd')
#' ]
#' spectra[
#' i = N == 'BUV615',
#' j = c('N', 'S') := list('BUV615', 'Siglec F')
#' ]
#'
#' ##
#' plot_trace(spectra[N == "eFluor 450"])
#'
#' ## this particular control is low-event count/rare
#' ## impacted by AF -- the spectra is off...
#' plot_trace(spectra[N == "BUV615"])
#'
#' ## rerun with a small tweak
#' spectra <- spectracle(
#' raw.reference.controls[c(1,3)],
#' top.expressing.override = c("BUV615 (Cells)" = 20)
#' )
#'
#' ## spectra is derived
#' plot_trace(spectra[N == "BUV615"])
#' }
#'
spectracle <- function(
    raw.reference.controls,
    top.expressing.override = NULL,
    print.timings= FALSE,
    ...
)
{
  ## defensive coding can go here to check for raw, conformity, software, parameters/vars, etc.
  ## ...

  ## accepts a directory (containing .fcs files) or accepts .fcs file paths
  l <- length(raw.reference.controls)
  if(!all(grepl(".fcs", raw.reference.controls)) & l == 1){
    ## get paths to raw .fcs files if a directory;
    ## reference group controls
    ref.paths <- list.files(raw.reference.controls, full.names = T)
  }else if(all(grepl(".fcs", raw.reference.controls)) & l != 1){
    ref.paths <- raw.reference.controls
  }

  ## reference group -- raw, single-color controls;
  ## cells only for now

  tictoc::tic.clear()
  tictoc::tic.clearlog()
  tictoc::tic("read/concatenate flowstate")

  suppressMessages(
    ref <- flowstate::read.flowstate(
      fcs.file.paths = ref.paths,
      colnames.type = 'N',
      concatenate = T
    )
  )

  tictoc::toc(log = TRUE, quiet = TRUE)
  ##

  tictoc::tic("preprocess -- add keywords, remove: saturating events; doublets")

  ## software naming convention into keywords;
  ## Spectroflo only for now
  flowstate:::reference.group.keywords(ref)
  ## remove saturating events
  flowstate::select_nonsaturating(ref)
  ref <- subset(ref, select.nonsaturating)
  ref$data[, select.nonsaturating := NULL]
  ## remove doublet events
  flowstate::select_singlets(ref, quantiles = c(0.85, 0.975))#@params
  ref <- subset(ref, select.singlets)
  ref$data[, select.singlets := NULL]

  tictoc::toc(log = TRUE, quiet = TRUE)

  ## variables needed
  cols.by <- ref$data[, names(.SD), .SDcols = is.factor]
  detectors.pn <- ref$parameters[TYPE == "Raw_Fluorescence", N]
  cols.detector <- names(ref$data)[ref$data[, sapply(.SD, attr, 'N') %in% detectors.pn]]
  cols.scatter <- grep("[FS]SC", names(ref$data), value = T)
  cols.mdat <- names(ref$keywords)[names(ref$keywords) %in% c('$CYT', '$CYTSN', 'CREATOR', '$PROJ', '$DATE')]
  mdat <- ref$keywords[, unique(.SD), .SDcols = cols.mdat]

  tictoc::tic("autofluorescence -- characterizing")

  ## autofluorescence -- generalized; vector (mean)
  af.vec.mean <- ref$data[
    i = N == "AF" & tissue.type == "Cells",
    j = sapply(.SD, mean),
    .SDcols = cols.detector
  ]
  ## autofluorescence -- generalized; vector (median)
  ## for adding to spectra during final output
  af.vec.median <- ref$data[
    i = N == "AF" & tissue.type == "Cells",
    j = lapply(.SD, stats::median),
    .SDcols = cols.detector,
    by = cols.by
  ][, detector := names(which.max(.SD)), .SDcols = cols.detector]
  ## add dominant AF detector to [['data']]
  ## UPDATES BY REFERENCE
  ref$data[
    i = N == "AF" & tissue.type == "Cells",
    j = detector := names(which.max(af.vec.mean))
  ]

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::tic("removing AF vector (projection-based orthogonalization) and deriving peak detector")

  ## METHOD -- projection-based orthogonalization (googled this term for R code...)
  ## normalize the (nuisance) vector
  v_unit <- af.vec.mean / sqrt(sum(af.vec.mean^2))
  v_unit_t <- t(v_unit)
  ## apply METHOD to [['data']];
  ## UPDATES BY REFERENCE
  ref$data[
    i = N != 'AF',
    j = detector := {
      ## the internal bits are vectorized here (I think...);
      ## possible speed gain with another projection method; lm/residuals?
      ## matrix operation -- .SD as matrix
      mat <- as.matrix(.SD)
      ## calculate the projection of each row onto vector and subtract it
      ## matrix multiplication -- remove AF/nuisance vector
      mat_cleaned <- mat - (mat %*% v_unit) %*% v_unit_t
      ## peak detector derived from 'mat_cleaned'
      ## UPDATES BY REFERENCE in [['data']]
      detector.peak <- names(which.max(apply(mat_cleaned, 2, mean)))
      ##
    },
    .SDcols = cols.detector,
    by = cols.by
  ][, detector := factor(detector)]
  ## update 'cols.by'
  cols.by <- ref$data[, names(.SD), .SDcols = is.factor]

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::tic("ordering peak detector vectors to select top expressing ('spectral') events")

  ## add 'spectral.events' logical to [['data']]
  ## UPDATES BY REFERENCE
  ref$data[, spectral.events := FALSE]
  ref$data[
    i = N == "AF" & tissue.type == "Cells",
    j = spectral.events := TRUE
  ]
  ##
  ref$data[
    i = N != "AF",
    j = spectral.events := {
      ## basic vector sorting for top (peak detector) expressing events
      detector <- as.character(.BY$detector)
      ## index the top (peak detector) expressing events -- 1000
      i.top <- order(.SD[[detector]], decreasing = T)[1:1000]
      ## cosine similarity against af.vector
      i.cs <- cosine.similarity.mat(
        x = as.matrix(.SD[i.top]),
        af.vec.median[, unlist(.SD), .SDcols = is.numeric]
      )
      ## set top n expressing events
      n <- 200
      if(!is.null(top.expressing.override)){
        id <- as.character(.BY$sample.id)
        if(names(top.expressing.override) %in% id){
          n <- top.expressing.override[[id]]
        }
      }
      ## re-index -- lowest cosine score
      i.top <- i.top[order(i.cs)][1:n]
      ## index and set the logical
      spectral.events[i.top] <- TRUE
      ## return the logical
      spectral.events
    },
    .SDcols = cols.detector,
    by = cols.by
  ]
  ## subset ref to retain only spectral events
  ref.spectral <- subset(ref, spectral.events)
  ## CONDITIONAL: if 'ref' is not needed (no return), rm and gc
  rm(ref) ; invisible(gc())
  ## update peak detector; some fluors/peak detectors are fully resolved at this point
  ref.spectral$data[
    i = N != "AF" & tissue.type == "Cells",
    j = detector := {
      names(which.max(lapply(.SD, mean)))
    },
    .SDcols = cols.detector,
    by = cols.by
  ]
  ref.spectral$data[, detector := factor(detector, levels = cols.detector)]

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::tic("AF scatter matching using nearest neighbors")

  ## scatter match to subtract AF; per-nearest two cells/neighbors
  ## input/data matrix for FNN::knnx.index()
  ## subset AF ; [['data']] only; 'cols.scatter' only
  af.scatter <- (
    subset(ref.spectral, N == "AF" & tissue.type == 'Cells')
    [['data']][, .SD, .SDcols = cols.scatter]
  )
  ## subset AF ; [['data']] only; 'cols.detector' only
  af.detectors <- (
    subset(ref.spectral, N == "AF" & tissue.type == 'Cells')
    [['data']][, .SD, .SDcols = cols.detector]
  )
  ## some data.table abuse below...
  ## trying to capture the derivatives/intermediates so that they can be returned if needed
  ## ...
  k <- 2#@params
  cols.nni <- paste0('nni.', seq(k))
  ## 1) scatter match to 'af.scatter' to create an index 'nni'
  ref.spectral$data[
    i = N != "AF",
    j = (cols.nni) := {
      ## 1) nearest neighbors scatter match between spectral events and AF/unstained
      nni <- FNN::knnx.index(
        data = af.scatter,
        query = .SD,
        k = k
      )
      ## UPDATE BY REFERENCE the indices that match to 'af.scatter'
      as.data.frame(nni)
    },
    .SDcols = cols.scatter,
    by = cols.by
  ]

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::tic("AF medians matching")

  ## 2) median for each group/row of points in 'nni' to derive scatter-matched 'af.detectors' vectors
  ## this is the slowest part, at least according to my timings -- there is a faster way, I'm sure;
  ## matrix operations will probably be faster here
  af.medians <- ref.spectral$data[
    i = N != "AF",
    j = {
      ## this seems quite a hack to match up the points...but it works
      ## data.table optimized loop using data.table::set to index the row groupings in 'af.detectors'
      for(i in seq(.N)){
        data.table::set(
          x = af.detectors,
          i = c(t(.SD[i])),
          j = 'group',
          value = i
        )
      }
      ## medians by group
      af.medians <- af.detectors[
        i = !is.na(group),
        j = lapply(.SD, stats::median),
        keyby = group
      ]
      ## NULL out the temporary groupings
      af.detectors[, group := NULL]
      ## return the medians
      af.medians
    },
    .SDcols = cols.nni,
    by = cols.by
  ]

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::tic("AF medians subtraction and deriving of spectra")

  ## 3) subtract out the matched AF vectors (rows); medians for linear vector; normalize
  spectra <- ref.spectral$data[
    i = N != "AF",
    j = {
      ## 3 CONDITIONAL) af.medians (matched to sample) can have fewer rows due to duplicate nni matches;
      ## get the group index to drop rows in .SD so they match up for subtraction
      i.group <- af.medians[i = sample.id == .BY$sample.id, group]
      spectral.vec.linear <- sapply(
        (.SD[i.group] - af.medians[
          i = sample.id == .BY$sample.id,
          j = .SD,
          .SDcols = cols.detector
        ]),
        stats::median
      )
      spectral.vec.norm <- spectral.vec.linear / max(spectral.vec.linear)
      spectral.vec.norm[spectral.vec.norm < 0] <- 0
      as.list(spectral.vec.norm)
    },
    .SDcols = cols.detector,
    by = cols.by
  ]

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::tic("Preparing final output -- normalized spectra")

  ## prepare spectra for final output
  ## exhaustive metadata can/should be added for reproducibility/transparency
  spectra <- rbind(spectra, af.vec.median)
  spectra[
    ,
    laser := factor(
      gsub("\\d+", "", detector),
      levels = unique(gsub("\\d+", "", cols.detector))
    )
  ]
  data.table::setcolorder(spectra, c(cols.by, 'laser'))
  spectra <- cbind(deriving.function = 'spectracle', mdat, spectra)
  data.table::setorder(spectra, laser, detector)
  ##AF to last position
  spectra[N != "AF", ord := seq(.N)]
  spectra[N == "AF", ord := spectra[, .N]]
  data.table::setorder(spectra, ord)[, ord := NULL]
  ##

  tictoc::toc(log = TRUE, quiet = TRUE)
  log <- tictoc::tic.log(format = TRUE)
  if(print.timings){
    print(unlist(log))
  }
  ##
  invisible(spectra)
}

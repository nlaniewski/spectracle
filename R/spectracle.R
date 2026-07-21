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
#' @param filter.top.expressing a character vector -- default `NULL`; to modulate the automated selection of 'spectral.events', the supplied character vector should match to unique .fcs sample name(s). If defined, a conditional filtering will take place to return only the 0.99 (cosine) most similar 'spectral.events'.
#' @param top.expressing.override a named numeric vector -- default `NULL`; to override the automated selection of 'spectral.events', the name(s) of the supplied numeric vector should match to unique .fcs sample name(s). Example as follows: c("BUV615 (Cells)" = 20); for 'BUV615', only 20 top expressing events will be used to characterize/derive the final spectra (see example below).
#' @param print.timings logical -- default `FALSE`; if `TRUE`, a [tictoc][tictoc::tic.log] will print to the console, detailing execution times.
#' @param unstained.ids a named character vector -- default `NULL`; in the odd case where the unstained sample(s) are ambiguously named, `unstained.ids` can be defined.  The named character vector should take the following form (representative examples of ambiguously named samples):
#' \itemize{
#'   \item `c("Unstained (Beads)" = "BioLegend Positive Beads (Beads)")`
#'   \item `c("Unstained (Cells)" = "PBMC (Cells)")`
#' }
#' @param ... not defined; placeholder.
#'
#' @returns
#' A [data.table][data.table::data.table] of normalized `[0,1]` spectra; an additional *generalized* autofluoresence signature is appended to the table.
#' @export
#'
#' @examples
#'
#' ## due to large file sizes, the .fcs files used in this example were first
#' ## fully processed -- spectra was derived and saved; load here
#'
#' spectra.files <- list.files(
#' system.file("extdata/prepared_spectra", package = "spectracle"),
#' full.names = TRUE
#' )
#' names(spectra.files) <- sub(".rds", "", basename(spectra.files))
#'
#' spectra.examples <- lapply(spectra.files, readRDS)
#'
#' ## expt1
#'
#' spectra <- spectra.examples$spectra_expt1_nofilter
#' plot_trace(spectra[N == "eFluor 450"])
#'
#' ## this particular control is low-event count/rare
#' ## impacted by AF -- the spectra is off...
#' plot_trace(spectra[N == "BUV615"])
#'
#' ## ...rerun with a small tweak using:
#' ## `filter.top.expressing = c("BUV615 (Cells)")`
#' spectra <- spectra.examples$spectra_expt1_filtered
#'
#' ## spectra is derived
#' plot_trace(spectra[N == "BUV615"])
#'
#'
#' ## expt2 -- beads and cells
#' spectra <- spectra <- spectra.examples$spectra_expt2
#'
#' plot_trace(spectra)
#'
#' spectra[, .(tissue.type, N, hash.md5)]
#'
spectracle <- function(
    raw.reference.controls,
    unstained.ids = NULL,
    filter.top.expressing = NULL,
    top.expressing.override = NULL,
    print.timings= FALSE,
    ...
)
{
  ## defensive coding can go here to check for raw, conformity, software, parameters/vars, etc.
  ## ...

  ## accepts a directory (containing .fcs files) or accepts .fcs file paths
  l <- length(raw.reference.controls)
  if(!all(grepl(".fcs", raw.reference.controls)) && l == 1){
    ## single directory path: enumerate .fcs files within it
    ref.paths <- list.files(
      raw.reference.controls,
      full.names = TRUE,
      pattern = ".fcs"
    )
  }else if(all(grepl(".fcs", raw.reference.controls)) && l != 1){
    ref.paths <- raw.reference.controls
  }

  ## reference group -- raw, single-color controls;
  ## cells only for now

  tictoc::tic.clear()
  tictoc::tic.clearlog()

  tictoc::tic("Making a Spectacle of Your Spectra with spectracle")

  tictoc::tic("read/concatenate flowstate")

  suppressMessages(
    ref <- flowstate::read.flowstate(
      fcs.file.paths = ref.paths,
      colnames.type = 'N',
      concatenate = T
    )
  )
  ## name-fix
  for(id in names(unstained.ids)){
    for(ii in c('keywords', 'data')){
      ref[[ii]][
        i = grep(unstained.ids[[id]], sample.id, fixed = T),
        j = sample.id := id
      ]
    }
  }

  tictoc::toc(log = TRUE, quiet = TRUE)
  ##

  tictoc::tic("preprocess -- add keywords, remove: saturating events; doublets")

  ## software naming convention into keywords;
  ## Spectroflo only for now
  reference.group.keywords(ref)
  ## variables needed
  vars <- get.vars(ref)
  ## remove saturating events
  flowstate::select_nonsaturating(ref)
  ref <- subset(ref, select.nonsaturating)
  ref$data[, select.nonsaturating := NULL]
  ## remove doublet events
  flowstate::select_singlets(ref, quantiles = c(0.85, 0.975))#@params
  ref <- subset(ref, select.singlets)
  ref$data[, select.singlets := NULL]
  ## bead-specific preprocessing
  if("Beads" %in% ref$data[, levels(tissue.type)]){
    ref$data[
      i = tissue.type != "Beads",
      j = select.beads := TRUE
    ]
    ref$data[
      i = tissue.type == "Beads",
      j = select.beads := {
        bounds.fsc <- peak.bounds(FSC_A, height.threshold = 0.1)
        bounds.ssc <- peak.bounds(SSC_A, height.threshold = 0.1)
        data.table::`%between%`(FSC_A, bounds.fsc) & data.table::`%between%`(SSC_A, bounds.ssc)
      },
      by = c(vars$cols.by)
    ]
    ref <- subset(ref, select.beads)
    ref$data[, select.beads := NULL]
  }

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::tic("autofluorescence -- characterizing")

  ## autofluorescence -- generalized (for now);
  ## mean and median vectors (mean)
  af.signatures <- .af.signatures(ref)

  ## add dominant AF detector to [['data']]
  ## UPDATES BY REFERENCE
  ref$data[
    i = N == "AF",
    j = detector := af.signatures[
      i = tissue.type == .BY$tissue.type & vector.type == 'mean'][['detector']],
    by = c(vars$cols.by)
  ]

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::tic("removing AF vector (projection-based orthogonalization) and deriving peak detector")

  ### START peak detector

  ## for 'Cells'
  if('Cells' %in% ref$data[, levels(tissue.type)]){
    ## METHOD -- projection-based orthogonalization (googled this term for R code...)
    ## normalize the (intrusive/nuisance) vector
    v <- af.signatures[
      i = tissue.type == "Cells" & vector.type == 'mean',
      j = unlist(.SD),
      .SDcols = is.numeric
    ]
    v_unit <- v / sqrt(sum(v^2))
    v_unit_t <- t(v_unit)
  }
  ## apply METHOD to [['data']];
  ## UPDATES BY REFERENCE
  ref$data[
    i = N != 'AF',# & tissue.type == "Cells",
    j = detector := {
      if(.BY$tissue.type == 'Cells'){
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
      }else if(.BY$tissue.type == 'Beads'){
        names(which.max(sapply(.SD, mean)))
      }
    },
    .SDcols = vars$detectors,
    by = c(vars$cols.by)
  ]

  ## END peak detector

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::tic("ordering peak detector vectors to select top expressing ('spectral') events")

  ## START 'spectral.events'

  ## add 'spectral.events' logical to [['data']]
  ## UPDATES BY REFERENCE
  ref$data[, spectral.events := FALSE]
  ## set 'AF' as TRUE
  ref$data[
    i = N == "AF",
    j = spectral.events := TRUE
  ]
  ## for 'Cells'
  if('Cells' %in% ref$data[, levels(tissue.type)]){
    ## 'AF' median vector for inside the data.table loop -- 'Cells'
    v.median <- af.signatures[
      i = tissue.type == "Cells" & vector.type == 'median',
      j = unlist(.SD),
      .SDcols = is.numeric
    ]
  }
  ##
  ref$data[
    i = N != "AF",# & tissue.type == "Cells",
    j = spectral.events := {
      if(.BY$tissue.type == 'Cells'){
        ## basic vector sorting for top (peak detector) expressing events
        detector <- as.character(.BY$detector)
        ## index the top (peak detector) expressing events -- 1000
        i.top <- order(.SD[[detector]], decreasing = T)[1:1000]
        ## cosine similarity against af.vector
        i.cs <- cosine.similarity.mat(
          x = as.matrix(.SD[i.top]),
          reference.vector = v.median
        )
        ## set top n expressing events
        n <- 200
        if (!is.null(top.expressing.override)) {
          id <- as.character(.BY$sample.id)
          if (id %in% names(top.expressing.override)) {
            n <- top.expressing.override[[id]]
          }
        }
        ## re-index -- lowest cosine score
        i.top <- i.top[order(i.cs)][1:n]

        if (!is.null(filter.top.expressing)) {
          id <- as.character(.BY$sample.id)
          if (id %in% filter.top.expressing) {
            ## a final (final?) filter
            i.max <- which.max(.SD[i.top][[detector]])
            i.cs <- cosine.similarity.mat(
              x = as.matrix(.SD[i.top]),
              reference.vector = unlist(.SD[i.top][i.max])
            )
            i.top <- i.top[which(i.cs > 0.99)]
          }
        }

        # n <- length(i.top)
        # r <- range(.SD[i.top])
        # colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(n)
        # colors <- rev(colors)
        # for(i in seq(n)){
        #   plot_spectral.trace.base(
        #     unlist(.SD[i.top[i]]),
        #     ylim = r,
        #     add.lines = isTRUE(i != 1),
        #     col = colors[i],
        #     sub = sprintf("Detector (peak): %s", detector)
        #     # ...
        #   )
        # }

      }else if(.BY$tissue.type == 'Beads'){
        ## basic vector sorting for top (peak detector) expressing events
        detector <- as.character(.BY$detector)
        ## index the top (peak detector) expressing events -- 200
        i.top <- order(.SD[[detector]], decreasing = T)[1:200]
      }
      ## index and set the logical
      spectral.events[i.top] <- TRUE
      ## return the logical
      spectral.events
    },
    .SDcols = vars$detectors,
    by = c(vars$cols.by, 'detector')
  ]

  ## subset ref to retain only spectral events
  ref.spectral <- subset(ref, spectral.events)
  ref.spectral$data[, spectral.events := NULL]
  ## CONDITIONAL: if 'ref' is not needed (no return), rm and gc
  rm(ref) ; invisible(gc())

  ## END 'spectral.events'

  ## update peak detector; some fluors/peak detectors are fully resolved at this point -- 'Cells'
  ref.spectral$data[
    i = N != "AF" & tissue.type == "Cells",
    j = detector := {
      names(which.max(lapply(.SD, mean)))
    },
    .SDcols = vars$detectors,
    by = c(vars$cols.by, 'detector')
  ]

  tictoc::toc(log = TRUE, quiet = TRUE)

  ## for 'Cells' only
  if('Cells' %in% ref.spectral$data[, levels(tissue.type)]){

    tictoc::tic("AF scatter matching using nearest neighbors")

    ## prepare AF subsets for eventual scatter-matching/subtraction
    af.sub       <- subset(ref.spectral, N == "AF" & tissue.type == "Cells")[["data"]]
    af.scatter   <- af.sub[, .SD, .SDcols = vars$scatter]
    af.detectors <- af.sub[, .SD, .SDcols = vars$detectors]
    rm(af.sub)

    k        <- 2L  #@params
    cols.nni <- paste0("nni.", seq_len(k))

    ## single-pass KNN across all non-AF events
    all.spectral.scatter <- ref.spectral$data[
      i = N != "AF" & tissue.type == 'Cells',
      j = .SD,
      .SDcols = vars$scatter
    ]# this might create a copy in memory; need to check mem address
    all.nni <- FNN::knnx.index(
      data  = as.matrix(af.scatter),
      query = as.matrix(all.spectral.scatter),
      k     = k
    )

    ref.spectral$data[
      i = N != "AF" & tissue.type == 'Cells',
      (cols.nni) := as.data.frame(all.nni)
    ]
    rm(all.spectral.scatter, all.nni)

    tictoc::toc(log = TRUE, quiet = TRUE)

    tictoc::tic("AF medians matching")

    ## 2) median for each group/row of points in 'nni' to derive scatter-matched 'af.detectors' vectors
    ## replaced with vectorised matrix indexing -- for k=2 the median of two values equals their mean
    af.det.mat <- as.matrix(af.detectors[, .SD, .SDcols = vars$detectors])
    n.det    <- ncol(af.det.mat)
    af.medians <- ref.spectral$data[
      i = N != "AF" & tissue.type == 'Cells',
      j = {
        nni.mat  <- as.matrix(.SD)# .N x k: each row = neighbor indices
        af.med   <- matrix(0, .N, n.det)
        for (ki in seq_len(k)) {
          af.med <- af.med + af.det.mat[nni.mat[, ki], , drop = FALSE]
        }
        af.med <- af.med / k# mean for k = 2; k is fixed for now
        colnames(af.med) <- vars$detectors
        cbind(
          data.table::data.table(group = seq_len(.N)),
          data.table::as.data.table(af.med))
      },
      .SDcols = cols.nni,
      by = c(vars$cols.by, 'detector')
    ]

    tictoc::toc(log = TRUE, quiet = TRUE)
  }

  tictoc::tic("AF medians subtraction and deriving of spectra")

  ## START AF/background subtraction; normalized spectra

  ## 3) subtract out the matched AF vectors (rows); medians for linear vector; normalize
  spectra <- ref.spectral$data[
    ,#i = N != "AF", #& tissue.type == 'Cells',
    j = {
      if(.BY$N != "AF"){
        if(.BY$tissue.type == 'Cells'){
          ## 3 CONDITIONAL) af.medians (matched to sample) can have fewer rows due to duplicate nni matches;
          ## get the group index to drop rows in .SD so they match up for subtraction
          i.group <- af.medians[i = sample.id == .BY$sample.id, group]
          spectral.vec.linear <- sapply(
            (.SD[i.group] - af.medians[
              i = sample.id == .BY$sample.id,
              j = .SD,
              .SDcols = vars$detectors
            ]),
            stats::median
          )
        }else if(.BY$tissue.type == "Beads"){
          vec <- af.signatures[
            i = tissue.type == "Beads" & vector.type == 'median',
            j = unlist(.SD),
            .SDcols = vars$detectors
          ]
          ##
          spectral.vec.linear <- sapply((.SD - vec), stats::median)
        }
        spectral.vec.norm <- spectral.vec.linear / max(spectral.vec.linear)
        spectral.vec.norm[spectral.vec.norm < 0] <- 0
        as.list(spectral.vec.norm)
      }else if(.BY$N == "AF"){
        spectral.vec.linear <- sapply(.SD, stats::median)
        spectral.vec.norm <- spectral.vec.linear / max(spectral.vec.linear)
        as.list(spectral.vec.norm)
      }
    },
    .SDcols = vars$detectors,
    by = c(vars$cols.by, 'detector')
  ]

  ## END AF/background subtraction; normalized spectra

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::tic("Preparing final output -- normalized spectra")

  ## prepare spectra for final output
  ## exhaustive metadata can/should be added for reproducibility/transparency
  # spectra <- rbind(
  #   spectra,
  #   af.signatures[
  #     i = tissue.type %in% spectra[, levels(tissue.type)] & vector.type == 'median',
  #     j = !'vector.type'
  #   ]
  # )
  spectra[
    ,
    laser := factor(
      gsub("\\d+", "", detector),
      levels = unique(gsub("\\d+", "", vars$detectors))
    )
  ]
  data.table::setcolorder(spectra, c(vars$cols.by, 'detector', 'laser'))
  ##
  pkg <- "spectracle"
  deriving.function <- sprintf("%s_%s", pkg, utils::packageVersion(pkg))
  spectra <- cbind(deriving.function, vars$mdat[spectra, on = 'sample.id'])
  ##
  data.table::setorder(spectra, laser, detector)
  ## AF to last position
  spectra[N != "AF", ord := seq(.N)]
  spectra[N == "AF", ord := max(spectra[['ord']], na.rm = T) + seq(.N)]
  data.table::setorder(spectra, ord)[, ord := NULL]
  ## add hashes
  spectra[
    ,
    hash.md5 := apply(.SD, 1, digest::digest, algo = "md5"),
    .SDcols = vars$detectors
  ]

  tictoc::toc(log = TRUE, quiet = TRUE)

  tictoc::toc(log = TRUE, quiet = TRUE)

  log <- tictoc::tic.log(format = TRUE)

  if(print.timings){
    print(unlist(log))
  }
  ##
  invisible(spectra)
}

## helper function: cosine similarity;
## scored against a reference vector
cosine.similarity.mat <- function(x, reference.vector){
  ## conformable?
  if(!all(names(reference.vector) %in% colnames(x))){
    stop('non-conformable')
  }
  ##
  .prod <- x %*% reference.vector
  .mat <- sqrt(rowSums(x^2))
  .vec <- sqrt(sum(reference.vector^2))
  ##
  cs <- (.prod / (.mat * .vec))
  ##
  return(cs)
}

## helper function: plot a spectral trace
plot_spectral.trace.base <- function(x, add.lines = F, ...){
  ## plot spectral trace
  if(!add.lines){
    plot(
      x,
      type = 'l',
      xlab = "Detectors",
      ylab = ifelse(max(x) == 1, "Expression [0,1]", "Expression"),
      xaxt = "n",
      ...
    )
    graphics::axis(
      1,
      at = seq(length(x)),
      labels = sub("-A$", "", names(x)),
      las = 2,
      cex.axis = 0.75
    )
  }else{
    graphics::lines(x, ...)
  }
}

#' @title Plot a Spectral Trace
#'
#' @param spectra [data.table][data.table::data.table] -- the return of [spectracle].
#'
#' @returns
#' A plot is printed to the active device.
#' @export
#'
plot_trace <- function(spectra){
  ##
  spectra[
    ,
    j = {
      plot_spectral.trace.base(
        x = unlist(.SD),
        main = ifelse(
          .BY$N == "AF",
          paste0(.BY$detector, '::', .BY$N, '::', .BY$tissue.type),
          paste0(.BY$detector, '::', .BY$N, '::', .BY$S, '::', .BY$tissue.type)
        ),
        col = "purple"
      )
      graphics::abline(v = which(levels(detector) %in% detector), lty = 'dashed')
      graphics::mtext(paste(vars$mdat, collapse = "::"), side = 1, line = 4, adj = 0, cex = 0.75)
    },
    .SDcols = is.numeric,
    by = c(spectra[, names(.SD), .SDcols = is.factor])
  ]
  ##
  invisible()
}

## helper function: quantile filtering of cosine similarity scores
## data-driven/programmatic method of mitigating the influence of AF;
## isolates peak detector/representative spectral vectors
## plot to see visualize representation of the effect of filtering
cosine.similarity.filter <- function(mat, method = c('median', 'mean'), plot = F, ...){
  ## stop/error message
  stopifnot(
    "mat needs a 'cosine.sim' column -- not found." = 'cosine.sim' %in% colnames(mat)
  )
  cols <- grep('cosine.sim', colnames(mat), value = T, invert = T)
  ## bins (quantile probs) for grouping cosine similarity ('cosine.sim') values
  bins <- seq(0.1, 1, length.out = 19)
  bins <- c(seq(0, 0.01, length.out = 9), bins/10, bins)
  bins <- sort(unique(bins))
  ## quantile breaks and groups from bins
  breaks <- stats::quantile(mat[, 'cosine.sim'], probs = bins)
  groups <- as.integer(cut(mat[, 'cosine.sim'], breaks = breaks, include.lowest = TRUE))
  ## which statistic to summarize groups
  res <- switch(
    match.arg(method),
    median = collapse::fmedian(mat[, cols], g = groups),
    mean = collapse::fmean(mat[, cols], g = groups)
  )

  ## !!! conditionals to handle edge-cases
  ## conditional for if a detector appears only once (artifact/unstable?) -- see eFlour 450::TCRgd (V11)
  if(any(table(max.col(res)) == 1)){
    drop.i <- as.numeric(names(which.min(table(max.col(res)))))
    ## if the detector appears only once but is in the first bin (least AF-like), retain -- see BUV605::Siglec F
    if(which(max.col(res) == drop.i) != 1){
      res <- res[!max.col(res) == drop.i,]
    }
  }
  ## !!!

  ## index the max for the whole res matrix
  i.max <- which(res == max(res), arr.ind = TRUE)
  ## peak detector
  detector.peak <- colnames(res)[i.max[, 2]]

  ## visualize
  if(plot){
    r <- range(res)
    colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(nrow(res))
    colors <- rev(colors)
    for(i in seq(nrow(res))){
      plot_spectral.trace.base(
        res[i, ],
        ylim = r,
        add.lines = isTRUE(i != 1),
        col = colors[i],
        sub = sprintf("Detector (peak): %s", detector.peak)
        # ...
      )
    }
    graphics::abline(v = which(colnames(res) %in% detector.peak), lty = "dashed")
  }
  ##
  return(
    list(
      csfm = res,
      detector.peak = detector.peak
    )
  )
}

cosine.filter.spectral.events <- function(x, threshold = 0.999){
  ##
  mat <- as.matrix(x)
  ## Normalize each row to unit length (Euclidean norm = 1)
  X_norm <- x / sqrt(rowSums(x^2))
  ## Matrix multiplication -- cosine similarity
  sim_matrix <- X_norm %*% t(X_norm)
  ## logical condition to identify matching pairs; affected by threshold
  threshold <- 0.999
  matching_rows <- which(rowSums(sim_matrix >= threshold & sim_matrix < 1) > 0)
  ## unique indicies
  unique(matching_rows)
}

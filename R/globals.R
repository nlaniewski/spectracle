.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
# dput(sort(c(...)))
# c(...) |> sort() |> dput()
utils::globalVariables(
  c(
    c("N", "TYPE", "detector", "group", "laser", "ord", "sample.id",
      "select.nonsaturating", "select.singlets", "tissue.type"),
    c(".mtext", "CREATOR", "hash.md5", "vars"),
    c('vector.type'),
    c('FSC_A', 'SSC_A', 'select.beads'),
    c('$CYT'),
    c('alias', 'S')
  )
)

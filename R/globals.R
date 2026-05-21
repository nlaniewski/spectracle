.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
# dput(sort(c(...)))
# c(...) |> sort() |> dput()
utils::globalVariables(
  c("N", "TYPE", "detector", "group", "laser", "ord", "sample.id",
    "select.nonsaturating", "select.singlets", "tissue.type")
)

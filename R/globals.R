.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
# dput(sort(c(...)))
# c(...) |> sort() |> dput()
utils::globalVariables(
  c(
    "$CYT", ".", ".mtext", "alias", "configuration", "CREATOR",
    "detector", "fluorochrome", "FSC_A", "group", "hash.md5", "laser",
    "N", "ord", "par", "S", "sample.id", "select.beads", "select.nonsaturating",
    "select.singlets", "SSC_A", "tissue.type", "TYPE", "V", "vars",
    "vector.type")
)

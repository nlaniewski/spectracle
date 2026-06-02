## subsampling raw reference control .fcs files is not sufficient for the purpose of examples;
## need the full (large file size) .fcs which often contain rare events that are needed to test/demonstrate functionality;
## to prevent a large final package/repo size, processing will take place here;
## the relatively much smaller derived spectra can then be saved as .rds objects for use in examples

## OMIP-069 -- human PBMC; 40 color (original panel)
raw.reference.controls <- system.file(
  "extdata/OMIP-069_Raw_Reference Group",
  package = "spectracle"
)
spectra <- spectracle(raw.reference.controls)
##
log.spectracle <- unlist(tictoc::tic.log())
spectracle.timings <- unlist(lapply(tictoc::tic.log(format = F), function(x)
  x$toc - x$tic))
##
saveRDS(
  list(
    raw.reference.controls = list.files(raw.reference.controls),
    timing.fcs.read = spectracle.timings[1],
    timing.spectracle.function = (
      spectracle.timings[length(spectracle.timings)] - spectracle.timings[1]
    ),
    spectra = spectra
  ),
  file = "inst/extdata/prepared_spectra/spectracle_OMIP-069_v1.rds"
)

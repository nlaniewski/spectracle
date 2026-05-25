## subsampling raw reference control .fcs files is not sufficient for the purpose of examples;
## need the full (large file size) .fcs;
## to prevent a large final package/repo size, processing will take place here;
## the relatively much smaller derived spectra can then be saved as .rds objects for use in examples

## expt1 -- mouse spleen challenging samples
raw.reference.controls <- system.file("extdata/expt1", package = "spectracle")
spectra <- spectracle(raw.reference.controls)
## name fix due to no marker '$PnS' names in these source files
## add empirically determined marker names
spectra[
i = N == '450',
j = c('N', 'S') := list('eFluor 450', 'TCRgd')
]
spectra[
i = N == 'BUV615',
j = c('N', 'S') := list('BUV615', 'Siglec F')
]
saveRDS(spectra, "inst/extdata/prepared_spectra/spectra_expt1_nofilter.rds")

spectra <- spectracle(raw.reference.controls, filter.top.expressing = "BUV615 (Cells)")
spectra[
  i = N == '450',
  j = c('N', 'S') := list('eFluor 450', 'TCRgd')
]
spectra[
  i = N == 'BUV615',
  j = c('N', 'S') := list('BUV615', 'Siglec F')
]
saveRDS(spectra, "inst/extdata/prepared_spectra/spectra_expt1_filtered.rds")

## expt2 -- OMIP-069 Expt4 bead and cell controls
raw.reference.controls <- system.file("extdata/expt2", package = "spectracle")
spectra <- spectracle(raw.reference.controls)
saveRDS(spectra, "inst/extdata/prepared_spectra/spectra_expt2.rds")
##

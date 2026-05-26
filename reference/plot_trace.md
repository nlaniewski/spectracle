# Plot a Spectral Trace

Plot a Spectral Trace

## Usage

``` r
plot_trace(spectra, benchmark = FALSE)
```

## Arguments

- spectra:

  [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) – the
  return of [spectracle](spectracle.md).

- benchmark:

  logical – default `FALSE`; if `TRUE`, derived spectra will be
  benchmarked against an internal reference library and the closest
  match (cosine similarity) will be reported as both a black, dashed
  trace and caption text.

## Value

A plot is printed to the active device.

## Examples

``` r

## due to large file sizes, the .fcs files used in this example were first
## fully processed -- spectra was derived and saved; load here

spectra.files <- list.files(
system.file("extdata/prepared_spectra", package = "spectracle"),
full.names = TRUE
)
names(spectra.files) <- sub(".rds", "", basename(spectra.files))

spectra.examples <- lapply(spectra.files, readRDS)

spectra <- spectra <- spectra.examples$spectra_expt2

plot_trace(spectra[sample.id == "CD16 BUV496 (Cells)"])
plot_trace(spectra[sample.id == "CD16 BUV496 (Cells)"], benchmark = TRUE)

```

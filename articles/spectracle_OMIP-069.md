# spectracle – OMIP-069

``` r

library(spectracle)
```

## Objective

Derive spectra using OMIP-069 source data.

## Source Data

Source data as presented in
[OMIP-069](https://pubmed.ncbi.nlm.nih.gov/32830910/) is downloaded from
its public repository:

- [figshare – OMIP-069](https://doi.org/10.6084/m9.figshare.25901635.v1)

Click to see how OMIP-069 source data was obtained.

``` r

## manually download or follow this script

## source directory; this path is user-dependent;
## the author of this article previously downloaded these source files to:
dir.source <- "W:/OMIP-069v2_flowstate/data_source/OMIP-069 40-Color Original Panel Experiment 1/"
## these source files were processed using spectracle();
## the results were saved to an .RDS solely for building this article

## BELOW IS THE METHOD USED TO DOWNLOAD/VERIFY SOURCE

## figshare name: Exp1_Original Panel_Donor1_FreshPBMCs_Fig6
## Experiment ($PROJ) name: OMIP-069 40-Color Original Panel Experiment 1
## Size: 1.17 GB
## MD5 checksum: 27a3a73c7669f4f4a0c8f32da8fa4dfd
## download link: https://ndownloader.figshare.com/files/46539499
## Citation: Bonilla, Diana Lucia (2024). 40-color commentary. 
## figshare. Preprint. https://doi.org/10.6084/m9.figshare.25901635.v1

## test for existence and integrity of source files (.fcs)
if(dir.exists(dir.source)){
  digests <- sapply(list.files(dir.source, full.names = T, recursive = T, pattern = ".fcs"),
                    digest::digest,
                    file = TRUE)
  digest <- digest::digest(as.vector(digests))
  if(digest != "f4c33799be919e6188be387c0b60ec34"){
    message("Digest of source files (.fcs) does not match historic value;\ndownloading source files...")
    download.source <- TRUE
  }else{
    download.source <- FALSE
  }
}else{
  download.source <- TRUE
}
if(download.source){
  temp <- tempfile(tmpdir = "data_temp/", fileext = ".zip")
  zip.link <- "https://ndownloader.figshare.com/files/46539499"
  curl::curl_download(
    url = zip.link,
    destfile = temp,
    quiet = FALSE,
    mode = "wb"
  )
  ## check the hash
  if (tools::md5sum(temp) != "27a3a73c7669f4f4a0c8f32da8fa4dfd") {
    stop("MD5 checksum (hash) of the downloaded .zip does not match the historic value.")
  }
  ## files contained in the .zip: 89 files;
  ## Raw, Unmixed, and SpectroFlo software files (.Expt, .WTML, and .UST)
  utils::unzip(temp, list = TRUE)
  ## unzip to "./data_temp"
  utils::unzip(temp, files = NULL, exdir = "./data_temp/")
  ## get the original 'Project' ('$PROJ') name to reconstruct directory structure
  proj <- "./data_temp/Raw/Reference Group/Unstained (Cells).fcs" |>
    flowstate:::readFCStext() |>
    flowstate:::keywords.to.data.table() |>
    _[][['$PROJ']]
  ## R scripts directory
  dir.analysis <- file.path("R", proj)
  if (!dir.exists(dir.analysis))
    dir.create(dir.analysis)
  ## create the source directory
  dir.out <- file.path("data_source/", proj)
  if (!dir.exists(dir.out))
    dir.create(dir.out)
  ## move the source .fcs files; Raw and Unmixed
  fcs.to.move <- list.files(
    "data_temp/",
    full.names = T,
    pattern = ".fcs",
    recursive = T
  )
  dirnames <- dirname(fcs.to.move) |> unique()
  sapply(dirnames, function(i) {
    dir.out.sub <- sub("data_temp", dir.out, i)
    if (!dir.exists(dir.out.sub))
      dir.create(dir.out.sub, recursive = T)
  }) |> invisible()
  file.rename(from = fcs.to.move,
              to = sub("data_temp", dir.out, fcs.to.move)) |> invisible()
  ## move the SpectroFlo experiment files
  dir.out.expt <- file.path(dir.out, "SpectroFlo_experiment")
  if (!dir.exists(dir.out.expt))
    dir.create(dir.out.expt)
  files.to.move <- list.files("data_temp/", full.names = T, pattern = "Expt|UST|WTML")
  file.rename(from = files.to.move,
              to = sub("data_temp", dir.out.expt, files.to.move)) |> invisible()
  ## cleanup
  unlink(temp)
  unlink("data_temp/Raw/", recursive = T)
  unlink("data_temp/Unmixed//", recursive = T)
}
##
```

### Source Files – Raw Reference Controls (Cells)

‘./Raw/Reference Group’ file paths

``` r

dir.source.raw <- grep("Raw.*Ref", list.dirs(dir.source), value = T)
raw.reference.controls <- list.files(
  dir.source.raw,
  full.names = T,
  pattern = ".fcs"
)
```

## Derive Spectra

[`spectracle()`](https://nlaniewski.github.io/spectracle/reference/spectracle.md)
is fully automated and by design has limited function arguments; simply
pass the directory where the source files are contained.

``` r

spectra <- spectracle(dir.source.raw)
```

### Spectra – OMIP-069

[`spectracle()`](https://nlaniewski.github.io/spectracle/reference/spectracle.md)
returns a standardized
[`data.table::data.table()`](https://rdrr.io/pkg/data.table/man/data.table.html)
containing relevant metadata:

- deriving function
- date-of-acquisition
- cytometer, serial number, software
- sample identifier (sample name, fluor, marker/stain, peak detector,
  laser)
- derived, normalized \[0,1\] spectra
- hash summaries for each spectra

``` r

spectra[]
```

### Spectra – Hashes

Derived using a fully data-driven/automated process, the resultant
spectra should reproduce given the same source data and function
version. A summary digest is generated from the individual hashes:

- Hash (md5) summary (historic): *bd830ba7113eca131f3c9c8229b6a543*

``` r

## hash (md5) summary for this instance
spectra[, digest::digest(hash.md5)]
#> [1] "bd830ba7113eca131f3c9c8229b6a543"
```

``` r

spectra[, .(N, S, detector, hash.md5)]
```

### Spectra – Namefix

As `spectra` is intended for use in unmixing raw data, the `N` and `S`
names will be used to form a conformant parameter name (`$PnN` and
`$PnS`); once derived and a `spectra` object is available in the
environment, these names can be edited.

``` r

spectra[
  i = sample.id == "LIVE DEAD Blue (Cells)",
  j = c('N', 'S') := list('LIVE DEAD Blue', 'Viability')
] |> invisible()
```

### Spectra – Traces

``` r

p <- plot_trace(spectra, plot.type = 'plotly')
```

#### UV

``` r

p$UV
```

#### V

``` r

p$V
```

#### B

``` r

p$B
```

#### YG

``` r

p$YG
```

#### R

``` r

p$R
```

#### Individual Traces – Benchmarked

A few choice individual traces (purple) – benchmarked against a
‘reference library’ (black).

##### UV7::BUV496::CD16

``` r

plot_trace(spectra[N == "BUV496"], benchmark = TRUE)
```

[![](spectracle_OMIP-069_files/figure-html/UV7-BUV496-CD16-1.png)](https://nlaniewski.github.io/spectracle/articles/spectracle_OMIP-069_files/figure-html/UV7-BUV496-CD16-1.png)

##### V7::BV510::CD3

``` r

plot_trace(spectra[N == "BV510"], benchmark = TRUE)
```

[![](spectracle_OMIP-069_files/figure-html/V7-BV510-CD3-1.png)](https://nlaniewski.github.io/spectracle/articles/spectracle_OMIP-069_files/figure-html/V7-BV510-CD3-1.png)

##### YG10::PE-Fire 810::TIGIT

``` r

plot_trace(spectra[N == "PE-Fire 810"], benchmark = TRUE)
```

[![](spectracle_OMIP-069_files/figure-html/YG10-PE-Fire810-TIGIT-1.png)](https://nlaniewski.github.io/spectracle/articles/spectracle_OMIP-069_files/figure-html/YG10-PE-Fire810-TIGIT-1.png)

### Spectra – Download

Download OMIP-069 Spectra

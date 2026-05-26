# spectracle – OMIP-109

``` r

library(spectracle)
```

## Objective

Derive spectra using OMIP-109 source data.

## Source Data

Source data as presented in
[OMIP-109](http://www.ncbi.nlm.nih.gov/pubmed/39466962) is downloaded
from its public repository:

- [figshare – OMIP-109](https://doi.org/10.6084/m9.figshare.25699221)

Click to see how OMIP-109 source data was obtained.

``` r

## manually download or follow this script

## source directory; this path is user-dependent;
## the author of this article previously downloaded these source files to:
dir.source <- "W:/OMIP-109_flowstate/data_source/OMIP-XXX 45-Color Experiment 1/"
## these source files were processed using spectracle();
## the results were saved to an .RDS solely for building this article

## BELOW IS THE METHOD USED TO DOWNLOAD/VERIFY SOURCE

## figshare name: 45-Color Full Spectrum Flow Cytometry Panel
## Experiment ($PROJ) name: OMIP-XXX 45-Color Experiment 1
## Size: 1.74 GB
## MD5 checksum: 623c49ca26c4584a3bcc44b6a3fa5843
## download link: https://ndownloader.figshare.com/files/45873810
## Citation: Jaimes, Maria (2024). 45-Color Full Spectrum Flow Cytometry Panel. 
## figshare. Dataset. https://doi.org/10.6084/m9.figshare.25699221.v1


## test for existence and integrity of source files (.fcs)
if(dir.exists(dir.source)){
  digests <- sapply(list.files(dir.source, full.names = T, recursive = T, pattern = ".fcs"),
                    digest::digest,
                    file = TRUE)
  digest <- digest::digest(as.vector(digests))
  if(digest != "fbe1796cfccead8c4096f942826f6613"){
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
  zip.link <- "https://ndownloader.figshare.com/files/45873810"
  curl::curl_download(
    url = zip.link,
    destfile = temp,
    quiet = FALSE,
    mode = "wb"
  )
  ## check the hash
  if (tools::md5sum(temp) != "623c49ca26c4584a3bcc44b6a3fa5843") {
    stop("MD5 checksum (hash) of the downloaded .zip does not match the historic value.")
  }
  ## files contained in the .zip: 100 files;
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

``` r

res <- readRDS(
  system.file(
    "extdata/prepared_spectra/spectracle_OMIP-109.rds",
    package = "spectracle"
  )
)
```

### Source Files – Raw Reference Controls (Cells)

‘./Raw/Reference Group’ file paths

``` r

# dir.source.raw <- grep("Raw.*Ref", list.dirs(dir.source), value = T)
# raw.reference.controls <- list.files(
#   dir.source.raw, 
#   full.names = T, 
#   pattern = ".fcs"
# )
##
data.frame(
  "./Raw/Reference Group/" = res$raw.reference.controls, 
  check.names = F
)
```

## Derive Spectra

[`spectracle()`](https://nlaniewski.github.io/spectracle/reference/spectracle.md)
is fully automated and by design has limited function arguments; simply
pass the directory where the source files are contained.

``` r

# spectra <- spectracle(dir.source.raw)
spectra <- res$spectra
```

### Spectra – OMIP-109

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

- Hash (md5) summary (historic): *cc8be5bd56a9f25337ce09d68ff80a45*

``` r

## hash (md5) summary for this instance
spectra[, digest::digest(hash.md5)]
#> [1] "cc8be5bd56a9f25337ce09d68ff80a45"
```

``` r

spectra[, .(N, S, detector, hash.md5)]
```

### Spectra – Traces

- [UV](#tabset-7-1)
- [V](#tabset-7-2)
- [B](#tabset-7-3)
- [YG](#tabset-7-4)
- [R](#tabset-7-5)
- [AF](#tabset-7-6)

&nbsp;

- - [CD11b Spark UV387 (Cells)](#tabset-1-1)
  - [CD45RA BUV395 (Cells)](#tabset-1-2)
  - [LIVE DEAD Blue (Cells)](#tabset-1-3)
  - [CD16 BUV496 (Cells)](#tabset-1-4)
  - [CCR5 BUV563 (Cells)](#tabset-1-5)
  - [CD314 BUV615 (Cells)](#tabset-1-6)
  - [CD39 BUV661 (Cells)](#tabset-1-7)
  - [CD38 BUV737 (Cells)](#tabset-1-8)
  - [CD8 BUV805 +BSB (Cells)](#tabset-1-9)

  ![](spectracle_OMIP-109_files/figure-html/plot-traces-1.png)

  ![](spectracle_OMIP-109_files/figure-html/plot-traces-2.png)

  ![](spectracle_OMIP-109_files/figure-html/plot-traces-3.png)

  ![](spectracle_OMIP-109_files/figure-html/plot-traces-4.png)

  ![](spectracle_OMIP-109_files/figure-html/plot-traces-5.png)

  ![](spectracle_OMIP-109_files/figure-html/plot-traces-6.png)

  ![](spectracle_OMIP-109_files/figure-html/plot-traces-7.png)

  ![](spectracle_OMIP-109_files/figure-html/plot-traces-8.png)

  ![](spectracle_OMIP-109_files/figure-html/plot-traces-9.png)

- [CCR7 BV421 (Cells)](#tabset-2-1)
- [CD123 Super Bright 436 (Cells)](#tabset-2-2)
- [CD11c eFluor 450 (Cells)](#tabset-2-3)
- [HLA-DR BV480 (Cells)](#tabset-2-4)
- [CD3 BV510 (Cells)](#tabset-2-5)
- [CD20 Spark Violet 538 (Cells)](#tabset-2-6)
- [IgM BV570 (Cells)](#tabset-2-7)
- [IgG BV605 (Cells)](#tabset-2-8)
- [CD28 BV650 (Cells)](#tabset-2-9)
- [CCR6 BV711 (Cells)](#tabset-2-10)
- [CXCR5 BV750 (Cells)](#tabset-2-11)
- [PD-1 BV785 (Cells)](#tabset-2-12)

![](spectracle_OMIP-109_files/figure-html/plot-traces-10.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-11.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-12.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-13.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-14.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-15.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-16.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-17.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-18.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-19.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-20.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-21.png)

- [TIM-3 BB515 (Cells)](#tabset-3-1)
- [CD57 cFluor B532 (Cells)](#tabset-3-2)
- [CD14 Spark Blue 550 (Cells)](#tabset-3-3)
- [CD2 PerCP-Cy5.5 (Cells)](#tabset-3-4)
- [TCRgd PerCP-Vio700 (Cells)](#tabset-3-5)
- [CD4 PerCP-Fire 806 (Cells)](#tabset-3-6)
- [DNAM-1 RealBlue 780 (Cells)](#tabset-3-7)

![](spectracle_OMIP-109_files/figure-html/plot-traces-22.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-23.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-24.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-25.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-26.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-27.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-28.png)

- [CD56 cFluor YG584 (Cells)](#tabset-4-1)
- [GRP56 PE (Cells)](#tabset-4-2)
- [CD24 cFluor YG610 (Cells)](#tabset-4-3)
- [CD337 PE-Dazzle594 (Cells)](#tabset-4-4)
- [CD103 PE-Fire 640 250ng (Cells)](#tabset-4-5)
- [CD45 PerCP (Cells)](#tabset-4-6)
- [CD95 PE-Cy5 (Cells)](#tabset-4-7)
- [CD25 cFluor BYG710 (Cells)](#tabset-4-8)
- [IgD cFluor BYG750 (Cells)](#tabset-4-9)
- [CXCR3 PE-Cy7 (Cells)](#tabset-4-10)
- [TIGIT PE-Fire 810 (Cells)](#tabset-4-11)

![](spectracle_OMIP-109_files/figure-html/plot-traces-29.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-30.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-31.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-32.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-33.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-34.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-35.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-36.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-37.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-38.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-39.png)

- [CD161 APC 500ng (Cells)](#tabset-5-1)
- [CD1c Alexa Fluor 647 (Cells)](#tabset-5-2)
- [CD19 Spark NIR685 (Cells)](#tabset-5-3)
- [CD127 cFluor R720 (Cells)](#tabset-5-4)
- [CD27 APC-H7 (Cells)](#tabset-5-5)
- [KLRG-1 APC-Fire 810 (Cells)](#tabset-5-6)

![](spectracle_OMIP-109_files/figure-html/plot-traces-40.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-41.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-42.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-43.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-44.png)

![](spectracle_OMIP-109_files/figure-html/plot-traces-45.png)

- [Unstained (Cells)](#tabset-6-1)

![](spectracle_OMIP-109_files/figure-html/plot-traces-46.png)

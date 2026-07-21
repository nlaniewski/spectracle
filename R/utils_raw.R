## helper function: check for conformity and if .fcs files are actually raw;
## stop otherwise to prevent reading of the files
## if tests pass, return software-specific metadata: identifiers and APD limits
## this all depends on if software/instrument adheres to FCS 3+ file standard
# check.fcs.files.raw <- function(fcs.file.paths){
#   ## list of keywords for each file
#   kw <- sapply(fcs.file.paths, AutoSpectral::readFCSheader)
#   kw <- lapply(kw, as.list)
#   ## list of parameters for each file; derived from keywords
#   parameters <- lapply(kw, parameters.to.data.table)
#   ## CREATOR (software) -- used to enforce expected naming convention/conformity
#   software <- unique(sapply(kw, '[[', 'CREATOR'))
#   if(length(software) != 1){
#     stop("Mixed software detected; are these .fcs files from the same instrument?")
#   }
#   ## software-specific tests
#   if(grepl("SpectroFlo", software)){
#     ## check TYPE
#     res <- sapply(parameters, function(.parameters){
#       .parameters[, any(grepl("Raw_Fluorescence", TYPE))]
#     })
#     ## are the files actually raw?
#     if(!all(res)){
#       stop(
#         sprintf("Non-conformant %s files: expect 'Raw_Fluorescence' data.types",
#                 software
#         )
#       )
#     }
#     ## do they have the same parameter names/number of detectors?
#     res <- sapply(parameters, function(.parameters){
#       .parameters[, .(N)]
#     })
#     if(length(unique(res)) != 1){
#       stop(
#         sprintf("Non-conformant %s files: differing names/number of parameters/detectors.",
#                 software
#         )
#       )
#     }
#     ## SpectroFlo identifier
#     ids <- as.vector(sapply(kw, '[[', 'TUBENAME'))
#     ## Cytek Aurora detector (APD) max -- linear value
#     val.saturating <- 4194304
#     ## detectors
#     detectors <- sapply(parameters, function(.parameters){
#       .parameters[TYPE == "Raw_Fluorescence", .(N)]
#     })
#     detectors <- unlist(unique(detectors))
#     ## scatter
#     scatter <- sapply(parameters, function(.parameters){
#       .parameters[grepl("scatter", TYPE, ignore.case = T), .(N)]
#     })
#     scatter <- unlist(unique(scatter))
#   }else if(grepl("some other software", software)){
#     message("some other software-specific code here to enforce naming convention/conformity")
#   }
#   ## return data-driven variables and proceed if files are raw...
#   return(
#     list(
#       identifiers = ids,
#       val.saturating = val.saturating,
#       detectors = detectors,
#       scatter = scatter,
#       keywords = kw,
#       parameters = parameters
#     )
#   )
# }

# raw.fluorescence.check <- function(fcs.file.paths){
#   invisible(sapply(fcs.file.paths, function(i){
#     parms <- parameters.to.data.table(readFCStext(i))
#     if(!'TYPE' %in% names(parms)){
#       stop("Are these .fcs files processed by flowstate? 'TYPE' is missing from [['parameters']].")
#     }else{
#       if(!"Raw_Fluorescence" %in% parms[['TYPE']]){
#         stop("'Raw_Fluorescence' value not found in 'TYPE' keyword; is this a raw reference control?")
#       }
#     }
#   }))
# }

## function to add metadata to [['keywords]] and [['data']]
reference.group.keywords <- function(flowstate){
  ## column identifier for [['data']] and [['keywords']]
  if(!'sample.id' %in% intersect(names(flowstate$data), names(flowstate$keywords))){
    stop("'sample.id' identifier not found.")
  }
  ## CREATOR/software
  software <- flowstate$keywords[, unique(CREATOR)]
  if(length(software) != 1){
    warning(
      sprintf("Mixed software: %s", paste(software, collapse = " ; "))
    )
    warning(
      "might not be able to reliably derive metadata"
    )
  }
  ## derive software-specific metadata based on established naming convention
  if(all(grepl("spectroflo", software, ignore.case = T))) reference.group.keywords.spectroflo(flowstate)
  ## prepend a new class for downstream workflow purposes
  data.table::setattr(flowstate, 'class', c("reference.group", class(flowstate)))
  ## return
  invisible(flowstate)
}

reference.group.keywords.spectroflo <- function(flowstate){
  ## add keyword metadata based on splitting sample.id;
  ## following SpectroFlo naming convention: marker(S) fluorophore(N) (type) -- literal spaces separating
  keywords.to.add <- c('tissue.type', 'N', 'S')
  ## test
  if(any(flowstate$data[, unique(sample.id)] != flowstate$keywords[, unique(sample.id)])){
    stop("Sample order does not match between [['data']] and [['keywords']].")
  }
  ## add keyword/value pairs based on gsub/regex; dependent on SpectrolFlo naming convention
  ## updates [['keywords']]
  flowstate$keywords[
    ,
    j = (keywords.to.add) := {
      ##
      if(!any(flowstate$keywords[, grepl("(.*)", sample.id)])){
        warning("Non-conformant SpectroFlo name: no (Beads) and/or (Cells) designation.")
      }
      ##
      tissue.type <- factor(gsub("^.*\\((.*?)\\).*$", "\\1", sample.id))
      res <- sub(" \\(.*$", "", sample.id)
      res[grep("unstained|negative", res, ignore.case = T)] <- "Unstained"
      ##
      marker <- factor(ifelse(
        tissue.type %in% c("Beads", "Cells") & grepl("Unstained", res),
        NA,
        sub(" +.*$", "", res)
      ))
      ##
      fluorophore <- factor(ifelse(
        tissue.type %in% c("Beads", "Cells") & is.na(marker),
        "AF",
        sub(".*? ", "", res)
      ))
      ##
      list(tissue.type, N = fluorophore, S = marker)
    }
  ]
  ## updates [['data']]
  flowstate::add.keywords.to.data(flowstate, keywords.to.add)
}

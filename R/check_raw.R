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


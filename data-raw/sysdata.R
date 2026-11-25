# Build the internal data stored in R/sysdata.rda.
# Run from the package root in a clean session: source("data-raw/sysdata.R")

# Nuclide parameters --------------------------------------------------------
# One row per parameter. Sources are given per row (`reference`, Rdpack keys
# from inst/REFERENCES.bib); `checked` is TRUE when the value was checked
# against the cited source during the v1 rewrite.
nuclide_parameters <- utils::read.csv(
  "data-raw/nuclide_parameters.csv",
  na.strings = "",
  colClasses = c(
    nuclide = "character", target = "character", parameter = "character",
    pathway = "character", element = "character", scaling = "character",
    value = "numeric", sd = "numeric", unit = "character",
    reference = "character", checked = "logical", note = "character"
  )
)

local({
  p <- nuclide_parameters
  stopifnot(
    all(p$parameter %in% c("half_life", "production", "attenuation")),
    all(p$scaling %in% c(NA, "St", "Lm", "Sf", "Sa")),
    all(is.na(p$value) | p$value >= 0),
    all(is.na(p$sd) | p$sd >= 0),
    !anyDuplicated(p[c("nuclide", "target", "parameter", "pathway",
                       "element", "scaling")]),
    # every value has a source, except stable-nuclide half-lives
    all(is.na(p$value) | !is.na(p$reference) | is.infinite(p$value))
  )
  # every reference key must exist in inst/REFERENCES.bib
  bib <- readLines("inst/REFERENCES.bib", encoding = "UTF-8")
  keys <- sub("^@[a-zA-Z]+\\{([^,]+),.*$", "\\1", grep("^@", bib, value = TRUE))
  used <- unique(unlist(strsplit(stats::na.omit(p$reference), ";")))
  missing <- setdiff(used, keys)
  if (length(missing)) {
    stop("Unknown reference keys: ", paste(missing, collapse = ", "))
  }
})

usethis::use_data(nuclide_parameters, internal = TRUE, overwrite = TRUE)

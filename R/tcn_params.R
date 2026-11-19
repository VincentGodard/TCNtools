#' Parameters for a nuclide, target and scaling scheme
#'
#' Builds the set of physical parameters needed by the production and
#' concentration functions of the package, for one nuclide in one target
#' and one scaling scheme: half-life and decay constant, sea level high
#' latitude (SLHL) reference production rates, and attenuation lengths.
#' Every value carries its source.
#'
#' Spallation production rates are the calibrated values of
#' \insertCite{borchers2016geological;textual}{TCNtools} (Table 7) for the
#' chosen scaling scheme. They depend on the scheme, so a rate must only be
#' used with the scaling factors of the same scheme. Their uncertainties for
#' 10Be and 26Al are those of \insertCite{phillips2016cronusearth;textual}{TCNtools}
#' (Tables 3 and 6). No uncertainty is given for 14C.
#'
#' Muon production rates are for a simple exponential approximation of
#' production at depth, with one term for stopped (negative) muons and one
#' for fast muons. They do not depend on the spallation scaling scheme.
#' 10Be and 26Al values are from \insertCite{braucher2011production;textual}{TCNtools}
#' (Table 6, sea level), 14C values from
#' \insertCite{lupker2015depthdependence;textual}{TCNtools} (SLHL).
#'
#' Attenuation lengths are not nuclide specific: 160, 1500 and 4320 g/cm2
#' for spallation, stopped muons and fast muons, as used by
#' \insertCite{braucher2011production;textual}{TCNtools}.
#'
#' 21Ne is not calibrated by Borchers et al. (2016). Its spallation rates
#' are those of the CRONUScalc program
#' \insertCite{marrero2016cosmogenic}{TCNtools}, which applies the
#' 21Ne/10Be production ratio of \insertCite{balco2009al26;textual}{TCNtools}
#' to the Borchers et al. (2016) 10Be rates. Independent estimates agree
#' within a few percent (e.g. 17.0 +/- 1.1 at/g/yr for St, 2 sigma, from the
#' SPICE calibration). The 21Ne fast-muon rate (0.203 at/g/yr at sea level,
#' no stopped-muon production) is that of
#' \insertCite{balco2019cosmogenic;textual}{TCNtools}, needed for long
#' exposure or exhumation histories where production at depth matters.
#' CRONUScalc models no 21Ne muons, so its spallation rate already includes
#' them: adding the muon term overestimates surface production by about
#' 1.2%, well within the uncertainty of the rate.
#'
#' If a value is missing (`NA`), the function warns. Values can be
#' supplied or replaced through `...`.
#'
#' @param nuclide Character scalar. One of `"Be10"`, `"Al26"`, `"C14"`,
#'   `"Ne21"`.
#' @param scaling Character scalar. Scaling scheme of the spallation
#'   production rate, following the notation of
#'   \insertCite{borchers2016geological;textual}{TCNtools}: `"St"` (Lal 1991,
#'   Stone 2000), `"Lm"` (time-dependent Lal/Stone), `"Sf"` (Lifton-Sato,
#'   total flux) or `"Sa"` (Lifton-Sato, nuclide specific). No default: the
#'   scheme must be chosen explicitly.
#' @param target Character scalar. Target mineral. Only `"quartz"` is
#'   available.
#' @param ... Named overrides of individual parameters. Names are those
#'   listed in the `parameter` column of `as.data.frame()` on the result
#'   (e.g. `half_life`, `production_spallation`, `production_muon_stop`,
#'   `attenuation_spallation`). Each value is a numeric of length 1 (value)
#'   or 2 (value and standard deviation), in the units of that parameter.
#'   An overridden parameter gets the reference `"user"`.
#'
#' @return An object of class `tcn_params`, a list with elements
#' \describe{
#'   \item{`nuclide`}{list with `name`, `half_life` (yr) and `lambda`
#'     (decay constant, 1/yr), the last two as `c(value, sd)`. `lambda` is
#'     0 for a stable nuclide.}
#'   \item{`target`}{target mineral.}
#'   \item{`scaling`}{scaling scheme.}
#'   \item{`production`}{data frame of SLHL production rates (at/g/yr), one
#'     row per pathway (`spallation`, `muon_stop`, `muon_fast`).}
#'   \item{`attenuation`}{data frame of attenuation lengths (g/cm2), one row
#'     per pathway.}
#'   \item{`source`}{data frame with one row per parameter: value, standard
#'     deviation, unit, reference key, whether the value was checked against
#'     its source, and notes.}
#' }
#'
#' @section Vectorisation:
#' All arguments are scalars: one call describes one nuclide.
#'
#' @seealso [tcn_params_list()] for several nuclides with one scheme,
#'   [tcn_params_table()] for the table of all values and their sources.
#'
#' @references
#' \insertAllCited{}
#'
#' Half-lives: \insertRef{chmeleff2010determination}{TCNtools},
#' \insertRef{korschinek2010new}{TCNtools},
#' \insertRef{nishiizumi2004preparation}{TCNtools},
#' \insertRef{kutschera2019halflife}{TCNtools}.
#'
#' @examples
#' be <- tcn_params("Be10", scaling = "St")
#' be
#' be$production
#' be$nuclide$lambda
#'
#' # override a value, e.g. a local production rate calibration
#' tcn_params("Be10", scaling = "St", production_spallation = c(4.2, 0.2))
#'
#' # all parameters with their sources
#' as.data.frame(tcn_params("C14", scaling = "Sa"))
#' @export
tcn_params <- function(nuclide, scaling, ..., target = "quartz") {
  tab <- nuclide_parameters
  check_choice(nuclide, sort(unique(stats::na.omit(tab$nuclide))))
  check_choice(scaling, scaling_schemes)
  check_choice(target, sort(unique(stats::na.omit(tab$target))))

  keep <- ((tab$nuclide %in% nuclide & tab$target %in% target) |
             is.na(tab$nuclide)) &
    (is.na(tab$scaling) | tab$scaling %in% scaling)
  rows <- tab[keep, ]
  rows$key <- ifelse(rows$parameter == "half_life", "half_life",
                     paste(rows$parameter, rows$pathway, sep = "_"))
  rows <- apply_overrides(rows, list(...))

  missing <- rows$key[is.na(rows$value)]
  if (length(missing) > 0) {
    warning(sprintf(
      "No value available yet for %s %s: %s. Supply them through `...`.",
      nuclide, target, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }

  hl <- rows[rows$parameter == "half_life", ]
  prod <- rows[rows$parameter == "production", ]
  att <- rows[rows$parameter == "attenuation", ]

  x <- structure(
    list(
      nuclide = list(
        name = nuclide,
        half_life = c(value = hl$value, sd = hl$sd),
        lambda = c(value = log(2) / hl$value,
                   sd = log(2) * hl$sd / hl$value^2)
      ),
      target = target,
      scaling = scaling,
      production = data.frame(pathway = prod$pathway, value = prod$value,
                              sd = prod$sd, unit = prod$unit),
      attenuation = data.frame(pathway = att$pathway, value = att$value,
                               sd = att$sd, unit = att$unit),
      source = data.frame(parameter = rows$key, value = rows$value,
                          sd = rows$sd, unit = rows$unit,
                          reference = rows$reference, checked = rows$checked,
                          note = rows$note)
    ),
    class = "tcn_params"
  )
  validate_tcn_params(x)
}

#' Parameters for several nuclides
#'
#' Calls [tcn_params()] once per nuclide, with the same scaling scheme and
#' target, and returns the results as a named list. Use it for multi-nuclide
#' calculations (e.g. 26Al/10Be), so that every nuclide uses the same scheme.
#'
#' @param nuclides Character vector of distinct nuclide names, each one of
#'   the values accepted by [tcn_params()].
#' @param scaling Character scalar. Scaling scheme, as in [tcn_params()].
#' @param target Character scalar. Target mineral, as in [tcn_params()].
#' @param ... Named overrides passed to [tcn_params()] for every nuclide.
#'   Only parameters shared by all nuclides are accepted (the attenuation
#'   lengths `attenuation_spallation`, `attenuation_muon_stop`,
#'   `attenuation_muon_fast`): a nuclide-specific value such as
#'   `production_spallation` or `half_life` would otherwise be applied to
#'   every nuclide. For a nuclide-specific override, call [tcn_params()] for
#'   that nuclide and replace the list element (see examples).
#'
#' @return A named list of `tcn_params` objects, one per nuclide, in the
#'   order of `nuclides`.
#'
#' @section Vectorisation:
#' `nuclides` is a vector, `scaling` and `target` are scalars.
#'
#' @seealso [tcn_params()] for one nuclide, [tcn_params_table()] for the
#'   table of all values and their sources.
#'
#' @examples
#' p <- tcn_params_list(c("Be10", "Al26"), scaling = "St")
#' names(p)
#' p$Al26$production
#'
#' # shared override: a different spallation attenuation length for all
#' p <- tcn_params_list(c("Be10", "Al26"), scaling = "St",
#'                      attenuation_spallation = 150)
#' p$Be10$attenuation
#'
#' # nuclide-specific override (e.g. a local 10Be production rate):
#' # `...` does not accept it, since it would apply to every nuclide.
#' # Build the list, then rebuild that nuclide and replace it in the list.
#' p <- tcn_params_list(c("Be10", "Al26"), scaling = "St")
#' p$Be10 <- tcn_params("Be10", scaling = "St",
#'                      production_spallation = c(4.2, 0.2))
#' p$Be10$production   # 4.2, reference "user"
#' p$Al26$production   # unchanged
#' @export
tcn_params_list <- function(nuclides, scaling, ..., target = "quartz") {
  if (!is.character(nuclides) || length(nuclides) == 0 || anyNA(nuclides)) {
    stop("`nuclides` must be a non-empty character vector without NA.",
         call. = FALSE)
  }
  if (anyDuplicated(nuclides)) {
    stop(sprintf("`nuclides` must not contain duplicates: %s.",
                 paste(unique(nuclides[duplicated(nuclides)]),
                       collapse = ", ")), call. = FALSE)
  }
  overrides <- names(list(...))
  tab <- nuclide_parameters
  shared <- unique(paste(tab$parameter, tab$pathway, sep = "_")[is.na(tab$nuclide)])
  specific <- setdiff(overrides[nzchar(overrides)], shared)
  if (length(specific) > 0) {
    stop(sprintf(paste0(
      "`...` only accepts parameters shared by all nuclides (%s), not: %s.\n",
      "To change it for one nuclide, rebuild that nuclide and replace it ",
      "in the list, e.g.\n",
      "  p <- tcn_params_list(%s, scaling = \"%s\")\n",
      "  p$%s <- tcn_params(\"%s\", scaling = \"%s\", %s = ...)"),
      paste(shared, collapse = ", "), paste(specific, collapse = ", "),
      deparse(nuclides), scaling,
      nuclides[1], nuclides[1], scaling, specific[1]),
      call. = FALSE)
  }
  out <- lapply(nuclides, function(n) {
    tcn_params(n, scaling = scaling, ..., target = target)
  })
  names(out) <- nuclides
  out
}

#' @export
print.tcn_params <- function(x, ...) {
  cat(sprintf("<tcn_params> %s in %s, scaling scheme %s\n",
              x$nuclide$name, x$target, x$scaling))
  s <- x$source
  out <- data.frame(
    parameter = s$parameter,
    value = formatC(s$value, digits = 4, format = "g"),
    sd = ifelse(is.na(s$sd), "", formatC(s$sd, digits = 3, format = "g")),
    unit = s$unit,
    reference = short_reference(s$reference)
  )
  print(out, row.names = FALSE, right = FALSE)
  cat(sprintf("decay constant: %s 1/yr\n",
              formatC(x$nuclide$lambda[["value"]], digits = 4, format = "g")))
  invisible(x)
}

# "chmeleff2010determination;korschinek2010new" -> "chmeleff2010, korschinek2010"
short_reference <- function(reference) {
  vapply(strsplit(reference, ";"), function(keys) {
    if (all(is.na(keys))) {
      return("")
    }
    paste(sub("^([a-z]+[0-9]{4}).*$", "\\1", keys), collapse = ", ")
  }, character(1))
}

#' @export
as.data.frame.tcn_params <- function(x, ...) {
  cbind(nuclide = x$nuclide$name, target = x$target, scaling = x$scaling,
        x$source)
}

# Scaling schemes, notation of Borchers et al. (2016)
scaling_schemes <- c("St", "Lm", "Sf", "Sa")

check_choice <- function(x, choices) {
  arg <- deparse(substitute(x))
  if (!is.character(x) || length(x) != 1 || is.na(x) || !x %in% choices) {
    stop(sprintf("`%s` must be one of %s.", arg,
                 paste0("\"", choices, "\"", collapse = ", ")),
         call. = FALSE)
  }
  invisible(x)
}

apply_overrides <- function(rows, overrides) {
  if (length(overrides) == 0) {
    return(rows)
  }
  nm <- names(overrides)
  if (is.null(nm) || any(nm == "")) {
    stop("All parameters in `...` must be named.", call. = FALSE)
  }
  unknown <- setdiff(nm, rows$key)
  if (length(unknown) > 0) {
    stop(sprintf("Unknown parameter(s) in `...`: %s. Valid names: %s.",
                 paste(unknown, collapse = ", "),
                 paste(rows$key, collapse = ", ")),
         call. = FALSE)
  }
  for (k in nm) {
    v <- overrides[[k]]
    if (!is.numeric(v) || !length(v) %in% 1:2 || anyNA(v) || any(v < 0)) {
      stop(sprintf(paste("`%s` must be a non-negative numeric of length 1",
                         "(value) or 2 (value, sd)."), k),
           call. = FALSE)
    }
    i <- match(k, rows$key)
    rows$value[i] <- v[1]
    rows$sd[i] <- if (length(v) == 2) v[2] else NA_real_
    rows$reference[i] <- "user"
    rows$checked[i] <- NA
    rows$note[i] <- "User-supplied value"
  }
  rows
}

validate_tcn_params <- function(x) {
  stopifnot(
    inherits(x, "tcn_params"),
    x$scaling %in% scaling_schemes,
    is.na(x$nuclide$half_life[["value"]]) || x$nuclide$half_life[["value"]] > 0,
    all(is.na(x$production$value) | x$production$value >= 0),
    all(is.na(x$attenuation$value) | x$attenuation$value > 0),
    setequal(x$production$pathway, c("spallation", "muon_stop", "muon_fast")),
    setequal(x$attenuation$pathway, c("spallation", "muon_stop", "muon_fast"))
  )
  x
}

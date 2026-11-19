#' Table of all nuclide parameters and their sources
#'
#' Returns the parameter table used by [tcn_params()]: one row per value,
#' with its uncertainty, unit, reference and notes. Use it to see or report
#' every value the package uses (e.g. `write.csv()` for a paper supplement,
#' or `knitr::kable()` in a report), optionally restricted to some nuclides
#' and scaling schemes.
#'
#' Rows that do not depend on the nuclide (attenuation lengths) have
#' `nuclide` and `target` set to `NA`. Rows that do not depend on the
#' scaling scheme (half-lives, muon rates, attenuation lengths) have
#' `scaling` set to `NA`. Both kinds are always kept by the filters, since
#' they apply to every nuclide and scheme.
#'
#' @param nuclides Character vector of nuclide names (`"Be10"`, `"Al26"`,
#'   `"C14"`, `"Ne21"`), or `NULL` (default) for all.
#' @param scaling Character vector of scaling schemes (`"St"`, `"Lm"`,
#'   `"Sf"`, `"Sa"`), or `NULL` (default) for all.
#' @param target Character vector of target minerals, or `NULL` (default)
#'   for all. Only `"quartz"` is available.
#'
#' @return A data frame with columns `nuclide`, `target`, `parameter`
#'   (`half_life`, `production` or `attenuation`), `pathway` (`spallation`,
#'   `muon_stop`, `muon_fast`), `element` (target element, reserved for
#'   chemistry-dependent nuclides), `scaling`, `value`, `sd` (one standard
#'   deviation, `NA` when none is available), `unit`, `reference` (keys of
#'   `inst/REFERENCES.bib`, separated by `;`), `checked` (value checked
#'   against its source during the v1 rewrite) and `note`. Units: half-lives
#'   in yr, production rates in at/g/yr (sea level, high latitude),
#'   attenuation lengths in g/cm2.
#'
#' @section Vectorisation:
#' All arguments are vectors of choices, combined with AND across arguments.
#'
#' @seealso [tcn_params()] for the parameters of one nuclide and scheme, as
#'   used by the calculation functions.
#'
#' @examples
#' # everything
#' tab <- tcn_parameters()
#' nrow(tab)
#'
#' # spallation rates of 10Be and 26Al in the St scheme
#' tab <- tcn_parameters(c("Be10", "Al26"), scaling = "St")
#' sp <- tab$parameter == "production" & tab$pathway %in% "spallation"
#' tab[sp, c("nuclide", "value", "sd", "reference")]
#'
#' # export, e.g. for a supplementary table
#' \dontrun{
#' write.csv(tcn_parameters(scaling = "Sa"), "tcn_parameters_Sa.csv",
#'           row.names = FALSE)
#' }
#' @export
tcn_parameters <- function(nuclides = NULL, scaling = NULL, target = NULL) {
  tab <- nuclide_parameters
  check_choices(nuclides, sort(unique(stats::na.omit(tab$nuclide))))
  check_choices(scaling, scaling_schemes)
  check_choices(target, sort(unique(stats::na.omit(tab$target))))

  keep <- rep(TRUE, nrow(tab))
  if (!is.null(nuclides)) keep <- keep & (is.na(tab$nuclide) | tab$nuclide %in% nuclides)
  if (!is.null(scaling)) keep <- keep & (is.na(tab$scaling) | tab$scaling %in% scaling)
  if (!is.null(target)) keep <- keep & (is.na(tab$target) | tab$target %in% target)

  out <- tab[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# NULL, or a non-empty character vector of valid choices without NA.
check_choices <- function(x, choices) {
  arg <- deparse(substitute(x))
  if (is.null(x)) {
    return(invisible(x))
  }
  if (!is.character(x) || length(x) == 0 || anyNA(x) || !all(x %in% choices)) {
    stop(sprintf("`%s` must be NULL or contain only %s.", arg,
                 paste0("\"", choices, "\"", collapse = ", ")),
         call. = FALSE)
  }
  invisible(x)
}

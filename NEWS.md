# TCNtools 0.9.0.9000

* Complete rewrite started on branch `v1`. The previous version is frozen
  at tag `v0.1.0`.
* New package infrastructure: roxygen-managed NAMESPACE with explicit
  exports, testthat 3rd edition, Rdpack references, GitHub Actions CI.
* Added `CITATION.cff`, a pkgdown workflow (deploys to `gh-pages` from
  `main` only) and `data-raw/` conventions for building package data.
* New `tcn_params()`: parameters for one nuclide (10Be, 26Al, 14C, 21Ne),
  target (quartz) and scaling scheme (St, Lm, Sf, Sa), with the source of
  every value. Spallation rates from Borchers et al. (2016), muon rates
  from Braucher et al. (2011) and Lupker et al. (2015). Half-lives: 10Be
  1.387 Ma (Chmeleff et al. 2010, Korschinek et al. 2010), 26Al 0.705 Ma
  (Nishiizumi 2004), 14C 5700 yr (Kutschera 2019). 21Ne spallation rates
  from CRONUScalc (Marrero et al. 2016, ratio of Balco and Shuster 2009),
  21Ne fast-muon rate from Balco et al. (2019). Uncertainties of the 10Be and 26Al rates from
  Phillips et al. (2016).
* New `tcn_params_list()`: `tcn_params()` for several nuclides with one
  scaling scheme, as a named list. Overrides are limited to parameters
  shared by all nuclides (attenuation lengths).
* New `tcn_params_table()`: the full table of parameters with their sources,
  optionally filtered by nuclide, scaling scheme and target.
* New `atm_pressure()`: mean atmospheric pressure from site elevation, with
  the standard atmosphere of Stone (2000) (`model = "stone2000"`) or the
  ERA-40 climatology used by CRONUS-Earth and CREp (`model = "era40"`,
  Uppala et al. 2005, Marrero et al. 2016), now vectorised over sites. The
  model has no default. Unlike v0.1.0, the ERA-40 grid is credited to its
  actual source rather than to Balco et al. (2008), which used NCEP-NCAR
  fields.

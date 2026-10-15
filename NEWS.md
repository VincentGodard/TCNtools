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
  from Braucher et al. (2011) and Lupker et al. (2015). The 26Al half-life
  and the 21Ne production rates are not yet filled in.

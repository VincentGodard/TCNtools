# data-raw

Scripts and source files used to build the package data. This directory is
excluded from the built package (`.Rbuildignore`).

Conventions:

- One script per dataset or group of related objects, named after the object
  it builds (e.g. `era40.R`, `geomagnetic.R`).
- Raw source files go in a subdirectory per source (e.g. `data-raw/era40/`),
  with their origin recorded at the top of the script: publication, URL or
  calculator version, retrieval date, and any transformation applied.
- Internal objects (not exported) are assembled by `sysdata.R`, which sources
  the individual scripts and calls `usethis::use_data(..., internal = TRUE)`.
  Exported datasets use `usethis::use_data()` and are documented in
  `R/data.R`.
- Scripts must run from a clean session from the package root, e.g.
  `source("data-raw/sysdata.R")`.

The v0.1.0 data were built from MATLAB `.mat` files in `dev/data/` by
`dev/data/packaging_data_internal.R` (see tag `v0.1.0`). They are migrated
here layer by layer, starting with layer 1.

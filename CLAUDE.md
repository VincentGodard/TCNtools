# TCNtools: v1 rewrite

R package for terrestrial cosmogenic nuclide (TCN) calculations.
Author: Vincent Godard (CEREGE, Aix-Marseille Université).
Repository: https://github.com/VincentGodard/TCNtools (staying on GitHub).

Working language for this package (code, documentation, commit messages,
issues, this file): English.

## Context and decisions

- Complete rewrite of the package: everything is rewritten and checked starting
  from the elementary building blocks, with full documentation, tests and
  examples for every function.
- The historical state is frozen at tag `v0.1.0`. TCNbook, the student lab and
  the Binder link must install that version
  (`remotes::install_github("VincentGodard/TCNtools@v0.1.0")`) until v1 is
  released.
- The rewrite happens on branch `v1`, started from an empty skeleton.
  The old code is browsed through a worktree at `../TCNtools-legacy`
  (`git worktree add ../TCNtools-legacy v0.1.0`) or with
  `git show v0.1.0:R/<file>.R`. Do not copy old code in bulk: each function is
  brought back, reviewed and tested individually.
- One branch per layer (`layer-NN-<name>`, e.g. `layer-04a-scaling-st-lm`),
  merged into `v1` through a pull request.
- Final switch: merge `v1` into `main`, tag `v1.0.0`, GitHub release (Zenodo
  DOI), then move GitHub Pages from `main /docs` to `gh-pages` (the pkgdown
  workflow already deploys to `gh-pages` from `main` only), and update TCNbook.

## Problems found in v0.1.0 (not to be repeated)

- No tests.
- `exportPattern("^[[:alpha:]]+")`: everything was exported, internal functions
  included (`likelihood1`, `prior1`, `proposal1`, `run_mcmc1`,
  `fun_solv_ss_eros_eul`, `LZ`, `E2R`, `Rv0`, `phi_v_slhl`, `d2r`, `r2d`...).
- Short and inconsistent argument names (`p`, `S`, `L`, `prm`, `Lambda`, `h`,
  `w`; `ero` vs `E`; `rho` vs `rhob`), implicit units.
- `docs/` (built pkgdown site) committed to the repository; `dev/`,
  `install.R`, `runtime.txt`, `index.md` missing from `.Rbuildignore`.
- `Depends: R (>= 2.10)`, outdated RoxygenNote.
- Dependencies: `akima` was in Imports but is only used in the
  `depth_profiles` vignette (drop it; if 2D interpolation is ever needed, use
  `interp`, since akima's ACM non-commercial licence is at odds with GPL-3).
  `pracma` only provides `trapz`, `interp2` and `flipud`, and `ellipse` is only
  used by `tnp_ellipse` (eigen-decomposition): both can be replaced by base R.
- Internal data (`sysdata.rda`) built from `.mat` files in `dev/data/`: migrate
  to `data-raw/` with documented sources (conventions in `data-raw/README.md`).
- Inherited TODOs (`dev/notes.md`): vectorisation of `solv_conc_eul`, E = 0 case
  in `conc_soil_mixing`, quartz enrichment factor.

## Conventions for v1

- Explicit exports via roxygen (`@export`); internal functions are marked
  `@noRd` or `@keywords internal`.
- Units fixed and restated in every `@param`: mass depth in g/cm², time in
  years, concentrations in at/g, denudation rate in g/cm²/yr (state it if
  cm/yr is accepted).
- Naming: snake_case, prefixed families (`prod_*`, `scaling_*`, `conc_*`,
  `solve_*`, `tnp_*`), full argument names, same name for the same quantity
  everywhere.
- Parameters: a single validated constructor `tcn_params()` (design in
  "Layer 1 design" below) passed to all downstream functions, replacing the
  v0.1.0 `prm` / `Lambda` pair. Every value carries its provenance.
- Documented vectorisation contract for every function (recycled vs scalar
  arguments).
- Input validation in base R (`stopifnot` or explicit messages); minimal
  dependencies.
- References via Rdpack (`inst/REFERENCES.bib`, taken from `inst/TCNtools.bib`
  at tag v0.1.0) and `\insertRef{}`.

## Definition of done for a function

1. Complete roxygen documentation: description, units, `@references`,
   `@return`.
2. Runnable example.
3. testthat (3rd edition) tests: published reference values or values from
   independent calculators (CRONUS-Earth / Balco's calculator, CREp),
   analytical limits, round-trip tests for inversions.
4. `devtools::check()` with no errors or warnings.
5. Entry in `NEWS.md`.

## Roadmap (layers, in order)

| Layer | Content | Validation |
|---|---|---|
| 0. Infrastructure | DESCRIPTION, roxygen NAMESPACE, testthat 3, Rdpack, `data-raw/`, `.Rbuildignore`, `NEWS.md`, `CITATION.cff`, GitHub Actions (`check-standard`, `test-coverage`, pkgdown → gh-pages) | check passes on an empty package |
| 1. Constants and data | decay constants, attenuation lengths, SLHL production rates; `data-raw/` scripts from the `.mat` files | values from the cited papers |
| 2. Atmosphere / geometry | `atm_pressure` (ERA40 / Stone 2000), angle utilities (internal) | published reference values |
| 3. Geomagnetism | VDM series, Rc computation | point values vs CRONUS-Earth / CREp |
| 4a. Scaling: St, Lm | `scaling_st`, `scaling_lm` | site grid vs online calculator / CREp |
| 4b. Scaling: LSDn | `scaling_lsd`, with the Sato et al. neutron, proton and muon fluxes and the Reedy cross-sections as internal building blocks | site grid vs online calculator / CREp; reference fluxes (`Ref_LSD`) |
| 5. Production at depth | exponential spallation; muons (Balco 2017, models 1a/1b) | Balco 2017 tables; integral identities |
| 6. Forward concentration models | analytical (steady state and transient), then numerical Eulerian / Lagrangian | limits: E = 0 → saturation, t → ∞, pure decay; numerical vs analytical |
| 7. Inversion | steady-state denudation rate and exposure age with uncertainties (root finding) | round trip forward → inverse |
| 8. Two nuclides | `tnp_curves`, `tnp_ellipse`, burial | consistency with layer 6 |

Deferred to v1.x (one GitHub issue each, label `v1.x`): MCMC solver,
`depth_profile_mc`, soil mixing family, chemistry-dependent nuclides
(`tcn_composition()`, `prod_cl36_*`; see "Layer 1 design").

The `dev/test_*.R` scripts at tag v0.1.0 can serve as a starting point for
tests.

### Notes per layer

- **Layer 1:** see "Layer 1 design" below. Reference production rates depend
  on the scaling scheme they were calibrated with (v0.1.0 `prm` has a single
  set, e.g. Be10 Pspal = 4.01 at/g/yr, presumably for St). v0.1.0 `Lambda`
  holds Lspal = 160, Lstop = 1500, Lfast = 4320 g/cm² (exponential muon
  approximation, distinct from the Balco 2017 models of layer 5).
- **Layer 4b:** the Sato fluxes (about 1250 lines in v0.1.0:
  `sato_neutrons`, `sato_neutrons_low`, `sato_protons`, `sato_muons`) are
  required by LSDn, so they belong here as internal functions rather than in
  v1.x. Needs `Ref_LSD` and `XSectsReedyAll` from `data-raw/`.
- **Layer 6:** v0.1.0 has no closed-form concentration functions (only the
  numerical `solv_conc_eul` and `solv_conc_lag`), so the analytical `conc_*`
  functions are new code. They take a production profile from layer 5 (see
  "Layer 1 design"), not hard-coded spallation and muon terms.
- **Example data:** the v0.1.0 exported dataset `tcn_depth_profiles` (from
  `dev/data/data_external/profiles.ods`) can serve as example data for
  layer 7, even with `depth_profile_mc` deferred.

### Layer 1 design: `tcn_params()`

Decided 2026-10-15. The parameter object describes a nuclide × target ×
scaling scheme system, not a nuclide alone (hence not `tcn_nuclide`):
attenuation lengths depend on the material, production rates on the target
and the calibration.

```r
tcn_params(nuclide = "Be10", target = "quartz", scaling = "LSDn", ...)
```

- One nuclide per call (scalar contract); a helper
  `tcn_params_list(c("Be10", "Al26"), ...)` (name to confirm) returns a named
  list. `...` takes named overrides; an overridden value gets source `"user"`.
- Returned S3 object, class `tcn_params`, one sub-list per dependency:

  | Element | Content | Depends on |
  |---|---|---|
  | `$nuclide` | half-life (value, sd), lambda (derived) | nuclide only |
  | `$target` | `"quartz"` or a `tcn_composition` object | |
  | `$scaling` | `"St"`, `"Lm"` or `"LSDn"` | |
  | `$production` | table of production pathways (below) | nuclide × target × scaling |
  | `$attenuation` | shared defaults (spallation 160 g/cm², ...), overridable | material, particle type |
  | `$source` | provenance of every value (Rdpack key or `"user"`) | |

- `$production` is a table of pathways, not a single number: one row per
  pathway (spallation, stopped muons, fast muons, neutron capture) and target
  element, with reference rate, sd, unit, source. Be10 in quartz is the
  degenerate case (one spallation row on O + Si, plus muon rows); 36Cl adds
  spallation rows on Ca, K, Ti, Fe and a low-energy neutron capture row, so
  it fits without an API change.
- Muon parameters for both representations: exponential approximation
  (P, L for stopped and fast muons) and Balco 2017 parameters (layer 5).
  Whether to keep the exponential approximation is still open (see below).
- Methods: `print()` (compact table with references), `validate_tcn_params()`
  (positivity, units, known scheme and target), `as.data.frame()`.
- Storage: one table `data-raw/nuclide_parameters.csv` with columns
  `nuclide, target, parameter, pathway, element, scaling, value, sd, unit,
  reference` (`scaling = NA` for scheme-independent constants), built into
  an internal object read by the constructor.

Chemistry-dependent nuclides (36Cl; also 3He and 21Ne outside quartz),
deferred to v1.x but anticipated now:

- `target` accepts either a mineral name (fixed stoichiometry) or a
  `tcn_composition()` object (target elements Ca, K, Ti, Fe, Cl; bulk rock
  water, B, Gd, Sm, ...). Same parameter class for all downstream functions.
- The constructor only stores and validates numbers. Turning a composition
  and elemental rates into an effective production rate is done by
  dedicated `prod_cl36_*` functions.

Interface between layers 5 and 6: 36Cl neutron-capture production is not
exponential with depth near the surface. Layer 5 therefore returns a
**production profile** (sum of exponential terms, or a function of depth),
and layer 6 models take that profile as input. Quartz nuclides use the
analytical shortcuts (sums of exponentials); other profiles fall back to the
numerical solvers.

Still open before writing layer 1:

1. Production rates: proposed to enter all three schemes in layer 1 (they
   are calibrated literature values), with layers 4a/4b adding the
   consistency test against calibration sites. Not yet confirmed.
2. Nuclides in v1: Be10, Al26, C14 (as in v0.1.0); He3 and/or Ne21 in v1 or
   later?
3. Calibration set per scheme (CRONUS v3 global, CREp defaults, or
   nuclide-specific sets).
4. Al26 half-life: v0.1.0 `lambda` gives 0.708 Ma, Nishiizumi (2004) gives
   0.705 Ma; check against the source before reuse.
5. Keep the exponential muon approximation alongside Balco 2017?

### Where v0.1.0 functions go

| Layer | v0.1.0 material |
|---|---|
| 1 | `prm`, `Lambda`; `dev/data/packaging_data_internal.R` |
| 2 | `atm_pressure` (ERA40), `angdist`, `d2r`, `r2d` |
| 3 | `get_vdm`, `vdm2rc`, `rc_ndp`, `lsdrc`, `solar_modulation`; `GMDB`, `Pal_LSD` |
| 4a | `scaling_st`, `scaling_lm` |
| 4b | `scaling_lsd`, `sato_*`; `Ref_LSD`, `XSectsReedyAll` |
| 5 | `mu_model1a`, `mu_model1b`, `LZ`, `Rv0`, `E2R`, `phi_v_slhl` |
| 6 | `solv_conc_eul`, `solv_conc_lag` |
| 7 | `solv_ss_eros_eul` |
| 8 | `tnp_curves`, `tnp_ellipse`, `tnp_burial` |
| v1.x | `likelihood1`, `prior1`, `posterior1`, `proposal1`, `run_mcmc1`, `depth_profile_mc`, `soil_mixing.R` |

## Work cycle per layer

```bash
git switch v1 && git switch -c layer-NN-<name>
# write R/, tests/testthat/, examples, vignette if relevant
# devtools::document(); devtools::test(); devtools::check()
git push -u origin layer-NN-<name>   # then PR into v1
```

## Progress

- [x] Tag `v0.1.0` set and pushed; README and Binder point to it
- [ ] TCNbook and student lab point to `v0.1.0` (not verifiable from this repo)
- [x] Branch `v1` created, empty skeleton, legacy worktree
- [x] Layer 0: infrastructure
- [ ] Layer 1: constants and data
- [ ] Layer 2: atmosphere / geometry
- [ ] Layer 3: geomagnetism
- [ ] Layer 4a: scaling St, Lm
- [ ] Layer 4b: scaling LSDn
- [ ] Layer 5: production at depth
- [ ] Layer 6: forward models
- [ ] Layer 7: inversion
- [ ] Layer 8: two nuclides
- [ ] Switch `v1` → `main`, `v1.0.0`, Zenodo DOI, GitHub Pages, TCNbook
- [ ] Remote branch `test_holepunch` (2 commits, Binder test only: empty
      `.binder/Dockerfile` and a README line) can be deleted;
      `variable_erosion` is fully contained in `v0.1.0`

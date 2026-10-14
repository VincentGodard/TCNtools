# TCNtools — refonte v1

Package R de calculs sur les nucléides cosmogéniques terrestres (TCN).
Auteur : Vincent Godard (CEREGE, Aix-Marseille Université).
Dépôt : https://github.com/VincentGodard/TCNtools (on reste sur GitHub).

## Contexte et décisions prises

- Refonte complète du package : tout réécrire et vérifier en partant des briques
  élémentaires, avec documentation, tests et exemples complets pour chaque fonction.
- L'état historique est figé au tag `v0.1.0`. TCNbook, le TP étudiants et le lien
  Binder doivent installer cette version (`remotes::install_github("VincentGodard/TCNtools@v0.1.0")`)
  tant que la v1 n'est pas publiée.
- La refonte se fait sur la branche `v1`, partie d'un squelette vide.
  L'ancien code se consulte via un worktree : `git worktree add ../TCNtools-legacy v0.1.0`
  (ou `git show v0.1.0:R/<fichier>.R`). Ne pas recopier l'ancien code en bloc :
  chaque fonction est réintégrée, revue et testée individuellement.
- Une branche par couche (`couche-NN-<nom>`), fusionnée dans `v1` par pull request.
- Bascule finale : merge `v1` → `main`, tag `v1.0.0`, release GitHub (DOI Zenodo),
  puis passage de GitHub Pages de `main /docs` à `gh-pages`, et mise à jour de TCNbook.

## Problèmes identifiés dans v0.1.0 (à ne pas reproduire)

- Aucun test.
- `exportPattern("^[[:alpha:]]+")` : tout était exporté, y compris les fonctions internes
  (`likelihood1`, `prior1`, `proposal1`, `run_mcmc1`, `fun_solv_ss_eros_eul`, `LZ`, `E2R`, `Rv0`, `phi_v_slhl`, `d2r`, `r2d`...).
- Noms d'arguments courts et incohérents (`p`, `S`, `L`, `prm`, `Lambda`, `h`, `w` ; `ero` vs `E` ; `rho` vs `rhob`), unités implicites.
- `docs/` (site pkgdown construit) commité dans le dépôt ; `dev/`, `install.R`, `runtime.txt`, `index.md` absents de `.Rbuildignore`.
- `Depends: R (>= 2.10)`, RoxygenNote obsolète.
- `akima` dans Imports sans appel visible : vérifier ; si une interpolation 2D est nécessaire, utiliser `interp` (akima a une licence ACM non commerciale, incompatible en esprit avec GPL-3).
- Données internes (`sysdata.rda`) construites depuis des `.mat` dans `dev/data/` : migrer vers `data-raw/` avec sources documentées.
- TODO hérités (`dev/notes.md`) : vectorisation de `solv_conc_eul`, cas E = 0 dans `conc_soil_mixing`, facteur d'enrichissement en quartz.

## Conventions pour la v1

- Exports explicites via roxygen (`@export`) ; les fonctions internes sont marquées `@noRd` ou `@keywords internal`.
- Unités fixées et rappelées dans chaque `@param` : profondeur massique en g/cm², temps en années,
  concentrations en at/g, taux de dénudation en g/cm²/an (préciser si cm/an est accepté).
- Nommage : snake_case, familles préfixées (`prod_*`, `scaling_*`, `conc_*`, `solve_*`, `tnp_*`),
  noms d'arguments complets, même nom pour la même grandeur partout.
- Paramètres des nucléides : un constructeur validé (ex. `tcn_nuclide("Be10")`) renvoyant une liste S3 légère
  avec provenance (taux de production de référence, demi-vie, longueurs d'atténuation, source bibliographique).
- Contrat de vectorisation documenté pour chaque fonction (arguments recyclés vs scalaires).
- Validation des entrées en base R (`stopifnot` ou messages explicites) ; dépendances minimales.
- Références via Rdpack (`inst/REFERENCES.bib`, récupéré de `inst/TCNtools.bib` du tag v0.1.0) et `\insertRef{}`.
- Documentation et code en anglais.

## Définition de « fait » pour une fonction

1. Documentation roxygen complète : description, unités, `@references`, `@return`.
2. Exemple exécutable.
3. Tests testthat (3e édition) : valeurs de référence publiées ou issues de calculateurs indépendants
   (CRONUS-Earth / calculateur de Balco, CREp), limites analytiques, tests aller-retour pour les inversions.
4. `devtools::check()` sans erreur ni warning.
5. Entrée dans `NEWS.md`.

## Feuille de route (couches, dans l'ordre)

| Couche | Contenu | Validation |
|---|---|---|
| 0. Infrastructure | DESCRIPTION, NAMESPACE roxygen, testthat 3, Rdpack, `data-raw/`, `.Rbuildignore`, `NEWS.md`, `CITATION.cff`, GitHub Actions (`check-standard`, `test-coverage`, pkgdown → gh-pages) | check passe à vide |
| 1. Constantes et données | constantes de désintégration, longueurs d'atténuation, taux de production SLHL ; scripts `data-raw/` depuis les `.mat` | valeurs des articles cités |
| 2. Atmosphère / géométrie | `atm_pressure` (ERA40 / Stone 2000), utilitaires d'angles (internes) | valeurs de référence publiées |
| 3. Géomagnétisme | séries de VDM, calcul de Rc | valeurs ponctuelles vs CRONUS-Earth / CREp |
| 4. Scaling | St, Lm, puis LSDn | grille de sites vs calculateur en ligne / CREp |
| 5. Production en profondeur | spallation exponentielle ; muons (Balco 2017, modèles 1a/1b) | tables de Balco 2017 ; identités intégrales |
| 6. Modèles directs de concentration | analytique (stationnaire et transitoire), puis numérique eulérien / lagrangien | limites : E = 0 → saturation, t → ∞, décroissance pure ; numérique vs analytique |
| 7. Inversion | taux de dénudation stationnaire et âge d'exposition avec incertitudes (recherche de racine) | aller-retour direct → inverse |
| 8. Deux nucléides | `tnp_curves`, `tnp_ellipse`, enfouissement | cohérence avec la couche 6 |

Reportés en v1.x (une issue GitHub chacune, label `v1.x`) : solveur MCMC, `depth_profile_mc`,
famille mélange des sols, éventuellement les flux de Sato s'ils ne servent qu'en intermédiaire pour LSD.

Les scripts `dev/test_*.R` du tag v0.1.0 peuvent servir de point de départ pour les tests.

## Cycle de travail par couche

```bash
git switch v1 && git switch -c couche-NN-<nom>
# écrire R/, tests/testthat/, exemples, vignette si pertinent
# devtools::document(); devtools::test(); devtools::check()
git push -u origin couche-NN-<nom>   # puis PR vers v1
```

## Avancement

- [ ] Tag `v0.1.0` posé et poussé ; TCNbook / TP / Binder pointent dessus
- [ ] Branche `v1` créée, squelette vide, worktree legacy
- [ ] Couche 0 — infrastructure
- [ ] Couche 1 — constantes et données
- [ ] Couche 2 — atmosphère / géométrie
- [ ] Couche 3 — géomagnétisme
- [ ] Couche 4 — scaling
- [ ] Couche 5 — production en profondeur
- [ ] Couche 6 — modèles directs
- [ ] Couche 7 — inversion
- [ ] Couche 8 — deux nucléides
- [ ] Bascule `v1` → `main`, `v1.0.0`, DOI Zenodo, GitHub Pages, TCNbook

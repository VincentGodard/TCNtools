#!/usr/bin/env bash
# Build the pkgdown website into docs/ from the current sources.
#
# - pkgdown renders every .md file at the package root as a page, with no way
#   to exclude one. The build therefore runs on a temporary copy holding only
#   the files git knows about (tracked, or untracked but not ignored), minus
#   notes/ and CLAUDE.md, so that working notes never reach the website.
# - Rdpack inserts the references from the REFERENCES.bib of the *installed*
#   TCNtools. To leave the copy in the user library untouched (v0.1.0, still
#   used by TCNbook and the student lab), the sources are installed in a
#   private library, put first on the library path for the build only.
#
# Usage: pkgdown/build-site.sh
set -euo pipefail
cd "$(dirname "$0")/.."
repo="$PWD"
lib="${TCNTOOLS_SITE_LIB:-$HOME/.cache/R-TCNtools-site-lib}"
mkdir -p "$lib"

src="$(mktemp -d)"
trap 'rm -rf "$src"' EXIT
git ls-files -z --cached --others --exclude-standard | while IFS= read -r -d '' f; do
  case "$f" in notes/* | CLAUDE.md) continue ;; esac
  [ -e "$f" ] && cp --parents "$f" "$src"
done

R CMD INSTALL --no-test-load --no-docs -l "$lib" "$src" > "$lib/install.log" 2>&1 \
  || { echo "install failed, see $lib/install.log"; exit 1; }
export R_LIBS="$lib${R_LIBS:+:$R_LIBS}"

rm -rf "$repo/docs"
Rscript -e "pkg <- pkgdown::as_pkgdown('$src', override = list(destination = '$repo/docs'))" \
        -e "pkgdown::build_site(pkg, install = FALSE, new_process = FALSE, preview = FALSE)"

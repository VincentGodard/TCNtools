"""Convert the ERA-40 climatology grid from ERA40.mat to era40_grid.csv.

One-off conversion, kept for traceability. R cannot read .mat files without
an extra package, so the grid is stored as plain text and data-raw/era40.R
reads the CSV.

Source file: common/production/ERA40.mat of the CRONUScalc code,
https://bitbucket.org/cronusearth/cronus-calc (commit cdbc8f2, md5
bd9f751f13f05bb1edf743bca7078ba8, MATLAB 5.0 MAT-file created 2011-03-31),
byte-identical to dev/data/data_internal/ERA40.mat at tag v0.1.0 of
TCNtools. It holds the long-term mean sea
level pressure (hPa) and mean 1000 hPa temperature (K) from the ERA-40
reanalysis (Uppala et al. 2005) on a 2.5 degree grid, as used by the
CRONUS-Earth calculators (Marrero et al. 2016) and CREp (Martin et al. 2017).

Transformations:
- the longitude column at 360 degrees duplicates the one at 0 degrees (they
  differ by less than 1e-11) and is dropped: data-raw/era40.R adds it back;
- the grid is written in long format, latitude from -90 to 90, longitude
  from 0 to 357.5, with 12 significant digits.

Usage, from the package root (needs scipy):
    git show v0.1.0:dev/data/data_internal/ERA40.mat > /tmp/ERA40.mat
    python data-raw/era40/era40_from_mat.py /tmp/ERA40.mat data-raw/era40/era40_grid.csv
"""

import csv
import sys

import numpy as np
import scipy.io as sio


def main(mat_path, csv_path):
    m = sio.loadmat(mat_path)
    lat = m["ERA40lat"].ravel()
    lon = m["ERA40lon"].ravel()
    slp = m["meanP"]
    t1000 = m["meanT"]

    assert slp.shape == t1000.shape == (lat.size, lon.size)
    assert np.allclose(np.diff(lat), -2.5) and np.allclose(np.diff(lon), 2.5)
    assert lon[0] == 0 and lon[-1] == 360
    assert np.allclose(slp[:, 0], slp[:, -1], rtol=0, atol=1e-9)
    assert np.allclose(t1000[:, 0], t1000[:, -1], rtol=0, atol=1e-9)

    # latitude south to north, drop the 360 degree column
    order = np.argsort(lat)
    keep = slice(0, lon.size - 1)
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["latitude", "longitude", "slp_hpa", "t1000_k"])
        for i in order:
            for j in range(lon.size)[keep]:
                w.writerow([f"{lat[i]:g}", f"{lon[j]:g}",
                            f"{slp[i, j]:.12g}", f"{t1000[i, j]:.12g}"])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

"""Convert the geomagnetic history of the Lm scheme from pmag_consts.mat to CSV.

One-off conversion, kept for traceability. R cannot read .mat files without
an extra package, so the tables are stored as plain text and
data-raw/geomagnetic.R reads them.

Source file: common/production/pmag_consts.mat of the CRONUScalc code,
https://bitbucket.org/cronusearth/cronus-calc (commit cdbc8f2, md5
db7e74b494b2ee6f40127ce3622ca72a, MATLAB 5.0 MAT-file created 2012-09-14).
Its geomagnetic fields are identical to dev/data/data_internal/PMag_Sep12.mat
at tag v0.1.0 of TCNtools. It holds the geomagnetic framework of Lifton et
al. (2014), used by CRONUScalc (get_tdsf.m) for the Lm scaling scheme.

Only the fields used by the Lm scheme are kept:
- t_fineRc, MM0_KCL, lat_pp_KCL, lon_pp_KCL: 0 to 6960 yr, dipole moment
  ratio M/M0 and position of the geomagnetic north pole, from the DGRF
  (2010 to 1950), CALS3k.3 (0 to 3 ka) and CALS7k.2 (3 to 7 ka);
- t_M, M: 7060 yr to 10 Myr, dipole moment ratio M/M0 from GLOPIS-75
  (7 to 18 ka) and PADM2M (18 ka to 2 Ma), then a constant value.
Times are in years before 2010 (t = 0 is 2010 CE), M0 is the 2010 DGRF
dipole moment. The gridded cutoff rigidities (TTRc, IHRc), solar
modulation (S, SPhi) and cross sections are LSD inputs (layer 4b).

Transformations: none on the values (written with 12 significant digits).

Usage, from the package root (needs scipy):
    python data-raw/geomagnetic/geomag_from_mat.py /path/to/pmag_consts.mat data-raw/geomagnetic
"""

import csv
import pathlib
import sys

import numpy as np
import scipy.io as sio


def write(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main(mat_path, out_dir):
    out = pathlib.Path(out_dir)
    p = sio.loadmat(mat_path, squeeze_me=True, struct_as_record=False)["pmag_consts"]

    t_pole = np.asarray(p.t_fineRc, dtype=float)
    assert t_pole.size == p.MM0_KCL.size == p.lat_pp_KCL.size == p.lon_pp_KCL.size
    assert np.all(np.diff(t_pole) > 0) and t_pole[0] == 0
    write(out / "geomag_poles.csv",
          ["time", "moment_ratio", "pole_latitude", "pole_longitude"],
          [[f"{t:.10g}", f"{m:.12g}", f"{la:.12g}", f"{lo:.12g}"]
           for t, m, la, lo in zip(t_pole, p.MM0_KCL, p.lat_pp_KCL, p.lon_pp_KCL)])

    t_m = np.asarray(p.t_M, dtype=float)
    assert t_m.size == p.M.size and np.all(np.diff(t_m) > 0)
    assert t_m[0] > t_pole[-1]
    write(out / "geomag_moment.csv",
          ["time", "moment_ratio"],
          [[f"{t:.10g}", f"{m:.12g}"] for t, m in zip(t_m, p.M)])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

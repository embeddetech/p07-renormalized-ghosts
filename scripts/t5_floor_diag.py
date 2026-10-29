"""Floor diagnostic: is the lambda=0 drift Trotter or truncation?
Usage: python t5_floor_diag.py <tag> <chi> <dt>
Truncation-dominated => floor grows with step count (finer dt worse,
coarser better; larger chi better). Trotter-dominated => floor ~ dt^2.
"""

import sys
import time

import numpy as np
import t5_mps as m


def main():
    tag, chi, dt = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    ops = m.site_ops()
    t0 = time.time()
    ts, dE = m.run_case(ops, 16, -1.0, 0.0, chi, dt=dt)
    with open("dt01_%s.txt" % tag, "w") as fh:
        fh.write("floor diag L=16 lam=0 chi=%d dt=%.3f\n" % (chi, dt))
        fh.write("maxabs_dE = %.6e\n" % np.max(np.abs(dE)))
        fh.write("dE_end_per_site = %.6e\n" % (dE[-1] / 8))
        fh.write("walltime_s = %.0f\n" % (time.time() - t0))
        fh.write("done\n")


if __name__ == "__main__":
    main()

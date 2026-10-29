"""dt=0.01 production ladder, one case per process (the lab server).
Usage: python t5_prod_dt01.py <case>
Writes dt01_<case>.txt when done. Full-throttle authorized by user
2026-07-12 (Normal priority, ~3 BLAS threads per process, 8-way).
"""

import sys
import time

import numpy as np
import t5_mps as m


CASES = {
    #  name     L   sign  lam    chi
    "floor":   (16, -1.0, 0.00,  40),
    "twin":    (16, +1.0, 0.06,  40),
    "lam003":  (16, -1.0, 0.03,  40),
    "lam006":  (16, -1.0, 0.06,  40),
    "lam012":  (16, -1.0, 0.12,  40),
    "L24":     (24, -1.0, 0.06,  40),
    "L32":     (32, -1.0, 0.06,  40),
    "chi64":   (16, -1.0, 0.06,  64),
}


def main():
    cid = sys.argv[1]
    L, sgn, lam, chi = CASES[cid]
    ops = m.site_ops()
    t0 = time.time()
    ts, dE = m.run_case(ops, L, sgn, lam, chi, dt=0.01)
    interior = (3 * L // 4) - (L // 4)
    rate = np.polyfit(ts[len(ts) // 3:],
                      dE[len(ts) // 3:], 1)[0] / interior
    with open("dt01_%s.txt" % cid, "w") as fh:
        fh.write("case %s L=%d sign=%+.0f lam=%.2f chi=%d dt=0.01\n"
                 % (cid, L, sgn, lam, chi))
        fh.write("maxabs_dE = %.6e\n" % np.max(np.abs(dE)))
        fh.write("dE_end_per_site = %.6e\n" % (dE[-1] / interior))
        fh.write("rate_per_site = %.6e\n" % rate)
        fh.write("walltime_s = %.0f\n" % (time.time() - t0))
        fh.write("done\n")


if __name__ == "__main__":
    main()

"""General one-off MPS run. Usage:
python t5_run_one.py <tag> <L> <sign> <lam> <chi> <dt>
Writes diag_<tag>.txt.
"""
import sys
import time
import numpy as np
import t5_mps as m


def main():
    tag = sys.argv[1]
    L, sign, lam = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    chi, dt = int(sys.argv[5]), float(sys.argv[6])
    ops = m.site_ops()
    t0 = time.time()
    ts, dE = m.run_case(ops, L, sign, lam, chi, dt=dt)
    interior = (3 * L // 4) - (L // 4)
    with open("diag_%s.txt" % tag, "w") as fh:
        fh.write("tag %s L=%d sign=%+.0f lam=%.2f chi=%d dt=%.3f\n"
                 % (tag, L, sign, lam, chi, dt))
        fh.write("maxabs_dE = %.6e\n" % np.max(np.abs(dE)))
        fh.write("dE_end_per_site = %.6e\n" % (dE[-1] / interior))
        fh.write("walltime_s = %.0f\n" % (time.time() - t0))
        fh.write("done\n")


if __name__ == "__main__":
    main()

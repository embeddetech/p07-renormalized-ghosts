"""Production MPS runs (post-gate): E1 lambda=0 control, E3 normal
twin, E2 volume sweep L=16/24/32, E4 coupling scaling, plus a
bond-dimension spot check. Certified settings: 2nd-order Trotter
dt=0.02, annealed vacuum prep, chi=40 (chi=64 spot check).
Appends to mps_prod_results.txt line by line."""

import time
import numpy as np
import t5_mps as m


def log(fh, msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def main():
    ops = m.site_ops()
    with open("mps_prod_results.txt", "a") as fh:
        log(fh, "prod start (chi=%d nmax=%d lam=%.2f)"
            % (m.CHI, m.NMAX, m.LAM))

        # E1: lambda=0 control (error floor)
        ts, dE = m.run_case(ops, 16, -1.0, 0.0, m.CHI)
        log(fh, "E1 lambda=0 L=16: max|dE| = %.3e (error floor)"
            % np.max(np.abs(dE)))

        # E3: normal twin
        ts, dEn = m.run_case(ops, 16, +1.0, m.LAM, m.CHI)
        log(fh, "E3 normal twin L=16: max|dE| = %.3e" %
            np.max(np.abs(dEn)))

        # E2: ghost volume sweep
        for L in (16, 24, 32):
            t0 = time.time()
            ts, dEg = m.run_case(ops, L, -1.0, m.LAM, m.CHI)
            interior = (3 * L // 4) - (L // 4)
            rate = np.polyfit(ts[len(ts) // 3:],
                              dEg[len(ts) // 3:], 1)[0] / interior
            log(fh, "E2 L=%d: dE(12)/site = %.4e  rate/site = %.4e  "
                "[%.0fs]" % (L, dEg[-1] / interior, rate,
                             time.time() - t0))

        # E4: coupling scaling at L=16
        rates = {}
        for lam in (0.03, 0.06, 0.12):
            ts, dEg = m.run_case(ops, 16, -1.0, lam, m.CHI)
            r = np.polyfit(ts[len(ts) // 3:],
                           dEg[len(ts) // 3:], 1)[0]
            rates[lam] = r
            log(fh, "E4 lam=%.2f: rate = %.4e" % (lam, r))
        log(fh, "E4 ratios (expect ~4 per doubling): %.2f  %.2f"
            % (rates[0.06] / rates[0.03], rates[0.12] / rates[0.06]))

        # chi spot check
        ts, dEg = m.run_case(ops, 16, -1.0, m.LAM, 64)
        log(fh, "chi=64 spot: dE(12) = %.4e (vs chi=40 above)"
            % dEg[-1])
        log(fh, "prod done")


if __name__ == "__main__":
    main()

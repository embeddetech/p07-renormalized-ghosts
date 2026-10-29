"""
T5: the dimension ladder of the QUANTUM channel sum, measured from
channel tables (couplings analytic, verified against the exact
Hamiltonian's vacuum column at small M).

PRE-REGISTERED PREDICTIONS (written before the first run; note this
CORRECTS the author's own statement of one message earlier, which
conflated the two objects):
  P1. d=1, equal masses, EXACT-resonant shell sum S_res(Lambda):
      SATURATES (the coherent b-component of the stage-4 rate fit;
      sum_k 1/omega^4 converges).
  P2. d=1, golden-rule smoothed channel sum Gamma_FGR(Lambda):
      grows LOGARITHMICALLY -- the near-resonant channel density at
      E ~ 0 is the boost-orbit volume (eta_max ~ ln Lambda), and the
      d=1 reduced measure is UV-finite: Gamma ~ I_s ln Lambda. This
      is Garriga-Vilenkin's D=2 log in channel-table form
      (paper 2's closeout showed the same object with a causal cap).
  P3. d=2: Gamma_FGR(Lambda) grows as a POWER ~ Lambda^2
      (interpolating paper 2's energy ladder 2(d-1): d=1 -> 0=log,
      d=3 -> 4).
Controls: (C1) analytic couplings vs exact H vacuum column at
M = 4, 6 (max relative error); (C2) smoothing-width robustness
(eta = 0.05/0.1/0.2); (C3) volume doubling at fixed cutoff
(extensivity of Gamma per volume).

Channel table: 4-particle states phi(k1,k2) chi(q1,q2), total
momentum zero. Analytic coupling (from H_int = lambda int phi^2
chi^2, box normalization, volume V = L^d):
  g = lambda / (V sqrt(w1 w2 w3 w4)) * (1/sqrt2 per identical pair).
Energy defect (ghost sign): E = w(k1)+w(k2)-w(q1)-w(q2).

ASCII only; numpy (+ t5_channel_ladder import for the verification
control only).
"""

import numpy as np
from collections import defaultdict

LAM = 0.02
MASS = 1.0


def omega_k(kvec, m=MASS):
    return np.sqrt(m * m + np.sum(np.asarray(kvec, float) ** 2))


def mode_list(M, d):
    if d == 1:
        return [(n,) for n in range(-M, M + 1) if n != 0]
    return [(a, b) for a in range(-M, M + 1) for b in range(-M, M + 1)
            if (a, b) != (0, 0)]


def pair_table(M, d, L):
    """Unordered pairs of modes: arrays of (sum-momentum key,
    1/sqrt(w1 w2) coupling factor incl. identical-pair 1/sqrt2,
    energy w1+w2)."""
    mds = mode_list(M, d)
    kfac = 2 * np.pi / L
    w = {n: omega_k([kfac * c for c in n]) for n in mds}
    buckets = defaultdict(lambda: ([], []))   # P -> (coup, energy)
    for i, n1 in enumerate(mds):
        for n2 in mds[i:]:
            P = tuple(a + b for a, b in zip(n1, n2))
            c = 1.0 / np.sqrt(w[n1] * w[n2])
            if n1 == n2:
                c /= np.sqrt(2.0)
            e = w[n1] + w[n2]
            buckets[P][0].append(c)
            buckets[P][1].append(e)
    return {P: (np.array(c), np.array(e))
            for P, (c, e) in buckets.items()}


def channel_sums(M, d, L, etas, tol=1e-9):
    """Return (S_res, {eta: Gamma_FGR}) over all zero-momentum
    4-particle channels."""
    V = L ** d
    pt = pair_table(M, d, L)
    S_res = 0.0
    G = {eta: 0.0 for eta in etas}
    for P, (cphi, ephi) in pt.items():
        mP = tuple(-x for x in P)
        if mP not in pt:
            continue
        cchi, echi = pt[mP]
        # channel coupling g = lam/V * cphi*cchi ; defect E = ephi-echi
        gg = (LAM / V) ** 2 * np.outer(cphi ** 2, cchi ** 2)
        dE = ephi[:, None] - echi[None, :]
        S_res += gg[np.abs(dE) < tol].sum()
        for eta in etas:
            G[eta] += float(np.sum(2 * np.pi * gg *
                                   np.exp(-dE ** 2 / (2 * eta ** 2)) /
                                   (eta * np.sqrt(2 * np.pi))))
    return S_res, G


def verify_couplings():
    """C1: analytic g vs the exact Hamiltonian's vacuum column."""
    import t5_channel_ladder as t5
    print("\n[C1] analytic couplings vs exact vacuum column (d=1):")
    for M in (4, 6):
        t5.L = 8.0
        H, mds, states, index = t5.build_H(M, 4, LAM, ghost=True)
        vac = index[tuple([0] * len(mds))]
        colv = np.asarray(H[:, vac].todense()).ravel()
        kfac = 2 * np.pi / 8.0
        errs = []
        for j, st in enumerate(states):
            if j == vac or abs(colv[j]) < 1e-14:
                continue
            if sum(st) != 4:
                continue                      # 4-particle rows only
            phis = [(f, n) for i, (f, n) in enumerate(mds)
                    for _ in range(st[i]) if f == 0]
            chis = [(f, n) for i, (f, n) in enumerate(mds)
                    for _ in range(st[i]) if f == 1]
            if len(phis) != 2 or len(chis) != 2:
                continue
            ws = [omega_k([kfac * n]) for (f, n) in phis + chis]
            g = LAM / 8.0 / np.sqrt(np.prod(ws))
            if phis[0] == phis[1]:
                g /= np.sqrt(2)
            if chis[0] == chis[1]:
                g /= np.sqrt(2)
            errs.append(abs(abs(colv[j]) - g) / g)
        print("  M=%d : channels checked %d ; max rel err %.2e"
              % (M, len(errs), max(errs)))


def main():
    print("=" * 72)
    print("T5: dimension ladder of the quantum channel sum "
          "(channel tables)")
    print("=" * 72)

    verify_couplings()

    etas = (0.05, 0.1, 0.2)

    # ---------------- d=1 ladder ----------------------------------------
    print("\n[P1+P2] d=1, L=16: S_res (exact shell) and Gamma_FGR vs "
          "cutoff")
    print("  M     Lambda   S_res[e-9]   G(0.05)[e-6]  G(0.1)[e-6]  "
          "G(0.2)[e-6]")
    Ms1 = (8, 16, 32, 64, 128, 256)
    res1 = []
    for M in Ms1:
        S, G = channel_sums(M, 1, 16.0, etas)
        res1.append((M, S, G))
        print("  %3d  %7.1f  %10.4f  %11.4f  %11.4f  %11.4f"
              % (M, 2 * np.pi * M / 16.0, S * 1e9,
                 G[0.05] * 1e6, G[0.1] * 1e6, G[0.2] * 1e6))
    # fits: log vs power for the FGR sum (last 4 points), eta=0.1
    lam_c = np.array([2 * np.pi * M / 16.0 for (M, S, G) in res1])
    g1 = np.array([G[0.1] for (M, S, G) in res1])
    s1 = np.array([S for (M, S, G) in res1])
    loglin = np.polyfit(np.log(lam_c[-4:]), g1[-4:], 1)
    powfit = np.polyfit(np.log(lam_c[-4:]), np.log(g1[-4:]), 1)
    print("  S_res saturation: last-doubling change %.2f%%"
          % (100 * (s1[-1] / s1[-2] - 1)))
    print("  Gamma_FGR: linear-in-ln fit slope %.3e /ln ; power fit "
          "exponent %.3f" % (loglin[0], powfit[0]))
    print("  (P2 verdict: log growth <=> good linear-in-ln fit with "
          "power exponent << 1)")

    # ---------------- d=2 ladder ----------------------------------------
    print("\n[P3] d=2, L=16: Gamma_FGR vs cutoff")
    print("  M    Lambda   G(0.05)[e-6]  G(0.1)[e-6]  G(0.2)[e-6]")
    Ms2 = (2, 3, 4, 6, 8)
    res2 = []
    for M in Ms2:
        S, G = channel_sums(M, 2, 16.0, etas)
        res2.append((M, G))
        print("  %2d  %7.1f  %11.4f  %11.4f  %11.4f"
              % (M, 2 * np.pi * M / 16.0,
                 G[0.05] * 1e6, G[0.1] * 1e6, G[0.2] * 1e6))
    lam2 = np.array([2 * np.pi * M / 16.0 for (M, G) in res2])
    g2 = np.array([G[0.1] for (M, G) in res2])
    p2 = np.polyfit(np.log(lam2[-3:]), np.log(g2[-3:]), 1)[0]
    print("  power-fit exponent (last 3): %.3f  (P3 prediction ~ 2)"
          % p2)

    # ---------------- [C3] volume doubling ------------------------------
    print("\n[C3] volume doubling at fixed cutoff (d=1, Lambda ~ "
          "12.6):")
    for (M, Lbox) in ((32, 16.0), (64, 32.0)):
        S, G = channel_sums(M, 1, Lbox, (0.1,))
        print("  M=%3d L=%4.1f : Gamma_FGR/L = %.4e"
              % (M, Lbox, G[0.1] / Lbox))
    print("  (ratio ~ 1 => extensive rate)")

    print("\ndone.")


if __name__ == '__main__':
    main()

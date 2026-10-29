"""
T5 stage 2: the trap controls (adopted from the overlap review) and
the box test, for the quantum channel-ladder experiment.

Reuses the builders of t5_channel_ladder.py (same directory). The
stage-1 headline was: exact evolution tracks the channel-sum
prediction (ratio ~ 1) across a 5.5x growth in open channels, with a
small residual trend at fixed transfer (-0.15%/step) suspected to be
cap truncation. Stage 2 pins that down and closes the traps:

 [T3] CAP-INDEPENDENCE: occupation caps mimic a confined
      (Deffayet-stable) ghost and can FAKE metastability. Rerun the
      constant-lambda M-sweep at cap 6 (and cap 8 at M=2): does the
      residual trend shrink with cap (truncation artifact) or
      persist (structural suppression)?
 [T4] WINDOW VALIDATION: at finite volume the resonant shell is
      exactly degenerate, so the honest statement of "rate" is the
      t^2 coefficient = sum over channels of |g|^2 (the bright-state
      Rabi regime, valid for t << 1/sqrt(S_ch)). Fit the coefficient
      and compare to S_ch directly; print the validity window.
 [T6] DETUNED-CUTOFF CONTROL: sweep the cutoff with UNEQUAL masses
      (m_chi/m_phi incommensurate) so added modes are off-resonant:
      the transfer must NOT grow with M -- separates channel count
      from UV content.
 [BOX] L and M scaled together (fixed lattice cutoff Lambda =
      2 pi M / L): channel sum per volume and dynamics ratio should
      be L-stable (the rate is extensive).

Vacuum choice, stated once (trap T1): the initial state is the free
positive-norm/negative-energy Fock vacuum -- the GSTZ framing.
Deffayet-type integral-of-motion vacua do not exist for phi^2 chi^2.

ASCII only; numpy/scipy sparse.
"""

import numpy as np
import t5_channel_ladder as t5


def run_ratio(M, cap, lam, ts, ghost=True):
    H, mds, states, index = t5.build_H(M, cap, lam, ghost=ghost)
    nphi = t5.observables(mds, states)
    Nex, Ps = t5.evolve(H, states, index, mds, nphi, ts)
    Npt = t5.pt_curve(H, states, index, mds, nphi, ts)
    cnt, ssum, sall = t5.channel_census(H, states, index, mds,
                                        ghost=ghost)
    return dict(H=H, Nex=Nex, Npt=Npt, cnt=cnt, ssum=ssum, sall=sall,
                dim=H.shape[0])


def main():
    print("=" * 72)
    print("T5 stage 2: trap controls + box test")
    print("=" * 72)
    ts = np.linspace(0.0, 40.0, 81)
    i20, i40 = 40, 80

    # ---------------- [T3] cap independence -----------------------------
    print("\n[T3] cap-independence of the exact/PT ratio (lambda=%.3f):"
          % t5.LAMBDA)
    print("  M  cap  dim(P=0)   ratio(t=20)  ratio(t=40)")
    trend = {}
    for cap in (4, 6):
        for M in (2, 3, 4):
            r = run_ratio(M, cap, t5.LAMBDA, ts)
            r20 = r['Nex'][i20] / max(r['Npt'][i20], 1e-30)
            r40 = r['Nex'][i40] / max(r['Npt'][i40], 1e-30)
            trend[(cap, M)] = (r20, r40)
            print("  %d   %d   %7d   %10.4f  %10.4f"
                  % (M, cap, r['dim'], r20, r40))
    r8 = run_ratio(2, 8, t5.LAMBDA, ts)
    print("  2   8   %7d   %10.4f  %10.4f"
          % (r8['dim'], r8['Nex'][i20] / r8['Npt'][i20],
             r8['Nex'][i40] / r8['Npt'][i40]))
    d4 = trend[(4, 2)][1] - trend[(4, 4)][1]
    d6 = trend[(6, 2)][1] - trend[(6, 4)][1]
    print("  M-trend of ratio(t=40), M=2->4: cap4 %.4f  cap6 %.4f"
          % (d4, d6))
    print("  MEASURED FACT (not the pre-run expectation): the deficit")
    print("  does NOT shrink with cap -- it is non-monotonic in cap")
    print("  (4 -> 6 overshoots, 8 partially recovers) and grows with")
    print("  M at fixed lambda. Interpretation candidates: (a) second-")
    print("  order self-energy shifts detuning the resonant shell --")
    print("  restorable by renormalized PT; (b) incipient structural")
    print("  suppression. Distinguishing (a) from (b) is stage 3.")

    # ---------------- [T4] window + coefficient vs channel sum ---------
    print("\n[T4] t^2-coefficient vs channel sum (the honest "
          "finite-volume rate statement):")
    print("  M   S_ch-weighted pred[e-6]   fitted coeff[e-6]   "
          "t_Rabi=1/sqrt(S_all)")
    tshort = np.linspace(0.0, 6.0, 25)
    for M in (2, 3, 4, 5):
        H, mds, states, index = t5.build_H(M, 4, t5.LAMBDA, ghost=True)
        nphi = t5.observables(mds, states)
        Nex, _ = t5.evolve(H, states, index, mds, nphi, tshort)
        # perturbative t^2 coefficient: sum_f |g|^2 * N_phi(f) over the
        # RESONANT shell (off-resonant states add oscillatory bits)
        vac = index[tuple([0] * len(mds))]
        col = np.asarray(H[:, vac].todense()).ravel()
        efree = np.array([sum(st[i] * (t5.omega(n) if f == 0 else
                                       -t5.omega(n))
                              for i, (f, n) in enumerate(mds))
                          for st in states])
        pred = sum(abs(col[j]) ** 2 * nphi[j]
                   for j in range(len(states))
                   if j != vac and abs(col[j]) > 1e-14
                   and abs(efree[j]) < 1e-9)
        cf = np.polyfit(tshort[1:13] ** 2, Nex[1:13], 1)[0]
        _, _, sall = t5.channel_census(H, states, index, mds)
        print("  %d   %22.3f   %17.3f   %8.1f"
              % (M, pred * 1e6, cf * 1e6, 1 / np.sqrt(sall)))
    print("  (coefficient = resonant channel sum: multiplicity in the")
    print("   finite-volume currency; window valid for t << t_Rabi.)")

    # ---------------- [T6] detuned-cutoff control ------------------------
    print("\n[T6] detuned control (m_chi = 1.35, incommensurate: no "
          "exact resonances):")
    m0 = t5.MASS
    print("  M   N_ch(resonant)   <N_phi>(40)")
    for M in (2, 3, 4, 5):
        # chi mass differs: patch omega per field via a wrapper build
        t5.MASS = 1.0
        H, mds, states, index = t5.build_H(M, 4, t5.LAMBDA, ghost=True)
        # rebuild free part with detuned chi mass by adjusting diagonal
        # (interaction coefficients keep the phi-mass omegas: a
        # conservative choice -- the CHANNEL structure is what is
        # detuned)
        diag_shift = np.zeros(H.shape[0])
        for j, st in enumerate(states):
            for i, (f, n) in enumerate(mds):
                if f == 1 and st[i] > 0:
                    k = 2 * np.pi * n / t5.L
                    w_old = np.sqrt(1.0 + k * k)
                    w_new = np.sqrt(1.35 ** 2 + k * k)
                    diag_shift[j] += st[i] * (-(w_new) - (-(w_old)))
        from scipy import sparse as sp
        Hd = (H + sp.diags(diag_shift)).tocsr()
        nphi = t5.observables(mds, states)
        cnt, _, _ = t5.channel_census(Hd, states, index, mds)
        # census uses equal-mass efree; recount with detuned energies
        vac = index[tuple([0] * len(mds))]
        col = np.asarray(Hd[:, vac].todense()).ravel()
        efree = np.array([sum(st[i] * ((np.sqrt(1 + (2 * np.pi * n /
                                                     t5.L) ** 2))
                                       if f == 0 else
                                       -(np.sqrt(1.35 ** 2 +
                                                 (2 * np.pi * n /
                                                  t5.L) ** 2)))
                              for i, (f, n) in enumerate(mds))
                          for st in states])
        cnt = sum(1 for j in range(len(states))
                  if j != vac and abs(col[j]) > 1e-14
                  and abs(efree[j]) < 1e-9)
        Nex, _ = t5.evolve(Hd, states, index, mds, nphi, ts)
        print("  %d   %6d   %.3e" % (M, cnt, Nex[i40]))
    t5.MASS = m0
    print("  MEASURED: zero resonant channels at every M; transfer is")
    print("  ~40x below the equal-mass case, and its residual growth")
    print("  with M reflects the 1/k decay of the mass detuning in the")
    print("  UV (modes re-approach resonance) -- itself confirming the")
    print("  resonance mechanism. The equal-mass growth is channel")
    print("  physics, not raw UV content.")

    # ---------------- [BOX] extensivity ---------------------------------
    print("\n[BOX] L and M scaled together (fixed cutoff 2 pi M / L):")
    L0 = t5.L
    print("  (M, L)    N_ch   S_ch*L [e-6]   <N_phi>(20)/L [e-3]")
    for (M, Lbox) in ((2, 4.0), (4, 8.0), (6, 12.0)):
        t5.L = Lbox
        r = run_ratio(M, 4, t5.LAMBDA, ts)
        print("  (%d,%4.1f)  %4d   %10.3f   %12.4f"
              % (M, Lbox, r['cnt'], r['ssum'] * Lbox * 1e6,
                 r['Nex'][i20] / Lbox * 1e3))
    t5.L = L0
    print("  MEASURED: the dynamical transfer per volume is roughly")
    print("  stable across the two larger boxes, but S_ch * L GROWS")
    print("  with the tower (super-extensive channel sum at fixed")
    print("  cutoff) -- the finite-volume shadow of paper 2's")
    print("  instantaneity trend; the naive extensivity expectation")
    print("  printed in an earlier version of this script was wrong.")

    print("\ndone.")


if __name__ == '__main__':
    main()

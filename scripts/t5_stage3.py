"""
T5 stage 3: THE DECISIVE BLOCK -- is the beyond-leading-order deficit
(stage 2, [T3]) renormalizable level-shift physics (multiplicity
survives, verdict a) or incipient structural suppression (the GSTZ
escape, verdict b)?

METHOD (Schur / Feshbach effective Hamiltonian -- exact folding of the
deep Hilbert space): P = span{vacuum, ALL zero-free-energy states}
(the degenerate shell, coupled or not), Q = the rest. The static
effective Hamiltonian at the vacuum energy E0 = 0,

    H_eff = P H P - P H Q (Q H Q - E0)^(-1) Q H P,

is Hermitian and contains the second-order self-energy shifts,
shell-shell induced couplings, and vertex corrections to all orders
in the Q-space folding (static approximation only). Evolving H_eff
(small dense matrix) and comparing the SHELL-population observable
    N_shell(t) = sum_{f in shell} |<f|psi(t)>|^2 N_phi(f)
against the full exact evolution decides the verdict:
  - full  ~ H_eff  ==>  the deficit is self-energy/detuning physics;
    the channel-sum survives at the renormalized level (verdict a);
  - full << H_eff  ==>  suppression beyond few-level folding
    (verdict b) -- the refutation signal.

Also: [ii] lambda-scaling of the naive deficit at fixed (M=3, cap 6)
-- detuning predicts deficit growth consistent with (dE t)^2, dE ~
lambda^2; [iii] cap convergence at M=2 (caps 4/6/8/10) with the
Schur ratio replacing the naive one.

ASCII only; numpy/scipy.
"""

import numpy as np
from scipy.linalg import expm
import t5_channel_ladder as t5

TOL_RES = 1e-9


def free_energies(states, mds):
    return np.array([sum(st[i] * (t5.omega(n) if f == 0 else
                                  -t5.omega(n))
                         for i, (f, n) in enumerate(mds))
                     for st in states])


def shell_split(H, states, index, mds):
    """P = vacuum + all |E_free| < tol states; Q = rest."""
    efree = free_energies(states, mds)
    vac = index[tuple([0] * len(mds))]
    P_idx = [j for j in range(len(states)) if abs(efree[j]) < TOL_RES]
    if vac not in P_idx:
        P_idx = [vac] + P_idx
    Q_idx = [j for j in range(len(states)) if j not in set(P_idx)]
    return vac, P_idx, Q_idx, efree


def schur_heff(H, P_idx, Q_idx, E0=0.0):
    Hd = H.todense()
    Hpp = np.asarray(Hd[np.ix_(P_idx, P_idx)])
    Hpq = np.asarray(Hd[np.ix_(P_idx, Q_idx)])
    Hqq = np.asarray(Hd[np.ix_(Q_idx, Q_idx)])
    A = Hqq - E0 * np.eye(len(Q_idx))
    mineig = np.min(np.abs(np.linalg.eigvalsh(A)))
    Heff = Hpp - Hpq @ np.linalg.solve(A, Hpq.conj().T)
    Heff = (Heff + Heff.conj().T) / 2
    return Heff, mineig


def shell_population(psi_cols, P_idx, nphi, vacpos_in_P):
    """N_shell(t) from full-state columns restricted to the shell."""
    out = []
    for p in psi_cols:
        out.append(sum(abs(p[j]) ** 2 * nphi[j] for j in P_idx))
    return np.array(out)


def run_case(M, cap, lam, ts):
    H, mds, states, index = t5.build_H(M, cap, lam, ghost=True)
    nphi = t5.observables(mds, states)
    vac, P_idx, Q_idx, efree = shell_split(H, states, index, mds)
    # full exact evolution
    from scipy.sparse.linalg import expm_multiply
    psi0 = np.zeros(H.shape[0], dtype=complex)
    psi0[vac] = 1.0
    cols = expm_multiply(-1j * H, psi0, start=ts[0], stop=ts[-1],
                         num=len(ts), endpoint=True)
    N_full = shell_population(cols, P_idx, nphi, None)
    # naive PT (free energies) restricted to the shell
    colv = np.asarray(H[:, vac].todense()).ravel()
    N_naive = np.zeros(len(ts))
    for j in P_idx:
        if j == vac or abs(colv[j]) < 1e-14:
            continue
        for a, t in enumerate(ts):
            N_naive[a] += (abs(colv[j]) * t) ** 2 * nphi[j]
    # Schur effective evolution
    Heff, mineig = schur_heff(H, P_idx, Q_idx, E0=0.0)
    vp = P_idx.index(vac)
    phi0 = np.zeros(len(P_idx), dtype=complex)
    phi0[vp] = 1.0
    nphi_P = np.array([nphi[j] for j in P_idx])
    N_eff = []
    for t in ts:
        ph = expm(-1j * Heff * t) @ phi0
        N_eff.append(np.real(np.sum(np.abs(ph) ** 2 * nphi_P)))
    N_eff = np.array(N_eff)
    return dict(N_full=N_full, N_eff=N_eff, N_naive=N_naive,
                nP=len(P_idx), nQ=len(Q_idx), mineig=mineig,
                dim=H.shape[0])


def main():
    print("=" * 72)
    print("T5 stage 3: Schur-renormalized channel sum vs full quantum "
          "evolution")
    print("=" * 72)
    ts = np.linspace(0.0, 40.0, 41)
    i40 = 40

    # ---------------- [i] the verdict table -----------------------------
    print("\n[i] shell population at t=40: full exact vs Schur-effective")
    print("    vs naive channel sum   (lambda = %.3f)" % t5.LAMBDA)
    print("  M  cap  dim   |P|  minQgap   full/naive  full/Schur")
    for cap in (4, 6):
        for M in (2, 3, 4):
            r = run_case(M, cap, t5.LAMBDA, ts)
            rn = r['N_full'][i40] / max(r['N_naive'][i40], 1e-30)
            rs = r['N_full'][i40] / max(r['N_eff'][i40], 1e-30)
            print("  %d   %d  %5d  %4d  %7.3f   %9.4f  %9.4f"
                  % (M, cap, r['dim'], r['nP'], r['mineig'], rn, rs))
    print("  VERDICT KEY: full/Schur ~ 1 with full/naive < 1 =>")
    print("  the deficit is folded-space self-energy physics (verdict")
    print("  a: renormalized multiplicity); full/Schur << 1 => beyond-")
    print("  few-level suppression (verdict b: the refutation signal).")

    # ---------------- [ii] lambda scaling of the naive deficit ----------
    print("\n[ii] lambda-scaling of the naive deficit (M=3, cap 6, "
          "t=40):")
    print("  lambda   1 - full/naive   1 - full/Schur")
    for lam in (0.005, 0.01, 0.02, 0.04):
        r = run_case(3, 6, lam, ts)
        dn = 1 - r['N_full'][i40] / max(r['N_naive'][i40], 1e-30)
        ds = 1 - r['N_full'][i40] / max(r['N_eff'][i40], 1e-30)
        print("  %.3f   %13.4f   %13.4f" % (lam, dn, ds))
    print("  (detuning physics: naive deficit grows ~lambda^2 at small")
    print("   deficit -- the shell is exactly degenerate, so the")
    print("   splitting is FIRST order, dE ~ lambda, entering as (dE t)^2;")
    print("   the Schur column should stay ~0 at every lambda if (a).)")

    # ---------------- [iii] cap convergence at M=2 ----------------------
    print("\n[iii] cap convergence, M=2 (Schur ratio should be cap-"
          "stable if (a)):")
    print("  cap  dim    full/naive  full/Schur")
    for cap in (4, 6, 8, 10):
        r = run_case(2, cap, t5.LAMBDA, ts)
        rn = r['N_full'][i40] / max(r['N_naive'][i40], 1e-30)
        rs = r['N_full'][i40] / max(r['N_eff'][i40], 1e-30)
        print("  %2d  %6d  %9.4f  %9.4f" % (cap, r['dim'], rn, rs))

    print("\ndone.")


if __name__ == '__main__':
    main()

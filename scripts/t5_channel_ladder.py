"""
T5: the quantum channel-ladder test -- exact (Krylov) real-time
evolution of a genuinely quantum 1+1d phi^2 chi^2 ghost theory on
momentum modes, vs the channel-sum (multiplicity) prediction.

THE QUESTION (papers 2/4 load-bearing): does the quantum vacuum decay
rate scale with the number of OPEN channels (multiplicity horn), or is
the channel sum dynamically suppressed (the GSTZ resummed-suppression
escape)? Classical lattices cannot see the channel ladder (paper 2's
decomposition); this is the quantum test. REFUTATION of channel
scaling here = evidence FOR the suppression escape = the program's
central premise fails. That is the point of running it.

MODEL: phi (normal) and chi (GHOST: H0 has -omega b^dag b) on modes
k_n = 2 pi n / L, n in {+-1..+-M} (zero modes dropped), equal masses
(maximal resonance), H_int = lambda int dx phi^2 chi^2 expanded in
modes with momentum conservation; full Fock basis with total quanta
<= cap in the total-momentum-0 sector. Ghost sign makes 4-particle
states (phi pair + chi pair) DEGENERATE with the vacuum: E_f =
2 omega(k) - 2 omega(q) = 0 at |k| = |q| -- the open channels.

OBSERVABLES:
 [1] channel census vs M: number of exactly-resonant 4-particle
     states and their coupling-weighted sum S_ch = sum |<f|H|0>|^2.
 [2] exact <N_phi>(t) (sparse Krylov) vs first-order perturbation
     theory N_pt(t) = sum_f |g_f|^2 [2 sin(D_f t/2)/D_f]^2 N_phi(f):
     the ratio exact/PT at the channel-sum level. Multiplicity holds
     iff the ratio stays ~1 as M grows; structural suppression would
     show as the ratio FALLING with channel count.
 [3] coupling sweep at fixed M: where back-reaction sets in (ratio
     vs lambda) -- distinguishes ordinary unitarity saturation from
     channel-structural suppression (the former moves with lambda,
     the latter with M).
 [4] controls: lambda = 0 (nothing happens); NORMAL-sign twin (no
     resonant channels: <N> stays at virtual-pair level, no secular
     growth); cap 4 vs 6 convergence at small M.

ASCII only; numpy/scipy sparse. Runtime: minutes.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from itertools import combinations_with_replacement

# ----------------------------------------------------------------------
# model parameters
# ----------------------------------------------------------------------
MASS = 1.0
L = 8.0
LAMBDA = 0.02
CAP = 4          # total-quanta cap (convergence-checked at 6)


def modes(M):
    """Mode list: (field, n) with field 0=phi, 1=chi; n = +-1..+-M."""
    ns = [n for n in range(-M, M + 1) if n != 0]
    return [(f, n) for f in (0, 1) for n in ns]


def omega(n):
    k = 2 * np.pi * n / L
    return np.sqrt(MASS ** 2 + k ** 2)


# ----------------------------------------------------------------------
# basis: occupation tuples, total quanta <= CAP, total momentum 0
# ----------------------------------------------------------------------

def build_basis(M, cap):
    """FULL basis (all momenta), total quanta <= cap. Operators need
    nonzero-momentum intermediates; the P=0 physical block is
    projected out after H is assembled (H conserves P)."""
    mds = modes(M)
    nm = len(mds)
    states = []
    for s in range(cap + 1):
        for combo in combinations_with_replacement(range(nm), s):
            occ = [0] * nm
            for c in combo:
                occ[c] += 1
            states.append(tuple(occ))
    states = sorted(set(states))
    index = {st: i for i, st in enumerate(states)}
    return mds, states, index


def momentum_of(st, mds):
    return sum(st[i] * n for i, (f, n) in enumerate(mds))


def ladder(mds, states, index, mode_i, dag):
    """Sparse a_i or a_i^dag on the basis (states leaving the basis are
    dropped -- consistent for H preserving quanta<=cap subspace up to
    the cap boundary, which the cap-convergence control quantifies)."""
    rows, cols, vals = [], [], []
    for j, st in enumerate(states):
        n = st[mode_i]
        if dag:
            new = list(st)
            new[mode_i] += 1
            tnew = tuple(new)
            if tnew in index:
                rows.append(index[tnew])
                cols.append(j)
                vals.append(np.sqrt(n + 1))
        else:
            if n > 0:
                new = list(st)
                new[mode_i] -= 1
                rows.append(index[tuple(new)])
                cols.append(j)
                vals.append(np.sqrt(n))
    dim = len(states)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(dim, dim))


def build_H(M, cap, lam, ghost=True):
    mds, states, index = build_basis(M, cap)
    dim = len(states)
    nm = len(mds)
    # free part
    diag = np.zeros(dim)
    for j, st in enumerate(states):
        e = 0.0
        for i, (f, n) in enumerate(mds):
            w = omega(n)
            e += st[i] * (w if (f == 0 or not ghost) else -w)
        diag[j] = e
    H = sparse.diags(diag).tocsr()

    # interaction: lam * int phi^2 chi^2
    # phi(x) = sum_n (a_n e^{ikx} + a_n^dag e^{-ikx})/sqrt(2 L w_n)
    # elementary op e = (mode_i, dag): momentum p(e) = -n if dag else n
    # coefficient factor 1/sqrt(2 L w) each; integral gives L * delta.
    phi_ids = [i for i, (f, n) in enumerate(mds) if f == 0]
    chi_ids = [i for i, (f, n) in enumerate(mds) if f == 1]

    def elem(field_ids):
        out = []
        for i in field_ids:
            n = mds[i][1]
            for dag in (0, 1):
                p = -n if dag else n
                out.append((i, dag, p, 1.0 / np.sqrt(2 * L * omega(n))))
        return out

    ephi, echi = elem(phi_ids), elem(chi_ids)
    Amats = {}

    def op(i, dag):
        key = (i, dag)
        if key not in Amats:
            Amats[key] = ladder(mds, states, index, i, dag)
        return Amats[key]

    # pair products, grouped by total momentum
    def pair_products(elems):
        prods = {}
        for (i1, d1, p1, c1) in elems:
            for (i2, d2, p2, c2) in elems:
                # operator order: first factor acts second (left)
                Ppair = op(i1, d1) @ op(i2, d2)
                key = p1 + p2
                prods.setdefault(key, []).append((c1 * c2, Ppair))
        # sum within each momentum class
        summed = {}
        for key, lst in prods.items():
            Msum = lst[0][0] * lst[0][1]
            for c, P in lst[1:]:
                Msum = Msum + c * P
            summed[key] = Msum
        return summed

    Pphi = pair_products(ephi)
    Pchi = pair_products(echi)
    Hint = None
    for pk, Aphi in Pphi.items():
        if -pk in Pchi:
            term = Aphi @ Pchi[-pk]
            Hint = term if Hint is None else Hint + term
    H = H + lam * L * Hint
    # hermitize (cap boundary truncation makes Hint slightly non-herm)
    H = (H + H.conj().T) / 2
    # project to the physical total-momentum-0 block
    keep = [j for j, st in enumerate(states)
            if momentum_of(st, mds) == 0]
    H = H[keep][:, keep].tocsr()
    states0 = [states[j] for j in keep]
    index0 = {st: i for i, st in enumerate(states0)}
    return H, mds, states0, index0


def observables(mds, states):
    nphi = np.array([sum(st[i] for i, (f, n) in enumerate(mds) if f == 0)
                     for st in states], dtype=float)
    return nphi


def channel_census(H, states, index, mds, tol=1e-9, ghost=True):
    """Resonant 4-particle channels: states f with <f|H|0> != 0 and
    E_f (free) == 0; return count and coupling sum."""
    vac = index[tuple([0] * len(mds))]
    col = np.asarray(H[:, vac].todense()).ravel()
    sgn = -1.0 if ghost else 1.0
    efree = np.array([sum(st[i] * (omega(n) if f == 0 else sgn *
                                   omega(n))
                          for i, (f, n) in enumerate(mds))
                      for st in states])
    cnt, ssum, s_all = 0, 0.0, 0.0
    for j in range(len(states)):
        if j == vac or abs(col[j]) < 1e-14:
            continue
        s_all += abs(col[j]) ** 2
        if abs(efree[j]) < tol:
            cnt += 1
            ssum += abs(col[j]) ** 2
    return cnt, ssum, s_all


def pt_curve(H, states, index, mds, nphi, ts):
    """First-order PT: N_pt(t) = sum_f |g_f|^2 [2 sin(D t/2)/D]^2 N_f."""
    vac = index[tuple([0] * len(mds))]
    col = np.asarray(H[:, vac].todense()).ravel()
    efree = np.array([sum(st[i] * (omega(n) if f == 0 else -omega(n))
                          for i, (f, n) in enumerate(mds))
                      for st in states])
    out = np.zeros(len(ts))
    for j in range(len(states)):
        g = col[j]
        if j == vac or abs(g) < 1e-14:
            continue
        D = efree[j]
        for a, t in enumerate(ts):
            if abs(D) < 1e-12:
                amp2 = (abs(g) * t) ** 2
            else:
                amp2 = abs(g) ** 2 * (2 * np.sin(D * t / 2) / D) ** 2
            out[a] += amp2 * nphi[j]
    return out


def evolve(H, states, index, mds, nphi, ts):
    vac = index[tuple([0] * len(mds))]
    psi0 = np.zeros(H.shape[0], dtype=complex)
    psi0[vac] = 1.0
    res = expm_multiply(-1j * H, psi0, start=ts[0], stop=ts[-1],
                        num=len(ts), endpoint=True)
    Nexp = np.array([np.real(np.vdot(p, nphi * p)) for p in res])
    Psur = np.array([abs(p[vac]) ** 2 for p in res])
    return Nexp, Psur


def main():
    print("=" * 72)
    print("T5: quantum channel-ladder test (ED/Krylov, phi^2 chi^2 "
          "ghost theory)")
    print("=" * 72)
    print("m=%.1f L=%.1f lambda=%.3f cap=%d ; modes +-1..+-M both "
          "fields" % (MASS, L, LAMBDA, CAP))

    ts = np.linspace(0.0, 40.0, 81)

    # ---------------- [1]+[2]: M sweep ----------------------------------
    print("\n[1,2] channel census + exact-vs-PT ratio vs M:")
    print("  M  dim    N_ch  S_ch[e-6]  S_all[e-6]  "
          "ratio(t=20)  ratio(t=40)")
    for M in (2, 3, 4, 5):
        H, mds, states, index = build_H(M, CAP, LAMBDA, ghost=True)
        nphi = observables(mds, states)
        cnt, ssum, sall = channel_census(H, states, index, mds)
        Nex, Ps = evolve(H, states, index, mds, nphi, ts)
        Npt = pt_curve(H, states, index, mds, nphi, ts)
        i20, i40 = 40, 80
        r20 = Nex[i20] / max(Npt[i20], 1e-30)
        r40 = Nex[i40] / max(Npt[i40], 1e-30)
        print("  %d  %5d  %4d  %9.3f  %9.3f  %10.3f  %10.3f"
              % (M, H.shape[0], cnt, ssum * 1e6, sall * 1e6, r20, r40))
    print("  (ratio ~ 1: full quantum evolution tracks the channel-sum")
    print("   prediction -- multiplicity; ratio falling with M would be")
    print("   the suppression signature.)")

    # ---------------- [2b]: constant-transfer M sweep -------------------
    # isolate channel-STRUCTURAL effects from ordinary back-reaction:
    # rescale lambda per M so the total perturbative transfer at t=40
    # is the same for every M; if the residual exact/PT deficit was
    # back-reaction (total-transfer driven), the ratio is now flat in
    # M; a fall persisting here would be channel-structural
    # suppression (the GSTZ escape signature).
    print("\n[2b] constant-transfer M sweep (lambda rescaled so "
          "N_pt(40) is fixed):")
    # reference transfer: M=2 at LAMBDA
    H, mds, states, index = build_H(2, CAP, LAMBDA, ghost=True)
    nphi = observables(mds, states)
    Npt_ref = pt_curve(H, states, index, mds, nphi, ts)[-1]
    for M in (2, 3, 4, 5):
        H, mds, states, index = build_H(M, CAP, LAMBDA, ghost=True)
        nphi = observables(mds, states)
        Npt40 = pt_curve(H, states, index, mds, nphi, ts)[-1]
        lam_M = LAMBDA * np.sqrt(Npt_ref / Npt40)   # N_pt ~ lambda^2
        H, mds, states, index = build_H(M, CAP, lam_M, ghost=True)
        nphi = observables(mds, states)
        Nex, _ = evolve(H, states, index, mds, nphi, ts)
        Npt = pt_curve(H, states, index, mds, nphi, ts)
        print("  M=%d : lambda=%.4f  ratio(t=20) %.4f  ratio(t=40) "
              "%.4f" % (M, lam_M, Nex[40] / max(Npt[40], 1e-30),
                        Nex[80] / max(Npt[80], 1e-30)))

    # ---------------- [3]: coupling sweep at M=4 ------------------------
    print("\n[3] coupling sweep (M=4): exact/PT ratio at t=20,40")
    for lam in (0.005, 0.02, 0.08, 0.32):
        H, mds, states, index = build_H(4, CAP, lam, ghost=True)
        nphi = observables(mds, states)
        Nex, _ = evolve(H, states, index, mds, nphi, ts)
        Npt = pt_curve(H, states, index, mds, nphi, ts)
        print("  lambda=%.3f : ratio %.3f / %.3f" % (
            lam, Nex[40] / max(Npt[40], 1e-30),
            Nex[80] / max(Npt[80], 1e-30)))

    # ---------------- [4]: controls --------------------------------------
    print("\n[4] controls:")
    H0, mds, states, index = build_H(3, CAP, 0.0, ghost=True)
    nphi = observables(mds, states)
    Nex0, _ = evolve(H0, states, index, mds, nphi, ts)
    print("  lambda=0        : max <N_phi> = %.2e (expect 0)"
          % Nex0.max())
    Hn, mds, states, index = build_H(3, CAP, LAMBDA, ghost=False)
    nphi = observables(mds, states)
    cntn, ssn, _ = channel_census(Hn, states, index, mds, ghost=False)
    Nexn, _ = evolve(Hn, states, index, mds, nphi, ts)
    print("  NORMAL twin     : resonant channels = %d; max <N_phi> = "
          "%.2e (virtual only)" % (cntn, Nexn.max()))
    for cap in (4, 6):
        H, mds, states, index = build_H(2, cap, LAMBDA, ghost=True)
        nphi = observables(mds, states)
        Nex, _ = evolve(H, states, index, mds, nphi, ts)
        print("  cap=%d (M=2)    : <N_phi>(40) = %.6e  dim=%d"
              % (cap, Nex[-1], H.shape[0]))

    print("\ndone.")


if __name__ == '__main__':
    main()

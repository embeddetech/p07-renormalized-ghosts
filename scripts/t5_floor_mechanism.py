"""
t5_floor_mechanism.py -- diagnostics behind the accumulated-truncation floor.

Emits the numbers quoted in the paper's floor-mechanism paragraph, which
previously lived only in ad-hoc drafting runs. Each block prints the
CONFIGURATION it was measured at, including WHICH BOND is cut: two of these
are L=4 quantities and one is an L=16 quantity, and at L=4 the middle bond
(1|2, 81 Schmidt values) and the edge bonds (9 Schmidt values) give different
entropies, so "the entropy" is not well defined without saying which cut.

Entropies here are computed by DIRECT DENSE RESHAPING of the state vector, not
by canonicalizing an MPS. That is deliberate: the canonicalization route is
easy to get silently wrong (SVD the left-canonical Q factor and every state
returns log(bond dimension); forget to right-canonicalize the other half and
the singular values are not Schmidt values at all). At L=4 the dense route is
cheap and leaves nothing to get wrong. The L=16 block necessarily uses the MPS.

  [1] L=4, lambda=0: the prepared vacuum is an exact eigenstate of the
      real-time generator (residual ~1e-15), so its exact entanglement is
      STATIONARY -> S = 0.274049 at the middle bond, at every t.

  [2] L=4, lambda=0: the exact vacuum's discarded Schmidt weight beyond
      rank 40 at the middle bond -> ~2.2e-11. Far below the observed floor,
      so the floor is NOT static representability.

  [3] L=4, lambda=0: simulated (truncated) middle-bond entropy, D=40 vs the
      D=100 control -> D=40 drifts ~1e-3 over [0,12]; D=100 flat to ~5e-5.
      The exact state does not move, so at L=4 there is no climb.

  [4] L=16, D=40, dt=0.02 (the PRODUCTION protocol): the simulated state's
      middle-bond entropy CLIMBS, 0.249 at t=0 to ~0.89 at t=12. This is the
      accumulation the floor paragraph refers to. Slow; run with --l16.
      Note the L=16 initial state is itself variational, so its climb mixes
      Trotter injection with imperfect preparation -- which is why the L=4
      D=100 control above is what shows the climb is not intrinsic.

Usage:  python t5_floor_mechanism.py [--l16]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t5_mps as m  # noqa: E402


def mps_to_dense(mps):
    """Contract an MPS down to a state vector (only sane for small L)."""
    out = np.array(mps[0], dtype=complex)
    for t in mps[1:]:
        out = np.tensordot(out, np.array(t, dtype=complex), axes=([-1], [0]))
    return out.reshape(-1)


def bond_spectrum(psi, d2, L, bond):
    """Schmidt spectrum across `bond` (sites 0..bond-1 | bond..L-1), dense."""
    psi = np.asarray(psi, dtype=complex).reshape(-1)
    psi = psi / np.linalg.norm(psi)
    mat = psi.reshape(d2 ** bond, d2 ** (L - bond))
    s = np.linalg.svd(mat, compute_uv=False)
    return s / np.linalg.norm(s)


def entropy_of(s):
    p = s ** 2
    p = p[p > 1e-18]
    return float(-np.sum(p * np.log(p)))


def dense_aux_vacuum(ops, L):
    """The auxiliary ground state, exactly, as t5_mps_gate.py prepares it.

    Built from t5_mps's own term constructors so the model cannot drift from
    the one the rest of the paper uses.
    """
    d2 = ops['d2']
    hs = m.onsite_h(ops, +1.0, 0.0)
    corr = -0.5 * (ops['Xp'] @ ops['Xp'] + ops['Xc'] @ ops['Xc'])
    hb = m.bond_h(ops, +1.0)
    dim = d2 ** L
    H = np.zeros((dim, dim), dtype=complex)
    for i in range(L):
        h_i = hs + (corr if i in (0, L - 1) else 0)
        H += np.kron(np.kron(np.eye(d2 ** i), h_i), np.eye(d2 ** (L - 1 - i)))
    for i in range(L - 1):
        H += np.kron(np.kron(np.eye(d2 ** i), hb), np.eye(d2 ** (L - 2 - i)))
    H = (H + H.conj().T) / 2
    _, V = np.linalg.eigh(H)
    return V[:, 0].astype(complex)


def real_time_dense(ops, L):
    """The lambda=0 real-time generator on the same lattice (ghost sign)."""
    d2 = ops['d2']
    hs = m.onsite_h(ops, -1.0, 0.0)
    corr = -0.5 * (ops['Xp'] @ ops['Xp'] - ops['Xc'] @ ops['Xc'])
    hb = m.bond_h(ops, -1.0)
    dim = d2 ** L
    H = np.zeros((dim, dim), dtype=complex)
    for i in range(L):
        h_i = hs + (corr if i in (0, L - 1) else 0)
        H += np.kron(np.kron(np.eye(d2 ** i), h_i), np.eye(d2 ** (L - 1 - i)))
    for i in range(L - 1):
        H += np.kron(np.kron(np.eye(d2 ** i), hb), np.eye(d2 ** (L - 2 - i)))
    return (H + H.conj().T) / 2


def tebd_entropy_run(ops, mps, L, chi, dt, tmax, d2, bond, label, dense=True):
    gh = m.two_site_gates(ops, -1.0, 0.0, L, dt / 2, imag=False)
    gf = m.two_site_gates(ops, -1.0, 0.0, L, dt, imag=False)
    nst = int(round(tmax / dt))
    every = max(1, int(round(3.0 / dt)))

    def S_now():
        if dense:
            return entropy_of(bond_spectrum(mps_to_dense(mps), d2, L, bond))
        return m_entropy_mps(mps, bond)

    s0 = S_now()
    print("      t= 0.0  S=%.6f   [%s]" % (s0, label))
    last = s0
    for s in range(1, nst + 1):
        m.tebd_sweep2(mps, gh, gf, chi, imag=False)
        if s % every == 0:
            last = S_now()
            print("      t=%4.1f  S=%.6f" % (s * dt, last))
    return s0, last


def m_entropy_mps(mps, bond):
    """Schmidt entropy at `bond` for a large MPS: canonicalize BOTH sides."""
    ts = [np.array(t, dtype=complex) for t in mps]
    L = len(ts)
    for i in range(bond):                      # left-canonicalize 0..bond-1
        dl, d, dr = ts[i].shape
        Q, R = np.linalg.qr(ts[i].reshape(dl * d, dr))
        ts[i] = Q.reshape(dl, d, -1)
        ts[i + 1] = np.tensordot(R, ts[i + 1], axes=(1, 0))
    for i in range(L - 1, bond, -1):           # right-canonicalize L-1..bond+1
        dl, d, dr = ts[i].shape
        Q, R = np.linalg.qr(ts[i].reshape(dl, d * dr).conj().T)
        ts[i] = Q.conj().T.reshape(-1, d, dr)
        ts[i - 1] = np.tensordot(ts[i - 1], R.conj().T, axes=(2, 0))
    dl, d, dr = ts[bond].shape
    s = np.linalg.svd(ts[bond].reshape(dl, d * dr), compute_uv=False)
    return entropy_of(s / np.linalg.norm(s))


def main():
    do_l16 = "--l16" in sys.argv
    ops = m.site_ops()
    d2 = ops['d2']
    L4, BOND4 = 4, 2                    # middle bond of a 4-site chain

    print("=" * 70)
    print("[1] L=4, lambda=0, middle bond (1|2) -- exact eigenstate, S STATIONARY")
    print("=" * 70)
    psi = dense_aux_vacuum(ops, L4)
    H = real_time_dense(ops, L4)
    E = float(np.real(psi.conj() @ (H @ psi)))
    print("    real-time residual ||(H-E)psi|| = %.3e  (exact eigenstate)" % (
        np.linalg.norm(H @ psi - E * psi)))
    print("    entropy at every bond of the exact vacuum:")
    for b in range(1, L4):
        sb = bond_spectrum(psi, d2, L4, b)
        print("      bond %d|%d : S=%.6f   (%d Schmidt values)" % (
            b - 1, b, entropy_of(sb), len(sb)))
    print("    the MIDDLE bond is the one quoted in the paper.")
    # eigenstate => exp(-iHt)psi = exp(-iEt)psi; a global phase cannot change
    # any Schmidt spectrum, and the residual above is what establishes it.
    S0 = entropy_of(bond_spectrum(psi, d2, L4, BOND4))
    for t in (0.0, 3.0, 6.0, 9.0, 12.0):
        St = entropy_of(bond_spectrum(np.exp(-1j * E * t) * psi, d2, L4, BOND4))
        print("      t=%4.1f  S_exact=%.6f" % (t, St))
    print("    => stationary at %.6f" % S0)

    print()
    print("=" * 70)
    print("[2] L=4, lambda=0, middle bond -- discarded Schmidt weight past rank 40")
    print("=" * 70)
    s = bond_spectrum(psi, d2, L4, BOND4)
    print("    Schmidt values at the middle bond : %d" % len(s))
    print("    discarded weight beyond rank 40   : %.4e" % float(np.sum(s[40:] ** 2)))
    print("    => far below the observed floor; NOT static representability")

    print()
    print("=" * 70)
    print("[3] L=4, lambda=0, middle bond -- simulated D=40 vs the D=100 control")
    print("=" * 70)
    for chi in (40, 100):
        mps = m.dense_to_mps(psi.copy(), L4, d2)
        a, b = tebd_entropy_run(ops, mps, L4, chi, 0.02, 12.0, d2, BOND4,
                                "D=%d" % chi, dense=True)
        print("    D=%3d : %.6f -> %.6f   (drift %.3e)" % (chi, a, b, b - a))
    print("    => D=40 drifts ~1e-3; D=100 flat to ~5e-5: the climb is not intrinsic")

    if do_l16:
        print()
        print("=" * 70)
        print("[4] L=16, D=40, dt=0.02, lambda=0 -- PRODUCTION protocol, CLIMBS")
        print("=" * 70)
        L16 = 16
        mps = m.prepare_vacuum(ops, L16, 40)
        a, b = tebd_entropy_run(ops, mps, L16, 40, 0.02, 12.0, d2, L16 // 2,
                                "L=16 D=40", dense=False)
        print("    L=16 : %.6f -> %.6f" % (a, b))
        print("    => the simulated state CLIMBS. The L=16 initial state is")
        print("       itself variational, so this mixes Trotter injection with")
        print("       imperfect preparation; block [3]'s D=100 control is what")
        print("       shows the climb is not intrinsic.")
    else:
        print()
        print("[4] L=16 production climb: skipped (slow). Re-run with --l16.")

    print()
    print("done.")


if __name__ == "__main__":
    main()

"""
T5: 32-site MPS (TEBD) cross-check of the channel-ladder program, in
position space -- a different representation, different algorithm,
different boundary conditions from the momentum-space ED.

MODEL (lattice units a=1): two bosons per site,
  H = H_phi - H_chi + lambda sum_x phi_x^2 chi_x^2,
  H_f = sum_x [pi_x^2/2 + m^2 f_x^2/2] + sum_<xy> (f_x - f_y)^2/2,
open boundaries. The free vacuum is prepared as the ground state of
the sign-flipped auxiliary H_aux = H_phi + H_chi by imaginary-time
TEBD (it is the SAME state as the ghost theory's Fock vacuum: the
chi-sector Bogoliubov vacuum is shared by +-H_chi); the ghost sign
enters the REAL-time evolution only.

PRE-REGISTERED EXPECTATIONS:
  E1. lambda = 0 control: the prepared vacuum is an eigenstate; the
      energy-transfer observable stays at the combined Trotter+
      truncation error floor.
  E2. ghost sign, lambda > 0: energy flows phi-ward at a rate per
      volume that is L-INDEPENDENT (extensivity) for L = 16, 24, 32
      -- the position-space image of the ED box test.
  E3. normal twin (+H_chi): bounded virtual wiggle, no secular
      transfer.
  E4. rate per volume scales ~ lambda^2 (golden-rule regime).
VALIDATION GATES (must pass before E1-E4 are trusted):
  G1. TEBD vs exact diagonalization at L=4, nmax=2: <E_phi>(t)
      curves agree to Trotter accuracy.
  G2. bond-dimension and Trotter-step halving stability.

Observable: Delta E_phi(t) = <H_phi>(t) - <H_phi>(0) summed over the
INTERIOR half of the chain (edge effects excluded), divided by the
interior length: the quantum version of paper 2's classical-lattice
E4 observable.

ASCII only; numpy only (hand-rolled TEBD, ~simple Vidal-less form).
"""

import numpy as np

MASS = 1.0
NMAX = 2                    # bosons per site per field (local dim 9)
CHI = 40                    # bond dimension
LAM = 0.06
W0 = np.sqrt(MASS ** 2 + 2.0)   # reference frequency for site ops


# ----------------------------------------------------------------------
# local operators
# ----------------------------------------------------------------------

def site_ops():
    d = NMAX + 1
    a = np.diag(np.sqrt(np.arange(1, d)), 1)
    x = (a + a.T) / np.sqrt(2 * W0)
    p = 1j * (a.T - a) * np.sqrt(W0 / 2)
    n = a.T @ a
    I = np.eye(d)
    # two-field site: phi (x) chi
    X_phi = np.kron(x, I)
    X_chi = np.kron(I, x)
    P2_phi = np.kron(np.real(p @ p), I)
    P2_chi = np.kron(I, np.real(p @ p))
    N_phi = np.kron(n, I)
    N_chi = np.kron(I, n)
    Id = np.eye(d * d)
    return dict(Xp=X_phi, Xc=X_chi, Pp=P2_phi, Pc=P2_chi,
                Np=N_phi, Nc=N_chi, I=Id, d2=d * d)


def onsite_h(ops, sign_chi, lam):
    """On-site part: pi^2/2 + (m^2+2) f^2/2 per field (the +2 is the
    diagonal part of the NN (f_x - f_y)^2 terms, bulk counting) +
    interaction. Edge sites get (m^2+1) -- handled by caller via a
    correction term."""
    hp = 0.5 * ops['Pp'] + 0.5 * (MASS ** 2 + 2) * ops['Xp'] @ ops['Xp']
    hc = 0.5 * ops['Pc'] + 0.5 * (MASS ** 2 + 2) * ops['Xc'] @ ops['Xc']
    hint = lam * (ops['Xp'] @ ops['Xp']) @ (ops['Xc'] @ ops['Xc'])
    return hp + sign_chi * hc + hint


def bond_h(ops, sign_chi):
    """NN cross terms: -phi_x phi_{x+1} - sign * chi_x chi_{x+1}."""
    return (-np.kron(ops['Xp'], ops['Xp'])
            - sign_chi * np.kron(ops['Xc'], ops['Xc']))


def two_site_gates(ops, sign_chi, lam, L, dt, imag=False):
    """List of two-site gate matrices exp(-i h_bond dt) with on-site
    parts split half-left/half-right (edges get full shares), plus
    edge m^2 correction (+2 -> +1 at chain ends)."""
    d2 = ops['d2']
    hs = onsite_h(ops, sign_chi, lam)
    corr = -0.5 * (ops['Xp'] @ ops['Xp'] +
                   sign_chi * ops['Xc'] @ ops['Xc'])
    gates = []
    for i in range(L - 1):
        h = bond_h(ops, sign_chi)
        # on-site terms split half to each adjacent bond; END sites
        # touch only one bond, so edge bonds carry their FULL share
        wl = 1.0 if i == 0 else 0.5
        wr = 1.0 if i == L - 2 else 0.5
        left = wl * hs + (corr if i == 0 else 0)
        right = wr * hs + (corr if i == L - 2 else 0)
        H2 = h + np.kron(left, np.eye(d2)) + np.kron(np.eye(d2), right)
        z = dt if imag else 1j * dt
        # eigendecomposition exponential
        w, V = np.linalg.eigh((H2 + H2.conj().T) / 2)
        G = (V * np.exp(-z * w)) @ V.conj().T
        gates.append(G.reshape(d2, d2, d2, d2))
    return gates


# ----------------------------------------------------------------------
# minimal MPS machinery
# ----------------------------------------------------------------------

def product_state(L, d2):
    return [np.zeros((1, d2, 1), dtype=complex) for _ in range(L)]


def init_vacuum_guess(L, d2):
    mps = product_state(L, d2)
    for A in mps:
        A[0, 0, 0] = 1.0
    return mps


def apply_gate(mps, G, i, chi_max):
    A, B = mps[i], mps[i + 1]
    cl, d, _ = A.shape
    _, _, cr = B.shape
    theta = np.tensordot(A, B, axes=(2, 0))        # cl d d cr
    theta = np.tensordot(G, theta, axes=([2, 3], [1, 2]))
    theta = theta.transpose(2, 0, 1, 3).reshape(cl * d, d * cr)
    U, S, Vh = np.linalg.svd(theta, full_matrices=False)
    k = min(chi_max, np.sum(S > 1e-12))
    U, S, Vh = U[:, :k], S[:k], Vh[:k]
    mps[i] = U.reshape(cl, d, k)
    mps[i + 1] = (np.diag(S) @ Vh).reshape(k, d, cr)
    return S


def normalize(mps):
    nrm = norm_of(mps)
    mps[-1] = mps[-1] / nrm
    return nrm


def norm_of(mps):
    E = np.ones((1, 1))
    for A in mps:
        E = np.einsum('ab,aic,bid->cd', E, A.conj(), A)
    return np.sqrt(abs(E[0, 0]))


def expectation_sum2(mps, ops2_list, positions):
    """Sum of two-site operator expectations via one left-to-right
    sweep per position (simple O(L^2), fine at these sizes)."""
    total = 0.0
    L = len(mps)
    for (op2, i) in zip(ops2_list, positions):
        E = np.ones((1, 1), dtype=complex)
        j = 0
        while j < L:
            A = mps[j]
            if j == i:
                B = mps[j + 1]
                theta = np.tensordot(A, B, axes=(2, 0))  # cl d d cr
                th2 = np.tensordot(op2, theta,
                                   axes=([2, 3], [1, 2]))
                th2 = th2.transpose(2, 0, 1, 3)
                E = np.einsum('ab,aijc,bijd->cd', E, theta.conj(),
                              th2, optimize=True)
                j += 2
            else:
                E = np.einsum('ab,aic,bid->cd', E, A.conj(), A,
                              optimize=True)
                j += 1
        total += np.real(E[0, 0])
    return total


def tebd_sweep(mps, gates, chi_max, imag=False):
    """First-order sweep (kept for reference)."""
    L = len(mps)
    for i in range(0, L - 1, 2):
        apply_gate(mps, gates[i], i, chi_max)
    for i in range(1, L - 1, 2):
        apply_gate(mps, gates[i], i, chi_max)
    if imag:
        normalize(mps)


def tebd_sweep2(mps, g_half, g_full, chi_max, imag=False):
    """Second-order Trotter: even(dt/2), odd(dt), even(dt/2)."""
    L = len(mps)
    for i in range(0, L - 1, 2):
        apply_gate(mps, g_half[i], i, chi_max)
    for i in range(1, L - 1, 2):
        apply_gate(mps, g_full[i], i, chi_max)
    for i in range(0, L - 1, 2):
        apply_gate(mps, g_half[i], i, chi_max)
    if imag:
        normalize(mps)


def dense_to_mps(psi, L, d2):
    """Exact MPS decomposition of a dense state (small L only)."""
    mps = []
    rest = psi.reshape(1, -1)
    chi_l = 1
    for i in range(L - 1):
        rest = rest.reshape(chi_l * d2, d2 ** (L - 1 - i))
        U, S, Vh = np.linalg.svd(rest, full_matrices=False)
        k = np.sum(S > 1e-12)
        mps.append(U[:, :k].reshape(chi_l, d2, k))
        rest = (np.diag(S[:k]) @ Vh[:k])
        chi_l = k
    mps.append(rest.reshape(chi_l, d2, 1))
    return mps


# ----------------------------------------------------------------------
# observables: H_phi two-site terms over the interior
# ----------------------------------------------------------------------

def hphi_terms(ops, L):
    d2 = ops['d2']
    hp_site = 0.5 * ops['Pp'] + 0.5 * (MASS ** 2 + 2) * \
        ops['Xp'] @ ops['Xp']
    terms, pos = [], []
    lo, hi = L // 4, 3 * L // 4          # interior window
    for i in range(lo, hi):
        op2 = (-np.kron(ops['Xp'], ops['Xp'])
               + 0.5 * np.kron(hp_site, np.eye(d2))
               + 0.5 * np.kron(np.eye(d2), hp_site))
        terms.append(op2.reshape(d2, d2, d2, d2))
        pos.append(i)
    return terms, pos


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def prepare_vacuum(ops, L, chi_max):
    """Annealed second-order imaginary-time TEBD."""
    mps = init_vacuum_guess(L, ops['d2'])
    for (dtau, steps) in ((0.1, 100), (0.03, 150), (0.01, 200)):
        gh = two_site_gates(ops, +1.0, 0.0, L, dtau / 2, imag=True)
        gf = two_site_gates(ops, +1.0, 0.0, L, dtau, imag=True)
        for s in range(steps):
            tebd_sweep2(mps, gh, gf, chi_max, imag=True)
    return mps


def run_case(ops, L, sign_chi, lam, chi_max, dt=0.02, tmax=12.0,
             mps=None):
    if mps is None:
        mps = prepare_vacuum(ops, L, chi_max)
    terms, pos = hphi_terms(ops, L)
    E0 = expectation_sum2(mps, terms, pos)
    gh = two_site_gates(ops, sign_chi, lam, L, dt / 2, imag=False)
    gf = two_site_gates(ops, sign_chi, lam, L, dt, imag=False)
    ts, dE = [0.0], [0.0]
    nst = int(round(tmax / dt))
    every = max(1, int(round(1.0 / dt)))
    for s in range(1, nst + 1):
        tebd_sweep2(mps, gh, gf, chi_max, imag=False)
        if s % every == 0:
            ts.append(s * dt)
            dE.append(expectation_sum2(mps, terms, pos) - E0)
    return np.array(ts), np.array(dE)


def exact_check(ops, L=4, lam=LAM, dt=0.05, tmax=4.0):
    """G1: dense ED of the same OBC lattice model vs TEBD."""
    d2 = ops['d2']
    dim = d2 ** L
    hs = onsite_h(ops, -1.0, lam)
    corr = -0.5 * (ops['Xp'] @ ops['Xp'] - ops['Xc'] @ ops['Xc'])
    hb = bond_h(ops, -1.0)

    def kronN(mats):
        out = np.array([[1.0]])
        for m in mats:
            out = np.kron(out, m)
        return out

    H = np.zeros((dim, dim), dtype=complex)
    I = np.eye(d2)
    for i in range(L):
        h_i = hs + (corr if i in (0, L - 1) else 0)
        H += kronN([I] * i + [h_i] + [I] * (L - 1 - i))
    # bond terms via kron of the two-site matrix
    for i in range(L - 1):
        left = np.eye(d2 ** i)
        right = np.eye(d2 ** (L - 2 - i))
        H += np.kron(np.kron(left, hb), right)
    H = (H + H.conj().T) / 2
    # auxiliary (both-positive) ground state
    hs_aux = onsite_h(ops, +1.0, 0.0)
    corr_aux = -0.5 * (ops['Xp'] @ ops['Xp'] + ops['Xc'] @ ops['Xc'])
    hb_aux = bond_h(ops, +1.0)
    Haux = np.zeros((dim, dim), dtype=complex)
    for i in range(L):
        h_i = hs_aux + (corr_aux if i in (0, L - 1) else 0)
        Haux += np.kron(np.kron(np.eye(d2 ** i), h_i),
                        np.eye(d2 ** (L - 1 - i)))
    for i in range(L - 1):
        Haux += np.kron(np.kron(np.eye(d2 ** i), hb_aux),
                        np.eye(d2 ** (L - 2 - i)))
    Haux = (Haux + Haux.conj().T) / 2
    w, V = np.linalg.eigh(Haux)
    psi = V[:, 0].astype(complex)
    # observable: interior H_phi (sites 1..2 window for L=4)
    hp_site = 0.5 * ops['Pp'] + 0.5 * (MASS ** 2 + 2) * \
        ops['Xp'] @ ops['Xp']
    # SAME interior window as hphi_terms: bonds i in [L//4, 3L//4)
    obs = np.zeros((dim, dim), dtype=complex)
    op2 = (-np.kron(ops['Xp'], ops['Xp'])
           + 0.5 * np.kron(hp_site, np.eye(d2))
           + 0.5 * np.kron(np.eye(d2), hp_site))
    for i in range(L // 4, 3 * L // 4):
        obs += np.kron(np.kron(np.eye(d2 ** i), op2),
                       np.eye(d2 ** (L - 2 - i)))
    from scipy.linalg import expm
    U = expm(-1j * H * dt * 20)
    e0 = np.real(psi.conj() @ obs @ psi)
    vals = []
    ph = psi.copy()
    for s in range(4):
        ph = U @ ph
        vals.append(np.real(ph.conj() @ obs @ ph) - e0)
    return np.array(vals)


def main():
    np.random.seed(0)
    ops = site_ops()
    print("=" * 72)
    print("T5 MPS cross-check (TEBD, position space, OBC): nmax=%d "
          "chi=%d lam=%.2f" % (NMAX, CHI, LAM))
    print("=" * 72)

    # ---------------- G1: ED cross-check at L=4 -------------------------
    print("\n[G1] TEBD vs dense ED, L=4 (interior bond observable):")
    ed = exact_check(ops, L=4)
    mps = prepare_vacuum(ops, 4, CHI)
    terms, pos = hphi_terms(ops, 4)
    E0 = expectation_sum2(mps, terms, pos)
    gates = two_site_gates(ops, -1.0, LAM, 4, 0.05, imag=False)
    tb = []
    for s in range(1, 81):
        tebd_sweep(mps, gates, CHI, imag=False)
        if s % 20 == 0:
            tb.append(expectation_sum2(mps, terms, pos) - E0)
    tb = np.array(tb)
    for k in range(4):
        print("  t=%3.0f : ED %.5e   TEBD %.5e   diff %.1e"
              % ((k + 1) * 1.0, ed[k], tb[k], abs(ed[k] - tb[k])))

    # ---------------- E1: lambda = 0 control ----------------------------
    print("\n[E1] lambda=0 control (L=16): energy-transfer drift "
          "(error floor):")
    ts, dE = run_case(ops, 16, -1.0, 0.0, CHI)
    print("  max |dE_phi| over t<=12: %.2e" % np.max(np.abs(dE)))

    # ---------------- E3: normal twin ------------------------------------
    print("\n[E3] normal twin (L=16, +H_chi): bounded transfer?")
    ts, dEn = run_case(ops, 16, +1.0, LAM, CHI)
    print("  max |dE_phi|: %.3e (bounded virtual)" % np.max(np.abs(dEn)))

    # ---------------- E2: ghost, volume sweep ---------------------------
    print("\n[E2] ghost sign: interior transfer per site vs L:")
    for L in (16, 24, 32):
        ts, dEg = run_case(ops, L, -1.0, LAM, CHI)
        interior = (3 * L // 4) - (L // 4)
        rate = np.polyfit(ts[1:], dEg[1:], 1)[0] / interior
        print("  L=%2d : dE_phi(12)/site = %.4e ; fitted rate/site = "
              "%.4e" % (L, dEg[-1] / interior, rate))

    # ---------------- E4: coupling scaling -------------------------------
    print("\n[E4] coupling scaling (L=16): rate ~ lambda^2 ?")
    base = None
    for lam in (0.03, 0.06, 0.12):
        ts, dEg = run_case(ops, 16, -1.0, lam, CHI)
        r = np.polyfit(ts[1:], dEg[1:], 1)[0]
        if base is None:
            base = (lam, r)
        print("  lambda=%.2f : rate = %.4e  (x%.2f vs lambda-"
              "doubling^2 = 4)" % (lam, r,
                                   r / base[1] if base[1] else 0))

    print("\ndone.")


if __name__ == '__main__':
    main()

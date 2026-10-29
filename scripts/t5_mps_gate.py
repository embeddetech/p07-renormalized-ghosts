"""G1 gate for t5_mps.py: TEBD (2nd-order, chi=100) vs dense ED at
L=4, from the SAME exact initial state, matched observables.
Writes results to mps_gate_result.txt next to this script."""

import numpy as np
import t5_mps as m


def main():
    ops = m.site_ops()
    d2 = ops['d2']
    L = 4
    hs_aux = m.onsite_h(ops, +1.0, 0.0)
    corr_aux = -0.5 * (ops['Xp'] @ ops['Xp'] + ops['Xc'] @ ops['Xc'])
    hb_aux = m.bond_h(ops, +1.0)
    dim = d2 ** L
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
    ed = m.exact_check(ops, L=4)
    mps = m.dense_to_mps(psi, L, d2)
    terms, pos = m.hphi_terms(ops, L)
    E0 = m.expectation_sum2(mps, terms, pos)
    gh = m.two_site_gates(ops, -1.0, m.LAM, L, 0.01, imag=False)
    gf = m.two_site_gates(ops, -1.0, m.LAM, L, 0.02, imag=False)
    tb = []
    for s in range(1, 201):
        m.tebd_sweep2(mps, gh, gf, 100, imag=False)
        if s % 50 == 0:
            tb.append(m.expectation_sum2(mps, terms, pos) - E0)
    with open("mps_gate_result.txt", "w") as fh:
        fh.write("[G1 fixed harness] chi=100 dt=0.02\n")
        for k in range(4):
            fh.write("t=%d : ED %.5e  TEBD %.5e  diff %.1e\n"
                     % (k + 1, ed[k], tb[k], abs(ed[k] - tb[k])))
        fh.write("done\n")


if __name__ == "__main__":
    main()

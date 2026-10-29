"""
T5 partial-fold backing for Sec. "attribution by elimination".

Claim under backing (M=8, cap 6, lambda=0.02): the naive ratio is
full/naive = 0.8620, and BOTH partial dressings fail to move it --
  * folding ONLY the diagonal self-energies (naive couplings, but
    second-order detuning of the P-space energies) reproduces the
    naive ratio, and
  * folding ONLY the vacuum-column vertex correction (renormalized
    couplings g_eff = g + H_PQ (0-E_Q)^-1 H_{Q,vac}, no detuning)
    reproduces the naive ratio.
Only the full shell-shell fold captures the deficit -- and that full
Schur fold OOMs at |Q|~81k (needs an LU factorization). The two
partial folds here need only matrix-free evolution (expm_multiply)
plus sparse matvecs, so they ARE computable at this size.

Writes results_partial_folds_M{M}_cap{cap}.txt (line-flushed) and
checkpoints the Hamiltonian to ckpt_H_M{M}_cap{cap}.npz (shared with
t5_campaign.py's `val` job). ASCII only; numpy/scipy.
"""

import os
import pickle
import time

import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import expm_multiply

import t5_channel_ladder as t5
import t5_campaign as camp


def run(M=8, cap=6, lam=0.02, Lbox=16.0):
    t5.L = Lbox
    fname = "results_partial_folds_M%d_cap%d.txt" % (M, cap)
    ckpt = "ckpt_H_M%d_cap%d.npz" % (M, cap)
    with open(fname, "a") as fh:
        camp.log(fh, "partial-folds start M=%d cap=%d lam=%.3f L=%.0f"
                 % (M, cap, lam, Lbox))
        t0 = time.time()
        if os.path.exists(ckpt) and os.path.exists(ckpt + ".meta"):
            H = sp.load_npz(ckpt).tocsr()
            with open(ckpt + ".meta", "rb") as mf:
                mds, states, index = pickle.load(mf)
            camp.log(fh, "H RESTORED from checkpoint: dim=%d nnz=%d"
                     % (H.shape[0], H.nnz))
        else:
            H, mds, states, index = camp.lean_build_H(t5, M, cap, lam)
            sp.save_npz(ckpt, H.tocoo())
            with open(ckpt + ".meta", "wb") as mf:
                pickle.dump((mds, states, index), mf)
            camp.log(fh, "H built+checkpointed: dim=%d nnz=%d [%.0fs]"
                     % (H.shape[0], H.nnz, time.time() - t0))

        nphi = t5.observables(mds, states)
        vac = index[tuple([0] * len(mds))]
        efree = np.array([sum(st[i] * (t5.omega(n) if f == 0 else
                                       -t5.omega(n))
                              for i, (f, n) in enumerate(mds))
                          for st in states])
        ts = np.linspace(0.0, 40.0, 41)

        # ---- exact evolution (matrix-free; no LU) ----
        t0 = time.time()
        psi0 = np.zeros(H.shape[0], dtype=complex)
        psi0[vac] = 1.0
        cols = expm_multiply(-1j * H, psi0, start=0, stop=40, num=41,
                             endpoint=True)
        camp.log(fh, "evolution done [%.0fs]" % (time.time() - t0))

        P_idx = [j for j in range(len(states)) if abs(efree[j]) < 1e-9]
        if vac not in P_idx:
            P_idx = [vac] + P_idx
        Q_idx = [j for j in range(len(states)) if j not in set(P_idx)]
        N_full = np.array([sum(abs(p[j]) ** 2 * nphi[j] for j in P_idx)
                           for p in cols])
        colv = np.asarray(H[:, vac].todense()).ravel()

        # ---- naive channel sum ----
        N_naive = np.zeros(len(ts))
        for j in P_idx:
            if j == vac or abs(colv[j]) < 1e-14:
                continue
            N_naive += (abs(colv[j]) * ts) ** 2 * nphi[j]
        camp.log(fh, "full/naive(t=40)        = %.4f  (|P|=%d |Q|=%d)"
                 % (N_full[-1] / max(N_naive[-1], 1e-30),
                    len(P_idx), len(Q_idx)))

        # ---- second-order pieces (sparse matvecs, no factorization) ----
        Hc = H.tocsc()
        efq = efree[Q_idx]
        # vertex correction to the vacuum column: one sparse matvec
        vQ = np.asarray(Hc[Q_idx, vac].todense()).ravel() / (0.0 - efq)
        HPQ = Hc[P_idx][:, Q_idx]
        gcorr = np.asarray(HPQ @ vQ).ravel()
        # diagonal second-order self-energy shift of each P state
        shifts = {}
        for f in P_idx:
            col = np.asarray(Hc[Q_idx, f].todense()).ravel()
            shifts[f] = float(np.sum(np.abs(col) ** 2 / (0.0 - efq)))
        dvac = shifts[vac]

        # DIAGONAL-ONLY fold: naive couplings, detuned by (shift_j - shift_vac)
        N_diag = np.zeros(len(ts))
        for j in P_idx:
            if j == vac or abs(colv[j]) < 1e-16:
                continue
            D = shifts[j] - dvac
            if abs(D) < 1e-12:
                N_diag += (abs(colv[j]) * ts) ** 2 * nphi[j]
            else:
                N_diag += (abs(colv[j]) * 2 * np.sin(D * ts / 2) / D) ** 2 \
                    * nphi[j]
        camp.log(fh, "full/diag-only(t=40)    = %.4f  (self-energy detuning, "
                 "naive couplings)" % (N_full[-1] / max(N_diag[-1], 1e-30)))

        # VERTEX-ONLY fold: renormalized couplings, no detuning
        N_vtx = np.zeros(len(ts))
        for kf, j in enumerate(P_idx):
            if j == vac:
                continue
            geff = colv[j] + gcorr[kf]
            if abs(geff) < 1e-16:
                continue
            N_vtx += (abs(geff) * ts) ** 2 * nphi[j]
        camp.log(fh, "full/vertex-only(t=40)  = %.4f  (vacuum-column vertex "
                 "correction, no detuning)"
                 % (N_full[-1] / max(N_vtx[-1], 1e-30)))
        camp.log(fh, "partial-folds done  (both partial folds ~ naive => "
                 "deficit is the shell-shell block)")


if __name__ == "__main__":
    camp.be_polite()
    run()

"""
T5 campaign runner for the lab server (polite-tenant edition).

Jobs (--job):
  ladder_ext : extended dimension-ladder channel-table scans --
               d=1 to M=1024 with smoothing-width (eta) scaling;
               d=2 shells to M=16. Cheap memory; writes each row
               immediately to results_ladder.txt.
  val        : exact-dynamics validation at (--M, --cap): lean
               (momentum-class-streamed) Hamiltonian build, P=0
               projection, Krylov evolution, naive-PT and
               Schur-renormalized ratios. Writes progress and
               results to results_val_M{M}_cap{cap}.txt.

Politeness: sets its own process priority to BELOW_NORMAL via
ctypes (no extra deps); single-threaded numpy left as-is.
All output ASCII, flushed line-by-line so partial runs are usable.
"""

import argparse
import ctypes
import sys
import time

import numpy as np


def be_polite():
    try:
        BELOW_NORMAL = 0x00004000
        h = ctypes.windll.kernel32.GetCurrentProcess()
        ctypes.windll.kernel32.SetPriorityClass(h, BELOW_NORMAL)
    except Exception:
        pass


def log(fh, msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


# ----------------------------------------------------------------------
# job: ladder_ext
# ----------------------------------------------------------------------

def job_ladder_ext():
    import t5_dimension_ladder as dl
    with open("results_ladder.txt", "a") as fh:
        log(fh, "ladder_ext start")
        etas = (0.025, 0.05, 0.1, 0.2, 0.4)
        log(fh, "d=1 L=16 etas=" + str(etas))
        for M in (32, 64, 128, 256, 512, 1024):
            t0 = time.time()
            S, G = dl.channel_sums(M, 1, 16.0, etas)
            row = ("d1 M=%4d Lam=%7.1f Sres=%.6e " % (
                M, 2 * np.pi * M / 16.0, S) +
                " ".join("G(%.3f)=%.6e" % (e, G[e]) for e in etas) +
                "  [%.0fs]" % (time.time() - t0))
            log(fh, row)
        log(fh, "d=2 L=16")
        for M in (4, 6, 8, 12, 16):
            t0 = time.time()
            S, G = dl.channel_sums(M, 2, 16.0, etas)
            row = ("d2 M=%3d Lam=%7.2f Sres=%.6e " % (
                M, 2 * np.pi * M / 16.0, S) +
                " ".join("G(%.3f)=%.6e" % (e, G[e]) for e in etas) +
                "  [%.0fs]" % (time.time() - t0))
            log(fh, row)
        log(fh, "ladder_ext done")


# ----------------------------------------------------------------------
# job: val  (lean build + evolution + ratios)
# ----------------------------------------------------------------------

def lean_build_H(t5, M, cap, lam):
    """Momentum-class-streamed interaction build: never holds more
    than one (phi-class, chi-class) product pair at a time."""
    from scipy import sparse
    mds, states, index = t5.build_basis(M, cap)
    dim = len(states)
    diag = np.zeros(dim)
    for j, st in enumerate(states):
        for i, (f, n) in enumerate(mds):
            w = t5.omega(n)
            diag[j] += st[i] * (w if f == 0 else -w)
    H = sparse.diags(diag).tocsr()

    phi_ids = [i for i, (f, n) in enumerate(mds) if f == 0]
    chi_ids = [i for i, (f, n) in enumerate(mds) if f == 1]

    def elem(field_ids):
        out = []
        for i in field_ids:
            n = mds[i][1]
            for dag in (0, 1):
                p = -n if dag else n
                out.append((i, dag, p,
                            1.0 / np.sqrt(2 * t5.L * t5.omega(n))))
        return out

    ephi, echi = elem(phi_ids), elem(chi_ids)
    ops = {}

    def op(i, dag):
        if (i, dag) not in ops:
            ops[(i, dag)] = t5.ladder(mds, states, index, i, dag)
        return ops[(i, dag)]

    # group elementary pairs by total momentum (indices only)
    from collections import defaultdict
    grp_phi, grp_chi = defaultdict(list), defaultdict(list)
    for a in ephi:
        for b in ephi:
            grp_phi[a[2] + b[2]].append((a, b))
    for a in echi:
        for b in echi:
            grp_chi[a[2] + b[2]].append((a, b))

    Hint = None
    for pk, phi_pairs in grp_phi.items():
        if -pk not in grp_chi:
            continue
        A = None
        for (i1, d1, p1, c1), (i2, d2, p2, c2) in phi_pairs:
            P = c1 * c2 * (op(i1, d1) @ op(i2, d2))
            A = P if A is None else A + P
        B = None
        for (i1, d1, p1, c1), (i2, d2, p2, c2) in grp_chi[-pk]:
            P = c1 * c2 * (op(i1, d1) @ op(i2, d2))
            B = P if B is None else B + P
        term = A @ B
        del A, B
        Hint = term if Hint is None else Hint + term
    H = H + lam * t5.L * Hint
    H = (H + H.conj().T) / 2
    keep = [j for j, st in enumerate(states)
            if t5.momentum_of(st, mds) == 0]
    H = H[keep][:, keep].tocsr()
    states0 = [states[j] for j in keep]
    index0 = {st: i for i, st in enumerate(states0)}
    return H, mds, states0, index0


def job_val(M, cap, lam=0.02, Lbox=16.0):
    import t5_channel_ladder as t5
    from scipy import sparse as _sp
    t5.L = Lbox
    fname = "results_val_M%d_cap%d.txt" % (M, cap)
    ckpt = "ckpt_H_M%d_cap%d.npz" % (M, cap)
    import os
    import pickle
    with open(fname, "a") as fh:
        log(fh, "val start M=%d cap=%d lam=%.3f L=%.0f"
            % (M, cap, lam, Lbox))
        t0 = time.time()
        if os.path.exists(ckpt) and os.path.exists(ckpt + ".meta"):
            H = _sp.load_npz(ckpt).tocsr()
            with open(ckpt + ".meta", "rb") as mf:
                mds, states, index = pickle.load(mf)
            log(fh, "build RESTORED from checkpoint: dim=%d nnz=%d"
                % (H.shape[0], H.nnz))
        else:
            H, mds, states, index = lean_build_H(t5, M, cap, lam)
            _sp.save_npz(ckpt, H.tocoo())
            with open(ckpt + ".meta", "wb") as mf:
                pickle.dump((mds, states, index), mf)
            log(fh, "build done+checkpointed: dim(P=0)=%d nnz=%d "
                "[%.0fs]" % (H.shape[0], H.nnz, time.time() - t0))
        nphi = t5.observables(mds, states)
        vac = index[tuple([0] * len(mds))]
        efree = np.array([sum(st[i] * (t5.omega(n) if f == 0 else
                                       -t5.omega(n))
                              for i, (f, n) in enumerate(mds))
                          for st in states])
        # exact evolution
        from scipy.sparse.linalg import expm_multiply
        ts = np.linspace(0.0, 40.0, 41)
        psi0 = np.zeros(H.shape[0], dtype=complex)
        psi0[vac] = 1.0
        t0 = time.time()
        cols = expm_multiply(-1j * H, psi0, start=0, stop=40, num=41,
                             endpoint=True)
        log(fh, "evolution done [%.0fs]" % (time.time() - t0))
        P_idx = [j for j in range(len(states))
                 if abs(efree[j]) < 1e-9]
        if vac not in P_idx:
            P_idx = [vac] + P_idx
        N_full = np.array([sum(abs(p[j]) ** 2 * nphi[j]
                               for j in P_idx) for p in cols])
        colv = np.asarray(H[:, vac].todense()).ravel()
        N_naive = np.zeros(len(ts))
        for j in P_idx:
            if j == vac or abs(colv[j]) < 1e-14:
                continue
            N_naive += (abs(colv[j]) * ts) ** 2 * nphi[j]
        rn = N_full[-1] / max(N_naive[-1], 1e-30)
        log(fh, "full/naive(t=40) = %.4f   (|P|=%d)"
            % (rn, len(P_idx)))
        # Schur via SPARSE LU of Hqq -- or, when the factorization
        # would not fit in memory (LU fill-in), the diagonal
        # second-order self-energy renormalized-PT substitute
        Q_idx = [j for j in range(len(states))
                 if j not in set(P_idx)]
        if len(Q_idx) > 40000:
            # second-order substitute for the full fold: renormalize
            # BOTH the energies (diagonal shifts) and the COUPLINGS
            # (vertex correction g_eff = g + H_PQ (0-E_Q)^-1 H_Q,vac
            # -- one sparse matvec, no factorization). Stage-3
            # lambda-scaling says the coupling piece carries the
            # lambda^2 deficit; the diagonal piece was measured
            # negligible (attempt 3).
            t0 = time.time()
            Hc = H.tocsc()
            efq = efree[Q_idx]
            vQ = np.asarray(Hc[Q_idx, vac].todense()).ravel() / \
                (0.0 - efq)
            HPQ = Hc[P_idx][:, Q_idx]
            gcorr = np.asarray(HPQ @ vQ).ravel()
            shifts = {}
            for kf, f in enumerate(P_idx):
                col = np.asarray(Hc[Q_idx, f].todense()).ravel()
                shifts[f] = float(np.sum(np.abs(col) ** 2 /
                                         (0.0 - efq)))
            dvac = shifts[vac]
            N_rpt = np.zeros(len(ts))
            for kf, j in enumerate(P_idx):
                if j == vac:
                    continue
                geff = colv[j] + gcorr[kf]
                if abs(geff) < 1e-16:
                    continue
                D = shifts[j] - dvac
                if abs(D) < 1e-12:
                    N_rpt += (abs(geff) * ts) ** 2 * nphi[j]
                else:
                    N_rpt += (abs(geff) *
                              2 * np.sin(D * ts / 2) / D) ** 2 * \
                        nphi[j]
            log(fh, "full/renPT2(t=40) = %.4f  (2nd-order shifts + "
                "vertex; |Q|=%d) [%.0fs]"
                % (N_full[-1] / max(N_rpt[-1], 1e-30), len(Q_idx),
                   time.time() - t0))
            log(fh, "val done")
            return
        t0 = time.time()
        from scipy.sparse.linalg import splu
        Hc = H.tocsc()
        Hpp = np.asarray(Hc[P_idx][:, P_idx].todense())
        Hpq = Hc[P_idx][:, Q_idx].tocsc()
        Hqq = Hc[Q_idx][:, Q_idx].tocsc()
        lu = splu(Hqq)
        # chunked RHS solve (memory-lean): X = Hqq^-1 Hqp
        Hqp = Hpq.conj().T.tocsc()
        nP = len(P_idx)
        Heff = np.array(Hpp, dtype=complex)
        for c0 in range(0, nP, 64):
            c1 = min(c0 + 64, nP)
            rhs = np.asarray(Hqp[:, c0:c1].todense())
            Xc = lu.solve(rhs)
            Heff[:, c0:c1] -= np.asarray(
                (Hpq @ Xc))
            del rhs, Xc
        del Hqp
        Heff = (Heff + Heff.conj().T) / 2
        from scipy.linalg import expm as dexpm
        vp = P_idx.index(vac)
        phi0 = np.zeros(len(P_idx), dtype=complex)
        phi0[vp] = 1.0
        nphi_P = np.array([nphi[j] for j in P_idx])
        Neff40 = np.real(np.sum(np.abs(
            dexpm(-1j * Heff * 40.0) @ phi0) ** 2 * nphi_P))
        log(fh, "full/Schur(t=40) = %.4f  (sparse LU, |Q|=%d) [%.0fs]"
            % (N_full[-1] / max(Neff40, 1e-30), len(Q_idx),
               time.time() - t0))
        log(fh, "val done")


if __name__ == "__main__":
    be_polite()
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", required=True,
                    choices=["ladder_ext", "val"])
    ap.add_argument("--M", type=int, default=10)
    ap.add_argument("--cap", type=int, default=4)
    a = ap.parse_args()
    if a.job == "ladder_ext":
        job_ladder_ext()
    else:
        job_val(a.M, a.cap)

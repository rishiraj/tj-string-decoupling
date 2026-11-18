"""
Latent-frame analysis of the exact 2-hole (delta = 1/8) ground state of the 4x4 t-J model.

Coordinates:  r  = position of hole 2 relative to hole 1 (15 values)
              s~ = spin configuration of the 14 electrons, listed in some r-dependent order ("frame")
The state is a 15 x 3432 matrix psi(r, s~).  Its Schmidt spectrum measures how well
"pair wavefunction (x) spin background" factorises in that frame.
"""
import itertools, pickle, sys, time
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.optimize import minimize
from collections import deque
from tj_core import TJ, square_cluster, POP

L4 = 4; N = 16; NS = 14
def site(x, y): return (x % L4) + L4 * (y % L4)
def xy(s): return s % L4, s // L4
SUB = np.array([(s % L4 + s // L4) & 1 for s in range(N)])
NBR = {s: [site(xy(s)[0]+1, xy(s)[1]), site(xy(s)[0]-1, xy(s)[1]), site(xy(s)[0], xy(s)[1]+1), site(xy(s)[0], xy(s)[1]-1)] for s in range(N)}
R0 = 1                                            # reference pair: holes at 0 and (1,0)
SMASK = np.array([sum(1 << i for i in c) for c in itertools.combinations(range(NS), 7)], dtype=np.int64)
SMASK.sort(); NSP = len(SMASK)                    # 3432
SRANK = -np.ones(1 << NS, dtype=np.int64); SRANK[SMASK] = np.arange(NSP)
BITS = ((SMASK[:, None] >> np.arange(NS)[None, :]) & 1)       # (3432,14) 1 = up

model = TJ(N, square_cluster(L4, L4), 7, 7)
PERMS = dict(Tx=[site(xy(s)[0]+1, xy(s)[1]) for s in range(N)], Ty=[site(xy(s)[0], xy(s)[1]+1) for s in range(N)],
             C4=[site(-xy(s)[1], xy(s)[0]) for s in range(N)], Mx=[site(-xy(s)[0], xy(s)[1]) for s in range(N)])

def sector_ground_state(H, seed=0):
    """Lowest state with K=(0,0), C4=-1, Mx=+1  (d_{x2-y2} pair symmetry)."""
    ops = {k: model.op_site_perm(p) for k, p in PERMS.items()}
    def proj(v):
        for k, lam in (("Tx", 1), ("Ty", 1)):
            w = v.copy(); acc = v.copy()
            for _ in range(3): w = ops[k] @ w; acc += w
            v = acc / 4
        w = v.copy(); acc = v.copy(); ph = 1
        for _ in range(3): w = ops["C4"] @ w; ph *= -1; acc += ph * w
        v = acc / 4
        v = 0.5 * (v + ops["Mx"] @ v)
        return v
    v0 = proj(np.random.default_rng(seed).standard_normal(model.dim)); v0 /= np.linalg.norm(v0)
    w, v = sla.eigsh(H, k=1, which='SA', v0=v0, tol=1e-11)
    psi = proj(v[:, 0]); nrm = np.linalg.norm(psi)
    assert nrm > 0.999, f"symmetry leak {nrm}"
    psi /= nrm
    return psi @ (H @ psi), psi, ops

def build_frame_map():
    """L : (r, s~ in raster order, site-ordered fermion convention) -> K=0 state in full basis."""
    u, d = model.u, model.d
    occ = u | d
    cols_state, cols_f, conv = [], [], []
    for r in range(1, N):
        sel = np.where(((occ & 1) == 0) & (((occ >> r) & 1) == 0))[0]
        osites = [s for s in range(N) if s not in (0, r)]
        us, ds = u[sel], d[sel]
        sm = np.zeros(len(sel), dtype=np.int64)
        for l, s in enumerate(osites): sm |= ((us >> s) & 1) << l
        inv = np.zeros(len(sel), dtype=np.int64)          # pairs (i up, j down, j<i)
        for a, i in enumerate(osites):
            for j in osites[:a]:
                inv += ((us >> i) & 1) & ((ds >> j) & 1)
        cols_state.append(sel); cols_f.append((r - 1) * NSP + SRANK[sm]); conv.append(1 - 2 * (inv & 1))
    st = np.concatenate(cols_state); f = np.concatenate(cols_f); cv = np.concatenate(conv)
    rows, cols, vals = [], [], []
    p = list(range(N))
    for ax in range(4):
        q = p
        for ay in range(4):
            idx, sgn = model.apply_site_perm(q)
            rows.append(idx[st]); cols.append(f); vals.append(0.25 * sgn[st] * cv)
            q = [PERMS["Ty"][s] for s in q]
        p = [PERMS["Tx"][s] for s in p]
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(model.dim, 15 * NSP)).tocsr()

# ---------------------------------------------------------------- frames
def string_path(r):
    """shortest path of hole 2 from R0 to r that avoids the other hole at the origin."""
    prev = {R0: None}; dq = deque([R0])
    while dq:
        s = dq.popleft()
        if s == r: break
        for n in NBR[s]:
            if n != 0 and n not in prev: prev[n] = s; dq.append(n)
    path = [r]
    while prev[path[-1]] is not None: path.append(prev[path[-1]])
    return path[::-1]

def lab0(s): return s - 2                           # reference labels: sites 2..15 -> 0..13

def frame_labels(kind):
    """for each r: array lab[l_raster] = new label, and the lattice site of each raster slot."""
    out = {}
    for r in range(1, N):
        osites = [s for s in range(N) if s not in (0, r)]
        if kind == "raster":
            m = {s: l for l, s in enumerate(osites)}
        elif kind == "site":                        # spins stay on their lattice sites; the one on R0 takes r's label
            m = {s: lab0(s) for s in osites if s != R0}
            if r != R0: m[R0] = lab0(r)
        elif kind == "string":                      # spins are dragged one step along the hole's path
            m = {s: lab0(s) for s in osites if s != R0}
            p = string_path(r)
            for k in range(1, len(p)): m[p[k - 1]] = lab0(p[k])
        out[r] = (np.array([m[s] for s in osites]), np.array(osites))
    return out

def frame_transform(labels, marshall=False):
    """returns newidx[r-1, a], sgn[r-1, a]:  psi_X[r, newidx] = sgn * psi_raster[r, a]."""
    newidx = np.zeros((15, NSP), dtype=np.int64); sgn = np.ones((15, NSP))
    for r in range(1, N):
        lab, osites = labels[r]
        nm = (BITS << lab[None, :]).sum(1)
        newidx[r - 1] = SRANK[nm]
        if marshall:   # (-1)^{# down spins sitting on a sublattice different from that of their label's home site}
            mism = (SUB[osites] != SUB[lab + 2]).astype(np.int64)
            sgn[r - 1] = 1 - 2 * (((1 - BITS) * mism[None, :]).sum(1) & 1)
    return newidx, sgn

def to_frame(M, fr):
    newidx, sgn = fr; out = np.empty_like(M)
    for i in range(15): out[i, newidx[i]] = sgn[i] * M[i]
    return out

def frame_operator(fr):
    newidx, sgn = fr
    rows = (np.arange(15)[:, None] * NSP + newidx).ravel(); cols = np.arange(15 * NSP)
    return sp.coo_matrix((sgn.ravel(), (rows, cols)), shape=(15 * NSP,) * 2).tocsr()

def schmidt(M):
    s = np.linalg.svd(M, compute_uv=False); p = s**2 / (s**2).sum()
    return p, -(p[p > 1e-15] * np.log(p[p > 1e-15])).sum()

# ---------------------------------------------------------------- variational rank-k ansatz in a frame
def optimise_rank_k(Hf, Gf, M0, k, maxiter=4000):
    """minimise <v|H|v>/<v|G|v>, v = sum_a chi_a(r) Phi_a(s~)  (15k + 3432k amplitudes)."""
    U, s, Vt = np.linalg.svd(M0, full_matrices=False)
    x0 = np.concatenate([(U[:, :k] * np.sqrt(s[:k])).ravel(), (Vt[:k].T * np.sqrt(s[:k])).ravel()])
    def f(x):
        A = x[:15 * k].reshape(15, k); B = x[15 * k:].reshape(NSP, k)
        v = (A @ B.T).ravel(); Hv = Hf @ v; Gv = Gf @ v; den = v @ Gv; E = (v @ Hv) / den
        g = (2 * (Hv - E * Gv) / den).reshape(15, NSP)
        return E, np.concatenate([(g @ B).ravel(), (g.T @ A).ravel()])
    res = minimize(f, x0, jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter, maxfun=2 * maxiter, ftol=1e-13, gtol=1e-7))
    A = res.x[:15 * k].reshape(15, k); B = res.x[15 * k:].reshape(NSP, k)
    return res.fun, (A @ B.T)

def fidelity(v, psi_f, Gf):
    """|<Psi_exact | L v>|^2 / <v|G|v>, using psi_f = L^T Psi_exact."""
    v = v.ravel(); return (v @ psi_f.ravel())**2 / (v @ (Gf @ v))

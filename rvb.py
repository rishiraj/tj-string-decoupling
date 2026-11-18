"""Baseline: Gutzwiller-projected d-wave pair (RVB) state, 2 holes on the periodic 4x4 cluster.

On this cluster k=(0,0) is a filled nodal point, so a pure d-wave pair function gives an A1 state that is
orthogonal to the B1 ground state.  Fix: give the four nodal Fermi points (+-pi/2, +-pi/2) an s-like pair
amplitude lam, then project onto the K=0, B1 sector.

With k=0 filled the pair function takes four values: a1 (eps=-2t shell), a2 (antinodes, eps=0), a3 (eps=+2t
shell), all with d-wave signs, and lam (nodal Fermi points).  Only x=a2/a1, z=a3/a1, y=lam/a1 matter, and the
projected wavefunction is a polynomial in them with 22 J-independent component vectors.  We build those once
and then optimise (i) the BCS form a_k = D_k/(xi_k+E_k) over (Delta, mu, lam) and (ii) free (x, z, y).
"""
import numpy as np, sys, os
from scipy.optimize import minimize
from frames import model, N, sector_ground_state
from excite import projector, Tm, HJ, Sp

terms = [(p, s, q) for q in (1, 3) for p in range(3) for s in range(5) if 0 <= 6 - p - s - q <= 4]
c = lambda x, z, y: np.array([x**p * z**s * y**q for p, s, q in terms])
if not os.path.exists("rvb_phi.npy"):
    P = projector(-1)
    pos = {m: [i for i in range(N) if (m >> i) & 1] for m in set(model.u.tolist()) | set(model.d.tolist())}
    RU = np.array([pos[m] for m in model.u.tolist()]); RD = np.array([pos[m] for m in model.d.tolist()])
    DX = ((RU[:, :, None] % 4 - RD[:, None, :] % 4) % 4).astype(np.int8); DY = ((RU[:, :, None] // 4 - RD[:, None, :] // 4) % 4).astype(np.int8)
    del RU, RD
    ks = [(2 * np.pi * a / 4, 2 * np.pi * b / 4) for a in range(4) for b in range(4)]
    def amp(x, z, y):
        a = np.zeros((4, 4))
        for kx, ky in ks:
            ck = int(round(np.cos(kx) + np.cos(ky))); dk = int(round(np.cos(kx) - np.cos(ky)))
            if abs(ck) == 2: continue                               # k=0 handled exactly below; (pi,pi) empty
            ak = y if dk == 0 else np.sign(dk) * {1: 1.0, 0: x, -1: z}[ck]
            for xx in range(4):
                for yy in range(4): a[xx, yy] += ak * np.cos(kx * xx + ky * yy)
        M = a[DX, DY]
        return P(np.linalg.det(M + 1.0) - np.linalg.det(M))         # a_{k=0} -> infinity limit, then projection
    rng = np.random.default_rng(3); pts = rng.uniform(0.4, 1.6, size=(len(terms) + 5, 3))
    F = np.array([amp(*p) for p in pts])
    A = np.array([c(*p) for p in pts[:-1]])
    Phi = np.linalg.lstsq(A, F[:-1], rcond=None)[0]
    res = np.linalg.norm(c(*pts[-1]) @ Phi - F[-1]) / np.linalg.norm(F[-1])
    print(f"{len(terms)} components; held-out residual {res:.1e}; cond = {np.linalg.cond(A):.0f}", flush=True)
    np.save("rvb_phi.npy", Phi)
Phi = np.load("rvb_phi.npy")
S = Phi @ Phi.T; T5 = Phi @ (Tm @ Phi.T); J5 = Phi @ (HJ @ Phi.T)

def bcs_ratios(D, mu):
    a = lambda eps, dk: dk * D / ((eps - mu) + np.sqrt((eps - mu)**2 + (dk * D)**2))
    a1 = a(-2.0, 1.0); return a(0.0, 2.0) / a1, a(2.0, 1.0) / a1

for J in [float(s) for s in sys.argv[1:]] or [0.4]:
    H = Tm + J * HJ; E0, psi, _ = sector_ground_state(H); Hm = T5 + J * J5
    def e_free(p):
        v = c(*p); return (v @ Hm @ v) / (v @ S @ v)
    def e_bcs(p):
        if p[0] <= 1e-3: return 1e3
        x, z = bcs_ratios(p[0], p[1]); return e_free((x, z, p[2]))
    out = []
    for f, grid in ((e_bcs, [(D, mu, y) for D in (0.1, 0.3, 1, 3) for mu in (-1.9, -1, -0.3, 0.3, 1) for y in (-3, -1, -0.3, 0.3, 1, 3)]),
                    (e_free, [(x, z, y) for x in (0.1, 0.4, 1, 2) for z in (0.03, 0.2, 0.6, 1.5) for y in (-3, -1, -0.3, 0.3, 1, 3)])):
        cand = sorted((f(g), g) for g in grid)[:4]
        r = min((minimize(f, g, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-10, maxiter=3000)) for _, g in cand), key=lambda r: r.fun)
        xz = bcs_ratios(r.x[0], r.x[1]) + (r.x[2],) if f is e_bcs else tuple(r.x)
        w = c(*xz) @ Phi; w /= np.linalg.norm(w)
        out.append((r.x, r.fun, (w @ psi)**2, np.linalg.norm(Sp @ w)**2, w))
    (pb, eb, fb, sb, wb), (pf, ef, ff, sf, wf) = out
    print(f"J={J}: exact {E0:.5f}")
    print(f"   d-BCS form  (Delta={pb[0]:.2f}, mu={pb[1]:.2f}, lam/a1={pb[2]:.2f}): E={eb:.5f}  dE/N={(eb-E0)/16:.4f}  fidelity={fb:.4f}  S(S+1)={sb:.2f}")
    print(f"   free ratios (a2/a1={pf[0]:.2f}, a3/a1={pf[1]:.2f}, lam/a1={pf[2]:.2f}): E={ef:.5f}  dE/N={(ef-E0)/16:.4f}  fidelity={ff:.4f}  S(S+1)={sf:.2f}", flush=True)

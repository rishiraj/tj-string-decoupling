"""Sanity check of the RVB machinery at half filling: projected d-wave pair state vs exact 4x4 Heisenberg ground state."""
import numpy as np, scipy.sparse.linalg as sla
from scipy.optimize import minimize
from frames import N, PERMS
from strings import half
ops = {k: half.op_site_perm(p) for k, p in PERMS.items()}
def P(v, c4=1):
    for k in ("Tx", "Ty"):
        w = v.copy(); acc = v.copy()
        for _ in range(3): w = ops[k] @ w; acc += w
        v = acc / 4
    w = v.copy(); acc = v.copy(); ph = 1
    for _ in range(3): w = ops["C4"] @ w; ph *= c4; acc += ph * w
    v = acc / 4
    return 0.5 * (v + ops["Mx"] @ v)
HJ = half.hamiltonian(0.0, 1.0); E0, phi = sla.eigsh(HJ, k=1, which='SA', tol=1e-12); phi = phi[:, 0]; E0 = E0[0]
pos = {m: [i for i in range(N) if (m >> i) & 1] for m in set(half.u.tolist())}
RU = np.array([pos[m] for m in half.u.tolist()]); RD = np.array([pos[m] for m in half.d.tolist()])
DX = (RU[:, :, None] % 4 - RD[:, None, :] % 4) % 4; DY = (RU[:, :, None] // 4 - RD[:, None, :] // 4) % 4
ks = [(2 * np.pi * a / 4, 2 * np.pi * b / 4) for a in range(4) for b in range(4)]
def state(p, c4=1):
    x, z, y = p; a = np.zeros((4, 4))
    for kx, ky in ks:
        ck = int(round(np.cos(kx) + np.cos(ky))); dk = int(round(np.cos(kx) - np.cos(ky)))
        if abs(ck) == 2: continue
        ak = y if dk == 0 else np.sign(dk) * {1: 1.0, 0: x, -1: z}[ck]
        for xx in range(4):
            for yy in range(4): a[xx, yy] += ak * np.cos(kx * xx + ky * yy)
    M = a[DX, DY]; w = P(np.linalg.det(M + 1.0) - np.linalg.det(M), c4); n = np.linalg.norm(w)
    return w / n if n > 1e-12 else None
def e(p):
    w = state(p); return 1e3 if w is None else w @ (HJ @ w)
cand = sorted((e((x, z, y)), (x, z, y)) for x in (0.2, 0.6, 1, 2) for z in (0.05, 0.3, 1) for y in (-2, -0.5, 0.5, 2))[:3]
r = min((minimize(e, g, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-9, maxiter=800)) for _, g in cand), key=lambda r: r.fun)
w = state(r.x)
print(f"half filling: exact E/N (S.S part) = {(E0+8)/16:.5f};  projected pair state (a2/a1={r.x[0]:.2f}, a3/a1={r.x[1]:.2f}, lam/a1={r.x[2]:.2f}): "
      f"E/N = {(r.fun+8)/16:.5f}, fidelity = {(w@phi)**2:.4f}")

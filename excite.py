"""Low-energy excitations: vibrational (same pair symmetry) and rotational (s-wave pair) string states vs exact.
Exact states are resolved by momentum K=0, point group (B1 = d-wave, A1 = s-wave) and total spin S=0.
For Sz=0 states S^2 = S^- S^+, applied matrix-free (the explicit S^2 matrix ran out of memory)."""
import sys
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.special import ai_zeros
from tj_core import TJ, square_cluster, POP
from frames import model, site, xy, N, PERMS
from strings import half, pair_operator, krylov

Tm = model.hamiltonian(1.0, 0.0); HJ = model.hamiltonian(0.0, 1.0)
ops = {k: model.op_site_perm(p) for k, p in PERMS.items()}
up1 = TJ(N, square_cluster(4, 4), 8, 6)
def s_plus():
    u, d = model.u, model.d; rows, cols, vals = [], [], []
    for i in range(N):
        idx = np.where((d >> i) & 1 == 1)[0]
        s1 = 1 - 2 * ((POP[u[idx]] + POP[d[idx] & ((1 << i) - 1)]) & 1)     # c_{i,dn}
        s2 = 1 - 2 * (POP[u[idx] & ((1 << i) - 1)] & 1)                     # c^+_{i,up}
        rows.append(up1.index(u[idx] | (1 << i), d[idx] ^ (1 << i))); cols.append(idx); vals.append(s1 * s2)
    return sp.coo_matrix((np.concatenate(vals).astype(float), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(up1.dim, model.dim)).tocsr()
Sp = s_plus(); SpT = Sp.T.tocsr()

def projector(c4):
    def P(v):
        for k in ("Tx", "Ty"):
            w = v.copy(); acc = v.copy()
            for _ in range(3): w = ops[k] @ w; acc += w
            v = acc / 4
        w = v.copy(); acc = v.copy(); ph = 1
        for _ in range(3): w = ops["C4"] @ w; ph *= c4; acc += ph * w
        v = acc / 4
        return 0.5 * (v + ops["Mx"] @ v)
    return P

def exact_sector(H, c4, k=3):
    P = projector(c4)
    def mv(v):
        v = v.ravel(); pv = P(v)
        return H @ pv + 3.0 * (SpT @ (Sp @ pv)) + 25.0 * (v - pv)
    A = sla.LinearOperator((model.dim,) * 2, matvec=mv, dtype=float)
    v0 = P(np.random.default_rng(1).standard_normal(model.dim))
    w, v = sla.eigsh(A, k=k, which='SA', v0=v0, tol=1e-8, ncv=30)
    o = np.argsort(w); w = w[o]; v = v[:, o]
    s2 = [np.linalg.norm(Sp @ v[:, i])**2 for i in range(k)]
    return w - 3.0 * np.array(s2), np.array(s2), v

if __name__ == "__main__":
    _, phi = sla.eigsh(half.hamiltonian(0.0, 1.0), k=1, which='SA', tol=1e-12); phi = phi[:, 0]
    seeds = {"d": pair_operator([(1, 0), (0, 1)], [1, -1]) @ phi, "s": pair_operator([(1, 0), (0, 1)], [1, 1]) @ phi}
    basis = {}
    for k, s in seeds.items():
        print(f"seed {k}: norm^2 = {s@s:.4f}")
        V = krylov(Tm, s / np.linalg.norm(s), 12); basis[k] = (V, V.T @ (Tm @ V), V.T @ (HJ @ V))
    V, Tk, Vk = basis["d"]
    b = np.diag(Tk, 1)[:6].mean(); V0 = Vk[0, 0]; dV = np.diag(Vk) - V0
    sig = np.polyfit(np.arange(5), dV[:5], 1)[0]
    print(f"effective string model (d channel): b = {b:.3f} t, V0 = {V0:.3f} J, tension sigma = {sig:.3f} J (fit l<=4)")
    Vs_, Ts, Vs = basis["s"]
    print("s channel: b_l =", np.round(np.diag(Ts, 1)[:6], 3), " V_l - V_0(d) =", np.round(np.diag(Vs)[:6] - V0, 3))
    a1 = -ai_zeros(1)[0][0]
    for J in [float(x) for x in sys.argv[1:]] or [0.2, 0.4, 0.8]:
        H = Tm + J * HJ
        ed, sd, vd = exact_sector(H, -1); es, ss, vs = exact_sector(H, +1)
        wd, cd = np.linalg.eigh(Tk + J * Vk); ws, cs = np.linalg.eigh(Ts + J * Vs)
        airy = V0 * J - 2 * b + a1 * b**(1/3) * (sig * J)**(2/3)
        fd = [(vd[:, i] @ (V @ cd[:, i]))**2 for i in range(2)]; fs = (vs[:, 0] @ (Vs_ @ cs[:, 0]))**2
        print(f"J={J}: exact B1(d): E={np.round(ed,4)} S(S+1)={np.round(sd,3)} | exact A1(s): E={np.round(es,4)} S(S+1)={np.round(ss,3)}")
        print(f"       string  d: {np.round(wd[:3],4)}   s: {np.round(ws[:3],4)}")
        print(f"       d->s gap: exact {es[0]-ed[0]:.4f}  string {ws[0]-wd[0]:.4f} (s-state fidelity {fs:.3f}) | "
              f"d vibrational gap: exact {ed[1]-ed[0]:.4f}  string {wd[1]-wd[0]:.4f} (fidelities {fd[0]:.3f}, {fd[1]:.3f})")
        print(f"       Airy closed form: {airy:.4f}   string {wd[0]:.4f}   exact {ed[0]:.4f}", flush=True)

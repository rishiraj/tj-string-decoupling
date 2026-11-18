"""Physics checks for the string ansatz on the 4x4 cluster, 2 holes:
sector dimension, staggered magnetisation, and scaling of the pair energy with J/t."""
import numpy as np, scipy.sparse.linalg as sla
from scipy.optimize import curve_fit
from tj_core import TJ, square_cluster, POP
from frames import model, site, xy, N, PERMS, sector_ground_state
from strings import half, pair_operator, krylov

Tm = model.hamiltonian(1.0, 0.0); HJ = model.hamiltonian(0.0, 1.0)

def sector_dim(mod):
    """dimension of the K=0, B1 sector from characters."""
    tr = 0.0; c4 = PERMS["C4"]; mx = PERMS["Mx"]; pg = []
    for k in range(4):
        for mm in range(2):
            p = list(range(N))
            for _ in range(mm): p = [mx[s] for s in p]
            for _ in range(k): p = [c4[s] for s in p]
            pg.append((p, (-1) ** k))
    for a in range(4):
        for b in range(4):
            for p, chi in pg:
                q = [site(xy(s)[0] + a, xy(s)[1] + b) for s in p]
                idx, sgn = mod.apply_site_perm(q)
                tr += chi * sgn[idx == np.arange(mod.dim)].sum()
    return tr / (16 * 8)
d0 = sector_dim(model); d1 = sector_dim(TJ(N, square_cluster(4, 4), 8, 6))
print(f"sector (K=0, B1) dimension: Sz=0 -> {d0:.0f}, Sz=1 -> {d1:.0f}  => singlets in sector: {d0-d1:.0f}")

Eaf, phi = sla.eigsh(half.hamiltonian(0.0, 1.0), k=1, which='SA', tol=1e-12); phi = phi[:, 0]; Eaf = Eaf[0]
seed = pair_operator([(1, 0), (0, 1)], [1, -1]) @ phi; seed /= np.linalg.norm(seed)
V = krylov(Tm, seed, 8); Tk = V.T @ (Tm @ V); Vk = V.T @ (HJ @ V)

A = sum(1 << s for s in range(N) if (xy(s)[0] + xy(s)[1]) % 2 == 0); B = ((1 << N) - 1) ^ A
def m2(mod, psi):
    ms = 0.5 * (POP[mod.u & A] - POP[mod.u & B] - POP[mod.d & A] + POP[mod.d & B])
    return 3 * (psi**2 * ms**2).sum() / N**2
print(f"undoped AF vacuum: <m_s^2> = {m2(half, phi):.4f}")

Jf = np.array([0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]); Ex, Es = [], []
print("  J/t    E_exact    E_string(n=8)  dE/N     fidelity  m2_exact  m2_string")
for J in Jf:
    H = Tm + J * HJ; E0, psi, _ = sector_ground_state(H)
    ev, c = np.linalg.eigh(Tk + J * Vk); s8 = V @ c[:, 0]
    print(f"  {J:4.2f}  {E0:10.5f}  {ev[0]:10.5f}   {(ev[0]-E0)/16:7.5f}  {(s8@psi)**2:7.4f}   {m2(model, psi):.4f}    {m2(model, s8):.4f}", flush=True)
    Ex.append(E0 - Eaf * J); Es.append(ev[0] - Eaf * J)
f = lambda J, a, b, nu: a + b * J**nu
for name, y in (("exact", np.array(Ex)), ("string n=8", np.array(Es))):
    p, _ = curve_fit(f, Jf, y, p0=(-6, 5, 0.8))
    print(f"E_2h - E_0h = a + b J^nu  [{name}]: a={p[0]:.3f} t, b={p[1]:.3f}, nu={p[2]:.3f}")

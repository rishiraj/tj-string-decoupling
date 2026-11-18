"""
String-length ("frozen background") representation of the two-hole ground state.

  |Psi_n> = sum_{l=0}^{n} a_l  T^l  Delta_d |Phi_AF>

Phi_AF  : exact undoped Heisenberg ground state (the "AF vacuum", contains no free parameter)
Delta_d : nearest-neighbour d_{x2-y2} singlet pair annihilation operator
T       : projected hopping.  T^l generates all strings of length l (retraced and looped paths included).
Compare with the Lanczos space of the full H from the same seed (which also lets the background relax).
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, pickle, sys
from tj_core import TJ, square_cluster, POP
from frames import model, site, xy, N, PERMS

half = TJ(N, square_cluster(4, 4), 8, 8)

def pair_annihilate(i, j):
    """Delta_ij = c_{i,up} c_{j,dn} + c_{j,up} c_{i,dn}  : (8,8) sector -> (7,7) sector, sparse matrix."""
    u, d = half.u, half.d
    rows, cols, vals = [], [], []
    for (a, b) in ((i, j), (j, i)):                       # c_{a,up} c_{b,dn}
        m = ((u >> a) & 1 == 1) & ((d >> b) & 1 == 1)
        idx = np.where(m)[0]
        s_dn = 1 - 2 * ((POP[u[idx]] + POP[d[idx] & ((1 << b) - 1)]) & 1)
        s_up = 1 - 2 * (POP[u[idx] & ((1 << a) - 1)] & 1)
        rows.append(model.index(u[idx] ^ (1 << a), d[idx] ^ (1 << b))); cols.append(idx); vals.append(s_up * s_dn)
    return sp.coo_matrix((np.concatenate(vals).astype(float), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(model.dim, half.dim)).tocsr()

def pair_operator(shell, sym):
    """sum over sites of pair operators at separations in `shell`, with sign pattern `sym`."""
    D = None
    for (dx, dy), sg in zip(shell, sym):
        for s in range(N):
            x, y = xy(s); op = sg * pair_annihilate(s, site(x + dx, y + dy))
            D = op if D is None else D + op
    return D

def krylov(A, v0, n):
    """orthonormal basis of span{v0, A v0, ..., A^n v0} with full re-orthogonalisation."""
    V = [v0 / np.linalg.norm(v0)]
    for _ in range(n):
        w = A @ V[-1]
        for _ in range(2):
            for q in V: w -= (q @ w) * q
        nw = np.linalg.norm(w)
        if nw < 1e-10: break
        V.append(w / nw)
    return np.array(V).T

if __name__ == "__main__":
    nmax = 14
    Tm = model.hamiltonian(1.0, 0.0); HJ = model.hamiltonian(0.0, 1.0)
    Eaf, phi = sla.eigsh(half.hamiltonian(0.0, 1.0), k=1, which='SA', tol=1e-12); phi = phi[:, 0]
    print(f"AF vacuum: E/J = {Eaf[0]:.6f}")
    Dd = pair_operator([(1, 0), (0, 1)], [1, -1])
    seed = Dd @ phi; seed /= np.linalg.norm(seed)
    print(f"seed energy/J <seed|H_J|seed> = {seed @ (HJ @ seed):.5f},  <seed|T|seed> = {seed @ (Tm @ seed):.5f}")
    V = krylov(Tm, seed, nmax)
    Tk = V.T @ (Tm @ V); Vk = V.T @ (HJ @ V)
    print("string hopping amplitudes b_l = <l|T|l+1>:", np.round(np.diag(Tk, 1), 4))
    print("diag <l|T|l>:", np.round(np.diag(Tk), 4))
    print("string potential V_l = <l|H_J|l> - V_0 (units of J):", np.round(np.diag(Vk) - Vk[0, 0], 4))
    print("off-diagonal <l|H_J|l+1>:", np.round(np.diag(Vk, 1), 4), " <l|H_J|l+2>:", np.round(np.diag(Vk, 2), 4))
    gs = pickle.load(open("gs44_sector.pkl", "rb")) if len(sys.argv) > 1 else {}
    out = dict(Tk=Tk, Vk=Vk, rows=[])
    for J in (0.2, 0.3, 0.4, 0.5, 0.8):
        H = Tm + J * HJ
        from frames import sector_ground_state
        E0, psi, _ = sector_ground_state(H)
        Vl = krylov(H, seed, nmax); Hl = Vl.T @ (H @ Vl)
        Hk = Tk + J * Vk
        print(f"\nJ/t = {J}  (U/t ~ {4/J:.0f}):  E_exact = {E0:.5f}   bare-pair weight Z = |<exact|Delta_d|AF>|^2 = {(psi@seed)**2:.4f}")
        print("   n   E_string(n)   dE/N    fidelity |  E_Lanczos-H(n)  dE/N   fidelity")
        for n in (0, 1, 2, 3, 4, 6, 8, 10, 12, 14):
            w, c = np.linalg.eigh(Hk[:n+1, :n+1]); f1 = (psi @ (V[:, :n+1] @ c[:, 0]))**2
            w2, c2 = np.linalg.eigh(Hl[:n+1, :n+1]); f2 = (psi @ (Vl[:, :n+1] @ c2[:, 0]))**2
            print(f"  {n:2d}  {w[0]:10.5f}  {(w[0]-E0)/16:7.4f}  {f1:7.4f}  |  {w2[0]:10.5f}  {(w2[0]-E0)/16:7.4f}  {f2:7.4f}")
            out["rows"].append((J, n, E0, w[0], f1, w2[0], f2))
    pickle.dump(out, open("strings.pkl", "wb"))

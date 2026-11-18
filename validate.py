"""Validation of the t-J code + test of the Schrieffer-Wolff (U/t >> 1) step against the Hubbard model."""
import itertools, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from tj_core import TJ, square_cluster, POP, between_mask

def hubbard_T_and_D(N, bonds, nup, ndn):
    um = np.array([sum(1<<i for i in c) for c in itertools.combinations(range(N), nup)], dtype=np.int64)
    dm = np.array([sum(1<<i for i in c) for c in itertools.combinations(range(N), ndn)], dtype=np.int64)
    def hop(masks):
        lut = {m:k for k,m in enumerate(masks)}
        r,c,v = [],[],[]
        for (i,j) in bonds:
            bm = between_mask(i,j)
            for s,dd in ((i,j),(j,i)):
                for k,m in enumerate(masks):
                    if (m>>s)&1 and not (m>>dd)&1:
                        r.append(lut[(m^(1<<s))|(1<<dd)]); c.append(k); v.append(-(1-2*(POP[m&bm]&1)))
        return sp.coo_matrix((v,(r,c)),shape=(len(masks),)*2).tocsr()
    Tu, Td = hop(um), hop(dm)
    T = sp.kron(Tu, sp.identity(len(dm))) + sp.kron(sp.identity(len(um)), Td)
    U_, D_ = np.repeat(um, len(dm)), np.tile(dm, len(um))
    return T.tocsr(), POP[U_ & D_], U_, D_

# --- (1) half-filled 4x4 Heisenberg benchmark:  E = -0.701780 J per site (S.S) - J/4 per bond
m = TJ(16, square_cluster(4,4), 8, 8)
E = sla.eigsh(m.hamiltonian(0.0, 1.0), k=1, which='SA')[0][0]
print(f"(1) 4x4 Heisenberg: E_SS/N = {(E+8)/16:.6f}   (literature -0.701780)")

# --- (2) matrix identity checks on the sqrt(10) cluster against Hubbard-derived operators
N=10; bonds = sorted({(min(i,(i+s)%N), max(i,(i+s)%N)) for i in range(N) for s in (1,3)})
for (nu,nd) in ((5,5),(4,4)):
    T, docc, U_, D_ = hubbard_T_and_D(N, bonds, nu, nd)
    P0 = np.where(docc==0)[0]; P1 = np.where(docc==1)[0]
    tj = TJ(N, bonds, nu, nd)
    # map Hubbard no-doublon states -> tJ index
    loc = tj.index(U_[P0], D_[P0])
    Pm = sp.coo_matrix((np.ones(len(P0)), (loc, P0)), shape=(tj.dim, T.shape[0])).tocsr()
    Q1 = sp.coo_matrix((np.ones(len(P1)), (P1, P1)), shape=T.shape).tocsr()
    PTP = Pm @ T @ Pm.T
    H2  = -(Pm @ T @ Q1 @ T @ Pm.T)        # times 1/U ; contains J-term + 3-site term, J = 4/U
    Hkin = tj.hamiltonian(1.0, 0.0); HJ = tj.hamiltonian(0.0, 4.0)
    print(f"(2) N=10 ({nu},{nd}): |PTP - H_t| = {abs(PTP-Hkin).max():.1e};  "
          f"|SW2 - H_J| = {abs(H2-HJ).max():.2e}  ({'must vanish at half filling' if nu+nd==N else '3-site terms, nonzero by design'})")
    if nu+nd < N:
        print("    U/t     E_Hubbard     E_tJ(J=4t^2/U)   E_tJ+3site     err(tJ)/N    err(tJ+3s)/N")
        for U in (8, 12, 16, 24, 32, 48, 64):
            Hh = T + U*sp.diags(docc.astype(float))
            Eh = sla.eigsh(Hh, k=1, which='SA')[0][0]
            Et = np.linalg.eigvalsh((Hkin + HJ/U).toarray())[0]
            E3 = np.linalg.eigvalsh((PTP + H2/U).toarray())[0]
            print(f"    {U:4d}   {Eh:11.6f}   {Et:11.6f}      {E3:11.6f}    {(Et-Eh)/N:+.5f}     {(E3-Eh)/N:+.6f}")

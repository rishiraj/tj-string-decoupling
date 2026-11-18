"""Check the pair-determinant machinery against an independent construction (projected Slater determinant)."""
import numpy as np, scipy.sparse.linalg as sla
from frames import N
from strings import half
HJ = half.hamiltonian(0.0, 1.0)
pos = {m: [i for i in range(N) if (m >> i) & 1] for m in set(half.u.tolist())}
RU = np.array([pos[m] for m in half.u.tolist()]); RD = np.array([pos[m] for m in half.d.tolist()])
q = np.pi / 2
S = [(0, 0), (q, 0), (-q, 0), (0, q), (0, -q), (np.pi, 0), (q, q), (-q, -q)]
X = lambda R: R % 4
Y = lambda R: R // 4
# independent: det_up * det_dn with plane-wave orbitals
phi_up = np.stack([np.exp(1j * (kx * X(RU) + ky * Y(RU))) for kx, ky in S], axis=-1)      # (dim, 8 electrons, 8 orbitals)
phi_dn = np.stack([np.exp(-1j * (kx * X(RD) + ky * Y(RD))) for kx, ky in S], axis=-1)
slater = np.linalg.det(phi_up) * np.linalg.det(phi_dn)
# pair machinery: a_k = 1 on S
a = np.zeros((4, 4), dtype=complex)
for kx, ky in S:
    for x in range(4):
        for y in range(4): a[x, y] += np.exp(1j * (kx * x + ky * y))
DX = (X(RU)[:, :, None] - X(RD)[:, None, :]) % 4; DY = (Y(RU)[:, :, None] - Y(RD)[:, None, :]) % 4
pair = np.linalg.det(a[DX, DY])
r = pair / slater; ok = np.abs(slater) > 1e-9
print("pair/slater ratio: min", np.round(r[ok].min(), 6), "max", np.round(r[ok].max(), 6), " |pair| where slater=0:", np.abs(pair[~ok]).max() if (~ok).any() else 0)
v = slater / np.linalg.norm(slater)
print("projected Slater determinant energy E/N (S.S part):", ((v.conj() @ (HJ @ v)).real + 8) / 16)

"""
Exact diagonalisation of the t-J model (strong-coupling limit of the repulsive
Hubbard model) on small periodic clusters.

Convention: |u,d> = prod_{i in u, ascending} c^+_{i,up} prod_{j in d, ascending} c^+_{j,dn} |0>
u, d are integer bitmasks; no double occupancy means u & d == 0.

H = -t sum_{<ij>,s} P (c^+_{is} c_{js} + h.c.) P + J sum_{<ij>} (S_i.S_j - n_i n_j / 4)
"""
import itertools
import numpy as np
import scipy.sparse as sp

POP = np.array([bin(i).count("1") for i in range(1 << 16)], dtype=np.int64)


def square_cluster(Lx, Ly):
    """Sites s = x + Lx*y, periodic. Returns list of nearest-neighbour bonds (i<j, unique)."""
    bonds = set()
    for y in range(Ly):
        for x in range(Lx):
            s = x + Lx * y
            for dx, dy in ((1, 0), (0, 1)):
                s2 = (x + dx) % Lx + Lx * ((y + dy) % Ly)
                if s2 != s:
                    bonds.add((min(s, s2), max(s, s2)))
    return sorted(bonds)


def between_mask(a, b):
    lo, hi = min(a, b), max(a, b)
    return ((1 << hi) - 1) & ~((1 << (lo + 1)) - 1)


class TJ:
    def __init__(self, N, bonds, nup, ndn):
        self.N, self.bonds, self.nup, self.ndn = N, bonds, nup, ndn
        us, ds = [], []
        for up_sites in itertools.combinations(range(N), nup):
            u = sum(1 << i for i in up_sites)
            rest = [i for i in range(N) if not (u >> i) & 1]
            for dn_sites in itertools.combinations(rest, ndn):
                us.append(u)
                ds.append(sum(1 << i for i in dn_sites))
        u = np.array(us, dtype=np.int64)
        d = np.array(ds, dtype=np.int64)
        key = (u << N) | d
        order = np.argsort(key)
        self.u, self.d, self.key = u[order], d[order], key[order]
        self.dim = len(self.key)

    def index(self, u, d):
        k = (u << self.N) | d
        idx = np.searchsorted(self.key, k)
        idx = np.clip(idx, 0, self.dim - 1)
        assert np.all(self.key[idx] == k)
        return idx

    def hamiltonian(self, t, J):
        N, u, d = self.N, self.u, self.d
        rows, cols, vals = [], [], []
        diag = np.zeros(self.dim)
        allidx = np.arange(self.dim)
        for (i, j) in self.bonds:
            bi, bj = 1 << i, 1 << j
            bm = between_mask(i, j)
            ui, uj = (u & bi) != 0, (u & bj) != 0
            di, dj = (d & bi) != 0, (d & bj) != 0
            ni, nj = ui | di, uj | dj
            # --- hopping (both directions), up then down
            for (src, dst) in ((i, j), (j, i)):
                bs, bd = 1 << src, 1 << dst
                # up electron hops src -> dst, dst must be empty of everything
                m = ((u & bs) != 0) & ((u & bd) == 0) & ((d & bd) == 0)
                if m.any():
                    un = (u[m] ^ bs) | bd
                    sgn = 1 - 2 * (POP[u[m] & bm] & 1)
                    rows.append(self.index(un, d[m])); cols.append(allidx[m]); vals.append(-t * sgn)
                m = ((d & bs) != 0) & ((d & bd) == 0) & ((u & bd) == 0)
                if m.any():
                    dn = (d[m] ^ bs) | bd
                    sgn = 1 - 2 * (POP[d[m] & bm] & 1)
                    rows.append(self.index(u[m], dn)); cols.append(allidx[m]); vals.append(-t * sgn)
            # --- J (Sz Sz - n n /4): parallel -> 0, antiparallel -> -1/2
            anti = (ui & dj) | (di & uj)
            diag += np.where(anti, -0.5 * J, 0.0)
            # --- J/2 (S+_i S-_j + h.c.);  S+_i S-_j = -(c+_{iu} c_{ju}) (c+_{jd} c_{id})
            for (a, b) in ((i, j), (j, i)):  # S+_a S-_b : needs dn at a, up at b
                ba, bb = 1 << a, 1 << b
                m = ((d & ba) != 0) & ((u & bb) != 0)
                if m.any():
                    su = 1 - 2 * (POP[u[m] & bm] & 1)   # up hop b->a
                    sd = 1 - 2 * (POP[d[m] & bm] & 1)   # dn hop a->b
                    un = (u[m] ^ bb) | ba
                    dn = (d[m] ^ ba) | bb
                    rows.append(self.index(un, dn)); cols.append(allidx[m]); vals.append(-0.5 * J * su * sd)
        rows.append(allidx); cols.append(allidx); vals.append(diag)
        H = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(self.dim, self.dim)).tocsr()
        return H

    # ---------- lattice symmetries (site permutation p: site s -> p[s]) ----------
    def perm_sign_table(self, p, npart):
        """sign of re-sorting c^+_{p(i1)} c^+_{p(i2)} ... for every mask with npart bits."""
        N = self.N
        masks = np.array([sum(1 << i for i in c) for c in itertools.combinations(range(N), npart)],
                         dtype=np.int64)
        inv = np.zeros(len(masks), dtype=np.int64)
        for i in range(N):
            for j in range(i + 1, N):
                if p[i] > p[j]:
                    inv += ((masks >> i) & 1) & ((masks >> j) & 1)
        newmask = np.zeros(len(masks), dtype=np.int64)
        for i in range(N):
            newmask |= ((masks >> i) & 1) << p[i]
        tab_sign = np.zeros(1 << N, dtype=np.int64)
        tab_new = np.zeros(1 << N, dtype=np.int64)
        tab_sign[masks] = 1 - 2 * (inv & 1)
        tab_new[masks] = newmask
        return tab_sign, tab_new

    def apply_site_perm(self, p):
        """Return (new_index, sign) such that P|n> = sign[n] |new_index[n]>."""
        su, nu = self.perm_sign_table(p, self.nup)
        sd, nd = self.perm_sign_table(p, self.ndn)
        return self.index(nu[self.u], nd[self.d]), su[self.u] * sd[self.d]

    def op_site_perm(self, p):
        idx, sgn = self.apply_site_perm(p)
        return sp.coo_matrix((sgn.astype(float), (idx, np.arange(self.dim))),
                             shape=(self.dim, self.dim)).tocsr()

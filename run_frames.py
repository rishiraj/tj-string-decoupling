import sys, pickle, time, numpy as np, scipy.sparse.linalg as sla
from frames import *
J = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4
t0 = time.time()
H = model.hamiltonian(1.0, J)
E0, psi, ops = sector_ground_state(H)
Emin = sla.eigsh(H, k=1, which='SA', tol=1e-9)[0][0]
print(f"J/t={J}  (U/t ~ {4/J:.1f})   E0[K=0, d-wave sector] = {E0:.6f}   global min = {Emin:.6f}   E0/N = {E0/16:.5f}")
print("  quantum numbers: " + " ".join(f"{k}={psi@(O@psi):+.3f}" for k, O in ops.items()))
Lm = build_frame_map()
Hfr = (Lm.T @ H @ Lm).tocsr(); Gfr = (Lm.T @ Lm).tocsr()
pf = Lm.T @ psi
print(f"  frame map check: Rayleigh = {(pf@(Hfr@pf))/(pf@(Gfr@pf)):.6f}, |psi_f|^2 = {pf@pf:.4f}")
M = pf.reshape(15, NSP)
# hole-hole distribution
Pr = (M**2).sum(1); Pr /= Pr.sum()
d2 = {}
for r in range(1, N):
    x, y = xy(r); x = min(x, 4-x); y = min(y, 4-y); d2.setdefault(x*x+y*y, 0.0); d2[x*x+y*y] += Pr[r-1]
print("  pair distribution P(|r|^2): " + ", ".join(f"{k}:{v:.3f}" for k, v in sorted(d2.items())))
frames = {"raster (1D snake squeezed)": frame_transform(frame_labels("raster")),
          "site-anchored":              frame_transform(frame_labels("site")),
          "site-anchored + Marshall":   frame_transform(frame_labels("site"), True),
          "geometric string":           frame_transform(frame_labels("string")),
          "geometric string + phase-string sign": frame_transform(frame_labels("string"), True)}
res = {}
print(f"  {'frame':40s} {'lam1^2':>7s} {'lam2^2':>7s} {'S_ent':>6s} | {'E(k=1)':>9s} {'dE/N':>8s} {'fid':>6s} | {'E(k=2)':>9s} {'E(k=4)':>9s}")
for name, fr in frames.items():
    F = frame_operator(fr); Hf = (F @ Hfr @ F.T).tocsr(); Gf = (F @ Gfr @ F.T).tocsr(); Mf = to_frame(M, fr)
    p, S = schmidt(Mf)
    row = [p[0], p[1], S]; Ek = {}
    for k in (1, 2, 4):
        E, v = optimise_rank_k(Hf, Gf, Mf, k); Ek[k] = E
        if k == 1: fid = fidelity(v, Mf, Gf)
    print(f"  {name:40s} {p[0]:7.4f} {p[1]:7.4f} {S:6.3f} | {Ek[1]:9.5f} {(Ek[1]-E0)/16:8.5f} {fid:6.4f} | {Ek[2]:9.5f} {Ek[4]:9.5f}", flush=True)
    res[name] = dict(p=p, S=S, Ek=Ek, fid=fid)
pickle.dump(dict(J=J, E0=E0, res=res, Pr=Pr, M=M), open(f"frames_J{J}.pkl", "wb"))
print(f"  [{time.time()-t0:.0f}s]")

# String-length ansatz for two holes in the 4x4 t-J model: scripts

All scripts need numpy and scipy only. Run from this folder. Each finishes in under five minutes on one core with 3 GB of memory.

| Script | What it does |
|---|---|
| `tj_core.py` | t-J Hamiltonian with fermion signs, lattice symmetry operators |
| `validate.py` | Heisenberg benchmark; Hubbard vs t-J (with and without three-site terms) on a 10-site cluster at 20% doping |
| `frames.py`, `run_frames.py` | Schmidt analysis of "pair wavefunction x spin background" product states in fixed frames. `python run_frames.py 0.4` |
| `strings.py` | The string ansatz: sum_l a_l T^l Delta_d \|AF vacuum>, compared with full Lanczos from the same seed |
| `extras.py` | Sector dimension, staggered magnetisation, J-scaling of the pair energy |
| `excite.py` | Excitations: string model vs exact states resolved by symmetry and spin. `python excite.py 0.4` |
| `rvb.py` | Gutzwiller-projected d-wave pair (RVB) baseline. `python rvb.py 0.4` |
| `rvb_check.py`, `rvb_slater_check.py` | Checks of the RVB machinery at half filling |

## Known limitations

- The RVB baseline is mechanically correct (it matches an independently built projected Slater determinant), but on the periodic 4x4 cluster it is a poor state even at half filling (E/N = -0.39 against the exact -0.70). Its two-hole numbers are not a fair proxy for RVB quality on larger lattices.
- The string model reproduces the ground state but not the excited states (see `excite.py` output).
- Two holes on 16 sites cannot address long-range order.

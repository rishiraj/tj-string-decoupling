# A string-length representation of the two-hole ground state of the t–J model: pair–antiferromagnet decoupling tested by exact diagonalization

Nov 18, 2026 · @Rishiraj Acharya

[![DOI](https://zenodo.org/badge/1462319304.svg)](https://doi.org/10.5281/zenodo.24330930)

## Abstract

A nine-amplitude wavefunction reproduces the exact two-hole ground state of the 4×4 t–J model with 97% fidelity for 0.25 ≤ J/t ≤ 1. The wavefunction is a d-wave singlet pair removed from the undoped Heisenberg ground state and dressed by powers of the projected hopping operator. Its only variational freedom is the distribution of string lengths.

In this representation the antiferromagnetic background and the pair enter as separate factors. The energy error is 0.002–0.004t per site, and the suppression of the staggered magnetization by doping is reproduced without modifying the background. Product states of a pair wavefunction and a spin background in any fixed labelling of the spins fail, reaching at most 59% fidelity.

The representation has clear limits. It does not describe excited states: predicted gaps are up to 2.6 times too large, and the overlap of its first excited d-wave state with the exact one is 3% or less for J/t ≥ 0.4. The residual 3% of the ground state requires relaxation of the background. A two-hole, 16-site cluster cannot address long-range order, stripes or superconductivity, and the strong-coupling reduction from the Hubbard model costs 0.01–0.02t per site at U/t = 8–16. We report all negative results and release the code.

## 1. Introduction

The lightly hole-doped two-dimensional Hubbard model at strong coupling has no controlled analytic solution. It is the minimal model of the cuprate superconductors \[1–4\]. Large-scale numerical work finds striped and d-wave superconducting states in close competition \[5–7\]. Uniform and striped states differ by about 0.01t per site \[5\].

Two standard frameworks carry assumptions that fail in this regime. Fermi liquid theory assumes that the electron is the correct quasiparticle. Mean-field BCS theory assumes that pairs form from those quasiparticles. Both treat antiferromagnetism and superconductivity as separate order parameters defined in laboratory coordinates.

Several proposals replace these coordinates. They include resonating valence bonds \[1\], SO(5) superspin rotations \[8\], squeezed space in one dimension \[9, 25\] and phase strings \[24\]. The oldest is the string picture, in which a moving hole leaves a trail of displaced spins \[10–14\]. It has been developed into a quantitative theory of single holes and hole pairs and tested in cold-atom experiments \[15–18\].

This paper asks one narrow question. Is there a representation of the two-hole ground state in which the antiferromagnetic background and the d-wave pair enter as separate factors? We answer it by exact diagonalization on a 4×4 cluster.

We report three findings:

1. Product states of a pair wavefunction and a spin background fail in every fixed spin labelling we tried (Section 5).
2. A representation in string length succeeds for the ground state with nine amplitudes (Sections 6 and 7).
3. The same representation fails for excited states (Section 8).

We do not claim a solution of the doped Hubbard model. The ansatz is a Krylov-space form of the string picture \[12–18\], and its novelty is modest. The contribution is a quantitative test of how far pair–background factorization holds, including where it breaks.

## 2. Model and methods

### 2.1 Hamiltonians

The starting point is the repulsive Hubbard model on the square lattice with nearest-neighbour hopping t and on-site repulsion U.

```latex
H_{\mathrm{Hub}} = -t \sum_{\langle ij \rangle, \sigma} \left( c^{\dagger}_{i\sigma} c_{j\sigma} + \mathrm{h.c.} \right) + U \sum_i n_{i\uparrow} n_{i\downarrow}
```

For U/t ≫ 1 a canonical transformation removes doubly occupied sites \[19\]. To second order it gives the t–J model with J = 4t²/U, plus three-site hopping terms of the same order that vanish at half filling. P projects out double occupancy.

```latex
H_{tJ} = -t \sum_{\langle ij \rangle, \sigma} P \left( c^{\dagger}_{i\sigma} c_{j\sigma} + \mathrm{h.c.} \right) P + J \sum_{\langle ij \rangle} \left( \mathbf{S}_i \cdot \mathbf{S}_j - \tfrac{1}{4} n_i n_j \right) \equiv T + J\, H_J
```

All main results use the t–J model without three-site terms. Section 4 quantifies the error of this reduction against the Hubbard model.

### 2.2 Cluster and symmetry sector

The main cluster is the 4×4 square lattice with periodic boundaries and two holes, a doping of 1/8. The sector with zero total S\_z has 411,840 states. We obtain ground states by sparse Lanczos diagonalization \[21\].

The ground state has zero momentum, is a spin singlet, and is odd under 90° rotation, which is d\_{x²−y²} pair symmetry. This symmetry sector contains 469 singlets. At J/t = 0.4 its lowest energy equals the global minimum, −10.68392t.

The ground state is degenerate with states at momenta (π,0) and (0,π). We resolved this three-fold degeneracy explicitly at J/t = 0.2. It is an artifact of the 4×4 torus, which is equivalent to a 2×2×2×2 hypercube. We select the zero-momentum d-wave state throughout.

### 2.3 Validation

- The half-filled 4×4 cluster gives a Heisenberg energy of −0.701780 J per site, matching the known value.
- On a 10-site cluster, the hopping and exchange matrices agree element by element with operators derived numerically from the Hubbard model.
- Variational energies are reported as errors per site, (E − E₀)/16, and wavefunction quality as fidelity, the squared overlap with the exact ground state.

## 3. A mean-field decoupling identity

At mean-field level, d-wave pairing is unchanged in form by antiferromagnetic order. This motivates the search for a strong-coupling analogue. The result is analytic and was not checked numerically in this work.

Consider commensurate order with wavevector Q = (π,π) and staggered field m. Split the dispersion into parts odd and even under a shift by Q, ε\_k = ε⁻\_k + ε⁺\_k. The spin-density-wave quasiparticles α and β have energies E^±\_k and real coherence factors u\_k and v\_k.

```latex
c_{k\sigma} = u_k \alpha_{k\sigma} - \sigma v_k \beta_{k\sigma}, \qquad c_{k+Q,\sigma} = \sigma v_k \alpha_{k\sigma} + u_k \beta_{k\sigma}, \qquad E^{\pm}_k = \varepsilon^{+}_k \pm \sqrt{(\varepsilon^{-}_k)^2 + m^2}
```

Take a singlet gap that is odd under the shift, Δ\_{k+Q} = −Δ\_k, as for d\_{x²−y²} pairing. For k in the magnetic zone the pair operator then keeps its form, because the cross terms cancel and u² + v² = 1.

```latex
\Delta_k \left( c_{k\uparrow} c_{-k\downarrow} - c_{k+Q\uparrow} c_{-k-Q\downarrow} \right) = \Delta_k \left( \alpha_{k\uparrow} \alpha_{-k\downarrow} - \beta_{k\uparrow} \beta_{-k\downarrow} \right)
```

Pairing is therefore purely intraband with an undressed gap, and the excitation energies follow directly.

```latex
\omega^{\pm}_k = \sqrt{ \left( E^{\pm}_k - \mu \right)^2 + \Delta_k^2 }
```

The identity holds for any band structure. For a gap that is even under the shift the same algebra gives a factor u² − v² and interband terms, so the two orders do not separate. The picture of pairing between spin-density-wave quasiparticles goes back to the spin-bag proposal \[20\].

This is a weak-coupling statement. It ignores the constraint of no double occupancy, which dominates at U/t ≫ 1. The remaining sections test whether any comparable separation survives that constraint.

## 4. Accuracy of the strong-coupling reduction

The t–J model misses the Hubbard ground-state energy by 0.01–0.02t per site at U/t = 8–16 and 20% doping. This is the same size as the energy difference between competing states \[5\], so the reduction itself limits what any t–J result says about the Hubbard model.

We compared ground-state energies on a periodic 10-site cluster with 8 electrons. Three-site terms were obtained by projecting the second-order Hubbard operator numerically, with no hand-coded signs.

**Table 1.** Ground-state energy of the Hubbard model and the error per site of two effective models, 10 sites, 8 electrons.

| U/t | Hubbard energy (t) | t–J error per site (t) | t–J + three-site error per site (t) |
| --- | --- | --- | --- |
| 8 | −8.1171 | +0.0115 | −0.0630 |
| 12 | −7.3497 | +0.0234 | −0.0238 |
| 16 | −6.9380 | +0.0211 | −0.0120 |
| 24 | −6.5182 | +0.0097 | −0.0047 |
| 32 | −6.3832 | +0.0070 | −0.0020 |
| 48 | −6.2565 | +0.0048 | −0.0009 |
| 64 | −6.1928 | +0.0037 | −0.0005 |

The plain t–J error falls as 1/U at large U, a factor 1.9 between U/t = 32 and 64. With three-site terms it falls as 1/U², a factor 4.1 over the same range. These are the scalings expected from the order of the neglected terms \[19\].

Three-site terms make the energy worse for U/t ≤ 12. The t–J error is also not monotonic between U/t = 8 and 12. This is consistent with a change of ground-state symmetry on this cluster, which we did not resolve.

## 5. The two-hole ground state and the failure of fixed-frame product states

### 5.1 Hole–hole distribution

The two holes are most often found at distance √2, where a bare nearest-neighbour d-wave pair has a node. At J/t = 0.4 the probabilities by separation are 0.288 at distance 1, 0.351 at √2, 0.108 at 2, 0.219 at √5 and 0.034 at 2√2.

The weight on the diagonal cannot come from a rigid pair. It must come from hole motion that rearranges the spins, with the d-wave sign carried jointly by the holes and the background.

### 5.2 Product states in fixed frames

The simplest form of decoupling is a product of a pair wavefunction χ(r) and a spin wavefunction Φ. We use translation symmetry to place one hole at the origin, leaving the relative position r with 15 values. The 14 spins are listed in an order that may depend on r, which we call a frame, giving 3,432 spin configurations.

The ground state is then a 15 × 3,432 matrix. Its Schmidt spectrum bounds the fidelity of any product state in that frame. We also minimized the energy of symmetrized product states and of sums of k products directly.

We tested three labellings, two of them also with a sublattice sign:

- **Raster:** spins listed in row-by-row site order, the two-dimensional analogue of squeezed space \[9, 25\].
- **Site-anchored:** each spin keeps the label of its lattice site.
- **String:** spins are dragged one step along a shortest path of the hole from a reference position \[15\].
- **Signs:** a factor −1 for each down spin on the wrong sublattice for its label, the Marshall sign \[23\] or its string counterpart \[24\].

**Table 2.** Pair–background factorization of the exact ground state in fixed frames, J/t = 0.4. The maximum entanglement entropy is ln 15 = 2.71.

| Frame | Largest Schmidt weight | Entanglement entropy | Fidelity of best product state | Energy error per site, one product (t) | Energy error per site, four products (t) |
| --- | --- | --- | --- | --- | --- |
| Raster | 0.242 | 2.10 | 0.321 | 0.097 | 0.014 |
| Site-anchored | 0.397 | 1.86 | 0.586 | 0.074 | 0.017 |
| Site-anchored with Marshall sign | 0.208 | 2.42 | 0.115 | 0.122 | 0.040 |
| String | 0.389 | 1.79 | 0.560 | 0.084 | 0.014 |
| String with phase-string sign | 0.184 | 2.50 | 0.302 | 0.141 | 0.043 |

No frame factorizes the state. The best single product reaches 59% fidelity with an energy error of 0.074t per site. Four products still leave 0.014t per site.

The reason is visible in Section 5.1. A hole pair at separation √2 can be reached from a horizontal or a vertical nearest-neighbour pair, and the two routes enter with opposite d-wave sign. One spin labelling per hole position cannot represent this interference. The natural coordinate is the path, not its endpoint.

## 6. The string-length representation

### 6.1 Definition

We replace the hole position by the string length ℓ as the coordinate. The ansatz removes a d-wave singlet pair from the undoped ground state and applies powers of the projected hopping T.

```latex
|\Psi_n\rangle = \sum_{\ell=0}^{n} a_\ell \, T^{\ell} \, \Delta_d \, |\Phi_{\mathrm{AF}}\rangle, \qquad \Delta_d = \sum_i \left( \Delta_{i,i+\hat{x}} - \Delta_{i,i+\hat{y}} \right), \qquad \Delta_{ij} = c_{i\uparrow} c_{j\downarrow} - c_{i\downarrow} c_{j\uparrow}
```

Φ\_AF is the exact Heisenberg ground state of the undoped cluster and contains no free parameter. T^ℓ generates every pair of hole paths of total length ℓ, including retraced and looped paths, with their fermion signs. Path interference is therefore built in.

The three ingredients are separate. Antiferromagnetic correlations enter only through Φ\_AF, pairing symmetry only through Δ\_d, and their interplay only through the n + 1 amplitudes a\_ℓ.

We orthonormalize the states T^ℓ Δ\_d Φ\_AF into a basis |ℓ⟩ and diagonalize H in it. The basis does not depend on J, so one basis serves all couplings. Two matrix elements vanish exactly, ⟨ℓ|T|ℓ⟩ and ⟨ℓ|H\_J|ℓ+1⟩, because each hop changes the parity of the number of holes on one sublattice and H\_J conserves it.

### 6.2 Ground-state accuracy

With n = 8 the ansatz reaches 95–97% fidelity and an energy error of 0.002–0.004t per site for J/t between 0.2 and 0.8. Nine amplitudes describe a state in a sector of 469 singlets.

**Table 3.** Fidelity and energy error of the string ansatz. The last column is the error of Lanczos iteration with the full Hamiltonian from the same starting state and the same number of basis states.

| J/t | Bare pair fidelity, n = 0 | Fidelity, n = 4 | Fidelity, n = 8 | Energy error per site, n = 8 (t) | Full Lanczos error per site, n = 8 (t) |
| --- | --- | --- | --- | --- | --- |
| 0.2 | 0.069 | 0.753 | 0.951 | 0.0021 | 0.0013 |
| 0.3 | 0.113 | 0.875 | 0.967 | 0.0021 | 0.0005 |
| 0.4 | 0.153 | 0.923 | 0.968 | 0.0025 | 0.0002 |
| 0.5 | 0.193 | 0.945 | 0.968 | 0.0029 | 0.0001 |
| 0.8 | 0.304 | 0.968 | 0.971 | 0.0039 | 0.0001 |

The bare pair alone is a poor state, with 7–30% fidelity. Strings of length up to four recover most of the weight, and the gain stops near n = 8.

&#91;embedded content: Output of strings.py · 4×4 t–J cluster, two holes, J/t = 0.4\]

**Figure 1.** Energy error per site against the string cutoff n. Both methods start from the same bare d-wave pair and use n + 1 basis states.

The string ansatz tracks full Lanczos iteration up to n ≈ 4 and then saturates, while Lanczos keeps converging. The difference is relaxation of the spin background, which powers of T alone cannot produce. This relaxation accounts for the last 3% of the ground state.

The plateau fidelity is nearly independent of coupling, 0.965–0.973 for J/t from 0.25 to 1. At smaller J the n = 8 truncation is insufficient: at J/t = 0.1 the fidelity is 0.63 and the error 0.006t per site.

## 7. Effective string Hamiltonian, energy scaling and antiferromagnetic order

### 7.1 A one-dimensional effective model

In the string basis the two-hole problem becomes a particle on a half-line with nearly constant hopping and a confining potential. The Hamiltonian is H = T + J·V, where T is tridiagonal with elements b\_ℓ and V couples ℓ to ℓ and ℓ ± 2.

**Table 4.** Matrix elements of the effective model in the d-wave channel. The potential is measured from V₀ = −15.319 J; the undoped vacuum is at −19.228 J.

| ℓ | Hopping b\_ℓ, from ℓ to ℓ+1 (t) | Potential V\_ℓ − V₀ (J) | Exchange coupling, from ℓ to ℓ+2 (J) |
| --- | --- | --- | --- |
| 0 | 3.32 | 0.00 | −0.08 |
| 1 | 2.96 | 1.68 | 0.05 |
| 2 | 3.19 | 3.13 | −0.03 |
| 3 | 3.25 | 3.89 | 0.39 |
| 4 | 2.91 | 4.53 | 0.63 |
| 5 | 3.32 | 5.28 | 0.71 |
| 6 | 2.94 | 5.39 | 1.09 |
| 7 | 3.53 | 6.42 | 0.75 |
| 8 | 2.93 | 6.21 | 0.90 |

The hopping is 3.2 ± 0.2 t. The potential rises by about 1.1 J per unit length for ℓ ≤ 4 and then bends over, saturating near 7 J. Saturation is expected on 16 sites, where long strings have already disordered the whole cluster.

The model is close to tridiagonal only for short strings. The coupling between ℓ and ℓ + 2 is below 0.1 J for ℓ ≤ 2 and of order J beyond ℓ = 4.

### 7.2 Kinetic origin of the d-wave symmetry

The first hop of the pair is 5.5 times stronger in the d-wave channel than in the s-wave channel. The s-wave and d-wave starting states have the same norm and the same exchange energy. Their first hopping elements are b₀ = 3.32t for d-wave and 0.60t for s-wave.

Paths from horizontal and vertical pairs to the same final configuration add for d-wave and nearly cancel for s-wave. On this cluster the preference for d-wave pairing is therefore a kinetic interference effect, not an exchange effect.

### 7.3 Closed-form estimate

A continuum treatment with constant hopping b and linear potential σJℓ gives an Airy-function ground state \[10, 14\].

```latex
E_0 \simeq V_0 J - 2b + a_1 \, b^{1/3} (\sigma J)^{2/3}, \qquad a_1 = 2.338
```

With b = 3.16t and σ = 1.13 this gives −10.43t at J/t = 0.4, against the exact −10.68t, and −15.37t at J/t = 0.8, against −16.03t. The errors are 0.016t and 0.041t per site. The closed form is qualitative only, because the potential is not linear and the lattice is discrete.

### 7.4 Scaling with J

The energy of the pair relative to the undoped vacuum fits a + b·J^ν over ten couplings from J/t = 0.1 to 1. The exact energies give ν = 0.76 with a = −6.38t. The string ansatz at n = 8 gives ν = 0.79 with a = −6.22t.

The continuum string prediction is ν = 2/3, and a single hole on 32 sites gives ν = 0.72 by the same kind of fit \[27\]. An exponent between 2/3 and 1 is expected, because the energy also contains terms linear in J from broken bonds. The fit has three parameters and no error analysis, so ν is an effective exponent.

### 7.5 Suppression of antiferromagnetic order

Two holes cut the squared staggered magnetization roughly in half, and the ansatz reproduces this with an unmodified vacuum. We use m\_s² = 3⟨(Σ\_i (−1)^i S^z\_i / N)²⟩. The undoped value is 0.277.

**Table 5.** Squared staggered magnetization with two holes.

| J/t | Exact | String ansatz, n = 8 |
| --- | --- | --- |
| 0.2 | 0.098 | 0.111 |
| 0.4 | 0.130 | 0.138 |
| 0.8 | 0.157 | 0.166 |
| 1.0 | 0.165 | 0.174 |

The ansatz overestimates m\_s² by 6% for J/t ≥ 0.4 and by 13% at J/t = 0.2. In this representation the loss of order is geometric: strings move spins onto the wrong sublattice, while the correlations of the vacuum stay intact.

## 8. Excitations

The string model does not describe the excited states. We compared its higher eigenvalues with exact states resolved by momentum, point-group symmetry and total spin. All exact states listed are singlets at zero momentum.

**Table 6.** Excitation gaps above the d-wave ground state, exact and from the string model with 13 basis states. The overlap is between the string-model state and the corresponding exact state.

| J/t | s-wave pair gap, exact (t) | s-wave pair gap, string (t) | s-wave overlap | Next d-wave state gap, exact (t) | Next d-wave state gap, string (t) | Next d-wave overlap |
| --- | --- | --- | --- | --- | --- | --- |
| 0.2 | 0.174 | 0.309 | 0.41 | 0.263 | 0.575 | 0.33 |
| 0.4 | 0.626 | 0.877 | 0.30 | 0.533 | 1.266 | 0.03 |
| 0.8 | 1.510 | 1.635 | 0.82 | 0.849 | 2.196 | 0.004 |

The model overestimates the gap to the next d-wave state by a factor 2.2–2.6, and its excited state has little or no overlap with the exact one. The s-wave pair gap is within 8% only at J/t = 0.8.

The exact gap to the next d-wave state is 1.1–1.3 J at all three couplings. This scaling suggests a spin excitation of the background bound to the pair, which a frozen vacuum cannot contain. We did not test this interpretation directly.

The separation of pair and background therefore holds for the ground state only. A description of the low-energy spectrum must let the background respond.

## 9. A projected d-wave pair baseline, and why it is not a fair comparison here

We have no fair variational baseline on this cluster. The standard compact ansatz for this problem is the Gutzwiller-projected d-wave BCS state \[1, 22\]. Our implementation on the periodic 4×4 lattice gives a poor state even at half filling. Earlier exact work reports significant overlaps of optimized resonating-valence-bond states with exact states on clusters of this size \[26\]. The weakness therefore most likely lies in our construction, not in the class of states.

The difficulty is a shell effect. The point k = (0,0) is a filled nodal point, and the four points (±π/2, ±π/2) sit exactly at the Fermi level on the gap nodes. A pure d-wave pair function then gives a two-hole state with s-wave symmetry, orthogonal to the ground state.

We repaired the symmetry by giving the nodal Fermi points an s-like pair amplitude and projecting onto the ground-state sector. The projected state is then a polynomial in three amplitude ratios with 22 coupling-independent component vectors. We optimized over the BCS form and over free ratios.

- **Check of the construction:** the pair-determinant amplitudes equal those of an independently built projected Slater determinant to machine precision.
- **Half filling:** the best state has −0.39 J per site against the exact −0.70, and 1% fidelity.
- **Two holes, J/t = 0.4:** the best state has an energy error of 0.098t per site and 4% fidelity. The BCS form and the free-ratio form give the same optimum.

The string ansatz's error of 0.0025t per site should therefore not be read as an improvement over projected BCS states. A fair comparison needs a construction that reproduces the overlaps of Ref. \[26\], or a cluster or boundary condition with no nodal point at the Fermi level, for example a tilted cluster or twisted boundary conditions chosen to avoid such points.

## 10. Discussion

### 10.1 What the results establish

On this cluster the two-hole ground state is, to 97%, two holes joined by a fluctuating string in an unmodified antiferromagnet. Three observations support this reading.

- The background factor is the undoped ground state, with nothing adjusted.
- The loss of staggered order on doping follows from string geometry alone (Section 7.5).
- The d-wave symmetry is selected by interference in the first hop (Section 7.2).

In laboratory coordinates antiferromagnetism and pairing appear to compete, since adding the pair halves the staggered order. In string coordinates the same state is a product-like form in which neither factor is degraded. In this limited sense the two orders decouple.

### 10.2 What they do not establish

The map from the bare pair to the ground state is a polynomial in T, which is not unitary. The representation is a variational subspace, not a duality. It gives no transformed Hamiltonian in which the two sectors are exactly independent.

The separation fails for excited states (Section 8) and leaves 3% of the ground state unexplained. Both gaps point to the same missing ingredient, the dynamics of the background.

### 10.3 Relation to earlier work

The physical picture is not new. Strings were introduced for single holes in antiferromagnets \[10, 11, 13, 14\] and applied to hole pairs in variational calculations \[12\]. Recent work treats pairs of holes joined by geometric strings and predicts their binding and spectra \[17, 18\]. Exact diagonalization of two holes on the 4×4 cluster is also long established \[4, 26\].

Our construction is a Krylov space of the hopping operator from a symmetry-selected starting state. What it adds is quantitative: a direct comparison with exact states, a comparison with fixed-frame factorization, a coupling-independent effective Hamiltonian (Table 4), and an explicit account of where the picture fails.

### 10.4 A proposed experimental test

Quantum gas microscopes resolve spins and holes in single snapshots of the Fermi–Hubbard model and have already been used to look for string patterns in the doped model \[16\]. Two predictions of the present picture can be tested on snapshots that contain exactly two holes in a region of antiferromagnetic correlations.

1. The hole–hole distance distribution peaks at the diagonal separation √2, not at nearest neighbours.
2. Staggered spin correlations recover when the spins along a shortest path between the two holes are shifted back by one site.

Both predictions come from a 16-site ground state. Finite temperature and larger systems may change them, and we have not estimated the required signal-to-noise ratio.

## 11. Limitations

- **System size:** two holes on 16 sites cannot show long-range order, phase coherence or stripes. A period-8 stripe at 1/8 doping does not fit on the cluster. Exact diagonalization of four holes on 32 sites at J/t = 0.3 finds unbound holes, unlike the two-hole case \[28\], so two-hole results need not carry over to finite doping.
- **Single cluster:** there is no finite-size scaling. The ground-state degeneracy with (π,0) and (0,π) is specific to the 4×4 torus.
- **Scalability:** the ansatz uses the exact undoped ground state and exact applications of T. On large lattices both need approximations, such as a variational vacuum and sampling of strings. We did not test any.
- **Model:** we use the t–J model without three-site terms and without next-nearest-neighbour hopping. Section 4 shows an error of 0.01–0.02t per site against the Hubbard model at U/t = 8–16.
- **Baseline:** we did not obtain a fair comparison with projected BCS states (Section 9).
- **Excitations:** the representation fails for excited states (Section 8).
- **Statistics:** the exponent in Section 7.4 comes from a three-parameter fit without uncertainty estimates.
- **Analytic result:** the identity of Section 3 was derived by hand and not verified numerically.

## 12. Conclusion

The exact two-hole ground state of the 4×4 t–J model is captured to 97% by a d-wave pair removed from the undoped antiferromagnet and dressed by strings of length up to eight. Nine amplitudes suffice, the basis is independent of J, and the energy error is 0.002–0.004t per site.

In this representation the antiferromagnetic background and the pair are separate factors, and the apparent competition between them is a geometric effect of the strings. Product states in fixed spin labellings do not achieve this.

The separation is restricted to the ground state of a small cluster. Excited states, the last 3% of the ground state and every question about long-range order require the background to respond.

Three extensions would decide whether the picture is more than a small-cluster observation:

1. Add starting states built from excited states of the undoped antiferromagnet, and repeat the excitation test.
2. Repeat the calculation on clusters without nodal points at the Fermi level, with a fair projected-BCS baseline.
3. Replace the exact vacuum by a variational one and sample strings stochastically, to reach four or more holes and larger lattices.

## Code availability

All scripts are in the archive `hubbard_strings.zip`. They need only NumPy and SciPy, and each runs in under five minutes on one core with 4 GB of memory.

| Script | Produces |
| --- | --- |
| `tj_core.py` | Hamiltonian and symmetry operators |
| `validate.py` | Section 2.3 and Table 1 |
| `frames.py`, `run_frames.py` | Section 5 and Table 2 |
| `strings.py` | Table 3, Figure 1 and Table 4 |
| `extras.py` | Sector dimension, Section 7.4 and Table 5 |
| `excite.py` | Section 7.2, Section 7.3 and Table 6 |
| `rvb.py`, `rvb_check.py`, `rvb_slater_check.py` | Section 9 |

## Statement on the use of AI

The calculations, the code and the first draft of this manuscript were produced with an AI model (Claude, Anthropic) under the direction of the author. The author is responsible for checking the results and the references before submission.

## References

Entries 5, 7, 17, 18 and 26–28 were checked against the publisher or arXiv record on 18 November 2026 and are linked.

1. P. W. Anderson, The resonating valence bond state in La₂CuO₄ and superconductivity, Science 235, 1196 (1987).
2. F. C. Zhang and T. M. Rice, Effective Hamiltonian for the superconducting Cu oxides, Phys. Rev. B 37, 3759 (1988).
3. P. A. Lee, N. Nagaosa and X.-G. Wen, Doping a Mott insulator: physics of high-temperature superconductivity, Rev. Mod. Phys. 78, 17 (2006).
4. E. Dagotto, Correlated electrons in high-temperature superconductors, Rev. Mod. Phys. 66, 763 (1994).
5. B.-X. Zheng et al., [Stripe order in the underdoped region of the two-dimensional Hubbard model](https://doi.org/10.1126/science.aam7127), Science 358, 1155 (2017).
6. M. Qin et al., Absence of superconductivity in the pure two-dimensional Hubbard model, Phys. Rev. X 10, 031016 (2020).
7. H. Xu, C.-M. Chung, M. Qin, U. Schollwöck, S. R. White and S. Zhang, [Coexistence of superconductivity with partially filled stripes in the Hubbard model](https://doi.org/10.1126/science.adh7691), Science 384, eadh7691 (2024).
8. S.-C. Zhang, A unified theory based on SO(5) symmetry of superconductivity and antiferromagnetism, Science 275, 1089 (1997).
9. M. Ogata and H. Shiba, Bethe-ansatz wave function, momentum distribution, and spin correlation in the one-dimensional strongly correlated Hubbard model, Phys. Rev. B 41, 2326 (1990).
10. L. N. Bulaevskii, E. L. Nagaev and D. I. Khomskii, A new type of auto-localized state of a conduction electron in an antiferromagnetic semiconductor, Sov. Phys. JETP 27, 836 (1968).
11. W. F. Brinkman and T. M. Rice, Single-particle excitations in magnetic insulators, Phys. Rev. B 2, 1324 (1970).
12. S. A. Trugman, Interaction of holes in a Hubbard antiferromagnet and high-temperature superconductivity, Phys. Rev. B 37, 1597 (1988).
13. B. I. Shraiman and E. D. Siggia, Mobile vacancies in a quantum Heisenberg antiferromagnet, Phys. Rev. Lett. 61, 467 (1988).
14. C. L. Kane, P. A. Lee and N. Read, Motion of a single hole in a quantum antiferromagnet, Phys. Rev. B 39, 6880 (1989).
15. F. Grusdt et al., Parton theory of magnetic polarons: mesonic resonances and signatures in dynamics, Phys. Rev. X 8, 011046 (2018).
16. C. S. Chiu et al., String patterns in the doped Hubbard model, Science 365, 251 (2019).
17. F. Grusdt, E. Demler and A. Bohrdt, [Pairing of holes by confining strings in antiferromagnets](https://doi.org/10.21468/SciPostPhys.14.5.090), SciPost Phys. 14, 090 (2023).
18. A. Bohrdt, E. Demler and F. Grusdt, [Dichotomy of heavy and light pairs of holes in the t–J model](https://doi.org/10.1038/s41467-023-43453-2), Nat. Commun. 14, 8017 (2023).
19. A. H. MacDonald, S. M. Girvin and D. Yoshioka, t/U expansion for the Hubbard model, Phys. Rev. B 37, 9753 (1988).
20. J. R. Schrieffer, X.-G. Wen and S.-C. Zhang, Dynamic spin fluctuations and the bag mechanism of high-T\_c superconductivity, Phys. Rev. B 39, 11663 (1989).
21. C. Lanczos, An iteration method for the solution of the eigenvalue problem of linear differential and integral operators, J. Res. Natl. Bur. Stand. 45, 255 (1950).
22. C. Gros, Superconductivity in correlated wave functions, Phys. Rev. B 38, 931 (1988).
23. W. Marshall, Antiferromagnetism, Proc. R. Soc. Lond. A 232, 48 (1955).
24. D. N. Sheng, Y. C. Chen and Z. Y. Weng, Phase string effect in a doped antiferromagnet, Phys. Rev. Lett. 77, 5102 (1996).
25. H. V. Kruis, I. P. McCulloch, Z. Nussinov and J. Zaanen, Geometry and the hidden order of Luttinger liquids: the universality of squeezed space, Phys. Rev. B 70, 075109 (2004).
26. Y. Hasegawa and D. Poilblanc, [Hole dynamics in the t–J model: an exact diagonalization study](https://www.osti.gov/biblio/5283545), Phys. Rev. B 40, 9035 (1989).
27. P. W. Leung and R. J. Gooding, Dynamical properties of the single-hole t–J model on a 32-site square lattice, [arXiv:cond-mat/9509087](https://arxiv.org/pdf/cond-mat/9509087) (1995).
28. P. W. Leung, Hole correlation and antiferromagnetic order in the t–J model, [arXiv:cond-mat/0511568](https://arxiv.org/pdf/cond-mat/0511568) (2005).

## Citation

If you build upon this analysis or use the exact-diagonalization scripts, please cite:

```bibtex
@misc{acharya2026string,
  author       = {Acharya, Rishiraj},
  title        = {A string-length representation of the two-hole ground state of the t--J model: pair--antiferromagnet decoupling tested by exact diagonalization},
  month        = nov,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.24330930},
  url          = {[https://doi.org/10.5281/zenodo.24330930](https://doi.org/10.5281/zenodo.24330930)}
}

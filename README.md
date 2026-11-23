# Normal-Flux DEC Surface MHD on the Icosahedral Sphere

A self-contained **Fortran 2018** prototype simulating compressible ideal magnetohydrodynamics (MHD) coupled to a normal magnetic flux on a closed triangulated 2-sphere (`S²`) using **Discrete Exterior Calculus (DEC)** and midpoint transport.

> **Prototype Disclaimer:** This repository is an experimental numerical prototype, not a validated energy-preserving variational integrator. All code and related mathematics was created locally, no AI generated content was used to help with this project involving mathematics or code.

---

## 1. Core Idea: Why Discrete Exterior Calculus for MHD?

Standard finite-difference or finite-volume fluid solvers discretize vector components on grid cells, which gradually smears out topological invariants like magnetic divergence (`∇ · B = 0`) over long runs. 

Instead of storing vectors at grid centers, **Discrete Exterior Calculus (DEC)** assigns physical quantities directly to the geometric elements of a triangulated mesh (`K`) according to their differential form degree:

| Mesh Element | Cochain Space | Physical Quantity | Fortran Array |
| :--- | :--- | :--- | :--- |
| **Vertices (`V`)** | 0-cochain `C⁰(K)` | Projected scalar potential (`q_v`), dual vertex areas | `q_v(nv)`, `star0(nv)` |
| **Edges (`E`)** | 1-cochain `C¹(K)` | Fluid velocity circulation (`u_e`), magnetic vector potential (`A_e`), fluxes | `u(ne)`, `A(ne)`, `flux_m(ne)` |
| **Faces (`F`)** | 2-cochain `C²(K)` | Triangle mass (`m_f`), normal magnetic flux (`B_f`), triangle area (`A_f`) | `m(nf)`, `B(nf)`, `area(nf)` |

Because the magnetic flux `B` is defined as the discrete exterior derivative of the edge potential (`B = d₁ A`) and `d₁ · d₀ = 0` holds identically as an integer matrix product, **net magnetic flux stays zero to machine precision (`~1e-16`) at every time step without divergence cleaning.**

---

## 2. Repository Layout

    ├── mhd_variational_integrator.f90   # Mesh builder, DEC operators, Hodge stars, and midpoint solver
    ├── main_sphere_alfven.f90           # Initialization, time-stepping driver, and diagnostic output
    ├── main_tex.pdf                     # Full mathematical write-up (compiled PDF)
    ├── .gitignore                       # Ignores compiled object, module, and executable binaries
    ├── LICENSE                          # MIT License
    └── README.md

---

## 3. Code Architecture & Walkthrough

### A. `mhd_variational_integrator.f90`

This module contains the entire geometric engine and time-stepping solver in contiguous, `allocatable` Fortran 2018 arrays with zero external dependencies:

1. **Icosahedral Sphere Construction (`build_icosphere`)**
   - Generates the 12 golden-ratio vertices (`ϕ = (1 + √5) / 2`) normalized onto the unit sphere (`|V| = 12`), 30 undirected edges (`|E| = 30`), and 20 outward-oriented triangular faces (`|F| = 20`), satisfying Euler's formula `|V| - |E| + |F| = 2`.
   - Precomputes face areas (`area`), unit outward normals (`normal`), edge lengths (`edge_len`), unit edge tangents (`edge_tangent`), and signed incidence arrays (`d0` mapping vertices to edges, `d1` mapping edges to faces).

2. **Metric Hodge Star Operators (`build_hodge_stars`)**
   - **`star0` (`⋆₀`):** Barycentric dual cell area at each vertex (`∑ A_f / 3` over incident faces).
   - **`star1_bary` (`⋆₁ᵇᵃʳʸ`):** Diagonal barycentric edge Hodge star (`(1 / ℓ_e²) ∑ A_f / 3`), mapping primal velocity 1-cochains to dual momentum 1-cochains.
   - **`star1_circ` (`⋆₁ᶜⁱʳᶜ`):** Signed circumcentric dual-to-primal edge ratio (`ℓ⋆_e / ℓ_e`), computed from adjacent triangle circumcenters for metric diagnostics.
   - **`star2` (`⋆₂`):** Inverse triangle area (`1 / A_f`), converting integrated face quantities (`m_f`, `B_f`) into intensive densities (`ρ_f = m_f / A_f`, `b_f = B_f / A_f`).

3. **Discrete Lie Advection & Interior Product (`edge_contraction`)**
   - Computes the discrete interior product `ι_u B` and `ι_u m` across each primal edge `e` shared by faces `f₊` and `f₋`:
     `(ι_u B)_e = 0.5 * u_e * ℓ_e * (B_{f+} / A_{f+} + B_{f-} / A_{f-})`
   - Applying `d₁` to `ι_u m` yields conservative face-to-face mass transport (`dm/dt = -d₁(ι_u m)`), while updating `dA/dt = -ι_u B` advances the vector potential along the flow.

4. **Momentum Residual & Transverse Velocity Reconstruction (`compute_momentum_rhs`)**
   - Reconstructs a 3D tangent velocity vector `v_f` on each triangle from its three boundary edge circulations, then projects the adjacent face velocities onto the in-surface edge normal to obtain the transverse velocity `u⊥_e` needed for the nonlinear convective/vorticity term.
   - Combines barotropic thermal enthalpy (`c_s² ln(ρ_f / ρ_ref)`) and normal magnetic pressure (`b_f² / (2 μ₀)`) on faces, projects that scalar onto vertices via area-weighted averaging (`q_v`), and takes `-d₀ q_v` to compute pressure and Lorentz forces along edges.

5. **Damped Fixed-Point Midpoint Stepper (`advance_step`)**
   - Advances `(uⁿ, mⁿ, Aⁿ) → (uⁿ⁺¹, mⁿ⁺¹, Aⁿ⁺¹)` using midpoint states `uⁿ⁺¹/² = 0.5 * (uⁿ + uⁿ⁺¹)`.
   - Iterates up to **80 damped Picard steps** (damping factor `ω = 0.5`) until the edge velocity update difference drops below `1.0e-10`.
   - If the iteration does not converge within 80 steps, it falls back cleanly to a single explicit step (`un = u - dt * rhs / star1_bary`) and flags the unconverged iteration count in the diagnostics.

---

### B. `main_sphere_alfven.f90`

The main driver initializes a smooth spherical harmonic-like perturbation on the unit icosahedron and integrates it forward in time:

- **Initial Conditions:** Sets uniform reference density `ρ_ref = 1.0`, sound speed `c_s = 0.2`, magnetic permeability `μ₀ = 1.0`, a divergence-free tangent velocity perturbation (`u_e`) from a vertex streamfunction, and a non-trivial edge magnetic potential (`A_e`) so `B = d₁ A` starts non-zero with exact zero net flux (`∑ B_f = 0`).
- **Diagnostic Monitoring:** At every 10th step (`nout = 10` up to `nsteps = 100` with `dt = 1.0e-3`), prints:
  - `exact_B_error` (`||B - d₁ A||_∞`) — Exact closure error (verifies machine-zero topological constraint).
  - `total_flux` (`∑ B_f`) and `total_mass` (`∑ m_f`) — Global magnetic flux and total mass conservation.
  - `kinetic`, `magnetic`, `internal` — Kinetic, normal magnetic, and isothermal free energy breakdown.
  - `iterations` and `residual` — Nonlinear midpoint solver convergence health.

---

## 4. Build & Run Instructions

Requires `gfortran` (GCC 9+ with Fortran 2018 support). From the root directory:

    gfortran -std=f2018 -O2 -Wall -Wextra -fcheck=all \
      mhd_variational_integrator.f90 main_sphere_alfven.f90 \
      -o sphere_mhd

    ./sphere_mhd

---

## 5. Known Limitations & Next Steps

- **Coarse Base Mesh:** Currently runs on the unsubdivided 20-face icosahedron; adding recursive midpoint/Loop geodesic refinement (`80`, `320`, `1280` faces) is the next step for spatial convergence studies.
- **Normal-Flux Physics Only:** Because `B ∈ C²(K)` is a primal 2-form on a 2D surface, it captures magnetic pressure (`b_f² / 2μ₀`) normal to `S²` rather than tangential Alfvén tension along the surface.
- **Diagonal Hodge Stars:** Uses diagonal barycentric `⋆₁` approximations rather than full Galerkin/Whitney 1-form mass matrices.
- **Midpoint Picard Solver:** Uses a practical damped fixed-point midpoint update with an explicit fallback rather than a full Jacobian Newton-Krylov DEL solve.

---

## 6. References & Foundational Literature

The discrete differential geometry and variational fluid transport structures used in this prototype build on the following foundational papers:

1. **Discrete Exterior Calculus (DEC) Foundations**
   - Hirani, A. N. (2003). *Discrete Exterior Calculus*. Ph.D. Thesis, California Institute of Technology.
   - Desbrun, M., Hirani, A. N., Leok, M., & Marsden, J. E. (2005). *Discrete Exterior Calculus*. arXiv:math/0508341.
   *(Defines primal simplicial cochains, discrete coboundary operators `d₀, d₁`, diagonal barycentric/circumcentric Hodge stars `⋆ₖ`, and discrete interior products `ι_u`.)*

2. **Variational Integrators & Discrete Euler-Poincaré Fluid Mechanics**
   - Pavlov, D., Mullen, P., Tong, Y., Kanso, E., Marsden, J. E., & Desbrun, M. (2011). *Structure-Preserving Discretization of Incompressible Fluids*. Physica D: Nonlinear Phenomena, 240(6), 443–458.
   - Gawlik, E. S., Mullen, P., Pavlov, D., Marsden, J. E., & Desbrun, M. (2011). *Geometric, Variational Discretization of Continuum Theories*. Physica D: Nonlinear Phenomena, 240(21), 1724–1760.
   *(Introduces discrete Lie advection of differential forms and variational Euler-Poincaré updates on simplicial meshes for compressible and incompressible fluids.)*

3. **Structure-Preserving & Exterior-Calculus Magnetohydrodynamics (MHD)**
   - Gawlik, E. S., & Gay-Balmaz, F. (2021). *A Variational Finite Element Discretization of Compressible Magnetohydrodynamics*. Foundations of Computational Mathematics, 21, 961–1001.
   - Kraus, M., Kormann, K., Morrison, P. J., & Sonnendrücker, E. (2017). *GEMPIC: Geometric Electromagnetic Particle-in-Cell Methods*. Journal of Plasma Physics, 83(4).
   *(Establishes exact magnetic flux preservation via potential formulations `B = d₁ A` and variational coupling between fluid advection and Lorentz forces.)*

## Credits
- lol me

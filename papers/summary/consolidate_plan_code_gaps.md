## Part A — Gaps Identified in Stage 1 (What Must Be Addressed in Stage 2)

Having reviewed the consolidated summary thoroughly, here are the structural and algorithmic gaps that a "golden reference pseudo-codebase" must close:

**Architectural gaps:**
The Stage 1 Algorithm 1 is monolithic — a single 227-line procedure mixing grid construction, operator assembly, time-stepping logic, linear algebra, boundary handling, monitoring events, and Greek extraction. A production-quality reference needs modular decomposition with explicit interfaces, data structures, and swappable strategy components.

**Missing algorithmic modules:**

| Module | Stage 1 status | What's needed |
|---|---|---|
| TR-BDF2 time stepper | Mentioned, not specified | Full two-stage pseudocode with coefficient formulas |
| American constraint (LCP) | Mentioned as "hook" | Policy iteration + penalty method inner loops |
| Non-uniform / graded meshes | Mentioned | Sinh-transform formula, mesh-grading algorithm, stencil coefficients on non-uniform grids |
| Keller box (simultaneous Delta) | Mentioned by Duffy | Staggered-grid setup, system assembly, extraction |
| All Greeks beyond Δ,Γ | Absent | Theta (from PDE residual), Vega (bump-and-reprice or sensitivity PDE), Rho |
| Local vol interpolation | Mentioned as "hook" | Surface lookup/interpolation inside operator assembly |
| Error / convergence monitor | Absent | Richardson extrapolation, oscillation detection, positivity check |
| Adaptive time-stepping | Mentioned (Phase 2) | Step-size controller with damping-aware logic |
| Multi-product dispatch | Absent | How the same engine serves calls, puts, digitals, barriers, Americans |

**Robustness infrastructure gaps:**

Stage 1 Algorithm 1 has no runtime diagnostics — no M-matrix verification, no oscillation detection, no convergence monitoring. A production reference must include these as first-class modules, not afterthoughts.

**Data structure gaps:**

Stage 1 works with raw arrays. A reference pseudo-codebase needs named data containers (grid objects, operator objects, solution snapshots, Greek bundles) so that the module interfaces are self-documenting.

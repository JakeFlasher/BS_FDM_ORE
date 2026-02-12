*Thinking...*

> The user wants me to draft a multi-round prompt strategy for implementing the Milev-Tagliani Crank-Nicolson variant scheme from their paper into QuantLib v1.23. Let me analyze what needs to be done:
> 
> 1. The paper describes two nonstandard FD schemes for Black-Scholes with discontinuous payoffs and low volatility:
>    - Duffy's exponentially fitted implicit scheme
>    - Milev-Tagliani's Crank-Nicolson variant
> 
> 2. The key mathematical elements from the paper:
>    - The CN variant modifies the reaction term (-rV) discretization using a 6-node stencil

# Multi-Round Implementation Prompt Strategy

Below is a complete, ordered prompt sequence. Each round targets a specific code module, passes forward any earlier modifications, and stays within context limits.

---

## Execution Order & Rationale

```
Round 1 ─► schemes/          (new scheme classes — zero dependencies on other changes)
Round 2 ─► operators/        (modified spatial operators — depends on Round 1 interfaces)
Round 3 ─► solvers/          (scheme registration — depends on Rounds 1–2)
Round 4 ─► root FD files     (old framework parallel implementation — independent but informed by 1–3)
Round 5 ─► pricingengines/   (user-facing engines — depends on Rounds 1–3)
Round 6 ─► utilities/ + meshers/ + stepconditions/ (supporting infrastructure — depends on all above)
```

---

## ROUND 1 PROMPT — `schemes/` directory

```markdown
# Task

You are implementing two nonstandard finite-difference schemes from the paper
"Low Volatility Options and Numerical Diffusion of Finite Difference Schemes"
(Milev & Tagliani, 2010) into QuantLib v1.23's new finite-difference framework.

The two schemes are:
1. **Duffy's implicit exponentially fitted scheme** (paper Section 2)
2. **Milev-Tagliani Crank-Nicolson variant** (paper Section 3)

In this round you will create NEW scheme classes inside the `schemes/` directory
only. You must output each changed or new file **in its entirety** — no diffs,
no "rest unchanged" placeholders.

## Mathematical specifications to implement exactly

### Scheme A — Exponentially Fitted Implicit (Duffy)

The fitted operator uses the standard implicit time discretization but replaces
the diffusion coefficient with a fitting factor:

$$
\rho_j^{n+1} = \frac{\mu_j^{n+1} h}{2}
  \coth\!\left(\frac{\mu_j^{n+1} h}{2\,\sigma_j^{n+1}}\right)
$$

where for Black-Scholes: μ(S,t) = rS, σ(S,t) = ½σ²S².

The resulting system is A U^{n+1} = U^n with tridiagonal A having entries:

- a_{j,j-1} = -k(ρ_j/h² - μ_j/(2h))
- a_{j,j}   = k(2ρ_j/h² - b_j + 1/k)
- a_{j,j+1} = -k(ρ_j/h² + μ_j/(2h))

The scheme is fully implicit, uniformly stable, oscillation-free, with
convergence |V(S_j,t_n) - U_j^n| ≤ c(h+k), c independent of σ.

In the limit σ→0, the fitting factor becomes |μh/2|, recovering first-order
upwind. The artificial diffusion is (1/2)μ·h·V_{SS} = (1/2)rSΔS·V_{SS}.

### Scheme B — Crank-Nicolson Variant (Milev-Tagliani)

This is standard Crank-Nicolson EXCEPT for the discretization of the reaction
term -rV. The mid-time value V(t + Δt/2) is replaced by a 6-node stencil:

$$
V\!\bigl(t+\tfrac{\Delta t}{2}\bigr)
= \omega_1(U^n_{j-1}+U^n_{j+1})
+ \bigl(\tfrac{1}{2}-2\omega_1\bigr)U^n_j
+ \omega_2(U^{n+1}_{j-1}+U^{n+1}_{j+1})
+ \bigl(\tfrac{1}{2}-2\omega_2\bigr)U^{n+1}_j
$$

The system is P U^{n+1} = N U^n where:

P = tridiag{
  rω₂ + (r/4)(S_j/ΔS) - (σS_j/(2ΔS))² ;
  1/Δt + (1/2)(σS_j/ΔS)² + r(1/2 - 2ω₂) ;
  rω₂ - (r/4)(S_j/ΔS) - (σS_j/(2ΔS))²
}

N = tridiag{
  -rω₁ - (r/4)(S_j/ΔS) + (σS_j/(2ΔS))² ;
  1/Δt - (1/2)(σS_j/ΔS)² - r(1/2 - 2ω₁) ;
  -rω₁ + (r/4)(S_j/ΔS) + (σS_j/(2ΔS))²
}

Parameters ensuring positivity (M-matrix for P, non-negative N):

$$
\omega_1 = \omega_2 = -\frac{r}{16\sigma^2}
$$

Time step restriction:

$$
\Delta t < \frac{1}{r(\frac{1}{2}-2\omega_1) + \frac{1}{2}(\sigma M)^2}
$$

where M is the number of spatial nodes.

The scheme has O(ΔS², Δt²) accuracy. Its artificial diffusion for small σ is:
(1/8)(r/σ · ΔS)² · V_{SS}.

## What to produce

Create the following new files, plus modify `all.hpp` to include them:

1. **`schemes/exponentiallyfittedscheme.hpp`** and
   **`schemes/exponentiallyfittedscheme.cpp`**
   — Implements Duffy's exponentially fitted implicit scheme.
   Must follow the same class structure as `ImplicitEulerScheme` (same typedefs,
   same method signatures: `step()`, `setStep()`).
   Must accept an `ext::shared_ptr<FdmLinearOpComposite>` and a `bc_set`.
   The fitting factor computation must be done inside `step()` based on the
   operator's mesh locations and the current local drift/diffusion values.

2. **`schemes/milevtaglianischeme.hpp`** and
   **`schemes/milevtaglianischeme.cpp`**
   — Implements the CN variant with modified reaction term.
   Must accept additional parameters: risk-free rate `r`, volatility `σ`,
   and spatial step `ΔS` (or equivalently the mesher) in its constructor.
   The ω₁, ω₂ parameters are computed internally.
   Must validate the time-step restriction in `setStep()`.
   Must follow the same class interface as `CrankNicolsonScheme`.

3. **Updated `schemes/all.hpp`** — Include the two new headers.

## Existing code provided

[ATTACH: schemes/cranknicolsonscheme.hpp, cranknicolsonscheme.cpp,
         impliciteulerscheme.hpp, impliciteulerscheme.cpp,
         expliciteulerscheme.hpp, expliciteulerscheme.cpp,
         boundaryconditionschemehelper.hpp,
         douglasscheme.hpp, douglasscheme.cpp,
         all.hpp]
[ATTACH: operators/fdmlinearopcomposite.hpp, triplebandlinearop.hpp]
[ATTACH: the paper PDF and paper summary sections 2–3]

## Constraints

- Output COMPLETE file contents for every file you create or modify.
- Use the exact QuantLib coding style: `ext::shared_ptr`, `Disposable<>`,
  `QL_REQUIRE`, `Real`, `Size`, `Time`, copyright headers, include guards.
- The schemes must work as drop-in replacements in `FdmBackwardSolver::rollback()`.
- Do NOT modify any files outside the `schemes/` directory in this round.
```

---

## ROUND 2 PROMPT — `operators/` directory

```markdown
# Task — Round 2: Operators

You are continuing the implementation of two nonstandard FD schemes from
Milev & Tagliani (2010). In Round 1 we created:
- `ExponentiallyFittedScheme` in schemes/
- `MilevTaglianiScheme` in schemes/

In this round you modify the `operators/` directory to support these schemes.

## What needs to change

### For the exponentially fitted scheme

The fitting factor ρ replaces the standard second-derivative coefficient.
Create a new operator class `FdmFittedBlackScholesOp` that:

1. In `setTime(t1, t2)`, computes μ_j = rS_j and σ_j = ½σ²S_j² at each node
2. Computes ρ_j = (μ_j·h/2)·coth(μ_j·h / (2·σ_j))
3. Assembles the tridiagonal system using ρ_j instead of σ_j for the
   second-derivative term, and uses the standard centered first-derivative
4. For σ_j → 0 (or |μ_j·h/(2σ_j)| > threshold), switches to the upwind
   limit ρ_j = |μ_j·h|/2

### For the Milev-Tagliani CN variant

Create `FdmMilevTaglianiBlackScholesOp` that:

1. Stores r, σ, and computes ω₁ = ω₂ = -r/(16σ²)
2. Provides methods `applyP()` and `applyN()` that compute P·x and N·x
   using the specific tridiagonal matrices from the paper
3. Provides `solveP()` that solves P·x = b using the Thomas algorithm
4. The P matrix must be verified as an M-matrix (positive diagonal,
   non-positive off-diagonals, diagonal dominance)

Alternatively, if modifying `FdmBlackScholesOp` is cleaner, add a mode
flag or subclass. The key requirement is that the reaction term -rV is
discretized using the 6-node stencil, not the standard CN approach.

## Files to produce (complete contents)

1. **`operators/fdmfittedblackscholesop.hpp`** and `.cpp` — New fitted BS operator
2. **`operators/fdmmilevtaglianiblackscholesop.hpp`** and `.cpp` — New CN variant operator
3. **Updated `operators/all.hpp`**

You may also need to modify:
4. **`operators/triplebandlinearop.hpp`** and `.cpp` — if the existing
   `solve_splitting()` / `axpyb()` / `apply()` methods are insufficient
   for the new stencil structures.

## Existing code provided

[ATTACH: All operators/*.hpp and *.cpp files from the provided codebase]
[ATTACH: The Round 1 output (new scheme files) for interface reference]
[ATTACH: Paper sections 2–3 and the paper summary sections 3.2–3.3]

## Constraints

- The new operators must implement the `FdmLinearOpComposite` interface.
- They must work with `FdmMesherComposite` and `FdmBlackScholesMesher`.
- The fitting factor computation must handle edge cases:
  when μ_j = 0, ρ_j = σ_j (standard centered);
  when σ_j → 0, use the upwind limit.
- Use `coth(x) = cosh(x)/sinh(x)` or equivalently `(exp(2x)+1)/(exp(2x)-1)`.
- Output COMPLETE file contents. No diffs.
```

---

## ROUND 3 PROMPT — `solvers/` directory

```markdown
# Task — Round 3: Solvers

Continuing the implementation. Rounds 1–2 produced:
- `ExponentiallyFittedScheme` and `MilevTaglianiScheme` (schemes/)
- `FdmFittedBlackScholesOp` and `FdmMilevTaglianiBlackScholesOp` (operators/)

In this round you modify the `solvers/` directory to register and use these
new schemes.

## What needs to change

### `fdmbackwardsolver.hpp` / `fdmbackwardsolver.cpp`

1. Add two new enum values to `FdmSchemeDesc::FdmSchemeType`:
   - `ExponentiallyFittedType`
   - `MilevTaglianiType`

2. Add static factory methods:
   - `FdmSchemeDesc::ExponentiallyFitted()`
   - `FdmSchemeDesc::MilevTagliani()`

3. In `FdmBackwardSolver::rollback()`, add switch cases for the two new types.
   For `ExponentiallyFittedType`: create `ExponentiallyFittedScheme`,
   wrap in `FiniteDifferenceModel<>`, call `rollback()`.
   For `MilevTaglianiType`: create `MilevTaglianiScheme`,
   wrap in `FiniteDifferenceModel<>`, call `rollback()`.

4. The `MilevTaglianiType` should use the Milev-Tagliani operator
   instead of the standard `FdmBlackScholesOp`. This may require the
   solver to create the operator differently. Consider whether
   `FdmBackwardSolver` needs to accept an alternative operator factory,
   or whether the scheme itself should encapsulate operator creation.

### `fdmblackscholessolver.hpp` / `fdmblackscholessolver.cpp`

Modify `performCalculations()` so that when the scheme type is
`MilevTaglianiType` or `ExponentiallyFittedType`, it creates the
appropriate operator instead of `FdmBlackScholesOp`.

### `fdmsolverdesc.hpp`

If additional fields are needed in `FdmSolverDesc` for the new schemes
(e.g., a flag for which operator to use), add them here.

## Files to produce (complete contents)

1. **`solvers/fdmbackwardsolver.hpp`** — Updated with new scheme types
2. **`solvers/fdmbackwardsolver.cpp`** — Updated with new switch cases
3. **`solvers/fdmblackscholessolver.hpp`** — Updated if needed
4. **`solvers/fdmblackscholessolver.cpp`** — Updated with operator selection
5. **`solvers/fdmsolverdesc.hpp`** — Updated if needed
6. **`solvers/fdm1dimsolver.cpp`** — Updated if the solver chain changes

## Existing code provided

[ATTACH: All solvers/*.hpp and *.cpp files]
[ATTACH: Round 1 output (scheme files)]
[ATTACH: Round 2 output (operator files)]
[ATTACH: Paper summary sections 3.2–3.3, 4.1–4.2 (algorithms)]

## Constraints

- The new schemes must be selectable via `FdmSchemeDesc` without
  modifying any code outside the `solvers/` and previously modified directories.
- Damping steps: For `ExponentiallyFittedType`, no separate damping is
  needed (the scheme is already non-oscillatory). For `MilevTaglianiType`,
  the scheme itself is non-oscillatory, so damping is optional.
- Backward compatibility: existing scheme types must continue to work
  identically.
- Output COMPLETE file contents.
```

---

## ROUND 4 PROMPT — Root-level `methods/finitedifferences/` files

```markdown
# Task — Round 4: Old Framework (Root-Level FD Files)

QuantLib v1.23 maintains a legacy finite-difference framework used by older
engines (e.g., `FDVanillaEngine`). This framework uses `TridiagonalOperator`,
`MixedScheme<>`, and `FiniteDifferenceModel<>`.

In this round, implement the same two schemes in the OLD framework for
backward compatibility.

## What to create/modify

### New file: `exponentiallyfittedscheme_old.hpp` (or similar name)

A template class analogous to `CrankNicolson<Operator>` but implementing
the exponentially fitted implicit scheme. It should:

1. Accept the same `operator_type` and `bc_set` as `ImplicitEuler<>`
2. In `setStep(dt)`, build the fitted operator using the fitting factor
3. In `step(a, t)`, solve the implicit system A·U^{n+1} = U^n

The fitting factor requires access to the grid (from the operator) and
the drift/diffusion coefficients. The old `BSMOperator` stores these
in its tridiagonal entries, so the fitted scheme may need to reconstruct
μ and σ from the operator's lower/mid/upper diagonals, or accept them
as constructor parameters.

### New file: `milevtaglianischeme_old.hpp`

A template class implementing the CN variant. Must:

1. Build the P and N matrices from the paper's formulas
2. In `step(a, t)`: compute N·a, then solve P·x = N·a
3. The ω parameters require r and σ, which must be passed to the constructor

### Modified: `bsmoperator.hpp` / `bsmoperator.cpp`

Add a new constructor or factory method that creates a fitted BSM operator
with the exponential fitting factor. The existing constructor uses:
```
pd = -(sigma2/dx-nu)/(2*dx)
pu = -(sigma2/dx+nu)/(2*dx)
pm = sigma2/(dx*dx)+r
```

The fitted version replaces `sigma2` with `2*ρ` in the second-derivative
terms, where ρ is the fitting factor.

### Modified: `fdtypedefs.hpp`

Add typedefs for the new scheme instantiations.

### Modified: `finitedifferencemodel.hpp`

No changes should be needed if the new schemes follow the same `step()`
/ `setStep()` interface, but verify.

## Files to produce (complete contents)

1. **`exponentiallyfittedscheme_old.hpp`** — New fitted scheme for old framework
2. **`milevtaglianischeme_old.hpp`** — New CN variant for old framework
3. **`bsmoperator.hpp`** and **`bsmoperator.cpp`** — Updated with fitted operator
4. **`fdtypedefs.hpp`** — Updated with new typedefs
5. **`all.hpp`** — Updated includes

## Existing code provided

[ATTACH: All root-level FD .hpp files: cranknicolson.hpp, mixedscheme.hpp,
         impliciteuler.hpp, expliciteuler.hpp, bsmoperator.hpp/cpp,
         tridiagonaloperator.hpp/cpp, finitedifferencemodel.hpp,
         boundarycondition.hpp/cpp, pde.hpp, pdebsm.hpp,
         fdtypedefs.hpp, operatortraits.hpp, stepcondition.hpp, all.hpp]
[ATTACH: Paper sections 2–3]

## Constraints

- The new old-framework classes must be usable as drop-in replacements
  for `CrankNicolson<TridiagonalOperator>` in `FiniteDifferenceModel<>`.
- The tridiagonal solve in the old framework uses
  `TridiagonalOperator::solveFor()` — the fitted operators must be
  compatible with this.
- Output COMPLETE file contents.
```

---

## ROUND 5 PROMPT — `pricingengines/`

```markdown
# Task — Round 5: Pricing Engines

In this round you update the user-facing pricing engines to expose the
new schemes and create specialized engines for the paper's test cases.

## What to create/modify

### Modified: `pricingengines/vanilla/fdblackscholesvanillaengine.hpp` / `.cpp`

Add the ability to select `FdmSchemeDesc::ExponentiallyFitted()` or
`FdmSchemeDesc::MilevTagliani()` as the scheme. The existing engine
already accepts an `FdmSchemeDesc`; the new scheme descriptors should
work transparently.

However, for `MilevTaglianiType`, the engine must pass r and σ to the
scheme/operator. Add these as optional parameters or extract them from
the process.

Add a `MakeFdBlackScholesVanillaEngine` builder method for the new schemes.

### Modified: `pricingengines/barrier/fdblackscholesbarrierengine.hpp` / `.cpp`

Update to support the new schemes for barrier option pricing. The paper
specifically tests with truncated call options and discrete double barrier
knock-out calls.

For discrete barriers, the engine must:
1. Apply the indicator function 1_{[L,U]}(S_j) at each monitoring date
2. This resets V to zero outside the corridor

### New: `pricingengines/vanilla/fdblackscholesvanillaengine_fitted.hpp` / `.cpp`
(Optional — only if a separate engine is cleaner)

A convenience engine that defaults to the exponentially fitted scheme.

## Files to produce (complete contents)

1. **`vanilla/fdblackscholesvanillaengine.hpp`** — Updated
2. **`vanilla/fdblackscholesvanillaengine.cpp`** — Updated
3. **`barrier/fdblackscholesbarrierengine.hpp`** — Updated
4. **`barrier/fdblackscholesbarrierengine.cpp`** — Updated

Also update the old-framework engine if applicable:
5. **`vanilla/fdvanillaengine.hpp`** — Updated for old framework
6. **`vanilla/fdvanillaengine.cpp`** — Updated

## Existing code provided

[ATTACH: pricingengines/vanilla/fdblackscholesvanillaengine.hpp/cpp,
         pricingengines/barrier/fdblackscholesbarrierengine.hpp/cpp,
         pricingengines/barrier/fdblackscholesrebateengine.hpp/cpp,
         pricingengines/vanilla/fdvanillaengine.hpp/cpp]
[ATTACH: Round 1–3 outputs (scheme, operator, solver files)]
[ATTACH: Paper section 4 (numerical examples) and definitions 4.1, 4.2]

## Constraints

- Backward compatibility: existing code using `FdmSchemeDesc::Douglas()`
  or `FdmSchemeDesc::CrankNicolson()` must continue to work unchanged.
- The engine must correctly handle the time-step restriction for the
  Milev-Tagliani scheme: Δt < 1/(r(1/2-2ω₁) + ½(σM)²).
  If the user-requested time steps violate this, either warn or auto-adjust.
- Output COMPLETE file contents.
```

---

## ROUND 6 PROMPT — `utilities/` + `meshers/` + `stepconditions/`

```markdown
# Task — Round 6: Utilities, Meshers, and Step Conditions

Final round. Update supporting infrastructure for the new schemes.

## What to modify

### `stepconditions/fdmstepconditioncomposite.cpp`

Ensure `vanillaComposite()` works correctly with discrete monitoring.
For the paper's discrete double barrier knock-out call (Definition 4.2):

V(S, t_i) = V(S, t_i^-) · 1_{[L,U]}(S)

at monitoring dates t_1, ..., t_F.

Create or verify a step condition that applies this indicator multiplication.
The existing `FdmBermudanStepCondition` is similar but for exercise.
A new `FdmDiscreteBarrierStepCondition` may be needed.

### `utilities/fdmdirichletboundary.cpp`

Verify that Dirichlet boundary conditions work correctly with the new
operators. The fitted operator and the CN variant both modify the
tridiagonal structure; the boundary condition application
(`applyBeforeApplying`, `applyAfterSolving`) must remain compatible.

### `meshers/fdmblackscholesmesher.cpp`

For low volatility problems (σ² ≪ r), the standard mesher may produce
inadequate grids. Consider adding a mode that:

1. Places more grid points near discontinuities (barrier levels, strike)
2. Uses a finer grid when σ is very small (the paper uses ΔS = 0.05 and
   ΔS = 0.025 in examples, with S_max = 140)

### `utilities/fdminnervaluecalculator.cpp`

For the truncated call payoff (Definition 4.1):

f[S(T)] = { S(T) - K  if S(T) ∈ [K,U],  0 otherwise }

Verify that `FdmLogInnerValue` correctly handles this discontinuous payoff.
The cell-averaging version `FdmCellAveragingInnerValue` may help smooth
the initial condition.

## Files to produce (complete contents)

1. **`stepconditions/fdmdiscretebarriercondition.hpp`** — New (if needed)
2. **`stepconditions/fdmdiscretebarriercondition.cpp`** — New (if needed)
3. **`stepconditions/fdmstepconditioncomposite.hpp`** — Updated
4. **`stepconditions/fdmstepconditioncomposite.cpp`** — Updated
5. **`stepconditions/all.hpp`** — Updated
6. **`utilities/fdmdirichletboundary.cpp`** — Updated if needed
7. **`meshers/fdmblackscholesmesher.cpp`** — Updated if needed

## Existing code provided

[ATTACH: All stepconditions/*.hpp/cpp, utilities/*.hpp/cpp, meshers/*.hpp/cpp]
[ATTACH: Round 1–5 outputs]
[ATTACH: Paper section 4, definitions 4.1–4.2, conditions (11)–(13)]

## Constraints

- The discrete barrier condition must correctly handle the indicator
  function at monitoring dates, setting values to zero outside [L,U].
- For the paper's examples: K=100, σ=0.001, T=1, r=0.05, L=95, U=110
  with daily monitoring (250 dates/year) or weekly (50 dates/year).
- The mesher should support the paper's test parameters:
  S_max=140, ΔS=0.05 down to ΔS=0.01.
- Output COMPLETE file contents.
```

---

## Verification Round Prompt (Optional Round 7)

```markdown
# Task — Round 7: Verification and Integration Test Scaffolding

Create a self-contained test file that reproduces the paper's numerical
examples. This is NOT a QuantLib test suite file, but a standalone
verification program.

## Test cases from the paper

### Example 4.1: Truncated Call Option
- r=0.05, σ=0.001, T=5/12, U=70, K=50, S_max=140
- Compare: Standard CN (should oscillate), Fully implicit (should oscillate),
  Duffy fitted (no oscillations, some diffusion), CN variant (no oscillations,
  some diffusion)
- Vary ΔS: 0.05, 0.025, 0.01
- Vary r: 0.01 (small), 0.5 (large)

### Example 4.2: Discrete Double Barrier Knock-Out Call
- K=100, σ=0.001, T=1, r=0.05, L=95 (or L=90 per Fig.2), U=110
- ΔS=0.025, Δt=0.001

### Verification criteria
- Standard CN and fully implicit: negative prices near barriers (spurious oscillations)
- Duffy fitted: positive prices, smooth, increasing diffusion with r
- CN variant: positive prices, smooth, controllable diffusion via ΔS
- All schemes: discrete maximum principle ||U^{n+1}||_∞ ≤ ||U^n||_∞

## Produce

A single **`test_milev_tagliani.cpp`** file that:
1. Sets up each test case
2. Runs all four schemes
3. Prints results in a table format
4. Flags oscillations (any negative price) and positivity violations

[ATTACH: All Round 1–6 outputs]
[ATTACH: Complete paper]
```

---

## Practical Notes for Each Round

When submitting each round to Claude Opus 4.6, attach:

| Round | Paper | Summary | Source Code to Attach | Previous Round Output |
|-------|-------|---------|----------------------|----------------------|
| 1 | Full PDF | Sections 2–3 | `schemes/*.hpp/cpp`, `operators/fdmlinearopcomposite.hpp`, `operators/triplebandlinearop.hpp` | — |
| 2 | Sections 2–3 | Sections 3.2–3.3 | All `operators/*.hpp/cpp` | Round 1 files |
| 3 | Sections 2–4 | Sections 3.2–4.2 | All `solvers/*.hpp/cpp` | Rounds 1–2 files |
| 4 | Full PDF | Sections 2–3 | All root-level FD `*.hpp/cpp` | Summary of Rounds 1–3 interfaces |
| 5 | Section 4 | Sections 4.1, 5.1–5.2 | Relevant `pricingengines/*.hpp/cpp` | Rounds 1–3 files |
| 6 | Section 4 | Definitions 4.1–4.2 | `utilities/*`, `meshers/*`, `stepconditions/*` | All prior rounds |

Each round's output becomes input to subsequent rounds, ensuring consistency across the full implementation.

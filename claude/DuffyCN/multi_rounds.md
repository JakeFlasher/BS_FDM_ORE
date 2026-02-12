*Thinking...*

> The user wants me to draft a multi-round prompting strategy for Claude Opus 4.6 to implement an improved Crank-Nicolson finite difference method in QuantLib v1.23 based on the Duffy paper's recommendations. Let me analyze what needs to be done:
> 
> 1. The paper recommends:
>    - Exponentially fitted spatial discretization (replacing centered differences)
>    - The fitting factor ρ = (μh/2σ)·coth(μh/2σ)
>    - Implicit Euler time stepping with the fitted scheme (not CN time stepping)
>    - Keller box scheme as an alternative
>    - Rannacher damping (implicit Euler startup steps)
>    
> 2. I need to structure this in rounds to avoid context overflow, starting with the mathematical core and working outward: first the exponentially fitted operators and scheme classes, then the solver infrastructure that uses them, followed by the pricing engines, and finally supporting utilities like meshers and boundary conditions.
> 
> For the solver layer, I'll need to modify the backward solver to dispatch to the new scheme, update the Black-Scholes solver and the 1D solver, and extend the scheme descriptor. Then I'll tackle the legacy framework by updating the mixed scheme and Crank-Nicolson implementations along with their operators, before moving to the pricing engines and utilities.
> 
> But I should reconsider the approach—the user wants the implementation to follow the mathematical formulations from the paper exactly, so I need to make sure each round respects that constraint.
> 
> Looking at what the paper actually proposes: exponentially fitted spatial discretization with implicit Euler time stepping (not Crank-Nicolson), or alternatively the Keller box scheme. The key is the fitting factor that accounts for the drift term, which modifies how the discrete spatial operator is constructed. So rather than just improving the existing CN scheme, I'm replacing the spatial discretization with the exponentially fitted version and potentially switching the time stepping method entirely.
> 
> Now I'm setting up the discrete formulation with the fitted operator applied at the new time level, using implicit Euler for time integration. I need to structure the implementation plan carefully—breaking it into rounds where each produces complete, working files rather than incremental changes. The first round should establish the mathematical core with the new fitted operators, then subsequent rounds can build out the solver infrastructure and test cases. I need to continue mapping out the remaining implementation rounds—solver infrastructure, pricing engines, and the root-level framework files that tie everything together.
> 
> Actually, I'm reconsidering the approach here. Rather than implementing directly, I should be drafting prompts strategically. The real challenge is ordering the rounds to minimize cross-dependencies, so each round produces self-contained code that the next one builds on.
> 
> The optimal structure starts with the core mathematical operators since they have no dependencies on other modified files, then moves to the scheme layer which builds on those operators, followed by the solver layer that depends on both.
> 
> Now I need to draft the actual prompts themselves—each one should be self-contained with the specific files to modify, the mathematical formulations to implement, context from previous rounds, and clear output format instructions.

# Multi-Round Prompting Plan: Implementing Duffy's Improved Crank-Nicolson in QuantLib v1.23

## Implementation Priority & Dependency Chain

```
Round 1: operators/          ← mathematical core (no upstream dependencies)
    ↓
Round 2: schemes/ + root FD  ← time-stepping logic (depends on Round 1)
    ↓
Round 3: solvers/            ← orchestration (depends on Rounds 1-2)
    ↓
Round 4: pricingengines/     ← user-facing API (depends on Rounds 1-3)
    ↓
Round 5: utilities/ + meshers/ ← supporting infrastructure (depends on Rounds 1-3)
```

---

## ROUND 1 PROMPT — Spatial Operators (Mathematical Core)

Attach: `methods.xml` (operators section), the paper PDF, `critique_CN.md` (Sections 3.4–3.5 only), and `hpp_structure_of_quantlib.md` (Section 2.3 only).

---

```
You are implementing Daniel J. Duffy's exponentially fitted finite-difference
spatial discretization from "A Critique of the Crank Nicolson Scheme" (Wilmott
2004) into QuantLib v1.23's new finite-differences framework.

This is ROUND 1 of 5. You are modifying ONLY the spatial operator layer
(methods/finitedifferences/operators/). Later rounds will modify schemes,
solvers, and engines to use your operators.

═══════════════════════════════════════════════════════════════════════════
MATHEMATICAL SPECIFICATION (implement exactly)
═══════════════════════════════════════════════════════════════════════════

The Black-Scholes PDE in log-space x = ln(S) is a convection-diffusion
equation:

  Lu = -∂u/∂t + σ(x,t) ∂²u/∂x² + μ(x,t) ∂u/∂x + b(x,t)u = f(x,t)

where for Black-Scholes:
  σ(x,t) = ½ vol²          (diffusion coefficient)
  μ(x,t) = r - q - ½ vol²  (drift/convection coefficient)
  b(x,t) = -r              (reaction/discount coefficient)

STANDARD centered differences (EXISTING in QuantLib) use:
  ∂²u/∂x² ≈ D₊D₋ Uⱼ = (Uⱼ₊₁ - 2Uⱼ + Uⱼ₋₁)/h²
  ∂u/∂x   ≈ D₀ Uⱼ   = (Uⱼ₊₁ - Uⱼ₋₁)/(2h)

FITTED discretization (TO IMPLEMENT) replaces the second derivative with a
fitted version:
  ρⱼ D₊D₋ Uⱼ + μⱼ D₀ Uⱼ + bⱼ Uⱼ

where the FITTING FACTOR is:
  ρⱼ = (μⱼ h / 2) · coth(μⱼ h / (2 σⱼ))

and coth(x) = (e^x + e^(-x)) / (e^x - e^(-x)) = (e^(2x) + 1) / (e^(2x) - 1)

CRITICAL LIMITING BEHAVIORS (must be handled numerically):
  • When σⱼ → 0 (pure convection):
      ρⱼ → +μⱼh/2  if μⱼ > 0    (forward upwind)
      ρⱼ → -μⱼh/2  if μⱼ < 0    (backward upwind)
  • When μⱼ → 0 (pure diffusion):
      ρⱼ → σⱼ                    (standard centered scheme recovered)
  • Use: lim_{x→0} x·coth(x) = 1

MONOTONICITY REQUIREMENT — the resulting tridiagonal coefficients must satisfy:
  aⱼ,ⱼ₋₁ = ρⱼ/h² - μⱼ/(2h) > 0   (ALWAYS, by construction of ρ)
  aⱼ,ⱼ   = -2ρⱼ/h² + bⱼ < 0      (ALWAYS, since ρ > 0 and b ≤ 0)
  aⱼ,ⱼ₊₁ = ρⱼ/h² + μⱼ/(2h) > 0   (ALWAYS, by construction of ρ)

═══════════════════════════════════════════════════════════════════════════
FILES TO PRODUCE (output complete files, not diffs)
═══════════════════════════════════════════════════════════════════════════

1. NEW FILE: operators/fdmfittedblackscholesop.hpp
   • Create class FdmFittedBlackScholesOp : public FdmLinearOpComposite
   • Same interface as FdmBlackScholesOp (size, setTime, apply, apply_mixed,
     apply_direction, solve_splitting, preconditioner, toMatrixDecomp)
   • Constructor takes same arguments as FdmBlackScholesOp PLUS an optional
     bool useFitting = true
   • Internally stores dxMap_ (FirstDerivativeOp), dxxMap_ (SecondDerivativeOp)
   • In setTime(): computes fitting factor ρⱼ for each grid point, then
     assembles the fitted operator using axpyb on TripleBandLinearOp

2. NEW FILE: operators/fdmfittedblackscholesop.cpp
   • Implement the fitting factor computation with numerical safeguards:
     - If |μⱼh/(2σⱼ)| < 1e-4, use Taylor expansion: coth(x) ≈ 1/x + x/3 - x³/45
     - If σⱼ < 1e-12, use the upwind limit directly
     - If |μⱼ| < 1e-12, use ρⱼ = σⱼ (pure diffusion limit)
   • The setTime() method must:
     a. Get r, q from term structures (same as FdmBlackScholesOp)
     b. Get vol² (local vol or constant vol, same as FdmBlackScholesOp)
     c. Compute μⱼ = r - q - ½vol² and σⱼ = ½vol² at each grid point
     d. Compute ρⱼ for each grid point
     e. Assemble: mapT_.axpyb(μ_array, dxMap_, dxxMap_.mult(ρ_array), Array(1, b))
        where b = -r

3. MODIFIED FILE: operators/fdmblackscholesop.hpp
   • Add: #include "fdmfittedblackscholesop.hpp" at top
   • No other changes to this file (preserve backward compatibility)

4. MODIFIED FILE: operators/fdmblackscholesop.cpp
   • No changes (preserve existing behavior for unfitted case)

5. UPDATE: operators/all.hpp
   • Add: #include <ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp>

═══════════════════════════════════════════════════════════════════════════
DESIGN CONSTRAINTS
═══════════════════════════════════════════════════════════════════════════

• PRESERVE all existing QuantLib v1.23 interfaces and behavior
• The fitted operator must work with ALL existing time-stepping schemes
  (implicit Euler, Crank-Nicolson, Douglas, etc.) since it replaces only
  the spatial discretization
• Follow QuantLib coding conventions: ext::shared_ptr, Disposable<Array>,
  QL_REQUIRE for preconditions, namespace QuantLib
• Use Real (not double), Size (not size_t), Time (not double) per QuantLib types
• Handle both constant-vol and local-vol cases (same logic as FdmBlackScholesOp)
• The fitting factor computation MUST handle non-uniform grids (use
  mesher->dminus and mesher->dplus, not a constant h)

═══════════════════════════════════════════════════════════════════════════
ATTACHED CODE CONTEXT
═══════════════════════════════════════════════════════════════════════════

[Attach: the operators section from methods.xml, specifically:
 - fdmblackscholesop.hpp + .cpp (as template)
 - triplebandlinearop.hpp + .cpp (for axpyb interface)
 - firstderivativeop.hpp + .cpp (for dxMap_)
 - secondderivativeop.hpp + .cpp (for dxxMap_)
 - fdmlinearopcomposite.hpp (interface to implement)
 - fdmlinearoplayout.hpp (for iteration)]

Output every file in full. Do not abbreviate or use "... rest unchanged ...".
```

---

## ROUND 2 PROMPT — Time-Stepping Schemes

Attach: `schemes/` from `methods.xml`, root-level FD files from `finitedifferences.xml`, Round 1 output files, `critique_CN.md` (Sections 3.2–3.3, 3.5), and `hpp_structure_of_quantlib.md` (Sections 2.1–2.2).

---

```
You are implementing Daniel J. Duffy's recommendations for improved time
stepping from "A Critique of the Crank Nicolson Scheme" (Wilmott 2004) into
QuantLib v1.23.

This is ROUND 2 of 5. In Round 1, we created FdmFittedBlackScholesOp (a new
spatial operator with exponential fitting). Now you are modifying the
TIME-STEPPING layer to:

(A) Add a new scheme type "FittedImplicitEuler" that pairs naturally with the
    fitted spatial operator (the paper's primary recommendation)
(B) Add Rannacher-style damping as a first-class option (a few fully-implicit
    steps followed by CN — already partially supported via dampingSteps, but
    we formalize it as a scheme variant)
(C) Modify the legacy CN framework to optionally use fitted spatial operators

═══════════════════════════════════════════════════════════════════════════
MATHEMATICAL SPECIFICATION
═══════════════════════════════════════════════════════════════════════════

PAPER'S FULLY DISCRETE FITTED SCHEME (equation 16):
  -(Uⱼⁿ⁺¹ - Uⱼⁿ)/k + ρⱼⁿ⁺¹ D₊D₋ Uⱼⁿ⁺¹ + μⱼⁿ⁺¹ D₀ Uⱼⁿ⁺¹ + bⱼⁿ⁺¹ Uⱼⁿ⁺¹ = fⱼⁿ⁺¹

This is IMPLICIT EULER in time with the fitted spatial operator. Written as:
  Uⱼⁿ⁺¹ - k·L_fitted·Uⱼⁿ⁺¹ = Uⱼⁿ + k·fⱼⁿ⁺¹

Convergence: |u(xⱼ,tₙ) - Uⱼⁿ| ≤ M(h+k) with M INDEPENDENT of h, k, and σ.
(First order in time, but with Richardson extrapolation → second order)

RANNACHER STARTUP (practical recommendation from literature):
  • Use d fully-implicit steps (d typically 2-4) at the start
  • Then switch to Crank-Nicolson for remaining steps
  • This damps oscillations from non-smooth initial data (payoff kinks)
  • QuantLib already has dampingSteps parameter — we ensure it works with
    the fitted operator

RICHARDSON EXTRAPOLATION for the fitted implicit scheme:
  • Run with step k → get U(k)
  • Run with step k/2 → get U(k/2)  
  • Extrapolate: U* = 2·U(k/2) - U(k) → second order in time
  • This recovers the accuracy advantage of CN without CN's oscillations

═══════════════════════════════════════════════════════════════════════════
FILES TO PRODUCE
═══════════════════════════════════════════════════════════════════════════

1. MODIFIED: solvers/fdmbackwardsolver.hpp
   • Add to FdmSchemeDesc::FdmSchemeType enum:
     FittedImplicitEulerType, RannacherCNType
   • Add static factory methods:
     FdmSchemeDesc::FittedImplicitEuler()
     FdmSchemeDesc::RannacherCrankNicolson(Size rannacherSteps = 2)
   • The mu field can encode rannacherSteps for the Rannacher variant

2. MODIFIED: solvers/fdmbackwardsolver.cpp
   • Add cases in FdmBackwardSolver::rollback() for:
     - FittedImplicitEulerType: use ImplicitEulerScheme (the existing one
       works since the fitted spatial operator handles the fitting)
     - RannacherCNType: explicit Rannacher startup logic:
       * First mu steps: ImplicitEulerScheme
       * Remaining steps: CrankNicolsonScheme
       (Note: this is different from the existing dampingSteps which always
        uses ImplicitEuler for damping regardless of main scheme. Here we
        make it a first-class scheme choice.)

3. MODIFIED: schemes/cranknicolsonscheme.hpp
   • No interface changes, but add a comment documenting that when paired
     with FdmFittedBlackScholesOp, the scheme gains monotonicity in space

4. MODIFIED: schemes/cranknicolsonscheme.cpp
   • No functional changes (the CN scheme is agnostic to which spatial
     operator it uses — the fitting happens in the operator)

5. MODIFIED (legacy): cranknicolson.hpp
   • Add a comment block explaining the fitted alternative
   • No functional changes to preserve backward compatibility

6. MODIFIED (legacy): mixedscheme.hpp
   • No functional changes

7. NEW FILE: schemes/richardsonextrapolationscheme.hpp
   • Template class RichardsonExtrapolationScheme<BaseScheme>
   • step() method: runs two sub-steps at half dt, extrapolates
   • This provides second-order time accuracy when BaseScheme = ImplicitEuler

8. NEW FILE: schemes/richardsonextrapolationscheme.cpp
   • (or header-only if template)

═══════════════════════════════════════════════════════════════════════════
ROUND 1 OUTPUT (for context — the fitted operator you depend on)
═══════════════════════════════════════════════════════════════════════════

[Paste: FdmFittedBlackScholesOp header and key method signatures from Round 1]

═══════════════════════════════════════════════════════════════════════════
ATTACHED CODE CONTEXT
═══════════════════════════════════════════════════════════════════════════

[Attach: schemes/ section from methods.xml:
 - cranknicolsonscheme.hpp + .cpp
 - impliciteulerscheme.hpp + .cpp
 - expliciteulerscheme.hpp + .cpp
 - boundaryconditionschemehelper.hpp
 And root-level:
 - cranknicolson.hpp
 - mixedscheme.hpp  
 - finitedifferencemodel.hpp
 And solvers:
 - fdmbackwardsolver.hpp + .cpp]

Output every file in full.
```

---

## ROUND 3 PROMPT — Solver Layer

Attach: `solvers/` from `methods.xml`, Round 1+2 output headers, `hpp_structure_of_quantlib.md` (Section 2.4).

---

```
You are implementing the solver layer changes for Duffy's improved
Crank-Nicolson scheme in QuantLib v1.23.

This is ROUND 3 of 5. Rounds 1-2 created:
  • FdmFittedBlackScholesOp (fitted spatial operator)
  • New FdmSchemeDesc types (FittedImplicitEulerType, RannacherCNType)
  • Modified FdmBackwardSolver::rollback() with new scheme dispatch
  • RichardsonExtrapolationScheme template

Now you modify the SOLVER layer so that:
  (A) FdmBlackScholesSolver can optionally use the fitted operator
  (B) Fdm1DimSolver properly supports the new scheme types
  (C) FdmSolverDesc carries any new configuration needed

═══════════════════════════════════════════════════════════════════════════
FILES TO PRODUCE
═══════════════════════════════════════════════════════════════════════════

1. MODIFIED: solvers/fdmsolverdesc.hpp
   • Add optional field: bool useExponentialFitting = false
   • Add optional field: bool useRichardsonExtrapolation = false

2. MODIFIED: solvers/fdmblackscholessolver.hpp
   • Add constructor parameter: bool useExponentialFitting = false
   • Store as member

3. MODIFIED: solvers/fdmblackscholessolver.cpp
   • In performCalculations():
     - If useExponentialFitting_: create FdmFittedBlackScholesOp
     - Else: create FdmBlackScholesOp (existing behavior)
     - Pass to Fdm1DimSolver as before

4. MODIFIED: solvers/fdm1dimsolver.hpp
   • No interface changes needed (it's agnostic to operator type)

5. MODIFIED: solvers/fdm1dimsolver.cpp
   • No changes needed (rollback is driven by scheme selection)

6. MODIFIED: solvers/fdmbackwardsolver.hpp (if not fully done in Round 2)
   • Ensure FdmSchemeDesc static factories are complete

7. MODIFIED: solvers/fdmbackwardsolver.cpp (if not fully done in Round 2)
   • Ensure rollback() handles Richardson extrapolation:
     If schemeDesc specifies Richardson:
       - Run rollback with steps N → get result1
       - Run rollback with steps 2N → get result2  
       - Return 2*result2 - result1

═══════════════════════════════════════════════════════════════════════════
ROUND 1-2 OUTPUT (headers only, for context)
═══════════════════════════════════════════════════════════════════════════

[Paste: FdmFittedBlackScholesOp.hpp, modified fdmbackwardsolver.hpp from
 previous rounds]

═══════════════════════════════════════════════════════════════════════════
ATTACHED CODE CONTEXT
═══════════════════════════════════════════════════════════════════════════

[Attach: solvers/ section from methods.xml:
 - fdmblackscholessolver.hpp + .cpp
 - fdm1dimsolver.hpp + .cpp
 - fdmsolverdesc.hpp
 - fdmbackwardsolver.hpp + .cpp (original, for reference)]

Output every file in full.
```

---

## ROUND 4 PROMPT — Pricing Engines

Attach: relevant engines from `pricingengines.xml`, Round 1-3 output headers, `hpp_structure_of_quantlib.md` (Section 2.8).

---

```
You are modifying QuantLib v1.23 pricing engines to expose the improved
Crank-Nicolson (Duffy exponentially fitted) scheme to end users.

This is ROUND 4 of 5. Rounds 1-3 created:
  • FdmFittedBlackScholesOp (fitted spatial operator)
  • New FdmSchemeDesc types + Richardson extrapolation support
  • Modified FdmBlackScholesSolver with useExponentialFitting option

Now you modify PRICING ENGINES to expose new options.

═══════════════════════════════════════════════════════════════════════════
FILES TO PRODUCE
═══════════════════════════════════════════════════════════════════════════

1. MODIFIED: pricingengines/vanilla/fdblackscholesvanillaengine.hpp
   • Add to constructor and Make builder:
     bool useExponentialFitting = false
   • Store as member

2. MODIFIED: pricingengines/vanilla/fdblackscholesvanillaengine.cpp
   • In calculate():
     - Pass useExponentialFitting_ to FdmBlackScholesSolver constructor
     - If using FittedImplicitEuler scheme, suggest dampingSteps = 0
       (the fitted scheme handles oscillation control inherently)

3. MODIFIED: pricingengines/barrier/fdblackscholesbarrierengine.hpp
   • Add: bool useExponentialFitting = false

4. MODIFIED: pricingengines/barrier/fdblackscholesbarrierengine.cpp
   • Pass useExponentialFitting_ to FdmBlackScholesSolver

5. MODIFIED: pricingengines/barrier/fdblackscholesrebateengine.hpp + .cpp
   • Same pattern as barrier engine

═══════════════════════════════════════════════════════════════════════════
USAGE EXAMPLE (for validation — include as comments in the engine)
═══════════════════════════════════════════════════════════════════════════

// Standard CN (existing behavior, unchanged):
auto engine1 = make_shared<FdBlackScholesVanillaEngine>(
    process, 100, 100, 0, FdmSchemeDesc::CrankNicolson());

// Duffy's fitted implicit scheme (new):
auto engine2 = make_shared<FdBlackScholesVanillaEngine>(
    process, 100, 100, 0, FdmSchemeDesc::FittedImplicitEuler(),
    false, -Null<Real>(),
    FdBlackScholesVanillaEngine::Spot,
    true  /* useExponentialFitting */);

// Rannacher CN with fitting (new):
auto engine3 = make_shared<FdBlackScholesVanillaEngine>(
    process, 100, 100, 2 /* Rannacher steps */,
    FdmSchemeDesc::CrankNicolson(),
    false, -Null<Real>(),
    FdBlackScholesVanillaEngine::Spot,
    true  /* useExponentialFitting */);

═══════════════════════════════════════════════════════════════════════════
ROUND 1-3 OUTPUT (headers only)
═══════════════════════════════════════════════════════════════════════════

[Paste: key headers from Rounds 1-3]

═══════════════════════════════════════════════════════════════════════════
ATTACHED CODE CONTEXT
═══════════════════════════════════════════════════════════════════════════

[Attach from pricingengines.xml:
 - vanilla/fdblackscholesvanillaengine.hpp + .cpp
 - barrier/fdblackscholesbarrierengine.hpp + .cpp
 - barrier/fdblackscholesrebateengine.hpp + .cpp]

Output every file in full.
```

---

## ROUND 5 PROMPT — Utilities, Meshers & Legacy BSM Operator

Attach: `utilities/` and `meshers/` from `methods.xml`, legacy FD files from `finitedifferences.xml`, Round 1-4 output headers.

---

```
You are completing the implementation of Duffy's improved Crank-Nicolson
in QuantLib v1.23 by modifying supporting utilities and the legacy framework.

This is ROUND 5 of 5 (final round). Rounds 1-4 created the fitted spatial
operator, new scheme types, solver modifications, and engine API changes.

Now you handle:
  (A) Legacy BSM operator with optional exponential fitting
  (B) Mesh quality improvements for the fitted scheme
  (C) Inner value calculator adjustments for payoff smoothing
  (D) Update all.hpp includes

═══════════════════════════════════════════════════════════════════════════
MATHEMATICAL CONTEXT FOR LEGACY OPERATOR
═══════════════════════════════════════════════════════════════════════════

The LEGACY BSMOperator uses on a UNIFORM log-grid with spacing dx:
  pd = -(σ²/dx - ν)/(2·dx)     where ν = r - q - σ²/2
  pu = -(σ²/dx + ν)/(2·dx)
  pm = σ²/(dx·dx) + r

With EXPONENTIAL FITTING, replace σ² with the fitting factor:
  ρ = (ν·dx/2) · coth(ν·dx/(2·σ²/2))   [note σ² appears as 2·(½σ²)]
    = (ν·dx) · coth(ν·dx/σ²) / 2

Then:
  pd_fitted = -(ρ/dx - ν/2)/dx = -ρ/(dx²) + ν/(2·dx)  → ALWAYS ≤ 0
  pu_fitted = -(ρ/dx + ν/2)/dx = -ρ/(dx²) - ν/(2·dx)  → ALWAYS ≤ 0  
  pm_fitted = 2ρ/(dx²) + r

═══════════════════════════════════════════════════════════════════════════
FILES TO PRODUCE
═══════════════════════════════════════════════════════════════════════════

1. NEW FILE: bsmfittedoperator.hpp
   • class BSMFittedOperator : public TridiagonalOperator
   • Same constructors as BSMOperator but computes ρ-based stencils
   • Include numerical safeguards for the coth computation

2. NEW FILE: bsmfittedoperator.cpp
   • Implement both uniform-dx and array-grid constructors
   • For the array-grid version, compute ρ at each grid point using
     local values of ν and σ² from the LogGrid

3. MODIFIED: bsmoperator.hpp
   • Add #include for new fitted operator (for discoverability)

4. MODIFIED: pricingengines/vanilla/fdvanillaengine.hpp + .cpp
   • Add option to use BSMFittedOperator instead of BSMOperator
   • In initializeOperator(): branch on useFitting_ flag

5. MODIFIED: meshers/fdmblackscholesmesher.cpp
   • When exponential fitting is in use, ensure the grid concentrates
     points near the strike (where payoff kink causes CN oscillations)
   • The existing concentrating mesher already does this — just ensure
     the default cPoint is set to (strike, 0.1) when fitting is on

6. MODIFIED: utilities/fdminnervaluecalculator.cpp
   • In FdmCellAveragingInnerValue::avgInnerValueCalc():
     This already implements cell-averaging via Simpson integration
     to smooth the payoff — verify it works correctly with the
     fitted scheme (no changes expected, but confirm)

7. UPDATE: all.hpp files
   • methods/finitedifferences/operators/all.hpp: add fdmfittedblackscholesop
   • methods/finitedifferences/all.hpp: add bsmfittedoperator
   • methods/finitedifferences/schemes/all.hpp: add richardsonextrapolationscheme
   
═══════════════════════════════════════════════════════════════════════════
ATTACHED CODE CONTEXT
═══════════════════════════════════════════════════════════════════════════

[Attach:
 - bsmoperator.hpp + .cpp
 - tridiagonaloperator.hpp + .cpp
 - fdvanillaengine.hpp + .cpp (from pricingengines.xml)
 - meshers/fdmblackscholesmesher.hpp + .cpp
 - utilities/fdminnervaluecalculator.hpp + .cpp
 - all.hpp files]

Output every file in full.
```

---

## Operational Notes for Running This Plan

**Passing context between rounds:** After each round, extract the *headers only* (`.hpp` files) from the output and paste them into the next round's prompt as "ROUND N OUTPUT." This provides the interface contracts without overwhelming the context window with implementation details.

**Verification checkpoints:** After each round, verify that the output compiles conceptually against the QuantLib interfaces by checking that all `#include` paths are valid, all base class pure virtuals are implemented, and all `ext::shared_ptr` / `Disposable<Array>` types are used correctly.

**Context budget per round:** Each round should stay under ~80K tokens of input. The largest round is Round 1 (operators section ≈ 30K tokens of context + 5K of prompt). If a round exceeds limits, split the attached code: provide headers for reference files and full source only for files being modified.

**If output is truncated:** Reissue the prompt for just the truncated file, stating "Continue outputting [filename] from where you left off. Begin at line [N]."

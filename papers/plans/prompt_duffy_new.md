*Thinking...*

> The user wants me to draft a multi-round prompt strategy for implementing the Crank-Nicolson improvements described in the Duffy paper into QuantLib v1.23. Let me analyze the paper and the code structure to determine:
> 
> 1. What changes are needed based on the paper
> 2. The priority order of changes
> 3. How to split the work across multiple prompts to avoid context overflow
> 
> From the paper, the key improvements are:
> 1. Exponentially fitted finite difference schemes (replacing centered differences with fitted coefficients using coth formula)
> 2. The Keller box scheme as an alternative
> 3. Better handling of non-smooth initial conditions
> 4. Better stability analysis (discrete maximum principle instead of von Neumann)
> 
> 
> 5. Proper handling of boundary conditions
> 6. Graceful degradation for small volatility (upwind schemes in the limit)
> 
> From the code structure analysis, the modification order should be:
> 
> **Round 1: Core operator modifications (exponential fitting)**
> - `methods/finitedifferences/operators/` - Modify spatial discretization operators to use exponential fitting
> - This is the foundation - all other changes depend on correct spatial operators
> 
> **Round 2: Scheme modifications**
> - `methods/finitedifferences/schemes/` - Modify or create new CN scheme with Rannacher damping built-in
> - Root-level `methods/finitedifferences/` files - Modify old framework
> 
> **Round 3: Solver and boundary condition modifications**
> - `methods/finitedifferences/solvers/` - Update solver orchestration
> - `methods/finitedifferences/utilities/` - Update boundary conditions
> 
> **Round 4: Mesher and engine modifications**
> - `methods/finitedifferences/meshers/` - Grid improvements
> - `pricingengines/vanilla/` - Engine configuration updates
> 
> **Round 5: Integration and testing**
> - Wire everything together
> - Test cases
> 
> The key mathematical implementations from the paper are the exponentially fitted scheme using ρ_j^h = (μ_j h / 2) coth(μ_j h / (2σ_j)) as a fitting factor, the fully discrete implicit scheme with fitted spatial operators, limiting behavior for small volatility, the Keller box scheme as an alternative, and discrete maximum principle properties ensuring monotonicity. I need to design 5-6 focused rounds where each prompt specifies relevant paper sections, code files to include, expected outputs, and how to pass results forward. Each round should be self-contained enough that a fresh Claude instance can execute it while maintaining the implementation strategy of creating new exponentially fitted operators, modifying schemes with Rannacher damping, and optionally adding the Keller box scheme as an alternative.

# Multi-Round Prompt Strategy: Implementing Duffy's Exponentially Fitted FD Scheme in QuantLib v1.23

Below is a complete prompt plan organized in **6 rounds**, ordered by dependency. Each round produces modified files that become inputs to subsequent rounds.

---

## Execution Order and Rationale

The paper's core contribution is replacing centered-difference spatial discretization with exponentially fitted coefficients, then using an implicit time-stepping scheme (not CN) that preserves monotonicity. The dependency chain is:

```
Round 1: New fitted spatial operators     ← foundation, no dependencies
Round 2: New scheme + modified CN scheme  ← depends on Round 1 operators
Round 3: Solver orchestration updates     ← depends on Round 2 schemes
Round 4: Old framework parallel changes   ← independent of Rounds 1-3, but logically after
Round 5: Engine + mesher wiring           ← depends on Rounds 1-3
Round 6: Verification & integration glue  ← depends on all above
```

---

## ROUND 1 — Fitted Spatial Operators

**Objective:** Create new exponentially fitted finite-difference operators that replace standard centered differences with the Duffy fitting factor ρ = (μh/2)·coth(μh/(2σ)).

**Files to include in prompt:**
- Paper PDF + critique_CN.md (full summary)
- hpp_structure_of_quantlib.md (architecture guide)
- From `methods.xml`: documents indexed 25 (firstderivativeop.cpp), 29 (secondderivativeop.cpp), 31 (triplebandlinearop.cpp), 14 (fdmblackscholesop.cpp), 21 (fdmlinearoplayout.cpp)
- From `methods.xml`: document 1 (boundarycondition.cpp) for BC interface reference

**Prompt text:**

```
You are implementing the exponentially fitted finite-difference scheme from Daniel Duffy's paper "A Critique of the Crank Nicolson Scheme" (Wilmott 2004) into QuantLib v1.23. This is Round 1 of 6: creating the fitted spatial operators.

## Paper Mathematics to Implement

The paper's key innovation is the exponentially fitted scheme. For the generic parabolic PDE:

Lu ≡ -∂u/∂t + σ(x,t)·∂²u/∂x² + μ(x,t)·∂u/∂x + b(x,t)·u = f(x,t)

The standard centered-difference scheme uses D+D- for ∂²/∂x² and D0 for ∂u/∂x. The fitted scheme replaces the second derivative coefficient with:

ρ_j^h = (μ_j·h / 2) · coth(μ_j·h / (2·σ_j))

where coth(x) = (e^x + e^{-x})/(e^x - e^{-x}).

The fitted discrete operator at each grid point j is:
  ρ_j^h · D+D- U_j + μ_j · D0 U_j + b_j · U_j = f_j

This produces tridiagonal coefficients:
  a_{j,j-1} = ρ_j^h/h² - μ_j/(2h)  > 0 always
  a_{j,j}   = -2ρ_j^h/h² + b_j      < 0 always  
  a_{j,j+1} = ρ_j^h/h² + μ_j/(2h)   > 0 always

These sign properties guarantee an M-matrix → monotonicity → no spurious oscillations.

**Limiting behavior (must be implemented):**
- When σ→0: lim ρ_j = |μ_j|·h/2, producing implicit upwind schemes
- When μ→0: lim(x·coth(x))=1, so ρ_j = σ_j, producing standard diffusion

For the Black-Scholes equation cast as:
  -∂V/∂t + σ²S²/2 · ∂²V/∂S² + rS · ∂V/∂S - rV = 0

After the log-transform x = ln(S), the convection coefficient is μ = r - q - σ²/2 and the diffusion coefficient is σ²/2.

## Your Task

Create the following NEW files and MODIFY existing files:

1. **NEW: `exponentiallyfittedop.hpp` / `exponentiallyfittedop.cpp`**
   A new TripleBandLinearOp subclass that computes the fitted spatial operator. It must:
   - Accept diffusion coefficient σ(x,t) and convection coefficient μ(x,t) arrays
   - Compute ρ_j = (μ_j·h/2)·coth(μ_j·h/(2·σ_j)) at each interior grid point
   - Handle the limiting cases: when |σ_j| < ε, use upwind; when |μ_j| < ε, use standard diffusion
   - Assemble the tridiagonal coefficients ensuring the M-matrix sign structure
   - Boundary points (j=0 and j=J) get zero entries (handled by BCs)

2. **NEW: `fdmfittedblackscholesop.hpp` / `fdmfittedblackscholesop.cpp`**
   A new FdmLinearOpComposite implementation for the Black-Scholes equation using exponential fitting. Similar interface to FdmBlackScholesOp but uses ExponentiallyFittedOp instead of separate FirstDerivativeOp + SecondDerivativeOp. The setTime() method must:
   - Extract r, q, σ from term structures
   - Compute μ = r - q - σ²/2 (convection) and d = σ²/2 (diffusion) at each grid point
   - Pass these to ExponentiallyFittedOp to build the fitted operator
   - Add the reaction term (-r) to the diagonal

3. **MODIFY: `triplebandlinearop.cpp`** — No changes needed to the Thomas algorithm itself, but verify the solve_splitting method works correctly with the new operator's coefficient structure.

Output every file in its entirety (complete compilable C++ files with all includes, namespaces, and method implementations). Use the QuantLib coding conventions visible in the attached source files.

[ATTACH: paper summary, architecture guide, and the 6 source files listed above]
```

---

## ROUND 2 — Scheme Modifications

**Objective:** Create a new `FittedImplicitScheme` and modify `CrankNicolsonScheme` to support Rannacher-style damping natively.

**Files to include:**
- Paper summary (abbreviated to Sections 4-6 + conclusions)
- Architecture guide (abbreviated to Sections 2.1-2.2)
- Round 1 output files (the new operator headers/sources)
- From `methods.xml`: documents 33 (cranknicolsonscheme.cpp), 35 (expliciteulerscheme.cpp), 37 (impliciteulerscheme.cpp), 32 (craigsneydscheme.cpp as template), 34 (douglasscheme.cpp as template)

**Prompt text:**

```
You are implementing Round 2 of 6: creating new time-stepping schemes based on Duffy's paper.

## Context from Round 1
In Round 1, we created:
- ExponentiallyFittedOp: a TripleBandLinearOp computing ρ_j = (μ_j·h/2)·coth(μ_j·h/(2·σ_j))
- FdmFittedBlackScholesOp: an FdmLinearOpComposite using the fitted operator

[ATTACH Round 1 output files here]

## Paper Mathematics for This Round

### Fitted Implicit Scheme (Paper equation 16):
The paper advocates a fully implicit scheme with the fitted operator:

  L_h^k U_j^n ≡ -(U_j^{n+1} - U_j^n)/k + ρ_j^{n+1}·D+D-·U_j^{n+1} + μ_j^{n+1}·D0·U_j^{n+1} + b_j^{n+1}·U_j^{n+1} = f_j^{n+1}

This is first-order in time but uniformly stable. Error bound (paper equation 18):
  |u(x_j,t_n) - U_j^n| ≤ M(h+k)  where M is independent of h, k, and σ.

### Richardson Extrapolation for Second-Order:
The paper recommends Richardson extrapolation to recover second-order accuracy:
- Solve with step k → solution U_k
- Solve with step k/2 → solution U_{k/2}  
- Extrapolated solution = 2·U_{k/2} - U_k

### Rannacher Startup for Standard CN:
As documented in the follow-up literature (Wade et al. 2007), use a few implicit Euler steps at startup to damp oscillations from non-smooth payoffs, then switch to CN.

## Your Task

1. **NEW: `fittedimplicitscheme.hpp` / `fittedimplicitscheme.cpp`**
   A new scheme class that:
   - Uses fully implicit Euler time-stepping with the fitted spatial operator
   - Has the same interface as CrankNicolsonScheme (step(), setStep())
   - Works with any FdmLinearOpComposite (but designed for FdmFittedBlackScholesOp)
   
2. **NEW: `fittedimplicitwithextrapolationscheme.hpp` / `fittedimplicitwithextrapolationscheme.cpp`**
   A wrapper scheme that:
   - Runs the fitted implicit scheme twice (full step and half step)
   - Applies Richardson extrapolation: result = 2·fine - coarse
   - Achieves second-order temporal accuracy while maintaining monotonicity
   
3. **MODIFY: `cranknicolsonscheme.hpp` / `cranknicolsonscheme.cpp`**
   Add a constructor parameter `Size rannacherSteps = 0` that:
   - If > 0, performs that many implicit Euler steps before switching to CN
   - The implicit steps use the same operator and boundary conditions
   - This addresses the paper's critique about CN oscillations with non-smooth payoffs

4. **MODIFY: `impliciteulerscheme.hpp` / `impliciteulerscheme.cpp`**
   Ensure the step(a, t, theta) method works correctly when theta=1.0 (full implicit) with the fitted operator's solve_splitting.

Output complete files. Follow QuantLib conventions from the attached source.

[ATTACH: paper summary sections, architecture sections, Round 1 outputs, and the 5 scheme source files]
```

---

## ROUND 3 — Solver Orchestration

**Objective:** Update `FdmBackwardSolver` to recognize new scheme types and update `FdmSchemeDesc` accordingly.

**Files to include:**
- Paper summary (abbreviated)
- Round 1 + Round 2 output file headers (just .hpp files to show interfaces)
- From `methods.xml`: document 44 (fdmbackwardsolver.cpp), 40 (fdm1dimsolver.cpp), 46 (fdmblackscholessolver.cpp)
- `finitedifferencemodel.hpp` content from the architecture guide

**Prompt text:**

```
You are implementing Round 3 of 6: updating the solver orchestration layer.

## Context from Previous Rounds
Round 1 created: ExponentiallyFittedOp, FdmFittedBlackScholesOp
Round 2 created: FittedImplicitScheme, FittedImplicitWithExtrapolationScheme, modified CrankNicolsonScheme with Rannacher steps

[ATTACH Round 1+2 header files showing interfaces]

## Your Task

1. **MODIFY: `fdmbackwardsolver.hpp` / `fdmbackwardsolver.cpp`**
   - Add new enum values to FdmSchemeDesc::FdmSchemeType:
     * FittedImplicitType
     * FittedImplicitExtrapolationType  
     * CrankNicolsonRannacherType
   - Add static factory methods:
     * FdmSchemeDesc::FittedImplicit()
     * FdmSchemeDesc::FittedImplicitExtrapolation()
     * FdmSchemeDesc::CrankNicolsonRannacher(Size rannacherSteps = 4)
   - Add cases to the rollback() switch statement for each new scheme type
   - For FittedImplicitType: use FittedImplicitScheme directly (no separate damping needed — the scheme is already monotone)
   - For FittedImplicitExtrapolationType: use FittedImplicitWithExtrapolationScheme
   - For CrankNicolsonRannacherType: use CrankNicolsonScheme with rannacherSteps from schemeDesc_.mu (repurpose mu parameter)

2. **MODIFY: `fdmblackscholessolver.hpp` / `fdmblackscholessolver.cpp`**
   - Add a boolean parameter `useFittedOperator` (default false)
   - When true, create FdmFittedBlackScholesOp instead of FdmBlackScholesOp
   - Pass through to Fdm1DimSolver

3. **MODIFY: `fdm1dimsolver.hpp` / `fdm1dimsolver.cpp`**
   - No interface changes needed; verify it works with new operator types

4. **VERIFY: `finitedifferencemodel.hpp`**
   - The rollbackImpl template should work unchanged with new scheme types since they implement the same step()/setStep() interface

Output complete modified files.

[ATTACH: architecture guide solver section, Round 1+2 headers, and the 3 solver source files]
```

---

## ROUND 4 — Old Framework Parallel Changes

**Objective:** Implement exponential fitting in the legacy `BSMOperator` / `MixedScheme` framework.

**Files to include:**
- Paper summary (Sections 4-5 on fitting factor)
- Round 1 output (ExponentiallyFittedOp for reference)
- Root-level FD files: `bsmoperator.cpp`, `tridiagonaloperator.cpp`, `mixedscheme.hpp`, `cranknicolson.hpp`, `boundarycondition.cpp`, `pde.hpp`/`pdebsm.hpp`, `operatortraits.hpp`
- `fdvanillaengine.cpp` from pricingengines

**Prompt text:**

```
You are implementing Round 4 of 6: parallel changes to the old (legacy) finite-difference framework.

QuantLib v1.23 has TWO FD frameworks:
1. New framework (methods/finitedifferences/schemes/ + operators/ + solvers/) — modified in Rounds 1-3
2. Old framework (methods/finitedifferences/ root files) — THIS ROUND

The old framework uses TridiagonalOperator + MixedScheme<> + FiniteDifferenceModel<>. We need to add exponential fitting here too.

## Paper Mathematics (same as Round 1)
Fitting factor: ρ_j = (μ_j·h/2)·coth(μ_j·h/(2·σ_j))
Tridiagonal coefficients with guaranteed M-matrix sign structure.

## Your Task

1. **NEW: `fittedbsmoperator.hpp` / `fittedbsmoperator.cpp`**
   A new TridiagonalOperator subclass analogous to BSMOperator but using exponential fitting:
   - Constructor takes grid, r, q, sigma (same signature pattern as BSMOperator)
   - Computes μ = r-q-σ²/2 and diffusion = σ²/2 at each grid point
   - Applies fitting factor ρ = (μ·h/2)·coth(μ·h/(2·diffusion))
   - Assembles tridiagonal with M-matrix signs
   - Handles limiting cases (small σ → upwind, small μ → standard)

2. **NEW: `fittedimplicit.hpp`**
   A typedef analogous to CrankNicolson<> but for fully implicit (θ=1):
   ```cpp
   template <class Operator>
   class FittedImplicit : public MixedScheme<Operator> {
       FittedImplicit(const operator_type& L, const bc_set& bcs)
       : MixedScheme<Operator>(L, 1.0, bcs) {}  // θ=1 for fully implicit
   };
   ```

3. **MODIFY: `bsmoperator.hpp` / `bsmoperator.cpp`**
   - No changes to existing code, but add a comment noting FittedBSMOperator as the preferred alternative

4. **MODIFY: `pde.hpp` / `pdebsm.hpp`**
   - Add a `PdeFittedBSM` class that works with FittedBSMOperator
   - Or extend PdeBSM to have a `fitted()` flag

5. **MODIFY: `fdvanillaengine.hpp` / `fdvanillaengine.cpp`**  
   - Add option to use FittedBSMOperator in initializeOperator()
   - Add `useFittedOperator_` flag, defaulting to false for backward compatibility

Output complete files.

[ATTACH: paper fitting sections, Round 1 output for reference, and all old-framework source files listed]
```

---

## ROUND 5 — Engine and Mesher Wiring

**Objective:** Update the user-facing pricing engines to expose the new scheme options.

**Files to include:**
- Round 1-3 output headers (interfaces only)
- From `pricingengines.xml`: documents 85 (fdblackscholesvanillaengine.cpp), 15 (fdblackscholesbarrierengine.cpp), 16 (fdblackscholesrebateengine.cpp)
- From `methods.xml`: documents 5 (fdmblackscholesmesher.cpp), 3 (concentrating1dmesher.cpp)

**Prompt text:**

```
You are implementing Round 5 of 6: wiring the new fitted schemes into user-facing engines.

## Context
Rounds 1-3 created:
- ExponentiallyFittedOp, FdmFittedBlackScholesOp (spatial operators)
- FittedImplicitScheme, FittedImplicitWithExtrapolationScheme (time stepping)
- Modified CrankNicolsonScheme with Rannacher steps
- Updated FdmBackwardSolver with new FdmSchemeDesc types
- Updated FdmBlackScholesSolver with useFittedOperator flag

[ATTACH Round 1-3 header files]

## Your Task

1. **MODIFY: `fdblackscholesvanillaengine.hpp` / `fdblackscholesvanillaengine.cpp`**
   - Add constructor parameter: `bool useFittedOperator = false`
   - Pass to FdmBlackScholesSolver
   - Update MakeFdBlackScholesVanillaEngine builder to support `.withFittedOperator(true)`
   - Default scheme when useFittedOperator=true should be FdmSchemeDesc::FittedImplicitExtrapolation()
   - Document that users can also pass FdmSchemeDesc::CrankNicolsonRannacher() for CN with damping

2. **MODIFY: `fdblackscholesbarrierengine.cpp`**
   - Add `useFittedOperator` parameter, pass through to solver
   - Barrier options are especially sensitive to CN oscillations (paper Section 9 mentions this)

3. **MODIFY: `fdblackscholesrebateengine.cpp`**
   - Same as barrier engine — add fitted operator support

4. **MODIFY: `fdmblackscholesmesher.cpp`**
   - When using fitted operators, the mesher can use a uniform grid (fitted scheme handles non-uniformity internally)
   - Add a comment noting that the concentration point is less critical with fitted operators since oscillations are suppressed

5. **VERIFY mesh quality:**
   - The paper notes that non-uniform meshes destroy CN's second-order accuracy
   - But the fitted implicit scheme has |u(x_j,t_n) - U_j^n| ≤ M(h+k) regardless of mesh structure
   - Document this advantage in comments

Output complete files.

[ATTACH: engine source files, mesher source files, Round 1-3 headers]
```

---

## ROUND 6 — Integration Verification and Keller Box Scheme

**Objective:** Add the Keller box scheme (paper Section 8.1) as an advanced alternative, and create a test/example file demonstrating all new capabilities.

**Files to include:**
- Paper summary Section 8.1 (Keller box)
- All previous round output headers
- A minimal set of existing test patterns from QuantLib

**Prompt text:**

```
You are implementing Round 6 of 6: the Keller box scheme and integration verification.

## Keller Box Scheme (Paper Section 8.1)

The Keller box scheme reduces the second-order BS PDE to a first-order system:
  a·∂u/∂x = v
  ∂v/∂x = ∂u/∂t - c·u - S

Then discretizes using box-averaging in both x and t:
  a_{j-1/2}^n · D_x^- u_j^n = v_{j-1/2}^n
  D_x^- v_j^{n-1/2} = D_t^- u_{j-1/2}^n - c_{j-1/2}^{n-1/2}·u_{j-1/2}^{n-1/2} - S_{j-1/2}^{n-1/2}

Properties: unconditionally stable, second-order for BOTH u and u_x (delta), supports piecewise smooth data, Richardson extrapolation gives two extra orders per extrapolation.

For BS: a = σ²/2, the system naturally produces both option price (u) and delta (v/a) simultaneously.

## Your Task

1. **NEW: `kellerboxscheme.hpp` / `kellerboxscheme.cpp`**
   A new scheme implementing the Keller box method:
   - Stores both u (price) and v (flux = a·∂u/∂x) arrays
   - step() method solves the coupled first-order system
   - Provides getU() and getV() for extracting both price and delta
   - Uses one-sided differences and box averages as described in paper equations (21)-(22)
   - Handles piecewise smooth BCs via mid-time enforcement (equation 24)

2. **NEW: `fdmkellerboxblackscholessolver.hpp` / `fdmkellerboxblackscholessolver.cpp`**
   A specialized solver using the Keller box scheme that:
   - Returns both valueAt(x) and derivativeX(x) with second-order accuracy
   - No need for numerical differentiation of the interpolated price for delta

3. **NEW: `example_fitted_cn.cpp`** (demonstration file)
   A comprehensive example showing:
   - European call pricing with standard CN vs fitted implicit vs Keller box
   - Comparison near the strike (where CN oscillates for delta/gamma)
   - Barrier option pricing showing improvement
   - Small volatility case showing graceful degradation
   - Print a comparison table of prices, deltas, gammas for each method

4. **Update `fdmbackwardsolver.cpp`** to add KellerBoxType to the scheme enum.

Output complete files.

[ATTACH: paper Section 8.1, all previous round headers, any needed existing solver patterns]
```

---

## Summary of All Rounds

| Round | Focus | New Files | Modified Files | Dependencies |
|-------|-------|-----------|---------------|--------------|
| 1 | Fitted spatial operators | 4 (2 hpp + 2 cpp) | 0 | None |
| 2 | Time-stepping schemes | 4-6 new + 2 modified | 2 (CN + implicit) | Round 1 |
| 3 | Solver orchestration | 0 | 3-4 | Rounds 1-2 |
| 4 | Old framework | 2-4 new + 2-3 modified | 3 | Independent (uses Round 1 math) |
| 5 | Engine wiring | 0 | 4-5 | Rounds 1-3 |
| 6 | Keller box + examples | 4-6 new | 1 | Rounds 1-5 |

**Key principle for each round:** Always attach the paper's mathematical formulations relevant to that round, the complete source files being modified, and the header files from previous rounds showing the interfaces created. Request complete file output, never diffs.

*Thinking...*

> The user wants me to draft a multi-round prompt strategy for implementing the Milev-Tagliani Crank-Nicolson variant scheme from the paper "Low Volatility Options and Numerical Diffusion of Finite Difference Schemes" into QuantLib v1.23's codebase.
> 
> Let me analyze what needs to be done:
> 
> 1. The paper describes two nonstandard FD schemes:
>    - Duffy's exponentially fitted implicit scheme
>    - Milev-Tagliani's Crank-Nicolson variant
> 
> 2. The key mathematical modifications are:
>    - The CN variant modifies the reaction term (-rV) discretization using a 6-node stencil
>    - Parameters ω1 = ω2 = -r/(16σ²)
>    - Time step restriction: Δt < 1/(r(1/2 - 2ω1) + 1/2(σM)²)
>    - The scheme produces tridiagonal matrices P and N where P is an M-matrix
>    - Positivity preserving under the parameter conditions
> 
> 3. The implementation needs to touch multiple parts of QuantLib:
>    - schemes/ - New scheme implementation
>    - operators/ - Modified operator for the reaction term
>    - solvers/ - Registration of new scheme
>    - meshers/ - Possibly for grid considerations
>    - pricingengines/ - For user-facing configuration
> 
> The critical insight is that this isn't just a parameter tweak to standard CN—it fundamentally changes how the P and N matrices are assembled because the reaction term uses a weighted 6-node stencil instead of standard averaging. I need to structure the prompts so each round builds on the previous one without exceeding context limits, starting with the mathematical core and then integrating upward through the solver chain to the pricing engines. I'm realizing the reaction term discretization is the key difference—it needs the ω-weighted stencil while the convection and diffusion parts stay standard. The scheme will be Black-Scholes specific since ω depends on r and σ, which is fine given the paper's scope.
> 
> For the multi-round structure, I should couple the scheme and operator modifications together in Round 1 since they're interdependent, then handle engine integration in Round 2, and save optional improvements like Duffy's method for Round 3. Let me now draft these prompts with precise mathematical specifications and clear code structure requirements. I need to work through the coordinate transformation carefully. The paper derives the Milev-Tagliani scheme in S-space with specific ω parameters, but QuantLib's framework uses log(S) space, so I need to re-derive those parameters for the transformed PDE where the convection coefficient becomes (r - q - σ²/2) and diffusion becomes σ²/2. The M-matrix conditions that guarantee positivity will yield different ω values in log-space, which is critical for maintaining the scheme's theoretical guarantees. /(4Δx) - D/(Δx²) < 0 if ω2 < 0 (and for typical parameters)
> 
> Actually, the conditions are complex. The paper derives ω1 = ω2 = -r/(16σ²) for S-space.
> 
> For the log(S) space implementation, we can either:
> 1. Derive the equivalent ω values for log(S) space
> 2. Work in S-space (less standard for QuantLib but matches the paper exactly)
> 3. Use the same ω formula and verify it works
> 
> I think option 3 is actually reasonable. The ω parameter controls the contribution of the reaction term to neighboring nodes. The formula ω = -r/(16σ²) ensures that the off-diagonal contributions from the reaction term are small enough to maintain M-matrix properties. This should work similarly in log(S) space since the reaction term is identical (-rV in both spaces).
> 
> Now let me check whether the M-matrix conditions hold in log(S) space with a uniform grid. I'm looking at the matrix entries for the lower diagonal, diagonal, and upper diagonal terms to verify they satisfy the necessary sign constraints for stability.
> 
> The key insight is that with ω = -r/(16σ²), the off-diagonal entries should remain negative and the diagonal entry positive, which is what we need. For reasonable parameter choices and grid spacing, these conditions should hold.
> 
> Rather than getting bogged down in the mathematical details, I should move forward with writing the implementation prompts. The paper's choice of ω should be a solid starting point, and I can verify the stability numerically if needed.
> 
> Now I'm thinking about how to structure these prompts practically. Since the user will be running separate Claude instances for each round, each one needs sufficient context to work independently. I should map out which files and information each round requires, then write the prompts accordingly.
> 
> Looking at the methods.xml file, it's quite comprehensive—covering all the finite difference components, solvers, boundary conditions, and related utilities across 86 documents. That's potentially too much to include in a single prompt, so I need to be selective about what each round actually needs rather than dumping everything at once. The user will be providing XML files, so I should reference the specific file paths within those XMLs rather than asking them to extract manually. Now I'm thinking through the scope—I'll focus on implementing the Milev-Tagliani Crank-Nicolson variant as the primary approach, with Duffy's exponentially fitted scheme as an optional addition since both are covered in the paper. For the round structure, Round 1 should handle the new scheme class, modifications to FdmBlackScholesOp, and updates to FdmBackwardSolver, with attachments including the paper, summary, and relevant code structure documentation to stay within context limits.

# Multi-Round Prompt Strategy for Implementing Milev-Tagliani CN Variant in QuantLib v1.23

Below is a complete set of prompts designed for sequential use with separate Claude 4.6 Opus instances. Each round builds on the previous one's output.

---

## Round 1 — Core Mathematical Implementation (Scheme + Operator + Solver Registration)

**Priority: HIGHEST — everything else depends on this**

**Attachments to provide:**
1. The PDF of the paper (pages 2–8 especially)
2. `low_vola_FDM_diff.md` (paper summary)
3. `methods.xml` (the QuantLib FD source code)
4. `hpp_structure_of_quantlib.md` (code structure analysis)

**Prompt:**

```
You are implementing two nonstandard finite-difference schemes from the paper
"Low Volatility Options and Numerical Diffusion of Finite Difference Schemes"
by Milev & Tagliani (Serdica Math. J. 36, 2010, pp. 223–236) into QuantLib
v1.23's finite-difference framework.

The two schemes are:
(A) The Milev-Tagliani Crank-Nicolson variant (Section 3 of the paper)
(B) Duffy's implicit exponentially fitted scheme (Section 2 of the paper)

Both target the Black-Scholes PDE:
  -∂V/∂t + rS·∂V/∂S + (1/2)σ²S²·∂²V/∂S² − rV = 0

QuantLib's new FD framework works in log-space x = ln(S), so the PDE becomes:
  -∂V/∂t + (r − q − σ²/2)·∂V/∂x + (σ²/2)·∂²V/∂x² − rV = 0

══════════════════════════════════════════════════════
SCHEME A: Milev-Tagliani CN Variant (Second-order accurate)
══════════════════════════════════════════════════════

This is a standard Crank-Nicolson for the convection and diffusion terms,
but with a MODIFIED discretization of the reaction term (−rV).

Instead of the standard CN midpoint average:
  −r·V(t+Δt/2) ≈ −r·(U^n + U^{n+1})/2

The reaction term V(t+Δt/2) is replaced by a 6-node weighted average
(paper eq. 8):

  V(t+Δt/2) = ω₁(U^n_{j-1} + U^n_{j+1}) + (1/2 − 2ω₁)U^n_j
             + ω₂(U^{n+1}_{j-1} + U^{n+1}_{j+1}) + (1/2 − 2ω₂)U^{n+1}_j

Parameters (paper eq. 9):
  ω₁ = ω₂ = −r/(16σ²)

Time-step restriction (paper eq. 9):
  Δt < 1 / [ r·(1/2 − 2ω₁) + (1/2)·(σM)² ]
where M = number of spatial nodes.

This produces the system P·U^{n+1} = N·U^n with tridiagonal matrices
(paper p. 228, adapted to log-space with Δx spacing):

Let μ = r − q − σ²/2 (drift in log-space), D = σ²/2 (diffusion coeff).

P = tridiag{
  lower: r·ω₂ + μ/(4Δx) − D/Δx²,
  diag:  1/Δt + D/Δx² + r·(1/2 − 2ω₂),
  upper: r·ω₂ − μ/(4Δx) − D/Δx²
}

N = tridiag{
  lower: −r·ω₁ − μ/(4Δx) + D/Δx²,
  diag:  1/Δt − D/Δx² − r·(1/2 − 2ω₁),
  upper: −r·ω₁ + μ/(4Δx) + D/Δx²
}

Properties (proven in paper):
• P is an irreducible diagonally dominant M-matrix ⇒ P⁻¹ > 0
• N ≥ 0 (non-negative entries) under the parameter conditions
• Positivity preserving: U^{n+1} = P⁻¹N·U^n ≥ 0 when U⁰ ≥ 0
• Discrete maximum principle: ‖U^{n+1}‖∞ ≤ ‖U^n‖∞
• Discretization error: O(Δx², Δt²)
• Artificial diffusion: (1/8)(r·Δx/σ)² · ∂²V/∂x²

══════════════════════════════════════════════════════
SCHEME B: Duffy's Exponentially Fitted Implicit Scheme (First-order accurate)
══════════════════════════════════════════════════════

This is a fully implicit scheme with a fitting factor ρ that replaces the
standard centered-difference coefficient for the second derivative.

In the general PDE form (paper eq. 2):
  -∂V/∂t + μ(x,t)·∂V/∂x + σ(x,t)·∂²V/∂x² + b(x,t)·V = 0

The fitting factor (paper eq. 3):
  ρ^{n+1}_j = (μ^{n+1}_j·h/2) · coth(μ^{n+1}_j·h / (2σ^{n+1}_j))

where h = Δx (spatial step).

For Black-Scholes in log-space:
  μ = r − q − σ²/2,  σ_coeff = σ²/2,  b = −r

The fitted operator (paper p. 226):
  L^h_k U^n_j = −(U^{n+1}_j − U^n_j)/k + μ^{n+1}_j·(U^{n+1}_{j+1} − U^{n+1}_{j-1})/(2h)
               + ρ^{n+1}_j·δ²_x U^{n+1}_j/h² + b^{n+1}_j·U^{n+1}_j

Produces: A·U^{n+1} = U^n with tridiagonal A (paper eq. 4):
  A = tridiag{
    lower: (−ρ/h² + μ/(2h))·k,
    diag:  (2ρ/h² − b + 1/k)·k,
    upper: (−ρ/h² − μ/(2h))·k
  }

Properties (paper Section 2):
• A is M-matrix ⇒ A⁻¹ > 0
• ‖A⁻¹‖∞ ≤ 1/(1+kr) < 1
• Uniform convergence: |V(Sj,tn) − U^n_j| ≤ c(h+k), c independent of σ
• When σ→0, reduces to first-order upwind with artificial diffusion
  (1/2)·μ·h·∂²V/∂x²

══════════════════════════════════════════════════════
IMPLEMENTATION REQUIREMENTS
══════════════════════════════════════════════════════

You must output COMPLETE, COMPILABLE source files (not diffs). Every file you
modify or create must be output in its entirety.

Create/modify the following files:

1. NEW FILE: methods/finitedifferences/schemes/milevtaglianischeme.hpp
2. NEW FILE: methods/finitedifferences/schemes/milevtaglianischeme.cpp
3. NEW FILE: methods/finitedifferences/schemes/exponentiallyfittedscheme.hpp
4. NEW FILE: methods/finitedifferences/schemes/exponentiallyfittedscheme.cpp
5. MODIFIED: methods/finitedifferences/operators/fdmblackscholesop.hpp
   — Add methods to expose: (a) the risk-free rate r, (b) the dividend yield q,
     (c) the volatility σ at each grid point, (d) the mesher's Δx information.
     These are needed by the new schemes to compute ω and ρ parameters.
6. MODIFIED: methods/finitedifferences/operators/fdmblackscholesop.cpp
   — Implement the new accessor methods.
7. MODIFIED: methods/finitedifferences/solvers/fdmbackwardsolver.hpp
   — Add MilevTaglianiType and ExponentiallyFittedType to FdmSchemeDesc enum.
   — Add static factory methods FdmSchemeDesc::MilevTagliani() and
     FdmSchemeDesc::ExponentiallyFitted().
8. MODIFIED: methods/finitedifferences/solvers/fdmbackwardsolver.cpp
   — Add cases for the two new scheme types in the rollback() switch.

Design notes:
• Both new schemes need access to r, q, σ, Δx — obtain these from the operator
  or from the mesher. Add whatever accessor methods are needed.
• The MilevTaglianiScheme must compute ω = −r/(16σ²) and assemble P, N as
  TripleBandLinearOps, then solve P·U^{n+1} = N·U^n using the Thomas algorithm
  (TripleBandLinearOp::solve_splitting).
• The ExponentiallyFittedScheme must compute the fitting factor ρ at each grid
  point and assemble A, then solve A·U^{n+1} = U^n.
• Both schemes must apply boundary conditions correctly (call bcSet_ methods).
• Follow the same coding conventions as existing schemes (see
  cranknicolsonscheme.cpp, impliciteulerscheme.cpp for patterns).
• Use `ext::shared_ptr`, `Disposable<Array>`, `QL_REQUIRE` macros, etc.
  exactly as in the existing codebase.

The reference files in the attached methods.xml are organized by source path.
Pay particular attention to:
• schemes/cranknicolsonscheme.cpp — pattern for a two-part scheme
• schemes/impliciteulerscheme.cpp — pattern for implicit solving
• operators/fdmblackscholesop.cpp — the BS operator you must modify
• operators/triplebandlinearop.cpp — the tridiagonal operator/solver
• solvers/fdmbackwardsolver.cpp — where scheme dispatch happens

Output each file with a clear header showing the full path.
```

---

## Round 2 — Pricing Engine Integration

**Attachments to provide:**
1. `low_vola_FDM_diff.md` (paper summary — for mathematical context)
2. `pricingengines.xml` (QuantLib pricing engine source code)
3. `hpp_structure_of_quantlib.md` (code structure analysis)
4. **All Round 1 output files** (paste the complete files from Round 1)

**Prompt:**

```
You are continuing the implementation of two nonstandard FD schemes
(Milev-Tagliani CN variant and Duffy's exponentially fitted scheme) into
QuantLib v1.23. Round 1 has already produced the core scheme implementations
and modified the operator and solver registration.

The Round 1 outputs (attached above) created:
• milevtaglianischeme.hpp/.cpp — The Milev-Tagliani CN variant scheme
• exponentiallyfittedscheme.hpp/.cpp — Duffy's exponentially fitted scheme
• Modified fdmblackscholesop.hpp/.cpp — With accessors for r, q, σ, Δx
• Modified fdmbackwardsolver.hpp/.cpp — With new scheme types registered

Now you must integrate these schemes with the user-facing pricing engines
so that users can select them when pricing options.

══════════════════════════════════════════════════════
REQUIRED MODIFICATIONS
══════════════════════════════════════════════════════

1. MODIFIED: pricingengines/vanilla/fdblackscholesvanillaengine.hpp
   — Add support for selecting MilevTagliani or ExponentiallyFitted scheme
     via the MakeFdBlackScholesVanillaEngine builder pattern.

2. MODIFIED: pricingengines/vanilla/fdblackscholesvanillaengine.cpp
   — Pass through the scheme selection to FdmSolverDesc/FdmSchemeDesc.

3. MODIFIED: pricingengines/barrier/fdblackscholesbarrierengine.hpp
4. MODIFIED: pricingengines/barrier/fdblackscholesbarrierengine.cpp
   — Same integration for barrier options.

5. MODIFIED: pricingengines/barrier/fdblackscholesrebateengine.hpp
6. MODIFIED: pricingengines/barrier/fdblackscholesrebateengine.cpp
   — Same integration for rebate engine used by barrier in-barriers.

Key design requirements:
• The default scheme should remain unchanged (Douglas or whatever was default).
• Users should be able to write:
    auto engine = MakeFdBlackScholesVanillaEngine(process)
        .withFdmSchemeDesc(FdmSchemeDesc::MilevTagliani())
        .withTGrid(100)
        .withXGrid(200);
• The MilevTagliani scheme should automatically use damping steps
  (implicit Euler) for the first few steps when payoff has discontinuities,
  just as the existing framework does for standard CN.
• The ExponentiallyFitted scheme is already fully implicit, so it may not
  need damping steps, but preserve the option.

Look at the existing fdblackscholesvanillaengine.cpp calculate() method
and its MakeFdBlackScholesVanillaEngine builder pattern in
pricingengines.xml for the pattern to follow.

Output each modified file COMPLETELY (not diffs). Include the full file
with a header showing the path.
```

---

## Round 3 — Old Framework Modifications + Boundary/Utility Enhancements

**Attachments to provide:**
1. `low_vola_FDM_diff.md` (paper summary)
2. `methods.xml` (for the old-framework files)
3. `hpp_structure_of_quantlib.md`
4. **Key Round 1 output files** (at minimum: the scheme .hpp headers, and the modified fdmbackwardsolver.hpp showing the new enum values)

**Prompt:**

```
You are continuing the implementation of the Milev-Tagliani CN variant and
Duffy's exponentially fitted scheme into QuantLib v1.23.

Rounds 1-2 implemented the new-framework versions. Now implement support in
the OLD (legacy) finite-difference framework used by FDVanillaEngine and
related engines.

══════════════════════════════════════════════════════
BACKGROUND: OLD FRAMEWORK ARCHITECTURE
══════════════════════════════════════════════════════

The old framework (in methods/finitedifferences/ root files) uses:
• TridiagonalOperator for the spatial operator (bsmoperator.cpp)
• MixedScheme<Operator> template for θ-schemes (mixedscheme.hpp)
• CrankNicolson<Operator> = MixedScheme with θ=0.5 (cranknicolson.hpp)
• FiniteDifferenceModel<Evolver> for time stepping (finitedifferencemodel.hpp)
• BoundaryCondition (Neumann/Dirichlet) for boundaries (boundarycondition.cpp)

The old-framework CN step in mixedscheme.hpp does:
  explicit part: a = (I − (1−θ)·Δt·L) · a
  implicit part: solve (I + θ·Δt·L) · x = a

══════════════════════════════════════════════════════
REQUIRED MODIFICATIONS
══════════════════════════════════════════════════════

1. NEW FILE: methods/finitedifferences/milevtaglianischeme_old.hpp
   — Template class MilevTaglianiScheme<Operator> for the old framework.
   — Must accept a TridiagonalOperator and boundary conditions.
   — At each step, must reconstruct P and N matrices using the formulas:
     ω = −r/(16σ²), then assemble P and N, solve P·x = N·a.
   — Constructor should take the operator L, boundary conditions,
     AND the risk-free rate r and volatility σ needed for ω computation.

2. NEW FILE: methods/finitedifferences/exponentiallyfittedscheme_old.hpp
   — Template class ExponentiallyFittedScheme<Operator> for the old framework.
   — Fully implicit with fitting factor ρ = (μh/2)·coth(μh/(2σ_coeff)).
   — Constructor needs drift μ, diffusion coefficient σ_coeff, and step h.

3. MODIFIED: methods/finitedifferences/bsmoperator.hpp
4. MODIFIED: methods/finitedifferences/bsmoperator.cpp
   — Add accessor methods for r, q, σ (volatility), and dx (grid spacing)
     so that the new scheme classes can compute ω and ρ parameters.

5. MODIFIED: pricingengines/vanilla/fdvanillaengine.hpp
6. MODIFIED: pricingengines/vanilla/fdvanillaengine.cpp
   — The old FDVanillaEngine base class. If possible, add a method or
     configuration option to select the Milev-Tagliani or exponentially
     fitted scheme instead of standard CrankNicolson.

The key reference files in methods.xml are:
• finitedifferences/boundarycondition.cpp (document index 1)
• finitedifferences/bsmoperator.cpp (document index 2)
• Look for cranknicolson.hpp, mixedscheme.hpp, finitedifferencemodel.hpp,
  operatortraits.hpp, tridiagonaloperator.cpp in the root FD directory
  (these are header-only or in the root — some may be in the hpp_structure doc).

Output each file COMPLETELY with full path headers.
```

---

## Round 4 — Verification Test Suite

**Attachments to provide:**
1. The PDF of the paper (especially Section 4 — Numerical Results, Examples 4.1–4.5)
2. `low_vola_FDM_diff.md` (paper summary)
3. **All output files from Rounds 1–3** (or at minimum the key headers)

**Prompt:**

```
You are writing a comprehensive test suite to verify the implementation of the
Milev-Tagliani CN variant and Duffy's exponentially fitted scheme in QuantLib
v1.23. The implementation was completed in Rounds 1–3 (attached).

══════════════════════════════════════════════════════
TEST CASES FROM THE PAPER
══════════════════════════════════════════════════════

The paper provides specific numerical examples. Implement tests that reproduce
these results and verify the key properties (positivity, oscillation-freedom,
convergence).

TEST 1: Truncated Call Option (Paper Example 4.1, Definition 4.1)
  Parameters: r=0.05, σ=0.001, T=5/12, U=70 (upper barrier), K=50 (strike)
  Payoff: f(S) = S−K if S∈[K,U], else 0
  Grid: Smax=140, ΔS=0.05, Δt=0.01
  Verify:
  (a) Standard CN produces spurious oscillations near S=70 (Fig. 1)
  (b) Milev-Tagliani variant produces oscillation-free, positive solution
  (c) Duffy's scheme produces oscillation-free positive solution
  (d) Both nonstandard schemes show numerical diffusion that decreases with
      smaller ΔS

TEST 2: Truncated Call with different r values (Paper Examples 4.3–4.4)
  Same as Test 1 but with:
  (a) r=0.01 — small r, less numerical diffusion expected
  (b) r=0.5  — large r, more numerical diffusion expected
  Verify both schemes produce positive solutions; verify diffusion increases
  with r.

TEST 3: Grid refinement (Paper Example 4.5)
  Parameters: r=0.5, σ=0.001, T=5/12, U=70, K=50, Smax=140
  Compare ΔS=0.05 vs ΔS=0.025 vs ΔS=0.01
  Verify numerical diffusion diminishes with smaller ΔS.
  At ΔS=0.01, Δt=0.001, verify Duffy and MT schemes are nearly
  indistinguishable and accurate.

TEST 4: Discrete Double Barrier Knock-Out Call (Paper Example 4.2, Def 4.2)
  Parameters: K=100, σ=0.001, T=1, r=0.05, L=95, U=110
  (Paper Fig.2 caption uses L=90; use both to check)
  ΔS=0.025, Δt=0.001
  Initial condition: V(S,0) = (S−K)⁺ · 1_{[L,U]}(S)
  Monitoring: V(S,tᵢ) = V(S,tᵢ⁻) · 1_{[L,U]}(S)
  Verify: Standard implicit/CN produce oscillations; nonstandard schemes don't.

TEST 5: Positivity preservation
  Use very low volatility σ=0.0001 with various r values.
  Verify solutions remain non-negative at all grid points for both schemes.

TEST 6: Convergence order
  Use a European vanilla call where the analytic Black-Scholes solution is known.
  Parameters: S₀=100, K=100, r=0.05, σ=0.20, T=1.0
  Measure error vs grid refinement.
  Verify: MT scheme shows O(Δx², Δt²) convergence.
  Verify: Duffy scheme shows O(Δx, Δt) convergence.

TEST 7: Comparison with analytic solution for low-volatility regime
  Parameters: S₀=60, K=50, r=0.05, σ=0.001, T=5/12
  Truncated call payoff.
  Compare numerical solution to the analytical value (which is approximately
  the discounted payoff shifted by drift).
  Verify both nonstandard schemes are more accurate than standard CN.

══════════════════════════════════════════════════════
OUTPUT REQUIREMENTS
══════════════════════════════════════════════════════

Create a single complete test file:
  test-suite/lowvolatilityfdm.cpp

Follow QuantLib's test-suite conventions:
• Use BOOST_AUTO_TEST_SUITE / BOOST_AUTO_TEST_CASE macros
• Use QL_EPSILON-level tolerances for exact properties (positivity)
• Use appropriate tolerances for convergence tests
• Print informative messages on failure
• Include setup of market data (yield curves, vol surfaces, etc.)

Also create:
  test-suite/lowvolatilityfdm.hpp
  — Declaration of the test suite

And modify:
  test-suite/quantlibtestsuite.cpp (or equivalent)
  — Register the new test suite

Output all files completely.
```

---

## Usage Instructions

**Execution order is strict:** Round 1 → Round 2 → Round 3 → Round 4.

Before each round, verify the previous round's output compiles conceptually (check that includes, class names, and method signatures are consistent).

Between rounds, if you spot inconsistencies (e.g., a method name in the scheme doesn't match what the solver expects), fix it before proceeding to the next round by asking the current-round Claude to revise.

**Context management tips:**

For Round 1 (heaviest), attach only `methods.xml`, the paper PDF, `low_vola_FDM_diff.md`, and `hpp_structure_of_quantlib.md`. Do NOT attach `pricingengines.xml`, `math.xml`, `processes.xml`, or `termstructures.xml` — they aren't needed for Round 1.

For Round 2, attach only `pricingengines.xml` plus the Round 1 outputs and the summary doc. Drop `methods.xml` since its relevant parts are now captured in Round 1's output.

For Rounds 3–4, trim aggressively — include only the specific files referenced in the prompt plus previous outputs.

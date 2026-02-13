
### SHARED PREAMBLE

This preamble is included verbatim at the start of every round. It is approximately 1,500 tokens and provides the stable context that all rounds need.

```xml
<system_context>
You are an expert C++ developer implementing improvements to QuantLib v1.23's
finite difference framework for Black-Scholes option pricing. You have deep
familiarity with QuantLib's coding conventions, class hierarchy, and the modern
Fdm* framework (not the deprecated legacy MixedScheme framework).

SOURCE-OF-TRUTH / NO-HALLUCINATION RULES (mandatory):
- Use ONLY the provided source files and prior-round outputs as factual reference for:
  * class/method signatures
  * include paths
  * base-class requirements and pure-virtual methods
  * coding conventions and header styles
  * sign conventions and operator semantics
- If any referenced file content is missing (i.e., an [INSERT FULL CONTENT ...]
  placeholder is not filled), DO NOT guess or reconstruct.
  STOP and ask the user to provide the missing file(s).
- If a requirement is ambiguous (e.g., coordinate definition ln(S) vs ln(S/K),
  whether applyTo is called at every step vs stopping times, or whether a scheme
  is copied by value), you MUST:
  (a) search the provided files for the relevant implementation detail, and
  (b) implement the behavior that matches the source.
  If still ambiguous, STOP and ask a precise clarification question.

INTERNAL VERIFICATION PROTOCOL (mandatory; do before coding; do not print):
- Verify what the 1D BS mesher stores as spatial coordinate x (ln(S), ln(S/K),
  or other) by reading the provided FdmBlackScholesMesher and/or by inspecting
  how FdmBlackScholesOp interprets mesher->location().
- Verify how FirstDerivativeOp / SecondDerivativeOp incorporate nonuniform mesh
  spacing (which mesher spacing functions they use).
- Verify how FiniteDifferenceModel handles StepCondition::applyTo:
  * is applyTo called at every time step, or only at stopping times?
  * is the evolver copied by value inside the model?
- Verify TripleBandLinearOp sign conventions before implementing any M-matrix
  diagnostics (i.e., which diagonal coefficients must be nonnegative vs
  nonpositive for the operator representation actually stored in mapT_).

CODING CONVENTIONS (mandatory):
- Namespace: all classes in `namespace QuantLib { }`
- Smart pointers: use `ext::shared_ptr`, NOT `std::shared_ptr`
- Return arrays: use `Array` (QuantLib uses move semantics; match existing
  method signatures — if similar methods return `Disposable<Array>`, do so too)
- Assertions: use `QL_REQUIRE(condition, message)` for preconditions,
  `QL_ENSURE(condition, message)` for postconditions
- Include guards: `#ifndef quantlib_<path_underscored>_hpp` /
  `#define quantlib_<path_underscored>_hpp`
- Copyright header: match the style of adjacent files in the same directory
- Includes: use angle brackets `<ql/...>` for QuantLib headers
- No `using namespace std;` in headers
- Follow existing naming: `camelCase` for methods, `camelCase_` for private
  members, `CamelCase` for classes

OUTPUT RULES:
- Produce COMPLETE files, never diffs or partial snippets
- Only output files that are new or modified in this round
- Each file must compile against unmodified QuantLib v1.23 headers (plus any
  files produced in prior rounds, which will be provided if needed)
- Include all necessary #include directives

</system_context>

<mathematical_context>
PDE in log-coordinate x and time-to-maturity τ (τ=0 at payoff, τ=T today):

  -u_τ + a(x,τ)·u_xx + b(x,τ)·u_x − r·u = 0

where a(x,τ) = σ²/2,  b(x,τ) = (r−q) − σ²/2.

IMPORTANT IMPLEMENTATION NOTE (coordinate convention):
- The notation above uses a generic log-coordinate x. The actual QuantLib
  implementation may use x=ln(S) or x=ln(S/K) (or a variant).
- In code, you MUST follow the convention used by the provided QuantLib v1.23
  sources (FdmBlackScholesMesher / FdmBlackScholesOp), especially for:
  * barrier comparisons and alignment targets
  * any conversion between S-space and mesh coordinate space
- Do not “assume ln(S/K)” just because it appears in this preamble.

EXPONENTIAL FITTING:
  θ_j = b_j·h / (2·a_j)       (local Péclet parameter)
  ρ_j = θ_j·coth(θ_j)         (fitting factor; ρ ≈ 1 + θ²/3 for |θ| < 1e-8)

Fitted stencil at interior node j (uniform grid spacing h):
  ℓ_j = a_j·ρ_j/h² − b_j/(2h)     (lower off-diagonal)
  d_j = −2·a_j·ρ_j/h² − r          (main diagonal)
  u_j⁺ = a_j·ρ_j/h² + b_j/(2h)    (upper off-diagonal)

NONUNIFORM GRID NOTE:
- The proof that exponential fitting guarantees ℓ_j ≥ 0 and u_j⁺ ≥ 0 is
  strictly for uniform spacing h.
- When using a nonuniform 1D mesher (e.g., sinh), define a per-node effective
  spacing consistent with QuantLib’s derivative operators, e.g.:
    h_j := 0.5*(dplus_j + dminus_j)
  where dplus/dminus are the local forward/backward spacings used by the
  derivative operators. This is required when computing θ_j.

IMPORTANT: Exponential fitting GUARANTEES ℓ_j ≥ 0 and u_j⁺ ≥ 0 for all
finite a_j > 0 (uniform-grid proof). Upwind fallback is only needed as a
degenerate guard when a_j ≈ 0.

This stencil decomposes as: b_j·(∂/∂x) + a_j·ρ_j·(∂²/∂x²) − r·I,
which means it can be assembled via TripleBandLinearOp::axpyb():
  mapT_.axpyb(b_array, firstDerivOp, secondDerivOp.mult(aρ_array), Array(1,-r))

RANNACHER-SMOOTHED CN: After each discontinuity event (payoff at τ=0,
monitoring projection), perform 2 implicit Euler half-steps at dt/2, then
resume Crank-Nicolson.

GREEKS: Δ = u_x/S,  Γ = (u_xx − u_x)/S²  where x is a log-coordinate.

SINH MESH: x(ξ) = x_center + c·sinh(α·(ξ − ξ₀)) for ξ ∈ [0,1], with ξ₀
found by bisection to satisfy endpoint constraints x(0)=xMin, x(1)=xMax.
</mathematical_context>
```

---

### ROUND 0 (OPTIONAL) — Preflight Fact-Check Questionnaire (no code)

Use this optional round if you want the model to audit the provided QuantLib v1.23
sources and resolve the key ambiguities BEFORE any code is generated.

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 0 (OPTIONAL): Preflight mode — do NOT write code.

Answer these questions with citations to specific identifiers/lines in the provided sources:

1) What does the 1D BS mesher store as spatial coordinate x?
   - ln(S), ln(S/K), or something else?
   - Show where this is defined (constructor / method / comment).
   - If the needed file is not provided, STOP and ask for it.

2) How do FirstDerivativeOp / SecondDerivativeOp incorporate nonuniform mesh spacing?
   - Identify which mesher spacing functions they use (dplus/dminus, locations differences, etc.)

3) Does FiniteDifferenceModel call StepCondition::applyTo at every time step,
   or only at stopping times? Does it call applyTo at the initial 'from' time?

4) Does FiniteDifferenceModel copy the evolver by value?

5) Are FdmSchemeDesc fields const or non-const in v1.23?

6) Confirm the correct sign expectations for “M-matrix satisfied” checks in TripleBandLinearOp
   for the BS operator representation stored in mapT_ (i.e., which diagonals should be
   nonnegative/nonpositive in the actual internal convention).

If any answer cannot be determined from the provided sources, say exactly what additional file is needed.
</task>

<output_specification>
Output only a short preflight report (no code).
</output_specification>
```

---

### ROUND 1 — New Fitted Operator + Sinh Mesher + Diagnostics Utility

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 1 OF 7: Create three new standalone class pairs (.hpp + .cpp). These
depend only on existing QuantLib v1.23 base classes, not on any other new files.

OUTPUT: 6 files total (3 headers + 3 implementations).
</task>

<source_files>
The following existing QuantLib headers are provided as reference for
interfaces, patterns, and conventions. Study them before writing code.

FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.hpp
[INSERT FULL CONTENT OF fdmblackscholesop.hpp]

FILE: ql/methods/finitedifferences/operators/triplebandlinearop.hpp
[INSERT FULL CONTENT OF triplebandlinearop.hpp]

FILE: ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp
[INSERT FULL CONTENT OF fdmlinearopcomposite.hpp]

FILE: ql/methods/finitedifferences/operators/firstderivativeop.hpp
[INSERT FULL CONTENT OF firstderivativeop.hpp]

FILE: ql/methods/finitedifferences/operators/secondderivativeop.hpp
[INSERT FULL CONTENT OF secondderivativeop.hpp]

FILE: ql/methods/finitedifferences/meshers/fdm1dmesher.hpp
[INSERT FULL CONTENT OF fdm1dmesher.hpp]

FILE: ql/methods/finitedifferences/meshers/concentrating1dmesher.hpp
[INSERT FULL CONTENT OF concentrating1dmesher.hpp]

FILE: ql/methods/finitedifferences/meshers/uniform1dmesher.hpp
[INSERT FULL CONTENT OF uniform1dmesher.hpp]

FILE: ql/methods/finitedifferences/meshers/fdmmesher.hpp
[INSERT FULL CONTENT OF fdmmesher.hpp]

FILE: ql/math/array.hpp
[INSERT FULL CONTENT OF array.hpp — or at minimum the class declaration]
</source_files>

<implementation_guidance>

=== FILE PAIR 1: FdmFittedBlackScholesOp ===

Path: ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp
      ql/methods/finitedifferences/operators/fdmfittedblackscholesop.cpp

This class implements FdmLinearOpComposite for a 1D Black-Scholes PDE with
exponentially fitted spatial discretization.

PRE-CODING FACT CHECKS (mandatory; do not print):
- Confirm from the provided sources what the mesher coordinate means (ln(S) vs ln(S/K))
  and mirror FdmBlackScholesOp’s convention (do not assume).
- Confirm how FirstDerivativeOp/SecondDerivativeOp handle nonuniform meshes and use
  the same local spacing definition when computing θ_i.

DESIGN DECISION (already made — do not deviate):
The operator stores FirstDerivativeOp dxMap_, SecondDerivativeOp dxxMap_, and
TripleBandLinearOp mapT_ — the same member layout as FdmBlackScholesOp. The
setTime() method computes per-node fitting factors and assembles mapT_ using
the public axpyb() method.

NONUNIFORM SPACING REQUIREMENT:
In setTime(), compute a per-node effective spacing h_i consistent with the
derivative operators, e.g.
  h_i = 0.5*(dplus_i + dminus_i)
where dplus_i/dminus_i are the forward/backward spacings along 'direction'.
If the mesher interface differs, derive h_i from locations in a way consistent
with FirstDerivativeOp/SecondDerivativeOp implementations.

  void FdmFittedBlackScholesOp::setTime(Time t1, Time t2) {
      // For each layout point i:
      //   0. Compute local spacing h_i consistent with derivative ops.
      //   1. Get sigma from volTS_ or localVol_ (mirror FdmBlackScholesOp logic)
      //   2. Compute a_i = 0.5*sigma*sigma
      //   3. Compute b_i and r_effective exactly as in FdmBlackScholesOp,
      //      including quanto adjustments if present
      //      (typically b_i = (r - q) - a_i in log-space)
      //   4. Clamp a_i to a small positive value (e.g. 1e-20) before dividing
      //   5. Compute theta_i = b_i*h_i/(2*a_i), rho_i = fittingFactor(theta_i)
      //   6. Store convection[i] = b_i
      //            fittedDiffusion[i] = a_i*rho_i
      // Then:
      //   mapT_.axpyb(convection, dxMap_,
      //              dxxMap_.mult(fittedDiffusion), Array(1, -r_effective));
  }

The static fittingFactor method:
  static Real fittingFactor(Real theta) {
      if (std::fabs(theta) < 1e-8)
          return 1.0 + theta*theta/3.0;
      return theta / std::tanh(theta);  // = theta * coth(theta)
  }

For degenerate cases (a_i ≈ 0): clamp a_i to a small positive value
(e.g., 1e-20) before computing theta. This is the ONLY fallback needed —
do not implement per-node upwind switching.

The mMatrixSatisfied() diagnostic method scans mapT_ coefficients after
setTime() and verifies all off-diagonals have the correct sign.

CRITICAL SIGN-CONVENTION REQUIREMENT:
- You MUST verify the sign convention of TripleBandLinearOp coefficients in the
  provided sources before deciding what “correct sign” means.
- Implement mMatrixSatisfied() and mMatrixViolationCount() consistent with
  the actual operator representation stored in mapT_ (not with assumptions).

Constructor signature:
  FdmFittedBlackScholesOp(
      ext::shared_ptr<FdmMesher> mesher,
      ext::shared_ptr<GeneralizedBlackScholesProcess> process,
      Real strike,
      bool localVol = false,
      Real illegalLocalVolOverwrite = -Null<Real>(),
      Size direction = 0,
      ext::shared_ptr<FdmQuantoHelper> quantoHelper
          = ext::shared_ptr<FdmQuantoHelper>());

Additional public methods beyond FdmLinearOpComposite:
  bool mMatrixSatisfied() const;   // true if all off-diags have correct sign
  Size mMatrixViolationCount() const;  // number of nodes with wrong sign (should be 0)


=== FILE PAIR 2: FdmSinhConcentrating1dMesher ===

Path: ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.hpp
      ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.cpp

Inherits from Fdm1dMesher. Constructor takes:
  FdmSinhConcentrating1dMesher(
      Real xMin, Real xMax, Size size,
      Real xCenter,
      Real alpha = 3.0,
      const std::vector<Real>& alignTargets = std::vector<Real>());

Algorithm:
1. If alpha == 0 (or sufficiently close): build a uniform grid in [xMin, xMax].
2. Else:
   a) Build uniform ξ_j = j/(size-1) for j = 0..size-1
   b) Solve for ξ₀ via bisection for asymmetric cases (xCenter != midpoint)
      using the two endpoint constraints:
        xMin = xCenter + c*sinh(alpha*(0 - ξ0))
        xMax = xCenter + c*sinh(alpha*(1 - ξ0))
      with unknowns (c, ξ0). Use bisection on ξ0 and compute c from one
      endpoint equation inside the objective. Document the objective function
      clearly in comments.
      (For symmetric xCenter == (xMin+xMax)/2: ξ0 = 0.5, c = (xMax-xMin)/(2*sinh(alpha/2)).)
   c) Compute x_j = xCenter + c*sinh(alpha*(ξ_j - ξ₀))
   d) Numerical safety: force x[0] = xMin, x[last] = xMax

Alignment targets (FIXED ORDER + SINGLE SHIFT):
3. Compute minSpacing = min_j (x[j+1]-x[j]) after base grid generation.
4. For each alignment target:
   - Find nearest node index k to that target.
   - Compute shift = target - x[k].
5. Choose at most ONE global shift: the smallest |shift| across targets
   such that |shift| < minSpacing/2.
6. Apply the shift to the entire grid if chosen.
7. After shifting, re-pin endpoints exactly:
     x[0] = xMin, x[last] = xMax
   (This fixes the logical inconsistency of shifting after pinning endpoints.)

8. Populate locations_, dplus_, dminus_ from Fdm1dMesher base class
   following the same pattern as Concentrating1dMesher's constructor.

Study Concentrating1dMesher's constructor to understand the base class
initialization pattern (how to set locations_, dplus_, dminus_).


=== FILE PAIR 3: FdmDiagnostics ===

Path: ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp
      ql/methods/finitedifferences/utilities/fdmdiagnostics.cpp

Lightweight utility class (not a step condition, not an operator).

  struct FdmDiagnosticsReport {
      Real minValue;
      Size negativeCount;
      Real oscillationScore;
      Size mMatrixViolationCount;
      Size nanCount;
  };

  class FdmDiagnostics {
    public:
      enum Level { Off, Light, Full };

      explicit FdmDiagnostics(Level level = Off);

      FdmDiagnosticsReport checkSolution(const Array& u) const;

      static Real oscillationScore(const Array& u);
      // Count sign changes in Δu_j = u[j+1]-u[j], normalize by (size-1)

      static FdmDiagnosticsReport merge(
          const FdmDiagnosticsReport& a,
          const FdmDiagnosticsReport& b);
      // Takes worst-case of each field

      Level level() const;
    private:
      Level level_;
  };

oscillationScore algorithm: iterate j = 0..size-2, compute du = u[j+1]-u[j],
track sign changes (ignoring zero-crossings where |du| < 1e-15), return
(change_count) / max(1, size-2).

</implementation_guidance>

<constraints>
- Do NOT modify any existing QuantLib files in this round
- The fitted operator must use axpyb() to populate mapT_ — do NOT attempt
  to access protected TripleBandLinearOp members
- Use ext::shared_ptr and ext::make_shared throughout
- Each .hpp file must have proper include guards and copyright header
- Each .cpp file must include its own header first, then other QL headers
- The sinh mesher must handle the edge case alpha = 0 (uniform grid)
- The diagnostics class must be thread-safe (no mutable state in static methods)
</constraints>

<output_specification>
Produce exactly 6 files in this order:
1. fdmfittedblackscholesop.hpp
2. fdmfittedblackscholesop.cpp
3. fdmsinhconcentrating1dmesher.hpp
4. fdmsinhconcentrating1dmesher.cpp
5. fdmdiagnostics.hpp
6. fdmdiagnostics.cpp

Each file must be complete and self-contained.
Begin each file with the QuantLib copyright header matching the style of
adjacent files in the same directory.
</output_specification>

<quality_checklist>
Before outputting, verify:
□ FdmFittedBlackScholesOp implements ALL pure virtual methods of FdmLinearOpComposite
□ setTime() uses axpyb() — no protected member access
□ fittingFactor returns 1.0 + theta²/3 for small theta, theta/tanh(theta) otherwise
□ Fitted operator computes local h_i consistent with derivative operators on nonuniform meshes
□ Sinh mesher populates locations_, dplus_, dminus_ correctly
□ Sinh mesher handles alpha ≈ 0 gracefully (falls back to uniform)
□ Diagnostics oscillationScore is O(N) and allocation-free
□ All ext::shared_ptr, not std::shared_ptr
□ All QL_REQUIRE for preconditions
□ Include guards match QuantLib naming convention
□ If any needed source file content is missing, STOP and ask (do not reconstruct)
</quality_checklist>
```

---

### ROUND 2 — Barrier Projection Condition + Policy Iteration

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 2 OF 7: Create two new class pairs for discrete barrier projection
and American option LCP solving. These depend only on existing QuantLib
base classes.

OUTPUT: 4 files total (2 headers + 2 implementations).
</task>

<source_files>
FILE: ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/stepconditions/fdmamericanstepcondition.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/stepcondition.hpp
[INSERT FULL CONTENT — this is the base class StepCondition<Array>]

FILE: ql/methods/finitedifferences/meshers/fdmmesher.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/operators/triplebandlinearop.hpp
[INSERT FULL CONTENT]

FILE: ql/math/array.hpp
[INSERT FULL CONTENT or class declaration]
</source_files>

<implementation_guidance>

=== FILE PAIR 1: FdmBarrierProjectionCondition ===

Path: ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp
      ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.cpp

Inherits from StepCondition<Array>. At discrete monitoring times, zeroes out
grid values for nodes outside a corridor [lowerBarrier, upperBarrier].

Constructor:
  FdmBarrierProjectionCondition(
      std::vector<Time> monitoringTimes,
      Real lowerBarrier,          // in S-space (NOT log-space)
      Real upperBarrier,          // in S-space
      ext::shared_ptr<FdmMesher> mesher,
      Size direction = 0);

FACT-CHECK (mandatory; do not print):
- You MUST verify from provided sources what mesher->location(iter, direction)
  represents (ln(S) vs ln(S/K) or other).
- If this cannot be determined from the provided sources, STOP and ask for the
  missing relevant files (e.g., FdmBlackScholesMesher/FdmBlackScholesOp).

The constructor precomputes outsideIndices_: a sorted vector of layout indices
whose spatial coordinate (in the SAME coordinate space as mesher->location(iter, direction))
falls outside [lowerBarrier, upperBarrier].

Edge cases:
- monitoringTimes empty: constructor succeeds; applyTo is always a no-op
- lowerBarrier = 0: means no lower barrier
- upperBarrier = +inf: means no upper barrier

Coordinate conversion design:
- Store barriers internally in the same coordinate system as the mesher locations:
  * If x = ln(S):     lnLower_ = log(lowerBarrier), lnUpper_ = log(upperBarrier)
  * If x = ln(S/K):   you need log(lowerBarrier/strike) etc.
    BUT strike is not in the constructor; therefore if the mesher uses ln(S/K)
    and strike cannot be obtained from the provided sources in a correct way,
    STOP and ask for a design clarification (do not guess).
- Then compare directly against mesher->location(iter, direction).

  void applyTo(Array& a, Time t) const override;
  // If t matches any monitoring time within tolerance (1e-10), set a[i]=0
  // for all i in outsideIndices_. Otherwise no-op.

  const std::vector<Time>& monitoringTimes() const;
  // Accessor for FdmStepConditionComposite to register as stopping times


=== FILE PAIR 2: FdmPolicyIterationLCP ===

Path: ql/methods/finitedifferences/utilities/fdmpolicyiteration.hpp
      ql/methods/finitedifferences/utilities/fdmpolicyiteration.cpp

IMPORTANT DESIGN DECISION: This class operates on PLAIN ARRAYS, not on
TripleBandLinearOp directly. This avoids the protected-member-access problem.

The caller provides:
- The tridiagonal system coefficients as three Arrays (lower, diag, upper)
  and a right-hand side Array
- The exercise values (payoff) as an Array
- The class solves the LCP internally using its own Thomas implementation

Interface:
  class FdmPolicyIterationLCP {
    public:
      FdmPolicyIterationLCP(
          Size maxIterations = 50,
          Real tolerance = 1e-12);

      Array solve(
          const Array& lower,
          const Array& diag,
          const Array& upper,
          const Array& rhs,
          const Array& phi) const;

      Size lastIterationCount() const;

    private:
      Size maxIterations_;
      Real tolerance_;
      mutable Size lastIterations_;

      static Array thomasSolve(
          Array lower, Array diag, Array upper, Array rhs);
  };

Policy iteration algorithm:
1. Initial guess: solve A·u = rhs via Thomas, then u = max(u, phi)
2. Loop:
   a. Determine active set: active[i] = (u[i] <= phi[i] + tol)
   b. Build modified system: for active[i], set lower[i]=0, diag[i]=1,
      upper[i]=0, rhs[i]=phi[i]
   c. Solve modified system via Thomas
   d. If no active set changes and max|u_new - u| < tol, converge
3. Final projection: u[i] = max(u[i], phi[i])

The internal thomasSolve is a standard Thomas algorithm on Array copies
(to avoid modifying caller's data).

</implementation_guidance>

<constraints>
- FdmBarrierProjectionCondition must handle the case of empty monitoringTimes
  (constructor succeeds, applyTo is always a no-op)
- FdmBarrierProjectionCondition must handle lowerBarrier = 0 (no lower barrier)
  and upperBarrier = +inf (no upper barrier)
- FdmPolicyIterationLCP must NOT depend on TripleBandLinearOp — it works with
  plain Arrays only
- FdmPolicyIterationLCP::solve must handle the case where phi is everywhere
  below the unconstrained solution (no active nodes → return unconstrained solve)
- Both classes must be const-correct (solve is const with mutable counter)
</constraints>

<output_specification>
Produce exactly 4 files:
1. fdmbarrierprojectioncondition.hpp
2. fdmbarrierprojectioncondition.cpp
3. fdmpolicyiteration.hpp
4. fdmpolicyiteration.cpp
</output_specification>

<quality_checklist>
□ FdmBarrierProjectionCondition correctly inherits StepCondition<Array>
□ applyTo is a no-op when t does not match any monitoring time
□ outsideIndices_ is computed once in constructor, not per call
□ Coordinate conversion matches verified QuantLib source convention
□ Policy iteration converges in constant iterations for typical American options
□ Thomas solver handles n=1 edge case
□ All Arrays are passed by const reference where not modified
□ No TripleBandLinearOp dependency in policy iteration
□ If coordinate convention cannot be verified from provided files, STOP and ask
</quality_checklist>
```

---

### ROUND 3 — Scheme Modifications (Damping Restart + FdmSchemeDesc)

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 3 OF 7: Modify CrankNicolsonScheme to support monitoring-restart damping,
and modify FdmSchemeDesc to carry the new configuration parameter.

These changes are internal to the scheme layer and do not depend on Rounds 1-2.

OUTPUT: 3 files (1 modified header containing FdmSchemeDesc and FdmBackwardSolver
declaration, 1 modified CN scheme header, 1 modified CN scheme implementation).
NOTE: In this round, only FdmSchemeDesc changes inside fdmbackwardsolver.hpp —
the FdmBackwardSolver class body is NOT changed until Round 4.
</task>

<source_files>
FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
[INSERT FULL CONTENT — this file contains both FdmSchemeDesc struct AND
FdmBackwardSolver class declaration]

FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.cpp
[INSERT FULL CONTENT — if available; if unavailable STOP and ask for it; do not reconstruct]

FILE: ql/methods/finitedifferences/schemes/impliciteulerscheme.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/schemes/expliciteulerscheme.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/schemes/boundaryconditionschemehelper.hpp
[INSERT FULL CONTENT]
</source_files>

<implementation_guidance>

=== FILE 1: Modified FdmSchemeDesc (in fdmbackwardsolver.hpp) ===

FdmSchemeDesc gains a new field: Size monitoringDampingSteps.

APPROACH: Check whether the existing members (type, theta, mu) are const or
non-const in the provided header. Then:

If members are NON-CONST (likely): simply add the new field with a default
  value and update all static factories to set it to 0.

If members are CONST: add a second constructor that takes the additional
  parameter, and update all factories. Keep the old constructor for backward
  compatibility (initializing new field to 0).

Add a new static factory:
  static FdmSchemeDesc CrankNicolsonWithMonitoringDamping(
      Size monitoringDampingSteps = 2);
  // Returns {CrankNicolsonType, 0.5, 0.0, monitoringDampingSteps}

CRITICAL: Do NOT change the FdmBackwardSolver class declaration or
implementation in this round. Only FdmSchemeDesc changes. Leave the
FdmBackwardSolver class body exactly as-is.


=== FILE 2: Modified CrankNicolsonScheme ===

The header gains:
1. A new constructor parameter: Size dampingHalfSteps = 0
2. Two new public methods: notifyDiscontinuity(), isDamping()
3. Three new private members: dampingHalfSteps_, dampingRemaining_,
   inDampingPhase_

ADDITIONAL ROBUSTNESS REQUIREMENT (NEW):
- Add: QL_REQUIRE(dampingHalfSteps % 2 == 0, "dampingHalfSteps must be even");
  (0 allowed). This ensures damping always corresponds to an integer number of
  full steps when implemented as two half-steps per step() call.

Modified constructor signature:
  CrankNicolsonScheme(
      Real theta,
      const ext::shared_ptr<FdmLinearOpComposite>& map,
      const bc_set& bcSet = bc_set(),
      Real relTol = 1e-8,
      ImplicitEulerScheme::SolverType solverType
          = ImplicitEulerScheme::BiCGstab,
      Size dampingHalfSteps = 0);    // NEW — defaults to 0 (no restart damping)

New methods:
  void notifyDiscontinuity();
  // Sets inDampingPhase_ = true, dampingRemaining_ = dampingHalfSteps_
  // If dampingHalfSteps_ == 0, this is a no-op

  bool isDamping() const;
  // Returns inDampingPhase_

Modified step() method logic:
  void CrankNicolsonScheme::step(array_type& a, Time t) {
      if (inDampingPhase_ && dampingRemaining_ > 0) {
          // Two implicit Euler half-steps at dt_/2
          Time halfDt = dt_ * 0.5;
          implicit_->setStep(halfDt);

          // First half-step
          implicit_->step(a, t, 1.0);

          // Second half-step (at t - halfDt, covering remaining half)
          implicit_->step(a, t - halfDt, 1.0);

          // Restore full step for future use
          implicit_->setStep(dt_);

          dampingRemaining_ -= 2;
          if (dampingRemaining_ <= 0) {
              inDampingPhase_ = false;
          }
          return;
      }

      // Standard CN: explicit half + implicit half
      if (theta_ != 1.0)
          explicit_->step(a, t, 1.0 - theta_);
      if (theta_ != 0.0)
          implicit_->step(a, t, theta_);
  }

IMPORTANT TIMING NOTE: In QuantLib's backward-in-time convention, step(a, t)
advances FROM time t backward. The implicit scheme's step(a, t, theta) solves
(I - theta*dt*L(t-dt, t)) u_new = ... So two half-steps at dt/2 from time t
go: t → t-dt/2 → t-dt, covering the same interval as one full step.

Modified setStep():
  void CrankNicolsonScheme::setStep(Time dt) {
      dt_ = dt;
      explicit_->setStep(dt);
      implicit_->setStep(dt);
      // NOTE: do not reset damping state here — setStep is called by
      // FiniteDifferenceModel before each rollback segment
  }

</implementation_guidance>

<constraints>
- The default behavior when dampingHalfSteps=0 must be IDENTICAL to the
  original CrankNicolsonScheme — no behavioral change for existing code
- notifyDiscontinuity() when dampingHalfSteps_==0 must be a no-op
- Do NOT modify any other scheme files (Douglas, Hundsdorfer, etc.)
- Do NOT modify FdmBackwardSolver in this round
- All existing static factories on FdmSchemeDesc must still compile and
  produce the same behavior
</constraints>

<output_specification>
Produce exactly 3 files:
1. fdmbackwardsolver.hpp (modified — FdmSchemeDesc struct changed,
   FdmBackwardSolver class unchanged)
2. cranknicolsonscheme.hpp (modified)
3. cranknicolsonscheme.cpp (modified — complete implementation)
</output_specification>

<quality_checklist>
□ FdmSchemeDesc::CrankNicolson() returns monitoringDampingSteps = 0
□ All other FdmSchemeDesc factories return monitoringDampingSteps = 0
□ CrankNicolsonScheme with dampingHalfSteps=0 behaves identically to original
□ notifyDiscontinuity() with dampingHalfSteps_=0 is a no-op
□ QL_REQUIRE enforces dampingHalfSteps is even
□ step() during damping performs exactly 2 implicit Euler half-steps per call
□ implicit_->setStep is restored to dt_ after damping half-steps
□ numberOfIterations() still works correctly (returns implicit's count)
□ FdmBackwardSolver class declaration is byte-for-byte identical to original
□ If cranknicolsonscheme.cpp content is missing, STOP and ask (do not reconstruct)
</quality_checklist>
```

---

### ROUND 4 — Backward Solver + Black-Scholes Solver Modifications

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 4 OF 7: Modify FdmBackwardSolver to support monitoring-restart damping
via time-segment splitting, and modify FdmBlackScholesSolver to support the
fitted operator selection.

This round depends on:
- Round 1: FdmFittedBlackScholesOp (new operator)
- Round 3: Modified CrankNicolsonScheme (with notifyDiscontinuity) and
  modified FdmSchemeDesc (with monitoringDampingSteps field)

OUTPUT: 4 files (2 modified headers + 2 modified implementations).
</task>

<source_files>
FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
[INSERT THE ROUND 3 OUTPUT VERSION — with modified FdmSchemeDesc]

FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.cpp
[INSERT FULL CONTENT if available; if unavailable STOP and ask for it; do not reconstruct]

FILE: ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp
[INSERT FULL CONTENT if available; if unavailable STOP and ask for it; do not reconstruct]

FILE: ql/methods/finitedifferences/solvers/fdm1dimsolver.hpp
[INSERT FULL CONTENT — for reference only, NOT modified]

FILE: ql/methods/finitedifferences/finitedifferencemodel.hpp
[INSERT FULL CONTENT — for reference on how rollback works internally]

ROUND 1 OUTPUT (for reference):
FILE: ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp
[INSERT ROUND 1 OUTPUT HEADER]

ROUND 3 OUTPUT (for reference):
FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp
[INSERT ROUND 3 OUTPUT HEADER]
</source_files>

<implementation_guidance>

=== FILE 1: Modified FdmBackwardSolver ===

FACT-CHECKS (mandatory; do not print):
- Read finitedifferencemodel.hpp rollback logic to confirm:
  * if the evolver is copied by value
  * when StepCondition::applyTo is called (each step and/or stopping times)
  * whether applyTo is called at the initial 'from' time (important for segment boundaries)
- Confirm time ordering conventions (from > to, backward marching).

CRITICAL DESIGN FIX (monitoring restart damping coherence):
- Round 3 adds CrankNicolsonScheme::notifyDiscontinuity() implementing true
  Rannacher restart smoothing (2 implicit Euler half-steps at dt/2).
- Round 4 MUST use that mechanism for monitoring restart damping.
- Do NOT replace it with “full implicit Euler damping segments” that bypass
  notifyDiscontinuity(); that would make Round 3 unused and would deviate
  from the intended two-half-step restart smoothing.

Baseline behavior preservation:
- When schemeDesc_.monitoringDampingSteps == 0 OR schemeDesc_.type != CrankNicolsonType:
  rollback behavior must be IDENTICAL to the original implementation.
- When monitoringDampingSteps > 0 but there are NO stopping times inside the CN phase interval:
  behavior must be identical to monitoringDampingSteps==0 (to preserve numerical identity).

Original structure reminder:
  void FdmBackwardSolver::rollback(..., Size steps, Size dampingSteps) {
      // Phase 1: damping steps with implicit Euler
      // Phase 2: main scheme steps (CN/Douglas/Hundsdorfer/etc.)
  }

NEW STRUCTURE when monitoring restart damping is enabled (CN only):
1) Phase 1 (unchanged): if dampingSteps != 0, run the original implicit Euler
   damping phase exactly as before.

2) Phase 2 (CN only): run CN from dampingTo -> to, but split into segments at
   condition_->stoppingTimes() that lie STRICTLY inside (to, dampingTo) (use tolerance).
   Let CNFrom = dampingTo, CNTo = to, CNSteps = steps - dampingSteps.

3) Build segment boundaries in descending (backward) order:
     boundaries = {CNFrom, t_stop_1, t_stop_2, ..., CNTo}
   where t_stop_k are filtered/sorted stopping times within (CNTo, CNFrom).

4) Distribute CNSteps across segments proportionally to segment length, ensuring:
   - sum(segmentSteps) == CNSteps
   - each segment with positive length gets at least 1 step
   - if impossible (too many segments vs steps), QL_REQUIRE with a clear message

5) For each segment [segFrom, segTo] with segSteps:
   - Construct a CrankNicolsonScheme with dampingHalfSteps = schemeDesc_.monitoringDampingSteps.
   - If this segment begins immediately AFTER a discontinuity time (i.e., segFrom is a stopping time
     that has just been applied at the end of the previous segment), call
       cnEvolver.notifyDiscontinuity();
     BEFORE passing cnEvolver into FiniteDifferenceModel for this segment.
   - Then rollback this segment via FiniteDifferenceModel<CrankNicolsonScheme>.

IMPORTANT SEMANTICS NOTE:
- Segmenting must not cause double-application of step conditions at segment boundaries.
  Verify FiniteDifferenceModel’s behavior from source:
  * If applyTo is only called after each step at the new time t, segment boundaries are safe.
  * If applyTo is also called at the initial segFrom time, you MUST adjust the segmentation
    logic to prevent double-application (e.g., excluding segFrom from the time list used
    by that model call, or using tolerance-based filtering).

Other scheme types:
- The switch statement and behavior for Douglas/Hundsdorfer/etc. must remain unchanged.

You may factor out a helper rollbackSegment(...) to reuse the original two-phase
(damping + main scheme) logic, but DO NOT change the behavior of the non-CN schemes.


=== FILE 2: Modified FdmBlackScholesSolver ===

Simpler change: add bool useFittedOperator_ member and constructor parameter.

Constructor gains:
  bool useFittedOperator = false

In performCalculations():
  ext::shared_ptr<FdmLinearOpComposite> op;
  if (useFittedOperator_) {
      op = ext::make_shared<FdmFittedBlackScholesOp>(
          solverDesc_.mesher, process_.currentLink(), strike_,
          localVol_, illegalLocalVolOverwrite_, 0,
          quantoHelper_.empty() ? ext::shared_ptr<FdmQuantoHelper>()
                                : quantoHelper_.currentLink());
  } else {
      op = ext::make_shared<FdmBlackScholesOp>(/* same args */);
  }
  solver_ = ext::make_shared<Fdm1DimSolver>(solverDesc_, schemeDesc_, op);

Include the new header:
  #include <ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp>

</implementation_guidance>

<constraints>
- When monitoringDampingSteps == 0: rollback behavior must be IDENTICAL to
  the original implementation (byte-for-byte same numerical results)
- When monitoringDampingSteps > 0 but no stopping times exist in the
  condition's stopping times within the CN phase interval: behavior identical
  to monitoringDampingSteps==0
- Do NOT modify FiniteDifferenceModel
- The switch statement for other scheme types (Douglas, Hundsdorfer, etc.)
  must remain unchanged
- FdmBlackScholesSolver's existing constructor signature must remain valid
  (new parameter has default value)
- fdmfittedblackscholesop.hpp must be included in the .cpp, not the .hpp
- If any required .cpp source file content is missing, STOP and ask (do not reconstruct)
</constraints>

<output_specification>
Produce exactly 4 files:
1. fdmbackwardsolver.hpp (modified — FdmBackwardSolver updated; FdmSchemeDesc from Round 3 retained)
2. fdmbackwardsolver.cpp (modified — restructured rollback for CN monitoring restart damping)
3. fdmblackscholessolver.hpp (modified — new constructor parameter)
4. fdmblackscholessolver.cpp (modified — conditional operator creation)
</output_specification>

<quality_checklist>
□ Original rollback behavior preserved when monitoringDampingSteps == 0
□ CN phase segmented at stopping times; notifyDiscontinuity used for restart smoothing
□ Step counts distributed proportionally and sum to total
□ No double-application of discrete dividend / step conditions at segment boundaries
□ Other scheme types unchanged
□ FdmBlackScholesSolver compiles with both old and new constructor forms
□ FdmFittedBlackScholesOp header is #included in the .cpp, not the .hpp
</quality_checklist>
```

---

### ROUND 5 — Wiring: Composite, Mesher Factory, Header Registrations

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 5 OF 7: Wire the new components into existing QuantLib infrastructure.
Modify supporting classes and update header registrations.

OUTPUT: 6-8 files (modified composites, factories, and all.hpp headers).
</task>

<source_files>
FILE: ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.cpp
[INSERT FULL CONTENT if available; if unavailable STOP and ask for it; do not reconstruct]

FILE: ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/meshers/fdmblackscholesmesher.cpp
[INSERT FULL CONTENT if available; if unavailable STOP and ask for it; do not reconstruct]

FILE: ql/methods/finitedifferences/operators/all.hpp
[INSERT FULL CONTENT — typically a list of #includes]

FILE: ql/methods/finitedifferences/meshers/all.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/stepconditions/all.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/utilities/all.hpp
[INSERT FULL CONTENT]

ROUND 1 OUTPUTS (headers only, for #include):
[INSERT fdmfittedblackscholesop.hpp HEADER]
[INSERT fdmsinhconcentrating1dmesher.hpp HEADER]
[INSERT fdmdiagnostics.hpp HEADER]

ROUND 2 OUTPUTS (headers only):
[INSERT fdmbarrierprojectioncondition.hpp HEADER]
[INSERT fdmpolicyiteration.hpp HEADER]
</source_files>

<implementation_guidance>

=== FILE 1: Modified FdmStepConditionComposite ===

Add a new static factory method:

  static ext::shared_ptr<FdmStepConditionComposite>
  barrierMonitoredComposite(
      const DividendSchedule& cashFlow,
      const ext::shared_ptr<Exercise>& exercise,
      const ext::shared_ptr<FdmMesher>& mesher,
      const ext::shared_ptr<FdmInnerValueCalculator>& calculator,
      const Date& refDate,
      const DayCounter& dayCounter,
      const std::vector<Date>& monitoringDates,
      Real lowerBarrier,
      Real upperBarrier);

This factory:
1. Calls the existing vanillaComposite() to get the base condition set
   (dividends, exercise, snapshot)
2. Converts monitoringDates to Times using dayCounter and refDate
3. Creates an FdmBarrierProjectionCondition with these times and barriers
4. Adds it to the condition list and its monitoring times to stopping times
5. Ensures stopping times are sorted and unique (use tolerance for Time comparisons)
6. Returns the assembled composite

IMPLEMENTATION NOTE: Look at how vanillaComposite() builds its condition list
and stopping times, then replicate the pattern with the addition of the
barrier projection condition.

Edge-case requirement:
- Empty monitoringDates must return the same result as vanillaComposite().


=== FILE 2: Modified FdmBlackScholesMesher ===

Add a static factory method:

  static ext::shared_ptr<Fdm1dMesher> createSinhMesher(
      Size size,
      const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
      Time maturity,
      Real strike,
      Real sinhAlpha = 3.0,
      Real xMinConstraint = Null<Real>(),
      Real xMaxConstraint = Null<Real>(),
      Real eps = 0.0001,
      Real scaleFactor = 1.5,
      const std::vector<Real>& barrierAlignTargets = std::vector<Real>());

This factory:
1. Computes domain [xMin, xMax] using the same volatility-based heuristic as
   the existing constructor (using eps and scaleFactor with the process's vol)
2. Overrides with constraints if provided
3. Creates FdmSinhConcentrating1dMesher with:
   - xCenter = log(process->x0()) (spot in log-space, consistent with existing mesher convention)
   - alpha = sinhAlpha
   - alignment targets = log() of each barrier target converted into the SAME coordinate as x
     (FACT-CHECK: if the mesher coordinate is ln(S/K), then alignment targets must be log(S/K),
      i.e., log(target/strike). Do not guess; follow verified convention.)
4. Returns the mesher as Fdm1dMesher pointer


=== FILES 3-6: Updated all.hpp headers ===

Each all.hpp gets one or more new #include lines:

operators/all.hpp:
  #include <ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp>

meshers/all.hpp:
  #include <ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.hpp>

stepconditions/all.hpp:
  #include <ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp>

utilities/all.hpp:
  #include <ql/methods/finitedifferences/utilities/fdmpolicyiteration.hpp>
  #include <ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp>

</implementation_guidance>

<constraints>
- The existing vanillaComposite factory must remain unchanged
- The new barrierMonitoredComposite must handle empty monitoringDates gracefully
  (return same result as vanillaComposite)
- createSinhMesher must not modify the class's existing constructor or any other method
- all.hpp additions must be in alphabetical order within existing includes
- If any needed .cpp source file content is missing, STOP and ask (do not reconstruct)
</constraints>

<output_specification>
Produce the following files (only those that change):
1. fdmstepconditioncomposite.hpp (modified — new factory declaration)
2. fdmstepconditioncomposite.cpp (modified — new factory implementation)
3. fdmblackscholesmesher.hpp (modified — new static method)
4. fdmblackscholesmesher.cpp (modified — new static method implementation)
5. operators/all.hpp (modified)
6. meshers/all.hpp (modified)
7. stepconditions/all.hpp (modified)
8. utilities/all.hpp (modified)
</output_specification>

<quality_checklist>
□ barrierMonitoredComposite with empty monitoring dates = vanillaComposite behavior
□ createSinhMesher correctly computes log-space domain from process vol
□ Barrier align targets are converted consistently with verified coordinate convention
□ all.hpp entries are alphabetically sorted
□ No circular include dependencies introduced
</quality_checklist>
```

---

### ROUND 6 — Pricing Engine Modifications

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 6 OF 7: Modify the user-facing pricing engine to expose the new
capabilities. This is the integration point where all improvements become
accessible to user code.

IMPORTANT: The file pricingengines/vanilla/fdblackscholesvanillaengine.hpp
was NOT provided in the source materials. Before implementing:
1. Verify the file exists at this path in your QuantLib v1.23 source tree
2. If it exists and its full content is provided, apply changes as described below
3. If the path differs OR content is missing, STOP and ask for the correct path + contents.
   Do NOT reconstruct from memory.

OUTPUT: 2 files (modified engine header + implementation).
</task>

<source_files>
FILE: ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp
[INSERT FULL CONTENT — if unavailable STOP and ask; do not reconstruct]

FILE: ql/pricingengines/vanilla/fdblackscholesvanillaengine.cpp
[INSERT FULL CONTENT — if unavailable STOP and ask; do not reconstruct]

ROUND 4 OUTPUTS (for reference — solver interfaces):
[INSERT fdmblackscholessolver.hpp from Round 4]
[INSERT fdmbackwardsolver.hpp from Round 3/4]

ROUND 5 OUTPUTS (for reference — factories):
[INSERT fdmstepconditioncomposite.hpp from Round 5]
[INSERT fdmblackscholesmesher.hpp from Round 5]
</source_files>

<implementation_guidance>

The engine constructor gains two parameters with backward-compatible defaults:

  FdBlackScholesVanillaEngine(
      ext::shared_ptr<GeneralizedBlackScholesProcess> process,
      Size tGrid = 100,
      Size xGrid = 100,
      Size dampingSteps = 0,
      const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
      bool localVol = false,
      Real illegalLocalVolOverwrite = -Null<Real>(),
      bool useFittedOperator = false,     // NEW
      Real sinhAlpha = 0.0);              // NEW: 0 = standard mesher

New private members:
  bool useFittedOperator_;
  Real sinhAlpha_;

Modified calculate() method:

1. MESHER SELECTION:
   If sinhAlpha_ > 0:
     mesher = FdmBlackScholesMesher::createSinhMesher(
         xGrid_, process_, maturity, strike, sinhAlpha_, ...);
   Else:
     mesher = existing construction (unchanged)

2. SOLVER CONSTRUCTION:
   Pass useFittedOperator_ to FdmBlackScholesSolver constructor.

3. SCHEME SELECTION (no change needed):
   The user passes schemeDesc which may be
   FdmSchemeDesc::CrankNicolsonWithMonitoringDamping(2).

USER EXAMPLE (include as a comment in the engine header):

  // Price a European call with all improvements enabled:
  auto engine = ext::make_shared<FdBlackScholesVanillaEngine>(
      process,
      /* tGrid */ 800, /* xGrid */ 800,
      /* dampingSteps */ 2,
      FdmSchemeDesc::CrankNicolsonWithMonitoringDamping(2),
      /* localVol */ false,
      /* illegalLocalVolOverwrite */ -Null<Real>(),
      /* useFittedOperator */ true,
      /* sinhAlpha */ 3.0);
  option.setPricingEngine(engine);
  Real price = option.NPV();

</implementation_guidance>

<constraints>
- Existing code that constructs FdBlackScholesVanillaEngine with the original
  parameter set must compile and produce identical results
- useFittedOperator=false must produce bit-identical results to the original
- sinhAlpha=0.0 must produce bit-identical results to the original
- The #include for FdmFittedBlackScholesOp is in the .cpp only (not header)
- If the engine file path differs or contents are missing, STOP and ask (no reconstruction)
</constraints>

<output_specification>
Produce 2 files:
1. fdblackscholesvanillaengine.hpp (modified)
2. fdblackscholesvanillaengine.cpp (modified)
</output_specification>

<quality_checklist>
□ All default parameter values match the original engine's behavior
□ sinhAlpha=0 and useFittedOperator=false produce identical results to original
□ The #include for FdmFittedBlackScholesOp is in the .cpp only (not header)
□ The user example in comments compiles
□ No new public dependencies leak into the engine's header
</quality_checklist>
```

---

### ROUND 7 — Test Suite

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 7 OF 7: Create a comprehensive test file that validates all improvements
from Rounds 1-6. The tests serve as both validation and regression protection.

OUTPUT: 1 file (test implementation).

FACT-CHECK (mandatory; do not print):
- Inspect the provided existing QuantLib test file(s) to confirm:
  * which framework is used (Boost.Test vs custom macros),
  * naming conventions,
  * how tests are registered.
- If the needed existing test file(s) are not provided, STOP and ask for them.

</task>

<source_files>
Provide ALL output headers from Rounds 1-6 so the test can compile:
[INSERT all .hpp files from Rounds 1-6]

Also provide the QuantLib test infrastructure headers and at least one example test file:
[INSERT any available test-suite/*.hpp or test patterns showing how QuantLib
registers test cases — if using Boost.Test, show the BOOST_AUTO_TEST_SUITE
pattern; if using QuantLib's custom framework, show that]
</source_files>

<implementation_guidance>

File: test-suite/fdmimprovedcn.cpp
(or adapt naming to match QuantLib's test convention — check existing test
files for the pattern)

The file should define a test suite with the following test cases:

[... keep ALL original test specs T1–T7 exactly as in the legacy prompt ...]

NOTE ON TEST INFRASTRUCTURE: QuantLib's test suite typically uses either
Boost.Test macros or a custom QUANTLIB_TEST_CASE registration. Check existing
test files (e.g., test-suite/europeanoption.cpp) for the pattern and replicate.

</implementation_guidance>

<constraints>
- All tests must pass with the default improved-CN configuration
- Tests must not depend on external data files
- Each test should complete in < 10 seconds on a modern machine
- Use QL_CHECK_CLOSE, QL_CHECK_SMALL, or BOOST_CHECK_CLOSE as appropriate
- T4 (DKO) tolerance should be wider than T1 (vanilla) since no closed-form
  reference exists
- T6 may need the policy iteration class to be called directly (not through
  the engine) if the engine doesn't expose American option configuration
- If test framework patterns are not provided, STOP and ask (no reconstruction)
</constraints>

<output_specification>
Produce 1 file:
1. fdmimprovedcn.cpp (complete test file)

Include a comment block at the top listing:
- Which test cases map to which failure modes (F1-F8)
- What grid sizes are used
- Expected run time
</output_specification>

<quality_checklist>
□ All 7 tests are present and have meaningful assertions
□ Reference values are computed, not hardcoded (except for very stable constants)
□ Richardson convergence test (T1) checks the order, not just the value
□ Positivity test (T2) uses FdmDiagnostics from Round 1
□ M-matrix test (T3) runs both standard and fitted operators
□ DKO test (T4) uses monitoring-restart damping
□ American test (T6) verifies complementarity explicitly
□ All tests use QuantLib's existing infrastructure (Process, TermStructure, etc.)
□ Test file compiles against QuantLib v1.23 + all Round 1-6 outputs
□ If any required test infra files are missing, STOP and ask (no reconstruction)
</quality_checklist>
```

---

### Optional ROUND 8 — Legacy Framework (if needed)

```xml
[INSERT SHARED PREAMBLE, replacing the line about not using legacy framework with:]
This round specifically targets the DEPRECATED legacy MixedScheme/
FiniteDifferenceModel framework. These changes are optional and may not be
accepted upstream. Implement only if legacy engine support is required.

<task>
ROUND 8 (OPTIONAL): Create FittedBSMOperator for the legacy framework,
and add notifyDiscontinuity() support to MixedScheme.

OUTPUT: 4 files (2 new + 2 modified).
</task>

[... unchanged from legacy prompt; but if any required legacy sources are missing, STOP and ask ...]
```

---

## Usage Instructions

**Context management:** Each round's `[INSERT ...]` placeholder should be filled with the actual file content from the QuantLib v1.23 source tree. Include only the files listed — adding unnecessary files wastes context and can confuse the model.

**Carrying forward outputs:** Rounds 4-7 reference outputs from prior rounds. When executing Round 4, include the Round 1 and Round 3 output *headers* (not implementations) in the `<source_files>` section. When executing Round 7, include ALL output headers from Rounds 1-6.

**Verification between rounds:** After each round, verify the output compiles against QuantLib v1.23 before proceeding. Common issues to check: missing `#include` directives, wrong `ext::shared_ptr` vs `std::shared_ptr`, missing `override` keywords, `Disposable<>` vs plain return types.

**If a round fails:** Re-run with the compilation errors appended to the prompt inside an `<error_context>` tag.

**FdmSchemeDesc const-ness:** Before executing Round 3, check the actual `FdmSchemeDesc` definition in your copy of `fdmbackwardsolver.hpp`. If members ARE const, add a note to the Round 3 prompt: "Members are declared const. Use approach (b): add a new constructor overload." If they are NOT const, add: "Members are non-const. Simply add the new field with a default initializer."

**Hard stop rule reminder:** If any placeholder file content is missing, STOP and ask the user to provide it; do not reconstruct or hallucinate implementations.


```xml
<!-- ================================================================== -->
<!--  SHARED PREAMBLE  (v2.1 — revised 2026-02-16)                     -->
<!--  Include verbatim at the start of every round.                     -->
<!-- ================================================================== -->

<system_context>
You are an expert C++ developer implementing two positivity-preserving finite
difference schemes from Milev & Tagliani (2010) into QuantLib 1.42-dev's modern
Fdm* framework. The schemes address spurious oscillations and negative prices
that standard Crank–Nicolson produces when σ² ≪ r and payoffs are discontinuous.

PAPER REFERENCE (corrected formulations):
Scheme 1: Fully implicit with Duffy's exponentially fitted diffusion coefficient.
Scheme 2: Crank–Nicolson variant with a modified 6-node reaction-term stencil.
Both are described in the GOLDEN REFERENCE DOCUMENT (Round 1.5 output), whose
corrected equation chain (CE-1 through CE-25) is authoritative.

SOURCE-OF-TRUTH / NO-HALLUCINATION RULES (mandatory):
- Use ONLY the provided source files and prior-round outputs as factual reference
  for class/method signatures, include paths, base-class requirements, and
  coding conventions.
- If any referenced file content is missing (i.e., an [INSERT FULL CONTENT ...]
  placeholder is not filled), DO NOT guess or reconstruct.
  STOP and ask the user to provide the missing file(s).
- If a requirement is ambiguous, search the provided files for the relevant
  implementation detail and implement the behavior that matches the source.
  If still ambiguous, STOP and ask a precise clarification question.

INTERNAL VERIFICATION PROTOCOL (mandatory; do before coding; do not print):
- Verify what the 1D BS mesher stores as spatial coordinate x (ln(S), ln(S/K),
  or other) by reading FdmBlackScholesMesher and FdmBlackScholesOp.
- Verify how FirstDerivativeOp / SecondDerivativeOp incorporate nonuniform mesh
  spacing (which mesher spacing functions they use).
- Verify the TripleBandLinearOp sign conventions: which off-diagonal coefficients
  must be nonneg for the system matrix (I − θ·dt·L) to be an M-matrix.
- Verify how FdmBackwardSolver handles step conditions (rollback segmentation at
  stopping times, and whether applyTo is called at the initial 'from' time).
- Verify that ext::make_shared<T>(...) is available in QuantLib's ext namespace.
  If it is not aliased, fall back to ext::shared_ptr<T>(new T(...)).

CODING CONVENTIONS (mandatory):
- Namespace: all classes in `namespace QuantLib { }`
- Smart pointers: use `ext::shared_ptr`, NOT `std::shared_ptr`
- Factory calls: prefer `ext::make_shared<T>(...)` if available; verify from
  provided headers before use. Fallback: `ext::shared_ptr<T>(new T(...))`.
- Return arrays: use `Array`
- Assertions: `QL_REQUIRE(condition, message)` for preconditions,
  `QL_ENSURE(condition, message)` for postconditions
- Include guards: `#ifndef quantlib_<path_underscored>_hpp`
- Copyright header: match the style of adjacent files in the same directory
- Includes: use angle brackets `<ql/...>` for QuantLib headers
- No `using namespace std;` in headers
- Naming: `camelCase` for methods, `camelCase_` for private members,
  `CamelCase` for classes

OUTPUT RULES:
- Produce COMPLETE files, never diffs or partial snippets
- Only output files that are new or modified in this round
- Each file must compile against unmodified QuantLib 1.42-dev headers plus any
  files produced in prior rounds
- Include all necessary #include directives
</system_context>

<mathematical_context>
PAPER'S ORIGINAL PDE (S-space, time-to-expiry τ; CE-1):
  −∂V/∂τ + rS·∂V/∂S + ½σ²S²·∂²V/∂S² − rV = 0

QUANTLIB LOG-SPACE PDE (x = ln(S), backward time τ):
  −u_τ + b·u_x + a·u_xx − r·u = 0
  where a = σ²/2,  b = (r − q) − σ²/2

═══ SCHEME 1: EXPONENTIALLY FITTED OPERATOR (CE-7 through CE-15) ═══

In log-space, the Péclet parameter and fitting factor at interior node j are:
  Pe_j = b · h_j / (2·a)                   (local Péclet number)
  ρ_j = Pe_j · coth(Pe_j)                  (fitting factor; ρ → 1+Pe²/3 for |Pe|→0)

Fitted stencil (replacing the standard centered operator):
  ℓ_j = a·ρ_j/h² − b/(2h)     (lower off-diagonal)
  d_j = −2·a·ρ_j/h² − r        (main diagonal)
  u_j = a·ρ_j/h² + b/(2h)     (upper off-diagonal)

NON-NEGATIVITY PROOF (for M-matrix):
  The lower off-diagonal can be rewritten as:
    ℓ_j = (b/(2h))·(coth(Pe) − 1)
  and the upper as:
    u_j = (b/(2h))·(coth(Pe) + 1)

  For b > 0 (Pe > 0): coth(Pe) > 1, so coth(Pe)−1 > 0; b/(2h) > 0. ℓ_j > 0.
                       coth(Pe) + 1 > 0; u_j > 0.
  For b < 0 (Pe < 0): coth(Pe) < −1, so coth(Pe)−1 < −2 (negative); b/(2h) < 0.
                       Product of two negatives: ℓ_j > 0.
                       coth(Pe)+1 < 0 (negative); b/(2h) < 0. u_j > 0.
  For b = 0: ρ → 1, ℓ_j = a/h² > 0, u_j = a/h² > 0.

  On UNIFORM meshes this is guaranteed for any h > 0 and any parameter values.
  On NONUNIFORM meshes (e.g. Concentrating1dMesher), the guarantee holds when
  the mesh ratio d⁺/d⁻ is moderate (typically < 3:1), which is satisfied by
  QuantLib's standard meshers. The mMatrixViolationCount() diagnostic verifies
  this at runtime.

This ensures the system matrix (I − θ·dt·L) is an M-matrix.

Artificial diffusion (CE-15, low-σ limit): ½·r·S·ΔS·V_SS in S-space.

NOTE ON PAPER TYPO (p. 227): The upwind scheme formulas on p. 227 show a
denominator of 2h, but direct computation shows the correct denominator is h.
The subsequent consistency analysis (eq. 7) and numerical diffusion formula
½μh·V_SS are correct — they correspond to the actual limiting scheme with h.

NONUNIFORM GRID NOTE: For non-uniform meshes, compute a per-node effective
spacing consistent with QuantLib's derivative operators:
  h_j := 0.5*(dplus_j + dminus_j)
At boundary nodes (first/last), dplus or dminus is Null<Real>(); skip these
nodes in per-node computations (they are handled by boundary conditions).

═══ SCHEME 2: CN VARIANT WITH MODIFIED REACTION TERM (CE-16 through CE-23) ═══

The paper's CN variant replaces −r·u_j with the 6-node stencil:
  −r·[ω·u_{j−1} + (1−2ω)·u_j + ω·u_{j+1}]
at each time level, with ω₁ = ω₂ = ω (symmetric weighting).

PARAMETER CHOICE (CE-19):
  ω = −r/(16σ²)

EFFECTIVE-DIFFUSION REFORMULATION FOR LOG-SPACE ASSEMBLY:
  The paper's CN variant applies the FULL off-diagonal reaction weight ω at
  each time level (eq. 8: each level sums to ½, but ω appears unsplit).
  Standard CN time-stepping halves the operator. To compensate, the spatial
  operator L must contain 2|rω| on its off-diagonals so that after CN's
  ½-factor, the system matrix reproduces the paper's P and N.

  This is mathematically equivalent to using an enhanced diffusion coefficient:
    a_eff = σ²/2 + r²·h²/(8·σ²)
  with the standard reaction term −r, assembled via the standard axpyb() pattern.

  PROOF: With a_eff, standard CN produces system matrix P with:
    P_lower = −0.5·Δt·(a_eff/h² − b/(2h))
    Expanding the a_eff term:
      a_eff/h² = σ²/(2h²) + r²/(8σ²)
    So: P_lower = −Δt·σ²/(4h²) − Δt·r²/(16σ²) + Δt·b/(4h)
    Since rω = −r²/(16σ²) and a = σ²/2:
      P_lower = Δt·(rω − a/(2h²) + b/(4h))
    which matches the paper's CE-17 lower entry (after factoring out 1/Δt).

  IMPORTANT: The coefficient is r²h²/(8σ²), NOT r²h²/(16σ²).
  The factor-of-2 arises because the paper places FULL ω at each CN time level
  while standard CN halves the operator contribution. Without the doubling, we
  would get r²h²/(16σ²) — which is wrong by exactly a factor of 2.

TIME-STEP CONSTRAINT (CE-20, log-space form):
  For N ≥ 0 (explicit-side diagonal non-negative):
    1 − 0.5·Δt·(2·a_eff/h² + r) ≥ 0
  ⟹  Δt < 1 / [σ²/(2h²) + r²/(8σ²) + r/2]
  This is less restrictive than the S-space form (no growing (σM)² term)
  but still very tight for small σ. For σ=0.001, r=0.05, h=0.01:
  the r²/(8σ²) = 312.5 term dominates, giving Δt < 0.0032.

Artificial diffusion (CE-23): ⅛·(r·ΔS/σ)²·V_SS in S-space.

═══ DISCRETE BARRIER MONITORING (CE-6) ═══

At each monitoring date t_i, after the time-step solve:
  U_j ← U_j · 𝟙_{[L,U]}(S_j)
This re-introduces discontinuities that the schemes must handle.

═══ M-MATRIX DIAGNOSTIC (CE-10, Proposition 4.1) ═══

For the operator L stored in mapT_, the M-matrix condition for the system
matrix (I − θ·dt·L) requires:
  All off-diagonals of L (lower_[] and upper_[]) must be NON-NEGATIVE.
  (Because system off-diag = −θ·dt·(operator off-diag), and θ·dt > 0.)

To access protected lower_[]/upper_[] arrays of TripleBandLinearOp, construct
a temporary ModTripleBandLinearOp from mapT_ and use its public accessors.

═══ NOTE ON FittedCrankNicolson COMBINATION ═══

The paper uses Scheme 1 (fitted operator) with fully-implicit time stepping ONLY.
The FittedCrankNicolson factory (Scheme 1 operator + CN time stepping) is a
natural EXTENSION not present in the original paper. It may offer improved
temporal accuracy but its M-matrix and positivity properties differ from the
paper's proven guarantees. Use with appropriate testing.
</mathematical_context>
```

---

### ROUND 0

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 0 (OPTIONAL): Preflight mode — do NOT write code.

You are given the complete set of QuantLib 1.42-dev header AND implementation
files from the methods/finitedifferences/ directory tree. Perform a systematic
audit to answer the questions below and produce a file-level implementation
roadmap.

PART A: Source Audit Questions

Answer each with citations to specific identifiers/lines in the provided sources:

1) COORDINATE CONVENTION:
   a) What does FdmBlackScholesMesher store as the spatial coordinate x?
      Is it ln(S), ln(S/K), or something else?
   b) How does FdmBlackScholesOp interpret mesher->location() in its setTime()?
   c) If the mesher uses ln(S/K), how is the strike communicated?

2) DERIVATIVE OPERATORS:
   a) How do FirstDerivativeOp / SecondDerivativeOp compute their stencils?
   b) What spacing functions do they use (dplus/dminus from the mesher)?
   c) What is the exact stencil at interior nodes? At boundary nodes?

3) OPERATOR ASSEMBLY:
   a) What does TripleBandLinearOp::axpyb(a, op1, op2, rhs) compute?
   b) What are the exact index conventions for lower_[], diag_[], upper_[]?
   c) How does mult(Array) work?

4) TIME-STEPPING (MODERN FRAMEWORK):
   a) How does FdmBackwardSolver::rollback() handle stopping times during
      the CN phase? Does it segment the rollback at stopping times?
   b) Is applyTo called at the initial 'from' time?
   c) How are scheme instances (e.g., CrankNicolsonScheme) created inside
      rollback — by value, by pointer? Do they persist across segments?

5) SCHEME CONFIGURATION:
   a) Are FdmSchemeDesc members (type, theta, mu) const or non-const?
   b) How does FdmBackwardSolver's rollback method handle damping steps?
   c) How does CrankNicolsonScheme combine explicit and implicit sub-steps?
      (It delegates to ExplicitEulerScheme::step(a,t,theta) and
       ImplicitEulerScheme::step(a,t,theta) via friend access.)
   d) Does ImplicitEulerScheme have a public setStep(Time dt) method?
      Does calling it on the implicit sub-scheme affect the explicit
      sub-scheme's internal dt_?

6) SIGN CONVENTIONS:
   a) For the operator L stored in mapT_, what sign convention do the
      off-diagonals use? (lower_[i] multiplies u[i-neighbor_below])
   b) For the system matrix (I − θ·dt·L), what sign must the off-diagonals
      of L have for the system to be an M-matrix?
      (Answer: L off-diags must be ≥ 0 so that −θ·dt·L off-diags ≤ 0.)
   c) Verify with the standard FdmBlackScholesOp: are its off-diagonals
      always non-negative? Under what grid conditions?

7) SMART POINTERS:
   a) Is ext::make_shared<T>(...) available in QuantLib's ext namespace?
   b) If not, what is the correct fallback pattern?

PART B: Coordinate Translation

Using the verified coordinate convention from Part A:

8) Translate the paper's S-space fitting factor (CE-7):
     ρ_j = (μ_j·ΔS/2)·coth(μ_j·ΔS/(2·σ_d^(j)))
   to the log-space Péclet parameter Pe_j and fitting factor in terms of
   QuantLib's operator coefficients a = σ²/2 and b = r−q−σ²/2.

9) Translate the paper's CN variant parameter choice (CE-19):
     ω = −r/(16σ²)
   to the log-space effective diffusion coefficient. Show that the correct
   value for use with standard CN time-stepping is:
     a_eff = σ²/2 + r²·h²/(8·σ²)
   Derive this by:
   a) Writing the CN variant's L operator with distributed reaction
   b) Showing that standard CN halves the off-diagonal contribution
   c) Showing that doubling the off-diagonal in L (to compensate)
      gives the additional diffusion r²h²/(8σ²), NOT r²h²/(16σ²)
   d) Verifying the resulting P_lower matches the paper's eq. (CE-17)

10) Express the M-matrix condition for both schemes in terms of the
    TripleBandLinearOp coefficients (lower_[], upper_[]) after axpyb().

PART C: File-Level Implementation Roadmap

11) Based on the audit, produce a table listing:
    - Each new file to create (path, class name, purpose)
    - Each existing file to modify (path, what changes, which round)
    - Dependencies between rounds
    - Estimated output size per round

If any answer cannot be determined from the provided sources, say exactly
what additional file is needed.
</task>

<source_files>
FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.hpp
[INSERT FULL CONTENT OF fdmblackscholesop.hpp]

FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.cpp
[INSERT FULL CONTENT OF fdmblackscholesop.cpp — CRITICAL for setTime()/axpyb()]

FILE: ql/methods/finitedifferences/operators/triplebandlinearop.hpp
[INSERT FULL CONTENT OF triplebandlinearop.hpp]

FILE: ql/methods/finitedifferences/operators/triplebandlinearop.cpp
[INSERT FULL CONTENT OF triplebandlinearop.cpp — CRITICAL for axpyb() impl]

FILE: ql/methods/finitedifferences/operators/modtriplebandlinearop.hpp
[INSERT FULL CONTENT OF modtriplebandlinearop.hpp — needed for diagnostic]

FILE: ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp
[INSERT FULL CONTENT OF fdmlinearopcomposite.hpp]

FILE: ql/methods/finitedifferences/operators/firstderivativeop.hpp
[INSERT FULL CONTENT OF firstderivativeop.hpp]

FILE: ql/methods/finitedifferences/operators/firstderivativeop.cpp
[INSERT FULL CONTENT OF firstderivativeop.cpp]

FILE: ql/methods/finitedifferences/operators/secondderivativeop.hpp
[INSERT FULL CONTENT OF secondderivativeop.hpp]

FILE: ql/methods/finitedifferences/operators/secondderivativeop.cpp
[INSERT FULL CONTENT OF secondderivativeop.cpp]

FILE: ql/methods/finitedifferences/meshers/fdm1dmesher.hpp
[INSERT FULL CONTENT OF fdm1dmesher.hpp]

FILE: ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp
[INSERT FULL CONTENT OF fdmblackscholesmesher.hpp]

FILE: ql/methods/finitedifferences/meshers/fdmblackscholesmesher.cpp
[INSERT FULL CONTENT OF fdmblackscholesmesher.cpp — CRITICAL for coordinate]

FILE: ql/methods/finitedifferences/meshers/fdmmesher.hpp
[INSERT FULL CONTENT OF fdmmesher.hpp]

FILE: ql/methods/finitedifferences/meshers/fdmmeshercomposite.cpp
[INSERT FULL CONTENT — location()/dplus()/dminus() delegation]

FILE: ql/methods/finitedifferences/meshers/concentrating1dmesher.hpp
[INSERT FULL CONTENT OF concentrating1dmesher.hpp]

FILE: ql/methods/finitedifferences/meshers/uniform1dmesher.hpp
[INSERT FULL CONTENT OF uniform1dmesher.hpp]

FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp
[INSERT FULL CONTENT OF cranknicolsonscheme.hpp]

FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.cpp
[INSERT FULL CONTENT OF cranknicolsonscheme.cpp — CRITICAL for step() impl]

FILE: ql/methods/finitedifferences/schemes/impliciteulerscheme.hpp
[INSERT FULL CONTENT OF impliciteulerscheme.hpp]

FILE: ql/methods/finitedifferences/schemes/impliciteulerscheme.cpp
[INSERT FULL CONTENT — step(a,t,theta) protected method]

FILE: ql/methods/finitedifferences/schemes/expliciteulerscheme.hpp
[INSERT FULL CONTENT OF expliciteulerscheme.hpp]

FILE: ql/methods/finitedifferences/schemes/expliciteulerscheme.cpp
[INSERT FULL CONTENT — step(a,t,theta) protected method]

FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
[INSERT FULL CONTENT OF fdmbackwardsolver.hpp]

FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.cpp
[INSERT FULL CONTENT OF fdmbackwardsolver.cpp — CRITICAL for rollback()]

FILE: ql/methods/finitedifferences/finitedifferencemodel.hpp
[INSERT FULL CONTENT OF finitedifferencemodel.hpp]

FILE: ql/methods/finitedifferences/stepcondition.hpp
[INSERT FULL CONTENT OF stepcondition.hpp]

FILE: ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp
[INSERT FULL CONTENT OF fdmstepconditioncomposite.hpp]

FILE: ql/math/array.hpp
[INSERT FULL CONTENT or class declaration]
</source_files>

<output_specification>
Output a structured preflight report with:
- Part A: Numbered answers with source citations
- Part B: Translated formulas with explicit verification
- Part C: Implementation roadmap table

No code in this round.
</output_specification>
```

---

### ROUND 1

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 1 OF 7: Create the exponentially fitted spatial operator implementing
Scheme 1 from Milev & Tagliani (2010) in QuantLib's log-space framework, plus
a diagnostic utility for solution quality verification.

These depend only on existing QuantLib 1.42-dev base classes.

OUTPUT: 4 files total (2 headers + 2 implementations).
</task>

<source_files>
FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.hpp
[INSERT FULL CONTENT OF fdmblackscholesop.hpp]

FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.cpp
[INSERT FULL CONTENT OF fdmblackscholesop.cpp — CRITICAL: needed to see
how setTime() assembles mapT_ via axpyb()]

FILE: ql/methods/finitedifferences/operators/triplebandlinearop.hpp
[INSERT FULL CONTENT OF triplebandlinearop.hpp]

FILE: ql/methods/finitedifferences/operators/modtriplebandlinearop.hpp
[INSERT FULL CONTENT OF modtriplebandlinearop.hpp — needed for diagnostic]

FILE: ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp
[INSERT FULL CONTENT OF fdmlinearopcomposite.hpp]

FILE: ql/methods/finitedifferences/operators/firstderivativeop.hpp
[INSERT FULL CONTENT OF firstderivativeop.hpp]

FILE: ql/methods/finitedifferences/operators/secondderivativeop.hpp
[INSERT FULL CONTENT OF secondderivativeop.hpp]

FILE: ql/methods/finitedifferences/meshers/fdm1dmesher.hpp
[INSERT FULL CONTENT OF fdm1dmesher.hpp]

FILE: ql/methods/finitedifferences/meshers/fdmmesher.hpp
[INSERT FULL CONTENT OF fdmmesher.hpp]

FILE: ql/math/array.hpp
[INSERT FULL CONTENT or class declaration]
</source_files>

<implementation_guidance>

=== FILE PAIR 1: FdmFittedBlackScholesOp (Scheme 1 in log-space) ===

Path: ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp
      ql/methods/finitedifferences/operators/fdmfittedblackscholesop.cpp

This class implements FdmLinearOpComposite for a 1D Black–Scholes PDE with
exponentially fitted spatial discretization.

MATHEMATICAL BASIS (log-space adaptation of CE-7):
  Pe_j = b · h_j / (2·a)          (Péclet parameter)
  ρ_j = Pe_j · coth(Pe_j)         (fitting factor; ρ ≥ 1 for all Pe)
The fitted diffusion coefficient is a·ρ_j (replacing the physical a in the
second-derivative stencil), guaranteeing non-negative off-diagonals on
uniform meshes and typical nonuniform meshes.

Same member layout as FdmBlackScholesOp (dxMap_, dxxMap_, mapT_),
same constructor signature plus any additional members for fitting.

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

The setTime() method computes per-node fitting factors and assembles mapT_:

  void FdmFittedBlackScholesOp::setTime(Time t1, Time t2) {
      // 1. Get r, q, vol exactly as FdmBlackScholesOp does
      // 2. For each layout point i:
      //    a) Skip boundary nodes where dplus or dminus is Null<Real>()
      //       (use standard coefficients or zero for those nodes)
      //    b) Compute local effective spacing:
      //       h_i = 0.5*(mesher_->dplus(iter, direction_)
      //                 + mesher_->dminus(iter, direction_))
      //    c) Compute a_i = 0.5*vol*vol  (or from local vol if applicable)
      //    d) Compute b_i = (r - q) - a_i  (drift in log-space)
      //    e) Compute Pe_i = b_i * h_i / (2.0 * a_i)
      //       Guard: clamp a_i to max(a_i, 1e-20) before dividing
      //    f) Compute rho_i = fittingFactor(Pe_i)
      //    g) Store: fittedDiffusion[i] = a_i * rho_i
      //              convection[i] = b_i
      // 3. Assemble via axpyb:
      //    mapT_.axpyb(convection, dxMap_,
      //                dxxMap_.mult(fittedDiffusion), Array(1, -r));
  }

Static fitting factor:
  static Real fittingFactor(Real Pe) {
      if (std::fabs(Pe) < 1e-8)
          return 1.0 + Pe*Pe/3.0;     // Taylor: Pe·coth(Pe) ≈ 1 + Pe²/3
      if (std::fabs(Pe) > 300.0)
          return std::fabs(Pe);        // tanh saturates to ±1 in IEEE 754
      return Pe / std::tanh(Pe);       // = Pe·coth(Pe)
  }

  NOTE ON GUARD THRESHOLDS:
  - The small-Pe guard at 1e-8 is exact to machine precision (next Taylor
    term Pe⁴/45 ≈ 2e-33 at Pe=1e-8).
  - The large-Pe guard at 300 is conservative: std::tanh returns exactly ±1.0
    for |x| ≳ 19.1 in IEEE 754 double precision, so |Pe| > 20 would suffice.
    The value 300 is kept as a harmless margin.

M-MATRIX DIAGNOSTIC METHODS (use ModTripleBandLinearOp for access):

  bool mMatrixSatisfied() const;
  // Construct ModTripleBandLinearOp from mapT_ (uses public copy constructor).
  // Check that mod.lower(i) >= 0 and mod.upper(i) >= 0 for all interior i.

  Size mMatrixViolationCount() const;
  // Returns the number of nodes where any off-diagonal is negative.
  // For the fitted operator on uniform meshes, this should ALWAYS be 0.
  // On nonuniform meshes, nonzero indicates the mesh is too skewed.

  INCLUDE: #include <ql/methods/finitedifferences/operators/modtriplebandlinearop.hpp>

Pure virtual methods from FdmLinearOpComposite: implement ALL of them,
mirroring FdmBlackScholesOp's implementations for size(), setTime(),
apply(), apply_mixed(), apply_direction(), solve_splitting(), preconditioner().


=== FILE PAIR 2: FdmDiagnostics ===

Path: ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp
      ql/methods/finitedifferences/utilities/fdmdiagnostics.cpp

Lightweight utility for verifying solution quality.

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
      // Light: minValue, negativeCount, nanCount only
      // Full: all fields including oscillationScore

      static Real oscillationScore(const Array& u);
      // Count sign changes in Δu_j = u[j+1]−u[j], ignoring |Δu| < 1e-15.
      // Normalize by max(1, size−2). Score of 0 = monotone; >0.1 = oscillating.
      // Algorithm is O(N) with no heap allocation.

      static FdmDiagnosticsReport merge(
          const FdmDiagnosticsReport& a,
          const FdmDiagnosticsReport& b);

      Level level() const;
    private:
      Level level_;
  };

</implementation_guidance>

<constraints>
- Do NOT modify any existing QuantLib files in this round
- The fitted operator must use axpyb() to populate mapT_ — do NOT attempt
  to access protected TripleBandLinearOp members for assembly
- For the fitting factor: use Pe/tanh(Pe), NOT Pe*coth(Pe), to
  avoid implementing coth manually (tanh is in <cmath>)
- The M-matrix diagnostic must use ModTripleBandLinearOp(mapT_) for access
  to the protected lower_[]/upper_[] arrays
- At boundary nodes (where dplus or dminus is Null<Real>()), skip the
  per-node fitting computation; use unmodified diffusion coefficient
- Use ext::shared_ptr throughout
- Each .hpp must have proper include guards and copyright header
- Each .cpp must include its own header first
- Diagnostics class must be thread-safe (no mutable state in static methods)
</constraints>

<output_specification>
Produce exactly 4 files in this order:
1. fdmfittedblackscholesop.hpp
2. fdmfittedblackscholesop.cpp
3. fdmdiagnostics.hpp
4. fdmdiagnostics.cpp

Each file must be complete and self-contained.
</output_specification>

<quality_checklist>
□ FdmFittedBlackScholesOp implements ALL pure virtual methods of FdmLinearOpComposite
□ setTime() uses axpyb() — no protected member access on TripleBandLinearOp
□ fittingFactor returns 1.0 + Pe²/3 for small Pe, Pe/tanh(Pe) otherwise
□ Extra guard for |Pe| > 300 to prevent tanh returning exact ±1
□ Fitting factor uses a_i = max(0.5*σ², 1e-20) to prevent division by zero
□ Per-node h_i is computed consistently with derivative operators (dplus+dminus)/2
□ Boundary nodes (Null<Real> spacing) are skipped in fitting computation
□ mMatrixSatisfied() uses ModTripleBandLinearOp for protected member access
□ mMatrixViolationCount() should return 0 on uniform meshes (mathematical guarantee)
□ Diagnostics oscillationScore is O(N) and allocation-free
□ All ext::shared_ptr, not std::shared_ptr
□ Include guards match QuantLib naming convention
□ If any needed source file content is missing, STOP and ask
</quality_checklist>
```

---

### ROUND 2

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 2 OF 7: Create the Crank–Nicolson variant spatial operator implementing
Scheme 2 from Milev & Tagliani (2010) in QuantLib's log-space framework, and
create the discrete barrier projection step condition.

These depend only on existing QuantLib 1.42-dev base classes.

OUTPUT: 4 files total (2 headers + 2 implementations).
</task>

<source_files>
FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.cpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/operators/triplebandlinearop.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/operators/modtriplebandlinearop.hpp
[INSERT FULL CONTENT — for M-matrix diagnostic]

FILE: ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/stepcondition.hpp
[INSERT FULL CONTENT — base class StepCondition<Array>]

FILE: ql/methods/finitedifferences/meshers/fdmmesher.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/meshers/fdm1dmesher.hpp
[INSERT FULL CONTENT]

FILE: ql/math/array.hpp
[INSERT FULL CONTENT or class declaration]

ROUND 1 OUTPUT (for reference):
FILE: ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp
[INSERT ROUND 1 OUTPUT HEADER]
</source_files>

<implementation_guidance>

=== FILE PAIR 1: FdmCNVariantBlackScholesOp (Scheme 2 in log-space) ===

Path: ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.hpp
      ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.cpp

MATHEMATICAL BASIS:
The paper's CN variant (eq. 8) places FULL off-diagonal reaction weight ω at
each time level. Since QuantLib's standard CN halves the operator, the operator
must contain 2|rω| on off-diagonals to produce the paper's system matrices.

This is equivalent to using an enhanced diffusion coefficient:
  a_eff = σ²/2 + r²·h²/(8·σ²)
with the standard reaction term −r. This allows assembly via the standard
axpyb() pattern with NO special off-diagonal manipulation.

IMPORTANT: The coefficient is r²h²/(8σ²), NOT r²h²/(16σ²).
The factor-of-2 arises because eq.(8) uses full ω at each time level.

Constructor signature:
  FdmCNVariantBlackScholesOp(
      ext::shared_ptr<FdmMesher> mesher,
      ext::shared_ptr<GeneralizedBlackScholesProcess> process,
      Real strike,
      bool localVol = false,
      Real illegalLocalVolOverwrite = -Null<Real>(),
      Size direction = 0,
      ext::shared_ptr<FdmQuantoHelper> quantoHelper
          = ext::shared_ptr<FdmQuantoHelper>());

The setTime() method:
  void FdmCNVariantBlackScholesOp::setTime(Time t1, Time t2) {
      // 1. Get r, q, vol exactly as FdmBlackScholesOp does
      // 2. Compute omega = -r / (16.0 * vol * vol)
      // 3. For each layout point i:
      //    a) Skip boundary nodes where dplus/dminus is Null<Real>()
      //    b) Compute h_i = 0.5*(dplus + dminus)
      //    c) Compute a_i = 0.5*vol*vol
      //    d) Compute b_i = (r - q) - a_i
      //    e) Compute effectiveDiffusion[i] = a_i + r*r*h_i*h_i/(8.0*vol*vol)
      //       NOTE: This is a_eff = σ²/2 + r²h²/(8σ²)
      //    f) Store: convection[i] = b_i
      // 4. Assemble via axpyb (STANDARD pattern — no custom off-diagonal):
      //    mapT_.axpyb(convection, dxMap_,
      //                dxxMap_.mult(effectiveDiffusion), Array(1, -r));
  }

Additional public methods:
  Real omega() const;                    // returns −r/(16σ²) last computed
  Real timestepConstraint() const;       // returns max safe Δt from CE-20

  // timestepConstraint() implements the log-space form:
  //   Δt_max = 1.0 / (a_eff_max/h_min² + r/2)
  //   where h_min is the MINIMUM effective spacing across all interior nodes,
  //   and a_eff_max is the corresponding effective diffusion at that node.
  //   Return Δt_max.
  //
  //   CAUTION: On concentrated meshes (e.g. near the strike), h_min can be
  //   extremely small, making Δt_max impractically tight. Users should monitor
  //   how many time steps this constraint implies and consider whether the
  //   mesh concentration is compatible with the CN variant scheme.

  bool mMatrixSatisfied() const;         // checks off-diagonal signs via
                                         // ModTripleBandLinearOp
  Size mMatrixViolationCount() const;


=== FILE PAIR 2: FdmBarrierProjectionCondition ===

Path: ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp
      ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.cpp

Implements the discrete barrier monitoring from CE-6:
  U_j ← U_j · 𝟙_{[L,U]}(S_j)

Constructor:
  FdmBarrierProjectionCondition(
      std::vector<Time> monitoringTimes,
      Real lowerBarrier,          // in S-space (NOT log-space); 0 means no lower
      Real upperBarrier,          // in S-space; Null<Real>() means no upper
      ext::shared_ptr<FdmMesher> mesher,
      Size direction = 0);

COORDINATE CONVERSION:
- The FdmBlackScholesMesher header comment says "1-d mesher for the
  Black-Scholes process (in ln(S))". So locations are x = ln(S).
- VERIFY from provided .cpp source before implementing. If mesher uses
  a different coordinate, STOP and ask.
- If x = ln(S): lnLower = log(lowerBarrier), lnUpper = log(upperBarrier)
- Precompute outsideIndices_ in constructor: vector of layout indices where
  location(iter, direction) < lnLower - tol OR > lnUpper + tol,
  with tol = 1e-12 for floating-point tolerance.

  void applyTo(Array& a, Time t) const override;
  // If t matches any monitoring time (within tolerance 1e-10), set a[i]=0
  // for all i in outsideIndices_. Otherwise no-op.

  const std::vector<Time>& monitoringTimes() const;

Edge cases:
- Empty monitoringTimes: constructor succeeds; applyTo always no-op
- lowerBarrier = 0: no lower barrier (only upper knockout)
- upperBarrier = Null<Real>(): no upper barrier (only lower knockout)

</implementation_guidance>

<constraints>
- FdmCNVariantBlackScholesOp must implement ALL pure virtual methods of
  FdmLinearOpComposite
- Assembly uses ONLY axpyb() with the effective diffusion coefficient — NO
  off-diagonal manipulation, NO TripleBandLinearOp::add(), NO
  ModTripleBandLinearOp for assembly (only for diagnostic)
- The effective diffusion coefficient is a + r²h²/(8σ²), NOT /(16σ²)
- M-matrix diagnostic uses ModTripleBandLinearOp for read access
- FdmBarrierProjectionCondition correctly converts S-space barriers to
  ln(S) coordinates (VERIFIED from source, not assumed)
- outsideIndices_ computed once in constructor
- Use ext::shared_ptr throughout
</constraints>

<output_specification>
Produce exactly 4 files:
1. fdmcnvariantblackscholesop.hpp
2. fdmcnvariantblackscholesop.cpp
3. fdmbarrierprojectioncondition.hpp
4. fdmbarrierprojectioncondition.cpp
</output_specification>

<quality_checklist>
□ Assembly uses a_eff = σ²/2 + r²h²/(8σ²) with standard axpyb() pattern
□ timestepConstraint() returns 1/(a_eff_max/h_min² + r/2) using minimum h
□ timestepConstraint() documents that concentrated meshes may produce very tight limits
□ mMatrixSatisfied() uses ModTripleBandLinearOp for access
□ Boundary nodes with Null<Real> spacing are skipped
□ FdmBarrierProjectionCondition converts S-space barriers to ln(S)
□ applyTo is a no-op when t does not match any monitoring time
□ outsideIndices_ computed once in constructor with tolerance
□ All edge cases handled (empty times, single barrier, etc.)
□ If coordinate convention cannot be verified, implementation STOPs and asks
</quality_checklist>
```

---

### ROUND 3

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 3 OF 7: Modify FdmSchemeDesc to carry new configuration parameters and
modify CrankNicolsonScheme to support Rannacher-style damping restart after
discrete monitoring events.

OUTPUT: 3 files (1 modified header, 1 modified scheme header, 1 modified
scheme implementation).
NOTE: Only FdmSchemeDesc changes inside fdmbackwardsolver.hpp —
the FdmBackwardSolver class body is NOT changed until Round 4.
</task>

<source_files>
FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
[INSERT FULL CONTENT — contains both FdmSchemeDesc struct AND FdmBackwardSolver]

FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.cpp
[INSERT FULL CONTENT — CRITICAL: step() implementation needed]

FILE: ql/methods/finitedifferences/schemes/impliciteulerscheme.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/schemes/expliciteulerscheme.hpp
[INSERT FULL CONTENT]
</source_files>

<implementation_guidance>

=== FILE 1: Modified FdmSchemeDesc (in fdmbackwardsolver.hpp) ===

FdmSchemeDesc has ALL CONST members. New fields must also be const.

APPROACH: Replace the existing 3-parameter constructor with a 5-parameter
constructor that has defaults for backward compatibility:

  struct FdmSchemeDesc {
      enum FdmSchemeType { HundsdorferType, DouglasType,
                           CraigSneydType, ModifiedCraigSneydType,
                           ImplicitEulerType, ExplicitEulerType,
                           MethodOfLinesType, TrBDF2Type,
                           CrankNicolsonType };

      FdmSchemeDesc(FdmSchemeType type, Real theta, Real mu,
                    Size monitoringDampingSteps = 0,
                    Size operatorType = 0);

      const FdmSchemeType type;
      const Real theta, mu;
      const Size monitoringDampingSteps;  // NEW: Rannacher half-steps after monitoring
      const Size operatorType;             // NEW: 0=standard, 1=fitted, 2=CN variant

      // ALL existing static factories — produce identical behavior
      static FdmSchemeDesc Douglas();
      static FdmSchemeDesc CrankNicolson();
      static FdmSchemeDesc ImplicitEuler();
      static FdmSchemeDesc ExplicitEuler();
      static FdmSchemeDesc CraigSneyd();
      static FdmSchemeDesc ModifiedCraigSneyd();
      static FdmSchemeDesc Hundsdorfer();
      static FdmSchemeDesc ModifiedHundsdorfer();
      static FdmSchemeDesc MethodOfLines(Real eps=0.001, Real relInitStepSize=0.01);
      static FdmSchemeDesc TrBDF2();

      // NEW static factories
      static FdmSchemeDesc FittedImplicit();
      // Returns {ImplicitEulerType, 1.0, 0.0, 0, 1}

      static FdmSchemeDesc FittedCrankNicolson();
      // Returns {CrankNicolsonType, 0.5, 0.0, 0, 1}
      // NOTE: This is an EXTENSION beyond the paper. The paper uses Scheme 1
      // (fitted operator) only with fully-implicit time stepping. CN time
      // stepping with the fitted operator may offer improved temporal accuracy
      // but its positivity properties are NOT proven by the paper.

      static FdmSchemeDesc CNVariant(Size monitoringDampingSteps = 0);
      // Returns {CrankNicolsonType, 0.5, 0.0, monitoringDampingSteps, 2}

      static FdmSchemeDesc CrankNicolsonWithDamping(Size monitoringDampingSteps = 2);
      // Returns {CrankNicolsonType, 0.5, 0.0, monitoringDampingSteps, 0}
  };

All EXISTING static factories pass 0,0 for the new parameters (via defaults).
This is fully backward-compatible.

CRITICAL: Do NOT change FdmBackwardSolver class body in this round.


=== FILES 2-3: Modified CrankNicolsonScheme ===

Add Rannacher monitoring-restart capability.

PREREQUISITE VERIFICATION (mandatory; do before coding; do not print):
  Inspect the provided ImplicitEulerScheme header to verify:
  1. ImplicitEulerScheme has a public setStep(Time dt) method.
  2. CrankNicolsonScheme is declared as a friend of ImplicitEulerScheme
     (or the three-argument step(a,t,theta) is accessible).
  3. Calling implicit_->setStep(halfDt) does NOT affect explicit_->dt_.
  If any of these cannot be confirmed, STOP and ask.

Modified constructor:
  CrankNicolsonScheme(
      Real theta,
      const ext::shared_ptr<FdmLinearOpComposite>& map,
      const bc_set& bcSet = bc_set(),
      Real relTol = 1e-8,
      ImplicitEulerScheme::SolverType solverType
          = ImplicitEulerScheme::BiCGstab,
      Size dampingHalfSteps = 0);    // NEW — must be even; 0 = no restart

New methods:
  void notifyDiscontinuity();
  // Sets inDampingPhase_ = true, dampingRemaining_ = dampingHalfSteps_
  // If dampingHalfSteps_ == 0, this is a no-op

  bool isDamping() const;

QL_REQUIRE: dampingHalfSteps must be even (0 allowed).

Modified step() method:
  void CrankNicolsonScheme::step(array_type& a, Time t) {
      if (inDampingPhase_ && dampingRemaining_ > 0) {
          // Two implicit Euler half-steps at dt/2
          Time halfDt = dt_ * 0.5;
          implicit_->setStep(halfDt);
          implicit_->step(a, t, 1.0);              // t → t−dt/2
          implicit_->step(a, t - halfDt, 1.0);     // t−dt/2 → t−dt
          implicit_->setStep(dt_);                  // restore
          dampingRemaining_ -= 2;
          if (dampingRemaining_ <= 0) inDampingPhase_ = false;
          return;
      }
      // Standard CN (unchanged):
      if (theta_ != 1.0) explicit_->step(a, t, 1.0 - theta_);
      if (theta_ != 0.0) implicit_->step(a, t, theta_);
  }

Default behavior (dampingHalfSteps=0) must be BIT-IDENTICAL to original.

</implementation_guidance>

<constraints>
- Default behavior must be IDENTICAL to original CrankNicolsonScheme
- notifyDiscontinuity() when dampingHalfSteps_==0 must be a no-op
- Do NOT modify any other scheme files
- Do NOT modify FdmBackwardSolver class in this round
- All existing FdmSchemeDesc factories produce identical behavior
- The 5-parameter constructor replaces (not supplements) the 3-parameter one
- Verify ImplicitEulerScheme::setStep() accessibility from source before coding
</constraints>

<output_specification>
Produce exactly 3 files:
1. fdmbackwardsolver.hpp (FdmSchemeDesc modified, FdmBackwardSolver unchanged)
2. cranknicolsonscheme.hpp (modified)
3. cranknicolsonscheme.cpp (modified — complete implementation)
</output_specification>

<quality_checklist>
□ FdmSchemeDesc constructor is 5-param with defaults; old 3-param call sites compile
□ All existing FdmSchemeDesc factories return operatorType=0, monitoringDampingSteps=0
□ New factories produce correct operatorType values
□ FittedCrankNicolson factory includes comment noting it is an extension beyond paper
□ CrankNicolsonScheme with dampingHalfSteps=0 behaves identically to original
□ QL_REQUIRE enforces dampingHalfSteps is even
□ step() during damping performs exactly 2 implicit Euler half-steps per call
□ implicit_->setStep() accessibility verified from source (not assumed)
□ implicit_->setStep restored to dt_ after damping half-steps
□ FdmBackwardSolver class declaration unchanged from input
</quality_checklist>
```

---

### ROUND 4

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 4 OF 7: Modify FdmBackwardSolver to support monitoring-restart damping
and operator type selection. Modify FdmBlackScholesSolver to create the
appropriate operator based on the scheme description.

Depends on Rounds 1, 2, 3.

OUTPUT: 4 files (2 modified headers + 2 modified implementations).
</task>

<source_files>
FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
[INSERT ROUND 3 OUTPUT VERSION — with modified FdmSchemeDesc]

FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.cpp
[INSERT FULL CONTENT — CRITICAL for rollback() implementation]

FILE: ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp
[INSERT FULL CONTENT — CRITICAL for performCalculations()]

FILE: ql/methods/finitedifferences/solvers/fdm1dimsolver.hpp
[INSERT FULL CONTENT — for reference only, NOT modified]

FILE: ql/methods/finitedifferences/finitedifferencemodel.hpp
[INSERT FULL CONTENT — for reference on rollback internals]

ROUND 1 OUTPUT (header only):
FILE: ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp
[INSERT ROUND 1 HEADER]

ROUND 2 OUTPUT (header only):
FILE: ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.hpp
[INSERT ROUND 2 HEADER]

ROUND 3 OUTPUT (headers only):
FILE: cranknicolsonscheme.hpp
[INSERT ROUND 3 HEADER]
</source_files>

<implementation_guidance>

=== FILE 1-2: Modified FdmBackwardSolver ===

When schemeDesc_.monitoringDampingSteps > 0 AND schemeDesc_.type == CrankNicolsonType:
1. Phase 1 (unchanged): if dampingSteps != 0, run original implicit damping.
2. Phase 2 (CN): segment the CN phase at stopping times.
   At each segment boundary that is a stopping time, apply damping restart:
   - Call cnEvolver.notifyDiscontinuity() to trigger Rannacher restart
   - Continue rollback; the scheme's step() handles the damping internally

   CONSERVATIVE POLICY (M1 fix): Apply notifyDiscontinuity() at ALL stopping
   times, not just monitoring times. The rollback has no inherent way to
   distinguish monitoring times from exercise or dividend stopping times.
   Extra damping at non-monitoring stopping times adds ~2 implicit half-steps
   per occurrence — a negligible cost that prevents missed restarts at
   monitoring discontinuities. Add a comment:
     // Conservative: apply damping restart at every stopping time.
     // Extra damping at non-monitoring times is harmless (~2 half-steps).
     // Missing damping at monitoring times can cause spurious oscillations.

When monitoringDampingSteps == 0: behavior IDENTICAL to original.
When type != CrankNicolsonType: behavior IDENTICAL to original.


=== FILE 3-4: Modified FdmBlackScholesSolver ===

Add operator selection based on schemeDesc_.operatorType:

In performCalculations():
  ext::shared_ptr<FdmLinearOpComposite> op;
  switch (schemeDesc_.operatorType) {
      case 0:  // standard
          op = ext::make_shared<FdmBlackScholesOp>(...);
          break;
      case 1:  // Scheme 1: exponentially fitted
          op = ext::make_shared<FdmFittedBlackScholesOp>(...);
          break;
      case 2:  // Scheme 2: CN variant
          op = ext::make_shared<FdmCNVariantBlackScholesOp>(...);
          break;
      default:
          QL_FAIL("Unknown operator type " << schemeDesc_.operatorType);
  }

  NOTE: If ext::make_shared is not available (verify from Round 0 or source),
  use ext::shared_ptr<T>(new T(...)) instead.

Include new headers in .cpp only (NOT in .hpp):
  #include <ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp>
  #include <ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.hpp>

</implementation_guidance>

<constraints>
- When monitoringDampingSteps==0: rollback IDENTICAL to original
- When operatorType==0: solver IDENTICAL to original
- Do NOT modify FiniteDifferenceModel
- New operator headers #included in .cpp only
- If any .cpp source file content is missing, STOP and ask
- Damping applied at ALL stopping times when enabled (conservative policy)
</constraints>

<output_specification>
Produce exactly 4 files:
1. fdmbackwardsolver.hpp (from Round 3, with FdmBackwardSolver changes if needed)
2. fdmbackwardsolver.cpp (modified rollback)
3. fdmblackscholessolver.hpp (modified if needed)
4. fdmblackscholessolver.cpp (modified — conditional operator creation)
</output_specification>

<quality_checklist>
□ Original rollback behavior preserved when monitoringDampingSteps==0
□ Original solver behavior preserved when operatorType==0
□ CN phase notifies discontinuity at ALL stopping times when damping enabled
□ Comment explains conservative damping policy
□ Other scheme types (Douglas, Hundsdorfer) unchanged
□ New headers #included in .cpp only
□ If any required .cpp source is missing, STOP and ask
</quality_checklist>
```

---

### ROUND 5

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 5 OF 7: Modify the user-facing pricing engine to expose new capabilities,
wire the barrier projection condition into the step condition composite, and
update header registrations.

OUTPUT: 6-8 files.
</task>

<source_files>
FILE: ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp
[INSERT FULL CONTENT — if unavailable STOP and ask]

FILE: ql/pricingengines/vanilla/fdblackscholesvanillaengine.cpp
[INSERT FULL CONTENT — if unavailable STOP and ask]

FILE: ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.cpp
[INSERT FULL CONTENT — if unavailable STOP and ask]

FILE: ql/methods/finitedifferences/operators/all.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/stepconditions/all.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/utilities/all.hpp
[INSERT FULL CONTENT]

ROUND 1-4 OUTPUTS (headers only, for #include paths):
[INSERT all new .hpp headers from Rounds 1-4]
</source_files>

<implementation_guidance>

=== FILE 1-2: Modified FdBlackScholesVanillaEngine ===

The operatorType and monitoringDampingSteps are carried inside FdmSchemeDesc.
Users select the scheme via the static factories:
  FdmSchemeDesc::FittedImplicit()
  FdmSchemeDesc::FittedCrankNicolson()
  FdmSchemeDesc::CNVariant()
  FdmSchemeDesc::CrankNicolsonWithDamping(2)

Existing constructor signature must still compile.


=== FILE 3-4: Modified FdmStepConditionComposite ===

Add a new static factory for barrier-monitored options:

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

WIRING REQUIREMENTS (m3 fix — make explicit):
  This factory must:
  1. Convert monitoringDates to monitoringTimes (using dayCounter and refDate).
  2. Create FdmBarrierProjectionCondition with these times + barriers.
  3. Insert ALL monitoringTimes into the stoppingTimes set (so the solver
     segments the rollback at monitoring dates and applies the projection).
  4. Associate the FdmBarrierProjectionCondition with those stopping times
     in the conditions list.
  5. Also include dividend and exercise stopping times (as vanillaComposite does).

Empty monitoringDates must return same result as vanillaComposite().

#include for FdmBarrierProjectionCondition in the .cpp file.


=== FILES 5-7: Updated all.hpp headers ===

operators/all.hpp: add
  #include <ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp>
  #include <ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.hpp>

stepconditions/all.hpp: add
  #include <ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp>

utilities/all.hpp: add
  #include <ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp>

All additions in alphabetical order.

</implementation_guidance>

<constraints>
- Existing constructor signatures must still compile
- FdmSchemeDesc::Douglas() default produces identical results to original
- all.hpp additions alphabetically sorted
- No circular include dependencies
- If engine file path or content missing, STOP and ask
- barrierMonitoredComposite must explicitly insert monitoring times as stopping times
</constraints>

<output_specification>
Produce the following files:
1. fdblackscholesvanillaengine.hpp (modified)
2. fdblackscholesvanillaengine.cpp (modified)
3. fdmstepconditioncomposite.hpp (modified)
4. fdmstepconditioncomposite.cpp (modified)
5. operators/all.hpp (modified)
6. stepconditions/all.hpp (modified)
7. utilities/all.hpp (modified)
</output_specification>

<quality_checklist>
□ Existing code compiles without changes
□ Default schemeDesc produces identical results to original
□ barrierMonitoredComposite with empty monitoring dates = vanillaComposite
□ barrierMonitoredComposite explicitly registers monitoring times as stopping times
□ all.hpp entries alphabetically sorted
□ No new public dependencies in engine .hpp
</quality_checklist>
```

---

### ROUND 6

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 6 OF 7: Create a comprehensive test file that validates implementations
against the paper's numerical examples and theoretical guarantees.

OUTPUT: 1 file (test implementation).

FACT-CHECK (mandatory; do not print):
- Inspect any provided test file to confirm framework (Boost.Test, custom
  macros) and naming/registration conventions.
- If no test file is provided, STOP and ask for one.
</task>

<source_files>
All output headers from Rounds 1-5:
[INSERT all .hpp files from Rounds 1-5]

QuantLib test infrastructure:
[INSERT at least one existing test file showing registration pattern,
e.g., test-suite/fdmlinearop.cpp or test-suite/europeanoption.cpp]
</source_files>

<implementation_guidance>

File: test-suite/fdmpositivitypreserving.cpp

=== TEST T1: Fitted Operator M-Matrix Guarantee ===
Setup: FdmFittedBlackScholesOp with parameters:
  σ = {0.001, 0.01, 0.1, 0.5, 1.0}, r = {0.01, 0.05, 0.1, 0.5}
  Grid sizes: 50, 200, 800. Both uniform and concentrated meshes.
Assertion: mMatrixViolationCount() == 0 for ALL combinations on uniform meshes.


=== TEST T2: Positivity Preservation (Scheme 1 vs Standard CN) ===
Setup: Truncated call (Definition 4.1):
  K=50, U=70, T=5/12, r=0.05, σ=0.001, S_max=140
  Log-space grid: 800 nodes, Δt=0.01

GRID RESOLUTION NOTE (M3 fix): The paper uses ~2800 S-space nodes (ΔS=0.05,
S_max=140). This test uses 800 log-space nodes, which is a QUALITATIVELY
different mesh. This test validates positivity/oscillation properties, NOT
exact reproduction of the paper's figures. For approximate reproduction of
the paper's S-space resolution, see the optional T2b variant below.

Run with:
  (a) Standard FdmBlackScholesOp + CrankNicolson
  (b) FdmFittedBlackScholesOp + ImplicitEuler (Scheme 1)

Assertions:
  (a) Standard CN: negativeCount > 0 OR oscillationScore > 0.05
  (b) Fitted: min(u) ≥ 0, negativeCount == 0, oscillationScore < 0.01


=== TEST T2b (OPTIONAL): High-Resolution Positivity Check ===
Same setup as T2, but with 4000 log-space nodes and Δt=0.005.
Same assertions, wider tolerance (oscillationScore < 0.05 for scheme 1).
Purpose: Approximate the paper's S-space grid resolution in log-space.


=== TEST T3: Positivity Preservation (Scheme 2 — CN Variant) ===
Same truncated call as T2. Run with FdmCNVariantBlackScholesOp + CN.
Assertions: min(u) ≥ −1e-10, oscillationScore < 0.01.
Additional: document that timestepConstraint() < Δt used.


=== TEST T4: Convergence Rate (Scheme 1) ===
European call (smooth payoff): S₀=100, K=100, T=1.0, r=0.05, q=0.02, σ=0.2
Reference: Black-Scholes closed form.
Grid sizes N = 50, 100, 200, 400.
Assertion: Richardson error ratio approaches 2.0 (first-order spatial).


=== TEST T5: Convergence Rate (Scheme 2) ===
Same European call. Grid sizes N = 50, 100, 200, 400.
Assertion: Richardson error ratio approaches 4.0 (second-order spatial).


=== TEST T6: Discrete Double Barrier Knock-Out ===
K=100, L=95, U=110, T=1, r=0.05, σ=0.001, monthly monitoring (12 dates).
Use FdmSchemeDesc::FittedImplicit() and FdmSchemeDesc::CNVariant(2).
Both: min(u) ≥ 0, oscillationScore < 0.01.
Standard CN without damping should show oscillations near barriers.


=== TEST T7: Artificial Diffusion Comparison ===
Truncated call at multiple grid sizes (200, 400, 800, 1600 log-space nodes).
Verify: error ratio ~2 for Scheme 1, ~4 for Scheme 2.


=== TEST T8: Backward Compatibility ===
European call with FdmSchemeDesc::Douglas() and FdmSchemeDesc::CrankNicolson().
NPV must match stored reference values to machine epsilon.

</implementation_guidance>

<constraints>
- All tests must not depend on external data files
- Each test should complete in < 10 seconds
- T2/T3: use 0.0 for strict positivity on Scheme 1; -1e-10 for Scheme 2
- T8 uses hardcoded reference values from unmodified QuantLib
- T2 and T2b document that log-space grids are NOT S-space equivalents
</constraints>

<output_specification>
Produce 1 file:
1. fdmpositivitypreserving.cpp (complete test file)
</output_specification>

<quality_checklist>
□ All 8+ tests present with meaningful assertions
□ T1 checks M-matrix guarantee across many parameter combinations
□ T2 demonstrates the problem (CN fails) and the fix (Scheme 1 works)
□ T2 documents grid resolution difference from the paper
□ T3 validates Scheme 2 and documents CE-20 constraint violation
□ T4/T5 check convergence rates via Richardson ratios
□ T6 uses barrier monitoring with both schemes
□ T7 quantifies artificial diffusion differences
□ T8 ensures backward compatibility
□ All tests use QuantLib infrastructure (Process, TermStructure, etc.)
□ Log-space grid sizes specified directly (not S-space equivalents)
□ If test framework patterns are missing, STOP and ask
</quality_checklist>
```

---

### ROUND 7

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 7 (OPTIONAL): Create a standalone S-space implementation of both schemes
exactly as described in the paper, for validation against the log-space
implementations from Rounds 1-6.

OUTPUT: 2 files (1 header + 1 implementation).
</task>

<source_files>
ROUND 1.5 OUTPUT (Golden Reference Document):
[INSERT the complete Golden Reference Document from round_1.5_output.md,
specifically §3 (Mathematical Specification) and §5 (Monolithic Pseudocode)]

ROUND 1.4 OUTPUT (Corrected Equation Chain):
[INSERT the Corrected Equation Chain CE-1 through CE-25]

QuantLib headers (minimal):
FILE: ql/math/array.hpp
[INSERT class declaration]
</source_files>

<implementation_guidance>

=== FILES: FdmSSpaceReferenceSolver ===

Path: ql/methods/finitedifferences/utilities/fdmsspacereferencesolver.hpp
      ql/methods/finitedifferences/utilities/fdmsspacereferencesolver.cpp

A standalone solver implementing BOTH schemes in S-space exactly as in the
Golden Reference Document. Does NOT use QuantLib's FDM framework.

  class FdmSSpaceReferenceSolver {
    public:
      enum Scheme { ExponentiallyFitted, CrankNicolsonVariant };

      struct Result {
          Array prices;
          Array grid;
          Size timeSteps;
          bool positivityPreserved;
          Real minPrice;
          Real oscillationScore;
      };

      FdmSSpaceReferenceSolver(
          Real r, Real sigma, Real K,
          Real T, Real Smax,
          Real deltaS, Real deltaT,
          Scheme scheme);

      Result solveTruncatedCall(Real upperCutoff) const;

      Result solveDiscreteBarrierCall(
          Real lowerBarrier, Real upperBarrier,
          const std::vector<Time>& monitoringTimes) const;

    private:
      static Array thomasSolve(
          const Array& sub, const Array& diag,
          const Array& sup, const Array& rhs);

      static Real computeRho(Real mu_j, Real sigma_d_j, Real deltaS);

      Real r_, sigma_, K_, T_, Smax_, deltaS_, deltaT_;
      Scheme scheme_;
  };

IMPLEMENTATION NOTES:
- Follow the Golden Reference §5 pseudocode line-by-line
- Scheme 1: assemble A per CE-9, solve A·U^{n+1} = U^n
- Scheme 2: assemble P and N per CE-17/CE-18 with ω = -r/(16σ²),
  solve P·U^{n+1} = N·U^n
- Apply monitoring projection per CE-6

</implementation_guidance>

<constraints>
- Does NOT use QuantLib's FDM framework
- Uses only ql/math/array.hpp for Array storage
- Thread-safe (no mutable state)
</constraints>

<output_specification>
Produce exactly 2 files:
1. fdmsspacereferencesolver.hpp
2. fdmsspacereferencesolver.cpp
</output_specification>

<quality_checklist>
□ Thomas algorithm matches §3(g) of Golden Reference
□ Fitting factor matches CE-7 with all numerical guards
□ Scheme 1 matrix A matches CE-9 exactly
□ Scheme 2 matrices P, N match CE-17/CE-18 with ω from CE-19
□ Monitoring projection matches CE-6
□ Handles CE-20 constraint violation gracefully (warns, does not crash)
</quality_checklist>
```

---

## Usage Instructions

**Prerequisite: Gather QuantLib 1.42-dev sources.** Before starting, collect the actual file contents for every `[INSERT ...]` placeholder. The critical `.cpp` files are listed per round.

**Execution order:**

| Step | Round | Depends On | Key Output |
|:---|:---|:---|:---|
| 1 | Round 0 (optional) | None | Audit report + roadmap |
| 2 | Round 1 | None | FdmFittedBlackScholesOp + FdmDiagnostics |
| 3 | Round 2 | Round 0 or 1 (coordinate verification) | FdmCNVariantBlackScholesOp + FdmBarrierProjectionCondition |
| 4 | Round 3 | None | Modified CrankNicolsonScheme + FdmSchemeDesc |
| 5 | Round 4 | R1, R2, R3 | Modified FdmBackwardSolver + FdmBlackScholesSolver |
| 6 | Round 5 | R1-R4 | Modified engine + wiring + all.hpp |
| 7 | Round 6 | R1-R5 | Test suite |
| 8 | Round 7 (optional) | None | S-space reference solver |

Rounds 1 and 3 are independent and can run in parallel. **Round 2 requires coordinate verification** (from Round 0 or Round 1 output confirming x = ln(S)). Rounds 4–6 are strictly sequential.

**Carrying forward outputs:** When executing Round 4+, include prior rounds' output *headers* (not implementations) in the `<source_files>` section.

**Verification between rounds:** After each round, verify compilation against QuantLib 1.42-dev. Common issues: missing `#include`, `ext::shared_ptr` vs `std::shared_ptr`, missing `override`.

**If a round fails:** Re-run with compilation errors appended inside an `<error_context>` tag.

**Key implementation notes:**

| Topic | Guidance |
|:---|:---|
| **FdmSchemeDesc const-ness** | Members ARE const. The 5-parameter constructor with defaults (Round 3) is the correct approach. |
| **TripleBandLinearOp API** | `add()` and `mult()` exist per the header. `ModTripleBandLinearOp` provides public read/write accessors to protected arrays. For Scheme 2, use the effective-diffusion approach (no `add()` or `ModTripleBandLinearOp` needed for assembly — only for diagnostic). |
| **ext::make_shared** | Verify availability in Round 0. If unavailable, use `ext::shared_ptr<T>(new T(...))`. |
| **Damping policy** | Conservative: damp at ALL stopping times when `monitoringDampingSteps > 0` (Round 4). Cannot distinguish monitoring from other stopping times inside rollback. |
| **FittedCrankNicolson** | An extension NOT from the paper. Note this in documentation and test with care. |
| **Hard stop rule** | If any placeholder file content is missing, the model must STOP and ask — never reconstruct from memory. |

**Key mathematical references across rounds:**

| Concept | Paper Reference | Corrected Value | Round |
|:---|:---|:---|:---|
| Fitting factor | §2, Eq.(3) | Pe·coth(Pe) in log-space | R1 |
| Effective diffusion (CN variant) | §3, Eq.(8)+(9) | a + r²h²/(8σ²) — **NOT** /(16σ²) | R2 |
| Time-step constraint (log-space) | §3, Eq.(9) | 1/(σ²/(2h²) + r²/(8σ²) + r/2) | R2 |
| M-matrix guarantee | §2, p.226 | Uniform meshes: unconditional. Nonuniform: verify at runtime | R1, R6 |
| Parameter ω | §3, Eq.(9) | −r/(16σ²) (S-space value, correct) | R2 |
| Monitoring projection | §4, Eq.(13) | U_j ← U_j · 𝟙_{[L,U]}(S_j) | R2, R6 |
| Paper typo | §2, p.227 | Upwind denominator should be h, not 2h | Reference only |


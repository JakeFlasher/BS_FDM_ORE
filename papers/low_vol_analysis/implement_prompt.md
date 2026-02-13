*Thinking...*

> The user wants me to help them draft a detailed, multi-round prompting framework for implementing the algorithms from the Milev & Tagliani paper into QuantLib v1.23. They've provided:
> 
> 1. An existing multi-round plan template (`consolidated_multi_round_plan_v2.md`) - this is the format/template to follow
> 2. The QuantLib HPP structure document (`hpp_structure_of_quantlib.md`) - detailed code reading of QuantLib v1.23
> 3. The paper summary (`paper_summary_low_vola_FDM_diff.md`) - summary of the Milev & Tagliani 2010 paper
> 
> 
> 4. Rounds 1.1-1.5 covering mathematical extraction, audit, gap resolution, verification, and golden reference documentation
> 
> The existing template implements log-space exponentially fitted schemes with sinh mesh concentration and policy iteration for American options, but the Milev & Tagliani paper operates in S-space with different exponential fitting and scheme approaches that need to be carefully distinguished and extracted.
> 
> The core challenge is bridging the paper's S-space uniform grid formulation with QuantLib's log-space FDM framework, which requires understanding how the paper's exponentially fitted implicit scheme and Crank-Nicolson variant with 6-node stencil discretization map into the existing template's architecture.
> 
> I need to design a multi-round framework that systematically builds from mathematical foundations through QuantLib implementation, where each round addresses a specific layer: translating the paper's schemes to log-space coordinates, implementing the fitted operator within QuantLib's composite operator pattern, handling the modified reaction term discretization, and finally integrating discrete barrier monitoring and diagnostic checks. Scheme 2 actually modifies the time discretization itself, not just the spatial operator—it changes how the reaction term gets discretized within the Crank-Nicolson framework by using different P and N matrices, which likely requires a custom scheme class rather than just a new operator.
> 
> I'm reconsidering the implementation strategy: the exponentially fitted operator (Scheme 1) is the most practical starting point since it fits within QuantLib's existing architecture, followed by discrete barrier monitoring, the CN variant, diagnostics for M-matrix and positivity checks, and Rannacher smoothing integration.
> 
> I'll structure this in rounds: first auditing the codebase and mapping the paper's math to QuantLib's conventions, then building the fitted operator and diagnostics as standalone components, followed by the CN variant scheme, discrete barrier projection, and finally Rannacher smoothing modifications.
> 
> But looking more carefully at the paper, Scheme 2 is quite specialized and only works in S-space, which complicates things in QuantLib's log-space framework.
> 
> The exponentially fitted approach is more general and maps naturally to what's already in the template, so I should focus on that while translating the paper's S-space formulations to log-space and implementing the M-matrix analysis and barrier monitoring that are unique to the paper.
> 
> The user wants me to design a framework specifically for this paper's algorithms, not just adapt the existing template, so I need to plan rounds that faithfully implement the paper's methods while leveraging QuantLib's infrastructure—starting with a preflight round to establish the mathematical bridge between S-space and log-space, then building the fitted operator, diagnostics, and test suite around that translation. The issue is that QuantLib's standard Crank-Nicolson applies the same operator L to both the explicit and implicit steps, but this variant needs different operators (N for explicit, P for implicit). I could either create a new scheme class that manages two separate operators, or exploit the fact that when ω₁ = ω₂, the modification becomes symmetric and might be encodable directly into a single operator that produces the combined behavior through the standard split.
> 
> Looking at the math more carefully, the reaction term gets spatially smoothed—it becomes a weighted average across neighboring grid points rather than just the center point. This translates to off-diagonal contributions in the tridiagonal matrix that QuantLib builds during `setTime()`. Since ω is negative, these off-diagonals help preserve the M-matrix property. So both schemes can actually be implemented as modified `FdmLinearOpComposite` classes: Scheme 1 tweaks the diffusion coefficient while using backward Euler, and Scheme 2 modifies the reaction discretization while using standard Crank-Nicolson.
> 
> The implementation roadmap should be: preflight work to translate the paper's formulas to log-space and audit QuantLib conventions, then build the two fitted operators with diagnostics, add barrier monitoring logic, integrate into the scheme and solver layers, wire everything through the engine, and finally validate against the paper's numerical examples. The test suite validates the paper's examples like truncated calls and discrete knockouts, runs convergence tests via Richardson extrapolation, and verifies M-matrix and positivity properties across all three schemes. For the CN variant, I need to reconsider the architecture: since the paper uses different reaction terms ω₂ for the implicit side and ω₁ for the explicit side, but they're equal in the optimal choice, I can implement this as a single modified operator that works with the standard CN scheme rather than requiring separate L_P and L_N operators. I'm going back to the standard CN discretization to clarify the sign conventions: the PDE is -∂V/∂t + LV = 0, which discretizes to (I/Δt - L/2)U^{n+1} = (I/Δt + L/2)U^n, where L includes convection, diffusion, and the reaction term b = -r. For the CN variant, instead of splitting the reaction term evenly as ±r/2 on the diagonal, the discretization spreads it across neighboring nodes using weights ω₁ and ω₂, creating tridiagonal contributions to both P and N operators. For the Crank-Nicolson variant, I need to replace the scalar reaction term with a tridiagonal contribution where the sub- and super-diagonals get -rω while accounting for QuantLib's sign convention where the operator L already includes -r on the diagonal.
> 
> The paper's modification is derived specifically for S-space with uniform grids and the ω parameter depends on that structure, so translating it directly to log-space would require reworking the completing-the-square derivation to account for how uniform log-spacing maps to non-uniform S-spacing.
> 
> For the implementation, I should keep Scheme 1 in log-space where the fitting factor concept works naturally, implement Scheme 2 in its native S-space with a separate mesher, and provide comparison tests to show how both approaches integrate with QuantLib's framework.
> 
> The actual execution will follow a structured progression: first, a preflight round to audit QuantLib conventions and translate the paper's mathematics, then implement the fitted operator with diagnostic tools, followed by the CN variant operator in its appropriate space.
> 
> Next comes barrier monitoring infrastructure with step conditions, and finally modifications to the core scheme and solver components to support operator selection and monitoring-based damping. I'm realizing the user wants me to design the framework structure itself, not just follow an existing template—they want the output formatted according to their template but with rounds I determine based on the paper's requirements. Let me think through the logical progression: starting with a comprehensive code audit in Round 0 to map the QuantLib structure, then building the two schemes incrementally, adding monitoring and validation layers, integrating everything into the engine, and finally testing against the paper's examples. I should aim for 6-7 rounds total to keep it manageable while covering all the mathematical formulations, operator implementations, solver modifications, and test reproduction. I'm realizing the CN variant has a fundamental difference in log-space where coefficients are constant, making the completing-the-square argument inapplicable—the M-matrix condition becomes much simpler than in S-space. Given this complexity, I should prioritize the exponentially fitted operator as the main contribution in Round 1, treat the CN variant as secondary (noting it works best in S-space), and optionally include a direct S-space implementation for validation. Let me now write out the full multi-round framework following the template structure with a shared preamble covering the PDE formulations, fitting factor formulas, M-matrix conditions, and implementation conventions, then detail each round's specific guidance and file requirements. The CN variant modifies the reaction term by spreading it across a three-point stencil with weight ω, which adds ±rω to the off-diagonals and adjusts the diagonal coefficient, fundamentally changing how the scheme handles the reaction component compared to standard CN. Now I'm working out the M-matrix stability condition for the Crank-Nicolson variant, which requires the sub-diagonal element to be non-positive. When the coefficient a/h² - b/(2h) is negative (the problematic case), I need to choose ω sufficiently negative so that |ω| ≥ (b/(2h) - a/h²)/r, and expressing this in terms of the log-space parameters gives a condition that depends on both the grid spacing h and the model parameters.
> 
> The paper's S-space formula ω = -r/(16σ²) comes from a completing-the-square argument that doesn't directly translate to log-space, so the optimal ω in log-space would actually be different—specifically, ω* = σ²/(2rh²) - (r-q-σ²/2)/(2rh), which scales with h rather than being constant. For small volatility and typical grids, this is dominated by the -1/(2h) term, quite different from the S-space constant.
> 
> The CN variant modification effectively adds artificial diffusion to the PDE, and I can implement both schemes as modified spatial operators in QuantLib's framework—the fitted operator changes the diffusion coefficient directly, while the CN variant adds off-diagonal contributions to the reaction term. I'll provide both the paper's S-space formula as a default for ω and the log-space optimal choice as an alternative, then move forward with writing the actual implementation framework. I'm considering the appropriate scope and detail level for the preamble—it needs to be comprehensive enough for accurate code generation while staying self-contained, so I'll aim for around 2000 tokens to cover the paper's mathematical context, key formulas, and QuantLib conventions. Each round will then include task descriptions, source file references with placeholders, and implementation guidance with verification checklists to ensure clarity and correctness across the seven-round structure. I'm deciding what mathematical content belongs in the preamble versus individual rounds—the preamble will cover the core PDE formulations in both coordinate systems, fitting factors, and parameter choices, while round-specific sections handle implementation details and verification procedures. I'm also settling on file and class names for the two schemes and ready to write the complete framework. I'm verifying the naming conventions align with QuantLib's existing patterns, then setting up the mathematical framework by contrasting the paper's formulation in stock space against QuantLib's log-space representation, where I'll define the Péclet parameter and fitting factor for the first scheme's discretization. The operator in QuantLib is assembled as L = b·∂_x + a·∂_xx - r, and for the implicit scheme (I - θ·dt·L) to be an M-matrix, the off-diagonals of L itself must be non-negative. The exponential fitting ensures ℓ ≥ 0 and u ≥ 0 unconditionally, while the standard scheme requires h ≤ 2a/b, and the CN variant imposes an additional constraint on the weighting parameter ω. Now I'm ready to implement this. Let me verify the TripleBandLinearOp indexing: lower_[i] is the coefficient for u_{i-1}, diag_[i] for u_i, and upper_[i] for u_{i+1}, so for the discretized operator L, I need lower_[i] = a/h² - b/(2h) to be non-negative for the M-matrix property to hold. The mMatrixSatisfied() function checks the off-diagonal signs in the system matrix, which is the actual M-matrix diagnostic I need. Now I'm ready to write the complete framework following the template structure, but I need to account for non-uniform meshes since QuantLib uses log-space grids rather than uniform S-space grids like the paper assumes, which affects how the fitting factor and Crank-Nicolson variant scale with grid spacing.

# Multi-Round Prompting Framework: Positivity-Preserving FD Schemes for Black–Scholes

## Implementing Milev & Tagliani (2010) in QuantLib v1.23

This framework implements both schemes from the paper "Low Volatility Options and Numerical Diffusion of Finite Difference Schemes" (Serdica Math. J. 36, 2010, pp. 223–236) within QuantLib v1.23's modern FDM infrastructure. It follows the **corrected** mathematical formulations from the Golden Reference Document (Round 1.5), incorporating all errata resolutions from Rounds 1.1–1.4.

---

### SHARED PREAMBLE

This preamble is included verbatim at the start of every round. It provides the stable mathematical and coding context.

```xml
<system_context>
You are an expert C++ developer implementing two positivity-preserving finite
difference schemes from Milev & Tagliani (2010) into QuantLib v1.23's modern
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
- Verify how FiniteDifferenceModel handles StepCondition::applyTo (every step
  vs stopping times, and whether at the initial 'from' time).

CODING CONVENTIONS (mandatory):
- Namespace: all classes in `namespace QuantLib { }`
- Smart pointers: use `ext::shared_ptr`, NOT `std::shared_ptr`
- Return arrays: use `Array`; if similar methods return `Disposable<Array>`, match
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
- Each file must compile against unmodified QuantLib v1.23 headers plus any
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
  θ_j = b · h_j / (2·a)                   (local Péclet)
  ρ_j = θ_j · coth(θ_j)                   (fitting factor; ρ → 1+θ²/3 for |θ|<1e-8)

Fitted stencil (replacing the standard centered operator):
  ℓ_j = a·ρ_j/h² − b/(2h)     (lower off-diagonal; GUARANTEED ≥ 0)
  d_j = −2·a·ρ_j/h² − r        (main diagonal)
  u_j = a·ρ_j/h² + b/(2h)     (upper off-diagonal; GUARANTEED ≥ 0)

PROOF OF NON-NEGATIVITY (critical for M-matrix):
  ℓ_j = (b/(2h))·(coth(θ) − 1) ≥ 0  for all θ (b>0 ⟹ coth>1; b<0 ⟹ both negative)
  u_j = (b/(2h))·(coth(θ) + 1) ≥ 0  for all θ (same sign argument)

This GUARANTEES the system matrix (I − θ·dt·L) is an M-matrix for any h, dt.

Artificial diffusion (CE-15, low-σ limit): ½·r·S·ΔS·V_SS in S-space,
equivalently ½·|b|·h·u_xx in log-space when |θ| ≫ 1.

NONUNIFORM GRID NOTE: For non-uniform meshes, compute a per-node effective
spacing consistent with QuantLib's derivative operators:
  h_j := 0.5*(dplus_j + dminus_j)

═══ SCHEME 2: CN VARIANT WITH MODIFIED REACTION TERM (CE-16 through CE-23) ═══

The standard −r·u_j is replaced by the 6-node stencil:
  −r·[ω·u_{j−1} + (1−2ω)·u_j + ω·u_{j+1}]

at each time level, with ω₁ = ω₂ = ω (symmetric weighting).

In log-space, this adds off-diagonal reaction contributions to the operator:
  ℓ'_j = a/h² − b/(2h) − rω   (modified lower)
  d'_j = −2a/h² − r(1−2ω)     (modified diagonal)
  u'_j = a/h² + b/(2h) − rω   (modified upper)

PARAMETER CHOICE (CE-19, adapted to log-space):
  S-space optimal: ω = −r/(16σ²)   (from completing-the-square in S_j)
  Log-space: since coefficients are constant, ω must satisfy:
    ℓ'_j ≥ 0  ⟹  ω ≤ (a/h² − b/(2h))/r
    u'_j ≥ 0  ⟹  ω ≤ (a/h² + b/(2h))/r
  The binding constraint (tighter) determines the maximum |ω|.

TIME-STEP CONSTRAINT (CE-20): For Scheme 2 with CN time-stepping, the
diagonal of the explicit-side matrix N must be non-negative:
  Δt < 1 / [r(½ − 2ω) + ½(σM)²]   (S-space form)

Artificial diffusion (CE-23): ⅛·(r·ΔS/σ)²·V_SS in S-space.

═══ DISCRETE BARRIER MONITORING (CE-6) ═══

At each monitoring date t_i, after the time-step solve:
  U_j ← U_j · 𝟙_{[L,U]}(S_j)
This re-introduces discontinuities that the schemes must handle.

═══ M-MATRIX DIAGNOSTIC (CE-10, Proposition 4.1) ═══

For the system matrix (I − θ·dt·L) stored as mapT_:
  M-matrix satisfied ⟺ all off-diagonals of mapT_ are NON-NEGATIVE
  (because system off-diag = −θ·dt·(operator off-diag), and θ·dt > 0)
</mathematical_context>
```

---

### ROUND 0 (OPTIONAL) — Preflight Source Audit + Coordinate Translation

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 0 (OPTIONAL): Preflight mode — do NOT write code.

You are given the complete set of QuantLib v1.23 header files from the
methods/finitedifferences/ directory tree. Perform a systematic audit to
answer the questions below and produce a file-level implementation roadmap.

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

4) TIME-STEPPING:
   a) Does FiniteDifferenceModel call StepCondition::applyTo at every time
      step, or only at stopping times?
   b) Is applyTo called at the initial 'from' time?
   c) Is the evolver copied by value inside the model?

5) SCHEME CONFIGURATION:
   a) Are FdmSchemeDesc members (type, theta, mu) const or non-const?
   b) How does FdmBackwardSolver's rollback method handle damping steps?
   c) How does CrankNicolsonScheme combine explicit and implicit sub-steps?

6) SIGN CONVENTIONS:
   a) For the operator L stored in mapT_, what sign convention do the
      off-diagonals use?
   b) For the system matrix (I − θ·dt·L), what sign must the off-diagonals
      of L have for the system to be an M-matrix?
   c) Verify with the standard FdmBlackScholesOp: are its off-diagonals
      always non-negative? Under what grid conditions?

PART B: Coordinate Translation

Using the verified coordinate convention from Part A:

7) Translate the paper's S-space fitting factor (CE-7):
     ρ_j = (μ_j·ΔS/2)·coth(μ_j·ΔS/(2·σ_d^(j)))
   to the log-space Péclet parameter θ_j and fitting factor in terms of
   QuantLib's operator coefficients a = σ²/2 and b = r−q−σ²/2.

8) Translate the paper's CN variant parameter choice (CE-19):
     ω = −r/(16σ²)
   to log-space and identify whether this choice guarantees ℓ'_j ≥ 0
   for QuantLib's grid spacings.

9) Express the M-matrix condition for both schemes in terms of the
   TripleBandLinearOp coefficients (lower_[], upper_[]) after axpyb().

PART C: File-Level Implementation Roadmap

10) Based on the audit, produce a table listing:
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
[INSERT FULL CONTENT OF fdmblackscholesop.cpp]

FILE: ql/methods/finitedifferences/operators/triplebandlinearop.hpp
[INSERT FULL CONTENT OF triplebandlinearop.hpp]

FILE: ql/methods/finitedifferences/operators/triplebandlinearop.cpp
[INSERT FULL CONTENT OF triplebandlinearop.cpp]

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
[INSERT FULL CONTENT OF fdmblackscholesmesher.cpp]

FILE: ql/methods/finitedifferences/meshers/fdmmesher.hpp
[INSERT FULL CONTENT OF fdmmesher.hpp]

FILE: ql/methods/finitedifferences/meshers/concentrating1dmesher.hpp
[INSERT FULL CONTENT OF concentrating1dmesher.hpp]

FILE: ql/methods/finitedifferences/meshers/uniform1dmesher.hpp
[INSERT FULL CONTENT OF uniform1dmesher.hpp]

FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp
[INSERT FULL CONTENT OF cranknicolsonscheme.hpp]

FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.cpp
[INSERT FULL CONTENT OF cranknicolsonscheme.cpp]

FILE: ql/methods/finitedifferences/schemes/impliciteulerscheme.hpp
[INSERT FULL CONTENT OF impliciteulerscheme.hpp]

FILE: ql/methods/finitedifferences/schemes/expliciteulerscheme.hpp
[INSERT FULL CONTENT OF expliciteulerscheme.hpp]

FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
[INSERT FULL CONTENT OF fdmbackwardsolver.hpp]

FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.cpp
[INSERT FULL CONTENT OF fdmbackwardsolver.cpp]

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

### ROUND 1 — Exponentially Fitted Operator (Scheme 1) + Diagnostics

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 1 OF 7: Create the exponentially fitted spatial operator implementing
Scheme 1 from Milev & Tagliani (2010) in QuantLib's log-space framework, plus
a diagnostic utility for solution quality verification.

These depend only on existing QuantLib v1.23 base classes.

OUTPUT: 4 files total (2 headers + 2 implementations).
</task>

<source_files>
The following existing QuantLib headers are provided as reference.

FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.hpp
[INSERT FULL CONTENT OF fdmblackscholesop.hpp]

FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.cpp
[INSERT FULL CONTENT OF fdmblackscholesop.cpp — CRITICAL: needed to see
how setTime() assembles mapT_ via axpyb()]

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
exponentially fitted spatial discretization, following the corrected equation
chain CE-7 through CE-12 from the Golden Reference Document.

MATHEMATICAL BASIS (log-space adaptation of CE-7):
The paper's S-space fitting factor ρ = (μΔS/2)·coth(μΔS/(2σ_d)) translates
to log-space as follows. In log-space with PDE coefficient a = σ²/2 and
b = (r−q) − σ²/2:
  θ_j = b · h_j / (2·a)          (Péclet parameter)
  ρ_j = θ_j · coth(θ_j)          (fitting factor; ρ ≥ 1 for all θ)
The fitted diffusion coefficient is a·ρ_j (replacing the physical a in the
second-derivative stencil), guaranteeing non-negative off-diagonals.

PRE-CODING FACT CHECKS (mandatory; do not print):
- Confirm from FdmBlackScholesOp::setTime() how the standard operator is
  assembled and mirror the same structure with the fitted modification.
- Confirm how FirstDerivativeOp and SecondDerivativeOp handle nonuniform
  meshes, and derive h_j from the mesher's dplus/dminus.

DESIGN: Same member layout as FdmBlackScholesOp (dxMap_, dxxMap_, mapT_),
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
      //    a) Compute local effective spacing:
      //       h_i = 0.5*(mesher_->dplus(iter, direction_)
      //                 + mesher_->dminus(iter, direction_))
      //       (consistent with derivative operators on nonuniform meshes)
      //    b) Compute a_i = 0.5*vol*vol  (or from local vol if applicable)
      //    c) Compute b_i = (r - q) - a_i  (drift in log-space)
      //    d) Compute theta_i = b_i * h_i / (2.0 * a_i)
      //       Guard: clamp a_i to max(a_i, 1e-20) before dividing
      //    e) Compute rho_i = fittingFactor(theta_i)
      //    f) Store: fittedDiffusion[i] = a_i * rho_i
      //              convection[i] = b_i
      // 3. Assemble via axpyb:
      //    mapT_.axpyb(convection, dxMap_,
      //                dxxMap_.mult(fittedDiffusion), Array(1, -r));
  }

Static fitting factor (CE-7 numerical guard):
  static Real fittingFactor(Real theta) {
      if (std::fabs(theta) < 1e-8)
          return 1.0 + theta*theta/3.0;    // Taylor: θ·coth(θ) ≈ 1 + θ²/3
      return theta / std::tanh(theta);      // = θ·coth(θ)
  }

  // Note: for |θ| > 20, coth(θ) ≈ sign(θ), so ρ ≈ |θ|.
  // std::tanh handles this correctly (no overflow), so no extra guard needed.

M-MATRIX DIAGNOSTIC METHODS:

  bool mMatrixSatisfied() const;
  // Returns true iff ALL off-diagonals of mapT_ are ≥ 0.
  // VERIFY the sign convention from the provided TripleBandLinearOp source
  // before implementing: which arrays (lower_, upper_) store which coefficients.

  Size mMatrixViolationCount() const;
  // Returns the number of nodes where any off-diagonal is negative.
  // For the fitted operator, this should ALWAYS be 0 (mathematical guarantee).

IMPORTANT: The exponential fitting guarantees non-negative off-diagonals
for ALL grid spacings h > 0 and ALL parameter values (CE-10, verified in
Round 1.4 §3a). If mMatrixViolationCount() > 0, it indicates a bug.

Pure virtual methods from FdmLinearOpComposite: implement ALL of them,
mirroring FdmBlackScholesOp's implementations for size(), setTime(),
apply(), apply_mixed(), apply_direction(), solve_splitting(), preconditioner().


=== FILE PAIR 2: FdmDiagnostics ===

Path: ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp
      ql/methods/finitedifferences/utilities/fdmdiagnostics.cpp

Lightweight utility for verifying solution quality, motivated by the paper's
analysis of positivity preservation (CE-12) and oscillation freedom (CE-10).

  struct FdmDiagnosticsReport {
      Real minValue;                // min over all nodes (negative → positivity failure)
      Size negativeCount;           // number of nodes with u[i] < 0
      Real oscillationScore;        // sign-change frequency in Δu
      Size mMatrixViolationCount;   // from operator diagnostic (0 for fitted)
      Size nanCount;                // number of NaN/Inf values
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
      // Takes worst-case of each field (min of minValue, max of counts/scores)

      Level level() const;
    private:
      Level level_;
  };

</implementation_guidance>

<constraints>
- Do NOT modify any existing QuantLib files in this round
- The fitted operator must use axpyb() to populate mapT_ — do NOT attempt
  to access protected TripleBandLinearOp members
- For the fitting factor: use theta/tanh(theta), NOT theta*coth(theta), to
  avoid implementing coth manually (tanh is in <cmath>)
- The M-matrix diagnostic must verify the sign convention from the actual
  TripleBandLinearOp source before checking signs
- Use ext::shared_ptr and ext::make_shared throughout
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
□ fittingFactor returns 1.0 + θ²/3 for small θ, θ/tanh(θ) otherwise
□ Fitting factor uses a_i = max(0.5*σ², 1e-20) to prevent division by zero
□ Per-node h_i is computed consistently with derivative operators (dplus+dminus)/2
□ mMatrixSatisfied() checks the correct sign (verified from TripleBandLinearOp source)
□ mMatrixViolationCount() should return 0 for the fitted operator (mathematical guarantee)
□ Diagnostics oscillationScore is O(N) and allocation-free
□ All ext::shared_ptr, not std::shared_ptr
□ Include guards match QuantLib naming convention
□ If any needed source file content is missing, STOP and ask (do not reconstruct)
</quality_checklist>
```

---

### ROUND 2 — CN Variant Operator (Scheme 2) + Barrier Monitoring

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 2 OF 7: Create the Crank–Nicolson variant spatial operator implementing
Scheme 2 from Milev & Tagliani (2010) in QuantLib's log-space framework, and
create the discrete barrier projection step condition.

These depend only on existing QuantLib v1.23 base classes.

OUTPUT: 4 files total (2 headers + 2 implementations).
</task>

<source_files>
FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/operators/fdmblackscholesop.cpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/operators/triplebandlinearop.hpp
[INSERT FULL CONTENT]

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
[INSERT ROUND 1 OUTPUT HEADER — for fittingFactor() reuse if desired]
</source_files>

<implementation_guidance>

=== FILE PAIR 1: FdmCNVariantBlackScholesOp (Scheme 2 in log-space) ===

Path: ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.hpp
      ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.cpp

This class implements the Crank–Nicolson variant from §3 of Milev & Tagliani
(2010), adapted to QuantLib's log-space framework. The modification replaces
the standard reaction term −r·u_j with the 6-node stencil:
  −r·[ω·u_{j−1} + (1−2ω)·u_j + ω·u_{j+1}]
where ω < 0 is chosen to ensure the operator's off-diagonals are non-negative
(M-matrix property) under milder grid restrictions than standard centered.

MATHEMATICAL BASIS (CE-16 through CE-21, adapted to log-space):

In log-space, the modified operator stencil at node j is:
  ℓ'_j = a/h² − b/(2h) − r·ω     (lower off-diagonal)
  d'_j = −2a/h² − r·(1−2ω)        (main diagonal)
  u'_j = a/h² + b/(2h) − r·ω      (upper off-diagonal)

For M-matrix: ℓ'_j ≥ 0 and u'_j ≥ 0. Since ω < 0, the term −r·ω > 0,
which ADDS to the off-diagonals (helping satisfy the M-matrix condition).

PARAMETER CHOICE (two modes):
Mode A (paper's S-space formula, CE-19):
  ω_paper = −r / (16·σ²)
  This is the optimal S-space value; in log-space it provides EXTRA margin
  beyond what's needed, at the cost of more artificial diffusion.

Mode B (log-space optimal):
  ω_log = −max(0, b/(2h) − a/h²) / r    (per-node, just enough to fix sign)
  This minimizes artificial diffusion but varies across nodes.

DESIGN DECISION: Use Mode A as the default (consistent with the paper) and
store ω as a constructor parameter with default value. The user can override.

Constructor signature:
  FdmCNVariantBlackScholesOp(
      ext::shared_ptr<FdmMesher> mesher,
      ext::shared_ptr<GeneralizedBlackScholesProcess> process,
      Real strike,
      Real omega = Null<Real>(),     // Null → auto-compute from paper's formula
      bool localVol = false,
      Real illegalLocalVolOverwrite = -Null<Real>(),
      Size direction = 0,
      ext::shared_ptr<FdmQuantoHelper> quantoHelper
          = ext::shared_ptr<FdmQuantoHelper>());

The setTime() method:
  void FdmCNVariantBlackScholesOp::setTime(Time t1, Time t2) {
      // 1. Get r, q, vol exactly as FdmBlackScholesOp does
      // 2. Compute omega:
      //    If omega_ was not provided (Null): omega = -r/(16*vol*vol)
      //    Else: use the user-provided value
      // 3. For each layout point i:
      //    a) Compute a_i = 0.5*vol*vol, b_i = (r-q) - a_i
      //    b) convection[i] = b_i
      //       modifiedDiffusion[i] = a_i
      // 4. Assemble the STANDARD convection+diffusion part:
      //    mapT_.axpyb(convection, dxMap_,
      //                dxxMap_.mult(modifiedDiffusion),
      //                Array(1, -r*(1.0 - 2.0*omega)));
      //    NOTE: The diagonal reaction term uses −r(1−2ω) instead of −r.
      //
      // 5. Add the off-diagonal reaction contributions:
      //    This requires accessing mapT_ to add −r·ω to the off-diagonals.
      //    APPROACH: Use axpyb() to build the base operator, then construct
      //    a separate TripleBandLinearOp for the reaction off-diagonal and
      //    combine using add().
      //
      //    reactionOffDiag = TripleBandLinearOp(direction, mesher_->layout())
      //    For each node: lower = -r*omega, diag = 0, upper = -r*omega
      //    mapT_ = mapT_ + reactionOffDiag    (using TripleBandLinearOp::add)

  ALTERNATIVE ASSEMBLY (if TripleBandLinearOp::add is not available):
  Build the entire operator manually using a dedicated TripleBandLinearOp
  without relying on axpyb(). Set lower/diag/upper node-by-node using the
  full stencil formulas ℓ'_j, d'_j, u'_j. Verify whether TripleBandLinearOp
  provides node-level write access from the provided source.

Additional public methods:
  Real omega() const;                    // returns the active ω value
  Real timestepConstraint() const;       // returns max safe Δt from CE-20
  bool mMatrixSatisfied() const;         // checks off-diagonal signs
  Size mMatrixViolationCount() const;    // count of violations

  // timestepConstraint() implements CE-20 adapted to log-space:
  //   For each node j: bound_j = 1.0 / [r*(1-2ω) + (σ·S_j/ΔS)²]
  //   In log-space: bound = 1.0 / [r*(1/2 - 2ω) + σ²·M²/2]
  //   where M is related to the grid extent.
  //   Return the minimum bound across all nodes.


=== FILE PAIR 2: FdmBarrierProjectionCondition ===

Path: ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp
      ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.cpp

Implements the discrete barrier monitoring from CE-6:
  U_j ← U_j · 𝟙_{[L,U]}(S_j)
at specified monitoring times, zeroing out grid values outside the barrier corridor.

Constructor:
  FdmBarrierProjectionCondition(
      std::vector<Time> monitoringTimes,
      Real lowerBarrier,          // in S-space (NOT log-space); 0 means no lower
      Real upperBarrier,          // in S-space; Null<Real>() means no upper
      ext::shared_ptr<FdmMesher> mesher,
      Size direction = 0);

COORDINATE CONVERSION (CRITICAL):
- The constructor must convert S-space barriers to the mesher's coordinate:
  * VERIFY from provided sources what mesher->location(iter, direction) returns.
  * If x = ln(S): lnLower = log(lowerBarrier), lnUpper = log(upperBarrier)
  * If x = ln(S/K): need strike info — if not available, STOP and ask.
- Precompute outsideIndices_ in constructor: vector of layout indices where
  the spatial coordinate falls outside [lnLower, lnUpper].

  void applyTo(Array& a, Time t) const override;
  // If t matches any monitoring time (within tolerance 1e-10), set a[i]=0
  // for all i in outsideIndices_. Otherwise no-op.

  const std::vector<Time>& monitoringTimes() const;
  // Accessor for registering as stopping times in FdmStepConditionComposite

Edge cases:
- Empty monitoringTimes: constructor succeeds; applyTo always no-op
- lowerBarrier = 0: no lower barrier (only upper knockout)
- upperBarrier = Null<Real>(): no upper barrier (only lower knockout)

</implementation_guidance>

<constraints>
- FdmCNVariantBlackScholesOp must implement ALL pure virtual methods of
  FdmLinearOpComposite
- The off-diagonal reaction contributions must be added without accessing
  protected TripleBandLinearOp members — use only public API (axpyb, add,
  mult, or node-level setters if available)
- If TripleBandLinearOp does not expose an add() method or node-level write
  access in the provided source, document the limitation and propose an
  alternative assembly strategy; STOP and ask if no clean approach exists
- FdmBarrierProjectionCondition must NOT depend on knowing the strike unless
  the mesher coordinate requires it — verify from source first
- outsideIndices_ must be computed once in constructor, not per applyTo call
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
□ FdmCNVariantBlackScholesOp assembles the modified reaction term correctly
□ Default ω = −r/(16σ²) matches CE-19
□ timestepConstraint() returns a positive Real consistent with CE-20
□ mMatrixSatisfied() checks correct sign convention for off-diagonals
□ FdmBarrierProjectionCondition correctly converts S-space barriers to
  the mesher's coordinate system (VERIFIED from source, not assumed)
□ applyTo is a no-op when t does not match any monitoring time
□ outsideIndices_ computed once in constructor
□ All edge cases handled (empty times, single barrier, etc.)
□ If coordinate convention cannot be verified, implementation STOPs and asks
</quality_checklist>
```

---

### ROUND 3 — Scheme Modifications (Rannacher Restart + FdmSchemeDesc)

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 3 OF 7: Modify CrankNicolsonScheme to support Rannacher-style damping
restart after discrete monitoring events (where CE-6 re-introduces discontinuities),
and modify FdmSchemeDesc to carry the new configuration parameters including
scheme selection for Scheme 1 vs Scheme 2.

These changes are internal to the scheme layer and do not depend on Rounds 1-2.

OUTPUT: 3 files (1 modified FdmSchemeDesc/FdmBackwardSolver header,
1 modified CN scheme header, 1 modified CN scheme implementation).
NOTE: Only FdmSchemeDesc changes inside fdmbackwardsolver.hpp —
the FdmBackwardSolver class body is NOT changed until Round 4.
</task>

<source_files>
FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
[INSERT FULL CONTENT — contains both FdmSchemeDesc struct AND FdmBackwardSolver]

FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/schemes/cranknicolsonscheme.cpp
[INSERT FULL CONTENT — if unavailable STOP and ask; do not reconstruct]

FILE: ql/methods/finitedifferences/schemes/impliciteulerscheme.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/schemes/expliciteulerscheme.hpp
[INSERT FULL CONTENT]
</source_files>

<implementation_guidance>

=== FILE 1: Modified FdmSchemeDesc (in fdmbackwardsolver.hpp) ===

FdmSchemeDesc gains two new fields:
  Size monitoringDampingSteps;    // Rannacher half-steps after each monitoring event
  Size operatorType;              // 0 = standard, 1 = fitted (Scheme 1), 2 = CN variant (Scheme 2)

APPROACH: Check whether existing members are const or non-const. Add fields
accordingly (see template Round 3 for the two approaches).

New static factories:
  static FdmSchemeDesc FittedImplicit();
  // Returns {ImplicitEulerType, 1.0, 0.0, 0, 1}
  // Scheme 1 with fully implicit time-stepping (paper's recommendation)

  static FdmSchemeDesc FittedCrankNicolson();
  // Returns {CrankNicolsonType, 0.5, 0.0, 0, 1}
  // Scheme 1 with CN time-stepping (not in paper, but useful for 2nd-order accuracy)

  static FdmSchemeDesc CNVariant(Size monitoringDampingSteps = 0);
  // Returns {CrankNicolsonType, 0.5, 0.0, monitoringDampingSteps, 2}
  // Scheme 2 with CN time-stepping (paper's formulation)

  static FdmSchemeDesc CrankNicolsonWithDamping(Size monitoringDampingSteps = 2);
  // Returns {CrankNicolsonType, 0.5, 0.0, monitoringDampingSteps, 0}
  // Standard operator with monitoring restart damping

All EXISTING static factories must still compile and produce identical behavior
(operatorType=0, monitoringDampingSteps=0).

CRITICAL: Do NOT change FdmBackwardSolver class body in this round.


=== FILES 2-3: Modified CrankNicolsonScheme ===

Add Rannacher monitoring-restart capability. After a discrete monitoring event
re-introduces a discontinuity, perform implicit Euler damping half-steps before
resuming CN to suppress the spurious oscillations documented in CE-25 / C25.

Modified constructor:
  CrankNicolsonScheme(
      Real theta,
      const ext::shared_ptr<FdmLinearOpComposite>& map,
      const bc_set& bcSet = bc_set(),
      Real relTol = 1e-8,
      ImplicitEulerScheme::SolverType solverType
          = ImplicitEulerScheme::BiCGstab,
      Size dampingHalfSteps = 0);    // NEW — must be even; 0 = no restart damping

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
</constraints>

<output_specification>
Produce exactly 3 files:
1. fdmbackwardsolver.hpp (FdmSchemeDesc modified, FdmBackwardSolver unchanged)
2. cranknicolsonscheme.hpp (modified)
3. cranknicolsonscheme.cpp (modified — complete implementation)
</output_specification>

<quality_checklist>
□ All existing FdmSchemeDesc factories return operatorType=0, monitoringDampingSteps=0
□ New factories produce correct operatorType values
□ CrankNicolsonScheme with dampingHalfSteps=0 behaves identically to original
□ QL_REQUIRE enforces dampingHalfSteps is even
□ step() during damping performs exactly 2 implicit Euler half-steps per call
□ implicit_->setStep restored to dt_ after damping half-steps
□ FdmBackwardSolver class declaration unchanged from input
</quality_checklist>
```

---

### ROUND 4 — Solver Modifications (Backward Solver + BS Solver)

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 4 OF 7: Modify FdmBackwardSolver to support monitoring-restart damping
and operator type selection. Modify FdmBlackScholesSolver to create the
appropriate operator (standard, fitted, or CN variant) based on the scheme
description.

This round depends on:
- Round 1: FdmFittedBlackScholesOp
- Round 2: FdmCNVariantBlackScholesOp
- Round 3: Modified CrankNicolsonScheme + FdmSchemeDesc

OUTPUT: 4 files (2 modified headers + 2 modified implementations).
</task>

<source_files>
FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
[INSERT ROUND 3 OUTPUT VERSION — with modified FdmSchemeDesc]

FILE: ql/methods/finitedifferences/solvers/fdmbackwardsolver.cpp
[INSERT FULL CONTENT if available; if unavailable STOP and ask]

FILE: ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp
[INSERT FULL CONTENT if available; if unavailable STOP and ask]

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

=== FILE 1: Modified FdmBackwardSolver ===

FACT-CHECKS (mandatory; do not print):
- Read finitedifferencemodel.hpp to confirm: evolver copied by value? applyTo
  called when? Time ordering conventions?

The rollback method gains monitoring-restart damping for the CN case:

When schemeDesc_.monitoringDampingSteps > 0 AND schemeDesc_.type == CrankNicolsonType:
1. Phase 1 (unchanged): if dampingSteps != 0, run original implicit damping.
2. Phase 2 (CN): segment the CN phase at stopping times.
   At each segment boundary that is a monitoring stopping time:
   - Call cnEvolver.notifyDiscontinuity() to trigger Rannacher restart
   - Continue rollback; the scheme's internal step() handles the damping

When monitoringDampingSteps == 0: behavior IDENTICAL to original.
When type != CrankNicolsonType: behavior IDENTICAL to original.

Segmentation approach:
- Collect stopping times from condition_ that lie within the CN phase interval
- Split the CN rollback into segments at these times
- For each segment, create FiniteDifferenceModel<CrankNicolsonScheme> and rollback
- After each segment that ends at a stopping time, the step condition's applyTo
  will have been called (by the model); then notify the next segment's evolver.

See the template's Round 4 for detailed segment distribution logic.
The key difference: use CrankNicolsonScheme::notifyDiscontinuity() from Round 3.


=== FILE 2: Modified FdmBlackScholesSolver ===

Add operator selection based on schemeDesc_.operatorType:

In performCalculations():
  ext::shared_ptr<FdmLinearOpComposite> op;
  switch (schemeDesc_.operatorType) {
      case 0:  // standard
          op = ext::make_shared<FdmBlackScholesOp>(
              solverDesc_.mesher, process_.currentLink(), strike_,
              localVol_, illegalLocalVolOverwrite_, 0, quantoHelper);
          break;
      case 1:  // Scheme 1: exponentially fitted
          op = ext::make_shared<FdmFittedBlackScholesOp>(
              solverDesc_.mesher, process_.currentLink(), strike_,
              localVol_, illegalLocalVolOverwrite_, 0, quantoHelper);
          break;
      case 2:  // Scheme 2: CN variant
          op = ext::make_shared<FdmCNVariantBlackScholesOp>(
              solverDesc_.mesher, process_.currentLink(), strike_,
              Null<Real>(),  // auto-compute ω
              localVol_, illegalLocalVolOverwrite_, 0, quantoHelper);
          break;
      default:
          QL_FAIL("Unknown operator type " << schemeDesc_.operatorType);
  }
  solver_ = ext::make_shared<Fdm1DimSolver>(solverDesc_, schemeDesc_, op);

Include the new headers in the .cpp (NOT in .hpp):
  #include <ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp>
  #include <ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.hpp>

</implementation_guidance>

<constraints>
- When monitoringDampingSteps==0: rollback IDENTICAL to original
- When operatorType==0: solver IDENTICAL to original
- Do NOT modify FiniteDifferenceModel
- New operator headers #included in .cpp only (not .hpp)
- If any .cpp source file content is missing, STOP and ask
</constraints>

<output_specification>
Produce exactly 4 files:
1. fdmbackwardsolver.hpp (from Round 3, possibly with FdmBackwardSolver changes)
2. fdmbackwardsolver.cpp (modified rollback)
3. fdmblackscholessolver.hpp (modified — no new public API needed if operatorType
   is in FdmSchemeDesc; if constructor changes are needed, add with defaults)
4. fdmblackscholessolver.cpp (modified — conditional operator creation)
</output_specification>

<quality_checklist>
□ Original rollback behavior preserved when monitoringDampingSteps==0
□ Original solver behavior preserved when operatorType==0
□ CN phase segmented at stopping times; notifyDiscontinuity used for restart
□ Step counts distributed proportionally across segments
□ Other scheme types (Douglas, Hundsdorfer) unchanged
□ New headers #included in .cpp only
□ If any required .cpp source is missing, STOP and ask
</quality_checklist>
```

---

### ROUND 5 — Engine Modifications + Wiring

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 5 OF 7: Modify the user-facing pricing engine to expose the new
capabilities, wire the barrier projection condition into the step condition
composite, and update header registrations.

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
[INSERT FULL CONTENT if available; if unavailable STOP and ask]

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

=== FILE 1: Modified FdBlackScholesVanillaEngine ===

The engine constructor gains backward-compatible parameters:

  FdBlackScholesVanillaEngine(
      ext::shared_ptr<GeneralizedBlackScholesProcess> process,
      Size tGrid = 100, Size xGrid = 100,
      Size dampingSteps = 0,
      const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
      bool localVol = false,
      Real illegalLocalVolOverwrite = -Null<Real>());
      // EXISTING signature — must still compile

  // NEW overload or additional default parameters:
  // The operatorType and monitoringDampingSteps are carried inside FdmSchemeDesc.
  // Users select the scheme via:
  //   FdmSchemeDesc::FittedImplicit()          — Scheme 1 + backward Euler
  //   FdmSchemeDesc::FittedCrankNicolson()     — Scheme 1 + CN
  //   FdmSchemeDesc::CNVariant()               — Scheme 2 + CN
  //   FdmSchemeDesc::CrankNicolsonWithDamping(2)  — standard + monitoring restart

USER EXAMPLE (as comment in header):
  // Price with Scheme 1 (exponentially fitted + fully implicit):
  auto engine = ext::make_shared<FdBlackScholesVanillaEngine>(
      process, 800, 800, 2, FdmSchemeDesc::FittedImplicit());
  //
  // Price with Scheme 2 (CN variant with monitoring restart):
  auto engine = ext::make_shared<FdBlackScholesVanillaEngine>(
      process, 800, 800, 2, FdmSchemeDesc::CNVariant(2));


=== FILE 2: Modified FdmStepConditionComposite ===

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

This factory:
1. Calls vanillaComposite() for the base condition set
2. Converts monitoringDates to Times
3. Creates FdmBarrierProjectionCondition with barriers and times
4. Adds it to the condition list; adds monitoring times to stopping times
5. Returns the assembled composite

Empty monitoringDates must return same result as vanillaComposite().


=== FILES 3-5: Updated all.hpp headers ===

operators/all.hpp: add
  #include <ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp>
  #include <ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.hpp>

stepconditions/all.hpp: add
  #include <ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp>

utilities/all.hpp: add
  #include <ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp>

All additions in alphabetical order within existing includes.

</implementation_guidance>

<constraints>
- Existing constructor signature must still compile
- FdmSchemeDesc::Douglas() default produces identical results to original
- all.hpp additions alphabetically sorted
- No circular include dependencies
- If engine file path or content missing, STOP and ask
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
□ all.hpp entries alphabetically sorted
□ No new public dependencies in engine .hpp
</quality_checklist>
```

---

### ROUND 6 — Test Suite: Paper Reproduction + Validation

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 6 OF 7: Create a comprehensive test file that validates the
implementations against the paper's numerical examples and theoretical
guarantees. Tests serve as both validation and regression protection.

OUTPUT: 1 file (test implementation).

FACT-CHECK (mandatory; do not print):
- Inspect any provided test file(s) to confirm the framework used (Boost.Test,
  custom macros, etc.) and naming/registration conventions.
- If no test file is provided, STOP and ask for one.
</task>

<source_files>
All output headers from Rounds 1-5:
[INSERT all .hpp files from Rounds 1-5]

QuantLib test infrastructure:
[INSERT at least one existing test file showing registration pattern,
e.g., test-suite/europeanoption.cpp or similar]
</source_files>

<implementation_guidance>

File: test-suite/fdmpositivitypreserving.cpp
(adapt naming to match QuantLib's test conventions)

=== TEST T1: Fitted Operator M-Matrix Guarantee ===
Purpose: Verify the mathematical guarantee from CE-10 that the fitted
operator always produces non-negative off-diagonals.

Setup: Create FdmFittedBlackScholesOp with a range of parameters:
  - σ = {0.001, 0.01, 0.1, 0.5, 1.0}
  - r = {0.01, 0.05, 0.1, 0.5}
  - Grid sizes: 50, 200, 800
  - Both uniform and concentrated meshes
For each combination: call setTime(0, T), then check mMatrixSatisfied().

Assertion: mMatrixViolationCount() == 0 for ALL parameter combinations.
(This is a mathematical guarantee, not a numerical tolerance check.)


=== TEST T2: Positivity Preservation (Scheme 1 vs Standard CN) ===
Purpose: Reproduce the paper's central finding (C25, C26) that standard CN
produces negative prices for low-σ with discontinuous payoffs, while Scheme 1
preserves positivity.

Setup: Truncated call (paper's Definition 4.1):
  K=50, U=70 (upper cutoff), T=5/12, r=0.05, σ=0.001, S_max=140
  Grid: ΔS=0.05 equivalent in log-space (~2800 nodes), Δt=0.01

Run with:
  (a) Standard FdmBlackScholesOp + CrankNicolson
  (b) FdmFittedBlackScholesOp + ImplicitEuler (Scheme 1)

Assertions:
  (a) Standard CN: FdmDiagnostics::checkSolution reports negativeCount > 0
      OR oscillationScore > 0.05
  (b) Fitted: min(u) ≥ 0, negativeCount == 0, oscillationScore < 0.01


=== TEST T3: Positivity Preservation (Scheme 2 — CN Variant) ===
Purpose: Verify the CN variant (CE-17/CE-18) produces positive, non-oscillating
solutions for the same problematic setup as T2.

Setup: Same truncated call as T2.
Run with FdmCNVariantBlackScholesOp + CrankNicolson (Scheme 2).

Assertions: min(u) ≥ 0 (or ≥ −tolerance if time-step constraint CE-20 is
violated), oscillationScore < 0.01.

Additional check: verify timestepConstraint() < Δt used (documenting that
the paper's examples violate CE-20, per Round 1.3 §G7.3).


=== TEST T4: Convergence Rate (Scheme 1) ===
Purpose: Verify the convergence rate from CE-13: O(h+k) uniform in σ, and
O(h²+k) for fixed σ > 0.

Setup: European call (smooth payoff, to isolate spatial convergence):
  S₀=100, K=100, T=1.0, r=0.05, q=0.02, σ=0.2
  Reference: Black-Scholes closed form

Run Scheme 1 (fitted + implicit Euler) at grid sizes N = 50, 100, 200, 400.
Compute errors at S₀.

Assertion: Richardson error ratio between successive refinements should
approach 4.0 (for O(h²) spatial) or 2.0 (for O(h) when σ is very small).


=== TEST T5: Convergence Rate (Scheme 2) ===
Purpose: Verify O(h²+k²) from CE-23 for the CN variant.

Setup: Same European call as T4.
Run Scheme 2 (CN variant + CN time-stepping) at multiple grid sizes.

Assertion: Richardson error ratio should approach 4.0 (second-order spatial).


=== TEST T6: Discrete Double Barrier Knock-Out ===
Purpose: Reproduce the paper's Example 4.2 (discrete DKO call) and verify
that barrier projection + monitoring restart damping produces non-oscillating,
positive prices.

Setup: K=100, L=95, U=110, T=1, r=0.05, σ=0.001
  Monthly monitoring (12 dates)
  Use FdmSchemeDesc::FittedImplicit() for Scheme 1
  Use FdmSchemeDesc::CNVariant(2) for Scheme 2

Assertions:
  Both: min(u) ≥ 0, oscillationScore < 0.01
  Values should be "small" (the option is frequently knocked out with
  tight barriers and low vol)

Also test: Standard CN without damping should show oscillations near barriers
(oscillationScore > threshold).


=== TEST T7: Artificial Diffusion Comparison ===
Purpose: Verify the paper's prediction (CE-25) that Scheme 1 has artificial
diffusion O(ΔS) while Scheme 2 has O(ΔS²).

Setup: Truncated call as in T2, but run at multiple ΔS values.
Compute max |u_fitted − u_reference| and max |u_cnvariant − u_reference|
where the reference is a very fine grid solution.

Assertion: The error ratio (coarse/fine) for Scheme 1 should approach 2
(first-order) and for Scheme 2 should approach 4 (second-order), at least
when σ is not too small.


=== TEST T8: Backward Compatibility ===
Purpose: Verify that existing code using default FdmSchemeDesc produces
bit-identical results before and after our modifications.

Setup: Standard European call with FdmSchemeDesc::Douglas() and
FdmSchemeDesc::CrankNicolson(). Compare against stored reference values.

Assertion: NPV matches to machine epsilon.

</implementation_guidance>

<constraints>
- All tests must pass with the improved-CN configuration
- Tests must not depend on external data files
- Each test should complete in < 10 seconds on modern hardware
- T2/T3 tolerances: use 0.0 for strict positivity checks on Scheme 1;
  use small negative tolerance (-1e-10) for Scheme 2 if CE-20 is violated
- T8 uses stored reference values — compute them once from unmodified QuantLib
  and hardcode in the test
- Use appropriate test framework macros (Boost.Test or QuantLib custom)
</constraints>

<output_specification>
Produce 1 file:
1. fdmpositivitypreserving.cpp (complete test file)

Include a comment block at the top listing:
- Test-to-paper-result mapping
- Grid sizes used
- Expected run time
</output_specification>

<quality_checklist>
□ All 8 tests present with meaningful assertions
□ T1 checks M-matrix guarantee across many parameter combinations
□ T2 demonstrates the problem (standard CN fails) and the fix (Scheme 1 works)
□ T3 validates Scheme 2 and documents CE-20 constraint violation
□ T4/T5 check convergence rates via Richardson ratios
□ T6 uses barrier monitoring with both schemes
□ T7 quantifies artificial diffusion difference between schemes
□ T8 ensures backward compatibility
□ All tests use QuantLib infrastructure (Process, TermStructure, etc.)
□ If test framework patterns are missing, STOP and ask
</quality_checklist>
```

---

### ROUND 7 (OPTIONAL) — Direct S-Space Reference Implementation

```xml
[INSERT SHARED PREAMBLE HERE]

<task>
ROUND 7 (OPTIONAL): Create a standalone S-space implementation of both schemes
exactly as described in the paper, for validation against the log-space
implementations from Rounds 1-2.

This implementation uses the paper's native S-space coordinates with uniform
grids and follows the pseudocode from the Golden Reference Document (§5)
verbatim. It does NOT integrate into QuantLib's FDM framework — it is a
self-contained validation tool.

OUTPUT: 2 files (1 header + 1 implementation).
</task>

<source_files>
ROUND 1.5 OUTPUT (Golden Reference Document):
[INSERT the complete Golden Reference Document from round_1.5_output.md,
specifically §3 (Mathematical Specification) and §5 (Monolithic Pseudocode)]

ROUND 1.4 OUTPUT (Corrected Equation Chain):
[INSERT the Corrected Equation Chain CE-1 through CE-25 from round_1.4_output.md]

QuantLib headers (minimal):
FILE: ql/math/array.hpp
[INSERT class declaration]
</source_files>

<implementation_guidance>

=== FILES: FdmSSpaceReferenceSolver ===

Path: ql/methods/finitedifferences/utilities/fdmsspacereferencesolver.hpp
      ql/methods/finitedifferences/utilities/fdmsspacereferencesolver.cpp

A standalone solver that implements BOTH schemes in S-space exactly as
specified in the Golden Reference Document.

  class FdmSSpaceReferenceSolver {
    public:
      enum Scheme { ExponentiallyFitted, CrankNicolsonVariant };

      struct Result {
          Array prices;           // V[0..M] at t=T (present)
          Array grid;             // S[0..M]
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

      // Truncated call (Definition 4.1)
      Result solveTruncatedCall(Real upperCutoff) const;

      // Discrete DKO call (Definition 4.2)
      Result solveDiscreteBarrierCall(
          Real lowerBarrier, Real upperBarrier,
          const std::vector<Time>& monitoringTimes) const;

    private:
      // Direct Thomas algorithm (CE-24/§3(g))
      static Array thomasSolve(
          const Array& sub, const Array& diag,
          const Array& sup, const Array& rhs);

      // Fitting factor computation (CE-7)
      static Real computeRho(Real mu_j, Real sigma_d_j, Real deltaS);

      Real r_, sigma_, K_, T_, Smax_, deltaS_, deltaT_;
      Scheme scheme_;
  };

IMPLEMENTATION NOTES:
- Follow the pseudocode from §5 of the Golden Reference Document line-by-line
- For Scheme 1 (ExponentiallyFitted): assemble matrix A per CE-9, solve
  A·U^{n+1} = U^n at each time step using Thomas algorithm
- For Scheme 2 (CrankNicolsonVariant): assemble P and N per CE-17/CE-18
  with ω = -r/(16σ²), solve P·U^{n+1} = N·U^n
- Use the numerical guards for coth from CE-7 (§3(e) of Golden Reference)
- Apply monitoring projection per CE-6 at monitoring dates
- Track positivity (min price) and oscillations throughout

This solver is used in test T9 (add to Round 6 test file):

TEST T9: Cross-validate log-space vs S-space implementations.
Run both S-space and log-space implementations for the same parameters.
For smooth payoffs (European call), they should agree to O(h²) accuracy.
For discontinuous payoffs (truncated call), both should show the same
qualitative behavior (positive, non-oscillating for fitted scheme).

</implementation_guidance>

<constraints>
- This solver does NOT use QuantLib's FDM framework (no FdmMesher, no
  FdmLinearOpComposite, no schemes)
- It uses only ql/math/array.hpp for Array storage
- It must handle the edge cases from the Golden Reference (CE-20 constraint
  violation, boundary conditions, monitoring dates on grid)
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
□ Positivity tracking throughout time-stepping
□ Handles CE-20 constraint violation gracefully (warns, does not crash)
□ Row sum verification: A has row sum 1+rΔt (CE-9); P has 1/Δt+r/2, N has 1/Δt-r/2
</quality_checklist>
```

---

## Usage Instructions

**Prerequisite: Gather QuantLib v1.23 sources.** Before starting, collect the actual file contents for every `[INSERT ...]` placeholder. The files needed per round are listed in each round's `<source_files>` section. Do NOT provide unnecessary files — this wastes context.

**Execution order:**

| Step | Round | Depends On | Key Output |
|:---|:---|:---|:---|
| 1 | Round 0 (optional) | None | Audit report + roadmap |
| 2 | Round 1 | None | FdmFittedBlackScholesOp + FdmDiagnostics |
| 3 | Round 2 | None | FdmCNVariantBlackScholesOp + FdmBarrierProjectionCondition |
| 4 | Round 3 | None | Modified CrankNicolsonScheme + FdmSchemeDesc |
| 5 | Round 4 | R1, R2, R3 | Modified FdmBackwardSolver + FdmBlackScholesSolver |
| 6 | Round 5 | R1-R4 | Modified engine + wiring + all.hpp |
| 7 | Round 6 | R1-R5 | Test suite |
| 8 | Round 7 (optional) | None | S-space reference solver |

Rounds 1, 2, and 3 are independent and can be executed in any order or in parallel. Rounds 4-6 are sequential and depend on prior outputs.

**Carrying forward outputs:** When executing Round 4+, include prior rounds' output *headers* (not implementations) in the `<source_files>` section.

**Verification between rounds:** After each round, verify the output compiles against QuantLib v1.23. Common issues: missing `#include` directives, `ext::shared_ptr` vs `std::shared_ptr`, missing `override` keywords, `Disposable<>` vs plain return types.

**If a round fails:** Re-run with compilation errors appended inside an `<error_context>` tag.

**FdmSchemeDesc const-ness:** Before Round 3, check the actual definition. If members are const, note this in the prompt for the constructor overload approach.

**TripleBandLinearOp API:** Before Round 2, check whether `add()` or node-level setters exist. If not, note this in the prompt so the CN variant can use an alternative assembly strategy.

**Hard stop rule:** If any placeholder file content is missing, the model must STOP and ask — never reconstruct from memory.

**Key mathematical references across rounds:**

| Concept | Paper Reference | Golden Ref | Round |
|:---|:---|:---|:---|
| Fitting factor | §2, E(3) | CE-7 | R1 |
| Fitted matrix A | §2, E(4)/U7 | CE-9 | R1, R7 |
| M-matrix guarantee | §2, C3-C4 | CE-10 | R1, R6(T1) |
| Positivity | §2, C6 | CE-12 | R6(T2) |
| CN variant P, N | §3, U18/U19 | CE-17/CE-18 | R2, R7 |
| Parameter ω | §3, E(9) | CE-19 | R2 |
| Time-step constraint | §3, E(9) | CE-20 | R2, R6(T3) |
| Monitoring projection | §4, E(13) | CE-6 | R2, R6(T6) |
| Artificial diffusion (Scheme 1) | §2, E(7) | CE-15 | R6(T7) |
| Artificial diffusion (Scheme 2) | §3, E(10) | CE-23 | R6(T7) |
| Errata corrections | — | §6 of Round 1.5 | All rounds |

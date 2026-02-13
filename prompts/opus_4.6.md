*Thinking...*

> The user wants me to revise and rewrite their prompting framework for a complex mathematical research project. Let me understand what they're asking:
> 
> 1. They have 2 seed papers (the PDFs attached - Milev & Tagliani 2010, and Duffy 2004)
> 2. They have related papers as review notes
> 3. They want to use Claude Opus 4.6 for this research
> 4. The task has 2 stages:
>    - Stage 1: Devise a novel SOTA algorithm that improves on the 2 seed papers
>    - Stage 2: Create implementation roadmaps using pseudocode and workflow diagrams based on the algorithm and QuantLib framework
> 
> 
> 5. They want the framework revised to follow Anthropic's official Opus 4.6 guidelines, support multi-round prompting for each stage, and prioritize mathematical rigor. The attached QuantLib header files provide context about the existing finite difference, lattice, and Monte Carlo methods framework. The current framework has several structural problems: it's too monolithic and tries to accomplish everything in one output, over-constrains the format with ASCII diagrams and line counts, conflates mathematical derivation with software architecture, and doesn't align with how Claude 4.6 performs best on focused tasks rather than massive combined prompts. Additionally, the stage labeling is confusing since what's called "Stage 1" references a consolidated summary that may already exist, and Stage 2 appears to duplicate Stage 1's prompt structure.
> 
> I should restructure this into genuinely distinct multi-round workflows where Stage 1 progresses through critical analysis of seed papers, literature surveying for SOTA methods, synthesizing a novel algorithm, and formalizing the mathematics with proofs, while Stage 2 separately handles architecture mapping to QuantLib, module-by-module pseudocode development, testing frameworks, and integration guidance. Each round should be self-contained and achievable, with outputs feeding naturally into subsequent rounds, clear artifact naming for cross-referencing, and freedom for Claude to communicate naturally rather than forcing rigid formatting constraints. I'm examining the QuantLib finite difference implementation to understand the distinction between the deprecated older framework with basic schemes like Crank-Nicolson and the newer multi-dimensional framework that uses composite operators, advanced schemes like Douglas and Craig-Sneyd, and specialized boundary conditions for exotic options.

# Revised Multi-Round Prompting Framework for Numerical PDE Research

Below is a complete redesign of your two-stage framework. I've restructured it based on several observations from your materials and what works well with this model.

---

## Design Rationale

Your original framework has three structural problems that will undermine output quality:

**Problem 1: Monolithic prompts.** Each stage tries to extract an entire document in one shot. Claude 4.6 Opus produces its best analytical work when given focused, bounded tasks with clear deliverables — not when asked to fill a 15-section template simultaneously. The token budget inevitably forces shallow treatment of later sections.

**Problem 2: Premature architecture.** Stage 1 jumps directly to "devise a SOTA method" without first establishing what the actual deficiencies are across the literature. The seed papers (Milev–Tagliani 2010 on numerical diffusion, Duffy 2004 on CN pathologies) identify problems but propose different remedies with different trade-offs. A novel method needs a precise gap analysis first.

**Problem 3: QuantLib mapping is under-specified.** Your Stage 2 references QuantLib but the prompt doesn't engage with the actual framework structure. The attached headers reveal QuantLib has *two* FD frameworks — the old deprecated one (`MixedScheme`, `CrankNicolson` template) and the new one (`FdmLinearOpComposite`, `FdmBackwardSolver`, scheme classes). Any implementation roadmap must target the new framework specifically, mapping to its operator/mesher/solver/scheme decomposition.

The revised framework below splits each stage into focused rounds with explicit artifact hand-offs.

---

## Stage 1: Devise a Novel Algorithm (4 Rounds)

### Round 1.1 — Critical Analysis of Seed Papers

```markdown
You are conducting mathematical research on finite difference methods for
Black-Scholes PDE pricing. Two seed papers are attached:

  Paper A: Milev & Tagliani (2010), "Low Volatility Options and Numerical
           Diffusion of Finite Difference Schemes" — Serdica Math. J. 36
  Paper B: Duffy (2004), "A Critique of the Crank-Nicolson Scheme" —
           Wilmott Magazine

Analyze both papers and produce a structured critical assessment with the
following sections. Be mathematically precise; cite equations by their
paper-internal numbering (e.g., "Paper A, Eq. (3)").

SECTION 1 — PROBLEM IDENTIFICATION
For each paper, extract:
  (a) The precise failure modes identified (spurious oscillations, numerical
      diffusion, positivity violation, etc.). State the mathematical
      conditions under which each failure occurs (e.g., σ² ≪ r, non-smooth
      initial data).
  (b) The PDE formulation used (coordinate system, domain, boundary
      conditions). Note any differences between the two papers.
  (c) The diagnostic criteria used to detect failure (visual, norm-based,
      spectral, etc.).

SECTION 2 — PROPOSED REMEDIES
For each paper, extract:
  (a) The finite difference scheme(s) proposed as improvements.
  (b) The theoretical properties claimed (stability, convergence order,
      positivity, M-matrix structure). Reproduce the key result statements
      with their conditions.
  (c) The artificial diffusion expressions derived. For Paper A, this means
      the ½rSΔS and ⅛(rΔS/σ)² terms. Explain what PDE each scheme
      effectively solves.

SECTION 3 — LIMITATIONS AND GAPS
For each paper, identify:
  (a) What the paper does NOT address (e.g., American exercise, Greeks
      accuracy, adaptive methods, multi-factor).
  (b) Where the proposed remedy introduces new problems (e.g., exponential
      fitting introduces O(h) diffusion; CN-variant requires restrictive
      Δt bound from Eq. (9) of Paper A).
  (c) Numerical evidence quality — are the test cases sufficient? What
      parameter regimes are untested?

SECTION 4 — CROSS-PAPER COMPARISON
  (a) Where do the papers agree on diagnosis?
  (b) Where do they disagree or offer incompatible remedies?
  (c) What is the precise relationship between Duffy's exponentially fitted
      scheme (Paper B §4-5) and Milev-Tagliani's presentation of the same
      scheme (Paper A §2)?

Label your output as ARTIFACT: seed_paper_analysis
```

### Round 1.2 — Literature Survey Synthesis

```markdown
Continuing the research from the previous round (reference: seed_paper_analysis).

Below are review notes from related papers in the field. [PASTE REVIEW NOTES HERE]

Using the seed paper analysis and these review notes, produce a structured
survey organized by TECHNIQUE CATEGORY rather than by paper. The categories
are:

CATEGORY A — SPATIAL DISCRETIZATION
Cover: central differences, exponential fitting (Il'in, Duffy), upwind schemes,
compact/Hermite schemes, non-uniform mesh strategies (Tavella-Randall grading,
sinh concentration). For each technique, state: convergence order, positivity
properties, behavior as σ→0, computational cost relative to standard central
differences.

CATEGORY B — TEMPORAL DISCRETIZATION
Cover: fully implicit Euler, Crank-Nicolson (standard), Rannacher smoothing
(implicit startup steps), TR-BDF2, exponential time integration, Richardson
extrapolation in time. For each: order of accuracy, L-stability (yes/no),
damping ratio for high-frequency modes, behavior with non-smooth data.

CATEGORY C — PAYOFF AND INITIAL CONDITION TREATMENT
Cover: direct pointwise sampling, cell averaging, projection onto FE basis,
smoothing by analytic sub-stepping (Rannacher), implicit regularization.
For each: effect on convergence order, interaction with spatial scheme.

CATEGORY D — AMERICAN EXERCISE AND CONSTRAINTS
Cover: penalty methods, policy (Howard) iteration, PSOR, operator splitting.
State convergence properties and per-step cost.

CATEGORY E — GREEKS COMPUTATION
Cover: direct differentiation of discrete solution, PDE-based (Theta from
residual), bump-and-reprice, simultaneous price-delta schemes (Keller box).
State: which Greeks achieve full scheme order vs. lose an order.

For each category, conclude with a GAP STATEMENT: what combination of
properties is NOT achieved by any single existing method?

Label your output as ARTIFACT: literature_survey
```

### Round 1.3 — Novel Algorithm Design

```markdown
Continuing the research (reference: seed_paper_analysis, literature_survey).

Based on the identified gaps, design a novel finite difference algorithm for
the one-factor Black-Scholes PDE that specifically addresses the combined
failure scenario: discontinuous payoff AND low volatility (σ² ≪ r), while
maintaining second-order accuracy for smooth problems.

Your algorithm must integrate techniques from across the surveyed categories
into a coherent scheme. Produce the following:

SECTION 1 — DESIGN CHOICES AND JUSTIFICATION
For each of the five technique categories (A-E from the survey), state which
approach you select and why it is preferred over alternatives in the context
of the combined failure scenario. Be specific about trade-offs.

SECTION 2 — ALGORITHM SPECIFICATION
Present the complete algorithm as a mathematical procedure. Use the following
structure:
  2.1 — Coordinate transformation (state the change of variables)
  2.2 — Spatial grid construction (uniform or graded; specify parameters)
  2.3 — Payoff initialization (smoothing method with formulas)
  2.4 — Spatial operator assembly (full expressions for tridiagonal
         coefficients, fitting factor, M-matrix conditions)
  2.5 — Time stepping scheme (state machine: which method at which steps,
         transition conditions)
  2.6 — Boundary enforcement (per time step)
  2.7 — American constraint enforcement (if applicable; specify method)
  2.8 — Greeks extraction (formulas in the chosen coordinate system)

For each subsection, state the convergence order and the key conditions
(on grid parameters, time step, etc.) required for it to hold.

SECTION 3 — THEORETICAL PROPERTIES
State and justify (proof sketch level, not full proofs):
  (a) Positivity preservation conditions
  (b) Discrete maximum principle conditions
  (c) Convergence order for price, Delta, Gamma — distinguishing the
      smooth-data case from the non-smooth case
  (d) Behavior in the limiting regimes: σ→0, r→0, K=S₀ (ATM)

SECTION 4 — RELATIONSHIP TO EXISTING METHODS
Precisely characterize how your algorithm reduces to known methods under
parameter specialization:
  - With fitting factor = 1: reduces to ___
  - With Rannacher steps = 0: reduces to ___
  - With cell averaging disabled: reduces to ___
  - With σ large relative to r: behaves like ___

Label your output as ARTIFACT: algorithm_design
```

### Round 1.4 — Monolithic Pseudocode and Error Analysis

```markdown
Continuing the research (reference: algorithm_design).

Produce two deliverables:

DELIVERABLE 1 — COMPLETE PSEUDOCODE
Write the full algorithm from algorithm_design §2 as a single self-contained
pseudocode procedure. Use C-like syntax with typed variables. Target 150-250
lines. Requirements:
  - All array indices explicit
  - All special cases handled (boundary nodes, first/last time steps,
    monitoring dates for barrier options)
  - Named constants for tuning parameters with recommended default values
  - Comments referencing algorithm_design §X.Y for each logical block
  - Input: a parameter struct (S₀, K, T, r, q, σ, option_type, exercise_type,
    barriers, monitoring_dates, grid_params)
  - Output: price at S₀, plus Delta, Gamma, Theta

DELIVERABLE 2 — TRUNCATION ERROR ANALYSIS
For the spatial and temporal discretization combined:
  (a) Derive (at proof-sketch level) the leading truncation error terms for
      the fitted spatial operator on a uniform log-price grid
  (b) Derive the leading error terms for the time stepping scheme in both
      the damped (implicit Euler) and normal (CN or TR-BDF2) phases
  (c) State the combined error as O(h^p + k^q) with explicit p, q
  (d) Identify the cross-terms or coupling conditions (e.g., when does the
      spatial fitting factor's O(h) artificial diffusion dominate the O(h²)
      truncation error of central differences?)
  (e) Give the modified equation (the PDE actually solved to leading order)
      for your scheme, analogous to Paper A Eq. (7) and (10)

Label your output as ARTIFACT: pseudocode_and_error_analysis
```

---

## Stage 2: QuantLib Implementation Roadmap (4 Rounds)

### Round 2.1 — QuantLib Framework Mapping

```markdown
You are designing an implementation of a custom FDM pricing engine within
the QuantLib C++ library (version ≥1.36, new FD framework).

INPUTS YOU HAVE:
  1. The algorithm specification (reference: algorithm_design)
  2. The pseudocode (reference: pseudocode_and_error_analysis)
  3. The QuantLib FDM header files (attached as methods.xml) — specifically
     the new framework under finitedifferences/operators/,
     finitedifferences/schemes/, finitedifferences/meshers/, and
     finitedifferences/solvers/

NOTE: QuantLib has TWO FD frameworks. The OLD framework (CrankNicolson<>,
ImplicitEuler<>, MixedScheme<>, TridiagonalOperator — all marked
[[deprecated]]) must NOT be targeted. Target the NEW framework built on:
  - FdmLinearOpComposite (operator interface)
  - FdmMesher / Fdm1dMesher (spatial grid)
  - FdmBackwardSolver + FdmSchemeDesc (time stepping)
  - Fdm1DimSolver / Fdm2DimSolver (top-level solvers)
  - BoundaryCondition<FdmLinearOp> (boundary conditions)
  - StepCondition<Array> (American exercise, barriers)
  - TripleBandLinearOp (tridiagonal spatial operators)

Produce:

SECTION 1 — COMPONENT MAPPING TABLE
Map each logical component from the algorithm to the QuantLib class it should
extend, wrap, or replace:

| Algorithm Component | QuantLib Base Class/Interface | Action |
|---|---|---|
| Spatial grid | Fdm1dMesher | Extend: new mesher class |
| Spatial operator | TripleBandLinearOp / FdmLinearOpComposite | ... |
| Time stepper | FdmSchemeDesc + scheme class | ... |
| ... | ... | ... |

For "Action", use one of: REUSE (use existing class as-is), EXTEND (subclass),
COMPOSE (wrap existing classes), REPLACE (new class implementing existing
interface), NEW (no existing interface fits).

SECTION 2 — NEW CLASSES REQUIRED
For each class requiring EXTEND, REPLACE, or NEW:
  (a) Proposed class name (following QuantLib naming conventions)
  (b) Parent class / interface implemented
  (c) Key methods to implement (signature sketch, not full pseudocode)
  (d) Which QuantLib header to model after (closest existing analog)

SECTION 3 — DATA FLOW THROUGH QUANTLIB
Trace the pricing call path through QuantLib's existing solver infrastructure:
  FdmBlackScholesSolver → Fdm1DimSolver → FdmBackwardSolver →
  [scheme].step() → [operator].apply() / solve_splitting()

Show where your custom components plug in and what existing flow they modify.
Use an ASCII diagram showing the call chain with your new classes highlighted.

SECTION 4 — INTEGRATION CONSTRAINTS
List concrete constraints imposed by QuantLib's architecture:
  (a) Thread safety requirements (QuantLib uses LazyObject pattern)
  (b) Memory ownership model (ext::shared_ptr everywhere)
  (c) Time direction convention (QuantLib solves backward: maturity→0)
  (d) Coordinate convention (QuantLib's FdmBlackScholesMesher works in
      ln(S) space — confirm compatibility with algorithm_design §2.1)
  (e) Any incompatibilities between the algorithm's requirements and
      QuantLib's existing interfaces

Label your output as ARTIFACT: quantlib_mapping
```

### Round 2.2 — Module Pseudocode (Core Operator and Mesher)

```markdown
Continuing implementation design (reference: quantlib_mapping,
algorithm_design, pseudocode_and_error_analysis).

Produce detailed pseudocode for the two most mathematically critical new
QuantLib components. These are the SPATIAL components — the mesher and
the fitted operator — which contain the core numerical innovation.

MODULE 1 — FITTED MESHER (extends Fdm1dMesher)
  Purpose: Builds a 1D log-price grid with optional sinh-concentration
  near strike/barriers.

  Pseudocode requirements:
  - Constructor logic: compute grid points, dplus/dminus arrays
  - Two strategies: UniformLogPrice, SinhGraded
  - For SinhGraded: specify the concentration formula, the parameter
    selection heuristic, and barrier-node alignment logic
  - Show how the output (locations_, dplus_, dminus_ arrays) is populated
  - 30-50 lines

MODULE 2 — FITTED SPATIAL OPERATOR (extends FdmLinearOpComposite)
  Purpose: Assembles the exponentially fitted tridiagonal operator for
  the 1D Black-Scholes PDE in log-price coordinates.

  Pseudocode requirements:
  - setTime(t1, t2): compute vol, drift, discount at current time
  - apply(r): multiply operator by array r (explicit application)
  - solve_splitting(direction, r, s): solve (I - s·L)·x = r via Thomas
  - Internal: fitting factor computation per algorithm_design §2.4
  - Internal: M-matrix verification, fallback to upwind if violated
  - Show the tridiagonal coefficient assembly in full detail (this is the
    key numerical content — do not abbreviate)
  - Show how TripleBandLinearOp's lower_/diag_/upper_ arrays are set
  - 50-80 lines

For both modules, include:
  - Preconditions (what must be true about inputs)
  - Postconditions (what the module guarantees about outputs)
  - Complexity (time and space)

Label your output as ARTIFACT: core_module_pseudocode
```

### Round 2.3 — Module Pseudocode (Time Stepping and Constraints)

```markdown
Continuing implementation design (reference: quantlib_mapping,
algorithm_design, core_module_pseudocode).

Produce detailed pseudocode for the time-stepping and constraint-enforcement
components.

MODULE 3 — RANNACHER-SMOOTHED CN SCHEME
  Purpose: Time stepper implementing Rannacher startup (m_damp implicit
  Euler half-steps) followed by standard Crank-Nicolson.

  QuantLib integration: This should work with FdmBackwardSolver's scheme
  dispatch. Examine how existing schemes (ImplicitEulerScheme,
  CrankNicolsonScheme, DouglasScheme) implement the step(a, t) and
  setStep(dt) interface.

  Pseudocode requirements:
  - State machine with states: RANNACHER_DAMPING, NORMAL_CN,
    POST_MONITOR_DAMPING (for discrete barrier monitoring dates)
  - step(a, t): the full logic including state transitions
  - Show how implicit Euler half-steps compose (two solves per full step
    during damping phase)
  - Show the CN step as the average of explicit and implicit applications
  - Handle monitoring dates: detect proximity, trigger re-damping
  - 40-60 lines

MODULE 4 — TR-BDF2 SCHEME (alternative time stepper)
  Purpose: L-stable second-order scheme as alternative to Rannacher-CN.

  Pseudocode requirements:
  - Two-stage update: trapezoidal stage then BDF2 stage
  - γ = 2 - √2 parameter
  - Show both tridiagonal solves
  - 25-40 lines

MODULE 5 — AMERICAN EXERCISE STEP CONDITION (extends StepCondition<Array>)
  Purpose: Enforces early exercise after each time step.

  Two strategy variants:
  (a) POLICY ITERATION: Show the iteration loop with active-set tracking.
      applyTo(a, t): iterate until active set stabilizes.
      15-30 lines.
  (b) PENALTY METHOD: Show the penalized system modification.
      15-25 lines.

MODULE 6 — DISCRETE BARRIER STEP CONDITION (extends StepCondition<Array>)
  Purpose: At monitoring dates, zeros the solution outside [L, U] corridor.

  Pseudocode requirements:
  - applyTo(a, t): check if t is a monitoring date, apply corridor projection
  - Signal to time stepper that re-damping is needed (via a shared flag or
    callback — specify the mechanism compatible with QuantLib's
    FdmStepConditionComposite)
  - 15-25 lines

Label your output as ARTIFACT: stepping_module_pseudocode
```

### Round 2.4 — Integration, Testing, and Workflow Diagrams

```markdown
Continuing implementation design (reference: quantlib_mapping,
core_module_pseudocode, stepping_module_pseudocode, algorithm_design).

Produce three final deliverables:

DELIVERABLE 1 — ORCHESTRATOR / SOLVER CLASS
Pseudocode for the top-level solver class (analogous to FdmBlackScholesSolver)
that wires all custom modules together.

Show the complete solve lifecycle:
  1. Construct mesher (Module 1)
  2. Construct operator (Module 2)
  3. Construct boundary conditions
  4. Construct step conditions (Module 5 and/or 6)
  5. Construct time stepper (Module 3 or 4)
  6. Initialize solution array (payoff processing per algorithm_design §2.3)
  7. Call FdmBackwardSolver::rollback (or equivalent)
  8. Extract price and Greeks from solution

Show how product type (European/American/Barrier) dispatches to different
strategy selections. Use a configuration table:

| Product | Mesher Strategy | Scheme | StepConditions | Boundary |
|---|---|---|---|---|

30-50 lines of pseudocode.

DELIVERABLE 2 — GREEKS EXTRACTION
Pseudocode for extracting Delta, Gamma, Theta from the solved array.
  - Delta and Gamma: from the log-price solution array using the formulas
    from algorithm_design §2.8, with interpolation to off-grid S₀
  - Theta: from PDE residual at the solution point
  - Vega, Rho: bump-and-reprice protocol (specify bump size, one-sided vs
    central, what gets reconstructed per bump)
  20-35 lines.

DELIVERABLE 3 — VALIDATION TEST SUITE SPECIFICATION
Define 6-8 test cases in a table:

| ID | Product | S₀ | K | T | r | q | σ | Reference Value | Source | Convergence Order | Failure Mode Tested |

Include:
  - At least one smooth case (EU vanilla, moderate σ) → verify O(h²+k²)
  - At least one low-vol case (σ=0.001, r=0.05) → verify no artificial diffusion
  - At least one discontinuous payoff (digital or truncated call)
  - At least one American option
  - At least one discrete-monitored barrier
  - At least one ATM case testing Greeks accuracy

For each test, specify how the reference value is obtained (BS closed-form,
high-resolution numerical, published benchmark).

Provide the convergence test procedure as a workflow diagram (ASCII):
  Run at J=100,N=50 → Run at J=200,N=100 → Run at J=400,N=200 →
  Compute error ratios → Check order → PASS if ratio ∈ [3.5, 4.5] for O(h²)

Label your output as ARTIFACT: integration_and_testing
```

---

## Usage Notes

**Artifact chaining.** Each round's output is labeled with an `ARTIFACT:` tag. In subsequent rounds, reference prior artifacts by name. Claude maintains context within a conversation, so you don't need to re-paste prior outputs — just reference them. If you start a new conversation, paste the relevant prior artifacts into the prompt.

**Handling the review notes.** In Round 1.2, you'll need to provide the review notes from the follow-up papers. Format them as a structured list with author, year, key contribution, and relevant equations. Don't paste full papers — distill to the technically critical content.

**Adapting to intermediate results.** If Round 1.3 produces an algorithm that differs significantly from what you expected (e.g., it recommends TR-BDF2 over Rannacher-CN as the primary scheme), adjust the subsequent round prompts accordingly. The framework is designed so each round's prompt can be fine-tuned based on what the previous round produced.

**QuantLib version specificity.** The attached headers show classes deprecated in version 1.42. The new framework has been stable since approximately version 1.20. If you're targeting a specific version, add a note to Round 2.1 specifying it. The key interfaces (`FdmLinearOpComposite`, `FdmBackwardSolver`, `FdmSchemeDesc`) have been stable across recent versions.

**What I removed from your original framework and why:**

The original prompt demanded ASCII box-drawing diagrams with minimum line counts, specific table column layouts, and exhaustive coverage of extension hooks (Heston 2D, multi-asset, jump-diffusion). These are aspirational but competing for token budget with the core technical content. The revised framework focuses the token budget on the novel algorithm and its implementation, which is where mathematical rigor matters most. Extensions can be a separate Stage 3 if needed.

The original also specified a `SYSTEM_ROLE` asserting Claude is a "principal software architect at a tier-1 quantitative finance institution." This kind of role-play instruction is unnecessary for Claude 4.6 Opus — it already engages deeply with technical material when the task is clearly specified. The revised framework gets the same quality by specifying the *task* precisely rather than the *persona*.

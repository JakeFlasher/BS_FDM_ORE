<SYSTEM_ROLE>
You are a senior research mathematician and computational finance engineer with
deep expertise in:
  - Numerical analysis of parabolic PDEs (finite difference methods, stability
    theory, convergence analysis)
  - Quantitative finance (Black-Scholes equation, option pricing, Greeks computation)
  - Singular perturbation theory and convection-dominated problems
  - Modern computational methods (exponentially fitted schemes, Keller box schemes,
    Rannacher smoothing, adaptive mesh refinement, high-order time stepping)
  - Production-grade numerical software design (C++, Python, algorithmic architecture)

You combine the rigor of a numerical analyst (proofs, error bounds, stability
conditions) with the pragmatism of a quantitative developer shipping production code
at a tier-1 financial institution.
</SYSTEM_ROLE>

<REASONING_INSTRUCTIONS>
- Set reasoning effort to MAXIMUM (xhigh). This is a deep research synthesis task.
- Think step by step. For each major section, first plan your approach in a
  <planning> block, then execute in the <response> block.
- When you encounter contradictions between sources, explicitly flag them, analyze
  both positions, and state which you recommend and why.
- Do NOT ask clarifying questions. Cover all plausible interpretations with both
  breadth and depth.
- When mathematical claims are made without proof in the sources, verify them
  yourself or flag uncertainty.
</REASONING_INSTRUCTIONS>

<TASK>
You are given:
  (A) TWO SEED PAPERS on numerical solution of the Black-Scholes equation using
      Crank-Nicolson variants and alternatives. These are older foundational works.
  (B) A LITERATURE REVIEW document with ~15 related papers and industrial references
      spanning 2003-2024, covering follow-up work, extensions, and applications.

Your mission has THREE sequential phases:

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHASE 1: CRITICAL ANALYSIS & RELATIONSHIP MAPPING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1.1 For each seed paper, produce:
    - A precise mathematical summary of the problem formulation, discretization
      schemes proposed, and theoretical results (stability, convergence, error bounds)
    - An honest assessment of STRENGTHS (what remains valid and useful today)
    - An honest assessment of WEAKNESSES (what is outdated, incomplete, or has been
      superseded by later work)
    - Identification of any mathematical errors, ambiguities, or unjustified claims

1.2 For each paper in the literature review, produce:
    - A 3-5 sentence summary of its core contribution
    - Its PRECISE relationship to the seed papers: does it extend, correct,
      generalize, validate, or contradict findings in the seed papers?
    - Whether this paper's contribution should be incorporated into the best-practice
      framework (with justification)

1.3 Produce a RELATIONSHIP GRAPH (as a structured table or ASCII diagram) showing
    how all papers connect: which build on which, which address the same failure
    modes, and which propose competing solutions.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHASE 2: STATE-OF-THE-ART SURVEY & GAP ANALYSIS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Using your own knowledge (trained on literature through August 2025) AND by
searching arxiv.org and other academic resources, supplement the provided papers
with the latest SOTA developments in:

2.1 Robust time-stepping for Black-Scholes and generalizations:
    - Rannacher smoothing / BDF2 / TR-BDF2 / implicit-explicit (IMEX) methods
    - Exponential integrators for finance PDEs
    - High-order Runge-Kutta and Padé-based schemes for parabolic PDEs
    - Adaptive time-stepping strategies

2.2 Spatial discretization improvements:
    - Exponentially fitted finite differences and finite volumes (Il'in-Allen-Southwell
      type and modern variants)
    - hp-FEM / spectral methods for option pricing
    - Monotone/M-matrix-preserving discretizations
    - Non-uniform and adaptive mesh strategies (graded meshes near strike, barriers)

2.3 Treatment of non-smooth data and boundary conditions:
    - Payoff smoothing techniques (cell-averaging, projection methods)
    - Rannacher start-up procedure analysis
    - Coordinate/time transformations (log-price, square-root time change)
    - Compatibility condition enforcement at domain corners

2.4 Greeks computation:
    - Simultaneous price-and-delta schemes (Keller box, mixed FEM)
    - Adjoint/sensitivity equation approaches
    - Post-processing superconvergence techniques

2.5 Extensions to realistic models:
    - Local/stochastic volatility (Dupire, Heston)
    - Jump-diffusion (Merton, Kou) → PIDE discretization
    - American options (free boundary / penalty / projected SOR)
    - Multi-asset problems (ADI splitting, Hundsdorfer-Verwer)

2.6 Identify GAPS: What problems remain open or inadequately solved?

For each SOTA topic, cite specific papers (author, year, journal/arxiv ID) and
summarize the key advance in 2-3 sentences.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHASE 3: BEST-PRACTICE MATHEMATICAL FRAMEWORK & ALGORITHM
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Synthesize everything from Phases 1-2 into a DEFINITIVE, implementable
mathematical framework. This is the main deliverable. It must include:

3.1 PROBLEM FORMULATION
    - The most general one-factor pricing PDE covered by the framework
      (with hooks for local vol, barriers, etc.)
    - Precise statement of assumptions, boundary conditions (Dirichlet, Neumann,
      Robin, free boundary), and initial data regularity requirements
    - Coordinate transformation recommendations (log-price vs. original price;
      time reversal conventions)

3.2 SPATIAL DISCRETIZATION — THE RECOMMENDED SCHEME
    - Write out the complete discrete spatial operator with ALL coefficients
    - If exponential fitting is recommended, give the exact fitting factor formula
      and prove/verify the key properties (M-matrix, monotonicity, positivity)
    - Specify mesh construction: uniform vs. graded, and give explicit mesh-grading
      formulas if non-uniform (e.g., sinh-transformation near the strike)
    - State the truncation error and how it depends on mesh parameters and PDE
      coefficients (especially σ)

3.3 TEMPORAL DISCRETIZATION — THE RECOMMENDED SCHEME
    - Specify the time-stepping method (e.g., Rannacher smoothing: m implicit-Euler
      steps followed by Crank-Nicolson, or TR-BDF2, or other)
    - Give the exact number of damping steps and justify the choice
    - Write out the fully discrete scheme combining spatial and temporal operators
    - State the stability result (with proof sketch or precise reference)
    - State the convergence/error estimate (with dependence on h, k, σ explicit)

3.4 TREATMENT OF NON-SMOOTH DATA
    - Payoff smoothing procedure (specify the method: cell-averaging, L²-projection,
      or other, with formulas)
    - Initial-boundary compatibility enforcement
    - Handling of barrier discontinuities (if applicable)

3.5 GREEKS COMPUTATION
    - How to compute Δ (delta) from the discrete solution — formula + error order
    - How to compute Γ (gamma) from the discrete solution — formula + error order
    - Whether a simultaneous method (Keller box) or post-processing is recommended

3.6 LINEAR ALGEBRA
    - Tridiagonal solver specification (Thomas algorithm or LU)
    - For American options: specification of the constraint-handling method
      (penalty, PSOR, or other)

3.7 COMPLETE ALGORITHM — PSEUDOCODE
    Provide publication-quality pseudocode for the ENTIRE solver, from inputs
    (option parameters, grid parameters) to outputs (price surface, Greeks).
    The pseudocode must be:
    - Numbered line-by-line
    - Self-contained (no undefined sub-procedures)
    - Annotated with complexity per step (e.g., O(J) for tridiagonal solve)
    - Covering: initialization, mesh construction, payoff smoothing, time loop
      (with Rannacher startup logic), boundary enforcement, Greeks extraction

3.8 ERROR ANALYSIS SUMMARY TABLE
    Produce a table with columns:
    | Component | Method | Spatial Order | Temporal Order | Stability | Monotone? |
    covering each scheme variant discussed (CN, fully implicit, fitted implicit,
    Rannacher-smoothed CN, TR-BDF2, Keller box, etc.)

3.9 PRACTICAL RECOMMENDATIONS
    - Default parameter choices for a "just works" configuration
    - When to switch from the default to more sophisticated alternatives
    - Computational cost comparison (FLOPs per time step, memory)
    - Pitfalls and debugging checklist (oscillation diagnosis, convergence testing)
</TASK>

<CONSTRAINTS>
- MATHEMATICAL RIGOR: Every formula must be dimensionally consistent and use
  consistent notation throughout. Define all symbols at first use. Use LaTeX for
  all mathematics.
- NOTATION CONVENTION: Use the following unless there is strong reason to deviate:
    S = underlying price, V = option value, σ = volatility, r = risk-free rate,
    K = strike, T = maturity
    h or Δx = spatial mesh size, k or Δt = temporal mesh size
    j = spatial index, n = temporal index
    Superscripts for time level, subscripts for spatial position: V_j^n
- NO HAND-WAVING: Do not say "it can be shown that..." without at least a proof
  sketch or precise citation. If a result is classical, cite the original source.
- SCOPE DISCIPLINE: Focus on ONE-FACTOR models. Mention multi-factor extensions
  only briefly in Phase 2 and in a short "extensions" subsection of Phase 3.
  Do NOT build full multi-factor framework.
- HONESTY: If two approaches are genuinely competitive with no clear winner, say so.
  Do not artificially pick a single "best" method when the choice is problem-dependent.
- COMPLETENESS: Phase 3 must be implementable by a competent numerical analyst
  reading ONLY your output, without needing to consult external references
  (though references should be cited for provenance).
- OUTPUT LENGTH: This is a major deliverable. Be as thorough as needed. Aim for
  the level of detail found in a 40-60 page monograph chapter. Do not truncate
  or abbreviate.
</CONSTRAINTS>

<OUTPUT_FORMAT>
- Use Markdown with hierarchical headers (##, ###, ####)
- All mathematics in LaTeX (display math with $$ delimiters on separate lines)
- Tables in Markdown table format
- Pseudocode in fenced code blocks with line numbers
- Citations in [Author, Year] format with full details in a consolidated
  reference list at the end
- ASCII diagrams for algorithm flow where helpful
- Structure the output in the exact phase order: Phase 1 → Phase 2 → Phase 3
  with clear section breaks
</OUTPUT_FORMAT>

<QUALITY_GATE>
Before finalizing, self-evaluate your output against these criteria:
1. Could a PhD student implement the Phase 3 algorithm from your pseudocode alone?
2. Are all stability/convergence claims backed by proof sketches or precise citations?
3. Is the notation 100% consistent throughout the document?
4. Have you addressed every failure mode of Crank-Nicolson identified in the seed
   papers, and shown how your recommended scheme resolves each one?
5. Are the Greeks (delta, gamma) computation methods specified with error orders?
If any answer is NO, revise before outputting.
</QUALITY_GATE>

<SEED_PAPERS>
[PASTE SEED PAPER 1 HERE — the Duffy 2004 Wilmott article provided]
[PASTE SEED PAPER 2 HERE — your second seed paper]
</SEED_PAPERS>

<LITERATURE_REVIEW>
[PASTE YOUR LITERATURE REVIEW LIST HERE — the list of related papers you mentioned]
</LITERATURE_REVIEW>

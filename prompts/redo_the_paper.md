*Thinking...*

> The user is shifting the goal from creating a novel algorithm to implementing an existing method from a single seed paper. The key challenges are:
> 
> 1. **Single seed paper implementation**: They want to implement exactly what's described in the paper, not create something new
> 2. **Mathematical gap-filling**: Papers often omit proofs, hide details, or skip steps that are crucial for implementation
> 3. **Verification of logic chains**: They want to ensure no weak logic or incorrect mathematical reasoning exists in the paper's presentation
> 4. **Supplementary materials**: They'll provide related papers/notes to fill in gaps the seed paper leaves
> 
> This is fundamentally different from the previous framework:
> - **Before**: Survey → Synthesize → Design novel method → Implement
> 
> 
> - **Now**: Deep-read single paper → Identify gaps/omissions → Fill gaps using supplementary materials → Verify mathematical correctness → Implement exactly as specified
> 
> The restructured rounds should focus on extracting every mathematical detail from the seed paper, pinpointing what's missing or underspecified, using supplementary materials to complete the picture, verifying each logical chain and formula, and finally producing implementation-ready pseudocode that faithfully represents the paper's method.
> 
> Since the user mentions providing only one seed paper going forward, I'm designing this as a generic framework for any single paper rather than targeting a specific one. The revised approach should have a mathematical extraction and verification stage as the foundation, followed by gap-filling and implementation stages. I'm realizing the user wants me to adapt the framework I created to specifically handle paper implementation rather than algorithm design, so I need to revise the GPT-5.2 Pro structure to emphasize extraction, gap-filling, and verification over innovation. The core challenge is that papers use compressed notation and skip steps that seem obvious to authors but aren't to implementers, so the framework needs to account for unpacking that implicit knowledge.

This is a fundamentally different task from "design a novel method." The core challenge shifts from *synthesis* to *forensic mathematical reconstruction* — you're reverse-engineering a paper's full reasoning chain, filling in what the authors left implicit, catching any errors, and producing an implementation-complete specification that is faithful to the paper but has no gaps.

Here's the revised framework.

---

## Structural Redesign Rationale

| Previous Goal (Novel Method) | New Goal (Faithful Implementation) |
|---|---|
| Survey literature broadly | Deep-read single paper forensically |
| Identify gaps in the *field* | Identify gaps in the *paper* |
| Synthesize techniques from multiple sources | Fill gaps using supplementary materials |
| Design novel algorithm | Reconstruct paper's algorithm exactly |
| Prove properties of *your* method | Verify properties *claimed* by the paper |
| Pseudocode for *your* design | Pseudocode faithful to *paper's* method |

The fundamental cognitive tasks are now: **Extract → Identify Gaps → Fill Gaps → Verify → Specify → Implement**. Stage 1 handles the first five (mathematical due diligence), Stage 2 handles the last (implementation).

---

## Stage 1: Mathematical Due Diligence (5 Rounds)

### Round 1.1 — Complete Mathematical Extraction

The goal is to produce a structured inventory of *everything* the paper states — equations, definitions, theorems, conditions, assumptions — with their dependency chain mapped.

```markdown
<context>
One seed paper is attached for faithful implementation:
  [Paper] [TITLE, AUTHORS, VENUE, YEAR — fill in when using]

<supplementary>
[Paste review notes from related papers, each labeled Note-1, Note-2, ...
 These are NOT alternative methods — they are reference material for
 filling gaps in the seed paper.]
</supplementary>
</context>

<task>
Perform a complete mathematical extraction of the seed paper.
Produce five deliverables:

1. NOTATION REGISTRY: Every symbol used, with:
   | Symbol | Meaning | First appearance (§/Eq) | Type (scalar/vector/operator/set) |
   Flag any symbol that is overloaded (used with different meanings)
   or undefined (used without introduction).

2. EQUATION CHAIN: Every numbered and unnumbered equation, in order
   of appearance, with:
   | Eq ID | Expression (reproduced in LaTeX) | Depends on (which prior Eqs) | Role |
   Role ∈ {definition, derived result, stated without proof, assumption,
   boundary/initial condition, discretization, algorithm step}

3. THEOREM/CLAIM INVENTORY: Every result the paper claims
   (theorems, lemmas, propositions, remarks presented as facts), with:
   | Claim ID | Statement (precise) | Proof provided? | Conditions required |

4. ASSUMPTION INVENTORY: Every assumption, both explicit and implicit.
   Explicit: stated by the paper ("we assume σ² > r").
   Implicit: unstated but required for a derivation to hold (e.g.,
   "the truncation S_max is large enough that boundary effects are
   negligible" — stated nowhere but assumed throughout).
   | Assumption ID | Statement | Explicit/Implicit | Used in (Eq/Claim refs) |

5. LOGICAL DEPENDENCY GRAPH: An ordered list showing the deductive
   structure: which results depend on which prior results/assumptions.
   Format as a DAG in text:
   A1, A2 → Eq(1) → Eq(2)
   A3, Eq(2) → Claim 1
   etc.
</task>

<constraints>
- Reproduce every equation exactly as stated in the paper, including
  any notation inconsistencies. Do NOT correct or improve at this stage.
- If the paper uses a result from a referenced work without restating
  it, log it in the Theorem inventory as "external, from [ref X]" and
  flag it as a gap.
- If a derivation step is unclear, log it but do NOT attempt to
  resolve it — that is Round 1.2.
- Cover the ENTIRE paper, including appendices and figure captions.
- Do not skip "obvious" equations — implementation requires every one.
</constraints>

<output_format>
- Five sections with Markdown tables as specified
- LaTeX notation for all equations
- Dependency graph as indented text with arrows
- Aim for exhaustive coverage over brevity
</output_format>
```

**Reasoning effort:** `high`

---

### Round 1.2 — Gap Identification & Classification

```markdown
<context>
Continuing from Round 1.1.
<prior_output>
[Paste Round 1.1 extraction in full — use complete paste, not summary,
 because every equation reference matters here]
</prior_output>

The seed paper is still attached for reference.
</context>

<task>
Systematically audit the paper's mathematical content for gaps,
errors, and ambiguities. For each issue found, classify and catalog it.

Audit procedure — for EACH equation transition (Eq(n) → Eq(n+1)):
  - Can Eq(n+1) be derived from Eq(n) and stated assumptions?
  - If yes: mark VERIFIED
  - If no: classify the gap (see categories below)

For EACH theorem/claim:
  - Is the proof complete? (all steps follow from prior results)
  - Are the conditions sufficient? (or is something stronger needed)
  - Are the conditions necessary? (or are they overly restrictive)

Gap categories:
  G1 — MISSING STEP: A derivation step is skipped. The result is
       likely correct but intermediate algebra is omitted.
  G2 — UNSTATED ASSUMPTION: The derivation requires an assumption
       not listed in the paper. State what assumption is needed.
  G3 — EXTERNAL DEPENDENCY: The result relies on a theorem from
       another paper that is not restated. Cite what is needed.
  G4 — SUSPECTED TYPO: A sign, factor, index, or coefficient appears
       wrong based on dimensional analysis or limiting-case checks.
       State what you suspect it should be and why.
  G5 — AMBIGUITY: The paper's statement admits multiple interpretations.
       State the alternatives.
  G6 — INCOMPLETENESS: The paper handles a special case but the
       general case (needed for implementation) is not addressed.
       Examples: boundary treatment only for Dirichlet but not Neumann;
       coefficient formulas only for constant σ but implementation
       needs σ(S,t); payoff formula only for calls but not puts.
  G7 — PROOF GAP: A step in a proof does not follow from the
       preceding steps, or the proof strategy has a logical flaw.

Produce:
1. GAP REGISTRY:
   | Gap ID | Location (§/Eq/Claim) | Category (G1–G7) | Description | Severity |
   Severity ∈ {critical: blocks implementation, significant: affects
   correctness, minor: cosmetic or easily resolved}

2. VERIFICATION LOG: For each equation transition, one line:
   | From → To | Status (VERIFIED / GAP) | Gap ID if applicable |

3. PROOF AUDIT: For each theorem/claim with a proof provided:
   - List the proof steps
   - Mark each step VERIFIED or flag the gap
   - Overall verdict: {sound, has fixable gaps, fundamentally flawed}

4. IMPLEMENTATION BLOCKERS: Subset of gaps that MUST be resolved
   before implementation can proceed, ordered by dependency
   (resolve upstream gaps first).
</task>

<constraints>
- Do NOT resolve gaps in this round. Identification only.
  Resolution is Round 1.3.
- Be conservative: if a step is non-obvious but you can see how it
  works after thought, mark it VERIFIED with a brief note, not as
  a gap.
- For suspected typos (G4), provide your reasoning (e.g., "dimensional
  analysis shows this term has units of 1/S², but surrounding terms
  have units of 1/S, suggesting a missing S factor").
- If the supplementary notes from related papers shed light on a gap,
  note this as "potentially resolvable via Note-X" but still log it.
- Every gap must have a unique ID (G1.1, G1.2, etc.) for reference
  in later rounds.
</constraints>

<output_format>
- Four sections with tables and structured text
- Verification log may be long — that is expected and desired
- Proof audits as numbered step lists
- Target: exhaustive, no length constraint
</output_format>
```

**Reasoning effort:** `xhigh`

---

### Round 1.3 — Gap Resolution

```markdown
<context>
Continuing from Rounds 1.1–1.2.
<prior_outputs>
[Paste Round 1.1 extraction and Round 1.2 gap registry in full]
</prior_outputs>

<supplementary>
[Same supplementary notes as Round 1.1. If the user has additional
 reference materials (textbooks, related papers) that specifically
 address identified gaps, add them here with labels.]
</supplementary>

The seed paper is still attached.
</context>

<task>
Resolve every gap in the Implementation Blockers list from Round 1.2,
plus as many other gaps as possible. For each gap, produce a
RESOLUTION ENTRY:

Format per gap:
  GAP [ID]: [brief description from Round 1.2]
  RESOLUTION TYPE: {derived, sourced, corrected, disambiguated, generalized}
  RESOLUTION:
    [Full mathematical derivation, citation, correction, or clarification.
     Show ALL intermediate steps — this is the gap-filling, so nothing
     should be left implicit.]
  VERIFICATION:
    [How you know the resolution is correct. Options:
     - Dimensional check
     - Limiting case check (e.g., as σ→0, as h→0)
     - Consistency with other equations in the paper
     - Agreement with supplementary reference [Note-X]
     - Independent re-derivation from first principles]
  STATUS: {resolved, partially resolved — state what remains,
           unresolvable — state what additional information is needed}

Specific instructions by gap category:

  G1 (missing step): Provide the full intermediate algebra.
  G2 (unstated assumption): State the assumption precisely, verify
      it is reasonable, and check whether the paper's results still
      hold under weaker alternatives.
  G3 (external dependency): Find and restate the needed result from
      the supplementary materials if available, or derive it. If the
      result cannot be sourced, state what is needed for the user to
      supply.
  G4 (suspected typo): Derive the correct expression from scratch.
      Show the derivation step-by-step so the correction is verifiable.
  G5 (ambiguity): Determine which interpretation is consistent with
      the rest of the paper. If both are consistent, state the
      implementation-relevant difference.
  G6 (incompleteness): Provide the general-case formula, derived
      by extending the paper's approach.
  G7 (proof gap): Provide the missing proof step, or identify a
      corrected proof strategy.
</task>

<constraints>
- Show ALL intermediate steps in derivations. The purpose of this
  round is to produce what the paper omitted.
- For G4 corrections: the corrected formula MUST be independently
  derived, not just "fixed by inspection." Provide the derivation.
- If a gap cannot be fully resolved with available materials, state
  precisely what is missing and mark as "partially resolved."
  Do not fabricate or speculate.
- Cross-reference: every resolution should cite which equations from
  Round 1.1 it connects to (e.g., "This fills the step between
  Eq(5) and Eq(6)").
- If a resolution changes the downstream dependency graph (e.g.,
  correcting a formula changes a later result), flag this explicitly.
</constraints>

<output_format>
- One subsection per gap, using the format above
- LaTeX for all mathematical content
- Derivations should be step-by-step, not compressed
- At the end: a summary table:
  | Gap ID | Status | Downstream impact (which later Eqs/Claims affected) |
</output_format>
```

**Reasoning effort:** `xhigh`

---

### Round 1.4 — End-to-End Verification

```markdown
<context>
Continuing from Rounds 1.1–1.3. The gap resolutions from Round 1.3
are now integrated.
<prior_outputs>
[Paste Round 1.1 equation chain, Round 1.2 gap registry,
 Round 1.3 resolutions]
</prior_outputs>
</context>

<task>
Before writing the implementation specification, perform a
complete end-to-end verification of the paper's method as corrected
and completed by Round 1.3.

Before starting: plan the verification by listing the 3–5 most
critical logical chains that must hold for the method to work.
Output this as a <planning> block.

Then verify:

1. CORRECTED EQUATION CHAIN: Reproduce the complete chain of
   equations from PDE to final discrete algorithm, incorporating
   Round 1.3 corrections and gap-fills. For each step, confirm it
   follows from the previous step. This is the AUTHORITATIVE chain
   — it supersedes the paper where corrections were made.

2. DIMENSIONAL CONSISTENCY: For every equation in the corrected
   chain, verify that both sides have the same dimensions/units.
   Report as a table:
   | Eq ID | LHS dimensions | RHS dimensions | Status |

3. LIMITING CASE VERIFICATION: Trace the method through each regime:
   (a) σ → 0: Does the scheme degrade gracefully (e.g., to upwind)?
       Reproduce the limiting coefficients explicitly.
   (b) r → 0, q → 0: Do the coefficients simplify correctly?
   (c) h → 0 with k fixed: Is consistency with the continuous PDE
       recovered?
   (d) Constant coefficients: Does the scheme reduce to the known
       textbook form?

4. PROPERTY VERIFICATION: For each property the paper claims
   (positivity, M-matrix, stability, convergence order, maximum
   principle):
   - Restate the claim with corrected conditions
   - Grade: {fully verified, verified with caveats, not verified}
   - If caveats: state what additional condition or proof step is needed

5. IMPLEMENTATION READINESS CHECKLIST:
   | Component | Fully specified? | Any remaining ambiguity? |
   Components: coordinate transform, grid construction, payoff init,
   operator assembly, time stepping, boundary conditions, Greeks
   extraction, parameter constraints.
</task>

<constraints>
- The corrected equation chain in §1 is the most important deliverable.
  It must be complete and self-contained — a reader should be able to
  implement from this chain alone.
- Where Round 1.3 corrected a formula, show BOTH the original and
  corrected versions, clearly labeled.
- For limiting cases, do the algebra explicitly — do not just state
  "it simplifies correctly."
- If any component is NOT fully specified after Rounds 1.1–1.3,
  flag it clearly. Do not paper over remaining gaps.
</constraints>

<output_format>
- Open with <planning> block
- Five numbered sections
- Corrected equation chain: numbered, with LaTeX, step-by-step
- Tables for dimensional check and readiness checklist
- Target: thorough rather than brief
</output_format>
```

**Reasoning effort:** `xhigh`

---

### Round 1.5 — Implementation-Ready Mathematical Specification

```markdown
<context>
Continuing from Rounds 1.1–1.4. Round 1.4 produced the verified
corrected equation chain and readiness checklist.
<prior_output>
[Paste Round 1.4 in full — especially the corrected equation chain
 and any remaining caveats]
</prior_output>
</context>

<task>
Produce the GOLDEN REFERENCE DOCUMENT: a single self-contained
mathematical specification of the paper's method, corrected and
completed, from which a developer can implement without consulting
the original paper.

Structure:

1. METHOD OVERVIEW: 1 paragraph stating what the method does, what
   PDE it solves, and what class of problems it targets.

2. NOTATION TABLE: Complete, consistent, no overloading.

3. MATHEMATICAL SPECIFICATION — sequential subsections:
   (a) Governing PDE with all coefficients explicit
   (b) Coordinate transformation (if any) with transformed PDE
   (c) Domain and boundary conditions (all cases)
   (d) Initial/payoff conditions (all product types the paper covers)
   (e) Spatial discretization: grid construction, operator coefficients
       (full formulas for every tridiagonal entry)
   (f) Temporal discretization: update equations for each phase/scheme
   (g) Linear system solution procedure
   (h) Any special treatment (fitting factors, smoothing, monitoring
       dates, exercise constraints)
   (i) Output extraction (price, Greeks if covered)

4. CONDITIONS AND CONSTRAINTS: All conditions on parameters
   (grid sizes, time steps, volatility/rate regimes) required for:
   - Well-definedness (matrix invertibility)
   - Positivity preservation
   - Stability
   - Stated convergence order

5. MONOLITHIC PSEUDOCODE: 100–250 numbered lines translating §3
   into algorithmic form. Requirements:
   - Typed variables (float, int, float[])
   - All indices explicit
   - Named constants with default values
   - Comments referencing §3 subsections: "// per §3(e)"
   - Handles all product types the paper addresses via branching
   - Input/output clearly defined

6. ERRATA RELATIVE TO ORIGINAL PAPER: A compact list of every
   correction made, with:
   | Location in paper | Original | Corrected | Reason |
</task>

<constraints>
- This document must be SELF-CONTAINED. A developer should never
  need to consult the original paper.
- Every formula must be explicit — no "see above" or "similarly."
- The specification must be FAITHFUL to the paper's method. Do not
  improve, optimize, or extend it. If the paper uses implicit Euler,
  specify implicit Euler — do not substitute Rannacher-CN.
- If any aspect remains genuinely ambiguous after Rounds 1.1–1.4,
  state the ambiguity and your best-judgment resolution, clearly
  labeled as such.
- The pseudocode is mathematical pseudocode, not code. No language-
  specific syntax.
</constraints>

<output_format>
- Six numbered sections as specified
- LaTeX for all formulas
- Pseudocode in fenced code block with line numbers
- Errata as Markdown table
- Target: 4000–6000 words
</output_format>
```

**Reasoning effort:** `high`

---

## Stage 2: Implementation (4 Rounds)

Stage 2 is largely the same as before but now references the golden reference document from Stage 1 Round 1.5 instead of a novel algorithm design.

### Round 2.1 — Architecture Mapping

```markdown
<context>
The implementation target is the method specified in the golden
reference document from Stage 1.
<golden_reference>
[Paste Round 1.5 output — the self-contained specification]
</golden_reference>

<quantlib_headers>
[Attach QuantLib FD framework headers if targeting QuantLib.
 If standalone implementation, omit this and replace the task
 below with standalone architecture design.]
</quantlib_headers>

If targeting QuantLib: use the NEW framework only (FdmLinearOpComposite,
FdmMesher, FdmBackwardSolver, FdmSchemeDesc, TripleBandLinearOp).
Do NOT use deprecated classes (CrankNicolson<>, MixedScheme<>, etc.).
</context>

<task>
Map the golden reference's components to the implementation target.

1. COMPONENT MAPPING TABLE:
   | Golden Ref Component (§ref) | Target Class/Module | Action |
   Action ∈ {REUSE, EXTEND, COMPOSE, REPLACE, NEW}

2. NEW/EXTENDED COMPONENT REGISTRY: For each non-REUSE:
   - Proposed name
   - Parent class/interface
   - Key method signatures (pseudocode)
   - Closest existing analog

3. DATA FLOW DIAGRAM: ASCII (40+ lines) showing all modules,
   data flow, hot/cold path distinction.

4. FAITHFULNESS CONSTRAINTS: List any points where the target
   framework's conventions differ from the golden reference's
   formulation (time direction, coordinate convention, storage
   layout). For each, specify the exact adaptation needed.
</task>

<constraints>
- The implementation must be FAITHFUL to the golden reference.
  If the framework's architecture pushes toward a different approach,
  flag the tension and propose the minimal adaptation.
- Do NOT write implementation code. Typed pseudocode only.
- Do NOT redesign the target framework's architecture.
</constraints>

<output_format>
- Tables in Markdown, diagram in fenced ASCII block
- Pseudocode signatures in fenced blocks
- Target: 2500–3500 words
</output_format>
```

**Reasoning effort:** `high`

---

### Round 2.2 — Core Module Pseudocode

```markdown
<context>
<prior_output>[Paste Round 2.1 mapping]</prior_output>
<golden_reference>[Paste Round 1.5 specification, especially §3 and §5]</golden_reference>
</context>

<task>
Detailed pseudocode for each module marked NEW or EXTEND in Round 2.1.

For EACH module provide:
  (a) PURPOSE — one sentence
  (b) INTERFACE — function signatures
  (c) INVARIANTS — preconditions, postconditions
  (d) ALGORITHM SKELETON — pseudocode, referencing golden reference
      sections: "// per GoldenRef §3(e)"
  (e) COMPLEXITY — time and space

Priority: the modules containing the paper's core numerical innovation
should receive the most detailed treatment (50–80 lines). Supporting
modules can be briefer (15–40 lines).

Faithfulness check: after each module, confirm that the pseudocode
implements exactly what the golden reference specifies — not a
variation or improvement.
</task>

<constraints>
- Reference golden reference sections instead of repeating formulas.
- Maximum 80 lines per module skeleton.
- C-like pseudocode. No language-specific syntax.
- No boilerplate.
- If the golden reference leaves an implementation choice open
  (e.g., Thomas algorithm vs. general tridiagonal solver), state
  your choice and note it is an implementation decision, not a
  mathematical one.
</constraints>

<output_format>
- One subsection per module, each containing (a)–(e)
- Pseudocode in fenced code blocks
- Target: 2500–3500 words
</output_format>
```

**Reasoning effort:** `high`

---

### Round 2.3 — Orchestrator & Integration

```markdown
<context>
<prior_outputs>[Paste Rounds 2.1 and 2.2]</prior_outputs>
<golden_reference>[Paste Round 1.5, especially §5 pseudocode]</golden_reference>
</context>

<task>
1. ORCHESTRATOR PSEUDOCODE (30–50 lines): Top-level driver wiring
   all modules. Show the complete solve lifecycle from input
   parameters to output price/Greeks.

2. CONFIGURATION TABLE: For each product type the paper addresses:
   | Product | Module settings | Default numerical parameters |

3. FAITHFULNESS CROSS-CHECK: A table mapping each line-range of
   the golden reference pseudocode (§5) to the corresponding module
   and line-range in the modular pseudocode. Confirm 1-to-1 coverage:
   | GoldenRef lines | Module | Module lines | Status |
   Status ∈ {exact match, adapted for framework, implementation choice}
</task>

<constraints>
- The faithfulness cross-check is mandatory. Every line of the golden
  reference pseudocode must map to something in the modular design.
- If any golden reference functionality is lost or changed in
  modularization, flag it explicitly.
</constraints>

<output_format>
- Three sections
- Pseudocode in fenced blocks, tables in Markdown
- Target: 1500–2500 words
</output_format>
```

**Reasoning effort:** `high`

---

### Round 2.4 — Testing & Validation

```markdown
<context>
<prior_outputs>[Summarize Rounds 2.1–2.3 key decisions]</prior_outputs>
<golden_reference>[Paste Round 1.5 §4 conditions and §6 errata]</golden_reference>
</context>

<task>
1. BENCHMARK TEST SUITE: Use the paper's own numerical examples
   as primary test cases. For each example in the paper:
   | Test ID | Paper reference (§/Fig/Table) | Parameters | Paper's reported value | Independent reference (if available) |

   Add 2–3 supplementary tests for edge cases not covered by the
   paper (informed by Round 1.4 limiting-case analysis).

2. CONVERGENCE TEST PROCEDURE:
   - ASCII flowchart for Richardson extrapolation convergence test
   - Expected convergence order (from golden reference §4)
   - PASS/FAIL criteria

3. REGRESSION TEST TABLE:
   | Test ID | Target value | Tolerance | What it validates |

4. DIAGNOSTIC CHECKS (per-run, not per-test):
   - Positivity check: any solution values < 0?
   - M-matrix check: are operator conditions satisfied?
   - Oscillation detector: sign-change count in solution differences
   - Condition number of tridiagonal system (optional)
</task>

<constraints>
- Primary validation: reproduce the paper's own reported results.
  If the implementation matches the paper's results, it is
  faithful. If it does not, the discrepancy must be investigated.
- Tolerances must account for the paper's reported grid parameters.
- If the paper does not report numerical values (only plots),
  note this and specify how to extract approximate values from
  the figures.
</constraints>

<output_format>
- Four sections, tables and flowchart
- Target: 1500–2500 words
</output_format>
```

**Reasoning effort:** `high`

---

## Summary: Previous vs. Revised Framework

| Aspect | Previous (Novel Method) | Revised (Paper Implementation) |
|---|---|---|
| Stage 1 goal | Design a new algorithm | Reconstruct and verify the paper's algorithm |
| Round 1.1 | Critical analysis of 2 papers | Complete extraction from 1 paper |
| Round 1.2 | SOTA literature survey | Gap identification and classification |
| Round 1.3 | Novel algorithm design | Gap resolution using supplementary materials |
| Round 1.4 | Pseudocode | End-to-end verification of corrected method |
| Round 1.5 | Verification | Implementation-ready golden reference |
| Supplementary materials serve to... | Inspire technique combinations | Fill specific identified gaps |
| Mathematical emphasis | Proof sketches for novel claims | Forensic audit of existing claims |
| Faithfulness constraint | N/A | Every round checks fidelity to the paper |
| Primary validation | Convergence to known solutions | Reproduce the paper's own reported results |
| xhigh rounds | Algorithm design (1.3) | Gap identification (1.2), gap resolution (1.3), verification (1.4) |

The key architectural difference: three rounds now operate at `xhigh` reasoning effort (1.2, 1.3, 1.4) instead of one, because the forensic work of finding gaps, filling them correctly, and verifying the complete chain demands deeper reasoning than designing something new — you're constrained to match the paper exactly while ensuring every hidden step is sound.

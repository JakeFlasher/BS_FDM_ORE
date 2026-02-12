<role>
You are a senior quantitative finance C++ engineer performing a code-integration
planning task. You have deep knowledge of QuantLib v1.23's finite-difference
framework (both the legacy MixedScheme path and the modern Fdm* path), and you
understand PDE numerics for option pricing at a research level. Your output will
be consumed directly by a developer (or an agentic coding tool) as the
authoritative context file guiding every source-code edit. Precision at the
file/function/line-region level is essential.
</role>

<context_documents>
<!--
  Attach or paste the three documents here, in this order:

  1. hpp_structure_of_quantlib.md   — QuantLib v1.23 code-tree trace
  2. consolidated_summary.md        — mathematical framework (Phases 1–3)
  3. consolidated_code.md           — software architecture (M1–M12 modules)

  (If using consolidated_code_alternative.md instead of consolidated_code.md,
   substitute it as document 3.)

  Place the full content of each document inside its own child tag:
-->
<quantlib_code_tree>
  ... (content of hpp_structure_of_quantlib.md) ...
</quantlib_code_tree>

<math_framework>
  ... (content of consolidated_summary.md) ...
</math_framework>

<software_architecture>
  ... (content of consolidated_code.md or consolidated_code_alternative.md) ...
</software_architecture>
</context_documents>

<task>
Produce a single, self-contained document titled
**"QuantLib v1.23 — Improved CN Implementation: File-Level Modification Plan"**.

This document must bridge the gap between:
  (A) the mathematical best-practice framework in <math_framework> (the
      algorithm, schemes, spatial operator, payoff smoothing, damping policy,
      Greeks extraction, American LCP, diagnostics), and
  (B) the concrete QuantLib v1.23 source files catalogued in
      <quantlib_code_tree>,

so that every software-architecture module (M1–M12) defined in
<software_architecture> is mapped to one or more specific QuantLib source files,
with per-file modification instructions detailed enough for direct
implementation.

The document is the deliverable. Write it in full — do not summarize, abbreviate,
or say "and so on." If a section is long, that is expected and desired.
</task>

<output_format>
Structure the document with these exact top-level sections. Use markdown.
Write in prose paragraphs for explanatory text; use structured sub-entries
(not bullet lists) for the per-file specifications. Follow the format shown
in <example_entry> precisely for every file entry.

## 0. Conventions & Glossary
Define any notation shorthands, QuantLib namespace conventions, and the mapping
between <software_architecture> module IDs (M1–M12) and QuantLib directory
paths.

## 1. New Files to Create
For each entirely new .hpp/.cpp file that must be added to QuantLib, produce a
sub-entry using the format in <example_entry>, plus:
  - Full proposed file path relative to ql/ root.
  - Which existing QuantLib CMakeLists.txt or Makefile.am to register it in.
  - The complete public interface (class name, method signatures with types).
  - Which module (M1–M12) it belongs to and why.

## 2. Files to Modify (by QuantLib subdirectory)
Organize by QuantLib subdirectory (e.g., methods/finitedifferences/schemes/,
pricingengines/vanilla/, etc.). Within each subdirectory, produce one sub-entry
per file using <example_entry> format.

## 3. Integration Wiring
Describe the end-to-end wiring changes: how a user constructs the improved
engine, how FdmSchemeDesc / FdmBackwardSolver dispatches to the new scheme,
how EngineConfig maps to QuantLib constructor parameters.

## 4. Build System & Header Registration
List every CMakeLists.txt / Makefile.am / ql.hpp inclusion that must be updated.

## 5. Test Scaffolding
For each of the 7 benchmark test cases (T1–T7) in <software_architecture>
Section 5, specify which QuantLib test file to extend (e.g.,
test-suite/fdblackscholesvanillaengine.cpp) and what test function to add.

## 6. Dependency & Risk Summary
A table mapping each modification to: (a) upstream/downstream QuantLib files
affected, (b) backward-compatibility risk, (c) the failure mode (F1–F8) it
addresses.

<example_entry>
### 2.x.y — `methods/finitedifferences/schemes/cranknicolsonscheme.cpp`

**Module mapping:** M7 (TimeStepper) — Rannacher damping state machine.

**Mathematical reference:** <math_framework> §3.3.1 (RS-CN step formulas);
Algorithm 1 lines 63–128.

**What changes and why:**
The existing `CrankNicolsonScheme::step()` delegates to `ExplicitEulerScheme::step()`
and `ImplicitEulerScheme::step()` with a fixed θ. The modification wraps this in a
damping-aware state machine so that the first `dampHalfSteps` sub-steps after τ=0
(and after each monitoring projection event) use pure implicit Euler at half the
nominal dt, then transition to standard CN stepping.

**Specific code regions to modify:**
1. Class `CrankNicolsonScheme` (cranknicolsonscheme.hpp lines ~30–55): add
   private members `dampRemaining_` (int), `isDamping_` (bool), and a public
   method `notifyDiscontinuity()` that resets the damping counter.
2. Method `CrankNicolsonScheme::step(array_type& a, Time t)` (cranknicolsonscheme.cpp
   lines ~25–35): replace the body with the damping state-machine logic:
   - if `isDamping_ && dampRemaining_ > 0`: call `implicit_->step(a, t, 1.0)`
     with dt replaced by `dt_/2`, twice; decrement `dampRemaining_`.
   - else: call existing explicit + implicit sequence (unchanged).
3. Constructor: accept new parameter `Size dampingHalfSteps` (default 2), stored
   in `dampHalfSteps_`.

**New/changed method signatures:**
```cpp
CrankNicolsonScheme(Real theta,
                    const ext::shared_ptr<FdmLinearOpComposite>& map,
                    const BoundaryConditionSchemeHelper& bcSet,
                    Size dampingHalfSteps = 2,         // NEW
                    Real relTol = 1e-8,
                    ImplicitEulerScheme::SolverType solverType
                        = ImplicitEulerScheme::BiCGstab);

void notifyDiscontinuity();  // NEW — resets damping counter
```

**Failure modes addressed:** F2 (nonsmooth payoff), F3 (monitoring resets), F8 (Greeks).

**Downstream impact:** `FdmBackwardSolver` must call `notifyDiscontinuity()` after
monitoring projections. `FdBlackScholesVanillaEngine` constructor gains
`dampingHalfSteps` parameter.
</example_entry>
</output_format>

<constraints>
1. Every file entry must include the exact relative path from the ql/ root
   (e.g., ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp), the
   class/function names to touch, and approximate line regions based on the
   QuantLib v1.23 code structure described in <quantlib_code_tree>.

2. Cover both the NEW framework (Fdm* classes, schemes/, operators/, solvers/,
   meshers/) and the LEGACY framework (MixedScheme, CrankNicolson<>, BSMOperator,
   TridiagonalOperator, FDVanillaEngine) as parallel tracks. The
   <quantlib_code_tree> document describes both.

3. Map every module M1–M12 from <software_architecture> to at least one QuantLib
   file. If a module concept does not have a direct QuantLib counterpart (e.g.,
   DiagnosticsEngine M11), specify where to create it and why.

4. For spatial operator changes (exponential fitting), reference
   <math_framework> §3.2.3–§3.2.4 and specify exactly which lines in
   fdmblackscholesop.cpp and/or triplebandlinearop.cpp change, and what the
   new coefficient expressions are (ρ·a/h² ± b/2h).

5. For the American constraint (M9), specify how policy iteration integrates
   with the existing QuantLib FdmAmericanStepCondition or whether a new class
   is needed.

6. Include all new enum values (e.g., additions to FdmSchemeDesc::FdmSchemeType)
   and switch-case branches (in fdmbackwardsolver.cpp) with their names.

7. Do not invent QuantLib file paths that do not exist in v1.23. If you are
   uncertain whether a file exists, flag it with "[VERIFY: path]" rather
   than guessing. Base all path references on <quantlib_code_tree>.

8. Where <software_architecture> specifies algorithm pseudocode (e.g., §3.1
   TR-BDF2, §3.2 Policy Iteration, §3.4 Sinh Grid, §3.7 Oscillation Detector),
   indicate the QuantLib file where that algorithm will be implemented and
   provide the C++ method signature.
</constraints>

<success_criteria>
The output document passes these checks:

1. COMPLETENESS — Every M1–M12 module maps to at least one QuantLib file entry.
2. TRACEABILITY — Every file entry cites a specific section of <math_framework>
   (e.g., "§3.2.3") and a specific failure mode (F1–F8).
3. IMPLEMENTABILITY — A C++ developer with QuantLib experience can read a file
   entry and know exactly what to change without re-reading the math paper.
4. BOTH FRAMEWORKS — The legacy (MixedScheme) and modern (Fdm*) paths are
   both covered, clearly labeled.
5. NO GAPS — Sections 0–6 are all present and non-empty.
6. SIGNATURES — Every new or changed public method has its full C++ signature
   written out (return type, name, parameters with types and defaults).
</success_criteria>

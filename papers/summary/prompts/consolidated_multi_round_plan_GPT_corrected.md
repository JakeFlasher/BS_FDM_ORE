## What I changed (rigorous prompt audit + fixes)

Your current `multi-round_prompt.md` is already strong (clear round structure, explicit constraints, QuantLib conventions). The main issues are **internal contradictions** and a few **underspecified / fragile** points that can cause Opus to implement the “wrong thing” even when it’s trying to follow instructions precisely.

Here are the highest-impact problems I found and how the replacement prompt fixes them:

1. **Log-coordinate inconsistency (x = ln(S/K) vs ln(S))**  
   - The shared math preamble says $$x=\ln(S/K)$$, but later steps (notably the barrier projection condition design) implicitly assume the mesher stores $$\ln(S)$$ because the condition signature has **no strike** parameter to convert barriers in S-space into log-moneyness.  
   - The fix: the new framework **requires an explicit, code-based verification** of what `FdmBlackScholesMesher` and `FdmBlackScholesOp` mean by “x”, and it provides a clean branch depending on what the source says—no guessing.

2. **Uniform-grid formulas used while adding a nonuniform sinh mesher**  
   - Exponential fitting in your preamble uses a single uniform $$h$$, but Round 1 introduces a nonuniform sinh mesher. Without a defined “local h”, Opus could pick an arbitrary spacing and silently break the intended monotonicity behavior.  
   - The fix: the new prompt defines **exactly how to compute a local effective spacing** from the mesher (`dplus/dminus`), and it explicitly states what properties are guaranteed (uniform-grid proof) vs what must be validated (nonuniform).

3. **Sinh mesher alignment step order is logically inconsistent**  
   - Your algorithm says “force endpoints” and then “shift entire grid” to align targets; shifting after forcing endpoints breaks the endpoint constraints.  
   - The fix: the replacement prompt changes the order and makes multi-target alignment behavior explicit (single best shift; then re-pin endpoints).

4. **Monitoring restart damping: Round 3 and Round 4 don’t agree**  
   - Round 3 adds `CrankNicolsonScheme::notifyDiscontinuity()` with **two implicit Euler half-steps** inside `step()`.  
   - Round 4 then chooses a fallback approach that uses **full implicit Euler steps** via `FdmBackwardSolver`, making the Round 3 restart logic effectively unused and not matching the Rannacher “two half-steps” spec.  
   - The fix: the replacement prompt standardizes on **one coherent approach**: restart smoothing is handled via the CN scheme’s `notifyDiscontinuity()` + segmenting at stopping times, so the “two half-steps” behavior is actually used.

5. **Missing/ambiguous source-file facts are handled informally**  
   - Anthropic’s official guidance for avoiding hallucinations recommends explicitly allowing “I don’t know,” restricting the model to provided documents, and using grounding/verification steps. ([docs.anthropic.com](https://docs.anthropic.com/en/docs/test-and-evaluate/strengthen-guardrails/reduce-hallucinations))  
   - The fix: the replacement prompt includes a **hard stop rule**: if a placeholder file isn’t provided, Opus must ask for it rather than reconstructing/hallucinating.

6. **Prompt structure not optimized for long-context performance**  
   - Anthropic recommends putting long documents **at the top** of the prompt and placing the query/instructions at the end. ([docs.anthropic.com](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/long-context-tips))  
   - The fix: each round template is reorganized so `<source_files>` come first, and the task request comes last.

All of the above improvements follow Anthropic’s official “prompt engineering overview” sequencing (clarity → examples → CoT → XML → role prompting → chaining → long context). ([docs.anthropic.com](https://docs.anthropic.com/en/docs/prompt-engineering))  
They also follow the Claude 4/4.6 best-practices emphasis on explicit requirements + context and being careful with examples/details. ([docs.anthropic.com](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/claude-4-best-practices))  
And they use XML-tag structuring as recommended. ([platform.claude.com](https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/use-xml-tags))  
Finally, they use role prompting in the system prompt as recommended. ([docs.anthropic.com](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/system-prompts))  

---

# Replacement framework: **Opus 4.6 multi-round prompt (QuantLib v1.23 Improved CN)**

Below is a **drop-in replacement** for your `multi-round_prompt.md`, designed for **Claude Opus 4.6**.

It is written as a **two-part template**:

1. **SYSTEM PROMPT** (stable across all rounds)  
2. **USER PROMPT** per round (Round 0 optional + Rounds 1–7 + Optional legacy Round 8)

It is formatted using XML tags, with the **source files placed first** in each round, per Anthropic’s long-context guidance. ([docs.anthropic.com](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/long-context-tips))  

Copy/paste as-is, then fill the `[INSERT FULL CONTENT ...]` placeholders.

---

## SYSTEM PROMPT (use as the Messages API `system` parameter)

```text
You are an expert C++ developer working on QuantLib v1.23.
Your task is to implement a staged set of improvements to QuantLib’s modern finite difference (Fdm*) framework for 1D Black-Scholes option pricing.

You MUST follow these constraints exactly:

(1) SOURCE-OF-TRUTH / NO-HALLUCINATION RULES (mandatory)
- Use ONLY the provided source files and prior-round outputs as factual reference for:
  * class/method signatures
  * include paths
  * base-class requirements
  * coding conventions and header styles
- If any referenced file content is missing, DO NOT guess or reconstruct.
  Stop and ask the user to provide the missing file(s).
- If a requirement is ambiguous (e.g., coordinate definition in mesher), you MUST:
  (a) search the provided files for the relevant implementation detail, and
  (b) implement the behavior that matches the source.
  If still ambiguous, stop and ask a precise clarification question.
- You are explicitly allowed to say “I don’t know” / “I need file X” when you lack source data.

(2) CODING CONVENTIONS (QuantLib style; mandatory)
- Namespace: all classes in namespace QuantLib { }
- Smart pointers: use ext::shared_ptr and ext::make_shared (NOT std::shared_ptr)
- Arrays: use QuantLib::Array; respect existing signatures (use Disposable<Array> where that’s the surrounding convention)
- Preconditions: QL_REQUIRE(condition, message)
- Postconditions: QL_ENSURE(condition, message)
- Headers:
  * include guards: #ifndef quantlib_<path_underscored>_hpp / #define ...
  * no “using namespace std;” in headers
  * includes: <ql/...> for QuantLib headers
  * match copyright header style used by adjacent files in the same directory
- Naming: CamelCase for classes, camelCase for methods, camelCase_ for private data members

(3) OUTPUT CONTRACT (mandatory)
- Output COMPLETE files, never diffs.
- Output ONLY the files requested for that round (no extras).
- Each file must compile against unmodified QuantLib v1.23 headers plus any prior-round outputs.

(4) INTERNAL VERIFICATION PROTOCOL (do before coding; do not print)
Before writing code each round, internally confirm:
- What “x” means in the relevant mesher/operator (ln(S) vs ln(S/K) etc.) by reading the provided files.
- What time variable convention the solver uses (rollback direction, t1/t2 meanings).
- Which methods are pure virtual and must be implemented.
- Which helper methods/fields are public vs protected (avoid forbidden access).

Do your reasoning step-by-step privately, but in the final output you MUST follow the round’s output format exactly.
```

---

## USER PROMPT TEMPLATES (one per round)

### Formatting note (important)
For each round, put **all large source files first**, then the task at the end. This follows Anthropic’s long-context guidance. ([docs.anthropic.com](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/long-context-tips))  

### Output format (used in every coding round)
Every file must be output as:

```text
FILE: <relative/path/from/ql/root>
```cpp
<full file contents>
```
```

No commentary before, between, or after files.

---

## ROUND 0 (OPTIONAL BUT RECOMMENDED) — Preflight fact-check questionnaire

Use this if you want Opus to *only* audit the provided QuantLib v1.23 sources and your plan before any code is generated.

```xml
<round id="0">
  <source_files>
    <!-- Put the REAL contents here -->
    <document path="ql/methods/finitedifferences/operators/fdmblackscholesop.hpp">
      [INSERT FULL CONTENT]
    </document>
    <document path="ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp">
      [INSERT FULL CONTENT]
    </document>
    <document path="ql/methods/finitedifferences/finitedifferencemodel.hpp">
      [INSERT FULL CONTENT]
    </document>
    <document path="ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp">
      [INSERT FULL CONTENT]
    </document>
    <document path="ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp">
      [INSERT FULL CONTENT]
    </document>
  </source_files>

  <task>
    You are in preflight mode. Do NOT write code.

    Answer these questions with citations to specific lines/identifiers in the provided sources:
    1) What does the 1D BS mesher store as spatial coordinate x? ln(S) or ln(S/K) or something else?
       - Show where (which constructor / method) it is defined.
    2) How do FirstDerivativeOp / SecondDerivativeOp incorporate nonuniform mesh spacing?
       - Identify which mesher spacing functions they use.
    3) Does FiniteDifferenceModel call StepCondition::applyTo at every time step, or only at stopping times?
    4) Does FiniteDifferenceModel copy the evolver by value?
    5) Are FdmSchemeDesc fields const or non-const in v1.23?
    6) Confirm the correct sign expectations for “M-matrix satisfied” checks in TripleBandLinearOp for the BS operator.
       (i.e., which diagonals should be nonnegative/nonpositive for the operator itself.)

    If any answer cannot be determined from the provided sources, say exactly what additional file is needed.
  </task>

  <output_spec>
    Output a single <preflight_report> block (no code).
  </output_spec>
</round>
```

---

## ROUND 1 — New fitted operator + sinh mesher + diagnostics

> This round creates **new files only** (no modifications).  
> The design is intentionally close to `FdmBlackScholesOp`.

```xml
<round id="1">
  <source_files>
    <!-- Put the full contents of these existing QuantLib headers here -->
    <document path="ql/methods/finitedifferences/operators/fdmblackscholesop.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/operators/triplebandlinearop.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/operators/firstderivativeop.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/operators/secondderivativeop.hpp">[INSERT FULL CONTENT]</document>

    <document path="ql/methods/finitedifferences/meshers/fdm1dmesher.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/meshers/concentrating1dmesher.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/meshers/uniform1dmesher.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/meshers/fdmmesher.hpp">[INSERT FULL CONTENT]</document>

    <document path="ql/math/array.hpp">[INSERT FULL CONTENT OR DECLARATION]</document>
  </source_files>

  <task>
    ROUND 1 OF 7 — Create 3 new standalone class pairs (.hpp + .cpp).
    Output exactly 6 files (3 headers + 3 implementations), in the specified order.

    IMPORTANT FACT-CHECKS (must do before coding; do not print):
    - Confirm from the provided sources whether the spatial coordinate stored in the mesher is ln(S) or ln(S/K).
      Your fitted-operator implementation MUST use the same convention as FdmBlackScholesOp.
    - Confirm how to obtain local mesh spacing along the chosen direction:
      prefer mesher->dplus/dminus (or equivalent) if available; otherwise derive effective h from the 1D mesher.

    FILE PAIR 1: FdmFittedBlackScholesOp
    Path:
      ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp
      ql/methods/finitedifferences/operators/fdmfittedblackscholesop.cpp

    Goal:
    - Implement a 1D Black-Scholes operator with exponential fitting for convection dominance.
    - Keep member layout parallel to FdmBlackScholesOp:
        FirstDerivativeOp dxMap_;
        SecondDerivativeOp dxxMap_;
        TripleBandLinearOp mapT_;

    Key algorithm (per node i):
    - sigma_i from volTS_ or localVol_ (mirror FdmBlackScholesOp logic)
    - a_i = 0.5*sigma_i^2
    - b_i = (r - q) - a_i   (log-space convection coefficient; include quanto adjustments exactly as in FdmBlackScholesOp)
    - Compute local effective spacing h_i:
        h_i = 0.5*(dplus_i + dminus_i)   (must use mesher spacing consistent with derivative ops)
      Guard: if h_i <= 0, QL_REQUIRE fail.
    - theta_i = b_i * h_i / (2*a_i)  (clamp a_i to a small positive value, e.g. 1e-20, before dividing)
    - rho_i = fittingFactor(theta_i):
        if |theta| < 1e-8: rho = 1 + theta^2/3
        else:              rho = theta / tanh(theta)   (theta * coth(theta))
    - Store:
        convection[i] = b_i
        fittedDiffusion[i] = a_i * rho_i
    - Assemble operator using TripleBandLinearOp::axpyb():
        mapT_.axpyb(convection, dxMap_,
                    dxxMap_.mult(fittedDiffusion),
                    Array(1, -r_effective));
      where r_effective matches FdmBlackScholesOp (incl. quanto helper if used there).

    Diagnostics:
    - mMatrixSatisfied(): true iff all off-diagonals of mapT_ satisfy the required sign pattern.
      (You MUST confirm the sign convention used by TripleBandLinearOp in the provided code;
       do not assume without checking.)
    - mMatrixViolationCount(): count of violations (should be 0).

    Constraints:
    - Do NOT access protected TripleBandLinearOp members.
    - Use axpyb() to populate mapT_.
    - No per-node upwind switching; only clamp a_i for degeneracy.

    FILE PAIR 2: FdmSinhConcentrating1dMesher
    Path:
      ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.hpp
      ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.cpp

    Inherit from Fdm1dMesher.

    Constructor:
      FdmSinhConcentrating1dMesher(
          Real xMin, Real xMax, Size size,
          Real xCenter,
          Real alpha = 3.0,
          const std::vector<Real>& alignTargets = std::vector<Real>());

    Requirements:
    - Handle alpha == 0 as uniform grid in [xMin, xMax].
    - For alpha > 0, implement:
        x(ξ) = xCenter + c*sinh(alpha*(ξ - ξ0)), ξ ∈ [0,1]
      Solve for ξ0 by bisection for asymmetric cases (xCenter not midpoint), and compute c accordingly.
      Make the bisection target function explicit in comments.
    - Alignment targets:
      * Choose at most one global shift (the smallest |shift| across targets) such that a nearest node snaps to a target
        IF |shift| < 0.5 * minSpacing.
      * Apply the shift to all nodes.
      * After shifting, re-pin endpoints exactly: x[0]=xMin, x[last]=xMax.
      * Do not attempt multiple sequential shifts that can undo previous alignment.
    - Populate locations_, dplus_, dminus_ following the same pattern used by Concentrating1dMesher.

    FILE PAIR 3: FdmDiagnostics
    Path:
      ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp
      ql/methods/finitedifferences/utilities/fdmdiagnostics.cpp

    Implement:
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
          static FdmDiagnosticsReport merge(const FdmDiagnosticsReport& a,
                                            const FdmDiagnosticsReport& b);
          Level level() const;
        private:
          Level level_;
      };

    oscillationScore:
    - O(N), allocation-free.
    - Count sign changes in du_j = u[j+1]-u[j], ignoring |du| < 1e-15.
    - Normalize by max(1, size-2).

    Thread safety:
    - static methods must not use shared mutable state.
  </task>

  <constraints>
    - Do NOT modify any existing QuantLib files in this round.
    - Each .cpp includes its own header first.
    - Use QuantLib conventions (ext::shared_ptr, QL_REQUIRE, etc.).
  </constraints>

  <output_spec>
    Output exactly 6 files in this order:
    1) ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp
    2) ql/methods/finitedifferences/operators/fdmfittedblackscholesop.cpp
    3) ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.hpp
    4) ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.cpp
    5) ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp
    6) ql/methods/finitedifferences/utilities/fdmdiagnostics.cpp
  </output_spec>
</round>
```

---

## ROUND 2 — Barrier projection step condition + policy iteration LCP

This round is mostly fine in your original prompt; the big fix is **coordinate conversion must be verified** from source.

```xml
<round id="2">
  <source_files>
    <document path="ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/stepconditions/fdmamericanstepcondition.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/pde/stepcondition.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/meshers/fdmmesher.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/operators/triplebandlinearop.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/math/array.hpp">[INSERT FULL CONTENT OR DECLARATION]</document>
  </source_files>

  <task>
    ROUND 2 OF 7 — Create 2 new class pairs (.hpp + .cpp). Output exactly 4 files.

    FILE PAIR 1: FdmBarrierProjectionCondition
    Path:
      ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp
      ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.cpp

    Inherits: StepCondition<Array>.

    Constructor:
      FdmBarrierProjectionCondition(
          std::vector<Time> monitoringTimes,
          Real lowerBarrier, Real upperBarrier,   // in S-space
          ext::shared_ptr<FdmMesher> mesher,
          Size direction = 0);

    Precompute outsideIndices_ ONCE in the constructor:
    - Sort monitoringTimes and store.
    - Convert barriers to the same coordinate space used by mesher->location(iter, direction):
      IMPORTANT: You MUST verify from provided sources whether that location is ln(S) or ln(S/K).
      If it’s ln(S): compare directly to ln(barrier).
      If it’s ln(S/K): you must determine how K is defined/accessed; if K is not accessible,
        STOP and ask for the relevant source file or adjust design (do not guess).

    Edge cases:
    - monitoringTimes empty: applyTo is always no-op.
    - lowerBarrier == 0: treat as no lower barrier.
    - upperBarrier == +infinity: treat as no upper barrier.

    applyTo(Array& a, Time t) const:
    - If t matches a monitoring time within tolerance 1e-10, set a[i]=0 for all outsideIndices_.

    Must provide:
    - const std::vector<Time>& monitoringTimes() const;

    FILE PAIR 2: FdmPolicyIterationLCP
    Path:
      ql/methods/finitedifferences/utilities/fdmpolicyiteration.hpp
      ql/methods/finitedifferences/utilities/fdmpolicyiteration.cpp

    Implement the policy-iteration solver operating on plain Arrays:
    - Solve LCP for American constraint using an internal Thomas solver on copies.
    - Handle n=1.
    - solve() is const; lastIterations_ is mutable.

    Keep algorithm as described in the original framework (active set, modify rows to identity, solve, converge).
  </task>

  <output_spec>
    Output exactly 4 files:
    1) ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp
    2) ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.cpp
    3) ql/methods/finitedifferences/utilities/fdmpolicyiteration.hpp
    4) ql/methods/finitedifferences/utilities/fdmpolicyiteration.cpp
  </output_spec>
</round>
```

---

## ROUND 3 — Scheme modifications (CN restart damping + SchemeDesc field)

This round’s content is kept, but the replacement framework makes Round 4 *actually use it*.

```xml
<round id="3">
  <source_files>
    <document path="ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/schemes/cranknicolsonscheme.cpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/schemes/impliciteulerscheme.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/schemes/expliciteulerscheme.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/schemes/boundaryconditionschemehelper.hpp">[INSERT FULL CONTENT]</document>
  </source_files>

  <task>
    ROUND 3 OF 7 — Modify FdmSchemeDesc and CrankNicolsonScheme.

    FACT-CHECK (must do before coding; do not print):
    - Determine whether FdmSchemeDesc members are const in v1.23 and choose the correct implementation strategy.

    FILE 1: ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
    - Add Size monitoringDampingSteps to FdmSchemeDesc (default 0).
    - Update factories to set monitoringDampingSteps=0.
    - Add new factory:
        CrankNicolsonWithMonitoringDamping(Size monitoringDampingSteps=2)
      returning {CrankNicolsonType, 0.5, 0.0, monitoringDampingSteps}.

    CRITICAL: The FdmBackwardSolver class declaration/body must remain byte-for-byte identical in this round.

    FILE 2-3: cranknicolsonscheme.hpp/.cpp
    - Add constructor parameter Size dampingHalfSteps = 0.
    - Add notifyDiscontinuity() and isDamping() const.
    - Implement damping as two implicit Euler half-steps at dt/2 when in damping phase.

    Additional requirements (added for robustness):
    - QL_REQUIRE(dampingHalfSteps % 2 == 0, "dampingHalfSteps must be even");
      (0 allowed)
  </task>

  <output_spec>
    Output exactly 3 files:
    1) ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
    2) ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp
    3) ql/methods/finitedifferences/schemes/cranknicolsonscheme.cpp
  </output_spec>
</round>
```

---

## ROUND 4 — Backward solver restart damping + fitted operator selection

Key fix here: **do not revert to “full-step implicit Euler damping segments”**; instead, **use the CN scheme’s notifyDiscontinuity** (from Round 3) so you actually get two half-steps at the segment start.

```xml
<round id="4">
  <source_files>
    <document path="ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp">[INSERT ROUND-3 OUTPUT VERSION]</document>
    <document path="ql/methods/finitedifferences/solvers/fdmbackwardsolver.cpp">[INSERT FULL CONTENT]</document>

    <document path="ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp">[INSERT FULL CONTENT]</document>

    <document path="ql/methods/finitedifferences/finitedifferencemodel.hpp">[INSERT FULL CONTENT]</document>

    <document path="ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp">[INSERT ROUND-3 OUTPUT VERSION]</document>
    <document path="ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp">[INSERT ROUND-1 OUTPUT HEADER]</document>
  </source_files>

  <task>
    ROUND 4 OF 7 — Modify:
    1) FdmBackwardSolver to support restart damping after stopping times (monitoring/discrete events),
       using CrankNicolsonScheme::notifyDiscontinuity() + segment splitting.
    2) FdmBlackScholesSolver to allow choosing FdmFittedBlackScholesOp.

    FACT-CHECKS (must do before coding; do not print):
    - Read FiniteDifferenceModel rollback logic to confirm:
      * how stoppingTimes are handled,
      * whether evolver is copied by value,
      * and whether StepCondition::applyTo can be applied at segment boundaries without double-applying discrete dividend events.
    - Confirm the time ordering conventions (from > to).

    A) FdmBackwardSolver changes:
    - Preserve original rollback behavior EXACTLY when:
        schemeDesc_.monitoringDampingSteps == 0
        OR schemeDesc_.type != CrankNicolsonType
    - When enabled (CN only):
      * Keep the existing initial dampingSteps behavior unchanged.
      * During the CN phase, split the time interval at condition_->stoppingTimes() that lie within the CN phase interval.
      * For each segment AFTER a stopping time, construct a CrankNicolsonScheme with dampingHalfSteps = schemeDesc_.monitoringDampingSteps,
        call notifyDiscontinuity() BEFORE passing it into FiniteDifferenceModel, then rollback that segment.
        (This avoids needing evolver state to persist across models and ensures 2 implicit Euler half-steps at dt/2.)

    Step allocation:
    - Distribute remaining steps across segments proportionally to segment length, ensuring:
      * sum(segmentSteps) == stepsRemaining
      * each segment with positive length has at least 1 step
      * if this is impossible, QL_REQUIRE with a clear message.

    B) FdmBlackScholesSolver changes:
    - Add bool useFittedOperator = false parameter (default false).
    - In performCalculations, instantiate either FdmBlackScholesOp or FdmFittedBlackScholesOp.
    - Include fdmfittedblackscholesop.hpp ONLY in the .cpp (not the header).
  </task>

  <output_spec>
    Output exactly 4 files:
    1) ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp
    2) ql/methods/finitedifferences/solvers/fdmbackwardsolver.cpp
    3) ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp
    4) ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp
  </output_spec>
</round>
```

---

## ROUND 5 — Wiring: composite factory + sinh mesher factory + all.hpp includes

(Kept similar; add dedup/sort requirements explicitly.)

```xml
<round id="5">
  <source_files>
    <document path="ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.cpp">[INSERT FULL CONTENT]</document>

    <document path="ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/meshers/fdmblackscholesmesher.cpp">[INSERT FULL CONTENT]</document>

    <document path="ql/methods/finitedifferences/operators/all.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/meshers/all.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/stepconditions/all.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/methods/finitedifferences/utilities/all.hpp">[INSERT FULL CONTENT]</document>

    <document path="ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp">[INSERT ROUND-1 OUTPUT HEADER]</document>
    <document path="ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.hpp">[INSERT ROUND-1 OUTPUT HEADER]</document>
    <document path="ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp">[INSERT ROUND-1 OUTPUT HEADER]</document>

    <document path="ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp">[INSERT ROUND-2 OUTPUT HEADER]</document>
    <document path="ql/methods/finitedifferences/utilities/fdmpolicyiteration.hpp">[INSERT ROUND-2 OUTPUT HEADER]</document>
  </source_files>

  <task>
    ROUND 5 OF 7 — Wire new components.

    1) Modify FdmStepConditionComposite:
       Add barrierMonitoredComposite(...) factory.
       - If monitoringDates is empty, return vanillaComposite unchanged.
       - Convert monitoringDates to Times and add to stopping times.
       - Ensure stopping times are sorted + unique (with tolerance).
       - Add FdmBarrierProjectionCondition to the condition list.

    2) Modify FdmBlackScholesMesher:
       Add createSinhMesher(...) static method.
       - Compute xMin/xMax exactly like existing heuristic code.
       - Apply constraints if provided.
       - Convert barrierAlignTargets (S-space) to log-space before passing to FdmSinhConcentrating1dMesher.

    3) Update all.hpp headers, alphabetically:
       operators/all.hpp: include fdmfittedblackscholesop.hpp
       meshers/all.hpp:   include fdmsinhconcentrating1dmesher.hpp
       stepconditions/all.hpp: include fdmbarrierprojectioncondition.hpp
       utilities/all.hpp: include fdmdiagnostics.hpp and fdmpolicyiteration.hpp
  </task>

  <output_spec>
    Output exactly these files (only those that change):
    1) ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp
    2) ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.cpp
    3) ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp
    4) ql/methods/finitedifferences/meshers/fdmblackscholesmesher.cpp
    5) ql/methods/finitedifferences/operators/all.hpp
    6) ql/methods/finitedifferences/meshers/all.hpp
    7) ql/methods/finitedifferences/stepconditions/all.hpp
    8) ql/methods/finitedifferences/utilities/all.hpp
  </output_spec>
</round>
```

---

## ROUND 6 — Pricing engine modifications (expose fitted operator + sinh mesher)

(Kept, but with a stronger “stop if missing” rule.)

```xml
<round id="6">
  <source_files>
    <document path="ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp">[INSERT FULL CONTENT]</document>
    <document path="ql/pricingengines/vanilla/fdblackscholesvanillaengine.cpp">[INSERT FULL CONTENT]</document>

    <document path="ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp">[INSERT ROUND-4 OUTPUT HEADER]</document>
    <document path="ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp">[INSERT ROUND-5 OUTPUT HEADER]</document>
  </source_files>

  <task>
    ROUND 6 OF 7 — Modify FD Black-Scholes vanilla engine to expose:
    - useFittedOperator (default false)
    - sinhAlpha (default 0.0 meaning original mesher)

    HARD REQUIREMENT:
    - If the engine file path differs in your tree, STOP and ask for the correct path + contents.
      Do NOT reconstruct from memory.

    Implement:
    - New constructor parameters with backward-compatible defaults.
    - Mesher selection in calculate():
        if sinhAlpha_ > 0: use FdmBlackScholesMesher::createSinhMesher(...)
        else: original behavior unchanged
    - Pass useFittedOperator_ through to FdmBlackScholesSolver

    Include:
    - Add a compilable usage example comment in the header (as in your original spec).
    - Keep new operator header includes in .cpp only, not the engine header.
  </task>

  <output_spec>
    Output exactly 2 files:
    1) ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp
    2) ql/pricingengines/vanilla/fdblackscholesvanillaengine.cpp
  </output_spec>
</round>
```

---

## ROUND 7 — Test suite

(Kept as-is, but add one more fact-check: confirm the test framework macros.)

```xml
<round id="7">
  <source_files>
    <!-- Provide ALL new headers from Rounds 1–6 plus needed QuantLib test infra -->
    <document path="test-suite/[SOME_EXISTING_TEST_FILE].cpp">[INSERT FULL CONTENT OF A SIMILAR TEST]</document>

    <!-- plus all new headers (Round outputs) -->
    [INSERT ALL OUTPUT HEADERS FROM ROUNDS 1–6]
  </source_files>

  <task>
    ROUND 7 OF 7 — Create a comprehensive test file for all improvements.

    FACT-CHECK:
    - Inspect the provided existing test file(s) to confirm:
      * which framework is used (Boost.Test vs custom macros),
      * naming conventions,
      * how tests are registered.

    Implement test-suite/fdmimprovedcn.cpp with tests T1–T7 (as in your original spec),
    including:
    - Richardson order check
    - positivity/oscillation diagnostics via FdmDiagnostics
    - fitted operator m-matrix diagnostic
    - discrete monitoring restart damping
    - American policy iteration complementarity

    Performance constraint:
    - Keep runtime < 10 seconds per test on modern hardware by choosing grids carefully.
  </task>

  <output_spec>
    Output exactly 1 file:
      test-suite/fdmimprovedcn.cpp
  </output_spec>
</round>
```

---

## OPTIONAL ROUND 8 — Legacy framework integration (only if you truly need it)

This keeps the legacy work clearly isolated, since your original system context says “not legacy” but your later plan mentions it. Keep it optional to avoid scope creep.

---

# Why this redesigned prompt should work better with Opus 4.6 (grounded in Anthropic guidance)

- **System prompt sets the role; user prompts contain the task** (role prompting best practice). ([docs.anthropic.com](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/system-prompts))  
- **Explicit instructions and context** are used everywhere (Claude 4.6 prompting best practices). ([docs.anthropic.com](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/claude-4-best-practices))  
- **XML tags** clearly separate documents, tasks, constraints, outputs. ([platform.claude.com](https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/use-xml-tags))  
- **Long source files placed first; tasks last** to improve long-context performance. ([docs.anthropic.com](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/long-context-tips))  
- **Hallucination guardrails**: permission to stop/ask for missing files + “no guessing,” aligned with Anthropic’s recommendations. ([docs.anthropic.com](https://docs.anthropic.com/en/docs/test-and-evaluate/strengthen-guardrails/reduce-hallucinations))  
- **Examples/details vigilance**: the new prompt eliminates “maybe check…” ambiguity and turns it into required preflight verification steps. ([docs.anthropic.com](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/claude-4-best-practices))  

---

## If you want, I can tailor this further to your exact workflow

Two quick questions that would let me tighten Round 4 (restart damping) and Round 2 (barrier coordinate conversion) even more:

1. In your QuantLib v1.23 tree, does `FdmBlackScholesMesher` store `x = ln(S)` or `x = ln(S/K)`?  
2. In `FiniteDifferenceModel::rollback`, is `StepCondition::applyTo` called on **every step**, or only at **stopping times**?

If you paste those two file bodies (or just the relevant functions), I can rewrite those round specs to be completely unambiguous and “mechanically followable” (zero hidden assumptions).

---
Learn more:
1. [https://docs.anthropic.com/en/docs/test-and-evaluate/strengthen-guardrails/reduce-hallucinations](https://docs.anthropic.com/en/docs/test-and-evaluate/strengthen-guardrails/reduce-hallucinations)
2. [https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/long-context-tips](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/long-context-tips)
3. [https://docs.anthropic.com/en/docs/prompt-engineering](https://docs.anthropic.com/en/docs/prompt-engineering)
4. [https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/claude-4-best-practices](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/claude-4-best-practices)
5. [https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/use-xml-tags](https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/use-xml-tags)
6. [https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/system-prompts](https://docs.anthropic.com/en/docs/build-with-claude/prompt-engineering/system-prompts)

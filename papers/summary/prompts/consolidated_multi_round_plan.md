*Thinking...*

> The user wants me to:
> 
> 1. Do a more detailed, comprehensive fact-checking of the file-level modification plan against the actual QuantLib v1.23 codebase structure
> 2. Use Anthropic's official prompting guide to draft a more effective, detailed, and comprehensive multi-round prompt strategy for Claude Opus 4.6
> 
> Let me first do thorough fact-checking, then draft the prompts following Anthropic's best practices.
> 
> Let me start by carefully examining what's in the attached documents and cross-referencing claims.
> 
> ## Detailed Fact-Checking
> 
> ### 1. File paths and existence verification
> 
> From the `hpp_structure_of_quantlib.md`, the document describes the QuantLib v1.23 file structure. Let me verify claims in `file_level_plan.md` against this.
> 
> 
> I'm verifying the file paths and struct definitions mentioned in the plan against the actual QuantLib structure documentation. The `cranknicolsonscheme.hpp` file checks out, and I'm now examining the `FdmSchemeDesc` struct definition more carefully—specifically whether its members are marked as const, since the plan makes specific claims about the member declarations that I need to validate against what's actually documented. Looking at the actual struct definition, the members are plain non-const fields, which contradicts my earlier fact-check—adding a new field would be straightforward with aggregate initialization. Now I'm checking whether the coefficient arrays in `TripleBandLinearOp` are protected or public...
> 
> the source shows they're stored as protected/private members, but the class exposes `axpyb()` as a public method to modify coefficients, and you can also construct it with explicit arrays. Let me trace through how `FdmBlackScholesOp` actually assembles the operator in `setTime()` by combining the derivative maps...
> 
> The key insight is that `axpyb()` computes a linear combination of operators, so the fitted operator approach won't work cleanly here since the fitting formula produces combined coefficients that don't decompose into separate first and second derivative components. I need to either find a way to construct the `TripleBandLinearOp` directly with the fitted coefficients, or use a different strategy altogether. her. After construction, its arrays are initialized but can be overwritten via the public interface.
> 
> Looking more closely at QuantLib v1.23, I see that `TripleBandLinearOp::axpyb` lets me combine scaled operators, and `SecondDerivativeOp::mult` scales each node's coefficients by a corresponding array value. So I could create a scaled second derivative operator where node j gets `a_j * ρ_j`, then combine it with the first derivative using `axpyb`. Let me trace through how `FdmBlackScholesOp::setTime()` actually constructs the operator...
> 
> Examining the operator construction... I see that for the fitted version, I need to modify the diffusion coefficient to be `a_j * ρ_j` instead of just `a_j`, while keeping the convection term `b_j * d/dx` and the reaction term `-r`. This can be implemented using axpyb by computing the fitted diffusion coefficients and passing them to scale the second derivative operator.
> 
> Setting up the coefficient arrays... I need to think through whether we can selectively apply upwind discretization to individual nodes that fail the M-matrix check, or if we need to bake the fallback logic into the coefficient computation upfront. The challenge is that `axpyb` sets all nodes at once, so if a node needs upwind instead of the fitted scheme, I'd either need to compute the coefficients with that decision already made per-node, or modify individual matrix entries afterward. The M-matrix condition depends on checking whether the off-diagonal coefficients stay non-negative, and if they don't at a particular node, that's where we'd need to switch strategies.
> 
> Now I'm working through what the upwind discretization looks like at each node—the formulas change depending on the sign of the convection coefficient, which means mixing upwind and centered differences across the domain gets complicated.
> 
> The real constraint is that the axpyb framework uses fixed derivative stencils, so I can't selectively apply upwind to just the nodes that need it. I'm weighing three approaches: stick with axpyb everywhere, manually patch the operator coefficients afterward, or build the entire tridiagonal system node-by-node from scratch. The third option is cleanest for a fitted operator but requires direct access to the internal arrays of TripleBandLinearOp, so I need to figure out if there's a public interface for that.
> 
> Looking more closely at the class, there's no public setter for individual coefficients, but I could subclass to gain protected access, or store coefficients externally and reconstruct the operator each time. Actually, there might be a simpler path I'm overlooking in how TripleBandLinearOp is structured.
> 
> I see now that I can multiply operators element-wise by node values and add them together, so the strategy would be to build the first and second derivative operators separately, then combine them with the fitted diffusion coefficients at each node. The challenge is that centered first derivatives won't guarantee an M-matrix at problematic nodes, so I'll likely need to use `axpyb` for the standard fitted case and handle edge cases separately.
> 
> Actually, there's a better approach: instead of patching individual nodes after the fact, I can compute per-node convection coefficients that already satisfy the M-matrix condition by adjusting the fitting factor ρ_j upward until both the lower and upper diagonal entries stay non-negative. This means finding the minimum ρ value that keeps the off-diagonals positive.
> 
> For the lower diagonal, I need ρ_j ≥ |θ_j| where θ_j depends on the sign of the convection coefficient b_j, and I'm verifying that the hyperbolic cotangent relationship θ·coth(θ) ≥ |θ| holds across all cases to ensure this bound is tight. Now I'm substituting the expression for ρ_j back into the left coefficient formula and simplifying—it factors nicely to (b_j/(2h))·(coth(θ_j) - 1), which is always positive regardless of the sign of b_j since the hyperbolic cotangent maintains the right inequality in both cases.
> 
> For the right coefficient, I'm doing the same substitution and getting (b_j/(2h))·(coth(θ_j) + 1), which also stays positive across all scenarios.
> 
> So exponential fitting actually guarantees the M-matrix property for off-diagonals automatically — the upwind fallback is really just a numerical safeguard for edge cases like when a_j approaches zero and the fitting computation becomes unstable. This is a significant simplification that means the per-node upwind logic in the plan is mostly redundant for typical cases.
> 
> Given this, I can use the axpyb approach directly: compute the fitted diffusion coefficients per node, then apply them through the existing map operations without needing any fallback logic in the standard path. This also means I don't need protected access to the TripleBandLinearOp internals at all.
> 
> Now I'm checking whether FdmBackwardSolver can actually query the stopping times from its condition composite...
> 
> Verifying the remaining claims about the FDM infrastructure: FdmStepConditionComposite does expose stopping times, FdmBlackScholesMesher can accommodate a static factory for the sinh mesher, FdmCellAveragingInnerValue already handles cell averaging via Simpson integration, and FdmAmericanStepCondition uses straightforward projection logic.
> 
> However, FdBlackScholesVanillaEngine wasn't included in the provided headers despite being referenced in the structure documentation, so I need to flag that path verification. The legacy framework classes carry deprecation warnings, which suggests that modifying them may not be the best use of effort—the prompts should probably steer away from those changes.
> 
> For notifyDiscontinuity() and ParallelEvolver, the core issue is that QuantLib's evolvers are templates without a shared base class, so I'd either need to add the method to each template or use SFINAE/if constexpr. Since CrankNicolsonScheme already has it in the plan and we're deprioritizing legacy work anyway, this becomes less urgent. Now I'm checking the solve_splitting signature matching for FdmLinearOpComposite.
> 
> solve_splitting against TripleBandLinearOp.solve_splitting — the fitted operator needs to delegate correctly by mapping the direction parameter and passing the right arguments. I'm also verifying the mathematical foundation: the fitting factor ρ = θ·coth(θ) is correct per Duffy and Milev-Tagliani, with a stable Taylor expansion for small θ values.
> 
> The exponential fitting itself guarantees the M-matrix property through non-negative off-diagonals, so the upwind fallback is just a numerical safeguard rather than mathematically required. Rannacher smoothing follows Giles-Carter with two implicit Euler half-steps at dt/2 each, preserving the time grid alignment. For the cell-averaged payoff, I'm checking the Gauss-Legendre quadrature implementation.
> 
> Now I'm verifying the architectural claims about how FdmBackwardSolver integrates the time-stepping scheme...
> 
> The solver implements a two-phase approach: it first applies optional damping steps using implicit Euler if a different scheme is specified, then switches to the main scheme (Crank-Nicolson or others) for the primary rollback computation, with each scheme wrapped in its corresponding FiniteDifferenceModel template. So `FiniteDifferenceModel` is actually a shared template that both frameworks use—the new framework wraps `CrankNicolsonScheme` in it while the legacy one wraps `CrankNicolson<TridiagonalOperator>`. This means any changes to `FiniteDifferenceModel::rollbackImpl()` will impact both paths, but the monitoring-restart logic sits at a higher level in `FdmBackwardSolver::rollback()`, so modifying that solver is the right approach.
> 
> For `CrankNicolsonScheme`, I can verify it takes theta, map, bcSet, relTol, and solverType based on the hpp structure, and adding a dampingHalfSteps parameter with a default value of 0 should work fine.
> 
> Now I'm thinking through what makes effective prompts for Claude based on Anthropic's guidance — being specific and clear, providing relevant context inline, using XML tags for structure, asking for step-by-step reasoning on complex tasks, including examples of expected output, keeping each prompt focused on one task, establishing an expert role, and prefilling response starts for code generation. Let me draft an improved multi-round prompt strategy that applies these principles, particularly emphasizing that the most critical instructions should appear at the beginning and end of the prompt.
> 
> For the fact-checking section, I need to be systematic and thorough, checking each claim methodically across different categories like the file structure of QuantLib, class interfaces, and mathematical correctness.
> 
> Now I'm going through the file structure claims to verify which files exist in the QuantLib library, starting with the finite difference schemes and operators. Continuing through the finite differences meshers and step conditions, verifying the remaining header files exist in the expected locations. utilities/fdminnervaluecalculator.hpp` — EXISTS ✓
> A20. `methods/finitedifferences/utilities/fdmdirichletboundary.hpp` — EXISTS ✓
> A21. `pricingengines/vanilla/fdblackscholesvanillaengine.hpp` — NOT PROVIDED in attached headers; EXISTS in QuantLib but exact interface unverified ⚠️
> A22. `methods/finitedifferences/bsmoperator.hpp` — EXISTS (legacy) ✓
> A23. `methods/finitedifferences/mixedscheme.hpp` — EXISTS (legacy, deprecated) ⚠️
> A24. `methods/finitedifferences/cranknicolson.hpp` — EXISTS (legacy, deprecated) ⚠️
> A25. `methods/finitedifferences/finitedifferencemodel.hpp` — EXISTS (shared) ✓
> A26. `methods/finitedifferences/tridiagonaloperator.hpp` — EXISTS ✓
> 
> Now I'm checking the class interface claims. The FdmSchemeDesc members are likely non-const despite what was claimed earlier—the aggregate initialization pattern works fine with mutable members and that's more typical for QuantLib's style, though I'd need to verify against the actual header. For TripleBandLinearOp, the coefficient arrays are protected as expected, but axpyb() gives public access to modify them through linear combination.
> 
> FdmBlackScholesOp correctly stores its derivative maps, CrankNicolsonScheme holds the explicit and implicit solvers as shared pointers, and FdmBackwardSolver creates scheme objects locally during rollback. The step method accessibility between CrankNicolsonScheme and ImplicitEulerScheme needs checking—it's calling `implicit_->step(a, t, theta_)` which implies either a public overload or a friend relationship.
> 
> FdmCellAveragingInnerValue provides the averaging interface, FdmStepConditionComposite has the stopping times and factory methods I need, but there's a signature mismatch between solve_splitting on FdmLinearOpComposite and TripleBandLinearOp that I'll need to bridge—the composite takes direction and array parameters while the band operator takes different arguments. The exponential fitting approach guarantees an M-matrix property, with the upwind fallback serving as a numerical safeguard.
> 
> Now checking the mathematical foundations: the fitting factor formula ρ = θ·coth(θ) follows Il'in and Duffy's work, the small-angle approximation ρ ≈ 1 + θ²/3 is a valid Taylor expansion, Rannacher smoothing uses two half-steps at dt/2 per Giles-Carter, the sinh mesh transformation is correct, and the Greeks calculations for delta and gamma are properly formulated for log-space coordinates.
> 
> The PDE coefficients are set correctly, Gauss-Legendre quadrature nodes and weights check out, and on the architectural side I can add dampingHalfSteps to the scheme constructor with a sensible default, implement notifyDiscontinuity() callbacks, and refactor the backward solver to handle multiple time segments—though this last piece gets tricky with managing damping steps that consume time from each segment and handling boundary conditions properly.
> 
> For the operator assembly, using axpyb for the fitted operator works as I've verified. Policy iteration can modify the tridiagonal system by creating modified coefficient copies, but I need to figure out how to access and change those coefficients without running into protected member restrictions. The barrier projection condition fits naturally as a StepCondition, the sinh mesher can inherit from the base mesher class, and backward compatibility looks solid since the new constructor parameter defaults to zero. The main question now is whether the new FdmSchemeDesc field will cause issues if existing members are const. The test suite needs to follow QuantLib's Boost.Test conventions for file naming and registration. For the algorithmic components, I'm verifying the core methods: the policy iteration active-set approach aligns with Reisinger-Witte, the tridiagonal solver is already available in QuantLib, Richardson extrapolation is mathematically sound, and the oscillation detection via sign changes in delta-u is a solid heuristic. The monitoring projection correctly zeros out values outside the barrier range for double knock-out options, and the corner compatibility constraint ensures boundary nodes match the payoff at time zero.
> 
> For implementing policy iteration with TripleBandLinearOp, I'm weighing two strategies: either bypass the operator entirely and work directly with coefficient arrays and a standalone Thomas solver, or try to adapt the existing solve_splitting method by modifying the RHS and operator—but that approach gets complicated when trying to enforce identity rows at active nodes. The simpler path is extracting coefficients into plain arrays, modifying them for active-set constraints, and using a direct tridiagonal solver.
> 
> Now I'm reconsidering where to apply this logic. Rather than trying to extract and reconstruct coefficients from TripleBandLinearOp, I could work at the level of the full implicit system that the time stepper already constructs, where the structure is more amenable to modification.
> 
> The challenge is that the solve happens deep inside ImplicitEulerScheme::step() through a call to TripleBandLinearOp::solve_splitting(), which builds and solves the tridiagonal system internally. For policy iteration to work correctly, I'd need to intercept and modify that solve rather than post-process afterward. This points toward wrapping the entire time step rather than trying to hook into the linear solve itself.
> 
> The cleanest approach seems to be having FdmPolicyIterationLCP wrap the step function and work directly with the operator's apply() method to build the full tridiagonal system explicitly, then modify rows for active nodes and solve using Thomas algorithm. I can extract the matrix coefficients using toMatrixDecomp() which gives sparse matrices—for 1D problems that's just one tridiagonal matrix.
> 
> This is getting complicated, so I'll simplify the guidance: policy iteration needs careful integration, and the most practical path is working with extracted matrix coefficients rather than trying to modify TripleBandLinearOp in place.
> 
> Now I'm thinking through how to structure the comprehensive prompts more carefully. I want to expand the fact-check section significantly, follow Anthropic's best practices for prompting, use proper XML structure throughout, weave mathematical context directly into the explanations, give clearer implementation guidance based on what I've verified, handle the axpyb insight for the fitted operator, tackle the policy iteration integration challenge, and deprioritize any legacy framework changes.
> 
> Looking at the round structure, the user's existing plan is actually reasonable, but I can refine it based on my findings. Round 1 should cover the independent pieces like the fitted operator, sinh mesher, and diagnostics. Round 2 handles step condition and policy iteration which don't depend on Round 1. Round 3 would address scheme modifications like CN damping and other enhancements.
> 
> The key insight is that each round should include only the source files it actually needs, with clear markers for what to include. The shared preamble needs to be concise but mathematically complete so all rounds have what they need.
> 
> For the prompting approach itself, I should follow Anthropic's recommendations by putting coding conventions and the role definition in the system prompt equivalent at the start, then structure the user prompt to build complexity progressively.
> 
> Now I'm organizing my response around what the user actually asked for: comprehensive fact-checking first, then using those findings to inform the multi-round prompts. I'll structure this as a fact-check report by category, then show how those findings shape the implementation priority and the actual prompts themselves, applying Anthropic's principles about XML structure, chain-of-thought reasoning, and explicit constraints.
> 
> I should also preserve the user's approach of using `[INSERT ...]` placeholders for source files—it's clean and lets them fill in actual code—but be more specific about which files matter and why, especially since Claude Opus has enough context to handle multiple source files per round without overloading. I'm realizing I should structure this response with multiple components: a fact-check report organized by severity, a corrected implementation strategy based on those findings, and detailed multi-round prompts that apply Anthropic's best practices like clear instructions, XML organization, and strategic placement of key information. Each prompt needs a solid structure with system context, specific instructions, and output format guidance to work effectively with Claude. I'm recognizing that extended thinking or chain-of-thought reasoning is valuable for complex code generation tasks, especially when working with QuantLib's intricate class hierarchies and mathematical formulas. Since each round involves substantial reasoning about design, interactions with existing code, and edge cases, structuring prompts to encourage step-by-step thinking before implementation makes sense. With Claude's large context window, I can include meaningful source code excerpts while staying focused on what's relevant to each task.
> 
> Now I'm reconsidering whether the seven-round structure is optimal—Round 1 packs three new file pairs together which feels heavy, and Round 4's backward solver is particularly complex. Let me evaluate if some rounds should be split or consolidated for better pacing and manageability.
> 
> Actually, the three classes in Round 1 are independent enough that a capable model should handle them together, though I could split it into separate rounds if context becomes an issue. I'll stick with seven rounds for now but keep that flexibility in mind.
> 
> I should finalize my response now rather than continuing to second-guess the structure. I have enough knowledge of Anthropic's prompting best practices to write this without needing to search for the latest guide.
> 
> The key principles I want to apply are being specific and direct, using concrete examples, structuring with XML tags, assigning Claude an expert role, leveraging the multi-round approach for complexity, prefilling responses with file headers, specifying complete files rather than diffs, and strategically placing critical information for long context windows.
> 
> Now I need to write out the comprehensive response with detailed prompts and examples. Given the scope—fact-checking plus seven detailed prompts with source code—I'll aim for thorough coverage while keeping each prompt to about 1-2 pages of instructions plus code placeholders, rather than padding to an excessive length.

# Comprehensive Fact-Check & Multi-Round Prompt Strategy

This response has two parts: a thorough fact-check of the `file_level_plan.md` against the actual QuantLib v1.23 codebase (as described in `hpp_structure_of_quantlib.md` and the attached headers), followed by a multi-round prompt strategy designed according to Anthropic's prompt engineering documentation.

---

## Part A — Comprehensive Fact-Check

I organize findings by category, with severity ratings: **✅ Verified**, **⚠️ Needs Correction**, **❌ Incorrect**, **❓ Unverifiable** (source not provided).

### A. File Paths and Existence

| # | Claimed Path | Status | Notes |
|---|---|---|---|
| A1 | `methods/finitedifferences/schemes/cranknicolsonscheme.hpp` | ✅ | Confirmed in hpp_structure §2.1 |
| A2 | `methods/finitedifferences/schemes/impliciteulerscheme.hpp` | ✅ | Confirmed |
| A3 | `methods/finitedifferences/schemes/expliciteulerscheme.hpp` | ✅ | Confirmed |
| A4 | `methods/finitedifferences/operators/fdmblackscholesop.hpp` | ✅ | Confirmed in hpp_structure §2.3 |
| A5 | `methods/finitedifferences/operators/triplebandlinearop.hpp` | ✅ | Confirmed |
| A6 | `methods/finitedifferences/operators/firstderivativeop.hpp` | ✅ | Confirmed |
| A7 | `methods/finitedifferences/operators/secondderivativeop.hpp` | ✅ | Confirmed |
| A8 | `methods/finitedifferences/solvers/fdmbackwardsolver.hpp` | ✅ | Confirmed in hpp_structure §2.4 |
| A9 | `methods/finitedifferences/solvers/fdmblackscholessolver.hpp` | ✅ | Confirmed |
| A10 | `methods/finitedifferences/solvers/fdm1dimsolver.hpp` | ✅ | Confirmed |
| A11 | `methods/finitedifferences/solvers/fdmsolverdesc.hpp` | ✅ | Confirmed |
| A12 | `methods/finitedifferences/meshers/fdmblackscholesmesher.hpp` | ✅ | Confirmed |
| A13 | `methods/finitedifferences/meshers/concentrating1dmesher.hpp` | ✅ | Confirmed |
| A14 | `methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp` | ✅ | Confirmed |
| A15 | `methods/finitedifferences/stepconditions/fdmamericanstepcondition.hpp` | ✅ | Confirmed |
| A16 | `methods/finitedifferences/utilities/fdminnervaluecalculator.hpp` | ✅ | Confirmed |
| A17 | `methods/finitedifferences/utilities/fdmdirichletboundary.hpp` | ✅ | Confirmed |
| A18 | `pricingengines/vanilla/fdblackscholesvanillaengine.hpp` | ❓ | **Not provided** in attached headers. hpp_structure §2.8 describes it but actual file absent. Plan marks `[VERIFY]`. |
| A19 | `methods/finitedifferences/bsmoperator.hpp` (legacy) | ✅ | Confirmed but carries `QL_DEPRECATED` |
| A20 | `methods/finitedifferences/mixedscheme.hpp` (legacy) | ✅ | Confirmed but deprecated |
| A21 | `methods/finitedifferences/cranknicolson.hpp` (legacy) | ✅ | Confirmed but deprecated |
| A22 | `methods/finitedifferences/finitedifferencemodel.hpp` | ✅ | Confirmed; shared between frameworks |
| A23 | `methods/finitedifferences/tridiagonaloperator.hpp` | ✅ | Confirmed |

**Finding:** All new-framework file paths verified. The critical unverifiable file is the pricing engine (A18). Legacy files exist but are deprecated. The plan should note that `finitedifferencemodel.hpp` is shared between both frameworks — changes to `rollbackImpl` there affect legacy engines too.

---

### B. Class Interfaces and Member Access

**B1 — `FdmSchemeDesc` members claimed as `const`**

The plan at §2.3.1 says the struct has `const FdmSchemeType type; const Real theta, mu;`. However, the hpp_structure shows aggregate initialization patterns (`return {CrankNicolsonType, 0.5, 0.0}`) that work with both const and non-const members. Examining QuantLib v1.23's actual convention: `FdmSchemeDesc` is a simple data struct. Looking at how factories return it by value with brace-initialization, and how `FdmBackwardSolver` stores `const FdmSchemeDesc schemeDesc_` (making the *copy* const, not the struct's members), the members themselves are most likely non-const. **Severity: ⚠️.** The previous prompt's claim that they "ARE const" is likely wrong, but without seeing the actual header, this cannot be definitively resolved. The safest implementation adds a new member with a default value and updates all factories.

**B2 — `TripleBandLinearOp` coefficient access**

The plan §1.1 states: "write directly into `mapT_` arrays `lower_[j]`, `diag_[j]`, `upper_[j]`." These are protected `unique_ptr<Real[]>` members. **Severity: ❌ as stated.** However, this can be resolved via the `axpyb()` public method. I prove this works for exponential fitting:

The existing `FdmBlackScholesOp::setTime()` assembles:
```cpp
mapT_.axpyb(convection_array, dxMap_, dxxMap_.mult(diffusion_array), discount_array);
```
This sets `mapT_ = convection · dxMap_ + dxxMap_scaled + discount · I`.

For the fitted operator, the desired stencil at node j is:
```
ℓ_j = a_j·ρ_j/h² − b_j/(2h)     u_j = a_j·ρ_j/h² + b_j/(2h)
```
This decomposes as `b_j · (centered-first-derivative) + a_j·ρ_j · (centered-second-derivative) − r · I`, which means we can use:
```cpp
Array fitted_diffusion(N);
for (j = 0; j < N; ++j) fitted_diffusion[j] = a_j * rho_j;
mapT_.axpyb(convection_coeff, dxMap_, dxxMap_.mult(fitted_diffusion), Array(1, -r));
```
This uses only public methods. **No protected access needed.**

**B3 — M-matrix guarantee from exponential fitting (mathematical proof)**

The plan's §1.1 describes a per-node upwind fallback when `ℓ_j < 0` or `u_j < 0`. This fallback is mathematically unnecessary for the exponential fitting case. Here is a complete proof:

Given ρ_j = θ_j·coth(θ_j) where θ_j = b_j·h/(2a_j):

For ℓ_j = a_j·ρ_j/h² − b_j/(2h):
- Substitute ρ_j: ℓ_j = (b_j/(2h))·coth(θ_j) − b_j/(2h) = (b_j/(2h))·(coth(θ_j) − 1)
- When b_j > 0: θ_j > 0, coth(θ_j) > 1, so ℓ_j > 0 ✓
- When b_j < 0: θ_j < 0, coth(θ_j) < −1, so (coth(θ_j) − 1) < −2, and b_j/(2h) < 0, giving ℓ_j > 0 ✓
- When b_j = 0: ρ_j = 1, ℓ_j = a_j/h² > 0 ✓

Similarly for u_j = (b_j/(2h))·(coth(θ_j) + 1) ≥ 0 by analogous argument.

**Severity: ⚠️.** The plan's upwind fallback mechanism is overengineered for the standard case. It should be retained only as a degenerate-case guard (a_j ≈ 0), not as a primary design feature. The `mMatrixSatisfied()` diagnostic should virtually always return `true` for exponentially fitted operators.

**B4 — `solve_splitting` signature mismatch**

Plan §1.1 says `FdmFittedBlackScholesOp` implements `FdmLinearOpComposite::solve_splitting(Size direction, const Array& r, Real s)`. But `TripleBandLinearOp::solve_splitting` has signature `(const Array& r, Real a, Real b)`. The bridge is straightforward — `FdmBlackScholesOp` already does this: `return mapT_.solve_splitting(r, s, 1.0)` when `direction == direction_`. **Severity: ✅** (correctly identified in hpp_structure, but plan should be explicit about the bridge).

**B5 — `ImplicitEulerScheme::step` accessibility**

The plan §2.1.1 has `CrankNicolsonScheme` calling `implicit_->step(a, t, 1.0)` during damping. The existing `CrankNicolsonScheme` already calls `implicit_->step(a, t, theta_)` in its normal operation, confirming this is accessible. **Severity: ✅.**

**B6 — `FdmCellAveragingInnerValue` already does cell averaging**

hpp_structure §2.7 confirms: "computes cell-averaged payoffs using Simpson integration." The plan's proposed Gauss-Legendre upgrade is a minor accuracy enhancement, not a critical change. **Severity: ✅** (plan correctly identifies this as existing, but overstates the importance of changing quadrature).

**B7 — Policy iteration with `TripleBandLinearOp`**

Plan §1.7 describes modifying a "copy of the `TripleBandLinearOp` coefficients." But `TripleBandLinearOp` doesn't expose public coefficient setters, and copying it doesn't give you write access to the copy's protected members either. **Severity: ⚠️.** The implementation must either:
(a) Use `toMatrixDecomp()` to extract a sparse matrix representation, then work with raw arrays, or
(b) Avoid modifying the operator entirely — instead, the policy iteration class extracts the tridiagonal coefficients from the system that `ImplicitEulerScheme` would build, operates on plain arrays, and uses its own Thomas solver.

Option (b) is cleaner and avoids fighting the encapsulation.

---

### C. Architectural Feasibility

**C1 — `FdmBackwardSolver::rollback()` segmentation**

The plan §2.3.1 describes splitting rollback into segments at monitoring times. The current implementation (per hpp_structure §2.4) creates a `FiniteDifferenceModel<Scheme>` and calls `model.rollback(rhs, from, to, steps, condition)`. For segmentation, the backward solver must call `model.rollback()` multiple times with different `[from, to]` ranges.

**Issue:** The step condition composite's `applyTo(a, t)` is called inside `FiniteDifferenceModel::rollbackImpl()` at stopping times. If we split the rollback at monitoring times, we need the projection to happen at the correct moment. If we split *before* the monitoring time, the projection happens at the start of the next segment. If we split *at* the monitoring time, the projection may be missed or doubled.

**Resolution:** The most robust approach is: for each segment `[from_seg, to_seg]` where `to_seg` is a monitoring time, run `model.rollback(rhs, from_seg, to_seg, steps_seg, condition)`. The step condition fires at `to_seg` (a stopping time), applying the projection. Then `notifyDiscontinuity()` is called, and the next segment starts with damping.

**Severity: ⚠️** — Feasible but requires careful handling of the interaction between `FiniteDifferenceModel`'s internal stopping-time logic and the outer segmentation. The plan should note this explicitly.

**C2 — Legacy framework modifications**

The plan proposes creating `FittedBSMOperator` (§1.3-1.4) and modifying `MixedScheme` (§2.7.2) and `FiniteDifferenceModel` (§2.7.3). Given that legacy FD classes are deprecated in v1.23, these modifications:
- May not be compiled if users disable deprecated code
- Risk breaking rarely-tested legacy engines
- Are unlikely to be accepted in QuantLib upstream

**Severity: ⚠️.** Recommend deprioritizing to optional/experimental status.

**C3 — `FdmSolverDesc` diagnostics field**

Plan §2.3.6 proposes adding `FdmDiagnostics::Level diagLevel` to `FdmSolverDesc`. Since `FdmSolverDesc` is a plain aggregate struct (no constructor), adding a field requires either:
- Providing a default value (C++14 in-class initializer, which QuantLib v1.23 supports)
- Or updating all construction sites

**Severity: ⚠️** — Minor, but all existing code that constructs `FdmSolverDesc` via brace-initialization must be updated or the new field must have a default.

---

### D. Mathematical Correctness

| # | Claim | Status | Notes |
|---|---|---|---|
| D1 | ρ = θ·coth(θ) for fitting factor | ✅ | Standard Il'in/Duffy result |
| D2 | Small-θ: ρ ≈ 1 + θ²/3 for \|θ\| < 1e-8 | ✅ | Taylor series correct; threshold reasonable |
| D3 | Stencil: ℓ = aρ/h² − b/(2h), u = aρ/h² + b/(2h), d = −2aρ/h² − r | ✅ | Consistent with consolidated_summary §3.2.4 |
| D4 | 2 implicit Euler half-steps at dt/2 for Rannacher | ✅ | Per Giles-Carter (2006) |
| D5 | Δ = u_x / S, Γ = (u_xx − u_x) / S² | ✅ | Correct for x = ln(S/K) |
| D6 | PDE coefficients a = σ²/2, b = (r−q) − σ²/2 | ✅ | Standard log-price form |
| D7 | Gauss-Legendre 3-point: ξ = {−√(3/5), 0, +√(3/5)}, w = {5/9, 8/9, 5/9} | ✅ | Exact for polynomials to degree 5 |
| D8 | Sinh mesh: x(ξ) = x_center + c·sinh(α·(ξ − ξ₀)) | ✅ | Standard concentration formula |
| D9 | Monitoring projection: u[j] ← u[j]·𝟙_{[L,U]}(S[j]) | ✅ | Correct for double knock-out |
| D10 | Policy iteration converges in O(J) iterations | ⚠️ | Per Reisinger-Witte it converges in few iterations (typically 2-5), not O(J). The plan correctly states maxIter=50 as a guard. |

---

### E. Consistency Between Documents

**E1 — Module M5 → fitted operator mapping**

The `consolidated_code.md` §2.M5 describes `SpatialOperator_assemble()` with M-matrix check and upwind fallback. The `file_level_plan.md` §1.1 maps this to `FdmFittedBlackScholesOp`. The mapping is consistent, except the plan overcomplicates the fallback (see B3 above).

**E2 — Module M7 → CrankNicolsonScheme damping**

Both `consolidated_code.md` §2.M7 and `file_level_plan.md` §2.1.1 describe the same damping state machine. Consistent. ✅

**E3 — Module M11 → diagnostics**

Both documents define the same diagnostic features. The `file_level_plan.md` provides a QuantLib-specific interface. Consistent. ✅

**E4 — Test cases**

Both `consolidated_code.md` §5.1 and `file_level_plan.md` §5 define tests T1–T7 with the same parameter sets. Minor differences: the plan adds a down-and-out rebate call (T5) while the code document has it as T5 with slightly different parameters. These are reconcilable.

**E5 — Discrepancy in American option time stepper**

`consolidated_code_alternative.md` §4.1 recommends `IMPLICIT_EULER` for American options, while `file_level_plan.md` §2.5.2 proposes policy iteration as an opt-in within the backward solver (compatible with any time stepper). The alternative code document's recommendation is the more conservative choice; the plan's approach is more flexible but requires careful integration.

---

### F. Critical Omissions in the Plan

**F1 — No guidance on `FiniteDifferenceModel::rollbackImpl` interaction.** The new-framework `FdmBackwardSolver` wraps scheme objects inside `FiniteDifferenceModel<Scheme>` and calls its `rollback`. The segmentation logic must account for how `FiniteDifferenceModel` internally handles stopping times, or risk double-application of conditions.

**F2 — No discussion of thread safety.** If `FdmFittedBlackScholesOp` stores mutable diagnostics counters (`mMatrixOk_`, `fallbackCount_`), these may cause issues in multithreaded pricing. Mark as `mutable` with appropriate documentation.

**F3 — No guidance on `ext::shared_ptr` vs `std::shared_ptr`.** QuantLib v1.23 uses `ext::shared_ptr` (aliased via `ql/shared_ptr.hpp`). All new code must use `ext::shared_ptr`, not `std::shared_ptr`.

**F4 — No `Disposable<>` convention clarity.** QuantLib v1.23 uses `Disposable<Array>` for return-value optimization (now largely unnecessary with move semantics but still present in the API). New methods returning `Array` should match the convention of similar existing methods.

**F5 — Missing `#include` dependency analysis.** The plan lists which `all.hpp` files to update but doesn't specify the include dependencies within each new file (e.g., `fdmfittedblackscholesop.hpp` needs `#include <ql/methods/finitedifferences/operators/triplebandlinearop.hpp>`, `#include <ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp>`, etc.).

---

## Part B — Corrected Implementation Priority

Based on the fact-check, here is the updated priority order with corrections applied:

| Round | Content | Key Corrections Applied |
|---|---|---|
| 1 | Fitted operator + sinh mesher + diagnostics | Use `axpyb` for fitted operator; simplify M-matrix fallback to degenerate-case guard only |
| 2 | Barrier projection condition + policy iteration | Policy iteration works on plain arrays extracted from system, not on `TripleBandLinearOp` directly |
| 3 | `CrankNicolsonScheme` damping + `FdmSchemeDesc` update | Don't assume const members; add constructor overload + update all static factories |
| 4 | `FdmBackwardSolver` + `FdmBlackScholesSolver` | Account for `FiniteDifferenceModel` stopping-time interaction in segmentation logic |
| 5 | Wiring: composite factory, mesher factory, headers | Minimal changes; use C++14 default member initializers where possible |
| 6 | Pricing engine | Verify file existence before modifying; provide fallback if path differs |
| 7 | Test suite | Follow QuantLib's Boost.Test conventions |

Legacy framework modifications (§1.3–1.4, §2.7.2–2.7.3) are deprioritized to an optional Round 8.

---

## Part C — Anthropic Prompt Engineering Principles Applied

The following Anthropic-recommended techniques are used in the prompts below:

**Technique 1 — XML tags for structure.** Each prompt uses `<system_context>`, `<mathematical_context>`, `<source_files>`, `<task>`, `<implementation_guidance>`, `<constraints>`, and `<output_specification>` tags to clearly delimit sections, following Anthropic's recommendation for long-context prompts.

**Technique 2 — Role assignment.** The system context establishes Claude as an expert C++ developer familiar with QuantLib, which Anthropic's docs recommend for domain-specific tasks.

**Technique 3 — Be specific and direct.** Each task is precisely scoped with exact file paths, class names, method signatures, and design decisions already made. Ambiguous choices are resolved in the prompt rather than left to the model.

**Technique 4 — Put critical instructions at start and end.** Each prompt begins with the most important constraints (output format, coding conventions) and ends with a quality checklist that reinforces them.

**Technique 5 — One concept per prompt.** Each round targets a single architectural layer (operators, conditions, schemes, solvers, wiring, engine, tests) to avoid cognitive overload.

**Technique 6 — Provide examples.** The `<source_files>` section provides existing QuantLib code that serves as a style example. The `<implementation_guidance>` section shows specific code patterns to follow.

**Technique 7 — Chain of thought.** Each prompt includes `<implementation_guidance>` that walks through design decisions, encouraging Claude to reason about the approach before generating code.

**Technique 8 — Explicit output format.** The `<output_specification>` tag specifies exactly which files to produce, that they must be complete (not diffs), and that they must compile against existing QuantLib headers.

**Technique 9 — Long context management.** Source files are included via `[INSERT ...]` placeholders so the user controls context size. Only files relevant to the current round are included.

**Technique 10 — Prefilling.** Each prompt ends with a suggested response start (copyright header) to guide format.

---

## Part D — Multi-Round Prompts

### SHARED PREAMBLE

This preamble is included verbatim at the start of every round. It is approximately 1,500 tokens and provides the stable context that all rounds need.

```xml
<system_context>
You are an expert C++ developer implementing improvements to QuantLib v1.23's
finite difference framework for Black-Scholes option pricing. You have deep
familiarity with QuantLib's coding conventions, class hierarchy, and the modern
Fdm* framework (not the deprecated legacy MixedScheme framework).

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
PDE in log-price x = ln(S/K), time-to-maturity τ (τ=0 at payoff, τ=T today):

  -u_τ + a(x,τ)·u_xx + b(x,τ)·u_x − r·u = 0

where a(x,τ) = σ²/2,  b(x,τ) = (r−q) − σ²/2.

EXPONENTIAL FITTING:
  θ_j = b_j·h / (2·a_j)       (local Péclet parameter)
  ρ_j = θ_j·coth(θ_j)         (fitting factor; ρ ≈ 1 + θ²/3 for |θ| < 1e-8)

Fitted stencil at interior node j (uniform grid spacing h):
  ℓ_j = a_j·ρ_j/h² − b_j/(2h)     (lower off-diagonal)
  d_j = −2·a_j·ρ_j/h² − r          (main diagonal)
  u_j⁺ = a_j·ρ_j/h² + b_j/(2h)    (upper off-diagonal)

IMPORTANT: Exponential fitting GUARANTEES ℓ_j ≥ 0 and u_j⁺ ≥ 0 for all
finite a_j > 0 (proof: ℓ_j = (b_j/(2h))·(coth(θ_j)−1) which is non-negative
for all signs of b_j). Upwind fallback is only needed as a degenerate guard
when a_j ≈ 0.

This stencil decomposes as: b_j·(∂/∂x) + a_j·ρ_j·(∂²/∂x²) − r·I,
which means it can be assembled via TripleBandLinearOp::axpyb():
  mapT_.axpyb(b_array, firstDerivOp, secondDerivOp.mult(aρ_array), Array(1,-r))

RANNACHER-SMOOTHED CN: After each discontinuity event (payoff at τ=0,
monitoring projection), perform 2 implicit Euler half-steps at dt/2, then
resume Crank-Nicolson.

GREEKS: Δ = u_x/S,  Γ = (u_xx − u_x)/S²  where x = ln(S/K).

SINH MESH: x(ξ) = x_center + c·sinh(α·(ξ − ξ₀)) for ξ ∈ [0,1], with ξ₀
found by bisection to satisfy endpoint constraints x(0)=xMin, x(1)=xMax.
</mathematical_context>
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

DESIGN DECISION (already made — do not deviate):
The operator stores FirstDerivativeOp dxMap_, SecondDerivativeOp dxxMap_, and
TripleBandLinearOp mapT_ — the same member layout as FdmBlackScholesOp. The
setTime() method computes per-node fitting factors and assembles mapT_ using
the public axpyb() method:

  void FdmFittedBlackScholesOp::setTime(Time t1, Time t2) {
      // For each layout point i:
      //   1. Get sigma from volTS_ or localVol_
      //   2. Compute a_i = 0.5*sigma*sigma, b_i = r-q-a_i (in log-space)
      //   3. Compute theta_i = b_i*h/(2*a_i), rho_i = fittingFactor(theta_i)
      //   4. Store b_i in convection array, a_i*rho_i in fittedDiffusion array
      // Then: mapT_.axpyb(convection, dxMap_,
      //                    dxxMap_.mult(fittedDiffusion), Array(1, -r));
  }

The static fittingFactor method:
  static Real fittingFactor(Real theta) {
      if (std::fabs(theta) < 1e-8)
          return 1.0 + theta*theta/3.0;
      return theta / std::tanh(theta);  // = theta * coth(theta)
  }

Note: coth(θ) = 1/tanh(θ), so θ·coth(θ) = θ/tanh(θ).

For degenerate cases (a_i ≈ 0): clamp a_i to a small positive value
(e.g., 1e-20) before computing theta. This is the ONLY fallback needed —
do not implement per-node upwind switching.

The mMatrixSatisfied() diagnostic method scans mapT_ coefficients after
setTime() and verifies all off-diagonals have the correct sign. This should
almost always return true (see mathematical proof in preamble).

Interface must match FdmLinearOpComposite exactly. Use FdmBlackScholesOp as
your template — the two classes have nearly identical structure, differing
only in setTime().

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
1. Build uniform ξ_j = j/(size-1) for j = 0..size-1
2. Find ξ₀ by bisection such that xMin = xCenter + c·sinh(α·(0 − ξ₀)) and
   xMax = xCenter + c·sinh(α·(1 − ξ₀)), where c = (xMax−xMin)/(2·sinh(α/2))
   for the symmetric case. For asymmetric (xCenter ≠ midpoint), bisect on ξ₀.
3. Compute x_j = xCenter + c·sinh(α·(ξ_j − ξ₀))
4. Force x[0] = xMin, x[size-1] = xMax (numerical safety)
5. For each alignment target: find nearest node, compute shift, apply shift
   to entire grid if shift < min(h)/2
6. Populate locations_, dplus_, dminus_ from Fdm1dMesher base class

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
□ Sinh mesher populates locations_, dplus_, dminus_ correctly
□ Sinh mesher handles alpha ≈ 0 gracefully (falls back to uniform)
□ Diagnostics oscillationScore is O(N) and allocation-free
□ All ext::shared_ptr, not std::shared_ptr
□ All QL_REQUIRE for preconditions
□ Include guards match QuantLib naming convention
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

FILE: ql/pde/stepcondition.hpp
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

The constructor precomputes outsideIndices_: a sorted vector of layout indices
whose spatial coordinate (in S-space, obtained via mesher->location(iter, direction)
converted from log-space via exp()) falls outside [lowerBarrier, upperBarrier].

  void applyTo(Array& a, Time t) const override;
  // If t matches any monitoring time within tolerance (1e-10), set a[i]=0
  // for all i in outsideIndices_. Otherwise no-op.

  const std::vector<Time>& monitoringTimes() const;
  // Accessor for FdmStepConditionComposite to register as stopping times

IMPORTANT: The mesher stores locations in LOG-space (x = ln(S/K)), so to
check against S-space barriers, convert: S_j = K * exp(x_j). However, the
mesher does not directly store K. Instead, use exp(mesher->location(iter, dir))
which gives S/K, then multiply by strike... Actually, check how
FdmBlackScholesMesher stores locations. The locations are in ln(S) space
(not ln(S/K)). Verify by examining how FdmBlackScholesOp reads coordinates.

Design choice: store barriers in log-space internally after converting in
the constructor:
  lnLower_ = std::log(lowerBarrier);
  lnUpper_ = std::log(upperBarrier);
Then compare directly against mesher->location(iter, direction).


=== FILE PAIR 2: FdmPolicyIterationLCP ===

Path: ql/methods/finitedifferences/utilities/fdmpolicyiteration.hpp
      ql/methods/finitedifferences/utilities/fdmpolicyiteration.cpp

IMPORTANT DESIGN DECISION: This class operates on PLAIN ARRAYS, not on
TripleBandLinearOp directly. This avoids the protected-member-access problem.

The caller (typically a modified backward solver or step condition) provides:
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

      // Solve: find u >= phi such that A·u >= rhs with complementarity
      // lower/diag/upper define tridiagonal A (size n)
      // rhs is right-hand side (size n)
      // phi is exercise values (size n)
      // Returns solution u
      Array solve(
          const Array& lower,     // a[1..n-1] (lower diagonal, size n-1 or n with a[0]=0)
          const Array& diag,      // b[0..n-1] (main diagonal, size n)
          const Array& upper,     // c[0..n-2] (upper diagonal, size n-1 or n with c[n-1]=0)
          const Array& rhs,       // right-hand side (size n)
          const Array& phi) const;  // exercise values (size n)

      Size lastIterationCount() const;

    private:
      Size maxIterations_;
      Real tolerance_;
      mutable Size lastIterations_;

      // Internal Thomas solver for modified systems
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
□ Policy iteration converges in constant iterations for typical American options
□ Thomas solver handles n=1 edge case
□ All Arrays are passed by const reference where not modified
□ No TripleBandLinearOp dependency in policy iteration
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
[INSERT FULL CONTENT — if available; otherwise note that you are reconstructing
the implementation from the header interface and the code description in
hpp_structure_of_quantlib.md section 2.1]

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

Verify that ImplicitEulerScheme::setStep() only stores dt_ and does not
trigger expensive recomputation. If it does, add a comment noting the
performance implication.

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
□ step() during damping performs exactly 2 implicit Euler half-steps per call
□ implicit_->setStep is restored to dt_ after damping half-steps
□ numberOfIterations() still works correctly (returns implicit's count)
□ FdmBackwardSolver class declaration is byte-for-byte identical to original
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
[INSERT FULL CONTENT if available — otherwise note reconstruction needed]

FILE: ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp
[INSERT FULL CONTENT if available]

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

The rollback() method must be restructured to support monitoring-restart
damping. The current structure (from hpp_structure document) is:

  void FdmBackwardSolver::rollback(..., Size steps, Size dampingSteps) {
      // Phase 1: damping steps with implicit Euler
      if (dampingSteps != 0) {
          ImplicitEulerScheme implicitEvolver(map_, bcSet_);
          FiniteDifferenceModel<ImplicitEulerScheme> dampingModel(implicitEvolver, ...);
          dampingModel.rollback(rhs, from, dampingTo, dampingSteps, *condition_);
      }
      // Phase 2: main scheme steps
      switch (schemeDesc_.type) {
          case CrankNicolsonType: {
              CrankNicolsonScheme cnEvolver(schemeDesc_.theta, map_, bcSet_);
              FiniteDifferenceModel<CrankNicolsonScheme> cnModel(cnEvolver, ...);
              cnModel.rollback(rhs, dampingTo, to, steps-dampingSteps, *condition_);
          } break;
          // ... other scheme types
      }
  }

NEW STRUCTURE when schemeDesc_.monitoringDampingSteps > 0 AND scheme is
CrankNicolsonType:

  1. Query condition_->stoppingTimes() to get monitoring times
  2. Filter to times within [to, from] (backward: from > to)
  3. Sort these times in descending order (backward marching)
  4. Add segment boundaries: {from, mon_1, mon_2, ..., mon_F, to}
  5. For each segment [segFrom, segTo]:
     a. If this is the first segment: use initial dampingSteps as before
     b. If this is after a monitoring event: the CN scheme will be in damping
        state (notifyDiscontinuity was called at end of previous segment)
     c. Create CN scheme with monitoringDampingSteps passed through
     d. Wrap in FiniteDifferenceModel and rollback this segment
     e. After rollback: if segTo is a monitoring time, call
        cnEvolver.notifyDiscontinuity() — BUT the evolver is local to the
        FiniteDifferenceModel, so this must happen BETWEEN segments

CRITICAL SUBTLETY: FiniteDifferenceModel<CrankNicolsonScheme>::rollback()
internally handles stopping times from the condition and calls
condition->applyTo(a, t) at those times. The monitoring projection is
applied by the condition inside the model's rollback. The
notifyDiscontinuity() call must happen AFTER the model's rollback returns
(i.e., after the projection has been applied at the segment boundary).

IMPLEMENTATION APPROACH: Create a single CrankNicolsonScheme for the entire
rollback. Between segments, call notifyDiscontinuity() on it. For each
segment, create a new FiniteDifferenceModel wrapping the SAME scheme object.

Wait — FiniteDifferenceModel takes the scheme by value (or reference?).
Check the template: if it copies the scheme, the damping state won't persist
across segments. If by reference, it will.

Looking at finitedifferencemodel.hpp: the evolver is stored as a member
`Evolver evolver_`. If the constructor takes it by value, modifications to
the external scheme won't be seen. If by reference, they will.

SAFEST APPROACH: Create the CN scheme, then for each segment:
  a. Set the scheme's step size for this segment
  b. Create FiniteDifferenceModel with the scheme (passed in)
  c. Call model.rollback()
  d. After return, call scheme.notifyDiscontinuity() if at monitoring boundary

If FiniteDifferenceModel copies the evolver, then we need an alternative:
use the scheme's built-in damping (from Round 3) which automatically does
Euler half-steps after notifyDiscontinuity(). But if the copy happens before
notifyDiscontinuity, the state is lost.

FALLBACK APPROACH (most robust): Don't rely on the scheme's internal
damping state at all. Instead, at each segment boundary after a monitoring
event, explicitly run dampingSteps of implicit Euler (just like the initial
damping), then run CN for the rest of the segment. This reuses the existing
two-phase pattern and avoids the copy-vs-reference issue entirely.

Choose this fallback approach. The modified rollback becomes:

  void rollback(rhs, from, to, steps, dampingSteps) {
      if (schemeDesc_.monitoringDampingSteps == 0 ||
          schemeDesc_.type != CrankNicolsonType) {
          // Original behavior unchanged
          originalRollback(rhs, from, to, steps, dampingSteps);
          return;
      }

      // Get monitoring times in [to, from]
      auto monTimes = filterMonitoringTimes(condition_->stoppingTimes(), to, from);

      // Build segment boundaries: from, mon_1, ..., mon_F, to
      auto segments = buildSegments(from, to, monTimes);

      // Distribute steps across segments proportional to time span
      auto stepsPerSeg = distributeSteps(segments, steps);

      for (size_t seg = 0; seg < segments.size()-1; ++seg) {
          Time segFrom = segments[seg];
          Time segTo = segments[seg+1];
          Size segSteps = stepsPerSeg[seg];
          Size segDamping = (seg == 0) ? dampingSteps
                                        : schemeDesc_.monitoringDampingSteps;

          // Reuse existing two-phase pattern per segment
          rollbackSegment(rhs, segFrom, segTo, segSteps, segDamping);
      }
  }

Where rollbackSegment encapsulates the original "damping then main" pattern.


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
- When monitoringDampingSteps > 0 but no monitoring times exist in the
  condition's stopping times: behavior identical to monitoringDampingSteps==0
- Do NOT modify FiniteDifferenceModel — work around its copy semantics
- The switch statement for other scheme types (Douglas, Hundsdorfer, etc.)
  must remain unchanged
- FdmBlackScholesSolver's existing constructor signature must remain valid
  (new parameter has default value)
</constraints>

<output_specification>
Produce exactly 4 files:
1. fdmbackwardsolver.hpp (modified — FdmBackwardSolver class gains
   rollbackSegment helper; FdmSchemeDesc from Round 3 retained)
2. fdmbackwardsolver.cpp (modified — restructured rollback)
3. fdmblackscholessolver.hpp (modified — new constructor parameter)
4. fdmblackscholessolver.cpp (modified — conditional operator creation)
</output_specification>

<quality_checklist>
□ Original rollback behavior preserved when monitoringDampingSteps == 0
□ Segment boundaries computed correctly in backward-time direction
□ Step counts distributed proportionally and sum to total
□ Each segment's damping + main stepping fills exactly that segment's time span
□ No double-application of monitoring projections
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
[INSERT FULL CONTENT if available]

FILE: ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp
[INSERT FULL CONTENT]

FILE: ql/methods/finitedifferences/meshers/fdmblackscholesmesher.cpp
[INSERT FULL CONTENT if available]

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
5. Returns the assembled composite

IMPLEMENTATION NOTE: Look at how vanillaComposite() builds its condition list
and stopping times, then replicate the pattern with the addition of the
barrier projection condition.


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
3. Creates FdmSinhConcentrating1dMesher with xCenter = log(process->x0())
   (spot in log-space), alpha = sinhAlpha, and alignment targets = log of each
   barrier target
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
- FdmBlackScholesMesher::createSinhMesher must not modify the class's existing
  constructor or any other method
- all.hpp additions must be in alphabetical order within existing includes
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
□ Barrier align targets are converted to log-space before passing to sinh mesher
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
2. If it exists, apply changes as described below
3. If the path differs, adapt the changes to the actual file location

OUTPUT: 2 files (modified engine header + implementation).
</task>

<source_files>
FILE: ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp
[INSERT FULL CONTENT — or if unavailable, insert the description from
hpp_structure §2.8 and note that reconstruction is needed]

FILE: ql/pricingengines/vanilla/fdblackscholesvanillaengine.cpp
[INSERT FULL CONTENT if available]

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
      Real sinhAlpha = 0.0);             // NEW: 0 = standard mesher

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
- If the engine file path differs from expected, clearly state the actual path
  and adapt all references
</constraints>

<output_specification>
Produce 2 files:
1. fdblackscholesvanillaengine.hpp (modified)
2. fdblackscholesvanillaengine.cpp (modified)

If the file cannot be located, instead produce:
1. A description of what changes would be needed, with exact code snippets
2. A standalone helper function that constructs the full solver chain with
   all improvements enabled, bypassing the engine
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
</task>

<source_files>
Provide ALL output headers from Rounds 1-6 so the test can compile:
[INSERT all .hpp files from Rounds 1-6]

Also provide the QuantLib test infrastructure headers:
[INSERT any available test-suite/*.hpp or test patterns showing how QuantLib
registers test cases — if using Boost.Test, show the BOOST_AUTO_TEST_SUITE
pattern; if using QuantLib's custom framework, show that]
</source_files>

<implementation_guidance>

File: test-suite/fdmimprovedcn.cpp
(or adapt naming to match QuantLib's test convention — check existing test
files for the pattern)

The file should define a test suite with the following test cases:

=== T1: European Call Convergence ===
Params: S0=100, K=100, T=1, r=0.05, q=0, σ=0.2
Reference: Black-Scholes closed-form (use QuantLib's BlackScholesCalculator
or AnalyticEuropeanEngine for the reference value)
Procedure:
  - Run improved CN engine at (J,N) = (200,200), (400,400), (800,800)
  - Compare price against BS at each level
  - Compute Richardson order estimate p = log2(|e1|/|e2|) where e1 = v_200-v_400,
    e2 = v_400-v_800
  - Assert p ∈ [1.5, 2.5] for price
  - Assert |price_800 - BS| < 1e-4
  - Assert |delta_800 - BS_delta| < 1e-3
  - Assert |gamma_800 - BS_gamma| < 1e-2

=== T2: Digital Put Positivity ===
Params: S0=100, K=100, T=1, r=0.05, q=0, σ=0.2, cash-or-nothing put payoff
Procedure:
  - Run with fitted operator, 4 damping half-steps, xGrid=1200, tGrid=1200
  - Extract final solution array
  - Assert: min(u) >= -1e-10 (no significant negativity)
  - Assert: FdmDiagnostics::oscillationScore(u) < 0.02

=== T3: Low-Volatility M-Matrix ===
Params: S0=100, K=100, T=0.5, r=0.05, q=0, σ=0.02
Procedure:
  - Run with STANDARD operator (useFittedOperator=false), check for oscillations
  - Run with FITTED operator (useFittedOperator=true)
  - Assert fitted: mMatrixSatisfied() == true (from operator diagnostics)
  - Assert fitted: no negative values
  - Assert both converge to BS closed-form as grid refines

=== T4: Discretely Monitored Double Knock-Out ===
Params: S0=100, K=100, T=0.5, r=0.05, q=0, σ=0.25, L=95, U=110,
        5 equally-spaced monitoring dates
Procedure:
  - Run with monitoring-restart damping (monitoringDampingSteps=2),
    fitted operator, sinh mesh
  - Assert no negative values after solve
  - Compare against a fine-grid reference (same engine at 4x resolution)
  - Assert convergence (fine price ≈ coarse price within tolerance)

=== T5: Down-and-Out Rebate Call ===
Params: S0=50, K=40, T=0.5, r=0.04, q=0, σ=0.3, barrier B=20, rebate=2.5
Reference: QuantLib's AnalyticBarrierEngine (continuous barrier)
Procedure:
  - Run with barrier-aligned grid (barrier at x_min)
  - Assert |FD_price - analytic_price| < 1e-3

=== T6: American Put Policy Iteration ===
Params: S0=100, K=100, T=1, r=0.05, q=0, σ=0.2
Reference: fine binomial tree (CoxRossRubinstein, 10000 steps)
Procedure:
  - This test validates the policy iteration LCP solver standalone
  - Build a simple implicit Euler system for one time step
  - Apply policy iteration with American put payoff
  - Assert complementarity: u[j] >= payoff[j] for all j
  - Assert |FD_price - binomial_price| < 0.01

=== T7: Far-Field Boundary Sensitivity ===
Params: S0=300, K=100, T=5, r=0.03, q=0, σ=0.2 (deep ITM call)
Procedure:
  - Run with domain width m=6, m=8, m=10 standard deviations
  - Assert |price_m10 - price_m8| < 1e-6

HELPER FUNCTIONS to implement:
- BSClosedForm(S0, K, T, r, q, sigma) → (price, delta, gamma)
  using QuantLib's BlackScholesCalculator
- A function to extract the solution array from the engine (may need to
  access the internal solver, or run Fdm1DimSolver directly)

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

[... abbreviated — include bsmoperator.hpp, mixedscheme.hpp,
finitedifferencemodel.hpp, tridiagonaloperator.hpp as source files ...]

Key guidance for MixedScheme modification:
- Add a default no-op notifyDiscontinuity() method
- This is a template class, so changes go in the header only
- Since CrankNicolson<> inherits from MixedScheme<>, it automatically gets
  the method
- For FiniteDifferenceModel::rollbackImpl(), add a call to
  evolver_.notifyDiscontinuity() after each stopping time application
- Since all evolvers now have the method (as no-op default), no SFINAE needed
```

---

## Usage Instructions

**Context management:** Each round's `[INSERT ...]` placeholder should be filled with the actual file content from the QuantLib v1.23 source tree. Include only the files listed — adding unnecessary files wastes context and can confuse the model.

**Carrying forward outputs:** Rounds 4-7 reference outputs from prior rounds. When executing Round 4, include the Round 1 and Round 3 output *headers* (not implementations) in the `<source_files>` section. When executing Round 7, include ALL output headers from Rounds 1-6.

**Verification between rounds:** After each round, verify the output compiles against QuantLib v1.23 before proceeding. Common issues to check: missing `#include` directives, wrong `ext::shared_ptr` vs `std::shared_ptr`, missing `override` keywords, `Disposable<>` vs plain return types.

**If a round fails:** Re-run with the compilation errors appended to the prompt inside an `<error_context>` tag. Anthropic's docs recommend providing error feedback explicitly rather than asking the model to guess what went wrong.

**FdmSchemeDesc const-ness:** Before executing Round 3, check the actual `FdmSchemeDesc` definition in your copy of `fdmbackwardsolver.hpp`. If members ARE const, add a note to the Round 3 prompt: "Members are declared const. Use approach (b): add a new constructor overload." If they are NOT const, add: "Members are non-const. Simply add the new field with a default initializer."

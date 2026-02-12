## Preflight Report

**Status: Partially answerable. Several questions can be resolved from web sources (Implementing QuantLib blog, GitHub tests); others require the actual `.cpp` source files for full confidence. The actual QuantLib v1.23 `.hpp`/`.cpp` files were NOT provided in this conversation — only architectural planning documents.**

---

### Question 1: What does the 1D BS mesher store as spatial coordinate?

**Answer: `x = ln(S)` (natural log of spot), NOT `ln(S/K)`.**

**Evidence (from web sources, not provided files):**

(a) The `GeneralizedBlackScholesProcess` documentation states it governs `d ln S(t)` and warns that "the internal calculations work on ln S".[[5]](https://github.com/lballabio/QuantLib/blob/master/ql/processes/blackscholesprocess.hpp)

(b) The `FdBlackScholesVanillaEngine::calculate()` creates the inner value calculator as `FdmLogInnerValue(payoff, mesher, 0)`, which recovers `S` via `exp(x)`[[2]](https://www.implementingquantlib.com/2016/09/chapter-7-new-framework.html) — this implies the mesher coordinate is `ln(S)`, since `exp(ln(S)) = S`.

(c) The `FdmBlackScholesMesher` constructor takes `strike` as a **concentration point** for grid clustering, not as a coordinate origin.

**⚠️ Residual uncertainty:** I do not have the actual `FdmBlackScholesMesher` constructor body. To be 100% certain, I need:
- **`ql/methods/finitedifferences/meshers/fdmblackscholesmesher.cpp`**

**Impact on implementation:** Barrier alignment targets must be expressed as `ln(B)`, not `ln(B/K)`. The Greeks formulas `Δ = u_x/S`, `Γ = (u_{xx} − u_x)/S²` are valid for any log-coordinate, but any S-to-x conversions in the fitted operator must use `ln(S)`.

---

### Question 2: How do FirstDerivativeOp / SecondDerivativeOp incorporate nonuniform mesh spacing?

**Answer: They use the mesher's `dplus` and `dminus` arrays for local spacing.**

**Evidence:**

(a) `FirstDerivativeOp` and `SecondDerivativeOp` both inherit from `TripleBandLinearOp` and take `(Size direction, const shared_ptr<FdmMesher>& mesher)` in their constructors.[[3]](https://www.implementingquantlib.com/2016/10/chapter-8-fd-basic-operators.html)

(b) The mesher base class stores `locations`, `dplus`, and `dminus` vectors[[2]](https://quantlibjl.readthedocs.io/en/latest/methods.html) (confirmed by the Julia port and QuantLib test code). The test code (search result 33-1) explicitly checks `mesher.dminus()` and `mesher.dplus()`.

(c) The nonuniform spacing formulas are from Bowen & Smith (2005), "Derivative formulae and errors for non-uniformly spaced points."[[3]](https://www.implementingquantlib.com/2016/10/chapter-8-fd-basic-operators.html)

**⚠️ Residual uncertainty:** The exact stencil formulas (whether standard 3-point nonuniform or symmetric 3-point) require the actual implementation files:
- **`ql/methods/finitedifferences/operators/firstderivativeop.cpp`**
- **`ql/methods/finitedifferences/operators/secondderivativeop.cpp`**

**Impact on implementation:** For the fitted operator, when computing the local Péclet parameter `θ_j = b_j·h_j/(2·a_j)`, the effective spacing `h_j` should be derived from the mesher's `dplus` and `dminus` values (e.g., `h_j = 0.5*(dplus_j + dminus_j)`) to match how the standard operators define their stencils.

---

### Question 3: Does FiniteDifferenceModel call `applyTo` at every time step, or only at stopping times?

**Answer: At EVERY time step. Also at stopping times (with adjusted sub-stepping). Also at the initial `from` time if it equals a stopping time.**

**Evidence (definitive — full `rollbackImpl` source found):**

The `rollbackImpl` loop shows that in the non-hit branch (`else`), after `evolver_.step(a,now)`, the method calls `condition->applyTo(a, next)` unconditionally at every step.[[1]](https://www.implementingquantlib.com/2015/08/chapter-8-part-5-of-n-finite-difference.html) When a stopping time is hit, the code subdivides the step around the stopping time and calls `applyTo` at the stopping time AND at the end of the subdivided step. At the very start, if the `from` time equals a stopping time, `applyTo` is called before the loop begins.

**Impact on implementation:** This means `FdmAmericanStepCondition::applyTo` (which does a pointwise max with the payoff) runs at every single time step, not just at stopping times. For the monitoring projection step condition we want to create, we must check inside `applyTo` whether the current time matches a monitoring date and only project at those times (do nothing otherwise). This is a lightweight check.

---

### Question 4: Does FiniteDifferenceModel copy the evolver by value?

**Answer: Yes, by value.**

**Evidence:**

The `FiniteDifferenceModel` class template stores `Evolver evolver_` as a private member (not a pointer or reference).[[1]](https://www.implementingquantlib.com/2015/08/chapter-8-part-5-of-n-finite-difference.html) The constructor takes `const Evolver& evolver` and copies it into the member.

**Impact on implementation:** If we add state to `CrankNicolson<>` or `MixedScheme<>` (like damping counters), that state will be copied into the model. Modifications to the evolver after constructing the model will NOT be visible to the model. This means the `notifyDiscontinuity()` approach must work on the COPY inside the model. The existing `rollbackImpl` code uses `evolver_` (the copy) directly, so calling `evolver_.notifyDiscontinuity()` from within `rollbackImpl` (after applying a step condition) will work correctly on the internal copy.

---

### Question 5: Are FdmSchemeDesc fields const or non-const?

**Answer: `const`.**

**Evidence:**

The `FdmSchemeDesc` struct declares `const FdmSchemeType type; const Real theta, mu;` — all three data fields are `const`.[[2]](https://www.implementingquantlib.com/2017/04/fd-solvers.html)

**Impact on implementation:** We CANNOT add a mutable field to `FdmSchemeDesc` and set it after construction. Any new field (e.g., `monitoringDampingSteps`) must also be `const` and set in the constructor. This means we need either a new constructor overload or a new static factory method. The existing `FdmSchemeDesc(FdmSchemeType type, Real theta, Real mu)` constructor must be augmented. Because all fields are const, the struct is effectively immutable — which is clean but requires updating all construction sites.

Note: The enum in the blog excerpt lists `HundsdorferType, DouglasType, CraigSneydType, ModifiedCraigSneydType, ImplicitEulerType, ExplicitEulerType`[[2]](https://www.implementingquantlib.com/2017/04/fd-solvers.html) — but this may be from an older version. QuantLib v1.23 likely also has `CrankNicolsonType`, `TrBDF2Type`, and `MethodOfLinesType`. The actual v1.23 header is needed:
- **`ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp`**

---

### Question 6: Sign conventions for M-matrix checks in `mapT_`

**Answer: For the BS operator stored in `mapT_`, the M-matrix condition requires `lower_[i] ≤ 0` and `upper_[i] ≤ 0` (nonpositive off-diagonals), with `diag_[i] > 0`.**

**Evidence and derivation:**

The `FdmBlackScholesOp` stores in `mapT_` the operator: `L_BS(t) = -(r−q−σ²/2)·∂/∂x − (σ²/2)·∂²/∂x² + r`.[[1]](https://www.implementingquantlib.com/2016/12/chapter-8-fd-operator-examples.html) Note the signs: the convection and diffusion terms carry negative signs; the reaction term `+r` is positive.

The `FdmBackwardSolver::rollback()` uses `ImplicitEulerScheme` for damping steps[[2]](https://www.implementingquantlib.com/2017/04/fd-solvers.html), which solves the system `(I + dt·mapT_)·u_new = u_old` (going backward from later to earlier calendar time with `dt > 0`).

The off-diagonals of the system matrix `(I + dt·mapT_)` are `dt·mapT_.lower[i]` and `dt·mapT_.upper[i]`. For M-matrix structure (nonpositive off-diagonals), we need:

- **`mapT_.lower[i] ≤ 0`** for all interior `i`
- **`mapT_.upper[i] ≤ 0`** for all interior `i`

For the standard centered-difference operator with `b = r−q−σ²/2`, `a = σ²/2`:

`mapT_.lower = b/(2h) − a/h²` — violates ≤ 0 when `bh/(2a) > 1` (Péclet > 1)

`mapT_.upper = −b/(2h) − a/h²` — violates ≤ 0 when `−bh/(2a) > 1` (Péclet < −1)

For the **fitted** operator with fitting factor `ρ`, the math framework's `L_{math}` has nonneg off-diags, and since `mapT_ = −L_{math}`, the fitted `mapT_` has nonpositive off-diags (**ℓ_j ≥ 0 in math framework ↔ mapT_.lower ≤ 0 in QuantLib**). This is the desired M-matrix property.

**⚠️ Critical implementation note:** The sign convention in `mapT_` is OPPOSITE to the math framework's `L_{math}`. When writing `assembleFittedOperator`, the coefficient stored in `mapT_.lower[i]` should be `−(a_j·ρ_j/h² − b_j/(2h))`, and `mapT_.upper[i]` should be `−(a_j·ρ_j/h² + b_j/(2h))`. The diagnostics check should verify `lower_[i] ≤ 0` and `upper_[i] ≤ 0`.

---

### Files Needed for Round 1 (Code Production)

To proceed with implementation, I need the following actual QuantLib v1.23 source files:

1. **`ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp`** and **`.cpp`** — confirm coordinate convention, constructor signature
2. **`ql/methods/finitedifferences/operators/fdmblackscholesop.hpp`** and **`.cpp`** — confirm `setTime`, `mapT_` assembly, `axpyb` usage, full `FdmLinearOpComposite` interface
3. **`ql/methods/finitedifferences/operators/triplebandlinearop.hpp`** — confirm `axpyb` signature, `lower_/diag_/upper_` access patterns, `solve_splitting` parameters
4. **`ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp`** and **`.cpp`** — confirm `FdmSchemeDesc` struct (full v1.23 enum), `rollback` implementation
5. **`ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp`** and **`.cpp`** — confirm constructor, `step()` body, interaction with `ImplicitEulerScheme`
6. **`ql/methods/finitedifferences/schemes/impliciteulerscheme.hpp`** — confirm `step()` and `setStep()` signatures
7. **`ql/methods/finitedifferences/operators/firstderivativeop.cpp`** and **`secondderivativeop.cpp`** — confirm nonuniform stencil formulas
8. **`ql/methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp`** — confirm `stoppingTimes()` accessor and condition application pattern
9. **`ql/methods/finitedifferences/solvers/fdm1dimsolver.hpp`** and **`.cpp`** — confirm Greeks interpolation interface
10. **`ql/methods/finitedifferences/utilities/fdminnervaluecalculator.hpp`** — confirm `FdmCellAveragingInnerValue` interface and quadrature

Please provide these files for Round 1.

---
Learn more:
1. [Chapter 8, part 11: Black-Scholes finite-difference operators — Implementing QuantLib](https://www.implementingquantlib.com/2016/12/chapter-8-fd-operator-examples.html)
2. [Chapter 8, part 5 of n: finite-difference models — Implementing QuantLib](https://www.implementingquantlib.com/2015/08/chapter-8-part-5-of-n-finite-difference.html)
3. [Chapter 8, part 7 of n: the new finite-difference framework — Implementing QuantLib](https://www.implementingquantlib.com/2016/09/chapter-7-new-framework.html)
4. [QuantLib Documentation](https://www.quantlib.org/docs.shtml)
5. [Chapter 8, part 13: finite-difference schemes and solvers — Implementing QuantLib](https://www.implementingquantlib.com/2017/04/fd-solvers.html)
6. [Pricing Methods — QuantLib.jl 0.0.1 documentation](https://quantlibjl.readthedocs.io/en/latest/methods.html)
7. [QuantLib: ql/pricingengines/barrier/fdhestonbarrierengine.cpp Source File](https://rkapl123.github.io/QLAnnotatedSource/d9/d30/fdhestonbarrierengine_8cpp_source.html)
8. [quantlib/QuantLib/ql/methods/finitedifferences/trbdf2.hpp at master · MattPD/quantlib](https://github.com/MattPD/quantlib/blob/master/QuantLib/ql/methods/finitedifferences/trbdf2.hpp)
9. [Chapter 8, part 10: basic finite-difference operators — Implementing QuantLib](https://www.implementingquantlib.com/2016/10/chapter-8-fd-basic-operators.html)
10. [QuantLib/fdmlinearop.cpp at master · lballabio/QuantLib](https://github.com/lballabio/QuantLib/blob/master/test-suite/fdmlinearop.cpp)
11. [FdmStepConditionComposite | quantlib.js](https://quantlib.js.org/docs/classes/_ql_methods_finitedifferences_stepconditions_fdmstepconditioncomposite_.fdmstepconditioncomposite.html)
12. [QuantLib: DMinus Class Reference](https://rkapl123.github.io/QLAnnotatedSource/dd/de6/class_quant_lib_1_1_d_minus.html)
13. [Stochastic Processes — QuantLib-Python Documentation 1.39 documentation](https://quantlib-python-docs.readthedocs.io/en/latest/stochastic_processes.html)
14. [QuantLib, a free/open-source library for quantitative finance](https://www.quantlib.org/)
15. [QuantLib/ql/processes/blackscholesprocess.hpp at master · lballabio/QuantLib](https://github.com/lballabio/QuantLib/blob/master/ql/processes/blackscholesprocess.hpp)
16. [Pricing Engines — QuantLib.jl 0.0.1 documentation](https://quantlibjl.readthedocs.io/en/latest/pricing_engines.html)
17. [QuantLib: C:/dev/QuantLib/ql/experimental/finitedifferences/fdklugeextouspreadengine.hpp Source File](https://rkapl123.github.io/QLAnnotatedSource/df/d40/fdklugeextouspreadengine_8hpp_source.html)
18. [Package ‘RQuantLib’ July 21, 2025 Title R Interface to the 'QuantLib' Library](https://cran.r-project.org/web/packages/RQuantLib/RQuantLib.pdf)
19. [Math Tools — QuantLib-Python Documentation 1.40 documentation](https://quantlib-python-docs.readthedocs.io/en/latest/mathTools.html)
20. [Quantlib libray en Python for Options | by Paulino Seoane | Medium](https://medium.com/@pjseoane/quantlib-libray-en-python-for-options-ff557d565e74)
21. [QuantLib - Wikipedia](https://en.wikipedia.org/wiki/QuantLib)
22. [Releases · lballabio/QuantLib](https://github.com/lballabio/QuantLib/releases)
23. [Read QuantLib Python Cookbook | Leanpub](https://leanpub.com/quantlibpythoncookbook/read)
24. [Option Model Handbook, Part III: European Option Pricing With QuantLib Python - G B](http://gouthamanbalaraman.com/blog/european-option-binomial-tree-quantlib-python.html)
25. [QuantLib download | SourceForge.net](https://sourceforge.net/projects/quantlib/)
26. [Using QuantLib interactively — Implementing QuantLib](https://www.implementingquantlib.com/2024/03/using-quantlib-interactively.html)
27. [QuantLib: Interpolation Class Reference](https://rkapl123.github.io/QLAnnotatedSource/d9/dbc/class_quant_lib_1_1_interpolation.html)
28. [Introduction to QuantLib Python - G B](http://gouthamanbalaraman.com/blog/quantlib-basics.html)
29. [GitHub - lballabio/QuantLib: The QuantLib C++ library](https://github.com/lballabio/QuantLib)
30. [QuantLib: Introduction](https://rkapl123.github.io/QLAnnotatedSource/index.html)
31. [Chapter 8, part 1 of n: the finite-differences framework. — Implementing QuantLib](https://www.implementingquantlib.com/2015/06/chapter-8-part-1-of-n-finite.html)
32. [Tutorials and Examples — Implementing QuantLib](https://www.implementingquantlib.com/p/tutorials.html)
33. [QuantLib · PyPI](https://pypi.org/project/QuantLib/)
34. [QuantLib/ql/methods/finitedifferences at master · lballabio/QuantLib](https://github.com/lballabio/QuantLib/tree/master/ql/methods/finitedifferences)
35. [QuantLib-Risks: Enhance Derivative Analysis - XAD Automatic Differentiation](https://auto-differentiation.github.io/quantlib-risks/practical-guide/)
36. [Basics — QuantLib-Python Documentation 1.40 documentation](https://quantlib-python-docs.readthedocs.io/en/latest/basics.html)
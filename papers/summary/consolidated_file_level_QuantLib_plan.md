*Thinking...*

> The user wants me to produce a comprehensive, self-contained document titled "QuantLib v1.23 — Improved CN Implementation: File-Level Modification Plan" that bridges the mathematical framework from `consolidated_summary.md` with the QuantLib v1.23 source code structure from `hpp_structure_of_quantlib.md`, guided by the software architecture from `consolidated_code.md` and `consolidated_code_alternative.md`.
> 
> This is a massive document. Let me carefully work through each section:
> 
> ## 0. Conventions & Glossary
> - Map M1-M12 to QuantLib directories
> - Define notation
> 
> ## 1. New Files to Create
> 
> 
> - Identify files that don't exist in QuantLib v1.23 but are needed
> 
> ## 2. Files to Modify (by QuantLib subdirectory)
> - Go through each subdirectory and specify changes
> 
> ## 3. Integration Wiring
> - End-to-end user flow
> 
> ## 4. Build System & Header Registration
> 
> ## 5. Test Scaffolding
> 
> ## 6. Dependency & Risk Summary
> 
> Now I need to examine what actually exists in QuantLib v1.23 and determine what needs to be created or modified. Looking at the finite difference framework paths, I'm focusing on the `ql/methods/finitedifference/` directory structure and identifying which new Fdm* classes and supporting infrastructure need to be added versus what can be integrated into existing files. Looking at the legacy framework components like CrankNicolson and MixedScheme templates, along with the FiniteDifferenceModel and supporting operators, then moving into the pricing engines that implement these methods for vanilla options.
> 
> Now I'm mapping the modular components M1 through M5 to their QuantLib equivalents—the mesh factories, coordinate transforms embedded in the mesher implementations, inner value calculators for payoff processing, volatility surface adapters tied to the BlackScholesProcess, and continuing through the remaining modules.
> 
> Continuing with M5 through M7, the spatial operators are spread across multiple files handling finite difference schemes and boundary conditions, with legacy implementations still present alongside the newer modular structure. The time-stepping logic is implemented through various schemes like Crank-Nicolson and implicit Euler, again with both current and legacy versions coexisting in the codebase.
> 
> Now I'm mapping out the remaining components: the tridiagonal solver embedded within the triple band operator, the American option constraint handling through step conditions, the Greeks calculation engine that interpolates and computes derivatives, a diagnostics module that doesn't currently exist in QuantLib, and the orchestrator that ties everything together in the backward solver.
> 
> For the new implementations, I need to consider whether to create an exponentially fitted spatial operator as a new class or modify the existing Black-Scholes operator, build out a dedicated diagnostics utility, potentially add policy iteration for the linear complementarity problem in American options, and decide if the sinh grid mesher should leverage the existing concentrating mesher or be built from scratch.
> 
> Looking at the existing QuantLib structure more carefully, I see there's both a legacy TRBDF2 framework in `trbdf2.hpp` and a newer scheme-based version in `schemes/trbdf2scheme.hpp`, so I need to be mindful of which framework I'm working within.
> 
> Now checking what scheme types are already supported in `FdmSchemeDesc` - it has an enum that includes HundsdorferType, DouglasType, CraigSneydType, and importantly, TrBDF2Type is already listed there, which means the infrastructure for this scheme is partially in place.
> 
> Looking at how `FdmBackwardSolver` implements Rannacher-style damping, it runs initial implicit Euler steps before switching to the main scheme, and I can see the damping model being set up with the appropriate finite difference model.
> 
> So QuantLib v1.23 already has the startup damping mechanism, but it's missing several key features: restarting damping based on solution monitoring, an exponentially fitted spatial operator to maintain the M-matrix property, a diagnostics engine for detecting positivity violations and oscillations, policy iteration for American option constraints instead of just pointwise maximum operations, and sinh-graded mesh refinement.
> 
> Now I need to focus on how to implement the fitted spatial operator. The critical piece is modifying `FdmBlackScholesOp::setTime()` where the finite difference stencils are assembled—I'll need to replace the standard centered difference coefficients with fitted ones that preserve monotonicity. The operator currently uses separate first and second derivative maps with standard discretization, so I need to understand how those are constructed and where to inject the fitting logic. I'm considering whether to create a dedicated `FdmFittedBlackScholesOp` class that incorporates fitted derivative operators, or modify the existing `setTime()` method to apply the fitting adjustments directly. The fitting factor depends on the local grid spacing and the PDE coefficients—specifically the diffusion and convection terms—and transforms the finite difference stencil coefficients accordingly. Option A feels cleaner since it avoids modifying existing code paths. The real constraint here is that fitting fundamentally changes the combined operator—the fitting factor ρ scales the diffusion coefficient and shifts the diffusion-convection balance—so I need to construct the tridiagonal coefficients directly rather than layering separate derivative operators. For the monitoring discontinuity restart, I'm looking at how `FdmBackwardSolver::rollback()` applies damping at the start but lacks a restart mechanism after monitoring events, and I need to trace through `FiniteDifferenceModel::rollbackImpl()` to understand where that restart logic should fit.
> 
> I think modifying `FdmBackwardSolver::rollback()` to split the backward evolution at monitoring times and reapply damping is cleaner than adding a callback to the CN scheme, since the solver already manages the damping-to-main-stepping transition. Now I'm ready to structure the implementation—I'll need new operator classes for the fitted Black-Scholes dynamics and supporting infrastructure.
> 
> For the sinh mesher, I should verify whether the existing `Concentrating1dMesher` already implements the Stage2 specification or if I need to create a specialized version that matches the exact formula.
> 
> On policy iteration, the current `FdmAmericanStepCondition` uses a simple operator splitting approach—solving the PDE then applying the exercise constraint pointwise. Policy iteration would be more robust but requires solving modified tridiagonal systems instead, so I'll need to build a new implementation for that.
> 
> The algorithm iterates between determining which nodes are in the exercise region and solving a modified linear system where active rows become identity equations with payoff values as the right-hand side. The cleanest approach in QuantLib would be creating a dedicated policy iteration class that wraps the linear complementarity problem and integrates into the backward solver. of the simple step condition
> 
> Actually, let me reconsider. The existing architecture has the step condition applied *after* the linear solve. Policy iteration needs to modify the linear system *during* the solve. So we need a different integration point.
> 
> For the new framework, the `ImplicitEulerScheme::step()` calls `map_->solve_splitting()` for 1D. We could create a wrapper that intercepts this solve and applies policy iteration.
> 
> Alternatively, we could create a new `FdmAmericanPolicyIterationSolver` that replaces `Fdm1DimSolver` for American options.
> 
> I think the cleanest approach within QuantLib's architecture is:
> 1. Keep `FdmAmericanStepCondition` for the simple projection approach
> 2. Add a new utility class `FdmPolicyIterationConstraint` that can be used within the backward solver
> 3. Modify `FdmBackwardSolver` to optionally use the policy iteration approach
> 
> Now I need to write out the full document. I'm also thinking through what enum values should be added to `FdmSchemeDesc::FdmSchemeType`. The current options include various schemes like Hundsdorfer, Douglas, Craig-Sneyd variants, implicit and explicit Euler, method of lines, TR-BDF2, and Crank-Nicolson. I might want to add a Rannacher CN variant that handles explicit monitoring-restart damping differently than the current approach, though I need to reconsider how the existing Crank-Nicolson type fits into this.
> 
> The current `rollback()` method in `FdmBackwardSolver` applies all damping steps first using implicit Euler, then runs the main scheme steps. But for monitoring restarts, we'd need to interleave damping and main stepping, which would require restructuring how the method works.
> 
> I can leverage the fact that `FdmStepConditionComposite` already tracks `stoppingTimes()` which includes monitoring times. The cleanest approach would be to modify `rollback()` to accept these discontinuity times, split the backward evolution at each one, and reapply damping after each discontinuity rather than just at the beginning.
> 
> Looking at the current `rollback()` signature, I could keep it as-is but have it internally consult the condition's stopping times to handle the damping restart logic. For the legacy framework, I need to check how `CrankNicolson<TridiagonalOperator>` fits into this picture.
> 
> The legacy path uses `FiniteDifferenceModel` with `CrankNicolson` and `BSMOperator`, and since `FiniteDifferenceModel` already respects stopping times and `AmericanCondition` handles projection, I'd mainly need to create a fitted operator variant and potentially adjust the first few time steps to use implicit Euler for damping.
> 
> Now I'm planning out the file structure—I'll need new headers and implementations for the fitted Black-Scholes operator in both the new framework and legacy framework, plus a diagnostics utility module. Now I'm mapping out the remaining modifications across the solver infrastructure—the backward solver needs the major damping and scheme dispatch work, the 1D solver should improve Greeks extraction, the Black-Scholes solver accepts a fitted operator option, and the descriptor gains diagnostics fields. The mesher integration brings in the sinh option, while the step conditions expose monitoring times for restart coordination. I'm realizing the scope is too broad and I need to narrow my focus. Let me reconsider what's actually required: a precise file-level modification plan that covers both the new Fdm framework and legacy MixedScheme, maps all M1-M12 modules, and includes exact paths, class names, and line regions. I should also verify which files actually exist in the codebase rather than inventing paths, and ensure I'm referencing the correct mathematical sections and specifying how policy iteration integrates with the American option pricing.
> 
> Now I'm checking the actual file structure from the attached QuantLib headers to confirm the paths in the finitedifferences directory. Continuing through the meshers directory with the remaining implementations for Black-Scholes variants, CEV processes, Heston variance, and composite meshing strategies, then moving into the operators directory to verify the finite difference method operator implementations for various pricing models. Continuing through the remaining operator headers, then moving into the schemes subdirectory where I'm verifying the boundary condition helpers and various finite difference schemes like Crank-Nicolson and Douglas methods. Continuing through the solvers directory with the remaining FDM implementations, then moving into the step conditions subdirectory where I'm verifying the American, Bermuda, and swing option conditions along with storage and averaging components. Continuing through the utilities directory, I'm verifying the remaining FDM-related headers for inner value calculations, boundary conditions, and various risk-neutral density calculators.
> 
> I notice the pricing engine headers aren't in the attached files, but the code tree documentation describes them in detail—they're located in `ql/pricingengines/` rather than the methods directory. Since only header files were provided and the actual implementations exist elsewhere in QuantLib, I'll work from the structural information already documented.
> 
> This response is going to be quite extensive given the user's request for a full, detailed document without abbreviations. I'll merge the module definitions from both consolidated code documents, preferring the alternative version where they diverge, and aim for comprehensive coverage of all implementable details.
> 
> I should note that some pricing engine files like the FD Black-Scholes vanilla engine are referenced from the code tree document but aren't in the attached materials, so I'll document them based on that analysis. Similarly, QuantLib's test-suite directory isn't included, but I can reference the standard naming conventions used there. Now I'll structure this carefully and begin writing the full document, estimating roughly 30 pages across all sections.
> 
> I'm realizing the key insight here: QuantLib v1.23 already has `CrankNicolsonScheme` wrapping explicit and implicit Euler schemes, and `FdmBackwardSolver` already implements Rannacher damping at startup. So the main additions I need to focus on are monitoring-restart damping modifications to `FdmBackwardSolver` and a new exponentially fitted operator for the Black-Scholes case.
> 
> I also need to enhance cell-averaged payoff handling with Gauss quadrature, add policy iteration for the linear complementarity problem, create diagnostics utilities, potentially add a sinh mesher, and ensure configuration propagates through engine constructors and solver descriptions. For the legacy framework, I'll need a fitted operator and damping restart capability in `FiniteDifferenceModel`.
> 
> Now I'm checking whether `FdmCellAveragingInnerValue` already exists in v1.23...
> 
> Checking the inner value calculator header... I can see that `FdmCellAveragingInnerValue` already implements cell averaging with `avgInnerValue()`, and the solver is already calling this method during initialization. The current implementation uses Simpson's rule for quadrature, but the math framework suggests 3-point Gauss-Legendre might be better—so the main question is whether we should enhance the quadrature method itself.
> 
> Looking at what's already in QuantLib v1.23, most of the cell-averaging infrastructure is already there. The real additions needed are just corner compatibility checks and optionally upgrading to Gauss-Legendre quadrature. I should figure out what actually needs to change versus what's already implemented.
> 
> So QuantLib has Rannacher damping, cell-averaged payoffs, log-price transforms, concentration meshers, TR-BDF2 schemes, American step conditions, and boundary handling all built in. What's missing is monitoring-restart damping after barrier projections and exponentially fitted spatial discretization.
> 
> Now I'm checking whether QuantLib already has corridor projection at monitoring times for barrier handling by looking at the Dirichlet boundary implementation.
> 
> I see that `FdmDirichletBoundary` sets values on one side, but for discrete barrier monitoring, I'd need something that zeros out nodes outside the corridor at specific times—essentially a step condition like `FdmBarrierProjectionCondition`. Looking through the existing step conditions like `FdmBermudanStepCondition` and `FdmSimpleSwingCondition`, I don't find a specific discrete barrier handler in v1.23, so this would need to be a custom implementation rather than relying on existing QuantLib infrastructure.
> 
> Now I'm ready to write the documentation. I'll be thorough but focused on the key implementation details. One thing I need to verify: the code tree mentions `fdvanillaengine.hpp` in the vanilla pricingengines folder, but I don't see it in the provided headers, so I'll mark that as needing verification.

# QuantLib v1.23 — Improved CN Implementation: File-Level Modification Plan

---

## 0. Conventions & Glossary

### 0.1 Path Convention

All file paths are relative to the `ql/` root directory of the QuantLib v1.23 source tree. For example, `methods/finitedifferences/schemes/cranknicolsonscheme.hpp` refers to `ql/methods/finitedifferences/schemes/cranknicolsonscheme.hpp`.

### 0.2 Framework Labels

QuantLib v1.23 contains two parallel finite-difference frameworks. Every file entry in this document is labeled with one of:

**[NEW-FW]** — The modern `Fdm*` class hierarchy rooted in `methods/finitedifferences/operators/`, `schemes/`, `solvers/`, `meshers/`, `stepconditions/`, and `utilities/`. This framework supports multi-dimensional problems, uses `FdmLinearOpComposite` as the operator interface, and dispatches time stepping through `FdmBackwardSolver` and `FdmSchemeDesc`.

**[LEGACY-FW]** — The older template-based framework using `MixedScheme<Operator>`, `CrankNicolson<TridiagonalOperator>`, `FiniteDifferenceModel<Evolver>`, `BSMOperator`, and `TridiagonalOperator`. This framework is restricted to 1D tridiagonal problems and is used by older pricing engines such as `FDVanillaEngine`.

### 0.3 Module-to-Directory Mapping

The software architecture defines modules M1–M12. Each maps to one or more QuantLib directories as follows.

**M1 (GridFactory)** maps to `methods/finitedifferences/meshers/`. The primary files are `fdmblackscholesmesher.hpp/.cpp`, `concentrating1dmesher.hpp/.cpp`, `uniform1dmesher.hpp/.cpp`, and `fdmmeshercomposite.hpp/.cpp` for the new framework, plus grid construction logic embedded in pricing engines for the legacy framework.

**M2 (CoordinateTransform)** has no dedicated QuantLib file. The log-price transform is implicit in `fdmblackscholesmesher.cpp` (which constructs the grid in log-space) and in how operators interpret the grid. Greeks conversion formulas live in `solvers/fdm1dimsolver.hpp/.cpp`.

**M3 (PayoffProcessor)** maps to `methods/finitedifferences/utilities/fdminnervaluecalculator.hpp/.cpp`. The class `FdmCellAveragingInnerValue` already provides Simpson-rule cell averaging via its `avgInnerValue()` method, used during initialization in `Fdm1DimSolver`.

**M4 (VolSurfaceAdapter)** maps to `processes/blackscholesprocess.hpp/.cpp` and the volatility term structure hierarchy under `termstructures/volatility/`. The operator classes call `process->blackVolatility()` or `process->localVolatility()`.

**M5 (SpatialOperator)** maps to `methods/finitedifferences/operators/` for the new framework (primarily `fdmblackscholesop.hpp/.cpp`, `triplebandlinearop.hpp/.cpp`, `firstderivativeop.hpp/.cpp`, `secondderivativeop.hpp/.cpp`) and to `methods/finitedifferences/bsmoperator.hpp/.cpp` for the legacy framework.

**M6 (BoundaryHandler)** maps to `methods/finitedifferences/utilities/fdmdirichletboundary.hpp/.cpp`, `fdmtimedepdirichletboundary.hpp/.cpp`, and `methods/finitedifferences/boundarycondition.hpp/.cpp` (legacy). Discrete monitoring corridor projection is a new capability requiring a new step condition.

**M7 (TimeStepper)** maps to `methods/finitedifferences/schemes/` (new framework: `cranknicolsonscheme.hpp/.cpp`, `impliciteulerscheme.hpp/.cpp`, `expliciteulerscheme.hpp/.cpp`, `trbdf2scheme.hpp`) and to `methods/finitedifferences/mixedscheme.hpp`, `cranknicolson.hpp`, `finitedifferencemodel.hpp` (legacy framework). The backward-solver driver is `methods/finitedifferences/solvers/fdmbackwardsolver.hpp/.cpp`.

**M8 (TridiagSolver)** maps to `methods/finitedifferences/operators/triplebandlinearop.hpp/.cpp` (the `solve_splitting` method implements Thomas) for the new framework, and to `methods/finitedifferences/tridiagonaloperator.hpp/.cpp` (`solveFor` method) for the legacy framework.

**M9 (AmericanConstraint)** maps to `methods/finitedifferences/stepconditions/fdmamericanstepcondition.hpp/.cpp` for the simple projection approach. Policy iteration requires a new utility class.

**M10 (GreeksEngine)** maps to `methods/finitedifferences/solvers/fdm1dimsolver.hpp/.cpp` (which provides `interpolateAt`, `derivativeX`, `derivativeXX`, `thetaAt`) and to the `fetchResults` methods of pricing engines.

**M11 (DiagnosticsEngine)** has no QuantLib counterpart. New files are required.

**M12 (Orchestrator)** maps to `methods/finitedifferences/solvers/fdmbackwardsolver.hpp/.cpp` for scheme dispatch and to `pricingengines/vanilla/fdblackscholesvanillaengine.hpp/.cpp` [VERIFY: exact path in ql/pricingengines/vanilla/] for user-facing configuration.

### 0.4 Mathematical Notation

Consistent with the math framework: `τ` is time-to-maturity (τ=0 at payoff, τ=T today), `x = ln(S/K)` is log-moneyness, `a(x,τ) = σ²/2` is the diffusion coefficient, `b(x,τ) = (r−q) − σ²/2` is the convection coefficient, `r` is risk-free rate, `q` is dividend yield. The fitting factor is `ρ_j = θ_j coth(θ_j)` where `θ_j = b_j h/(2 a_j)` and `h` is the local grid spacing.

### 0.5 Failure Mode Codes

F1 (convection dominance), F2 (nonsmooth payoff), F3 (monitoring resets), F4 (positivity violation), F5 (numerical diffusion), F6 (far-field truncation), F7 (American LCP), F8 (Greeks instability). These are defined in the math framework §1.3.1.

---

## 1. New Files to Create

### 1.1 — `methods/finitedifferences/operators/fdmfittedblackscholesop.hpp`

**Module mapping:** M5 (SpatialOperator) — exponentially fitted BS operator for the new framework.

**Mathematical reference:** Math framework §3.2.3 (fitting factor ρ), §3.2.4 (fitted stencil coefficients), §3.1.4 (log-price PDE coefficients).

**Register in:** `methods/finitedifferences/operators/all.hpp` (add `#include` line) and the corresponding `Makefile.am` / `CMakeLists.txt` for the `operators/` subdirectory.

**Why it is needed:** The existing `FdmBlackScholesOp` assembles the spatial operator by combining `FirstDerivativeOp` and `SecondDerivativeOp` with standard centered-difference stencils. These stencils fail the M-matrix sign condition when convection dominates (σ² ≪ r), causing oscillations (F1, F4). The fitted operator replaces the standard combination with node-by-node fitted coefficients that guarantee nonnegative off-diagonals and therefore M-matrix structure under appropriate conditions.

**Complete public interface:**

```cpp
namespace QuantLib {

    class FdmMesher;
    class GeneralizedBlackScholesProcess;
    class FdmQuantoHelper;

    class FdmFittedBlackScholesOp : public FdmLinearOpComposite {
      public:
        FdmFittedBlackScholesOp(
            const ext::shared_ptr<FdmMesher>& mesher,
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
            Real strike,
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            Size direction = 0,
            ext::shared_ptr<FdmQuantoHelper> quantoHelper
                = ext::shared_ptr<FdmQuantoHelper>(),
            bool enableMMatrixFallback = true);

        Size size() const override;
        void setTime(Time t1, Time t2) override;

        Disposable<Array> apply(const Array& r) const override;
        Disposable<Array> apply_mixed(const Array& r) const override;
        Disposable<Array> apply_direction(Size direction,
                                          const Array& r) const override;
        Disposable<Array> solve_splitting(Size direction,
                                          const Array& r, Real s) const override;
        Disposable<Array> preconditioner(const Array& r, Real s) const override;

        Disposable<std::vector<SparseMatrix>> toMatrixDecomp() const override;

        // Diagnostics accessors
        bool mMatrixSatisfied() const;
        Size fallbackNodeCount() const;

      private:
        const ext::shared_ptr<FdmMesher> mesher_;
        const ext::shared_ptr<YieldTermStructure> rTS_, qTS_;
        const ext::shared_ptr<BlackVolTermStructure> volTS_;
        const ext::shared_ptr<LocalVolTermStructure> localVol_;
        const Array x_;
        TripleBandLinearOp mapT_;
        const Real strike_;
        const Real illegalLocalVolOverwrite_;
        const Size direction_;
        const ext::shared_ptr<FdmQuantoHelper> quantoHelper_;
        const bool enableMMatrixFallback_;
        mutable bool mMatrixOk_;
        mutable Size fallbackCount_;

        void assembleFittedOperator(Time t1, Time t2);
        static Real fittingFactor(Real theta);
    };
}
```

**Implementation notes for `assembleFittedOperator`:** For each interior node `j`, compute `a_j = σ_j²/2`, `b_j = (r−q) − σ_j²/2`, `θ_j = b_j·h_j/(2·a_j)`, `ρ_j = fittingFactor(θ_j)`. Then write directly into `mapT_` arrays: `lower_[j] = a_j·ρ_j/h² − b_j/(2h)`, `diag_[j] = −2·a_j·ρ_j/h² − r`, `upper_[j] = a_j·ρ_j/h² + b_j/(2h)`. If `lower_[j] < 0` or `upper_[j] < 0` and `enableMMatrixFallback_` is true, recompute that node with pure upwind convection discretization and set `fallbackCount_++`. The static method `fittingFactor(θ)` returns `θ·coth(θ)` for `|θ|≥1e-8` and `1.0 + θ²/3` for `|θ|<1e-8`.

**Failure modes addressed:** F1 (M-matrix enforcement under convection dominance), F4 (positivity preservation), F5 (controlled numerical diffusion via fitting rather than crude upwinding).

---

### 1.2 — `methods/finitedifferences/operators/fdmfittedblackscholesop.cpp`

Implementation file for the class declared in §1.1. Contains `assembleFittedOperator`, `setTime`, `apply`, `solve_splitting`, `fittingFactor`, and the `toMatrixDecomp` method. The `solve_splitting` method delegates to `mapT_.solve_splitting()` which uses the Thomas algorithm already implemented in `TripleBandLinearOp`. Register in the same build targets as the header.

---

### 1.3 — `methods/finitedifferences/fittedbsmoperator.hpp`

**Module mapping:** M5 (SpatialOperator) — exponentially fitted BS operator for the legacy framework.

**Mathematical reference:** Math framework §3.2.3–§3.2.4; identical coefficient formulas to §1.1 but using `TridiagonalOperator` base class.

**Register in:** `methods/finitedifferences/all.hpp` and the root `Makefile.am`.

**Why it is needed:** The legacy `BSMOperator` uses centered differences (see `bsmoperator.hpp` which sets `pd`, `pm`, `pu` with standard stencils). For parity with the new framework's fitted operator, the legacy path needs a `FittedBSMOperator` that inherits from `TridiagonalOperator` and applies exponential fitting.

**Complete public interface:**

```cpp
namespace QuantLib {

    class FittedBSMOperator : public TridiagonalOperator {
      public:
        FittedBSMOperator(Size size, Real dx,
                          Rate r, Rate q, Volatility sigma,
                          bool enableFallbackUpwind = true);
        FittedBSMOperator(const Array& grid,
                          Rate r, Rate q, Volatility sigma,
                          bool enableFallbackUpwind = true);

        bool mMatrixSatisfied() const;
        Size fallbackNodeCount() const;

      private:
        bool mMatrixOk_;
        Size fallbackCount_;
        static Real fittingFactor(Real theta);
    };
}
```

**Implementation notes:** The constructor iterates over interior nodes, computes the same `θ_j`, `ρ_j` fitting as in §1.1, and calls `setMidRow(i, lower_i, diag_i, upper_i)` inherited from `TridiagonalOperator`. For boundary rows, uses `setFirstRow` and `setLastRow` with linear extrapolation coefficients matching the existing `BSMOperator` convention.

**Failure modes addressed:** F1, F4 (legacy framework path).

---

### 1.4 — `methods/finitedifferences/fittedbsmoperator.cpp`

Implementation of the class declared in §1.3. Register alongside the header.

---

### 1.5 — `methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp`

**Module mapping:** M6 (BoundaryHandler) — discrete monitoring corridor projection for the new framework.

**Mathematical reference:** Math framework §3.1.2 (corridor projection at monitoring times), Algorithm 1 lines 183–193.

**Register in:** `methods/finitedifferences/stepconditions/all.hpp` and the stepconditions `Makefile.am`.

**Why it is needed:** QuantLib v1.23 has `FdmDirichletBoundary` for continuous barrier conditions, but no step condition that zeroes out solution values outside a corridor `[L, U]` at discrete monitoring times. This is essential for discretely monitored double knock-out options (F3). The projection acts as a step condition applied at specific stopping times.

**Complete public interface:**

```cpp
namespace QuantLib {

    class FdmBarrierProjectionCondition : public StepCondition<Array> {
      public:
        FdmBarrierProjectionCondition(
            std::vector<Time> monitoringTimes,
            Real lowerBarrier,
            Real upperBarrier,
            ext::shared_ptr<FdmMesher> mesher,
            Size direction = 0);

        void applyTo(Array& a, Time t) const override;

        const std::vector<Time>& monitoringTimes() const;
        bool isMonitoringTime(Time t, Time tolerance = 1e-10) const;

      private:
        const std::vector<Time> monitoringTimes_;
        const Real lowerBarrier_, upperBarrier_;
        const ext::shared_ptr<FdmMesher> mesher_;
        const Size direction_;
        std::vector<Size> outsideIndices_;
    };
}
```

**Implementation notes:** The constructor precomputes `outsideIndices_` — the set of layout indices whose spatial coordinate falls outside `[L, U]`. In `applyTo`, if `t` matches a monitoring time (within tolerance), set `a[i] = 0` for all `i` in `outsideIndices_`. Otherwise do nothing. The `monitoringTimes()` accessor is used by `FdmStepConditionComposite` to register stopping times.

**Failure modes addressed:** F3 (discrete monitoring resets), F4 (maintains nonnegative values outside corridor as zero).

---

### 1.6 — `methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.cpp`

Implementation of the class declared in §1.5.

---

### 1.7 — `methods/finitedifferences/utilities/fdmpolicyiteration.hpp`

**Module mapping:** M9 (AmericanConstraint) — policy iteration LCP solver.

**Mathematical reference:** Software architecture §3.2 (Policy Iteration for American LCP).

**Register in:** `methods/finitedifferences/utilities/all.hpp` and the utilities `Makefile.am`.

**Why it is needed:** The existing `FdmAmericanStepCondition` applies a simple pointwise projection `a[i] = max(a[i], payoff[i])` after each time step (operator-splitting approach). Policy iteration solves the LCP more accurately by determining the active set (exercise region) and solving a modified linear system iteratively. This produces a more precise early-exercise boundary, particularly relevant for American options where Γ accuracy matters (F7, F8).

**Complete public interface:**

```cpp
namespace QuantLib {

    class FdmPolicyIterationLCP {
      public:
        FdmPolicyIterationLCP(
            Array exerciseValues,
            Size maxIterations = 50,
            Real tolerance = 1e-12);

        // Solve LCP: find u >= phi s.t. A u >= rhs where equality holds
        // on the continuation set. Takes a tridiagonal solve functor.
        Disposable<Array> solve(
            const TripleBandLinearOp& systemOp,
            const Array& rhs,
            Real scaleFactor,
            const ext::shared_ptr<FdmMesher>& mesher) const;

        // For legacy framework (tridiagonal)
        Disposable<Array> solve(
            const TridiagonalOperator& A,
            const Array& rhs) const;

        Size lastIterationCount() const;

      private:
        const Array exerciseValues_;
        const Size maxIterations_;
        const Real tolerance_;
        mutable Size lastIterations_;
    };
}
```

**Implementation notes:** The `solve` method implements the active-set iteration from the software architecture §3.2: start with an initial guess from the unconstrained solve, determine which nodes should exercise (`u[i] <= phi[i]`), build a modified system with identity rows at exercise nodes and payoff in the RHS, solve, check if the active set changed, and repeat until convergence. For the new framework, it modifies a copy of the `TripleBandLinearOp` coefficients; for the legacy framework, it modifies a copy of the `TridiagonalOperator`.

**Failure modes addressed:** F7 (robust American LCP), F8 (accurate exercise boundary for Greeks).

---

### 1.8 — `methods/finitedifferences/utilities/fdmpolicyiteration.cpp`

Implementation of the class declared in §1.7.

---

### 1.9 — `methods/finitedifferences/utilities/fdmdiagnostics.hpp`

**Module mapping:** M11 (DiagnosticsEngine) — runtime health monitoring.

**Mathematical reference:** Software architecture §3.7 (oscillation detector), §3.8 (Richardson estimator).

**Register in:** `methods/finitedifferences/utilities/all.hpp` and the utilities `Makefile.am`.

**Why it is needed:** QuantLib v1.23 has no built-in mechanism to detect positivity violations, oscillatory solutions, or M-matrix failures at runtime. During development and production monitoring, these diagnostics are essential for identifying when the numerical scheme is behaving poorly (all failure modes F1–F8).

**Complete public interface:**

```cpp
namespace QuantLib {

    struct FdmDiagnosticsReport {
        Real minValue;
        Size negativeCount;
        Real oscillationScore;
        std::vector<Size> oscillationIntervals;
        Size mMatrixViolationCount;
        Size fallbackUpwindCount;
        Size nanCount;
    };

    struct FdmConvergenceRecord {
        Real valueCoarse, valueFine, valueFiner;
        Real orderEstimate;
        Real richardsonExtrapolation;
    };

    class FdmDiagnostics {
      public:
        enum Level { Off, Light, Full };

        explicit FdmDiagnostics(Level level = Off);

        FdmDiagnosticsReport checkSolution(const Array& u) const;

        FdmDiagnosticsReport checkSolutionWithOperator(
            const Array& u,
            bool operatorMMatrixOk,
            Size operatorFallbackCount) const;

        static Real oscillationScore(const Array& u);
        static std::vector<Size> flagOscillationIntervals(const Array& u);

        static FdmConvergenceRecord richardsonEstimate(
            Real vCoarse, Real vFine, Real vFiner);

        static FdmConvergenceRecord richardsonEstimateAssumeOrder(
            Real vCoarse, Real vFine, Real assumedOrder);

        Level level() const;

      private:
        Level level_;
    };
}
```

**Implementation notes:** `checkSolution` scans the array for min value, negative count, NaN count (Light mode), and additionally runs the oscillation detector (Full mode). The oscillation detector counts sign changes in `Δu_j = u[j+1] − u[j]` and normalizes by grid size, following the algorithm in software architecture §3.7. The Richardson estimator implements §3.8: `p = log2(|v_h − v_{h/2}| / |v_{h/2} − v_{h/4}|)` and extrapolation `v* = v_{h/2} + (v_{h/2} − v_h)/(2^p − 1)`.

**Failure modes addressed:** All (F1–F8) via detection and reporting.

---

### 1.10 — `methods/finitedifferences/utilities/fdmdiagnostics.cpp`

Implementation of the class declared in §1.9.

---

### 1.11 — `methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.hpp`

**Module mapping:** M1 (GridFactory) — sinh-graded mesh with exact node alignment.

**Mathematical reference:** Software architecture §3.4 (Sinh-Graded Mesh Construction).

**Register in:** `methods/finitedifferences/meshers/all.hpp` and the meshers `Makefile.am`.

**Why it is needed:** The existing `Concentrating1dMesher` uses an `asinh`-based grading that works well in many cases, but the sinh-graded formulation with explicit barrier-node alignment (solving for ξ₀ via bisection so that barrier and strike coordinates land exactly on grid nodes) provides superior control for barrier option pricing (F3, F6). The new mesher implements the construction from software architecture §3.4 with an added alignment feature.

**Complete public interface:**

```cpp
namespace QuantLib {

    class FdmSinhConcentrating1dMesher : public Fdm1dMesher {
      public:
        FdmSinhConcentrating1dMesher(
            Real xMin, Real xMax, Size size,
            Real xCenter,
            Real alpha = 3.0,
            const std::vector<Real>& alignTargets
                = std::vector<Real>());

        Real alpha() const;
        Real xCenter() const;

      private:
        Real alpha_, xCenter_;
    };
}
```

**Implementation notes:** Follows the algorithm in software architecture §3.4 exactly: uniform `ξ_j ∈ [0,1]`, solve for `ξ₀` via bisection so that `x(0) = xMin` and `x(1) = xMax` with `x(ξ) = xCenter + c·sinh(α·(ξ − ξ₀))`. After construction, apply a global shift if any `alignTargets` entry is not within tolerance of a grid node (shift the whole grid to snap the nearest node to the target). Populate `locations_`, `dplus_`, `dminus_` arrays inherited from `Fdm1dMesher`.

**Failure modes addressed:** F6 (boundary alignment), F3 (barrier-node alignment for discrete monitoring).

---

### 1.12 — `methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.cpp`

Implementation of the class declared in §1.11.

---

## 2. Files to Modify (by QuantLib subdirectory)

### 2.1 `methods/finitedifferences/schemes/`

#### 2.1.1 — `methods/finitedifferences/schemes/cranknicolsonscheme.hpp`

**Module mapping:** M7 (TimeStepper) — monitoring-restart damping capability.

**Mathematical reference:** Math framework §3.3.1 (RS-CN step formulas, damping state machine); Algorithm 1 lines 63–128.

**What changes and why:** The existing `CrankNicolsonScheme` delegates to `ExplicitEulerScheme::step()` and `ImplicitEulerScheme::step()` with a fixed θ. It has no awareness of monitoring discontinuities that occur mid-solve. The modification adds a damping state machine: after a monitoring projection reintroduces a discontinuity, the scheme reverts to pure implicit Euler half-steps for a configurable number of sub-steps before resuming CN. This prevents the persistent oscillations that CN exhibits when stepping over a freshly introduced discontinuity (F2, F3).

**Specific code regions to modify:**

The class declaration (approximately lines 30–55 of the header) gains new private state members and a public notification method. The existing constructor signature gains an optional parameter.

**New/changed method signatures:**

```cpp
class CrankNicolsonScheme {
  public:
    CrankNicolsonScheme(
        Real theta,
        const ext::shared_ptr<FdmLinearOpComposite>& map,
        const bc_set& bcSet = bc_set(),
        Real relTol = 1e-8,
        ImplicitEulerScheme::SolverType solverType
            = ImplicitEulerScheme::BiCGstab,
        Size dampingHalfSteps = 0);           // NEW parameter

    void step(array_type& a, Time t);
    void setStep(Time dt);

    void notifyDiscontinuity();               // NEW — resets damping counter
    bool isDamping() const;                   // NEW — query damping state

    Size numberOfIterations() const;

  private:
    /* existing members unchanged */
    Size dampingHalfSteps_;                   // NEW
    Size dampingRemaining_;                   // NEW
    bool inDampingPhase_;                     // NEW
};
```

The `step` method body (in the `.cpp` file, approximately lines 25–40) changes: if `inDampingPhase_ && dampingRemaining_ > 0`, call `implicit_->step(a, t, 1.0)` with `dt_` set to half the nominal step (by calling `implicit_->setStep(dt_/2)` and stepping twice), then decrement `dampingRemaining_` by 2 and restore `implicit_->setStep(dt_)`. When `dampingRemaining_` reaches zero, set `inDampingPhase_ = false`. The `notifyDiscontinuity()` method sets `inDampingPhase_ = true` and `dampingRemaining_ = dampingHalfSteps_`.

**Failure modes addressed:** F2 (nonsmooth payoff re-injection), F3 (monitoring resets), F8 (Greeks stabilization after discontinuities).

**Downstream impact:** `FdmBackwardSolver` must call `notifyDiscontinuity()` after barrier projection events. The `FdmSchemeDesc` needs no enum change since `CrankNicolsonType` already exists; the damping parameter threads through the solver description.

---

#### 2.1.2 — `methods/finitedifferences/schemes/cranknicolsonscheme.cpp`

**Module mapping:** M7 (TimeStepper).

**Mathematical reference:** Math framework §3.3.1; Algorithm 1 lines 63–128.

**What changes and why:** The implementation of `step()` must incorporate the damping-phase conditional. The implementation of `notifyDiscontinuity()` must be added. The constructor must store the new parameter.

**Specific code regions to modify:**

The constructor body (approximately lines 10–20) must initialize the new members `dampingHalfSteps_`, `dampingRemaining_ = 0`, `inDampingPhase_ = false`.

The `step()` method body (approximately lines 25–40) is replaced with a conditional: when damping is active, perform two implicit Euler half-steps at `dt_/2` using `implicit_->setStep(dt_*0.5)` and `implicit_->step(a, t, 1.0)` called twice, then restore `implicit_->setStep(dt_)`. When not damping, execute the original body (explicit step with `1-θ`, implicit step with `θ`).

Add new method bodies for `notifyDiscontinuity()` and `isDamping()`.

**Failure modes addressed:** F2, F3, F8.

---

#### 2.1.3 — `methods/finitedifferences/schemes/impliciteulerscheme.hpp`

**Module mapping:** M7 (TimeStepper) — ensure half-step capability.

**Mathematical reference:** Math framework §3.3.1 (implicit Euler damping steps at k/2).

**What changes and why:** No interface change is strictly required; the existing `ImplicitEulerScheme` already supports arbitrary `dt` via `setStep()`. However, verify that calling `setStep()` multiple times within a single nominal step does not corrupt internal state. If the scheme caches operator decompositions that depend on `dt`, those caches must be invalidated on `setStep()` calls. Based on the existing code (lines ~30–50 of the `.cpp`), `setStep` only stores `dt_` and does not trigger expensive recomputation, so no modification is needed. Document this dependency.

**Specific code regions to modify:** None (verification only).

**Failure modes addressed:** F2, F3 (indirectly, by enabling damping sub-steps).

---

### 2.2 `methods/finitedifferences/operators/`

#### 2.2.1 — `methods/finitedifferences/operators/fdmblackscholesop.hpp`

**Module mapping:** M5 (SpatialOperator) — reference for the fitted replacement.

**Mathematical reference:** Math framework §3.2.4 (standard stencil that the fitted version replaces).

**What changes and why:** The class itself is NOT modified (preserving backward compatibility). However, a factory or configuration mechanism must allow users to select between `FdmBlackScholesOp` (standard centered differences) and `FdmFittedBlackScholesOp` (exponentially fitted). This selection occurs in `FdmBlackScholesSolver` (see §2.4.2). No changes to this file.

**Specific code regions to modify:** None.

---

#### 2.2.2 — `methods/finitedifferences/operators/triplebandlinearop.hpp`

**Module mapping:** M8 (TridiagSolver) — the Thomas algorithm lives here.

**Mathematical reference:** Math framework §3.6.1 (Thomas algorithm).

**What changes and why:** No modification required. The `solve_splitting` method already implements the Thomas algorithm and is used by both the standard and fitted operators. The fitted operator (`FdmFittedBlackScholesOp`) directly populates a `TripleBandLinearOp mapT_` member and delegates solves to it.

**Specific code regions to modify:** None.

---

### 2.3 `methods/finitedifferences/solvers/`

#### 2.3.1 — `methods/finitedifferences/solvers/fdmbackwardsolver.hpp`

**Module mapping:** M7 (TimeStepper) + M12 (Orchestrator) — scheme dispatch and monitoring-restart damping at the driver level.

**Mathematical reference:** Math framework §3.3.1 (damping restart after monitoring); Algorithm 1 lines 183–193 (monitoring projection); software architecture §2.M7 (damping state machine).

**What changes and why:** The existing `FdmBackwardSolver::rollback()` performs all damping steps at the start of the backward evolution, then does all main-scheme steps in a single block. For monitoring-restart damping, the rollback must be split at monitoring times: after each monitoring event triggers a barrier projection (which reintroduces a discontinuity), additional damping steps must be performed before resuming the main scheme. This requires the `rollback` method to be aware of monitoring times and to interleave damping and main-stepping segments.

Additionally, the `FdmSchemeDesc` struct gains a field for the number of monitoring-restart damping steps (distinct from the initial damping steps already in `FdmSolverDesc`).

**Specific code regions to modify:**

The `FdmSchemeDesc` struct (approximately lines 10–30 of the header) gains a new member:

```cpp
struct FdmSchemeDesc {
    /* existing members unchanged */
    FdmSchemeType type;
    Real theta, mu;
    Size monitoringDampingSteps;              // NEW — default 0

    // existing static factory methods get overloads:
    static FdmSchemeDesc CrankNicolson();     // unchanged
    static FdmSchemeDesc CrankNicolsonWithMonitoringDamping(
        Size monitoringDampingSteps = 2);     // NEW
    /* other factories unchanged */
};
```

The `FdmBackwardSolver` class declaration gains access to monitoring times:

```cpp
class FdmBackwardSolver {
  public:
    /* existing interface unchanged */
    void rollback(array_type& a,
                  Time from, Time to,
                  Size steps, Size dampingSteps);

  private:
    /* existing members */
    void rollbackSegment(array_type& a,
                         Time from, Time to,
                         Size steps, Size dampingSteps,
                         bool notifySchemeOfDiscontinuity);  // NEW helper
};
```

**New/changed method signatures (`.cpp` file):**

The `rollback` method (approximately lines 15–60 of the `.cpp`) is restructured. If `schemeDesc_.monitoringDampingSteps > 0`, it queries `condition_->stoppingTimes()` to identify monitoring times that fall within `[to, from]`. It then splits the backward evolution into segments between consecutive monitoring times. For each segment, it calls the existing damping + main-stepping logic. After each monitoring event (where the step condition applies a barrier projection), it invokes the scheme's `notifyDiscontinuity()` method (if the scheme is `CrankNicolsonScheme`) so that the next segment begins with damping sub-steps.

If `schemeDesc_.monitoringDampingSteps == 0`, the method behaves identically to the existing implementation (full backward compatibility).

**Failure modes addressed:** F3 (monitoring resets), F2 (re-injected discontinuities).

**Downstream impact:** `FdmBlackScholesSolver` and all pricing engines that construct `FdmSchemeDesc` can now pass monitoring-damping parameters.

---

#### 2.3.2 — `methods/finitedifferences/solvers/fdmbackwardsolver.cpp`

Implementation changes corresponding to §2.3.1. The major change is in the `rollback` body. The new helper `rollbackSegment` encapsulates the existing "damping then main" pattern for a single time segment.

---

#### 2.3.3 — `methods/finitedifferences/solvers/fdmblackscholessolver.hpp`

**Module mapping:** M12 (Orchestrator) — BS-specific solver that creates the operator.

**Mathematical reference:** Math framework §3.2 (spatial operator selection).

**What changes and why:** The existing `FdmBlackScholesSolver` hardcodes creation of `FdmBlackScholesOp` in its `performCalculations()` method. A boolean parameter must be added to select between the standard operator and the fitted operator.

**Specific code regions to modify:**

The constructor (approximately lines 10–20 of the header) gains a new parameter:

```cpp
class FdmBlackScholesSolver : public LazyObject {
  public:
    FdmBlackScholesSolver(
        Handle<GeneralizedBlackScholesProcess> process,
        Real strike,
        FdmSolverDesc solverDesc,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
        bool localVol = false,
        Real illegalLocalVolOverwrite = -Null<Real>(),
        Handle<FdmQuantoHelper> quantoHelper = Handle<FdmQuantoHelper>(),
        bool useFittedOperator = false);      // NEW

    /* existing methods unchanged */

  private:
    /* existing members */
    const bool useFittedOperator_;            // NEW
};
```

In the `.cpp` `performCalculations()` method (approximately lines 10–25), the operator creation is conditioned:

```cpp
ext::shared_ptr<FdmLinearOpComposite> op;
if (useFittedOperator_) {
    op = ext::make_shared<FdmFittedBlackScholesOp>(
        solverDesc_.mesher, process_.currentLink(), strike_,
        localVol_, illegalLocalVolOverwrite_, 0,
        quantoHelper_.empty()
            ? ext::shared_ptr<FdmQuantoHelper>()
            : quantoHelper_.currentLink());
} else {
    op = ext::make_shared<FdmBlackScholesOp>(
        solverDesc_.mesher, process_.currentLink(), strike_,
        localVol_, illegalLocalVolOverwrite_, 0,
        quantoHelper_.empty()
            ? ext::shared_ptr<FdmQuantoHelper>()
            : quantoHelper_.currentLink());
}
solver_ = ext::make_shared<Fdm1DimSolver>(solverDesc_, schemeDesc_, op);
```

**Failure modes addressed:** F1, F4 (via operator selection).

---

#### 2.3.4 — `methods/finitedifferences/solvers/fdmblackscholessolver.cpp`

Implementation changes corresponding to §2.3.3.

---

#### 2.3.5 — `methods/finitedifferences/solvers/fdm1dimsolver.hpp`

**Module mapping:** M10 (GreeksEngine) — Greeks extraction.

**Mathematical reference:** Math framework §3.5.1 (Delta/Gamma from log-derivatives), §3.5.2 (quadratic interpolation).

**What changes and why:** The existing `Fdm1DimSolver` already provides `interpolateAt(Real x)`, `derivativeX(Real x)`, `derivativeXX(Real x)`, and `thetaAt(Real x)`, which compute price, Delta-helper, Gamma-helper, and Theta respectively using cubic interpolation on the final solution. The Greeks conversion from log-derivatives to S-space derivatives is done externally (in the pricing engine's `fetchResults`). No structural change is needed here, but verify that the cubic interpolation uses at least quadratic-accurate derivative evaluation. The existing implementation uses `CubicInterpolation` which is adequate.

**Specific code regions to modify:** None required. This entry documents the existing M10 integration point for reference.

---

#### 2.3.6 — `methods/finitedifferences/solvers/fdmsolverdesc.hpp`

**Module mapping:** M12 (Orchestrator) — solver configuration.

**Mathematical reference:** Software architecture §1.3 (SolverConfig data structure).

**What changes and why:** The `FdmSolverDesc` struct is a plain aggregate with no modification needed for the core changes. The monitoring-restart damping configuration lives in `FdmSchemeDesc` (§2.3.1). The diagnostics level could optionally be added here as a configuration field for production use.

**Specific code regions to modify:**

Optionally add a diagnostics level field:

```cpp
struct FdmSolverDesc {
    /* existing fields unchanged */
    const ext::shared_ptr<FdmMesher> mesher;
    const FdmBoundaryConditionSet bcSet;
    const ext::shared_ptr<FdmStepConditionComposite> condition;
    const ext::shared_ptr<FdmInnerValueCalculator> calculator;
    const Time maturity;
    const Size timeSteps;
    const Size dampingSteps;
    FdmDiagnostics::Level diagLevel;          // NEW — optional
};
```

**Failure modes addressed:** All (via diagnostic configuration).

---

### 2.4 `methods/finitedifferences/meshers/`

#### 2.4.1 — `methods/finitedifferences/meshers/fdmblackscholesmesher.hpp`

**Module mapping:** M1 (GridFactory) — mesher selection.

**Mathematical reference:** Math framework §3.2.1 (grid construction), software architecture §3.4 (sinh grid).

**What changes and why:** The existing `FdmBlackScholesMesher` constructs a log-space grid using either `Uniform1dMesher` or `Concentrating1dMesher`. An additional factory method or constructor overload should allow the user to specify use of the new `FdmSinhConcentrating1dMesher` with explicit barrier alignment targets.

**Specific code regions to modify:**

Add a static factory method (approximately after line 35 of the header):

```cpp
class FdmBlackScholesMesher : public Fdm1dMesher {
  public:
    /* existing constructors unchanged */

    // NEW factory: sinh-graded mesher with barrier alignment
    static ext::shared_ptr<Fdm1dMesher> createSinhMesher(
        Size size,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Time maturity, Real strike,
        Real sinhAlpha = 3.0,
        Real xMinConstraint = Null<Real>(),
        Real xMaxConstraint = Null<Real>(),
        Real eps = 0.0001,
        Real scaleFactor = 1.5,
        const std::vector<Real>& barrierAlignTargets
            = std::vector<Real>());

    /* existing processHelper unchanged */
};
```

**Implementation notes:** The factory creates an `FdmSinhConcentrating1dMesher` with `xCenter = 0` (at-the-money in log-space), `alpha = sinhAlpha`, and `alignTargets` set to the log-barrier levels. The domain `[xMin, xMax]` is determined by the same volatility-based heuristic as the existing constructor (`eps` and `scaleFactor`).

**Failure modes addressed:** F6 (barrier alignment), F3 (grid quality for monitoring).

---

### 2.5 `methods/finitedifferences/stepconditions/`

#### 2.5.1 — `methods/finitedifferences/stepconditions/fdmstepconditioncomposite.hpp`

**Module mapping:** M6 (BoundaryHandler) + M7 (TimeStepper) — expose monitoring times for damping restart.

**Mathematical reference:** Software architecture §2.M6 (monitoring projection events trigger damping restart).

**What changes and why:** The existing `FdmStepConditionComposite` aggregates conditions and their stopping times. When the new `FdmBarrierProjectionCondition` (§1.5) is added to the composite, its monitoring times automatically become stopping times. No change to the composite class itself is needed for this integration; the stopping times are already accessible via `stoppingTimes()`. However, the `vanillaComposite` factory method should be extended with an overload or optional parameter to include barrier projection conditions.

**Specific code regions to modify:**

Add a new factory or extend `vanillaComposite` (approximately lines 30–40 of the header):

```cpp
class FdmStepConditionComposite : public StepCondition<Array> {
  public:
    /* existing interface unchanged */

    // NEW factory for discretely monitored barrier products
    static ext::shared_ptr<FdmStepConditionComposite>
    barrierMonitoredComposite(
        const DividendSchedule& schedule,
        const ext::shared_ptr<Exercise>& exercise,
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<FdmInnerValueCalculator>& calculator,
        const Date& refDate,
        const DayCounter& dayCounter,
        const std::vector<Date>& monitoringDates,
        Real lowerBarrier,
        Real upperBarrier);
};
```

**Failure modes addressed:** F3 (integration of monitoring projection into condition composite).

---

#### 2.5.2 — `methods/finitedifferences/stepconditions/fdmamericanstepcondition.hpp`

**Module mapping:** M9 (AmericanConstraint) — integration point for policy iteration.

**Mathematical reference:** Software architecture §3.2 (policy iteration), §2.M9 interface.

**What changes and why:** The existing class applies simple pointwise projection. For users who want policy iteration, a parallel class or a strategy flag is needed. The cleanest approach is to keep this class unchanged (it serves as the "simple projection" strategy) and have the pricing engine choose between using `FdmAmericanStepCondition` (projection) or using `FdmPolicyIterationLCP` (§1.7) within a modified backward solver.

**Specific code regions to modify:** None (this class remains as-is; the policy iteration alternative is a new utility class that integrates differently, at the solver level rather than the step-condition level).

---

### 2.6 `methods/finitedifferences/utilities/`

#### 2.6.1 — `methods/finitedifferences/utilities/fdminnervaluecalculator.hpp`

**Module mapping:** M3 (PayoffProcessor) — cell-averaged payoff initialization.

**Mathematical reference:** Math framework §3.4.1 (cell-averaged payoff with Gauss quadrature).

**What changes and why:** The existing `FdmCellAveragingInnerValue::avgInnerValue()` already performs cell averaging using (presumably) Simpson's rule. The math framework recommends 3-point Gauss-Legendre for higher accuracy on cells containing payoff kinks. This can be achieved by adding an optional quadrature selector or by enhancing the existing integration. However, since the existing Simpson approach is already production-grade for most payoffs and the infrastructure correctly uses `avgInnerValue` during initialization, only a minor enhancement is needed.

**Specific code regions to modify:**

Optionally add a quadrature kind enum and a constructor parameter to `FdmCellAveragingInnerValue` (approximately lines 50–65):

```cpp
class FdmCellAveragingInnerValue : public FdmInnerValueCalculator {
  public:
    enum QuadratureKind { Simpson, GaussLegendre3 };  // NEW

    FdmCellAveragingInnerValue(
        ext::shared_ptr<Payoff> payoff,
        ext::shared_ptr<FdmMesher> mesher,
        Size direction,
        ext::function<Real(Real)> gridMapping = identity<Real>(),
        QuadratureKind quad = Simpson);               // NEW parameter

    /* existing interface unchanged */

  private:
    /* existing members */
    QuadratureKind quad_;                             // NEW
};
```

In the `.cpp` implementation, the `avgInnerValueCalc` method (which performs the actual integration) is modified to use Gauss-Legendre nodes and weights when `quad_ == GaussLegendre3`: nodes `ξ = {-√(3/5), 0, +√(3/5)}`, weights `w = {5/9, 8/9, 5/9}`, evaluated over the cell interval `[x_{j-1/2}, x_{j+1/2}]`.

**Failure modes addressed:** F2 (improved payoff smoothing at kinks), F8 (cleaner Greeks initialization).

---

### 2.7 `methods/finitedifferences/` (root-level, legacy framework)

#### 2.7.1 — `methods/finitedifferences/bsmoperator.hpp`

**Module mapping:** M5 (SpatialOperator) — legacy framework standard operator.

**Mathematical reference:** Math framework §3.2.4 (standard stencil, for comparison with fitted version).

**What changes and why:** No modification. The existing `BSMOperator` is preserved. The new `FittedBSMOperator` (§1.3) is a separate class that serves as a drop-in replacement when the user selects fitted discretization in the legacy path.

**Specific code regions to modify:** None.

---

#### 2.7.2 — `methods/finitedifferences/mixedscheme.hpp`

**Module mapping:** M7 (TimeStepper) — legacy framework θ-scheme with optional damping restart.

**Mathematical reference:** Math framework §3.3.1 (damping sub-steps); Algorithm 1 lines 63–128.

**What changes and why:** The existing `MixedScheme<Operator>::step()` (approximately lines 85–115 of the header) performs a single explicit-then-implicit step. For monitoring-restart damping in the legacy path, we add a `notifyDiscontinuity()` method and a damping counter that, when active, forces `theta_` to 1.0 (pure implicit) with halved step size for a configurable number of sub-steps.

**Specific code regions to modify:**

Add new members and methods to the `MixedScheme` template class (approximately after line 60):

```cpp
template <class Operator>
class MixedScheme {
  public:
    /* existing constructors and methods */

    // NEW: monitoring-restart damping support
    void setDampingPolicy(Size halfSteps);
    void notifyDiscontinuity();
    bool isDamping() const;

  protected:
    /* existing members */
    Size dampHalfSteps_;       // NEW — total half-steps per damping episode
    Size dampRemaining_;       // NEW — remaining half-steps in current episode
};
```

The `step()` method (approximately lines 85–115) gains a preamble: if `dampRemaining_ > 0`, perform two implicit-only steps at `dt_/2` (by temporarily setting `theta_` to 1.0 and halving the operator scaling), then restore `theta_` and decrement `dampRemaining_` by 2.

**Failure modes addressed:** F2, F3 (legacy framework path).

---

#### 2.7.3 — `methods/finitedifferences/finitedifferencemodel.hpp`

**Module mapping:** M7 (TimeStepper) — legacy framework time-loop driver.

**Mathematical reference:** Math framework §3.3.1 (rollback with stopping-time handling and damping restart).

**What changes and why:** The existing `rollbackImpl` (approximately lines 55–105 of the header) handles stopping times but does not trigger damping restart after hitting one. When the evolver is a `MixedScheme` with damping capability (§2.7.2), the model should call `evolver_.notifyDiscontinuity()` after applying a step condition at a stopping time, so that subsequent steps include damping.

**Specific code regions to modify:**

In `rollbackImpl`, inside the stopping-time handling block (approximately lines 75–95), after the line that calls `condition->applyTo(a, stoppingTimes_[j])`, add:

```cpp
// Trigger damping restart in the evolver after a discontinuity event
evolver_.notifyDiscontinuity();
```

This requires that the `Evolver` type has a `notifyDiscontinuity()` method, which is satisfied by the modified `MixedScheme` (§2.7.2). For evolvers that do not have this method (e.g., `ParallelEvolver`), use SFINAE or a no-op default.

**Failure modes addressed:** F3 (monitoring-restart damping in legacy framework).

---

#### 2.7.4 — `methods/finitedifferences/tridiagonaloperator.hpp`

**Module mapping:** M8 (TridiagSolver) — legacy Thomas algorithm.

**Mathematical reference:** Math framework §3.6.1.

**What changes and why:** No modification needed. The `solveFor` method already implements a correct Thomas algorithm. The `FittedBSMOperator` (§1.3) inherits from `TridiagonalOperator` and uses `setMidRow` to set fitted coefficients; `solveFor` handles the rest.

**Specific code regions to modify:** None.

---

### 2.8 `pricingengines/vanilla/` [VERIFY: exact paths]

#### 2.8.1 — `pricingengines/vanilla/fdblackscholesvanillaengine.hpp`

**Module mapping:** M12 (Orchestrator) — user-facing configuration interface.

**Mathematical reference:** Software architecture §4.1 (product × configuration matrix), §4.2 (default parameters).

**What changes and why:** The existing constructor accepts `tGrid`, `xGrid`, `dampingSteps`, and `schemeDesc`. The improved CN engine needs additional parameters: whether to use the fitted spatial operator, the monitoring-restart damping count, the mesher type (uniform vs sinh), and optionally the diagnostics level.

**Specific code regions to modify:**

The constructor (approximately lines 15–30 of the header) is extended:

```cpp
class FdBlackScholesVanillaEngine : public ... {
  public:
    FdBlackScholesVanillaEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        Size tGrid = 100,
        Size xGrid = 100,
        Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
        bool localVol = false,
        Real illegalLocalVolOverwrite = -Null<Real>(),
        bool useFittedOperator = false,             // NEW
        Real sinhAlpha = 0.0);                      // NEW: 0 = uniform grid

    /* existing interface unchanged */

  private:
    /* existing members */
    bool useFittedOperator_;                        // NEW
    Real sinhAlpha_;                                // NEW
};
```

In the `.cpp` `calculate()` method, the mesher creation must be conditioned on `sinhAlpha_`: if `sinhAlpha_ > 0`, use `FdmBlackScholesMesher::createSinhMesher(...)` (§2.4.1); otherwise use the existing mesher construction. The operator flag `useFittedOperator_` must be forwarded to `FdmBlackScholesSolver`.

**Failure modes addressed:** F1, F4, F6 (via user-accessible configuration).

---

#### 2.8.2 — `pricingengines/vanilla/fdblackscholesvanillaengine.cpp`

Implementation changes corresponding to §2.8.1. The `calculate()` method gains conditional mesher and solver creation.

---

## 3. Integration Wiring

### 3.1 User Construction Path (New Framework)

A user prices an improved-CN European option as follows:

```cpp
// 1. Build process (unchanged)
auto process = ext::make_shared<BlackScholesMertonProcess>(
    Handle<Quote>(spot), dividendTS, riskFreeTS, volTS);

// 2. Select improved scheme with monitoring-restart damping
auto schemeDesc = FdmSchemeDesc::CrankNicolsonWithMonitoringDamping(
    /* monitoringDampingSteps = */ 2);

// 3. Construct engine with fitted operator and sinh grid
auto engine = ext::make_shared<FdBlackScholesVanillaEngine>(
    process,
    /* tGrid */ 800,
    /* xGrid */ 800,
    /* dampingSteps */ 2,       // initial Rannacher damping
    schemeDesc,
    /* localVol */ false,
    /* illegalLocalVolOverwrite */ -Null<Real>(),
    /* useFittedOperator */ true,
    /* sinhAlpha */ 3.0);

// 4. Price (unchanged)
option.setPricingEngine(engine);
Real npv = option.NPV();
```

### 3.2 Internal Call Chain

When `calculate()` fires, the chain is:

1. `FdBlackScholesVanillaEngine::calculate()` creates the mesher (sinh if `sinhAlpha_ > 0`), the cell-averaging inner value calculator, the step-condition composite (including `FdmBarrierProjectionCondition` if monitoring times exist), and the `FdmSolverDesc`.

2. It then creates `FdmBlackScholesSolver` with `useFittedOperator_ = true`.

3. `FdmBlackScholesSolver::performCalculations()` creates either `FdmFittedBlackScholesOp` or `FdmBlackScholesOp` based on the flag, then creates `Fdm1DimSolver`.

4. `Fdm1DimSolver::performCalculations()` creates `FdmBackwardSolver` with the selected `FdmSchemeDesc`.

5. `FdmBackwardSolver::rollback()` performs the modified backward evolution: for each segment between monitoring times, it executes initial damping (implicit Euler), then main stepping (CN), and after monitoring projections, calls `notifyDiscontinuity()` on the CN scheme to restart damping.

6. After rollback completes, `Fdm1DimSolver` interpolates results and the engine extracts Greeks.

### 3.3 Legacy Framework Path

For legacy engines using `FDVanillaEngine` or similar:

1. Replace `BSMOperator` with `FittedBSMOperator` in the engine's operator construction.

2. `FiniteDifferenceModel<CrankNicolson<TridiagonalOperator>>` drives time stepping. The modified `MixedScheme` (§2.7.2) supports damping restart. After stopping-time events, `FiniteDifferenceModel::rollbackImpl` calls `evolver_.notifyDiscontinuity()`.

3. No changes to `CrankNicolson<>` template itself are needed (it inherits from `MixedScheme`).

### 3.4 FdmSchemeDesc Dispatch

In `FdmBackwardSolver::rollback()`, the existing `switch(schemeDesc_.type)` for `CrankNicolsonType` now constructs the `CrankNicolsonScheme` with the `dampingHalfSteps` parameter from `schemeDesc_.monitoringDampingSteps`:

```cpp
case CrankNicolsonType: {
    CrankNicolsonScheme cnEvolver(
        schemeDesc_.theta, map_, bcSet_,
        schemeDesc_.monitoringDampingSteps,  // NEW parameter
        1e-8, ImplicitEulerScheme::BiCGstab);
    // ... rollback logic with monitoring-aware segmentation
} break;
```

No new `FdmSchemeType` enum value is required since the monitoring-damping capability is an extension of the existing `CrankNicolsonType`.

---

## 4. Build System & Header Registration

### 4.1 `methods/finitedifferences/operators/all.hpp`

Add the line:
```cpp
#include <ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp>
```

### 4.2 `methods/finitedifferences/all.hpp`

Add the line:
```cpp
#include <ql/methods/finitedifferences/fittedbsmoperator.hpp>
```

### 4.3 `methods/finitedifferences/stepconditions/all.hpp`

Add the line:
```cpp
#include <ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp>
```

### 4.4 `methods/finitedifferences/utilities/all.hpp`

Add the lines:
```cpp
#include <ql/methods/finitedifferences/utilities/fdmpolicyiteration.hpp>
#include <ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp>
```

### 4.5 `methods/finitedifferences/meshers/all.hpp`

Add the line:
```cpp
#include <ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.hpp>
```

### 4.6 Makefile.am / CMakeLists.txt

For each new `.hpp` and `.cpp` file listed in Section 1, add the file to the appropriate `Makefile.am` source list in its subdirectory. The relevant `Makefile.am` files are:

- `ql/methods/finitedifferences/operators/Makefile.am` — add `fdmfittedblackscholesop.hpp` and `fdmfittedblackscholesop.cpp`
- `ql/methods/finitedifferences/Makefile.am` — add `fittedbsmoperator.hpp` and `fittedbsmoperator.cpp`
- `ql/methods/finitedifferences/stepconditions/Makefile.am` — add `fdmbarrierprojectioncondition.hpp` and `fdmbarrierprojectioncondition.cpp`
- `ql/methods/finitedifferences/utilities/Makefile.am` — add `fdmpolicyiteration.hpp`, `fdmpolicyiteration.cpp`, `fdmdiagnostics.hpp`, `fdmdiagnostics.cpp`
- `ql/methods/finitedifferences/meshers/Makefile.am` — add `fdmsinhconcentrating1dmesher.hpp` and `fdmsinhconcentrating1dmesher.cpp`

If QuantLib v1.23 uses CMake (which it does alongside autotools), the corresponding `CMakeLists.txt` in each subdirectory must also be updated.

---

## 5. Test Scaffolding

### 5.1 Test File Location

All tests should be added to the QuantLib test suite. The primary test file for FD Black-Scholes pricing is [VERIFY: `test-suite/fdmblackscholes.cpp` or `test-suite/quantlibtestsuite.cpp`]. If a dedicated test file for the improved CN does not exist, create `test-suite/fdmimprovedcn.cpp` and register it in the test suite runner.

### 5.2 Test T1 — European Call (Baseline Greeks Convergence)

**Test function:** `testImprovedCN_EuropeanCallConvergence()`

**Parameters:** S₀=100, K=100, T=1, r=0.05, q=0, σ=0.2.

**Procedure:** Run the improved CN engine (fitted operator, Rannacher damping, cell-averaged payoff) at grid sizes (J,N) = (200,200), (400,400), (800,800). Compare price, delta, gamma against BS closed-form at each level. Verify that the Richardson order estimate for price is approximately 2.0 (within [1.5, 2.5]) and that delta and gamma errors decrease monotonically.

**Failure modes tested:** F2 (payoff kink at S=K), F8 (Greeks accuracy), F6 (domain sensitivity — repeat with different domain widths).

### 5.3 Test T2 — European Digital Put (Positivity and Damping)

**Test function:** `testImprovedCN_DigitalPutPositivity()`

**Parameters:** S₀=100, K=100, T=1, r=0.05, q=0, σ=0.2, digital payoff.

**Procedure:** Run with fitted operator and 4 damping half-steps. After solve, run `FdmDiagnostics::checkSolution()` on the final solution array. Assert `negativeCount == 0` and `oscillationScore < 0.02`.

**Failure modes tested:** F2 (discontinuous payoff), F4 (positivity).

### 5.4 Test T3 — Low Volatility Convection-Dominated (M-Matrix)

**Test function:** `testImprovedCN_LowVolNoOscillation()`

**Parameters:** S₀=100, K=100, T=0.5, r=0.05, q=0, σ=0.02.

**Procedure:** Run with both standard operator and fitted operator. For the standard operator, verify that diagnostics report M-matrix violations. For the fitted operator, verify `mMatrixSatisfied() == true` and `negativeCount == 0`. Compare prices against BS closed-form (both should converge, but fitted should be non-oscillatory at coarser grids).

**Failure modes tested:** F1 (convection dominance), F4 (positivity), F5 (numerical diffusion comparison).

### 5.5 Test T4 — Discretely Monitored Double Knock-Out

**Test function:** `testImprovedCN_DiscreteDKO()`

**Parameters:** S₀=100, K=100, T=0.5, r=0.05, q=0, σ=0.25, L=95, U=110, 5 monitoring dates.

**Procedure:** Run with monitoring-restart damping (2 half-steps after each monitoring event), fitted operator, sinh mesh with barrier alignment. Verify no negative values after any monitoring projection. Compare price against a high-resolution reference (fine PDE solve at 4× resolution).

**Failure modes tested:** F3 (monitoring resets), F2 (re-injected discontinuities), F8 (Greeks near barriers).

### 5.6 Test T5 — Down-and-Out Rebate Call

**Test function:** `testImprovedCN_DownOutRebateCall()`

**Parameters:** S₀=50, K=40, T=0.5, r=0.04, q=0, σ=0.3, B=20, Rebate=2.5.

**Procedure:** Run with barrier-aligned grid (barrier node exactly at ln(B/K)). Compare price against analytical continuous-barrier formula. Verify price within 1e-4 of analytical.

**Failure modes tested:** F6 (barrier boundary alignment), F4 (positivity near barrier).

### 5.7 Test T6 — American Put (Policy Iteration)

**Test function:** `testImprovedCN_AmericanPutPolicyIteration()`

**Parameters:** S₀=100, K=100, T=1, r=0.05, q=0, σ=0.2.

**Procedure:** Run with policy iteration LCP solver. Verify that `u[j] >= payoff[j]` for all j (complementarity). Verify convergence of price and delta against a high-resolution binomial tree reference.

**Failure modes tested:** F7 (American LCP accuracy), F4 (early-exercise constraint satisfaction).

### 5.8 Test T7 — Far-Field Boundary Sensitivity

**Test function:** `testImprovedCN_BoundarySensitivity()`

**Parameters:** S₀=300, K=100, T=5, r=0.03, q=0, σ=0.2 (deep ITM call).

**Procedure:** Run with domain widths m=6, m=8, m=10 (in standard deviations). Verify that the price at S₀ changes by less than 1e-6 between m=8 and m=10.

**Failure modes tested:** F6 (far-field truncation).

---

## 6. Dependency & Risk Summary

| Modification | Upstream Dependencies | Downstream Files Affected | Backward Compat. Risk | Failure Modes |
|---|---|---|---|---|
| `FdmFittedBlackScholesOp` (§1.1–1.2) | `FdmMesher`, `GeneralizedBlackScholesProcess`, `TripleBandLinearOp` | `FdmBlackScholesSolver` (§2.3.3), engine (§2.8.1) | **Low** — new class, no existing API change | F1, F4, F5 |
| `FittedBSMOperator` (§1.3–1.4) | `TridiagonalOperator` | Legacy FD engines | **Low** — new class, drop-in for `BSMOperator` | F1, F4 |
| `FdmBarrierProjectionCondition` (§1.5–1.6) | `FdmMesher`, `StepCondition<Array>` | `FdmStepConditionComposite` (§2.5.1), `FdmBackwardSolver` (§2.3.1) | **Low** — new class registered via composite | F3, F4 |
| `FdmPolicyIterationLCP` (§1.7–1.8) | `TripleBandLinearOp`, `TridiagonalOperator` | Engine `calculate()` methods, backward solver | **Low** — opt-in alternative to `FdmAmericanStepCondition` | F7, F8 |
| `FdmDiagnostics` (§1.9–1.10) | `Array` only | Engines (optional), test suite | **None** — purely additive utility | F1–F8 (detection) |
| `FdmSinhConcentrating1dMesher` (§1.11–1.12) | `Fdm1dMesher` | `FdmBlackScholesMesher` (§2.4.1), engines | **Low** — new mesher class | F3, F6 |
| `CrankNicolsonScheme` damping (§2.1.1–2.1.2) | `ImplicitEulerScheme` | `FdmBackwardSolver` | **Medium** — constructor signature extended; existing code using old constructor unaffected if new param defaults to 0 | F2, F3, F8 |
| `FdmBackwardSolver` monitoring restart (§2.3.1–2.3.2) | `FdmStepConditionComposite`, `CrankNicolsonScheme` | All pricing engines using FD | **Medium** — `rollback` behavior changes only when `monitoringDampingSteps > 0`; default 0 preserves old behavior | F3 |
| `FdmBlackScholesSolver` operator flag (§2.3.3–2.3.4) | `FdmFittedBlackScholesOp`, `FdmBlackScholesOp` | Pricing engines | **Low** — new constructor param with default `false` | F1, F4 |
| `FdmSchemeDesc` new field (§2.3.1) | None (struct) | `FdmBackwardSolver`, engines | **Medium** — struct gains field; existing construction must initialize it (default 0) | F3 |
| `FdmBlackScholesMesher` sinh factory (§2.4.1) | `FdmSinhConcentrating1dMesher` | Engines | **Low** — new static method, existing constructors unchanged | F6 |
| `FdmStepConditionComposite` barrier factory (§2.5.1) | `FdmBarrierProjectionCondition` | Engines for monitored barriers | **Low** — new static factory, existing factories unchanged | F3 |
| `FdmCellAveragingInnerValue` Gauss quad (§2.6.1) | None | All FD engines using cell averaging | **Low** — new optional param defaults to Simpson (existing behavior) | F2, F8 |
| `MixedScheme` damping (§2.7.2) | None | Legacy engines | **Low** — new methods with no-op defaults | F2, F3 |
| `FiniteDifferenceModel` restart (§2.7.3) | `MixedScheme` damping | Legacy engines using stopping times | **Medium** — behavior changes at stopping times if evolver supports damping; SFINAE guards needed for evolvers without the method | F3 |
| Engine constructor extension (§2.8.1–2.8.2) | `FdmBlackScholesSolver`, meshers | User code | **Low** — new params with backward-compatible defaults | F1, F4, F6 |

### Risk Mitigation Notes

The highest backward-compatibility risks are in `FdmBackwardSolver::rollback()` (§2.3.1) and `FdmSchemeDesc` (§2.3.1), because these are used by many engines. The mitigation is that all new parameters default to values that reproduce the existing behavior exactly: `monitoringDampingSteps = 0` means no monitoring-restart damping and the rollback proceeds identically to the current implementation. Existing user code that constructs `FdmSchemeDesc` via the static factories (e.g., `FdmSchemeDesc::Douglas()`) must be updated to initialize the new field, which is ensured by the struct constructor.

For `MixedScheme` (§2.7.2), the `notifyDiscontinuity()` call in `FiniteDifferenceModel` (§2.7.3) requires SFINAE or a compile-time check to avoid breaking evolvers that don't support this method. The safest approach is to add a default no-op `notifyDiscontinuity()` method at the `MixedScheme` level, making it available to all template instantiations including `CrankNicolson<>`, `ExplicitEuler<>`, and `ImplicitEuler<>`.

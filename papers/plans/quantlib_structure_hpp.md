*Thinking...*

> The user wants me to perform a detailed code-reading of QuantLib v1.23 to understand which subdirectories need modification to create a custom, improved Crank-Nicolson method for solving Black-Scholes equations numerically using finite differences.
> 
> They've provided:
> 1. A paper critique about CN schemes (critique_CN.md)
> 2. A detailed architectural analysis (quantlib_structure.md)
> 3. Source code for math utilities (math.xml)
> 4. Source code for finite difference methods (methods.xml)
> 5. Source code for pricing engines (pricingengines.xml)
> 
> 
> 6. Source code for stochastic processes (processes.xml)
> 7. Source code for term structures (termstructures.xml)
> 8. Source code for the finite differences framework headers and implementations (finitedifferences.xml)
> 
> Now I'm tracing through the execution path to map out how a finite difference Black-Scholes pricer works end-to-end. When a user instantiates an FD pricing engine like `FdBlackScholesVanillaEngine`, it kicks off a calculation chain that I need to follow through the codebase to identify all the subdirectories requiring modifications for the improved Crank-Nicolson implementation.
> 
> The engine's calculate method orchestrates several key components: it builds a mesher to discretize the underlying asset space, sets up an inner value calculator to evaluate payoffs at each node, configures step conditions to handle exercise features, establishes boundary conditions at the domain edges, constructs a solver descriptor with all the parameters, and finally instantiates the solver itself. Each of these components lives in its own subdirectory within the finite differences framework.
> 
> Once the solver is initialized, it creates the differential operator that encodes the Black-Scholes PDE, then delegates to a one-dimensional solver to handle the actual time-stepping. That solver in turn instantiates a backward solver and invokes its rollback method, which examines the scheme descriptor to determine which numerical scheme to use—for Crank-Nicolson it instantiates that scheme and wraps it in a finite difference model to execute the time steps.
> 
> The Crank-Nicolson implementation itself composes explicit and implicit Euler schemes, blending them with configurable weights during each time step. Meanwhile, the legacy framework represents Crank-Nicolson as a specialized mixed scheme with a fixed theta of 0.5, delegating the actual time-stepping logic to the mixed scheme handler.
> 
> Now I'm mapping out all the key files that need attention across the finite differences module—starting with the root-level implementations like the old Crank-Nicolson typedef, the mixed scheme that powers it, the model driver that orchestrates time stepping, and the supporting infrastructure for operators, boundary conditions, and step conditions. I'm continuing through the spatial operators with second derivatives and composite operator interfaces, then moving into the solver layer that orchestrates the backward solving process and handles Black-Scholes-specific logic, followed by the meshing components for grid generation. concentrating1dmesher.hpp/.cpp` - Concentrated mesher
>    - `fdmmeshercomposite.hpp/.cpp` - Composite mesher
>    - `fdm1dmesher.hpp` - Base 1D mesher
> 
> 6. `ql/methods/finitedifferences/stepconditions/` - Step conditions:
>    - `fdmstepconditioncomposite.hpp/.cpp` - Composite conditions
>    - `fdmamericanstepcondition.hpp/.cpp` - American exercise
>    - `fdmsnapshotcondition.hpp/.cpp` - Snapshot for theta
> 
> 7. `ql/methods/finitedifferences/utilities/` - Utilities:
>    - `fdmdirichletboundary.hpp/.cpp` - Dirichlet BC
>    - `fdminnervaluecalculator.hpp/.cpp` - Inner value calculation
>    - `fdmquantohelper.hpp/.cpp` - Quanto helper
>    - `fdmboundaryconditionset.hpp` - BC set typedef
> 
> The pricing engines and barrier option implementations will likely need modifications, along with the underlying process definitions for Black-Scholes dynamics.
> 
> Now I'm diving deeper into the actual code flow to trace exactly how the Crank-Nicolson scheme connects through the framework. The old implementation uses a template-based approach where CrankNicolson is defined as a MixedScheme with a specific theta parameter, and I need to understand how this maps to the new operator-based architecture. The Crank-Nicolson scheme splits the time step into explicit and implicit parts weighted by theta, then the explicit part applies the finite difference operator scaled by the time step and boundary conditions, while the implicit part handles the remaining portion through a linear solve.
> 
> For the implicit step, when there's only one spatial dimension it uses a direct tridiagonal solver, otherwise it falls back to iterative methods like BiCGstab or GMRES. The Black-Scholes operator computes interest and dividend rates from the term structures and constructs the finite difference matrix by combining the first and second derivative operators, which can then be applied to vectors or solved with a given time step.
> 
> Implementing the Thomas algorithm in the triple band linear operator to handle the tridiagonal solve with proper index reversal...

# Complete Code-Reading: QuantLib v1.23 Crank-Nicolson for Black-Scholes

## 1. High-Level Execution Trace

Before diving into subdirectories, let me trace the exact execution path when a user prices a European option using the FD Black-Scholes CN engine. This chain reveals every component you'd need to touch.

```
User code: option.setPricingEngine(FdBlackScholesVanillaEngine)
  │
  ▼
pricingengines/vanilla/fdblackscholesvanillaengine.cpp::calculate()
  ├── Creates FdmBlackScholesMesher           ← meshers/
  ├── Creates FdmLogInnerValue                ← utilities/
  ├── Creates FdmStepConditionComposite       ← stepconditions/
  ├── Creates FdmSolverDesc                   ← solvers/
  └── Creates FdmBlackScholesSolver           ← solvers/
        │
        ▼
      solvers/fdmblackscholessolver.cpp::performCalculations()
        ├── Creates FdmBlackScholesOp           ← operators/
        └── Creates Fdm1DimSolver               ← solvers/
              │
              ▼
            solvers/fdm1dimsolver.cpp::performCalculations()
              └── Creates FdmBackwardSolver       ← solvers/
                    │
                    ▼
                  solvers/fdmbackwardsolver.cpp::rollback()
                    ├── [damping] ImplicitEulerScheme    ← schemes/
                    └── [main]   CrankNicolsonScheme     ← schemes/
                          │        wrapped in FiniteDifferenceModel<>
                          │        from finitedifferencemodel.hpp
                          ▼
                        schemes/cranknicolsonscheme.cpp::step()
                          ├── ExplicitEulerScheme::step(a, t, 1-θ)
                          │     └── map_->apply(a)          ← operators/
                          └── ImplicitEulerScheme::step(a, t, θ)
                                └── map_->solve_splitting() ← operators/
                                      └── TripleBandLinearOp::solve_splitting()
                                            └── Thomas algorithm (tridiagonal solve)
```

There is also a parallel **old (legacy) framework** path used by `FDVanillaEngine` and related older engines:

```
cranknicolson.hpp: CrankNicolson<TridiagonalOperator> = MixedScheme<> with θ=0.5
  └── mixedscheme.hpp::step()
        ├── explicitPart_ = I - (1-θ)·dt·L  → applyTo(a)
        └── implicitPart_ = I + θ·dt·L      → solveFor(a)  [Thomas algorithm]
              wrapped in FiniteDifferenceModel<CrankNicolson<TridiagonalOperator>>
```

## 2. Subdirectory-by-Subdirectory Analysis

### 2.1 `methods/finitedifferences/schemes/` — **MUST MODIFY (Primary Target)**

This is where the new-framework CN implementation lives.

**`cranknicolsonscheme.hpp` / `cranknicolsonscheme.cpp`** — The CN scheme for the new multi-dimensional FD framework. Looking at the actual code:

```cpp
CrankNicolsonScheme::CrankNicolsonScheme(Real theta, ...) 
: theta_(theta),
  explicit_(make_shared<ExplicitEulerScheme>(map, bcSet)),
  implicit_(make_shared<ImplicitEulerScheme>(map, bcSet, relTol, solverType)) {}

void CrankNicolsonScheme::step(array_type& a, Time t) {
    if (theta_ != 1.0) explicit_->step(a, t, 1.0-theta_);  // explicit half
    if (theta_ != 0.0) implicit_->step(a, t, theta_);       // implicit half
}
```

This is the critical insertion point. Any improvement — Rannacher startup, exponential fitting, adaptive θ, Richardson extrapolation — would either modify this class or create a new scheme class alongside it.

**`expliciteulerscheme.hpp` / `expliciteulerscheme.cpp`** — Used by CN for the explicit half-step:

```cpp
void ExplicitEulerScheme::step(array_type& a, Time t, Real theta) {
    map_->setTime(std::max(0.0, t - dt_), t);
    bcSet_.applyBeforeApplying(*map_);
    a += (theta*dt_) * map_->apply(a);
    bcSet_.applyAfterApplying(a);
}
```

**`impliciteulerscheme.hpp` / `impliciteulerscheme.cpp`** — Used by CN for the implicit half-step. For 1D problems it uses a direct tridiagonal solve; for multi-dimensional problems it uses BiCGStab or GMRES:

```cpp
void ImplicitEulerScheme::step(array_type& a, Time t, Real theta) {
    map_->setTime(std::max(0.0, t-dt_), t);
    if (map_->size() == 1) {
        a = map_->solve_splitting(0, a, -theta*dt_);  // direct solve
    } else {
        // iterative BiCGStab or GMRES solve
    }
}
```

**`boundaryconditionschemehelper.hpp`** — Wraps boundary condition application for all schemes. If your improved CN changes how BCs are applied at each sub-step, this needs modification.

**Other scheme files for reference:** `douglasscheme`, `hundsdorferscheme`, `craigsneydscheme`, `modifiedcraigsneydscheme`, `methodoflinesscheme`, `trbdf2scheme` — these show the design patterns used for alternative ADI/splitting schemes and serve as templates if you want to build a new scheme from scratch.

---

### 2.2 `methods/finitedifferences/` (root-level files) — **MUST MODIFY (Old Framework)**

**`cranknicolson.hpp`** — The old-framework CN. It's a one-liner typedef:

```cpp
template <class Operator>
class CrankNicolson : public MixedScheme<Operator> {
    CrankNicolson(const operator_type& L, const bc_set& bcs)
    : MixedScheme<Operator>(L, 0.5, bcs) {}  // θ = 0.5
};
```

**`mixedscheme.hpp`** — The actual θ-scheme engine for the old framework. The `step()` method is where explicit/implicit splitting happens:

```cpp
template <class Operator>
void MixedScheme<Operator>::step(array_type& a, Time t) {
    if (theta_!=1.0) {  // explicit part
        if (L_.isTimeDependent()) { L_.setTime(t); explicitPart_ = I_-((1.0-theta_)*dt_)*L_; }
        for (i=0; i<bcs_.size(); i++) bcs_[i]->applyBeforeApplying(explicitPart_);
        a = explicitPart_.applyTo(a);
        for (i=0; i<bcs_.size(); i++) bcs_[i]->applyAfterApplying(a);
    }
    if (theta_!=0.0) {  // implicit part
        if (L_.isTimeDependent()) { L_.setTime(t-dt_); implicitPart_ = I_+(theta_*dt_)*L_; }
        for (i=0; i<bcs_.size(); i++) bcs_[i]->applyBeforeSolving(implicitPart_,a);
        implicitPart_.solveFor(a, a);
        for (i=0; i<bcs_.size(); i++) bcs_[i]->applyAfterSolving(a);
    }
}
```

**`finitedifferencemodel.hpp`** — The time-stepping loop driver. The `rollbackImpl` method handles stopping times and drives the evolver:

```cpp
void rollbackImpl(array_type& a, Time from, Time to, Size steps, ...) {
    Time dt = (from-to)/steps;
    evolver_.setStep(dt);
    for (Size i=0; i<steps; ++i, t -= dt) {
        // handle stopping times...
        evolver_.step(a, now);
        if (condition) condition->applyTo(a, next);
    }
}
```

If your improved CN requires adaptive time stepping or multi-step logic, this template needs revision.

**`tridiagonaloperator.hpp` / `tridiagonaloperator.cpp`** — The Thomas algorithm implementation for the old framework. The `solveFor()` method is the direct tridiagonal solver:

```cpp
void TridiagonalOperator::solveFor(const Array& rhs, Array& result) const {
    Real bet = diagonal_[0];
    result[0] = rhs[0]/bet;
    for (Size j=1; j<=n_-1; ++j) {
        temp_[j] = upperDiagonal_[j-1]/bet;
        bet = diagonal_[j]-lowerDiagonal_[j-1]*temp_[j];
        result[j] = (rhs[j] - lowerDiagonal_[j-1]*result[j-1])/bet;
    }
    // back substitution...
}
```

**`bsmoperator.hpp` / `bsmoperator.cpp`** — The BS spatial operator for the old framework. Constructs centered-difference stencils:

```cpp
BSMOperator::BSMOperator(Size size, Real dx, Rate r, Rate q, Volatility sigma) {
    Real sigma2 = sigma*sigma;
    Real nu = r-q-sigma2/2;
    Real pd = -(sigma2/dx-nu)/(2*dx);
    Real pu = -(sigma2/dx+nu)/(2*dx);
    Real pm = sigma2/(dx*dx)+r;
    setMidRows(pd,pm,pu);
}
```

If you want to implement exponential fitting (as recommended in the Duffy paper), you would replace these centered-difference coefficients with fitted coefficients using the `coth` formula.

**`boundarycondition.hpp` / `boundarycondition.cpp`** — Neumann and Dirichlet BCs for the old framework. The `NeumannBC` and `DirichletBC` classes modify the tridiagonal system before/after applying or solving.

**`pde.hpp` / `pdebsm.hpp`** — PDE coefficient definitions. `PdeBSM` extracts drift, diffusion, and discount from the BS process. `PdeOperator` / `GenericTimeSetter` use these to build time-dependent operators.

**`stepcondition.hpp`** — Base `StepCondition<Array>` interface. All step conditions (American exercise, etc.) derive from this.

**`operatortraits.hpp`** — Binds `TridiagonalOperator`, `Array`, and `BoundaryCondition` together into a trait set used by `MixedScheme` and `FiniteDifferenceModel`.

---

### 2.3 `methods/finitedifferences/operators/` — **MUST MODIFY**

**`fdmblackscholesop.hpp` / `fdmblackscholesop.cpp`** — The BS spatial operator for the new framework. The critical `setTime()` method assembles the PDE operator:

```cpp
void FdmBlackScholesOp::setTime(Time t1, Time t2) {
    const Rate r = rTS_->forwardRate(t1, t2, Continuous).rate();
    const Rate q = qTS_->forwardRate(t1, t2, Continuous).rate();
    // ...
    mapT_.axpyb(Array(1, r-q-0.5*v), dxMap_,
                dxxMap_.mult(0.5*Array(mesher_->layout()->size(), v)),
                Array(1, -r));
}
```

This assembles `L = (r-q-σ²/2)·∂/∂x + (σ²/2)·∂²/∂x² - r` using `dxMap_` (first derivative) and `dxxMap_` (second derivative). If you want to implement exponential fitting in the new framework, you would modify how these derivative operators are combined.

**`triplebandlinearop.hpp` / `triplebandlinearop.cpp`** — The new framework's tridiagonal operator. The `solve_splitting()` method contains the Thomas algorithm:

```cpp
Disposable<Array> TripleBandLinearOp::solve_splitting(const Array& r, Real a, Real b) const {
    // Thomson algorithm for: (b + a*diag)*x = r
    Size rim1 = reverseIndex_[0];
    Real bet = 1.0/(a*dptr[rim1]+b);
    retVal[reverseIndex_[0]] = r[rim1]*bet;
    for (Size j=1; j<=layout->size()-1; j++) {
        // forward elimination + back substitution
    }
}
```

**`firstderivativeop.hpp` / `firstderivativeop.cpp`** — Centered first derivative. At boundaries, uses upwinding (first order) or downwinding. Interior stencil:

```cpp
lower_[i] = -hp/zetam1;   // centered: -(h+)/(h-(h++h-))
diag_[i]  = (hp-hm)/zeta0;
upper_[i] = hm/zetap1;
```

**`secondderivativeop.hpp` / `secondderivativeop.cpp`** — Standard second derivative stencil. At boundaries, sets to zero (relying on boundary conditions):

```cpp
lower_[i] =  2.0/zetam1;
diag_[i]  = -2.0/zeta0;
upper_[i] =  2.0/zetap1;
```

**`fdmlinearopcomposite.hpp`** — Abstract interface that all composite operators must implement: `size()`, `setTime()`, `apply()`, `apply_mixed()`, `apply_direction()`, `solve_splitting()`, `preconditioner()`.

**`fdmlinearoplayout.hpp` / `fdmlinearoplayout.cpp`** — Grid layout management. The `neighbourhood()` method handles index arithmetic for multi-dimensional grids with reflective boundary handling.

**`ninepointlinearop.hpp` / `ninepointlinearop.cpp`** — Nine-point stencil for mixed derivatives (used in 2D problems like Heston). Not directly relevant for 1D BS but important if extending to 2D.

---

### 2.4 `methods/finitedifferences/solvers/` — **MUST MODIFY**

**`fdmbackwardsolver.hpp` / `fdmbackwardsolver.cpp`** — The scheme selector and backward evolution driver. This is where `FdmSchemeDesc` is defined and where the rollback dispatches to the appropriate scheme:

```cpp
void FdmBackwardSolver::rollback(...) {
    // Damping steps with implicit Euler
    if (dampingSteps != 0 && schemeDesc_.type != ImplicitEulerType) {
        ImplicitEulerScheme implicitEvolver(map_, bcSet_);
        FiniteDifferenceModel<ImplicitEulerScheme> dampingModel(...);
        dampingModel.rollback(rhs, from, dampingTo, dampingSteps, *condition_);
    }
    // Main stepping
    switch (schemeDesc_.type) {
      case CrankNicolsonType: {
          CrankNicolsonScheme cnEvolver(schemeDesc_.theta, map_, bcSet_);
          FiniteDifferenceModel<CrankNicolsonScheme> cnModel(...);
          cnModel.rollback(rhs, dampingTo, to, steps, *condition_);
      } break;
      // ... other scheme types
    }
}
```

The `FdmSchemeDesc` struct defines scheme configuration:

```cpp
struct FdmSchemeDesc {
    enum FdmSchemeType { HundsdorferType, DouglasType, CraigSneydType, ..., CrankNicolsonType };
    FdmSchemeType type; Real theta, mu;
    static FdmSchemeDesc CrankNicolson() { return {CrankNicolsonType, 0.5, 0.0}; }
};
```

If you add a new scheme type (e.g., `ImprovedCrankNicolsonType`), you need to add it to this enum and add a case to the switch statement.

**`fdmblackscholessolver.hpp` / `fdmblackscholessolver.cpp`** — BS-specific solver that creates `FdmBlackScholesOp` and `Fdm1DimSolver`:

```cpp
void FdmBlackScholesSolver::performCalculations() const {
    const auto op = make_shared<FdmBlackScholesOp>(
        solverDesc_.mesher, process_.currentLink(), strike_, localVol_, ...);
    solver_ = make_shared<Fdm1DimSolver>(solverDesc_, schemeDesc_, op);
}
```

**`fdm1dimsolver.hpp` / `fdm1dimsolver.cpp`** — 1D solver that performs the actual rollback and interpolates results:

```cpp
void Fdm1DimSolver::performCalculations() const {
    Array rhs(initialValues_.size());
    std::copy(initialValues_.begin(), initialValues_.end(), rhs.begin());
    FdmBackwardSolver(op_, solverDesc_.bcSet, conditions_, schemeDesc_)
        .rollback(rhs, solverDesc_.maturity, 0.0,
                  solverDesc_.timeSteps, solverDesc_.dampingSteps);
    // copy results and build cubic interpolation
}
```

**`fdmsolverdesc.hpp`** — The solver description struct bundling mesher, BCs, conditions, calculator, maturity, and step counts.

---

### 2.5 `methods/finitedifferences/meshers/` — **LIKELY MODIFY**

**`fdmblackscholesmesher.hpp` / `fdmblackscholesmesher.cpp`** — Generates the spatial grid in ln(S) space. Uses `Concentrating1dMesher` or `Uniform1dMesher` depending on whether a concentration point is specified:

```cpp
if (cPoint.first != Null<Real>() && log(cPoint.first) >= xMin && ...) {
    helper = make_shared<Concentrating1dMesher>(xMin, xMax, size,
        pair<Real,Real>(log(cPoint.first), cPoint.second));
} else {
    helper = make_shared<Uniform1dMesher>(xMin, xMax, size);
}
```

Grid quality directly affects CN accuracy. If your improvement involves adaptive mesh refinement or non-standard grids, modify here.

**`concentrating1dmesher.hpp` / `concentrating1dmesher.cpp`** — Generates grids concentrated around critical points (e.g., the strike). Uses `asinh`-based transformations.

**`fdm1dmesher.hpp`** — Base class for all 1D meshers. Stores `locations_`, `dplus_`, `dminus_`.

**`fdmmeshercomposite.hpp` / `fdmmeshercomposite.cpp`** — Assembles multi-dimensional meshes from 1D meshers.

---

### 2.6 `methods/finitedifferences/stepconditions/` — **POSSIBLY MODIFY**

**`fdmstepconditioncomposite.hpp` / `fdmstepconditioncomposite.cpp`** — Composite step condition that aggregates stopping times and conditions. The `vanillaComposite` factory creates the standard condition set:

```cpp
static shared_ptr<FdmStepConditionComposite> vanillaComposite(
    const DividendSchedule& cashFlow, const shared_ptr<Exercise>& exercise,
    const shared_ptr<FdmMesher>& mesher, ...);
```

**`fdmamericanstepcondition.hpp` / `fdmamericanstepcondition.cpp`** — American exercise: `a[i] = max(a[i], innerValue)`.

**`fdmsnapshotcondition.hpp` / `fdmsnapshotcondition.cpp`** — Captures solution values at a specific time for theta calculation.

---

### 2.7 `methods/finitedifferences/utilities/` — **LIKELY MODIFY**

**`fdmdirichletboundary.hpp` / `fdmdirichletboundary.cpp`** — Dirichlet boundary conditions for the new framework. Sets fixed values on boundary indices.

**`fdminnervaluecalculator.hpp` / `fdminnervaluecalculator.cpp`** — Payoff evaluation. `FdmLogInnerValue` evaluates `payoff(exp(x))` and `FdmCellAveragingInnerValue` computes cell-averaged payoffs using Simpson integration — this smoothing can help with non-smooth payoff issues identified in the Duffy paper.

**`fdmboundaryconditionset.hpp`** — Typedef for the boundary condition set.

**`fdmquantohelper.hpp` / `fdmquantohelper.cpp`** — Quanto adjustment helper.

---

### 2.8 `pricingengines/vanilla/` — **LIKELY MODIFY**

**`fdblackscholesvanillaengine.hpp` / `fdblackscholesvanillaengine.cpp`** — The user-facing engine. Creates the full solver chain. If your improved CN requires additional configuration parameters (e.g., Rannacher step count, fitting parameters), you'd add them here:

```cpp
FdBlackScholesVanillaEngine::FdBlackScholesVanillaEngine(
    shared_ptr<GeneralizedBlackScholesProcess> process,
    Size tGrid, Size xGrid, Size dampingSteps,
    const FdmSchemeDesc& schemeDesc, ...)
```

**`fdvanillaengine.hpp` / `fdvanillaengine.cpp`** — The old-framework FD vanilla engine. Uses `BSMOperator`, `SampledCurve`, and the old `CrankNicolson<>` evolver.

---

### 2.9 `pricingengines/barrier/` — **POSSIBLY MODIFY**

**`fdblackscholesbarrierengine.cpp`** and **`fdblackscholesrebateengine.cpp`** — These use the same FD solver chain as vanilla options. If your improved CN changes the solver interface or adds new configuration parameters, these engines need corresponding updates.

---

### 2.10 `processes/` — **POSSIBLY MODIFY**

**`blackscholesprocess.hpp` / `blackscholesprocess.cpp`** — `GeneralizedBlackScholesProcess` provides `riskFreeRate()`, `dividendYield()`, `blackVolatility()`, `localVolatility()`, and `x0()`. The FD operators read these to build PDE coefficients. If your improved CN needs higher-order derivatives of volatility or special coefficient access, the process interface may need extension.

---

### 2.11 `math/` — **POSSIBLY MODIFY**

**`sampledcurve.hpp` / `sampledcurve.cpp`** — Used by the old framework for grid representation and result interpolation. The `valueAtCenter()`, `firstDerivativeAtCenter()`, `secondDerivativeAtCenter()` methods extract option price, delta, and gamma.

**`matrixutilities/bicgstab.hpp` / `bicgstab.cpp`** — BiCGStab iterative solver used by `ImplicitEulerScheme` for multi-dimensional problems.

**`matrixutilities/gmres.hpp` / `gmres.cpp`** — GMRES iterative solver, alternative to BiCGStab.

**`matrixutilities/sparseilupreconditioner.hpp` / `sparseilupreconditioner.cpp`** — ILU preconditioner for sparse systems.

**`richardsonextrapolation.hpp` / `richardsonextrapolation.cpp`** — Already in QuantLib's math library. If your improved CN uses Richardson extrapolation for time-stepping accuracy, this utility is directly usable.

---

## 3. Dependency Map

```
                    ┌─────────────────────┐
                    │  pricingengines/     │
                    │  vanilla/ barrier/   │
                    │  asian/              │
                    └──────────┬──────────┘
                               │ creates solver chain
                               ▼
                    ┌─────────────────────┐
                    │  solvers/           │
                    │  fdmbackwardsolver  │◄── selects scheme
                    │  fdmblackscholes    │
                    │  fdm1dimsolver      │
                    └──────────┬──────────┘
                               │ uses
                    ┌──────────┴──────────┐
               ┌────┤                     ├────┐
               ▼    │                     │    ▼
    ┌──────────────┐│    ┌───────────┐    │┌──────────────┐
    │  schemes/    ││    │operators/ │    ││  meshers/    │
    │  CN scheme   ││    │fdmBS op  │    ││  BS mesher   │
    │  explicit    │◄────┤tripleband│    ││  concentr.   │
    │  implicit    ││    │1st/2nd   │    ││  composite   │
    └──────┬───────┘│    │deriv ops │    │└──────────────┘
           │        │    └───────┬──┘    │
           │        │            │       │
           ▼        │            ▼       │
    ┌──────────────┐│    ┌───────────┐   │
    │root-level FD ││    │utilities/ │   │
    │mixedscheme   ││    │boundary   │◄──┘
    │cranknicolson ││    │innervalue │
    │tridiagonal   ││    │dividend   │
    │bsmoperator   ││    └───────────┘
    │fdmodel       ││
    └──────────────┘│    ┌───────────┐
                    │    │stepconds/ │
                    ├───►│american   │
                    │    │bermudan   │
                    │    │snapshot   │
                    │    │composite  │
                    │    └───────────┘
                    │
                    │    ┌───────────┐
                    └───►│processes/ │
                         │BS process│
                         └───────────┘
```

## 4. Priority Classification for Modification

**Tier 1 — Must modify (core CN implementation):**

The `methods/finitedifferences/schemes/` directory contains the new-framework CN scheme implementation. The root-level `methods/finitedifferences/` files contain the old-framework implementation. Both `cranknicolsonscheme.cpp` (new) and `cranknicolson.hpp` + `mixedscheme.hpp` (old) are the primary targets. The `solvers/fdmbackwardsolver.cpp` is where scheme selection happens and new scheme types must be registered.

**Tier 2 — Likely modify (spatial discretization and solver orchestration):**

The `methods/finitedifferences/operators/` directory, particularly `fdmblackscholesop.cpp` and `triplebandlinearop.cpp`, controls how the spatial operator is assembled and solved. The `pricingengines/vanilla/` directory contains the engines that configure and invoke the solver chain. The `methods/finitedifferences/utilities/` directory handles boundary conditions and payoff evaluation. The `methods/finitedifferences/meshers/` directory controls grid generation.

**Tier 3 — Possibly modify (supporting infrastructure):**

The `processes/blackscholesprocess.cpp` provides market data to operators. The `math/` directory provides linear algebra and interpolation. The `methods/finitedifferences/stepconditions/` directory handles exercise features. The root-level files `finitedifferencemodel.hpp` and `tridiagonaloperator.cpp` provide the time-stepping loop and tridiagonal solver respectively.

**Tier 4 — Unlikely to modify:**

The directories `cashflows/`, `currencies/`, `indexes/`, `instruments/`, `legacy/`, `models/`, `patterns/`, `quotes/`, `time/`, and `utilities/` are not directly involved, though `instruments/` defines option types and `termstructures/` provides yield curves and volatility surfaces consumed by the operators.

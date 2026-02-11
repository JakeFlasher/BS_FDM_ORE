*Thinking...*

> The user wants to understand the QuantLib v1.23 codebase structure and identify which subdirectories need to be modified to create a custom, improved Crank-Nicolson method for solving Black-Scholes equations numerically using finite differences.
> 
> Let me analyze the code structure systematically, tracing the relevant components from the top level down through the directory tree.
> 
> The Crank-Nicolson method is a finite difference scheme used to solve PDEs (like the Black-Scholes PDE). In QuantLib, this involves several layers:
> 
> 1. **The FD scheme itself** - The Crank-Nicolson time-stepping scheme
> 
> 
> 2. **The FD operators** - Differential operators that discretize the spatial derivatives
> 3. **The FD meshers** - Grid/mesh generation for the spatial domain
> 4. **The FD solvers** - Orchestrate the solution process
> 5. **The boundary conditions** - Handle boundary conditions for the PDE
> 6. **The step conditions** - Handle early exercise, etc.
> 7. **The pricing engines** - Connect the FD machinery to specific instruments
> 8. **The processes** - Define the underlying stochastic process (Black-Scholes)
> 9. **The PDE definitions** - Define the specific PDE coefficients
> 
> Now I'm tracing through the actual code structure to see where each of these components lives and how they interconnect.
> 
> QuantLib actually has two separate finite difference frameworks running in parallel. The older one at the top level uses classic schemes like Crank-Nicolson (which is really just a mixed scheme with theta set to 0.5) and tridiagonal operators for 1D problems. The newer framework is organized into subdirectories with separate modules for schemes, operators, meshers, solvers, and step conditions, giving it a much more modular architecture. nitedifferencemodel.hpp` - the model that drives the time stepping
> - `bsmoperator.hpp` / `pde.hpp` / `pdebsm.hpp` - BS-specific PDE operators
> - `onefactoroperator.hpp` - one-factor PDE operator
> 
> The new framework uses:
> - `schemes/cranknicolsonscheme.hpp/.cpp` - the CN scheme for the new framework
> - `operators/` - various FD operators
> - `meshers/` - mesh generation
> - `solvers/` - FD solvers
> - `stepconditions/` - step conditions
> - `utilities/` - boundary conditions, inner value calculators, etc.
> 
> ### Directories to modify:
> 
> 1. **`methods/finitedifferences/`** - This is the PRIMARY directory
>    - `schemes/` - Contains the new CN scheme implementation
>    - `operators/` - FD operators for spatial discretization
>    - `meshers/` - Mesh and grid generation
>    - `solvers/` - FD solvers orchestrating the overall process
>    - `stepconditions/` - Step conditions for early exercise and other constraints
>    - `utilities/` - Boundary conditions and related helpers
> 
> 2. **`methods/finitedifferences/` top-level files** - The old framework files like `cranknicolson.hpp` and `mixedscheme.hpp` (which implements the theta-scheme that CN is based on), along with `finitedifferencemodel.hpp`, are the ones I need to focus on updating.
> 
> 3. **`pricingengines/`** - The engines that use FD methods:
>    - `vanilla/` - Contains vanilla option engines
>    - `barrier/` - Contains barrier option engines
>    - `asian/` - Contains Asian option engines (partially listed)
> 
> 4. **`processes/`** - The Black-Scholes process definition used by FD engines
> 
> 5. **`math/`** - Mathematical utilities including sampled curves, arrays, matrices, and matrix solvers for tridiagonal systems
> 
> 6. **`experimental/finitedifferences/`** - Extended FD operators and engines
> 
> 7. **Root `ql/` files** - Discretized asset definitions used by FD methods
> 
> Now I need to trace through the actual framework more carefully. The user wants to build a custom Crank-Nicolson method for solving Black-Scholes equations, so I should map out how the existing implementation chains together. The `CrankNicolson` class is built on `MixedScheme` with a θ parameter of 0.5, which itself relies on `OperatorTraits` to define the finite difference operators. Let me follow this dependency chain to understand what components I need to work with. Looking at the new framework structure, I can see the Crank-Nicolson scheme implementation for multi-dimensional finite difference, along with the operator hierarchy for handling derivatives and the backward solver that ties everything together. The framework also includes specialized solvers for Black-Scholes problems and 1D cases, mesh generation utilities, and boundary condition handling. For pricing engines, the old framework had vanilla option engines using the finite difference approach, with more specialized engines building on top of that foundation. I'm mapping out the directory structure for the finite differences methods—the old framework has everything at the root level with files like the Crank-Nicolson scheme and PDE operators, while the new framework organizes these into dedicated subdirectories for schemes, operators, solvers, mesh generation, step conditions, and utilities.
> 
> Now I'm seeing how the pricing engines layer on top of this, with separate directories for vanilla and barrier options that leverage the finite differences solvers, plus the supporting infrastructure like the Black-Scholes process, math utilities for curves and arrays, and root-level components for discretized assets and time grids.

# Code Architecture for a Custom Crank-Nicolson Scheme in QuantLib v1.23

This is a deep dive into the QuantLib v1.23 codebase to trace every component involved in the Crank-Nicolson (CN) finite difference method for Black-Scholes. I'll walk through the general architecture, then identify each relevant subdirectory and explain why it matters.

---

## 1. Overview: Two Parallel FD Frameworks

QuantLib v1.23 contains **two coexisting finite difference frameworks** within `methods/finitedifferences/`, and this is the single most important thing to understand before making any changes.

**The old (legacy) framework** lives in the top-level files of `methods/finitedifferences/` — files like `cranknicolson.hpp`, `mixedscheme.hpp`, `finitedifferencemodel.hpp`, `tridiagonaloperator.hpp`, `bsmoperator.hpp`, `pde.hpp`, and `pdebsm.hpp`. This framework is 1D-oriented, uses `TridiagonalOperator` as its workhorse, and defines `CrankNicolson<>` as a simple typedef of `MixedScheme<>` with θ = 0.5. It drives older pricing engines like `FDVanillaEngine` found in `pricingengines/vanilla/`.

**The new (multidimensional) framework** lives in the subdirectories of `methods/finitedifferences/` — namely `schemes/`, `operators/`, `meshers/`, `solvers/`, `stepconditions/`, and `utilities/`. Here, the Crank-Nicolson scheme is a standalone class in `schemes/cranknicolsonscheme.cpp/.hpp`. This framework supports multi-dimensional PDEs and drives the newer `Fd*` pricing engines like `FdBlackScholesVanillaEngine`.

Any "improved" CN implementation will need to decide which framework to target — or both. The following analysis covers both paths.

---

## 2. Subdirectories That Need Revision

### 2.1 `methods/finitedifferences/` — **The Core FD Directory (PRIMARY)**

This is the heart of the matter. It contains everything related to finite difference discretization and time-stepping.

**Top-level files (old framework):**

The old CN scheme is defined through a chain of abstractions. `cranknicolson.hpp` simply typedefs `MixedScheme` with θ = 0.5. The actual time-stepping logic — the splitting of implicit and explicit parts, the Thomas algorithm solve — lives in `mixedscheme.hpp`, which is a template parameterized on operator traits. The traits are defined in `operatortraits.hpp`, which binds together `TridiagonalOperator` (from `tridiagonaloperator.hpp/.cpp`), the `Array` value type, and `SampledCurve` as the boundary condition type. The model driver in `finitedifferencemodel.hpp` orchestrates the time evolution by calling the scheme's `step()` method at each time level. The spatial operator for Black-Scholes is built in `bsmoperator.hpp/.cpp` using finite difference stencils (`DPlus`, `DMinus`, `DPlusDMinus`, `DZero` from their respective headers) and PDE coefficients from `pdebsm.hpp` (which inherits the interface in `pde.hpp`). Boundary conditions are handled via `boundarycondition.hpp/.cpp`, and step conditions for features like early exercise are in `stepcondition.hpp` and `americancondition.hpp`.

If you want to improve the CN scheme in the old framework — for example by adding Richardson extrapolation, adaptive time-stepping, Rannacher startup steps, or improved boundary handling — you would modify `mixedscheme.hpp` and/or `cranknicolson.hpp`, potentially `finitedifferencemodel.hpp`, and possibly `tridiagonaloperator.hpp/.cpp` if your improvement changes how the linear system is assembled or solved.

**`schemes/` subdirectory (new framework):**

This contains `cranknicolsonscheme.cpp` and `cranknicolsonscheme.hpp`, which implement CN for the new multi-dimensional FD framework. The scheme here works with `FdmLinearOpComposite` rather than `TridiagonalOperator`, and it supports operator splitting via a θ parameter and an ADI-like structure. Also relevant in this directory are the adjacent schemes: `expliciteulerscheme`, `impliciteulerscheme`, `douglasscheme`, `hundsdorferscheme`, `craigsneydscheme`, `modifiedcraigsneydscheme`, `methodoflinesscheme`, and `trbdf2scheme`. Understanding these alternatives is useful if your improvement borrows ideas from higher-order or split-step methods. The `boundaryconditionschemehelper.hpp` file provides the boundary condition integration point for all schemes.

**`operators/` subdirectory:**

This contains the spatial discretization operators for the new framework. The key file for BS is `fdmblackscholesop.cpp/.hpp`, which assembles the BS differential operator using `FirstDerivativeOp`, `SecondDerivativeOp`, and `TripleBandLinearOp` (all in this directory). There's also `fdm2dblackscholesop` for 2D cases, and operators for other models (Heston, Bates, Hull-White, etc.). If your improved CN requires a different spatial discretization — perhaps higher-order stencils or a non-uniform coefficient treatment — these operators need modification. The `fdmlinearop.hpp` and `fdmlinearopcomposite.hpp` define the abstract interfaces that all operators must satisfy.

**`meshers/` subdirectory:**

Grid generation lives here. `fdmblackscholesmesher.cpp/.hpp` generates the spatial mesh for BS problems, typically using a concentration of points around the strike. `concentrating1dmesher.cpp/.hpp` provides the general concentrated 1D mesh. `fdmmeshercomposite.cpp/.hpp` assembles multi-dimensional meshes. If your CN improvement involves adaptive mesh refinement or a non-standard grid layout, this is where changes go.

**`solvers/` subdirectory:**

The solver classes orchestrate the full solution process by connecting scheme, operator, mesher, boundary conditions, and step conditions. `fdmbackwardsolver.cpp/.hpp` is the generic backward-in-time solver that calls the scheme's `step()` method. `fdmblackscholessolver.cpp/.hpp` is the BS-specific solver. `fdm1dimsolver.cpp/.hpp` handles 1D problems. `fdmndimsolver.hpp` is the N-dimensional template solver. `fdmsolverdesc.hpp` defines the solver description struct that bundles all configuration. These are the files where the overall time-stepping loop lives, so if your improvement involves, for instance, adaptive time-step selection or a modified iteration strategy, the solver layer needs changes.

**`stepconditions/` subdirectory:**

Step conditions implement features that interrupt the backward evolution, such as `fdmamericanstepcondition` (early exercise), `fdmbermudanstepcondition`, `fdmsnapshotcondition`, and `fdmstepconditioncomposite`. If your improved CN interacts differently with early-exercise boundaries or requires a modified projection step, these files are relevant.

**`utilities/` subdirectory:**

This contains boundary condition implementations (`fdmdirichletboundary`, `fdmdiscountdirichletboundary`, `fdmtimedepdirichletboundary`), inner value calculators (`fdminnervaluecalculator`), dividend handling (`fdmdividendhandler`), and quanto adjustments (`fdmquantohelper`). Boundary condition treatment is a common area for CN improvements — for example, implementing Crandall's method or higher-order boundary extrapolation.

---

### 2.2 `pricingengines/` — **Where FD Methods Meet Instruments**

The pricing engines are the user-facing entry points that connect the FD machinery to specific financial instruments.

**`pricingengines/vanilla/`** is the most important subdirectory here. It contains both old and new framework engines for vanilla options. On the old side, `fdvanillaengine.cpp/.hpp` is the base class, with `fdstepconditionengine.hpp`, `fdmultiperiodengine.hpp`, `fddividendengine.hpp`, and `fddividendshoutengine.hpp` as derived classes. On the new side, `fdblackscholesvanillaengine.cpp/.hpp` is the primary engine. There are also FD engines for Heston (`fdhestonvanillaengine`), Bates (`fdbatesvanillaengine`), SABR (`fdsabrvanillaengine`), and CIR (`fdcirvanillaengine`). The old-style `fdconditions.hpp` provides condition setup for the legacy framework.

**`pricingengines/barrier/`** contains `fdblackscholesbarrierengine.cpp/.hpp` and `fdblackscholesrebateengine.cpp/.hpp` for barrier options using the new FD framework, plus Heston-based barrier engines.

**`pricingengines/asian/`** contains `fdblackscholesasianengine.cpp/.hpp`.

**`pricingengines/basket/`** contains `fd2dblackscholesvanillaengine.cpp/.hpp` for 2D BS.

If you're creating a new CN scheme, you'll likely need to either modify these engines to accept your custom scheme configuration, or create new engine classes that use your improved scheme. At minimum, the engines that currently instantiate `CrankNicolsonScheme` (or the old `CrankNicolson<>`) will need to be updated to use your replacement.

---

### 2.3 `processes/` — **Stochastic Process Definitions**

The file `blackscholesprocess.cpp/.hpp` defines `GeneralizedBlackScholesProcess` and its variants (`BlackScholesProcess`, `BlackScholesMertonProcess`, `BlackProcess`, `GarmanKohlagenProcess`). These classes provide the drift, diffusion, and local volatility that the FD operators use to build PDE coefficients. If your improved CN requires additional information from the process — such as higher-order derivatives of volatility for Richardson extrapolation, or time-dependent coefficients handled in a specific way — you may need to extend the process interface.

The `eulerdiscretization.cpp/.hpp` and `endeulerdiscretization.cpp/.hpp` files provide time discretization for the SDE, which is related but distinct from the PDE discretization; they're less likely to need changes unless your improvement bridges Monte Carlo and FD approaches.

---

### 2.4 `math/` — **Mathematical Infrastructure**

Several components in `math/` underpin the FD framework. `array.hpp` defines the `Array` class used for the solution vector at each time step. `matrix.hpp/.cpp` defines the `Matrix` class. `sampledcurve.hpp/.cpp` provides the `SampledCurve` class used by the old framework for grid representation and interpolation of results. `comparison.hpp` provides `close_enough()` used throughout for floating-point comparison.

The `math/matrixutilities/` subdirectory contains sparse matrix and iterative solver infrastructure (`bicgstab.cpp/.hpp`, `gmres.cpp/.hpp`, `sparseilupreconditioner.cpp/.hpp`, `sparsematrix.hpp`) that the new framework's `TripleBandLinearOp` can use. If your improved CN requires a different linear solver — for example, a preconditioned Krylov method instead of the direct Thomas algorithm — this is where those solvers live.

The `math/interpolations/` subdirectory provides interpolation routines used when extracting results from the FD grid, including `cubicinterpolation.hpp` and `linearinterpolation.hpp`.

---

### 2.5 `experimental/finitedifferences/` — **Extended FD Operators and Engines**

This subdirectory extends the new FD framework with additional operators and engines for models like extended Ornstein-Uhlenbeck, ZABR, Kluge, and Heston SLV. Key files include `fdmblackscholesfwdop.cpp/.hpp` (the forward BS operator, useful for Dupire-style local vol), `fdmdupire1dop.cpp/.hpp`, and `modtriplebandlinearop.hpp` (a modified triple-band operator). If your improved CN needs a modified operator structure, `modtriplebandlinearop.hpp` is a precedent for how to do that within the existing architecture.

---

### 2.6 Root `ql/` Level — **Supporting Infrastructure**

Several root-level files provide infrastructure that the FD framework depends on. `discretizedasset.hpp/.cpp` defines `DiscretizedAsset` and `DiscretizedOption`, which provide the old framework's interface between lattice/FD methods and instruments — the `rollback()`, `partialRollback()`, `preAdjustValues()`, and `postAdjustValues()` methods. `numericalmethod.hpp` defines the `Lattice` base class with its `initialize()`, `rollback()`, `partialRollback()`, and `presentValue()` interface. `timegrid.hpp/.cpp` defines the `TimeGrid` class that manages the time discretization for both lattice and FD methods. These are unlikely to need heavy modification, but you should understand them because they define the contracts your improved scheme must respect.

---

### 2.7 `termstructures/` — **Market Data Input**

The FD operators need yield curves and volatility surfaces. The most relevant subdirectories are `termstructures/yield/` (yield term structures, especially `flatforward.hpp` and `zeroyieldstructure.hpp`) and `termstructures/volatility/equityfx/` (Black vol surfaces, local vol surfaces). You're unlikely to modify these, but the FD operators read from them extensively, so understanding the `BlackVolTermStructure` → `LocalVolSurface` chain matters for understanding how PDE coefficients are computed.

---

## 3. Summary of Subdirectories by Priority

**Must modify:**

The directory `methods/finitedifferences/` is the primary target. Within it, `schemes/` contains the CN scheme implementation for the new framework, while the root-level `cranknicolson.hpp` and `mixedscheme.hpp` contain the old framework implementation. The `solvers/` subdirectory controls the time-stepping loop. The `operators/` subdirectory (particularly `fdmblackscholesop`) defines spatial discretization for BS.

**Likely modify:**

The `pricingengines/` directory — especially `pricingengines/vanilla/` — contains the engines that instantiate and configure the FD solver chain. You'll need to update or create engines that use your improved scheme. The `methods/finitedifferences/utilities/` subdirectory handles boundary conditions, and `methods/finitedifferences/meshers/` handles grid generation, both of which are common targets for CN improvements.

**Possibly modify:**

The `processes/` directory (specifically `blackscholesprocess`) if you need additional process information. The `math/` directory if you need new linear solvers or interpolation methods. The `experimental/finitedifferences/` directory if you're extending to non-standard operators. The root-level `discretizedasset.hpp`, `timegrid.hpp`, and `numericalmethod.hpp` if your scheme changes the time-stepping interface.

**Unlikely to modify:**

The directories `cashflows/`, `currencies/`, `indexes/`, `instruments/`, `legacy/`, `models/`, `patterns/`, `quotes/`, `time/`, and `utilities/` are not directly involved in the FD scheme itself, though `instruments/` defines the option classes your engine will price, and `patterns/` provides the Observer/Observable framework that propagates market data changes.

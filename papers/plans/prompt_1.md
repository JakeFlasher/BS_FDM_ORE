*Thinking...*

> The user is asking about two things:
> 
> 1. The `experimental/finitedifferences` directory in QuantLib - what's in there and does it need modification?
> 2. What additional files/directories are needed for a custom test suite and benchmark?
> 
> Let me think about what's in QuantLib v1.23's experimental/finitedifferences and what a proper test suite would require.
> 
> In QuantLib, the `experimental/` directory contains code that's not yet promoted to the main library. For finite differences, there might be experimental schemes, operators, or solvers there. Let me also think about the test infrastructure.
> 
> For QuantLib's test suite, the tests are typically in `test-suite/` at the root level. 
> 
> Looking at the existing test files, I can see there are already tests for FDM linear operators, Heston model pricing, barrier options, American options, and European options - many of which use finite difference methods. There's also a speed benchmarking file in the test suite. For a custom benchmark suite, I'd need to set up additional infrastructure beyond what's already there.
> 
> Now I'm thinking through what's in the experimental finite differences directory and how it fits into the broader QuantLib structure. The experimental code likely contains newer or less-stable FD schemes, boundary condition implementations, and mesh generation utilities that are still being developed. I need to understand whether these experimental components would need modifications for the multi-round prompt approach, and how the test suite and benchmarks should be organized to cover both stable and experimental code paths. - Local vol forward PDE methods
> - Risk neutral density calculations
> - Some experimental schemes
> 
> But I don't have the exact file listing from the user's attachments. Let me focus on what I know and what the user needs.
> 
> For a proper test suite and benchmark, the user would need:
> 
> **Test Suite Files:**
> - `test-suite/fdmlinearop.cpp` (or similar) - existing tests for FDM operators
> - A new test file for the fitted scheme specifically
> - Tests comparing fitted vs unfitted results
> - Convergence rate tests
> - Oscillation detection tests (the key metric from Duffy's paper)
> - Greek accuracy tests (delta, gamma near strike)
> - Small volatility regime tests
> - Barrier option tests
> 
> **Benchmark Files:**
> - Performance comparison: fitted vs standard CN
> - Grid refinement studies
> - Time-stepping convergence studies
> 
> **Other infrastructure:**
> I'd also need to update the build configuration files to register the new test and benchmark files, and potentially add preprocessor flags or linker settings if the fitted scheme requires special compilation directives.





Good catches. Here's the complete picture of what's missing from the original plan.

---

## The `ql/experimental/finitedifferences/` Directory

In QuantLib v1.23, `ql/experimental/finitedifferences/` contains specialized engines and utilities for energy/commodity-style derivatives—things like the Kluge-ExtOU spread engine, jump model inner value calculators, and related infrastructure.[[7]](https://quantnet.com/threads/problems-with-installing-boost-1-55-0-and-quantlib-1-3-in-microsoft-visual-studio-2012.15461/) Other experimental code, such as the ZABR volatility model in `ql/experimental/volatility/zabr.cpp`, directly consumes components from both `ql/experimental/finitedifferences/` (like `glued1dmesher.hpp`) and the main `ql/methods/finitedifferences/` framework.[[10]](https://github.com/lballabio/QuantLib/blob/master/ql/methods/finitedifferences/cranknicolson.hpp)

The key files in `experimental/finitedifferences/` for v1.23 include forward-PDE operators (`fdmblackscholesfwdop`, `fdmsquarerootfwdop`, `fdmhestonfwdop`), specialized meshers (`glued1dmesher`), energy derivative engines (`fdextoujumpvanillaengine`, `fdklugeextouspreadengine`), dynamic programming VPP engines, and specialized inner value calculators. The Heston FD test file, for instance, includes `#include <ql/experimental/finitedifferences/fdmblackscholesfwdop.hpp>` and related experimental forward operators.[[2]](https://github.com/quantlibjs/test-suite)

**Impact on our implementation:** The experimental forward operators (`fdmblackscholesfwdop` etc.) are the *Fokker-Planck forward equation* versions of the BS operator—probability density evolvers rather than option price evolvers. Since the exponential fitting addresses the *backward* pricing PDE's convection-diffusion instability, the experimental forward operators are a **low-priority but not-zero** concern. If someone uses the fitted scheme for calibration via forward PDEs (local vol RND calculation), they'd eventually want fitted forward operators too. For the initial implementation, though, you can skip this directory.

---

## The `test-suite/` Directory — What You Need

The existing test suite in `test-suite/fdmlinearop.cpp` already contains tests for the FDM framework including layout, meshers, derivative operators, the tridiagonal solve, Heston barrier and American tests, BiCGStab, GMRES, and notably a test for "Crank-Nicolson with initial implicit damping steps for a digital option."[[2]](https://rkapl123.github.io/QLAnnotatedSource/d4/df6/test.html)

The existing `testCrankNicolsonWithDamping` test specifically validates CN with Rannacher-style damping on a digital option—exactly the failure mode Duffy's paper identifies.[[3]](https://deepwiki.com/lballabio/QuantLib/2.2-option-models)

The relevant existing test files are:

| File | What it tests | Relevance |
|---|---|---|
| `test-suite/fdmlinearop.cpp` | FDM operators, meshers, schemes, iterative solvers | **Primary**—add fitted operator tests here |
| `test-suite/europeanoption.cpp` | European vanilla options across all engines | **Primary**—add fitted FD BS engine tests |
| `test-suite/americanoption.cpp` | American options including FD engines | **High**—convergence/oscillation tests |
| `test-suite/barrieroption.cpp` | Barrier options including FD BS and Heston | **High**—barrier options are where CN fails most |
| `test-suite/dividendoption.cpp` | Dividend options with FD engines | **Medium**—discrete dividends + FD |
| `test-suite/fdheston.cpp` | FD Heston model tests | **Low**—unless extending fitted schemes to Heston |
| `test-suite/digitaloption.cpp` | Digital/binary options | **High**—non-smooth payoffs are the core issue |
| `test-suite/operators.cpp` | Legacy tridiagonal/BSM operators | **Medium**—legacy fitted BSM operator tests |

Beyond the source files you also need:

| File | Purpose |
|---|---|
| `test-suite/toplevelfixture.hpp` | Base fixture class all tests inherit from |
| `test-suite/utilities.hpp` / `utilities.cpp` | Helper macros and functions (`QUANTLIB_TEST_CASE`, `Flag`, etc.) |
| `test-suite/CMakeLists.txt` or `Makefile.am` | Build registration for new test files |
| `test-suite/quantlibbenchmark.cpp` | The benchmark runner (QuantLib's built-in speed test harness) |

---

## Complete Expanded Round Plan

The original 5 rounds should become **7 rounds** to cover experimental FD, tests, and benchmarks:

```
Round 1: operators/              ← fitted spatial operator (mathematical core)
Round 2: schemes/ + root FD      ← time-stepping logic  
Round 3: solvers/                ← solver orchestration
Round 4: pricingengines/         ← user-facing API
Round 5: utilities/ + meshers/   ← supporting infrastructure + legacy BSM
Round 6: test-suite/             ← correctness tests ← NEW
Round 7: benchmarks + experimental/ ← performance + forward PDE ← NEW
```

---

## ROUND 6 PROMPT — Test Suite

```
You are writing a comprehensive test suite for the exponentially fitted 
Crank-Nicolson implementation in QuantLib v1.23, following the exact failure 
modes identified in Duffy's "A Critique of the Crank Nicolson Scheme" (2004).

This is ROUND 6 of 7. Rounds 1-5 produced:
  • FdmFittedBlackScholesOp (fitted spatial operator)
  • New FdmSchemeDesc types (FittedImplicitEuler, RannacherCN)
  • RichardsonExtrapolationScheme
  • Modified FdmBlackScholesSolver with useExponentialFitting
  • Modified FdBlackScholesVanillaEngine exposing fitting option
  • BSMFittedOperator for legacy framework

═══════════════════════════════════════════════════════════════════════════
FILES TO PRODUCE
═══════════════════════════════════════════════════════════════════════════

1. NEW: test-suite/fittedfinitedifferences.cpp
   Complete Boost.Test file with these test cases:

   A. OPERATOR UNIT TESTS
   • testFittingFactorLimits
     - Verify ρ → σ as μ → 0 (pure diffusion limit)
     - Verify ρ → ±μh/2 as σ → 0 (upwind limit)
     - Verify ρ = (μh/2)·coth(μh/(2σ)) for normal parameters
     - Verify monotonicity: a_{j,j-1} > 0, a_{j,j} < 0, a_{j,j+1} > 0
   
   • testFittedOperatorApply
     - Apply FdmFittedBlackScholesOp to a known smooth function
     - Compare against analytic PDE residual
     - Verify second-order spatial accuracy via grid refinement
   
   • testFittedVsUnfittedOperator
     - For moderate vol (σ=20%), both operators should give similar results
     - For small vol (σ=1%), fitted operator should be dramatically better

   B. SCHEME-LEVEL CONVERGENCE TESTS
   • testFittedImplicitEulerConvergence
     - European call, K=100, S₀=100, r=5%, q=0%, σ=20%, T=1
     - Compare against analytic Black-Scholes at grid refinements
       (N=25,50,100,200 time steps; J=50,100,200,400 spatial points)
     - Verify O(h+k) convergence rate
     - Verify error constant M is independent of σ
   
   • testRannacherDampingForDigitalOption
     - Cash-or-nothing call, K=100, S₀=100
     - Compare: CN alone vs Rannacher(2 steps)+CN vs FittedImplicit
     - CN alone should show oscillations in price near strike
     - Rannacher+CN and FittedImplicit should be smooth

   • testRichardsonExtrapolation
     - Run fitted implicit with N steps → get price P₁
     - Run fitted implicit with 2N steps → get price P₂
     - Verify 2P₂ - P₁ converges at O(h+k²) rate

   C. OSCILLATION DETECTION TESTS (paper's key contribution)
   • testNoSpuriousOscillationsAtStrike
     - European call, scan price grid from S=80 to S=120
     - With standard CN: detect oscillation by checking sign changes
       in second difference of computed prices
     - With fitted scheme: verify monotone price profile
   
   • testNoSpuriousOscillationsSmallVol
     - σ=0.5% (convection-dominated regime)  
     - Standard CN should oscillate; fitted scheme should not
     - Use σ(t)=σ₀·exp(-α(T-t)) for time-dependent vol case
   
   • testDeltaGammaSmoothness
     - Compute Δ = ∂V/∂S and Γ = ∂²V/∂²S via FD on the price grid
     - Verify no oscillations in delta profile near strike
     - Compare delta accuracy against analytic BS delta

   D. FULL ENGINE-LEVEL TESTS
   • testFittedEngineVsAnalytic
     - European call/put at various moneynesses
     - Tolerance: 0.01% of analytic price for 200×200 grid
   
   • testFittedBarrierOption
     - Down-and-out call with barrier near strike
     - Compare against AnalyticBarrierEngine
   
   • testFittedAmericanOption
     - American put, compare against Barone-Adesi-Whaley
     - Verify early exercise boundary is captured correctly

   • testBackwardCompatibility
     - Verify that all EXISTING tests still pass when using
       standard (unfitted) operators with default parameters
     - Every engine with useExponentialFitting=false must give
       bit-identical results to the unmodified QuantLib

2. NEW: test-suite/fittedfinitedifferences.hpp
   • Header with test suite registration

3. MODIFIED: test-suite/CMakeLists.txt (or Makefile.am)
   • Add fittedfinitedifferences.cpp to the build

4. MODIFIED: test-suite/fdmlinearop.cpp
   • Add to existing suite():
     testFittedBlackScholesOperatorConsistency
     testCrankNicolsonDampingWithFitting (extend existing damping test)

═══════════════════════════════════════════════════════════════════════════
TEST DATA (reference values from analytic formulas)
═══════════════════════════════════════════════════════════════════════════

Use BlackCalculator/AnalyticEuropeanEngine as the reference for European 
options. Use AnalyticBarrierEngine for barriers. For American options, use
existing BaroneAdesiWhaley as an approximate reference with wider tolerance.

Standard test parameters (from QuantLib's existing test conventions):
  Spot: 100.0
  Strikes: {80, 90, 100, 110, 120}
  Rates: {0.01, 0.05, 0.10}
  Vols: {0.005, 0.05, 0.20, 0.50}  ← 0.005 is the critical small-vol case
  Maturities: {0.1, 0.5, 1.0, 2.0}
  Dividend yields: {0.0, 0.03}

═══════════════════════════════════════════════════════════════════════════
BOOST.TEST CONVENTIONS (match QuantLib style)
═══════════════════════════════════════════════════════════════════════════

#include "toplevelfixture.hpp"
#include "utilities.hpp"
using namespace QuantLib;
using namespace boost::unit_test_framework;

BOOST_FIXTURE_TEST_SUITE(QuantLibTests, TopLevelFixture)
BOOST_AUTO_TEST_SUITE(FittedFiniteDifferencesTests)
  BOOST_AUTO_TEST_CASE(testFittingFactorLimits) { ... }
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
```

---

## ROUND 7 PROMPT — Benchmark + Experimental Forward PDE

```
You are completing the implementation with performance benchmarks and 
optional experimental forward-PDE support.

This is ROUND 7 of 7 (final round).

═══════════════════════════════════════════════════════════════════════════
FILES TO PRODUCE
═══════════════════════════════════════════════════════════════════════════

1. NEW: test-suite/fittedfinancebenchmark.cpp
   Benchmark comparing:
   • Standard CN vs Fitted Implicit vs Rannacher+CN
   • Grid sizes: 50², 100², 200², 400²
   • Metrics: wall-clock time, max error vs analytic, oscillation count
   • Output as a formatted table to BOOST_TEST_MESSAGE

   Include timing for:
   • fitting factor computation overhead (coth evaluation)
   • tridiagonal solve time (should be identical)
   • total engine time

2. MODIFIED: ql/experimental/finitedifferences/fdmblackscholesfwdop.hpp/.cpp
   • Add optional exponential fitting to the forward operator
   • The forward operator discretizes: 
       ∂p/∂t = -∂/∂x[(r-q-½σ²)p] + ½∂²/∂x²[σ²p]
     which is the adjoint of the backward operator
   • Fitting applies to the diffusion term in the same way

3. UPDATE: ql/experimental/finitedifferences/all.hpp
   • Include any new experimental headers

4. MODIFIED: Makefile.am / CMakeLists.txt at ql/ level
   • Register new .cpp files in the build:
     ql/methods/finitedifferences/operators/fdmfittedblackscholesop.cpp
     ql/methods/finitedifferences/bsmfittedoperator.cpp
   • Register new test files
```

---

## Summary of ALL Files Across ALL Rounds

Here's the master checklist showing every file that needs to be produced or modified, with its round assignment:

| Round | Action | File |
|---|---|---|
| 1 | NEW | `operators/fdmfittedblackscholesop.hpp` |
| 1 | NEW | `operators/fdmfittedblackscholesop.cpp` |
| 1 | MOD | `operators/all.hpp` |
| 2 | MOD | `solvers/fdmbackwardsolver.hpp` |
| 2 | MOD | `solvers/fdmbackwardsolver.cpp` |
| 2 | NEW | `schemes/richardsonextrapolationscheme.hpp` |
| 2 | MOD | `schemes/all.hpp` |
| 3 | MOD | `solvers/fdmsolverdesc.hpp` |
| 3 | MOD | `solvers/fdmblackscholessolver.hpp` |
| 3 | MOD | `solvers/fdmblackscholessolver.cpp` |
| 4 | MOD | `pricingengines/vanilla/fdblackscholesvanillaengine.hpp` |
| 4 | MOD | `pricingengines/vanilla/fdblackscholesvanillaengine.cpp` |
| 4 | MOD | `pricingengines/barrier/fdblackscholesbarrierengine.hpp` |
| 4 | MOD | `pricingengines/barrier/fdblackscholesbarrierengine.cpp` |
| 4 | MOD | `pricingengines/barrier/fdblackscholesrebateengine.hpp` |
| 4 | MOD | `pricingengines/barrier/fdblackscholesrebateengine.cpp` |
| 5 | NEW | `bsmfittedoperator.hpp` |
| 5 | NEW | `bsmfittedoperator.cpp` |
| 5 | MOD | `pricingengines/vanilla/fdvanillaengine.hpp` |
| 5 | MOD | `pricingengines/vanilla/fdvanillaengine.cpp` |
| 5 | MOD | `all.hpp` (root FD) |
| 6 | NEW | `test-suite/fittedfinitedifferences.cpp` |
| 6 | NEW | `test-suite/fittedfinitedifferences.hpp` |
| 6 | MOD | `test-suite/fdmlinearop.cpp` |
| 6 | MOD | `test-suite/CMakeLists.txt` |
| 7 | NEW | `test-suite/fittedfinancebenchmark.cpp` |
| 7 | MOD | `experimental/finitedifferences/fdmblackscholesfwdop.hpp` |
| 7 | MOD | `experimental/finitedifferences/fdmblackscholesfwdop.cpp` |
| 7 | MOD | `CMakeLists.txt` / `Makefile.am` (top-level build) |

The `experimental/finitedifferences/` directory is intentionally last because it's the lowest priority—the forward PDE operators are used for calibration and density estimation rather than direct pricing, and the fitted scheme's benefits are most impactful on the backward pricing PDE.

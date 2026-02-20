## ROUND 0 — Preflight Audit Report (no code)

This report answers the requested audit questions **using only the provided QuantLib source files** and then translates the Milev–Tagliani / Duffy formulas into **QuantLib’s verified log-space convention**. Finally, it proposes a file-level implementation roadmap consistent with the preferred “spatial discretization selector” API.

---

# Part A — Source Audit Questions (with source citations)

> **Notation for citations:** I cite by **file path + identifier/function + key statement(s)**, since line numbers are not available in the pasted sources.

---

## 1) COORDINATE CONVENTION

### 1a) What does `FdmBlackScholesMesher` store as spatial coordinate \(x\)?

**Answer: it stores \(x = \ln(S)\)** (not \(\ln(S/K)\)).

**Evidence:**
- The header explicitly says “(in ln(S))”:  
  `ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp` — file brief: *“1-d mesher for the Black-Scholes process (in ln(S))”*.
- The implementation computes bounds as logs of forward minima/maxima and stores them directly as `locations_`:  
  `ql/methods/finitedifferences/meshers/fdmblackscholesmesher.cpp`, `FdmBlackScholesMesher::FdmBlackScholesMesher`:
  - `Real xMin = std::log(mi) - ...; Real xMax = std::log(ma) + ...;`
  - `helper = new Uniform1dMesher(xMin, xMax, size)` (or `Concentrating1dMesher(xMin, xMax, ...)`)
  - `locations_ = helper->locations();`

Also, downstream solver usage assumes the solver’s \(x\) is \(\ln(S)\):  
`ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp`, `FdmBlackScholesSolver::valueAt` returns `solver_->interpolateAt(std::log(s))`.

---

### 1b) How does `FdmBlackScholesOp` interpret `mesher->location()` / `mesher->locations()` in `setTime()`?

`FdmBlackScholesOp` **treats the mesher coordinate as log-spot \(x=\ln(S)\)** and converts to spot **only when it needs spot** (local-vol query).

**Evidence:**
- In the constructor, for the local-vol case it precomputes spot values as `Exp(mesher->locations(direction))`:  
  `ql/methods/finitedifferences/operators/fdmblackscholesop.cpp`, `FdmBlackScholesOp::FdmBlackScholesOp` initializer:  
  `x_((localVol) ? Array(Exp(mesher->locations(direction))) : Array())`  
  This is a direct “\(S=\exp(x)\)” conversion, proving the stored grid is \(x=\ln(S)\).
- In `setTime()` (local-vol branch), it calls `localVol_->localVol(..., x_[i], ...)` where `x_[i]` is the **spot** (not log):  
  same file, `FdmBlackScholesOp::setTime`, inside the `if (localVol_ != nullptr)` loop:  
  `v[i] = squared(localVol_->localVol(0.5*(t1+t2), x_[i], true));`

In the **constant-vol branch**, it does not use `mesher->location()` explicitly; it uses `strike_` to query forward variance:  
`volTS_->blackForwardVariance(t1, t2, strike_)/(t2-t1)` in `FdmBlackScholesOp::setTime`.

---

### 1c) If the mesher uses \(\ln(S/K)\), how is the strike communicated?

Not applicable here: **the mesher does not use \(\ln(S/K)\)**.

However, the strike is still communicated to:
- The mesher constructor (used for vol/bounds calculations):  
  `FdmBlackScholesMesher(Size,..., Time maturity, Real strike, ...)` in  
  `ql/methods/finitedifferences/meshers/fdmblackscholesmesher.hpp`.
- The operator (`FdmBlackScholesOp`) to query forward variance at the strike:  
  `FdmBlackScholesOp(..., Real strike, ...)` and `strike_` usage in  
  `ql/methods/finitedifferences/operators/fdmblackscholesop.hpp/.cpp`.

---

## 2) DERIVATIVE OPERATORS

### 2a) How do `FirstDerivativeOp` / `SecondDerivativeOp` compute their stencils?

They build **row-wise tridiagonal stencils** using the local mesh steps:
- `hm = mesher->dminus(iter, direction_)`
- `hp = mesher->dplus(iter, direction_)`

**Evidence:**
- `ql/methods/finitedifferences/operators/firstderivativeop.cpp`, `FirstDerivativeOp::FirstDerivativeOp`: reads `hm`, `hp` and constructs coefficients.
- `ql/methods/finitedifferences/operators/secondderivativeop.cpp`, `SecondDerivativeOp::SecondDerivativeOp`: same pattern.

Both use the same helper denominators:
- `zetam1 = hm*(hm+hp)`
- `zeta0  = hm*hp`
- `zetap1 = hp*(hm+hp)`

---

### 2b) What spacing functions do they use?

They use the **mesher’s** directional spacing:
- `FdmMesher::dplus(iter, direction)`
- `FdmMesher::dminus(iter, direction)`

**Evidence:**
- `FirstDerivativeOp` and `SecondDerivativeOp` constructors call `mesher->dminus(iter, direction_)` and `mesher->dplus(iter, direction_)`.  
  (`ql/methods/finitedifferences/operators/firstderivativeop.cpp`,  
  `ql/methods/finitedifferences/operators/secondderivativeop.cpp`)
- For composite meshers, those are forwarded to the underlying `Fdm1dMesher`:  
  `ql/methods/finitedifferences/meshers/fdmmeshercomposite.cpp`, `FdmMesherComposite::dplus/dminus`.

---

### 2c) Exact stencil (interior vs boundary)

#### First derivative \(u_x\)

Let node index be \(j\), with \(h_- = hm\), \(h_+ = hp\).

**Interior nodes** (`0 < j < last`):  
From `firstderivativeop.cpp`:
- `lower = -hp / (hm*(hm+hp))`
- `diag  = (hp-hm) / (hm*hp)`
- `upper =  hm / (hp*(hm+hp))`

So:
\[
u_x(x_j) \approx
\left(-\frac{h_+}{h_-(h_-+h_+)}\right)u_{j-1}
+\left(\frac{h_+ - h_-}{h_-h_+}\right)u_j
+\left(\frac{h_-}{h_+(h_-+h_+)}\right)u_{j+1}.
\]

**Lower boundary** (`j==0`): forward / upwind:
- `lower = 0`
- `diag = -(upper = 1/hp)`  
So \(u_x \approx (u_{1}-u_{0})/h_+\).  
Source: `ql/.../firstderivativeop.cpp`, branch `if (iter.coordinates()[direction_] == 0)`.

**Upper boundary** (`j==last`): backward / downwind:
- `upper = 0`
- `lower = -(diag = 1/hm)`  
So \(u_x \approx (u_{last}-u_{last-1})/h_-\).  
Source: same file, `else if (co == dim-1)`.

#### Second derivative \(u_{xx}\)

**Interior nodes** (`0 < j < last`):  
From `secondderivativeop.cpp`:
- `lower =  2/(hm*(hm+hp))`
- `diag  = -2/(hm*hp)`
- `upper =  2/(hp*(hm+hp))`

So:
\[
u_{xx}(x_j) \approx
\frac{2}{h_-(h_-+h_+)}u_{j-1}
-\frac{2}{h_-h_+}u_j
+\frac{2}{h_+(h_-+h_+)}u_{j+1}.
\]

**Boundary nodes** (`j==0` or `j==last`): all coefficients set to 0 (operator returns 0 there):  
`ql/.../secondderivativeop.cpp`, `if (co == 0 || co == dim-1) lower=diag=upper=0.0;`

---

## 3) OPERATOR ASSEMBLY

### 3a) What does `TripleBandLinearOp::axpyb(a, op1, op2, rhs)` compute?

It computes a **row-wise affine combination** of two operators plus a **diagonal shift**:

For each row \(i\):
- `lower[i] = y_lower[i] + s * x_lower[i]`
- `diag[i]  = y_diag[i]  + s * x_diag[i] + b_i`
- `upper[i] = y_upper[i] + s * x_upper[i]`

where \(s\) is taken from `a` (either a scalar `a[0]` if `a.size()==1`, or `a[i]` if `a.size()>1`), and \(b_i\) similarly from `b`.

**Evidence:**  
`ql/methods/finitedifferences/operators/triplebandlinearop.cpp`, `TripleBandLinearOp::axpyb` (all 3 cases of empty/nonempty `a`,`b` show the same pattern).

Conceptually:
\[
\text{this} \;\leftarrow\; y \;+\; \operatorname{diag}(a)\,x \;+\; \operatorname{diag}(b).
\]

This is exactly how `FdmBlackScholesOp` builds:
\[
L = b(x,t) \, D_x \;+\; a(x,t)\, D_{xx} \;-\; r\,I
\]
using `dxMap_`, `dxxMap_` and the diagonal `Array(1,-r)`.

---

### 3b) Exact index conventions for `lower_[]`, `diag_[]`, `upper_[]`

For each layout index \(i\), `apply()` computes:
\[
(\text{op}\cdot r)_i = \text{lower}[i]\cdot r[i0[i]] + \text{diag}[i]\cdot r[i] + \text{upper}[i]\cdot r[i2[i]].
\]

**Evidence:**  
`ql/methods/finitedifferences/operators/triplebandlinearop.cpp`, `TripleBandLinearOp::apply`:
`retVal[i] = r[i0ptr[i]]*lptr[i] + r[i]*dptr[i] + r[i2ptr[i]]*uptr[i];`

The neighbor indices `i0_[i]` and `i2_[i]` are set as:
- `i0_[i] = layout->neighbourhood(iter, direction, -1)`
- `i2_[i] = layout->neighbourhood(iter, direction, +1)`  
Source: `ql/.../triplebandlinearop.cpp`, constructor.

`FdmLinearOpLayout::neighbourhood` uses **reflection** outside the boundary (mirror), not clamping:
- if coordinate becomes `-1`, it is reflected to `+1`
- if coordinate becomes `dim`, it is reflected to `dim-2`  
Source: `ql/methods/finitedifferences/operators/fdmlinearoplayout.cpp`, `FdmLinearOpLayout::neighbourhood`.

This reflection is why derivative operators must explicitly zero-out unused bands at boundaries (which they do).

---

### 3c) How does `TripleBandLinearOp::mult(Array u)` work?

`mult(u)` scales **each row’s three coefficients** by `u[i]`:
- `ret.lower[i] = lower[i] * u[i]`
- `ret.diag[i]  = diag[i]  * u[i]`
- `ret.upper[i] = upper[i] * u[i]`

**Evidence:**  
`ql/methods/finitedifferences/operators/triplebandlinearop.cpp`, `TripleBandLinearOp::mult`.

This is a **left-multiplication by a diagonal matrix** \(\operatorname{diag}(u)\) in operator form, consistent with the comment in the header:
`// interpret u as the diagonal of a diagonal matrix, multiplied on LHS`  
(`ql/.../triplebandlinearop.hpp`).

---

## 4) TIME-STEPPING (MODERN FRAMEWORK)

### 4a) How does `FdmBackwardSolver::rollback()` handle stopping times? Does it segment rollback?

Yes. Stopping times are handled by `FiniteDifferenceModel::rollbackImpl()`, which **splits a time step** if it crosses any stopping time(s), performing sub-steps that land **exactly** on each stopping time and applying the step condition there.

**Evidence:**
- `ql/methods/finitedifferences/solvers/fdmbackwardsolver.cpp`, `FdmBackwardSolver::rollback`: for CN it constructs  
  `FiniteDifferenceModel<CrankNicolsonScheme> cnModel(..., condition_->stoppingTimes());` then `cnModel.rollback(...)`.
- `ql/methods/finitedifferences/finitedifferencemodel.hpp`, `FiniteDifferenceModel::rollbackImpl`:
  - loops over steps and computes `now` and `next`
  - scans stopping times and if `next <= stoppingTime < now` it:
    - sets a smaller step `evolver_.setStep(now-stoppingTime)`
    - calls `evolver_.step(a, now)`
    - calls `condition->applyTo(a, stoppingTime)`
  - then potentially completes the remainder to `next` and applies condition at `next`
  - resets `evolver_.setStep(dt)` afterwards

So rollback is explicitly **segmented** around stopping times.

---

### 4b) Is `applyTo` called at the initial `from` time?

**Only if** `from` equals the **largest** stopping time.

**Evidence:**  
`ql/.../finitedifferencemodel.hpp`, `rollbackImpl`:
```cpp
if(!stoppingTimes_.empty() && stoppingTimes_.back() == from) {
    if (condition) condition->applyTo(a,from);
}
```

So:
- If `from` is a stopping time and is the maximum stopping time, `applyTo(a, from)` is executed before the first PDE step.
- Otherwise, the first call is at the end of the first step (at `next`) or at an intermediate hit.

---

### 4c) How are scheme instances created inside rollback? Do they persist across segments?

- In `FdmBackwardSolver::rollback`, evolvers are created as **local objects by value** for most schemes:  
  e.g., `CrankNicolsonScheme cnEvolver(...)` (stack object).  
  Source: `ql/.../fdmbackwardsolver.cpp`, CN case.
- Each evolver is then stored **by value** inside `FiniteDifferenceModel` (moved into `evolver_`).  
  Source: `ql/.../finitedifferencemodel.hpp`, constructor `FiniteDifferenceModel(Evolver evolver, ...) : evolver_(std::move(evolver))`.
- Inside a single `FiniteDifferenceModel::rollbackImpl`, the evolver **persists across all segments**; only `setStep()` is called to change `dt_` for sub-steps.

For CN specifically:
- `CrankNicolsonScheme` owns `explicit_` and `implicit_` as `ext::shared_ptr` members created in its constructor via `ext::make_shared`.  
  Source: `ql/methods/finitedifferences/schemes/cranknicolsonscheme.cpp`, `CrankNicolsonScheme::CrankNicolsonScheme`.

Thus, within one rollback call:
- the scheme instance persists
- the internal explicit/implicit sub-schemes persist
- `setStep()` is called repeatedly (default step + sub-steps around stopping times)

---

## 5) SCHEME CONFIGURATION

### 5a) Are `FdmSchemeDesc` members const?

Yes—`type`, `theta`, and `mu` are declared `const`.

**Evidence:**  
`ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp`, `struct FdmSchemeDesc`:
```cpp
const FdmSchemeType type;
const Real theta, mu;
```

---

### 5b) How does `FdmBackwardSolver::rollback` handle damping steps?

If `dampingSteps != 0` and the chosen scheme is **not** `ImplicitEulerType`, it performs:
1) a **damping phase** using `ImplicitEulerScheme` from `from` down to `dampingTo`, for `dampingSteps`
2) the requested scheme from `dampingTo` down to `to`, for `steps`

If the scheme is already `ImplicitEulerType`, it does one rollback with `allSteps = steps + dampingSteps`.

**Evidence:**  
`ql/methods/finitedifferences/solvers/fdmbackwardsolver.cpp`, `FdmBackwardSolver::rollback`:
- computes `allSteps`, `dampingTo`
- runs damping model with `ImplicitEulerScheme` when applicable
- switch statement runs the main scheme afterwards (CN case uses `dampingTo` as the start time)

---

### 5c) How does `CrankNicolsonScheme` combine explicit and implicit sub-steps?

`CrankNicolsonScheme::step(a,t)` does:
- explicit Euler step with weight \(1-\theta\)
- implicit Euler step with weight \(\theta\)

**Evidence:**  
`ql/methods/finitedifferences/schemes/cranknicolsonscheme.cpp`, `CrankNicolsonScheme::step`:
```cpp
if (theta_ != 1.0) explicit_->step(a, t, 1.0-theta_);
if (theta_ != 0.0) implicit_->step(a, t, theta_);
```

For the standard CN setting `theta_=0.5`, this is the canonical:
\[
(I - \tfrac{1}{2}\Delta t L)u^{n+1} = (I + \tfrac{1}{2}\Delta t L)u^n
\]
form (implemented via “explicit then implicit”).

---

### 5d) Does `ImplicitEulerScheme` have a public `setStep(Time dt)`? Does setting it affect the explicit sub-scheme?

- Yes, `ImplicitEulerScheme::setStep(Time dt)` is public and sets its own `dt_`.  
  Evidence: `ql/.../impliciteulerscheme.hpp` declares `void setStep(Time dt);`, and `impliciteulerscheme.cpp` defines `dt_=dt;`.
- Calling it affects **only that `ImplicitEulerScheme` instance**. The explicit scheme has its own `dt_` member and `setStep`.  
  Evidence: `ql/.../expliciteulerscheme.hpp/.cpp` and `ql/.../impliciteulerscheme.hpp/.cpp` each store and set `dt_` independently.

`CrankNicolsonScheme::setStep` explicitly sets both:
`explicit_->setStep(dt_); implicit_->setStep(dt_);`  
Source: `ql/.../cranknicolsonscheme.cpp`, `setStep`.

---

## 6) SIGN CONVENTIONS

### 6a) For the operator \(L\) stored in `mapT_`, what sign convention do off-diagonals use?

`mapT_` is a `TripleBandLinearOp`. In `apply()`, its coefficients multiply neighbor values as:

- `lower_[i]` multiplies \(u\) at the **“minus one”** neighbor index `i0_[i]`
- `upper_[i]` multiplies \(u\) at the **“plus one”** neighbor index `i2_[i]`

So the row is:
\[
(Lu)_i = \text{lower}[i]\cdot u_{i0[i]} + \text{diag}[i]\cdot u_i + \text{upper}[i]\cdot u_{i2[i]}.
\]

**Evidence:**  
`ql/.../triplebandlinearop.cpp`, `TripleBandLinearOp::apply`.

---

### 6b) For system matrix \((I - \theta\,\Delta t\,L)\), what sign must off-diagonals of \(L\) have for an M-matrix?

In this framework, the **system matrix** used by implicit stepping is:
\[
A = I - \theta\,\Delta t\,L,
\]
so its off-diagonals are:
\[
A_{i,i\pm1} = -\theta\,\Delta t\,L_{i,i\pm1}.
\]

To make \(A\) an M-matrix (non-positive off-diagonals), we need:
\[
L_{i,i\pm1} \ge 0 \quad\Rightarrow\quad A_{i,i\pm1} \le 0
\]
because \(\theta\,\Delta t > 0\).

This aligns directly with the `TripleBandLinearOp::apply` convention above.

(You also typically need \(A\) diagonally dominant with positive diagonal; for the Black–Scholes operator the diagonal dominance is usually ensured once off-diagonals are non-negative, because the derivative stencils have zero row-sum and the reaction term contributes \(-r\) to the diagonal—see Part B.10.)

---

### 6c) Verify with standard `FdmBlackScholesOp`: are its off-diagonals always non-negative? Under what grid conditions?

**They are not guaranteed to be non-negative.** The drift term contribution (via `FirstDerivativeOp`) can make one off-diagonal negative unless diffusion dominates (i.e., small Péclet number / sufficiently fine grid).

**Evidence from stencil signs:**
- `SecondDerivativeOp` interior off-diagonals are strictly positive (`lower_>0`, `upper_>0`).  
  Source: `ql/.../secondderivativeop.cpp`, interior branch uses `2/zetam1`, `2/zetap1`.
- `FirstDerivativeOp` interior has `lower_<0`, `upper_>0` (for a standard increasing mesh).  
  Source: `ql/.../firstderivativeop.cpp`, interior branch sets `lower_ = -hp/zetam1` and `upper_ = hm/zetap1`.

**Evidence of how the BS operator combines them:**
- `FdmBlackScholesOp::setTime` assembles `mapT_` as:
  - drift coefficient \(b = r-q-0.5v\) multiplying `dxMap_`
  - diffusion coefficient \(a = 0.5v\) multiplying `dxxMap_`
  - plus diagonal shift `-r`  
  Source: `ql/.../fdmblackscholesop.cpp`, `mapT_.axpyb(Array(1, r-q-0.5*v), dxMap_, dxxMap_.mult(0.5*v), Array(1,-r))`.

**Uniform-grid condition (1D, constant coefficients):**  
For uniform \(h\), interior stencils reduce to:
- \(D_x\): lower \(=-\frac{1}{2h}\), upper \(=+\frac{1}{2h}\)
- \(D_{xx}\): lower \(=+\frac{1}{h^2}\), upper \(=+\frac{1}{h^2}\)

So the BS operator off-diagonals are:
\[
L_{\text{lower}} = \frac{a}{h^2} - \frac{b}{2h},\qquad
L_{\text{upper}} = \frac{a}{h^2} + \frac{b}{2h}.
\]

To have **both** \(\ge 0\), you need:
\[
\frac{a}{h^2} \ge \frac{|b|}{2h}
\quad\Longleftrightarrow\quad
h \le \frac{2a}{|b|}
\quad\Longleftrightarrow\quad
|Pe|=\left|\frac{b h}{2a}\right|\le 1.
\]

This is exactly the convection–diffusion “central differencing monotonicity” restriction; when violated, CN can oscillate even though it remains stable—consistent with the paper motivation.

On **nonuniform** meshes, the condition becomes **local** (node-dependent) and depends on both \(h_-\) and \(h_+\) via the stencils in Part A.2c.

---

## 7) SMART POINTERS

### 7a) Is `ext::make_shared<T>(...)` available?

Yes. It is used throughout the provided code.

**Evidence examples:**
- `ql/methods/finitedifferences/meshers/fdmblackscholesmesher.cpp`: constructs `QuantoTermStructure` with `ext::make_shared<QuantoTermStructure>(...)`.
- `ql/methods/finitedifferences/schemes/cranknicolsonscheme.cpp`: `explicit_(ext::make_shared<ExplicitEulerScheme>(...))`.
- `ql/methods/finitedifferences/solvers/fdmblackscholessolver.cpp`: `ext::make_shared<FdmBlackScholesOp>(...)`.

---

### 7b) If not, what fallback pattern is correct?

The codebase already uses the canonical QuantLib fallback:
\[
\texttt{ext::shared_ptr<T>(new T(args...))}
\]

**Evidence examples:**
- `ql/methods/finitedifferences/meshers/fdmblackscholesmesher.cpp`:  
  `helper = ext::shared_ptr<Fdm1dMesher>(new Uniform1dMesher(...));`
- `ql/pricingengines/barrier/fdblackscholesrebateengine.cpp`:  
  multiple `new` usages wrapped in `ext::shared_ptr` / `FdmBoundaryConditionSet::value_type(new ...)`.

---

## 8) GRID ALIGNMENT

### 8a) Does `FdmBlackScholesMesher` place a grid node exactly at \(\ln(K)\)?

**Not guaranteed**.

What it does:
- It can use `Concentrating1dMesher` around a critical point, but in the way it is called from `FdmBlackScholesMesher`, it does **not** request that the point be an **exact node**.

**Evidence:**
- `FdmBlackScholesMesher` chooses:
  - `Concentrating1dMesher(xMin, xMax, size, pair(log(cPoint.first), cPoint.second))`  
    Source: `ql/.../fdmblackscholesmesher.cpp`, constructor.
- The `Concentrating1dMesher` constructor signature shows `requireCPoint=false` by default:  
  `ql/methods/finitedifferences/meshers/concentrating1dmesher.hpp`:
  ```cpp
  Concentrating1dMesher(..., const std::pair<Real,Real>& cPoints, bool requireCPoint=false);
  ```
- `FdmBlackScholesMesher` does not pass `requireCPoint=true`, so it uses default `false`.

Therefore, you get **concentration around** \(\ln(K)\), but not a hard guarantee that \(\ln(K)\) is on the grid.

---

### 8b) Can the mesher be configured to place nodes at arbitrary points (e.g., barrier levels \(\ln(L)\), \(\ln(U)\))?

Yes, in two different senses:

1) **Domain boundary alignment (exact)**: `FdmBlackScholesMesher` supports `xMinConstraint` / `xMaxConstraint` that directly set the domain endpoints.

**Evidence:**
- `FdmBlackScholesMesher` constructor parameters include `Real xMinConstraint`, `Real xMaxConstraint`.  
  Source: `ql/.../fdmblackscholesmesher.hpp`.
- Implementation overrides `xMin`/`xMax` when constraints are not Null:  
  `ql/.../fdmblackscholesmesher.cpp`:
  - `if (xMinConstraint != Null<Real>()) xMin = xMinConstraint;`
  - `if (xMaxConstraint != Null<Real>()) xMax = xMaxConstraint;`
- Barrier engines use this to force the boundary at \(\ln(\text{barrier})\):  
  `ql/pricingengines/barrier/fdblackscholesbarrierengine.cpp`, `calculate()` sets:
  - `xMin = std::log(arguments_.barrier)` for down barriers
  - `xMax = std::log(arguments_.barrier)` for up barriers
  and passes them to `FdmBlackScholesMesher`.

2) **Interior point alignment (possible, but not via `FdmBlackScholesMesher` API)**: you can directly build a `Concentrating1dMesher` with multiple “required points” (see 8c) and wrap it in `FdmMesherComposite`.

---

### 8c) How does `Concentrating1dMesher` handle multiple concentration points?

The second constructor supports a vector of tuples:
\[
(\text{point},\ \text{density},\ \text{required\_bool})
\]
and it **ensures required points become part of the grid** by:
- building an ODE-based mapping \(x \mapsto y\)
- then constructing a **transform interpolation** that “snaps” required points onto grid coordinates (it uses Brent root finding to locate the transform parameter that yields the required physical point)

**Evidence:**  
`ql/methods/finitedifferences/meshers/concentrating1dmesher.hpp` declares:
```cpp
Concentrating1dMesher(..., const std::vector<std::tuple<Real, Real, bool>>& cPoints, Real tol=1e-8);
```
Implementation details in `ql/.../concentrating1dmesher.cpp` (second constructor):
- collects `points` and `betas`
- solves for scaling `a` via `Brent().solve(...)`
- solves ODE at all grid points
- for each `required` point, finds its location and inserts it into the transform map `w`
- builds `LinearInterpolation transform(u,z)`
- sets `locations_[i] = odeSolution(transform(i*dx))`

So it can enforce **multiple required alignment points**—but this functionality is not currently exposed by `FdmBlackScholesMesher`.

---

## 9) S_MAX SELECTION

### 9a) How does current `FdmBlackScholesSolver` determine domain bounds?

`FdmBlackScholesSolver` **does not determine bounds itself**; it takes a pre-built mesher from `FdmSolverDesc`.

**Evidence:**
- `FdmBlackScholesSolver` stores `FdmSolverDesc solverDesc_` and in `performCalculations()` it constructs the operator using `solverDesc_.mesher`.  
  Source: `ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp/.cpp`.
- The bounds are decided upstream when engines construct `FdmBlackScholesMesher`, which computes:
  - forward envelope \([mi, ma]\) (with dividends)
  - then
    \[
    x_{\min}=\ln(mi)-\sigma\sqrt{T}\,\Phi^{-1}(1-\varepsilon)\cdot \text{scaleFactor},\quad
    x_{\max}=\ln(ma)+\sigma\sqrt{T}\,\Phi^{-1}(1-\varepsilon)\cdot \text{scaleFactor}
    \]
  Source: `ql/.../fdmblackscholesmesher.cpp`, `sigmaSqrtT`, `normInvEps`, `xMin`, `xMax`.

So “\(S_{\max}\)” in S-space corresponds to \(x_{\max}\) in log-space, and it comes from `FdmBlackScholesMesher`.

---

### 9b) Is there a way to override automatic \(S_{\max}\) selection?

Yes.

Options visible in the provided code:
1) Provide `xMinConstraint` / `xMaxConstraint` to the mesher (exact override).  
   Source: `ql/.../fdmblackscholesmesher.hpp/.cpp` and barrier engines using it.
2) Use a different mesher altogether (since `FdmSolverDesc` takes `const ext::shared_ptr<FdmMesher> mesher;`).  
   Source: `ql/.../fdmsolverdesc.hpp`.

---

# Part B — Coordinate Translation & Formula Mapping (log-space, QuantLib conventions)

## Verified QuantLib PDE form

From `FdmBlackScholesOp::setTime`, the operator is assembled as:
- drift coefficient: \(b = r - q - \tfrac{1}{2}v\)
- diffusion coefficient: \(a = \tfrac{1}{2}v\)
- reaction term: \(-r\)

**Evidence:**  
`ql/.../fdmblackscholesop.cpp`, in the constant-vol case:
- `v = blackForwardVariance(...) / (t2-t1)`
- `mapT_.axpyb(Array(1, r - q - 0.5*v), dxMap_, dxxMap_.mult(0.5*Array(..., v)), Array(1, -r));`

Thus in log-space \(x=\ln(S)\), the backward-in-time pricing PDE is:
\[
-\frac{\partial u}{\partial t} + b\,\frac{\partial u}{\partial x} + a\,\frac{\partial^2 u}{\partial x^2} - r\,u = 0,
\]
with:
\[
a=\frac{\sigma^2}{2},\qquad b=(r-q)-\frac{\sigma^2}{2},
\]
where \(\sigma\) here is the (Black) volatility and \(\sigma^2=v\).

---

## 8) Translate paper’s S-space fitting factor (CE-7) to log-space \(Pe_j,\rho_j\)

### Starting point (generic convection–diffusion fitting)
Duffy-style exponential fitting for:
\[
-a_d\,u_{xx} + \mu\,u_x \quad\text{(or equivalently } a_d u_{xx} + \mu u_x \text{ depending on sign convention)}
\]
introduces the dimensionless parameter:
\[
Pe = \frac{\mu\,h}{2\,a_d},
\]
and fitting factor:
\[
\rho = Pe\coth(Pe),
\]
with the smooth limit \(\rho \to 1 + \frac{Pe^2}{3}\) as \(Pe\to 0\).

### Log-space identification for QuantLib
In QuantLib log-space:
- convection coefficient is \(\mu \equiv b\)
- diffusion coefficient is \(a_d \equiv a\)
- mesh step is \(h\) (uniform) or \(h_j\) (nonuniform)

So:
\[
Pe_j = \frac{b\,h_j}{2a},
\qquad
\rho_j = Pe_j\coth(Pe_j).
\]

### Nonuniform mesh note (QuantLib-consistent spacing)
QuantLib’s derivative operators use `dminus` and `dplus` at each node (Part A.2). For the fitted scheme on nonuniform meshes, a practical per-node effective step is:
\[
h_j := \tfrac{1}{2}\left(d^-_j + d^+_j\right),
\]
where \(d^-_j = \texttt{dminus}\) and \(d^+_j = \texttt{dplus}\), skipping boundary nodes where one of them is `Null<Real>()`.

That matches the *local-step* nature of QuantLib stencils (and avoids using a single global \(h\) on a concentrated mesh).

---

## 9) Translate CN-variant parameter choice (CE-19) to effective diffusion \(a_{\text{eff}}\)

### Given (paper)
\[
\omega = -\frac{r}{16\sigma^2}.
\]

The CN variant distributes the reaction term \(-r u\) to neighbors (a 6-node stencil across the two time levels), creating **off-diagonal reaction contributions**.

### Goal (QuantLib implementation via standard CN)
QuantLib’s standard CN does:
\[
P u^{n+1} = N u^n,\quad
P = I - \tfrac{1}{2}\Delta t\,L,\quad
N = I + \tfrac{1}{2}\Delta t\,L.
\]
(`CrankNicolsonScheme::step` = explicit then implicit, Part A.5c.)

So any “extra off-diagonal” we want in the **paper’s \(P\)** must come from \(L\) and then gets multiplied by \(-\tfrac{1}{2}\Delta t\).

### (a) Distributed reaction as an operator contribution
If the reaction \(-r u\) is replaced by:
\[
-r\left[\omega\,u_{j-1} + (1-2\omega)\,u_j + \omega\,u_{j+1}\right],
\]
then at each time level it contributes an **off-diagonal term**:
\[
L_{\text{react,off}} = -r\omega.
\]
Since \(\omega<0\), we have \(-r\omega>0\): it is diffusion-like and helps monotonicity.

### (b) Why standard CN halves the off-diagonal contribution
In standard CN:
\[
P_{\text{off}} = -\tfrac{1}{2}\Delta t\,L_{\text{off}}.
\]
But the paper’s formulation applies the full \(\omega\)-weight at each time level (not “halved again” beyond the CN midpoint averaging).

Therefore, to reproduce the paper’s \(P\) using **standard CN**, the operator must contain **twice** the desired off-diagonal contribution:
\[
L_{\text{off, add}} = 2(-r\omega).
\]

### (c) Convert the doubled off-diagonal into “extra diffusion”
On a uniform log-grid with spacing \(h\), the diffusion part contributes off-diagonal:
\[
L_{\text{diff,off}} = \frac{a}{h^2}.
\]
To add \(L_{\text{off, add}}\), we can add \(a_{\text{add}}\) such that:
\[
\frac{a_{\text{add}}}{h^2} = 2(-r\omega).
\]

Now plug \(\omega = -\frac{r}{16\sigma^2}\):
\[
2(-r\omega) = 2\left(-r\left(-\frac{r}{16\sigma^2}\right)\right)
= \frac{r^2}{8\sigma^2}.
\]
Hence:
\[
a_{\text{add}} = \frac{r^2 h^2}{8\sigma^2}.
\]

So the effective diffusion is:
\[
a_{\text{eff}} = a + a_{\text{add}}
= \frac{\sigma^2}{2} + \frac{r^2 h^2}{8\sigma^2}.
\]

This is the “**\(8\)** in the denominator” result (not \(16\)).

### (d) Verify \(P_{\text{lower}}\) matches the paper’s extra term
With standard CN:
\[
P_{\text{off, add}} = -\tfrac{1}{2}\Delta t\left(\frac{a_{\text{add}}}{h^2}\right)
= -\tfrac{1}{2}\Delta t\left(\frac{r^2}{8\sigma^2}\right)
= -\Delta t\,\frac{r^2}{16\sigma^2}.
\]
The paper’s distributed-reaction off-diagonal term in \(P\) is \(\Delta t\,(r\omega)\), and:
\[
\Delta t\,(r\omega)=\Delta t\left(r\left(-\frac{r}{16\sigma^2}\right)\right)
= -\Delta t\,\frac{r^2}{16\sigma^2}.
\]
They match exactly.

---

## 10) M-matrix condition in terms of `TripleBandLinearOp` coefficients after `axpyb()`

### Operator storage
`FdmBlackScholesOp` stores the spatial operator \(L\) inside `mapT_` (a `TripleBandLinearOp`).

`TripleBandLinearOp` defines the tridiagonal row \(i\) as:
- \(L_{i,i-1} \equiv \texttt{lower_[i]}\)
- \(L_{i,i}   \equiv \texttt{diag_[i]}\)
- \(L_{i,i+1} \equiv \texttt{upper_[i]}\)  
(see Part A.3b)

### System matrix for implicit stepping
Implicit stepping uses:
\[
A = I - \theta\,\Delta t\,L.
\]
So:
\[
A_{i,i\pm 1} = -\theta\,\Delta t\,L_{i,i\pm1}.
\]

### M-matrix off-diagonal condition (load-bearing, easy to check)
Since \(\theta\Delta t>0\), a **sufficient** and in practice decisive condition is:
\[
\boxed{L_{i,i-1} \ge 0 \ \text{and}\ L_{i,i+1} \ge 0\quad\forall i}
\]
i.e.:
\[
\boxed{\texttt{lower_[i] >= 0 \ \&\&\ upper_[i] >= 0}}
\]

This is exactly what your planned runtime diagnostic can test via `ModTripleBandLinearOp(mapT_)` (accessors `lower(i)`, `upper(i)` exist in `ql/.../modtriplebandlinearop.hpp`).

### Why diagonal dominance is then automatic (for the BS operator)
For QuantLib’s BS operator, both `FirstDerivativeOp` and `SecondDerivativeOp` have **row-sum = 0** (true for the shown boundary handling too), and the reaction contributes \(-r\) on the diagonal only. Therefore:
\[
\sum_j L_{i,j} = -r.
\]
Then:
\[
\sum_j A_{i,j} = 1 - \theta\Delta t(-r) = 1 + \theta\Delta t r > 0,
\]
and with \(A\) off-diagonals \(\le 0\), \(A\) is diagonally dominant with positive diagonal—supporting the M-matrix property.

(For CN positivity you also need the explicit-side matrix \(N\) to be nonnegative, which imposes a time-step restriction on the diagonal, but the **off-diagonal** nonnegativity condition is the same: \(L_{\text{off}}\ge 0\).)

---

## 11) Uniform vs nonuniform log-mesh: constancy of \(Pe\)

- On a **uniform** log grid (\(h_j \equiv h\)) with **constant** coefficients \(a,b\), the Péclet number is:
  \[
  Pe = \frac{b h}{2a},
  \]
  which is constant across all interior nodes.

- On a **nonuniform** grid (e.g., `Concentrating1dMesher`), `dplus/dminus` vary by node, hence \(h_j\) varies and:
  \[
  Pe_j = \frac{b h_j}{2a}
  \]
  varies per node.

- With **local volatility**, \(a=a(x,t)\) also varies per node and time, so \(Pe_j\) becomes both node- and time-dependent even on uniform grids.

**Impact on fitted-operator assembly:** you cannot compute a single global \(\rho\). You must compute \(\rho_j\) per interior node (and per time step if coefficients vary with time), consistent with QuantLib’s `setTime(t1,t2)` pattern.

---

## 12) Accuracy constraints in log-space; quantify for \(\sigma=0.001\), \(r=0.05\)

### Scheme 1 (exponentially fitted / upwind-limit diffusion)
In log-space, the dominant *numerical* diffusion term in the drift-dominated regime behaves like:
\[
a_{\text{num,1}} \sim \tfrac{1}{2}|b|\,h.
\]
To make it negligible vs physical diffusion \(a=\sigma^2/2\), require:
\[
\tfrac{1}{2}|b|h \ll a
\quad\Longrightarrow\quad
h \ll \frac{2a}{|b|}.
\]

With \(a=\sigma^2/2\) and \(b \approx r-q-\sigma^2/2\), for \(q=0\):
- \(\sigma=0.001 \Rightarrow \sigma^2 = 10^{-6}\)
- \(a = 5\times 10^{-7}\)
- \(b \approx 0.05 - 5\times 10^{-7} \approx 0.0499995\)

So:
\[
\frac{2a}{|b|} \approx \frac{10^{-6}}{0.05} = 2\times 10^{-5}.
\]
Hence **need \(h \ll 2\times 10^{-5}\)**.

### Scheme 2 (CN variant via effective diffusion)
Artificial diffusion term in log-space is:
\[
a_{\text{add}} = \frac{r^2 h^2}{8\sigma^2}.
\]
Require \(a_{\text{add}} \ll a=\sigma^2/2\):
\[
\frac{r^2 h^2}{8\sigma^2} \ll \frac{\sigma^2}{2}
\quad\Longrightarrow\quad
r^2 h^2 \ll 4\sigma^4
\quad\Longrightarrow\quad
h \ll \frac{2\sigma^2}{r}.
\]

For \(\sigma^2=10^{-6}\), \(r=0.05\):
\[
\frac{2\sigma^2}{r} = \frac{2\times 10^{-6}}{0.05} = 4\times 10^{-5}.
\]
Hence **need \(h \ll 4\times 10^{-5}\)**.

So for the very-low-vol test case, both schemes demand extremely fine log spacing to be “diffusion-free,” with Scheme 1 being slightly more restrictive in this parameter regime.

---

# Part C — File-Level Implementation Roadmap (preferred “spatial discretization selector” API)

You stated a preference for:

> **Add an explicit “spatial discretization” selector to `FdmBlackScholesSolver` (and engines), leaving `FdmSchemeDesc` strictly as time-stepping.**

Below is a practical multi-round roadmap consistent with the audited framework.

---

## Proposed public API object

Create a small config type, e.g.:

- `enum class FdmBlackScholesSpatialDiscretization { Standard, ExponentialFittingImplicit, CrankNicolsonVariant };`
- optionally a struct for per-scheme parameters/diagnostics toggles

**Rationale:** `FdmSchemeDesc` is already (and should remain) **time stepping only** (Part A.5a).

---

## Roadmap table

| ROUND | File (path) | New / Modify | Class / API | Purpose | Depends on |
|---:|---|---|---|---|---|
| 1 | `ql/methods/finitedifferences/solvers/fdmblackscholessolver.hpp/.cpp` | **Modify** | `FdmBlackScholesSolver` | Add spatial-discretization selector member + pass it into operator construction | Audit Part A confirms solver constructs `FdmBlackScholesOp` in `performCalculations()` |
| 1 | `ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp/.cpp` | **New** | `FdmFittedBlackScholesOp` (name TBD) | Implement **Scheme 1**: fitted stencil operator in log-space, compatible with `FdmLinearOpComposite` API; intended for **ImplicitEuler** stepping | Needs derivative/sign conventions (Part A.2, A.3, A.6) |
| 1 | `ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp/.cpp` | **Modify** | `FdBlackScholesVanillaEngine` + `MakeFdBlackScholesVanillaEngine` | Add parameter + builder method to select spatial discretization; thread through to `FdmBlackScholesSolver` | Depends on Round 1 solver API |
| 1 | `ql/pricingengines/barrier/fdblackscholesbarrierengine.hpp/.cpp` | **Modify** | `FdBlackScholesBarrierEngine` | Same threading as vanilla engine | Round 1 solver API |
| 1 | `ql/pricingengines/barrier/fdblackscholesrebateengine.hpp/.cpp` | **Modify** | `FdBlackScholesRebateEngine` | Same threading as barrier engine | Round 1 solver API |
| 2 | `ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.hpp/.cpp` | **New** | `FdmCnVariantBlackScholesOp` (name TBD) | Implement **Scheme 2** via **effective diffusion** \(a_{\text{eff}}\) assembly; meant to be used with standard `CrankNicolsonScheme` | Requires Part B.9 derivation + Part A.3 assembly semantics |
| 2 | `ql/pricingengines/vanilla/fdblackscholesshoutengine.hpp/.cpp` | **Modify** | `FdBlackScholesShoutEngine` | Thread spatial discretization into shout engine | Round 1 solver API |
| 2 | `ql/methods/finitedifferences/operators/*` (small helper header) | **New** | `FdmBlackScholesSpatialDiscretization` enum/struct | Centralize the selector type to avoid duplicating enum in engines/solver | None |
| 3 | (optional) `ql/methods/finitedifferences/operators/...` | **New** | M-matrix diagnostic helper | Add `mMatrixViolationCount()` by inspecting `mapT_` via `ModTripleBandLinearOp` (`modtriplebandlinearop.hpp`) | Part B.10 |
| 3 | (optional, outside given tree) `test-suite/...` | **New** | unit/regression tests | Verify positivity, no oscillations (esp. discontinuous payoff), and match baselines | Would require access to QuantLib test infrastructure (not provided here) |

### Estimated output size per ROUND (rough)
- **Round 1:** 4–6 files touched/added, ~800–1500 LOC total (new operator + API threading).
- **Round 2:** 2–4 files, ~500–1200 LOC (CN-variant operator + remaining engines).
- **Round 3 (optional):** diagnostics + tests; diagnostics ~100–300 LOC; tests can be 300–1500 LOC depending on coverage.

---

## Notes / constraints surfaced by the audit (implementation-impacting)

1) **Operator sign convention is “\(L\) goes into \((I-\theta\Delta t L)\)”**  
   The fitted/CN-variant operators must populate `mapT_` so that **`lower_[]` and `upper_[]` are non-negative** to make the system matrix an M-matrix (Part A.6).

2) **Stopping times are segmented automatically**  
   Discrete monitoring step conditions (barrier indicator) will be applied exactly at stopping times if implemented as an `FdmStepConditionComposite` component (Part A.4).

3) **Grid-point alignment is not guaranteed by current mesher defaults**  
   If strict barrier/strike alignment becomes a requirement, you’ll either:
   - extend `FdmBlackScholesMesher` to expose the multi-point `Concentrating1dMesher` constructor, or
   - build a custom mesher externally and pass it through `FdmSolverDesc`  
   (Part A.8).

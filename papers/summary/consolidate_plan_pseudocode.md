## SECTION 1: SYSTEM ARCHITECTURE OVERVIEW

### 1.1 Top-level module diagram (data + control flow; hot vs cold)

```text
Legend:
  [COLD] = called once / per pricing request
  [HOT]  = called per time step (or per damp substep)
  ─────▶  data flow (struct on arrow)
  ─ ─ ▶  control / events (flags, state)

┌──────────────────────────────────────────────────────────────────────────────────────────────┐
│ Inputs (external)                                                                              │
│  EngineConfig cfg        MarketData md         ProductSpec prod         ModelSpec model       │
│  (numerics+strategies)   (r,q,curves,spots)    (type, K,B,L/U,mon)      (PDE family flags)     │
└──────────────────────────────────────────────────────────────────────────────────────────────┘
                    │ cfg, md, prod, model
                    ▼
┌──────────────────────────────────────────────────────────────────────────────────────────────┐
│ M12 ┌──────────────────────────────────────────────────────────────────────────────────────┐ │
│     │ Orchestrator [COLD]                                                                  │ │
│     │ - validates inputs, selects strategies                                                │ │
│     │ - runs full solve lifecycle                                                           │ │
│     └──────────────────────────────────────────────────────────────────────────────────────┘ │
└──────────────────────────────────────────────────────────────────────────────────────────────┘
   │ StrategyBundle(sb)                       │ GridSpec(gs) + BoundarySpec(bc) + VolSurface(vol)
   │ (module strategy picks)                 │ (domain, J/N, mesh type, event times, BC style)
   ▼                                         ▼
┌───────────────────────────────┐        ┌───────────────────────────────┐
│ M2 CoordinateTransform [COLD] │        │ M4 VolSurfaceAdapter [HOT]     │
│  S ↔ x, Jacobians, maps       │        │  σ(S,τ) lookup/interp          │
└──────────────┬────────────────┘        └──────────────┬────────────────┘
               │ CoordinateTransform(xf)                 │ sigmaAt()
               ▼                                        │
┌──────────────────────────────────────────┐            │
│ M1 GridFactory [COLD]                    │            │
│  - SpatialGrid (uniform x or sinh-graded)│            │
│  - TemporalGrid (aligned to monitoring)  │            │
│  - barrier-node alignment                │            │
└──────────────┬───────────────────────────┘            │
               │ SpatialGrid(xg), TemporalGrid(tg)       │
               ▼                                        │
┌──────────────────────────────────────────┐            │
│ M3 PayoffProcessor [COLD]                │            │
│  - cell-averaged payoff init → Stage1 §3.4.1          │
│  - corner compatibility → Stage1 §3.4.2               │
└──────────────┬───────────────────────────┘            │
               │ SolutionState(st0) + payoffPhi[]        │
               ▼                                        │
┌──────────────────────────────────────────────────────────────────────────────────────────────┐
│ Time marching loop over TemporalGrid                                                          │
│                                                                                              │
│   for each step n: τ_n → τ_{n+1} with dt = tg.dt[n]                                           │
│                                                                                              │
│   [HOT path modules invoked inside the loop]                                                  │
└──────────────────────────────────────────────────────────────────────────────────────────────┘
               │ st (u[], τ), dt, τ_n, τ_{n+1}, prod, bc, vol, xf, xg
               ▼
┌───────────────────────────────┐          ┌───────────────────────────────┐
│ M6 BoundaryHandler [HOT]      │          │ M11 DiagnosticsEngine [HOT/*] │
│  - eval Dirichlet BCs         │          │  - positivity check           │
│  - apply corridor projection  │          │  - oscillation detector       │
│  - emits StepEvent            │          │  - M-matrix check             │
└──────────────┬────────────────┘          │  - step health summary        │
               │ BoundaryValues(bv),        └──────────────┬────────────────┘
               │ StepEvent(ev)                             │ DiagnosticsReport(dr_step)
               │                                           │
               ▼                                           │
┌──────────────────────────────────────────────────────────────────────────────────────────────┐
│ M7 TimeStepper [HOT]                                                                              │
│  Strategy: RannacherCN / TRBDF2 / PureImplicitEuler                                                │
│  - damping state machine (initial + post-monitor)                                                  │
│  - assembles systems via M5 and solves via M8 / M9                                                  │
│                                                                                                    │
│   ┌───────────────────────┐      OperatorCoeffs(op)      ┌───────────────────────┐               │
│   │ M5 SpatialOperator     │◀────────────────────────────▶│ M4 VolSurfaceAdapter  │               │
│   │   assemble tridiag Lh  │                              │   σ(S,τ)              │               │
│   └───────────┬───────────┘                              └───────────────────────┘               │
│               │ TridiagSystem(A,rhs)                                                          ┌───▼───────────────┐
│               ├──────────────────────────────────────────────────────────────────────────────▶│ M9 AmericanConstraint│
│               │                                                                               │   None/Policy/Pen  │
│               │                                                                               └───┬───────────────┘
│               │ TridiagSystem (modified) / u_next[]                                                │
│               ▼                                                                                    │
│   ┌───────────────────────┐                                                                       │
│   │ M8 TridiagSolver       │  (Thomas → Stage1 §3.6.1)                                             │
│   └───────────┬───────────┘                                                                       │
│               │ u_next[]                                                                           │
│               ▼                                                                                    │
│   updates SolutionState(st) + TimeStepperState(ts)  ─ ─▶ emits step diagnostics hooks               │
└───────────────────────────────────────────────────────────────────────────────────────────────────┘
               │ st_T (final), aggregated DiagnosticsReport(dr_all)
               ▼
┌───────────────────────────────┐
│ M10 GreeksEngine [COLD/*]     │
│  - Δ,Γ via log-derivatives → Stage1 §3.5.1/§3.5.2
│  - Θ via PDE residual (operator apply at τ=T)
│  - Vega/Rho via bump-and-reprice (calls M12)
└──────────────┬────────────────┘
               │ GreekBundle(gb)
               ▼
┌──────────────────────────────────────────────────────────────────────────────────────────────┐
│ Outputs                                                                                      │
│  PricingResult {price, Greeks, diagnostics, (optional) surfaces, convergence record}        │
└──────────────────────────────────────────────────────────────────────────────────────────────┘

Optional (offline / QA):
  M11.RichardsonEstimator (uses multiple PricingResult runs) ─────▶ ConvergenceRecord
```

---

### 1.2 Module registry

| Module ID | Name | Responsibility | Hot/Cold | Depends On |
|---|---|---|---|---|
| M1 | GridFactory | Build SpatialGrid/TemporalGrid; sinh grading; monitoring alignment; barrier-node alignment | Cold | M2, ProductSpec, EngineConfig |
| M2 | CoordinateTransform | S↔x mapping + Jacobians for Greeks | Cold | — |
| M3 | PayoffProcessor | Cell-averaged payoff init; payoff library; corner compatibility | Cold | M2, M6 (BC at τ=0), ProductSpec |
| M4 | VolSurfaceAdapter | Provide σ(S,τ) lookups (constant/local-vol grid) | Hot | VolSurface, SpatialGrid |
| M5 | SpatialOperator | Assemble fitted FD tridiagonal operator; M-matrix checks; fallback upwind | Hot | M4, SpatialGrid, ModelSpec, EngineConfig |
| M6 | BoundaryHandler | Evaluate BCs; apply discrete monitoring projection; emit discontinuity events | Hot | ProductSpec, BoundarySpec, SpatialGrid |
| M7 | TimeStepper | Time-marching (RS-CN, TR-BDF2, Euler); damping state machine; system assembly orchestration | Hot | M5, M6, M8, M9, TemporalGrid |
| M8 | TridiagSolver | Thomas algorithm solve for TridiagSystem | Hot | — |
| M9 | AmericanConstraint | LCP enforcement (policy iteration / penalty / none) | Hot (if American) | M8, ProductSpec |
| M10 | GreeksEngine | Compute Δ,Γ,Θ,Vega,Rho; interpolation; bump-and-reprice protocols | Cold/* | M2, M4, M5 (for Θ), M12 (bumps) |
| M11 | DiagnosticsEngine | Positivity, oscillation, M-matrix checks; Richardson estimator; step/summary reports | Hot/* | M5 (coeff inputs), SpatialGrid |
| M12 | Orchestrator | Wiring + dispatch; full solve lifecycle; returns PricingResult | Cold | M1–M11 |

`*` = can be run per-step (lightweight) or post-solve / offline (heavyweight).

---

### 1.3 Data structure registry

| Struct Name | Fields (name: type, brief) | Owned By | Passed To |
|---|---|---|---|
| EngineConfig | `TimeStepperType: Enum`, `SpatialOpType: Enum`, `AmericanType: Enum`, `MeshType: Enum`, `J: int`, `N: int`, `mDomain: float`, `dampHalfSteps: int`, `mMatrixEps: float`, `minDt: float`, `maxDt: float`, `adaptTol: float`, `bumpVol: float`, `bumpRate: float`, `diagLevel: Enum` | M12 | M1,M5,M7,M9,M10,M11 |
| GridSpec | `xMin: float`, `xMax: float`, `T: float`, `J: int`, `N: int`, `meshType: Enum`, `sinhAlpha: float`, `clusterX: float`, `alignTargetsX: float[]`, `monitoringTau: float[]` | M12 | M1 |
| SpatialGrid | `J: int`, `x: float[]`, `S: float[]`, `dxL: float[]` (x_j-x_{j-1}), `dxR: float[]` (x_{j+1}-x_j), `isUniform: bool`, `h: float` (if uniform), `xMin: float`, `xMax: float` | M1 | M3,M5,M6,M7,M10,M11 |
| TemporalGrid | `N: int`, `tau: float[]`, `dt: float[]`, `monitorIdx: int[]` (step indices hitting monitoring), `isPiecewiseUniform: bool` | M1 | M7,M6 |
| VolSurface | `type: Enum`, `constSigma: float`, `SNodes: float[]`, `tauNodes: float[]`, `sigmaGrid: float[][]`, `interp: Enum` | M12 | M4,M10 |
| ModelSpec | `pdeFamily: Enum` (BS/LocalVol), `r: float`, `q: float`, `allowDegenerateA: bool` | M12 | M5,M7,M10 |
| ProductSpec | `productType: Enum`, `payoffType: Enum`, `style: Enum` (EU/AM), `S0: float`, `K: float`, `T: float`, `barrierL: float`, `barrierU: float`, `rebate: float`, `rebateStyle: Enum`, `monitoringTau: float[]`, `corridorL: float`, `corridorU: float` | M12 | M1,M3,M6,M7,M9,M10 |
| BoundarySpec | `type: Enum`, `gL: FuncPtr(tau)->float`, `gR: FuncPtr(tau)->float`, `usesBarrierAsBoundary: bool`, `barrierSide: Enum`, `rebateFn: FuncPtr(tau)->float` | M12 | M6,M7,M3 |
| BoundaryValues | `left: float`, `right: float`, `tau: float` | M6 | M7,M3 |
| OperatorCoeffs | `n: int`, `tau: float`, `ell: float[]`, `diag: float[]`, `upp: float[]`, `flags: int[]` (per-node fallback markers) | M5 | M7,M11,M10(Θ) |
| TridiagSystem | `n: int`, `lower: float[]`, `diag: float[]`, `upper: float[]`, `rhs: float[]` | M7/M9 | M8 |
| SolutionState | `tau: float`, `u: float[]`, `uPrev: float[]` (optional), `event: StepEvent`, `stepIndex: int` | M7 | M6,M7,M9,M10,M11 |
| StepEvent | `isMonitoringHit: bool`, `projectionApplied: bool`, `needsDampingRestart: bool` | M6/M7 | M7,M11 |
| AmericanSpec | `enabled: bool`, `payoffPhi: float[]`, `maxIter: int`, `tol: float`, `penaltyLambda: float` | M12/M3 | M9,M7 |
| GreekBundle | `price: float`, `delta: float`, `gamma: float`, `theta: float`, `vega: float`, `rho: float` | M10 | M12 |
| DiagnosticsReport | `minU: float`, `negCount: int`, `oscScore: float`, `oscIntervals: int[]`, `mMatrixViolCount: int`, `fallbackUpwindCount: int`, `nanCount: int`, `notes: int[]` | M11 | M12 |
| ConvergenceRecord | `pPrice: float`, `pDelta: float`, `pGamma: float`, `richExtrapPrice: float`, `errEstPrice: float`, `gridInfo: int[]` | M11 | M12 |
| PricingResult | `gb: GreekBundle`, `diag: DiagnosticsReport`, `conv: ConvergenceRecord`, `uFinal: float[]` (optional), `grid: SpatialGrid` (optional) | M12 | external |

---

## SECTION 2: MODULE SPECIFICATIONS (one subsection per module)

### 2.M1 — GridFactory

- **(a) PURPOSE —** Build spatial/temporal grids consistent with product events and boundary alignment.

- **(b) INTERFACE —**
```pseudo
GridSpec GridFactory_makeGridSpec(EngineConfig cfg, ProductSpec prod, CoordinateTransform xf);

SpatialGrid GridFactory_buildSpatial(GridSpec gs, ProductSpec prod, CoordinateTransform xf);

TemporalGrid GridFactory_buildTemporal(GridSpec gs, ProductSpec prod);

int GridFactory_alignNodeIndex(SpatialGrid xg, float xTarget); // returns nearest node index
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `prod.T > 0`, `gs.J >= 3`, `gs.N >= 1`, `gs.xMin < gs.xMax`.
  - Post (SpatialGrid): `x[j]` strictly increasing; `S[j]=xf.xToS(x[j])`; if `isUniform` then `x[j]=xMin+j*h`.
  - Post (TemporalGrid): `tau[0]=0`, `tau[N]=prod.T`; all monitoring times in `prod.monitoringTau[]` appear exactly in `tau[]` (no floating tolerance reliance).
  - Barrier-node alignment (if requested): barrier x-targets become exact nodes within `mMatrixEps` tolerance.

- **(d) STRATEGY VARIANTS —**
  - `UniformLogX` (default): uniform `x` grid per → Stage1 §3.2.1.
  - `SinhGradedX`: nonuniform `x` clustered near `clusterX` (strike/barrier); uses §3.4 below.
  - `TemporalPiecewiseUniformAligned`: per-interval uniform dt with exact monitoring hits.

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
const int MIN_J = 3;

SpatialGrid GridFactory_buildSpatial(GridSpec gs, ProductSpec prod, CoordinateTransform xf) {
  assert(gs.J >= MIN_J);

  if (gs.meshType == UNIFORM_LOG_X) {
    // Build x[j] uniformly; per Stage1 §3.2.1 domain choice comes from gs.xMin/xMax
    xg = makeUniformGrid(gs.xMin, gs.xMax, gs.J);
  } else if (gs.meshType == SINH_GRADED_X) {
    // Build sinh-graded x; formula in §3.4 (Stage2)
    xg = buildSinhGrid(gs.xMin, gs.xMax, gs.J, gs.clusterX, gs.sinhAlpha);
  }

  // Barrier/feature alignment: shift grid so each xTarget lands on a node
  for each xTarget in gs.alignTargetsX {
    int j = GridFactory_alignNodeIndex(xg, xTarget);
    float shift = xTarget - xg.x[j];
    // shift whole grid (preserve spacing pattern); if violates domain policy, caller must widen gs
    for j=0..xg.J: xg.x[j] += shift;
  }

  // Populate S and dx
  for j=0..xg.J:
    xg.S[j] = xf.xToS(xg.x[j]);
  computeDxArrays(xg); // dxL[j], dxR[j], set isUniform/h if applicable
  return xg;
}

TemporalGrid GridFactory_buildTemporal(GridSpec gs, ProductSpec prod) {
  // Build event times: {0} ∪ monitoringTau ∪ {T}; ensure sorted unique
  float[] events = sortUnique([0] + prod.monitoringTau + [prod.T]);

  // Allocate steps per interval proportional to length (at least 1), sum to gs.N
  int[] nSeg = allocateStepsPerSegment(events, gs.N);

  tg.tau = [0]; tg.dt = [];
  for s=0..len(events)-2:
    float a = events[s], b = events[s+1];
    int m = nSeg[s];
    for i=1..m:
      float tau_i = a + (b-a) * (i/m);
      tg.dt.append(tau_i - tg.tau.last());
      tg.tau.append(tau_i);
  tg.monitorIdx = indicesWhere(tg.tau matches prod.monitoringTau exactly);
  return tg;
}
```

- **(f) COMPLEXITY —**
  - Spatial build: `O(J)` time, `O(J)` space.
  - Temporal build: `O(N + F)` time, `O(N)` space (`F` monitoring count).

---

### 2.M2 — CoordinateTransform

- **(a) PURPOSE —** Provide S↔x mappings and derivative Jacobians for Greeks conversion.

- **(b) INTERFACE —**
```pseudo
float CoordinateTransform_SToX(float S, float K);  // x = ln(S/K)
float CoordinateTransform_XToS(float x, float K);  // S = K * exp(x)

void CoordinateTransform_mapUxToGreeks(
  float S0, float u, float ux, float uxx,
  out float price, out float delta, out float gamma
); // per Stage1 §3.5.1
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `S > 0`, `K > 0`.
  - Post: mapping is invertible (`XToS(SToX(S)) == S` within float tolerance).

- **(d) STRATEGY VARIANTS —**
  - `LogStrikeCentered`: `x = ln(S/K)` (required by Stage1).
  - (No other variants in Stage2; extensions can add shifted log or scaled coords.)

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
float CoordinateTransform_SToX(float S, float K) { return log(S / K); }
float CoordinateTransform_XToS(float x, float K) { return K * exp(x); }

void CoordinateTransform_mapUxToGreeks(S0, u, ux, uxx, out price, out delta, out gamma) {
  // price = u; delta/gamma per Stage1 §3.5.1
  price = u;
  delta = ...; // per Stage1 §3.5.1
  gamma = ...; // per Stage1 §3.5.1
}
```

- **(f) COMPLEXITY —** `O(1)` time/space per call.

---

### 2.M3 — PayoffProcessor

- **(a) PURPOSE —** Produce initial condition vector (cell-averaged payoff) and enforce corner compatibility.

- **(b) INTERFACE —**
```pseudo
float Payoff_eval(ProductSpec prod, float S);

float[] PayoffProcessor_cellAverageInit(
  SpatialGrid xg, ProductSpec prod, BoundarySpec bc, CoordinateTransform xf
); // per Stage1 §3.4.1

void PayoffProcessor_enforceCornerCompatibility(
  inout float[] u0, SpatialGrid xg, BoundarySpec bc
); // per Stage1 §3.4.2
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `u0.length == J+1`.
  - Post: `u0[0]=bc.gL(0)`, `u0[J]=bc.gR(0)`; interior initialized by cell averaging.
  - Payoff discontinuities are handled by integration (no point-sampling at kinks).

- **(d) STRATEGY VARIANTS —**
  - `CellAverage3PointGauss` (default): 3-point Gauss-Legendre on each cell → Stage1 §3.7 lines 30–41.
  - `PointSample` (debug only): `u0[j]=Phi(S[j])` (not recommended).

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
float[] PayoffProcessor_cellAverageInit(xg, prod, bc, xf) {
  int J = xg.J;
  float[] u0 = new float[J+1];

  u0[0] = bc.gL(0);
  u0[J] = bc.gR(0);

  for j=1..J-1:
    // Integrate Phi(K*exp(x)) over [x_j - h/2, x_j + h/2]
    // Implementation per Stage1 §3.4.1 (3-pt Gauss)
    u0[j] = cellAverageGauss3(prod, xg, j, xf); // uses Payoff_eval internally
  endfor

  PayoffProcessor_enforceCornerCompatibility(u0, xg, bc);
  return u0;
}

void PayoffProcessor_enforceCornerCompatibility(inout u0, xg, bc) {
  // Clamp boundary nodes to BC at tau=0; per Stage1 §3.4.2
  u0[0]     = bc.gL(0);
  u0[xg.J]  = bc.gR(0);
  // Optional: if boundary-adjacent cell averaging used, re-average with clamped payoff (policy flag)
}
```

- **(f) COMPLEXITY —** `O(J)` time, `O(J)` space.

---

### 2.M4 — VolSurfaceAdapter

- **(a) PURPOSE —** Provide σ(S,τ) values for operator assembly.

- **(b) INTERFACE —**
```pseudo
float VolSurfaceAdapter_sigmaAt(VolSurface vol, float S, float tau);

VolSurface VolSurfaceAdapter_bumpAll(VolSurface vol, float dSigma); // for Vega bump
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `S > 0`, `0 <= tau <= T`.
  - Post: returns finite nonnegative σ; if out-of-grid for LocalVolGrid, applies boundary extrapolation rule (clamp).

- **(d) STRATEGY VARIANTS —**
  - `ConstantVol`: `σ = vol.constSigma`.
  - `LocalVolGridBilinear`: bilinear interpolation on `(SNodes, tauNodes)` grid.

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
float VolSurfaceAdapter_sigmaAt(vol, S, tau) {
  if (vol.type == CONST_VOL) return vol.constSigma;

  if (vol.type == LOCALVOL_GRID) {
    // clamp (S,tau) into grid bounds; then bilinear interpolate
    (i0,i1,wS) = bracket(vol.SNodes, S);
    (k0,k1,wT) = bracket(vol.tauNodes, tau);
    float s00 = vol.sigmaGrid[i0][k0];
    float s10 = vol.sigmaGrid[i1][k0];
    float s01 = vol.sigmaGrid[i0][k1];
    float s11 = vol.sigmaGrid[i1][k1];
    return bilerp(s00,s10,s01,s11,wS,wT);
  }

  assert(false);
}
```

- **(f) COMPLEXITY —** `O(1)` time per lookup (grid bracketing `O(log Ns)` if binary search; `O(1)` if cached), `O(1)` space.

---

### 2.M5 — SpatialOperator

- **(a) PURPOSE —** Assemble tridiagonal operator coefficients for the spatial discretization, with M-matrix enforcement.

- **(b) INTERFACE —**
```pseudo
OperatorCoeffs SpatialOperator_assemble(
  SpatialGrid xg, float tau, ModelSpec model, VolSurface vol, EngineConfig cfg
);

bool SpatialOperator_checkMMatrixPattern(OperatorCoeffs op, EngineConfig cfg);

void SpatialOperator_applyFallbackUpwindAtNode(
  inout OperatorCoeffs op, SpatialGrid xg, int j, float tau, ModelSpec model, VolSurface vol
);
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `xg.J >= 3`, `tau` in `[0,T]`, `op.n == J-1`.
  - Post: `op.ell[j], op.diag[j], op.upp[j]` defined for all interior nodes.
  - M-matrix runtime check (pattern-level): for implicit Euler matrix `I - k L_h`, require `ell>=0` and `upp>=0` per Stage1 §3.2.4; violations trigger fallback (per-node) if enabled by cfg.

- **(d) STRATEGY VARIANTS —**
  - `FittedFD` (default): fitting factor ρ → Stage1 §3.2.3; coefficients → Stage1 §3.2.4.
  - `CentralFD`: ρ ≡ 1 (comparison / debug).
  - `PureUpwind`: replaces convection discretization with upwind; used as fallback when M-matrix pattern fails.

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
OperatorCoeffs SpatialOperator_assemble(xg, tau, model, vol, cfg) {
  int J = xg.J;
  OperatorCoeffs op = allocOp(n = J-1, tau=tau);

  for j=1..J-1:
    float S = xg.S[j];
    float sig = VolSurfaceAdapter_sigmaAt(vol, S, tau);

    // Compute a(x,tau), b(x,tau) in log-space; per Stage1 §3.1.4
    // Compute fitting factor rho; per Stage1 §3.2.3
    // Compute ell/diag/upp for L_h; per Stage1 §3.2.4
    (op.ell[j], op.diag[j], op.upp[j]) = computeFittedStencil(Stage1_refs, xg, j, sig, model);

  endfor

  if (!SpatialOperator_checkMMatrixPattern(op, cfg) && cfg.SpatialOpType == FITTED_FD) {
    // Per-node fallback: only fix nodes that violate; keep rest fitted
    for j=1..J-1:
      if (op.ell[j] < -cfg.mMatrixEps || op.upp[j] < -cfg.mMatrixEps) {
        SpatialOperator_applyFallbackUpwindAtNode(op, xg, j, tau, model, vol);
        op.flags[j] |= FLAG_UPWIND_FALLBACK;
      }
    endfor
  }

  return op;
}

bool SpatialOperator_checkMMatrixPattern(op, cfg) {
  int viol = 0;
  for j=1..op.n:
    if (op.ell[j] < -cfg.mMatrixEps) viol++;
    if (op.upp[j] < -cfg.mMatrixEps) viol++;
  return (viol == 0);
}
```

- **(f) COMPLEXITY —**
  - Assembly: `O(J)` time, `O(J)` space for coefficient arrays.
  - Fallback scan: `O(J)`.

---

### 2.M6 — BoundaryHandler

- **(a) PURPOSE —** Evaluate boundary values and apply discrete monitoring projection, emitting discontinuity events.

- **(b) INTERFACE —**
```pseudo
BoundaryValues BoundaryHandler_eval(BoundarySpec bc, float tau);

StepEvent BoundaryHandler_applyMonitoringProjection(
  inout SolutionState st, SpatialGrid xg, ProductSpec prod, float tau
); // projection per Stage1 §3.1.2

void BoundaryHandler_applyDirichlet(inout SolutionState st, BoundaryValues bv);
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `st.u.length == J+1`.
  - Post: after `applyDirichlet`, boundary nodes match bc at `tau`.
  - Post: if monitoring projection applied, `st.u[j]=0` for `S[j] خارج corridor`; event flags set so TimeStepper restarts damping.

- **(d) STRATEGY VARIANTS —**
  - `DirichletVanilla`: vanilla call/put Dirichlet boundaries per Stage1 §3.1.2.
  - `BarrierKnockOut`: barrier boundary is Dirichlet `0` (or corridor truncation at boundary).
  - `BarrierRebate`: barrier boundary is Dirichlet `rebateFn(tau)` (rebate style encoded in bc).

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
BoundaryValues BoundaryHandler_eval(bc, tau) {
  BoundaryValues bv;
  bv.left  = bc.gL(tau);
  bv.right = bc.gR(tau);
  bv.tau = tau;
  return bv;
}

void BoundaryHandler_applyDirichlet(inout st, bv) {
  st.u[0]      = bv.left;
  st.u[st.u.length-1] = bv.right;
}

StepEvent BoundaryHandler_applyMonitoringProjection(inout st, xg, prod, tau) {
  StepEvent ev = {isMonitoringHit=false, projectionApplied=false, needsDampingRestart=false};

  if (!isMonitoringTime(prod.monitoringTau, tau)) return ev;

  ev.isMonitoringHit = true;

  // Corridor projection per Stage1 §3.1.2 (and Algorithm 1 lines 183–193)
  for j=0..xg.J:
    if (xg.S[j] < prod.corridorL || xg.S[j] > prod.corridorU) st.u[j] = 0.0;
  endfor

  ev.projectionApplied = true;
  ev.needsDampingRestart = true;
  return ev;
}
```

- **(f) COMPLEXITY —**
  - Boundary eval/apply: `O(1)`.
  - Projection: `O(J)` per monitoring event.

---

### 2.M7 — TimeStepper

- **(a) PURPOSE —** Advance the solution in time using selected time integration and damping policy.

- **(b) INTERFACE —**
```pseudo
enum TimeStepperState { DAMPING, NORMAL, POST_MONITOR_DAMPING };

struct TimeStepperCtx {
  TimeStepperState state;
  int remainingHalfSteps; // for Rannacher
  float gammaTRBDF2;      // for TR-BDF2
};

SolutionState TimeStepper_advance(
  inout TimeStepperCtx ctx,
  SolutionState st,
  SpatialGrid xg,
  TemporalGrid tg,
  int n, // step index (τ_n -> τ_{n+1})
  ProductSpec prod,
  ModelSpec model,
  VolSurface vol,
  BoundarySpec bc,
  AmericanSpec am,
  EngineConfig cfg,
  TridiagSolver solver
);
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `st.tau == tg.tau[n]`.
  - Post: returns `stNext.tau == tg.tau[n+1]`, with boundaries enforced at `τ_{n+1}`.
  - Damping contract (RS-CN): after any `StepEvent.needsDampingRestart`, perform `cfg.dampHalfSteps` implicit Euler half-steps (each `dt/2`) before resuming normal stepping → Stage1 §3.3.1.
  - TR-BDF2 contract: uses `γ = 2 - sqrt(2)` for L-stability; two solves per step → §3.1 (Stage2).

- **(d) STRATEGY VARIANTS —**
  - `RannacherCN` (default): Euler half-steps then CN → Stage1 §3.3.1.
  - `TRBDF2`: two-stage L-stable, γ=2-√2 → §3.1.
  - `PureImplicitEuler`: always implicit Euler (debug/robust).

- **(e) KEY ALGORITHM SKELETON —**
```text
State machine (ctx.state) for RS-CN:

┌──────────────┐   initial (n=0) or restart event   ┌─────────────────────┐
│   NORMAL     │◀───────────────────────────────────│ POST_MONITOR_DAMPING│
└──────┬───────┘                                    └─────────┬───────────┘
       │  after monitor projection event                        │ set remainingHalfSteps = cfg.dampHalfSteps
       │                                                        ▼
       │                                        halfsteps done? yes ─────▶ NORMAL
       ▼
┌──────────────┐
│   DAMPING    │  (used at start; same mechanics as POST_MONITOR_DAMPING)
└──────┬───────┘
       │ halfsteps done? yes
       └────────────────────────────────────────────▶ NORMAL
```

```pseudo
SolutionState TimeStepper_advance(inout ctx, st, xg, tg, n, prod, model, vol, bc, am, cfg, solver) {
  float tau_n   = tg.tau[n];
  float tau_np1 = tg.tau[n+1];
  float dt      = tg.dt[n];

  // Boundary at current level (defensive)
  BoundaryHandler_applyDirichlet(st, BoundaryHandler_eval(bc, tau_n));

  // Monitoring projection happens at τ_{n+1} boundary in this architecture
  // (caller may apply right after step; if applied, it sets ctx.state via event)

  if (cfg.TimeStepperType == RANNACHER_CN) {
    if (n == 0 && ctx.state != DAMPING) { ctx.state = DAMPING; ctx.remainingHalfSteps = cfg.dampHalfSteps; }

    if (ctx.state == DAMPING || ctx.state == POST_MONITOR_DAMPING) {
      // Run implicit Euler half-steps until we reach τ_{n+1}
      float half = 0.5 * dt;
      for i=1..2: // exactly two half-steps per full step (dt split)
        float tau_stage = tau_n + i * half;
        OperatorCoeffs op = SpatialOperator_assemble(xg, tau_stage, model, vol, cfg);
        TridiagSystem sys = buildImplicitEulerSystem(op, xg, st.u, bc, tau_stage, half); // → Stage1 §3.3.1
        float[] u_new = AmericanConstraint_solveOrNone(am, sys, solver);                 // → M9
        st.u = writeBackInteriorAndBC(u_new, st.u, bc, tau_stage);
      endfor

      ctx.remainingHalfSteps -= 2;
      if (ctx.remainingHalfSteps <= 0) ctx.state = NORMAL;

      st.tau = tau_np1;
      return st;
    }

    // NORMAL: one CN step
    OperatorCoeffs op0 = SpatialOperator_assemble(xg, tau_n,   model, vol, cfg);
    OperatorCoeffs op1 = SpatialOperator_assemble(xg, tau_np1, model, vol, cfg);
    TridiagSystem sysCN = buildCrankNicolsonSystem(op0, op1, xg, st.u, bc, tau_n, tau_np1, dt); // → Stage1 §3.3.1
    float[] u_next = AmericanConstraint_solveOrNone(am, sysCN, solver);
    st.u = writeBackInteriorAndBC(u_next, st.u, bc, tau_np1);
    st.tau = tau_np1;
    return st;
  }

  if (cfg.TimeStepperType == TR_BDF2) {
    st = TRBDF2_step(st, xg, bc, model, vol, am, cfg, dt, tau_n, tau_np1, solver); // → §3.1
    return st;
  }

  // Pure implicit Euler
  OperatorCoeffs opE = SpatialOperator_assemble(xg, tau_np1, model, vol, cfg);
  TridiagSystem sysE = buildImplicitEulerSystem(opE, xg, st.u, bc, tau_np1, dt);
  float[] uE = AmericanConstraint_solveOrNone(am, sysE, solver);
  st.u = writeBackInteriorAndBC(uE, st.u, bc, tau_np1);
  st.tau = tau_np1;
  return st;
}
```

- **(f) COMPLEXITY —**
  - RS-CN: `O(J)` per (half)step (assembly + Thomas), total `O(NJ)`; extra `O(J)` per monitoring restart.
  - TR-BDF2: `2 * O(J)` per step (two solves).

---

### 2.M8 — TridiagSolver

- **(a) PURPOSE —** Solve tridiagonal linear systems via Thomas algorithm.

- **(b) INTERFACE —**
```pseudo
float[] TridiagSolver_solve(TridiagSystem sys);

bool TridiagSolver_checkDiagDominance(TridiagSystem sys); // optional debug
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `sys.n >= 1`, `diag[i] != 0` for all i.
  - Post: returns `x[]` satisfying `A x = rhs` within solver tolerance.
  - In-place modification of `sys.diag`/`sys.rhs` allowed (caller must pass copy if needed).

- **(d) STRATEGY VARIANTS —**
  - `ThomasNoPivot` (default): per Stage1 §3.6.1.
  - `ThomasWithChecks` (debug): asserts no near-zero pivots, optional diag dominance test.

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
float[] TridiagSolver_solve(sys) {
  int n = sys.n;
  // forward elimination
  for i=1..n-1:
    float m = sys.lower[i] / sys.diag[i-1];
    sys.diag[i] = sys.diag[i] - m * sys.upper[i-1];
    sys.rhs[i]  = sys.rhs[i]  - m * sys.rhs[i-1];
  endfor
  // back substitution
  float[] x = new float[n];
  x[n-1] = sys.rhs[n-1] / sys.diag[n-1];
  for i=n-2..0:
    x[i] = (sys.rhs[i] - sys.upper[i] * x[i+1]) / sys.diag[i];
  endfor
  return x;
}
```

- **(f) COMPLEXITY —** `O(n)` time, `O(n)` space for solution.

---

### 2.M9 — AmericanConstraint

- **(a) PURPOSE —** Enforce early-exercise constraint (LCP) per time step.

- **(b) INTERFACE —**
```pseudo
float[] AmericanConstraint_solveOrNone(AmericanSpec am, TridiagSystem sys, TridiagSolver solver);

float[] AmericanConstraint_policyIteration(AmericanSpec am, TridiagSystem sys, TridiagSolver solver); // → §3.2

float[] AmericanConstraint_penaltyMethod(AmericanSpec am, TridiagSystem sys, TridiagSolver solver);   // → §3.3
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `am.enabled` implies `am.payoffPhi.length == sys.n` (interior nodes).
  - Post (American): `u[i] >= payoffPhi[i] - am.tol` (componentwise).
  - Post (European / None): returns linear solve result.

- **(d) STRATEGY VARIANTS —**
  - `None`: no-op; return `solver.solve(sys)`.
  - `PolicyIteration` (default): active-set iterations → §3.2.
  - `PenaltyMethod`: penalized nonlinear iteration → §3.3.

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
float[] AmericanConstraint_solveOrNone(am, sys, solver) {
  if (!am.enabled) return TridiagSolver_solve(sys);

  if (am.type == POLICY_ITERATION) return AmericanConstraint_policyIteration(am, sys, solver);
  if (am.type == PENALTY_METHOD)   return AmericanConstraint_penaltyMethod(am, sys, solver);

  assert(false);
}
```

- **(f) COMPLEXITY —**
  - None: `O(J)`.
  - Policy iteration: `O(I * J)` where `I` is iteration count.
  - Penalty: `O(I * J)`.

---

### 2.M10 — GreeksEngine

- **(a) PURPOSE —** Compute price and Greeks from final solution with robust interpolation and bump protocols.

- **(b) INTERFACE —**
```pseudo
GreekBundle GreeksEngine_computeAll(
  PricingResult base,
  EngineConfig cfg,
  ProductSpec prod,
  ModelSpec model,
  VolSurface vol,
  CoordinateTransform xf,
  Orchestrator engine // used only for bump-and-reprice
);

void GreeksEngine_quadInterpUxUxx(
  SpatialGrid xg, float[] u, float x0,
  out float u0, out float ux0, out float uxx0
); // per Stage1 §3.5.2
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: base solution corresponds to τ=T (“today” in Stage1 convention).
  - Post: `gb.price` equals interpolated `u(x0,T)`; Δ/Γ per Stage1 §3.5.1; Θ computed by PDE residual at τ=T (see below); Vega/Rho by central bumps.

- **(d) STRATEGY VARIANTS —**
  - `GreeksFromInterpLogDerivs` (default): Δ/Γ via `(u,ux,uxx)` in log-space → Stage1 §3.5.1–§3.5.2.
  - `BumpAll` (optional): bump-and-reprice Δ/Γ (slow; debugging only).

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
GreekBundle GreeksEngine_computeAll(base, cfg, prod, model, vol, xf, engine) {
  float x0 = xf.SToX(prod.S0, prod.K);

  // Interpolate u, ux, uxx at x0; per Stage1 §3.5.2
  (u0, ux0, uxx0) = GreeksEngine_quadInterpUxUxx(base.grid, base.uFinal, x0);

  GreekBundle gb;
  xf.mapUxToGreeks(prod.S0, u0, ux0, uxx0, out gb.price, out gb.delta, out gb.gamma); // Stage1 §3.5.1

  // Theta via PDE residual at (x0, tau=T): Theta_calendar = -u_tau ≈ -L(u)
  // Use operator apply with coefficients at tau=T (assembled once); no re-derivation here.
  gb.theta = computeThetaFromSpatialOperatorResidual(base, cfg, prod, model, vol, x0); // uses M5 + derivatives

  // Vega bump-and-reprice (central)
  float epsV = cfg.bumpVol;
  PricingResult upV = engine.priceOnly(cfg, prod, model, VolSurfaceAdapter_bumpAll(vol, +epsV));
  PricingResult dnV = engine.priceOnly(cfg, prod, model, VolSurfaceAdapter_bumpAll(vol, -epsV));
  gb.vega = (upV.gb.price - dnV.gb.price) / (2*epsV);

  // Rho bump-and-reprice (central)
  float epsR = cfg.bumpRate;
  PricingResult upR = engine.priceOnly(cfg, prod, modelWithRate(model, model.r + epsR), vol);
  PricingResult dnR = engine.priceOnly(cfg, prod, modelWithRate(model, model.r - epsR), vol);
  gb.rho = (upR.gb.price - dnR.gb.price) / (2*epsR);

  return gb;
}
```

- **(f) COMPLEXITY —**
  - Δ/Γ/Θ from one run: `O(1)` extra beyond base solve (plus operator assemble if Θ uses it).
  - Vega/Rho: +2 solves each → total `O(#bumps * N * J)`.

---

### 2.M11 — DiagnosticsEngine

- **(a) PURPOSE —** Provide runtime health checks and convergence estimates.

- **(b) INTERFACE —**
```pseudo
DiagnosticsReport DiagnosticsEngine_evalStep(
  SolutionState st, SpatialGrid xg, OperatorCoeffs op, EngineConfig cfg
);

float DiagnosticsEngine_minU(float[] u);
int   DiagnosticsEngine_countNeg(float[] u, float eps);

OscResult DiagnosticsEngine_detectOscillations(float[] u, SpatialGrid xg, EngineConfig cfg); // → §3.7

ConvergenceRecord DiagnosticsEngine_richardson(
  PricingResult coarse, PricingResult fine
); // → §3.8
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `u.length == J+1`.
  - Post: report fields are finite; if NaNs detected, increments `nanCount`.

- **(d) STRATEGY VARIANTS —**
  - `LightStepChecks`: min/neg/NaN + cheap osc score.
  - `FullStepChecks`: also M-matrix violations count and fallback counts.
  - `OfflineConvergence`: Richardson estimator on multi-run results.

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
DiagnosticsReport DiagnosticsEngine_evalStep(st, xg, op, cfg) {
  DiagnosticsReport dr;
  dr.minU = DiagnosticsEngine_minU(st.u);
  dr.negCount = DiagnosticsEngine_countNeg(st.u, eps = 0.0);
  dr.nanCount = countNaNs(st.u);

  if (cfg.diagLevel >= DIAG_OSC) {
    OscResult osc = DiagnosticsEngine_detectOscillations(st.u, xg, cfg); // §3.7
    dr.oscScore = osc.score;
    dr.oscIntervals = osc.flaggedIntervals;
  }

  if (cfg.diagLevel >= DIAG_M_MATRIX) {
    dr.mMatrixViolCount = countViolations(op.ell, op.upp, cfg.mMatrixEps);
    dr.fallbackUpwindCount = countFlags(op.flags, FLAG_UPWIND_FALLBACK);
  }

  return dr;
}
```

- **(f) COMPLEXITY —**
  - Step checks: `O(J)` time, `O(1)+O(#flags)` space.
  - Richardson: `O(1)` time once inputs exist (excludes cost of producing inputs).

---

### 2.M12 — Orchestrator (top-level driver)

- **(a) PURPOSE —** Wire modules, execute solve lifecycle, and return PricingResult for all supported products.

- **(b) INTERFACE —**
```pseudo
PricingResult Orchestrator_price(EngineConfig cfg, MarketData md, ProductSpec prod, ModelSpec model, VolSurface vol);

PricingResult Orchestrator_priceOnly(EngineConfig cfg, ProductSpec prod, ModelSpec model, VolSurface vol); // used by bumps
```

- **(c) INVARIANTS / CONTRACTS —**
  - Pre: `prod.S0>0`, `prod.K>0`, `prod.T>0`; corridor/barriers consistent.
  - Post: `PricingResult.gb.price` is the engine output at τ=T; diagnostics included as requested by cfg.

- **(d) STRATEGY VARIANTS —**
  - Dispatch occurs via `cfg` + `prod`:
    - SpatialOperator: `FittedFD/CentralFD/PureUpwind`
    - TimeStepper: `RannacherCN/TRBDF2/Euler`
    - AmericanConstraint: `None/PolicyIteration/PenaltyMethod`
    - BoundaryHandler: `Dirichlet/BarrierKnockOut/BarrierRebate`

- **(e) KEY ALGORITHM SKELETON —**
```pseudo
PricingResult Orchestrator_price(cfg, md, prod, model, vol) {
  // 1) Build transform + grid spec
  CoordinateTransform xf = makeLogTransform(K = prod.K);
  GridSpec gs = GridFactory_makeGridSpec(cfg, prod, xf);            // domain, alignment targets, monitoring
  SpatialGrid xg = GridFactory_buildSpatial(gs, prod, xf);
  TemporalGrid tg = GridFactory_buildTemporal(gs, prod);

  // 2) Boundary spec (Dirichlet functions)
  BoundarySpec bc = buildBoundarySpec(cfg, md, prod, xf);           // per Stage1 §3.1.2

  // 3) Payoff init (cell-averaged) + American payoff vector (interior)
  float[] u0 = PayoffProcessor_cellAverageInit(xg, prod, bc, xf);   // Stage1 §3.4.1/§3.4.2
  AmericanSpec am = buildAmericanSpec(cfg, prod, u0, xg);           // payoffPhi interior extracted

  // 4) Initialize state + timestep context
  SolutionState st = {tau=0.0, u=u0, stepIndex=0, event={...}};
  TimeStepperCtx ts = initTimeStepperCtx(cfg); // state=DAMPING for RS-CN, gamma for TR-BDF2

  DiagnosticsReport drAll = emptyDiagnostics();

  // 5) March forward in τ from 0 to T
  for n=0..tg.N-1:
    st.stepIndex = n;

    // Step: τ_n → τ_{n+1}
    st = TimeStepper_advance(ts, st, xg, tg, n, prod, model, vol, bc, am, cfg, solver=M8);

    // Apply monitoring projection exactly at τ_{n+1}; restart damping if needed
    StepEvent ev = BoundaryHandler_applyMonitoringProjection(st, xg, prod, tg.tau[n+1]);
    if (ev.needsDampingRestart) {
      ts.state = POST_MONITOR_DAMPING;
      ts.remainingHalfSteps = cfg.dampHalfSteps;
    }

    // Diagnostics (optional per cfg)
    if (cfg.diagLevel != DIAG_OFF) {
      OperatorCoeffs op = SpatialOperator_assemble(xg, st.tau, model, vol, cfg); // reuse allowed if cached
      DiagnosticsReport dr = DiagnosticsEngine_evalStep(st, xg, op, cfg);
      drAll = mergeDiagnostics(drAll, dr);
    }
  endfor

  // 6) Final price + Greeks
  PricingResult base;
  base.uFinal = st.u;
  base.grid = xg;
  base.diag = drAll;

  base.gb = GreeksEngine_computeAll(base, cfg, prod, model, vol, xf, engine=this);

  return base;
}
```

- **(f) COMPLEXITY —**
  - Base solve: `O(NJ)` time, `O(J)` memory.
  - With American: multiply by constraint iteration count.
  - With Vega/Rho bumps: multiply by number of bump solves.

---

## SECTION 3: CRITICAL ALGORITHMS NOT IN STAGE 1

### 3.1 TR-BDF2 Time Step (γ = 2 − √2)

```pseudo
SolutionState TRBDF2_step(
  SolutionState st_n,
  SpatialGrid xg,
  BoundarySpec bc,
  ModelSpec model,
  VolSurface vol,
  AmericanSpec am,
  EngineConfig cfg,
  float dt, float tau_n, float tau_np1,
  TridiagSolver solver
) {
  const float gamma = 2.0 - sqrt(2.0);     // per Bank et al. (1996); Bonaventura–Della Rocca (2015)
  float tau_g = tau_n + gamma * dt;

  // Stage 1 (TR on substep gamma*dt): CN-form on [tau_n, tau_g]
  OperatorCoeffs L_n = SpatialOperator_assemble(xg, tau_n, model, vol, cfg);
  OperatorCoeffs L_g = SpatialOperator_assemble(xg, tau_g, model, vol, cfg);

  TridiagSystem A1 = build_CN_LHS(L_g, dt1 = gamma*dt);             // (I - dt1/2 L_g)
  float[] rhs1 = build_CN_RHS(L_n, st_n.u, dt1 = gamma*dt);         // (I + dt1/2 L_n) u_n
  applyDirichletToSystem(A1, rhs1, bc, tau_g, xg);

  float[] u_g_int = AmericanConstraint_solveOrNone(am, withRHS(A1,rhs1), solver);
  float[] u_g = writeBackInteriorAndBC(u_g_int, st_n.u, bc, tau_g);

  // Stage 2 (variable-step BDF2 over remaining (1-gamma)dt using {u_n, u_g}):
  // (I - alpha*dt*L_{n+1}) u_{n+1} = c1*u_g + c0*u_n
  // alpha, c0, c1 are TR-BDF2 constants for this gamma (no derivation here)
  float alpha = (1.0 - gamma) / (2.0 - gamma);
  float c1    = 1.0 / (gamma * (2.0 - gamma));
  float c0    = -((1.0 - gamma)*(1.0 - gamma)) / (gamma * (2.0 - gamma));

  OperatorCoeffs L_np1 = SpatialOperator_assemble(xg, tau_np1, model, vol, cfg);

  TridiagSystem A2 = build_Implicit_LHS(L_np1, dt2 = alpha*dt);     // (I - alpha*dt*L_{n+1})
  float[] rhs2 = linearCombo(c1, u_g, c0, st_n.u);                  // interior combo
  applyDirichletToSystem(A2, rhs2, bc, tau_np1, xg);

  float[] u_np1_int = AmericanConstraint_solveOrNone(am, withRHS(A2,rhs2), solver);
  float[] u_np1 = writeBackInteriorAndBC(u_np1_int, st_n.u, bc, tau_np1);

  return {tau=tau_np1, u=u_np1, stepIndex=st_n.stepIndex+1, event=st_n.event};
}
```

---

### 3.2 Policy Iteration for American LCP (active-set)

```pseudo
float[] PolicyIteration_LCP(
  TridiagSystem A, float[] rhs,
  float[] phi, AmericanSpec am, TridiagSolver solver
) {
  int n = A.n;
  float[] u = maxVec(phi, initialGuessFromLinearSolve(A, rhs, solver)); // safe start
  bool[] exercise = new bool[n];

  for iter=0..am.maxIter-1:
    // Update active set (policy): enforce whichever constraint is tighter
    // LCP form: max(phi - u, A*u - rhs) = 0  (elementwise)
    float[] Au = tridiagMatVec(A, u);
    int changes = 0;
    for i=0..n-1:
      bool exNew = (phi[i] - u[i]) >= (Au[i] - rhs[i]);
      if (exNew != exercise[i]) changes++;
      exercise[i] = exNew;
    endfor
    if (changes == 0) break;

    // Build modified system: exercise rows -> identity; continuation rows -> original A
    TridiagSystem M = copyTridiag(A);
    float[] b = rhs.copy();
    for i=0..n-1:
      if (exercise[i]) {
        zeroRowTridiag(M, i);
        M.diag[i] = 1.0;
        b[i] = phi[i];
      }
    endfor

    float[] uNew = TridiagSolver_solve(withRHS(M,b));
    if (maxAbs(uNew - u) < am.tol) { u = uNew; break; }
    u = uNew;
  endfor

  // Ensure final projection (tolerance)
  for i=0..n-1: u[i] = max(u[i], phi[i]);
  return u;
}
```

---

### 3.3 Penalty Method for American LCP (semi-smooth Newton style)

```pseudo
float[] PenaltyMethod_LCP(
  TridiagSystem A, float[] rhs,
  float[] phi, AmericanSpec am, TridiagSolver solver
) {
  int n = A.n;
  float lambda = am.penaltyLambda;
  float[] u = TridiagSolver_solve(withRHS(copyTridiag(A), rhs.copy())); // start from European

  for iter=0..am.maxIter-1:
    // D_i = 1 if u_i < phi_i else 0
    float[] D = new float[n];
    for i=0..n-1: D[i] = (u[i] < phi[i]) ? 1.0 : 0.0;

    // Solve (A + lambda*diag(D)) uNew = rhs + lambda*diag(D)*phi
    TridiagSystem M = copyTridiag(A);
    float[] b = rhs.copy();
    for i=0..n-1:
      M.diag[i] += lambda * D[i];
      b[i]      += lambda * D[i] * phi[i];
    endfor

    float[] uNew = TridiagSolver_solve(withRHS(M,b));
    if (maxAbs(uNew - u) < am.tol) { u = uNew; break; }
    u = uNew;
  endfor

  for i=0..n-1: u[i] = max(u[i], phi[i]);
  return u;
}
```

---

### 3.4 Sinh-Graded Mesh Construction (x-grid)

```pseudo
SpatialGrid buildSinhGrid(float xMin, float xMax, int J, float xCenter, float alpha) {
  // ξ_j uniform in [0,1]; map to [-1,1]; sinh cluster around ξ=0.5
  float[] x = new float[J+1];

  // Choose scale c so endpoints match exactly:
  // xMin = xCenter + c*sinh(alpha*(ξ0-0.5)), xMax = xCenter + c*sinh(alpha*(ξJ-0.5))
  // With ξ0=0, ξJ=1 => symmetric: sinh(-alpha/2), sinh(+alpha/2)
  float denom = sinh(alpha * 0.5);
  float c = (xMax - xMin) / (2.0 * denom);

  for j=0..J:
    float xi = (float)j / (float)J;              // [0,1]
    float z  = alpha * (xi - 0.5);               // [-alpha/2, +alpha/2]
    x[j] = xCenter + c * sinh(z);
  endfor

  // Shift to ensure exact endpoints (numeric safety)
  float shift = xMin - x[0];
  for j=0..J: x[j] += shift;

  SpatialGrid xg = {J=J, x=x, xMin=x[0], xMax=x[J], isUniform=false};
  computeDxArrays(xg);
  return xg;
}
```

---

### 3.5 Keller Box Scheme (simultaneous price + delta)

```pseudo
// Variables: u(x,τ), p(x,τ)=u_x. PDE in log-space: -u_τ + a p_x + b p - r u = 0; and p - u_x = 0.
// Discretize on boxes (j-1/2,j+1/2) with staggered averages; Duffy reference (box scheme).
BlockTriSystem KellerBox_assembleStep(
  SpatialGrid xg, float tau_n, float tau_np1, float dt,
  float[] u_n, float[] p_n, ModelSpec model, VolSurface vol, EngineConfig cfg, BoundarySpec bc
) {
  // Assemble block-tridiagonal for unknowns (u_{j}^{n+1}, p_{j}^{n+1}) over interior nodes
  // Each row is 2x2 block; uses midpoint/box averages (no re-derivation here).
  BlockTriSystem B = allocBlockTri(nBlocks = xg.J-1, blockSize=2);

  for j=1..xg.J-1:
    // compute a,b at box midpoints; sigmaAt via M4
    // fill blocks: B.lower[j], B.diag[j], B.upper[j] (each 2x2), and rhs block (2)
    fillKellerBoxBlocks(B, j, Stage1_PDE_coeff_refs, xg, tau_n, tau_np1, dt, u_n, p_n, model, vol, bc);
  endfor
  applyBlockDirichlet(B, bc, tau_np1, xg);
  return B;
}

(u_np1, p_np1) KellerBox_solve(BlockTriSystem B) {
  // Block Thomas for 2x2 blocks
  return blockThomasSolve2x2(B);
}
```

---

### 3.6 Adaptive Time-Step Controller

```pseudo
float AdaptiveController_nextDt(
  float dt, float errEst, float tol, int order,
  bool inDampingPhase,
  float nextEventTau, float currentTau,
  EngineConfig cfg
) {
  if (inDampingPhase) return dt; // no adaptation during Rannacher/forced damping

  // PI-like controller (simplified): dt_new = dt * safety * (tol/err)^(1/(order+1))
  const float safety = 0.9;
  const float growMax = 2.0;
  const float shrinkMin = 0.5;

  float scale = (errEst <= 0) ? growMax : safety * pow(tol / errEst, 1.0 / (order + 1));
  scale = clamp(scale, shrinkMin, growMax);

  float dtNew = clamp(dt * scale, cfg.minDt, cfg.maxDt);

  // Event alignment: do not step past next monitoring time
  float remaining = nextEventTau - currentTau;
  if (dtNew > remaining) dtNew = remaining;

  return dtNew;
}
```

---

### 3.7 Oscillation Detector (sign changes of Δu)

```pseudo
struct OscResult { float score; int[] flaggedIntervals; };

OscResult detectOscillations(float[] u, SpatialGrid xg, EngineConfig cfg) {
  // Compute Δu_j = u_{j+1} - u_j on interior; count sign changes in Δu
  int J = xg.J;
  int changes = 0;
  int[] flags = [];

  float prev = u[1] - u[0];
  int prevSign = sign(prev, eps=0.0);

  int windowStart = 0;
  int windowChanges = 0;

  for j=1..J-1:
    float du = u[j+1] - u[j];
    int s = sign(du, eps=0.0);
    if (s != 0 && prevSign != 0 && s != prevSign) {
      changes++; windowChanges++;
    }
    if (s != 0) prevSign = s;

    // Window flagging (fixed node window; map to x-interval)
    if ((j - windowStart) >= cfg.oscWindowNodes) {
      if (windowChanges >= cfg.oscChangeThreshold) flags.append(windowStart);
      windowStart = j; windowChanges = 0;
    }
  endfor

  OscResult r;
  r.score = (float)changes / max(1, J-1);   // normalized severity
  r.flaggedIntervals = flags;
  return r;
}
```

---

### 3.8 Richardson Convergence Estimator

```pseudo
ConvergenceRecord RichardsonEstimate(PricingResult coarse, PricingResult fine) {
  // Assumes fine grid is (2J,2N) of coarse; uses scalar comparisons at S0
  float Pc = coarse.gb.price, Pf = fine.gb.price;
  float Dc = coarse.gb.delta, Df = fine.gb.delta;
  float Gc = coarse.gb.gamma, Gf = fine.gb.gamma;

  // If a third run (4J,4N) exists, prefer p from three-level ratios.
  // With two levels only, report p as "unknown" and provide err estimate ~ |Pf-Pc|.
  ConvergenceRecord cr;
  cr.pPrice = NAN; cr.pDelta = NAN; cr.pGamma = NAN;

  cr.errEstPrice = abs(Pf - Pc);
  cr.richExtrapPrice = Pf; // without p, no extrapolation; caller can upgrade with 3 levels

  cr.gridInfo = [coarse.grid.J, fine.grid.J, coarseTemporalN, fineTemporalN];
  return cr;
}
```

---

## SECTION 4: PRODUCT DISPATCH & CONFIGURATION

### 4.1 Product × configuration matrix

| Product | Spatial | TimeStepper | American | Boundary | PayoffSmooth |
|---|---|---|---|---|---|
| EU Call/Put | `FittedFD` | `RannacherCN` | `None` | `DirichletVanilla` | `CellAverage3PointGauss` |
| EU Digital | `FittedFD` + (fallback upwind on M-matrix viol) | `TRBDF2` (or RS-CN with more damping) | `None` | `DirichletVanilla` | `CellAverage3PointGauss` |
| AM Call/Put | `FittedFD` | `PureImplicitEuler` (or TRBDF2 if stable) | `PolicyIteration` | `DirichletVanilla` | `CellAverage3PointGauss` |
| Down-Out Barrier (cont. BC) | `FittedFD` | `RannacherCN` | `None` | `BarrierKnockOut` / `BarrierRebate` | `CellAverage3PointGauss` |
| Disc-Mon DKO (corridor projection) | `FittedFD` | `RannacherCN` with restart | `None` | `BarrierKnockOut` + projection | `CellAverage3PointGauss` |

---

### 4.2 Default numerical parameter recommendations

| Product | J (space) | N (time) | dampHalfSteps | mDomain (stddev) | Notes |
|---|---:|---:|---:|---:|---|
| EU Call/Put | 600 | 600 | 2 | 8 | baseline |
| EU Digital | 1200 | 1200 | 4 | 10 | more damping + finer grid |
| AM Put | 800 | 800 | 0–2 | 8 | prefer implicit Euler for LCP robustness |
| Down-Out Barrier | 800 | 1000 | 2 | 8–10 | align barrier node; widen domain if needed |
| Disc-Mon DKO | 1200 | 1500 | 2 (restart each monitor) | 10 | ensure monitoring-aligned TemporalGrid |

(Use `MeshType=SINH_GRADED_X` with `alpha≈2–4` when clustering is needed near strike/barriers.)

---

## SECTION 5: TESTING & VALIDATION FRAMEWORK

### 5.1 Benchmark test cases (5–7)

| Test ID | Params (S0,K,T,r,q,σ,barrier/mon) | Reference | Expected order (Price/Δ/Γ) | Failure modes hit (Stage1 §1.3.1) |
|---|---|---|---|---|
| T1 | EU Call: (100,100,1,0.05,0,0.2) | BS closed-form | ~2 / ~2 / ~2 (with damping) | F2,F8,F6 |
| T2 | EU Digital Call (cash-or-nothing), same market | analytic digital (or high-precision quad) | price ~2 after smoothing; Greeks slower | F2,F4,F8 |
| T3 | Low-vol convection: EU Call (100,100,0.5,0.05,0,0.02) | BS closed-form | price stable; avoid oscillations | F1,F4,F5,F8 |
| T4 | Down-and-out Call (continuous barrier), rebate=0 | analytic barrier (or QuantLib high-prec) | price convergent; Δ/Γ near barrier sensitive | F6,F8 |
| T5 | Discrete monitored double KO (corridor projection), 5 monitors | MC high-precision + fine PDE | qualitative: no negativity; convergence slower | F3,F2,F4,F8 |
| T6 | American Put (100,100,1,0.05,0,0.2) | binomial/PSOR high-prec | price monotone; Δ approx | F7,F4,F8 |
| T7 | Boundary sensitivity: EU Call, vary xMax/xMin | invariance check | price change < tol | F6 |

Mapping ensures coverage of **F1–F8** across the suite.

---

### 5.2 Convergence test procedure flowchart

```text
┌───────────────────────────────┐
│ Choose test (T1..T7)          │
└───────────────┬───────────────┘
                ▼
┌───────────────────────────────┐
│ Run solver at (J,N)           │
│  → PricingResult R1           │
└───────────────┬───────────────┘
                ▼
┌───────────────────────────────┐
│ Run solver at (2J,2N)         │
│  → PricingResult R2           │
└───────────────┬───────────────┘
                ▼
┌───────────────────────────────┐
│ Run solver at (4J,4N)         │
│  → PricingResult R3           │
└───────────────┬───────────────┘
                ▼
┌───────────────────────────────────────────┐
│ Compute ratios:                            │
│  e1=|R2-R1|, e2=|R3-R2|                    │
│  p≈log2(e1/e2)                             │
│ Check p against expected (price/delta/gamma)│
└───────────────┬───────────────────────────┘
                ▼
┌───────────────────────────────┐
│ PASS / FAIL + store regression │
└───────────────────────────────┘
```

---

### 5.3 Regression test table

| Test ID | Product | Params | Target | Tolerance | Tests Failure Mode |
|---|---|---|---|---:|---|
| T1 | EU Call | as above | price, Δ, Γ | 1e-6 / 1e-5 / 1e-4 | F2,F8 |
| T2 | EU Digital | as above | price ≥ 0, oscScore | price 1e-5; oscScore < 0.02 | F2,F4 |
| T3 | Low-vol call | as above | negCount=0 | 0 | F1,F4,F5 |
| T5 | Disc DKO | as above | negCount=0, stability | 0; oscScore < 0.05 | F3,F4,F8 |
| T6 | AM Put | as above | LCP residual | < 1e-8 | F7 |

---

## SECTION 6: EXTENSION HOOKS (brief, structural only)

### 6.1 Jump-Diffusion (Merton PIDE)

- **Modules gaining strategies**
  - M5 SpatialOperator: add `IntegralOperator` strategy (jump term).
  - M7 TimeStepper: add `IMEXSplit` time stepper (diffusion implicit, integral explicit).
  - M4 VolSurfaceAdapter: unchanged.

- **New data structures**
  - `JumpKernel {lambda: float, muJ: float, sigJ: float, trunc: float}`
  - `IntegralStencil {weights: float[], offsets: int[]}`

- **ASCII diff**
```text
M7 TimeStepper
  before:  uses M5(tridiag) ─▶ M8
  after:   uses M5(tridiag) ─▶ M8
                 + M5(integral apply) (explicit RHS add)
```

---

### 6.2 Stochastic Volatility (Heston 2D)

- **2D-aware modules**
  - M1: build 2D grid (x,v) + time grid.
  - M5: assemble 2D operator (incl mixed derivative).
  - M8: TridiagSolver replaced by `BandedSolver` / `ADI_LineSolver`.
  - M7: owns ADI splitting strategy.

- **ADI ownership**
  - M7 TimeStepper adds: `ADI_Douglas`, `ADI_MCS`, `ADI_HV` strategies (extension).

- **ASCII diff**
```text
M7(TimeStepper/ADI)
  line solves in x using SolverX  (tridiag)
  line solves in v using SolverV  (tridiag)
  mixed term handled per ADI variant
```

---

### 6.3 Multi-Asset (2–3 factors)

- **Implications**
  - M8: `TridiagSolver` → `BandedSolver` / iterative solver interface.
  - M5: operator splitting module added (dimension-wise assembly).
  - M11: diagnostics extended to multi-dim oscillation metrics.

- **ASCII diff**
```text
M5 SpatialOperator (1D) ─▶ becomes M5n OperatorAssembler (d dims)
M8 TridiagSolver         ─▶ becomes LinearSolver (banded/iterative)
M7 TimeStepper           ─▶ uses splitting strategy + solver interface
```

---

### Quality-gate checklist (self-verified)

1. **M5 implementable with Stage1 §3.2**: yes (assemble ρ, ℓ/d/u; M-matrix check; fallback hooks).
2. **M7 implementable for RS-CN + TR-BDF2**: yes (RS-CN references Stage1 §3.3.1; TR-BDF2 fully specified in §3.1).
3. **All data flows explicit**: yes (no globals; structs on arrows and interfaces).
4. **Orchestrator covers EU vanilla, American vanilla, discrete monitored barrier**: yes (Boundary projection + damping restart; AmericanSpec drives M9).
5. **Every strategy variant named/distinguishable**: yes (module variant lists).
6. **Tests cover F1–F8**: yes (suite mapping table in §5.1).

## 1. Paper Identity

- **Title:** *A Crank-Nicolson finite difference approach on the numerical estimation of rebate barrier option prices*  
- **Author(s):** Nneka Umeorah; Phillip Mashele  
- **Affiliation(s):**  
  - Unit for Business Mathematics and Informatics, North-West University (NWU), Potchefstroom, South Africa  
  - Centre for Business Mathematics and Informatics, North-West University (NWU), Potchefstroom, South Africa  
- **Venue:** *Cogent Economics & Finance*  
- **Year:** 2019 (online publication listed as April 2019 in the paper front matter)  
- **DOI:** 10.1080/23322039.2019.1598835  
- **Other identifiers:** No arXiv ID reported in the paper.

## 2. Problem Statement & Formulation

Barrier options are **path-dependent** derivatives whose value depends on whether the underlying asset price hits a barrier during the contract’s life. The paper targets **rebate knock-out barrier options**, where a **rebate** is paid to the option holder when the barrier is breached (paid either immediately at knock-out or paid at expiry, depending on contract specification).

The paper’s technical objective is to compute (numerically, on a discrete space–time grid) the price surface $$V(t,S)$$ for a non-dividend-paying **down-and-out (DO) call** with a rebate, under the Black–Scholes model, using a **Crank–Nicolson (CN) finite difference method**, and to (i) compare CN prices to extended Black–Scholes closed-form prices for continuously monitored barriers and (ii) analyze CN behavior for sensitivities (Greeks), with emphasis on **spurious oscillations** induced by discontinuities.

### Risk-neutral underlying dynamics (Black–Scholes)

The underlying price process $$S(t)$$ is assumed to follow geometric Brownian motion:
$$
dS(t) = S(t)\bigl(r\,dt + \sigma\,dW(t)\bigr),
$$
where $$r$$ is the risk-free rate, $$\sigma$$ is volatility, and $$W(t)$$ is a standard Brownian motion.

### Pricing PDE and conditions (down-and-out call with rebate)

Let $$V(t,S)$$ denote the option value at time $$t$$ when the underlying price is $$S$$. Under the Black–Scholes framework, the paper states that $$V(t,S)$$ satisfies:
$$
\frac{\partial V}{\partial t} + rS\frac{\partial V}{\partial S} + \frac{\sigma^2 S^2}{2}\frac{\partial^2 V}{\partial S^2} = rV.
$$

For the DO call payoff and boundary/terminal conditions, the paper uses:

- **Terminal condition (European call payoff):**
$$
V(T,S) = \max(S-K,0),
$$
where $$K$$ is strike.

- **Barrier boundary condition (rebate-at-knock-out form):**
$$
V(t,B) = R,
$$
where $$B$$ is the down barrier and $$R$$ is the rebate amount (constant in the paper’s main numerical examples).

- **Far-field boundary condition (large underlying):**
$$
V(t,\infty) = S - K e^{-r(T-t)}.
$$

The PDE domain is described as:
$$
\mathcal{D} = \{(t,S)\;:\; B \leq S \leq \infty,\; t \in [0,T]\}.
$$

#### Rebate paid at expiry (boundary-condition variant)

For the “rebate paid at expiry” variant, the paper states that the barrier boundary is modified by discounting the rebate to the evaluation time, i.e., conceptually:
$$
V(t,B) = R e^{-r(T-t)}.
$$
Ambiguity flag: the discrete-time version in Section 3.1 is rendered in the provided text with a sign that appears OCR-corrupted; the interpretation above is consistent with the paper’s own definition of an expiry-paid rebate as $$R e^{-r(T-t)}$$ and with the narrative comparison “expiry-paid rebate is cheaper because it is discounted.”

### Risk-neutral expectation / first-passage decomposition

The paper provides an integral representation (risk-neutral valuation) that separates “no barrier hit” terminal payoff from “barrier hit” rebate via the first-passage time density. In the notation used:

- $$f(\cdot)$$ is the density of the underlying at maturity,
- $$g(p;B;t,S)$$ is the first-passage time density of hitting $$B$$ starting from $$S$$ at time $$t$$,
- $$\tau_B = \inf\{u : S(u)=B\}$$ is the first hitting time.

The qualitative structure is:
- a discounted terminal payoff integral over paths that do not hit the barrier before expiry,
- plus a discounted integral of expected rebate over the hitting time distribution.

## 3. Core Methodology

The core computational method is a **Crank–Nicolson (CN) finite difference method** applied to the Black–Scholes PDE with barrier and far-field boundary conditions, solved by backward iteration from maturity $$T$$ to current time $$t=0$$. The CN discretization yields, at each time step, a tridiagonal linear system in the spatial grid values of $$V$$.

### Overall pipeline (paper’s end-to-end workflow)

```text
┌──────────────────────────────────────────┐
│ Inputs                                   │
│ (S0,K,B,R,r,σ,T,Smax; grid sizes N,M)    │
└───────────────────────┬──────────────────┘
                        ▼
┌──────────────────────────────────────────┐
│ Model + PDE specification                │
│ GBM + Black-Scholes PDE on S∈[B,Smax]    │
│ Terminal payoff + barrier + far boundary │
└───────────────────────┬──────────────────┘
                        ▼
┌──────────────────────────────────────────┐
│ Uniform grid construction                │
│ S_k = k·ΔS,  t_i = i·Δt                  │
│ ΔS = Smax/M,  Δt = T/N                   │
└───────────────────────┬──────────────────┘
                        ▼
┌──────────────────────────────────────────┐
│ CN discretization of PDE                 │
│ Tridiagonal system each step             │
│ W·V_i = X·V_{i+1} (+ boundary terms)     │
└───────────────────────┬──────────────────┘
                        ▼
┌──────────────────────────────────────────┐
│ Backward time-marching                   │
│ for i = N-1, …, 0: solve tridiagonal     │
└───────────────────────┬──────────────────┘
                        ▼
┌──────────────────────────────────────────┐
│ Outputs                                  │
│ Price V(0,S0) and Greeks via FD          │
│ Oscillation diagnostics + timestep control│
└──────────────────────────────────────────┘
```

### Discretization setup (space–time grid)

The paper discretizes a computational domain $$[0,S_{\max}]\times[0,T]$$ with a **uniform** grid:

- Asset grid:
$$
S_k = k\,\Delta S,\quad k=0,1,\dots,M,\quad M\Delta S = S_{\max}.
$$

- Time grid:
$$
t_i = i\,\Delta t,\quad i=0,1,\dots,N,\quad N\Delta t = T.
$$

Grid function values:
$$
V_{i,k} = V(t_i,S_k).
$$

### Terminal and boundary conditions in grid form (down-and-out call)

For a rebate knock-out barrier option “rebate paid at knock-out,” the paper writes:

- Terminal:
$$
V_{N,k} = \max(S_k - K,0).
$$

- Barrier boundary (at $$S=B$$):
$$
V_{i,k_B} = R,
$$
where $$k_B$$ is implicitly the index such that $$S_{k_B}=B$$.  
Ambiguity flag: the paper does not describe what is done when $$B$$ is not exactly a grid point; the discrete condition is written as if $$B$$ is a node.

- Far boundary (at $$S=S_{\max}$$):
$$
V_{i,M} = S_{\max} - K e^{-r(T-t_i)}.
$$

For a rebate paid at expiry, the paper states the barrier boundary becomes discounted:
$$
V_{i,k_B} = R e^{-r(T-t_i)}.
$$

### CN scheme for the Black–Scholes PDE (procedural form)

The paper derives the CN discretization by averaging a time-implicit (spatial derivatives at the next time level) and a time-explicit (spatial derivatives at the current time level) discretization, yielding a tridiagonal relation of the form (paper’s Equation (3.6)):

For interior nodes $$k=1,2,\dots,M-1$$ and time stepping backward:
$$
\lambda_k V_{i-1,k-1} + (1-\beta_k)V_{i-1,k} + \eta_k V_{i-1,k+1}
=
\lambda_k V_{i,k-1} + (1+\beta_k)V_{i,k} + \eta_k V_{i,k+1}.
$$

The paper defines coefficient scalars (as written):
$$
\lambda_k = \frac{\Delta t}{4}\bigl(rk-\sigma^2 k^2\bigr),\quad
\beta_k = \frac{\Delta t}{2}\bigl(\sigma^2 k^2 + r\bigr),\quad
\eta_k = \frac{\Delta t}{4}\bigl(rk+\sigma^2 k^2\bigr).
$$
Ambiguity flag: these expressions are in index-space $$k$$ rather than explicitly in terms of $$S_k$$ and $$\Delta S$$; the paper later rewrites the scheme in matrix form using $$S_k/\Delta S$$, which indicates an implicit nondimensionalization consistent with the standard Black–Scholes FD derivations.

### Matrix form and per-step solve

The paper states the CN system is tridiagonal at each time step and can be solved efficiently.

A generic CN solve step (consistent with the paper’s matrix discussion) can be expressed as:

1. Build tridiagonal matrix $$W$$ from the left-hand side coefficients in $$V_{i-1,\cdot}$$.
2. Build tridiagonal matrix $$X$$ from the right-hand side coefficients in $$V_{i,\cdot}$$.
3. Enforce boundary values at $$S=B$$ and $$S=S_{\max}$$ in the linear system (either by modifying the right-hand side or by eliminating boundary nodes).
4. Solve:
$$
W\,\mathbf{V}_{i-1} = X\,\mathbf{V}_{i}
$$
for $$\mathbf{V}_{i-1}$$, where $$\mathbf{V}_i$$ is the vector of grid values at time index $$i$$.

### Greeks computation and oscillation diagnosis (as used in the paper)

The paper’s sensitivity study focuses on **Greeks** (specifically delta, gamma, theta), computed numerically from the grid. It reports that even if the price surface appears smooth, the Greek surfaces can show **spurious oscillations** near discontinuities (strike, barrier, and monitoring dates in discretely monitored barriers).

The paper’s diagnostic setup includes:
- **Continuous monitoring**: oscillations mainly at strike due to payoff kink.
- **Discrete monitoring**: oscillations at strike and near barrier.

### Oscillation-mitigation via timestep restriction (paper’s approach)

The paper’s main stabilization lever is a “special timestep restriction” motivated by:
- characteristic grid diffusion reasoning (via a cited condition involving a characteristic diffusion time),
- enforcing positivity/maximum principle properties via M-matrix arguments (following Milev & Tagliani, 2010, as cited in the paper).

The paper’s stated restriction includes:
$$
\Delta t < \frac{2}{r + 2(\sigma m)^2},
$$
where $$m$$ is the number of spatial steps (paper’s notation).

It also highlights a regime condition:
$$
\sigma^2 > r,
$$
as part of the cited hypotheses supporting positivity and eigenvalue properties (see Section 4 below).

## 4. Theoretical Results

The paper’s “theory” consists of (i) a convergence argument for CN via stability and consistency on a transformed heat equation and (ii) conditions (positivity / eigenvalue properties) used to explain and mitigate spurious oscillations in Greeks.

No items are labeled “Theorem/Lemma/Proposition/Corollary” by the authors; the statements below reproduce the complete claims that are stated as such in the text (including named theorems invoked by the authors).

### Theorem (Lax Equivalence Theorem; invoked)

**Statement (as used by the paper):** Given a well-posed initial value problem, a consistent finite-difference method is convergent if and only if the method is stable.

**Proof strategy sketch:** The paper treats this as a standard result and uses it as a template: (i) show stability of the CN discretization for a canonical PDE (heat equation) and (ii) assert consistency via truncation-error order, concluding convergence.

### Proposition (Unconditional stability of CN for the heat equation; paper’s eigenvalue argument)

The paper transforms the Black–Scholes PDE to a heat equation:
$$
\frac{\partial u}{\partial \tau} = \frac{\partial^2 u}{\partial x^2},\quad -\infty < x < \infty,\quad \tau>0,
$$
under the stated change of variables:
$$
x = \log S + \bigl(r-\tfrac{1}{2}\sigma^2\bigr)(T-t),\qquad \tau = \tfrac{1}{2}\sigma^2(T-t),
$$
and a corresponding mapping $$u(x,\tau)$$ built from $$V(S,t)$$ (paper provides an explicit expression).

The CN discretization of the heat equation is written (paper’s Equation (3.8)):
$$
\frac{u_{i,j+1}-u_{i,j}}{\Delta\tau}
=
\frac{u_{i+1,j+1}-2u_{i,j+1}+u_{i-1,j+1}}{2\Delta x^2}
+
\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{2\Delta x^2}.
$$

Let
$$
\zeta = \frac{\Delta \tau}{2\Delta x^2}.
$$
The scheme is rearranged to:
$$
A\,\mathbf{u}^{j+1} = B\,\mathbf{u}^{j},
$$
with tridiagonal Toeplitz matrices:
$$
A = \operatorname{tridiag}(-\zeta,\; 1+2\zeta,\; -\zeta),\qquad
B = \operatorname{tridiag}(\zeta,\; 1-2\zeta,\; \zeta).
$$

The paper then uses the eigenvalues of the Toeplitz structure to state:
$$
\lambda_j(A) = 1 + 4\zeta \sin^2\!\Bigl(\frac{\pi j}{2M}\Bigr),\qquad j=1,\dots,M-1,
$$
and concludes that the stability condition
$$
\left|\frac{2}{\lambda_j(A)} - 1\right| < 1
\quad\Longleftrightarrow\quad
\left|\frac{1 - 4\zeta \sin^2(\pi j/(2M))}{1 + 4\zeta \sin^2(\pi j/(2M))}\right| < 1
$$
holds for all $$\zeta>0$$, so CN is stable for any $$\Delta \tau>0$$.

**Proof strategy sketch:** The argument is a spectral-radius bound for the linear iteration induced by CN. The key observation is that $$\lambda_j(A)>1$$ for $$\zeta>0$$, so the amplification factors remain inside the unit disk for all modes, yielding unconditional stability.

### Proposition (Consistency order and convergence of CN; truncation-error statement)

The paper states that CN exhibits “quadratic convergence in both space and time,” expressed as a local truncation error order:
$$
O\bigl((\Delta x)^2 + (\Delta \tau)^2\bigr).
$$

Combined with unconditional stability and the Lax Equivalence Theorem, the paper concludes that the CN discretization is convergent.

**Proof strategy sketch:** CN for diffusion-type PDEs is standardly second-order accurate in both time and space for smooth solutions. The paper’s logic is: truncation error establishes consistency; eigenvalue analysis establishes stability; Lax equivalence gives convergence.

### Proposition (Positivity / maximum principle / eigenvalue conditions used to suppress oscillations; as quoted from Milev & Tagliani-style conditions)

The paper rewrites the CN scheme in matrix form:
$$
W\mathbf{V}^{i-1} = X\mathbf{V}^{i},
\quad\text{so}\quad
\mathbf{V}^{i-1} = W^{-1}X\mathbf{V}^{i},
$$
and then states (citing Milev & Tagliani, 2010) that if the discretization satisfies the hypothesis:
$$
\frac{\sigma^2}{r} > 1,
$$
then the following properties are satisfied:

1. **Positivity:** $$V_{i-1,k} > 0$$, and under $$r<\sigma^2$$ one has an inequality of the form:
$$
\frac{rS_k}{\Delta S} - \Bigl(\frac{\sigma S_k}{\Delta S}\Bigr)^2 < 0.
$$
The matrix $$W$$ is described as an irreducible row-diagonally-dominant **M-matrix**, implying:
$$
W^{-1} > 0,\qquad \|W^{-1}\|_\infty \leq \frac{1}{1 + r\Delta t/2}.
$$

2. **Maximum principle:** With $$X \geq 0$$ and
$$
\|X\|_\infty = 1 - r\Delta t/2 > 0,
$$
one obtains a contraction-type bound:
$$
\|\mathbf{V}^{i-1}\|_\infty
=
\|W^{-1}X\mathbf{V}^{i}\|_\infty
\leq
\frac{1 - r\Delta t/2}{1 + r\Delta t/2}\,\|\mathbf{V}^{i}\|_\infty
\leq
\|\mathbf{V}^{i}\|_\infty.
$$

3. **Eigenvalue properties under timestep restriction:** Together with:
$$
\Delta t < \frac{2}{r + 2(\sigma m)^2},
$$
the paper states that eigenvalues of the iteration matrix are real and distinct, and by Gershgorin-type bounds the eigenvalues of $$W^{-1}X$$ lie inside:
$$
\left[\frac{1}{1 + r\Delta t/2},\; \frac{1}{1 + \Delta t\bigl(r/2 + (\sigma m)^2\bigr)}\right] \subset (0,1).
$$

The paper concludes: if the timestep restriction and the condition $$\sigma^2>r$$ are ignored, positivity is not ensured and negative/complex eigenvalues can occur, leading to spurious oscillations (especially visible in Greeks).

**Proof strategy sketch:** The logic is matrix-theoretic: (i) show $$W$$ is an M-matrix so $$W^{-1}$$ is nonnegative and bounded, (ii) show $$X$$ is elementwise nonnegative with norm bounded by a discount-like factor, (iii) conclude the iteration is monotone/contractive in the maximum norm, and (iv) use Gershgorin bounds plus timestep restriction to keep eigenvalues inside $$ (0,1)$$, which suppresses oscillatory modes.

### Condition set (Péclet-type conditions; stated as sufficient to avoid oscillations with central differencing)

The paper cites conditions (from Zvan, Forsyth, and Vetzal, 1998) under which central weighting avoids oscillations by ensuring positivity and applicability of the maximum principle. The paper gives two inequalities involving $$\Delta S_{i-1/2}$$, $$\Delta S_{i+1/2}$$, $$\Delta t$$, and local coefficients.

Ambiguity flag: the exact algebraic form of the second inequality is partially garbled in the provided text rendering; the paper’s conceptual use is: satisfy these conditions to ensure positivity and prevent spurious oscillations in CN when central differencing is used.

## 5. Experimental Evaluation

### Experimental configurations (datasets, baselines, metrics, hyperparameters)

The paper uses synthetic parameter sets (no empirical dataset) and evaluates numerical price convergence and Greek stability against analytic (extended Black–Scholes) prices and Monte Carlo estimators.

| Experiment | Contract(s) | Monitoring | Rebate timing | Baselines | Metrics reported | Key hyperparameters / settings |
|---|---|---|---|---|---|---|
| E1 (Table 1) | DO call, rebate | Continuous (per closed form) | At knock-out | Closed-form (extended BS, Eq. (2.7)) | CN price, computation time | $$S=50$$, $$S_{\max}=140$$, $$K=40$$, $$B=20$$, $$r=0.04$$, $$\sigma=0.3$$, $$T=0.5$$, $$R=2.5$$; CN grids $$N=M \in [150,500]$$ |
| E2 (Table 2) | DO call, rebate | Continuous (per closed form) | At knock-out | Closed-form (Eq. (2.7)); MCS; Antithetic MCS | Prices vs exact | Same as E1; CN uses $$N=M=400$$; MCS/AMCS use 100000 simulations with 400 time steps; varying $$S \in \{70,65,\dots,35\}$$ |
| E3 (Table 3) | DO call, rebate | Continuous (per closed form) | At knock-out | Closed-form (Eq. (2.7)) | CN price convergence vs exact | $$S=100$$, $$S_{\max}=260$$, $$K=100$$, $$B=60$$, $$r=0.08$$, $$\sigma=0.1$$, $$T=0.5$$, rebate $$R=0.04S$$; multiple grid aspect ratios and $$N=M$$ sweeps |
| E4 (Fig. 2–4) | DO barrier option (zero-rebate in Greek plots) | Continuous and discrete monitoring (5 dates shown) | N/A in figures (focus: oscillations) | CN without restriction vs CN with restriction | Greek smoothness; qualitative oscillation patterns | $$S=60$$, $$K=50$$, $$B=35$$, $$r=0.05$$, $$\sigma=0.2$$, $$T=0.75$$, $$M=150$$, $$N=25$$, $$S_{\max}=140$$ |
| E5 (Fig. 5–7, Table 4) | Vanilla call, ZRDO, RDO, RDOE | Continuous (per narrative) | Knock-out vs expiry | Closed-form comparisons (implied) | Comparative price curves; rebate-value dependence | $$B=120$$, $$K=125$$, $$r=0.06$$, $$\sigma=0.5$$, $$T=2$$; rebate rates 5%, 12%, 17% of $$S$$ |

### Key quantitative results (prices, convergence, runtime)

#### Table 1 (CN convergence and runtime; rebate at knock-out)

Parameters: $$S=50$$, $$S_{\max}=140$$, $$K=40$$, $$B=20$$, $$r=0.04$$, $$\sigma=0.3$$, $$T=0.5$$, rebate $$R=2.5$$. Exact price reported: 11.3777.

| $$N$$ | $$M$$ | CN value (CNV) | Computation time (seconds) |
|---:|---:|---:|---:|
| 150 | 150 | 11.4090 | 0.1722 |
| 200 | 200 | 11.3875 | 0.3572 |
| 250 | 250 | 11.3818 | 0.7278 |
| 300 | 300 | 11.3787 | 1.1481 |
| 350 | 350 | 11.3786 | 2.3930 |
| 400 | 400 | 11.3780 | 3.8790 |
| 450 | 450 | 11.3777 | 6.0094 |
| 500 | 500 | 11.3777 | 8.7626 |

Reported conclusion: increasing $$N$$ and $$M$$ together drives CNV toward the exact value; runtime increases monotonically with grid refinement.

#### Table 2 (CN vs Monte Carlo; varying spot)

Parameters as Table 1; CN uses $$N=M=400$$; MCS and AMCS use 100000 simulations and 400 time steps. Values reproduced:

| $$S$$ | Exact | CNV | MCS | AMCS |
|---:|---:|---:|---:|---:|
| 70 | 30.8026 | 30.5945 | 30.8960 | 30.8853 |
| 65 | 25.8226 | 25.7657 | 25.7487 | 25.7761 |
| 60 | 20.8777 | 20.8655 | 20.8328 | 20.8870 |
| 55 | 16.0225 | 16.0208 | 16.0336 | 16.0103 |
| 50 | 11.3777 | 11.3780 | 11.3847 | 11.3768 |
| 45 | 7.1737 | 7.1755 | 7.1558 | 7.1622 |
| 40 | 3.7590 | 3.7635 | 3.7883 | 3.7408 |
| 35 | 1.4876 | 1.4866 | 1.4655 | 1.4670 |

Reported conclusion: CN values are “fairly close” to exact values; AMCS improves over MCS via variance reduction.

#### Table 3 (grid aspect ratio effects; slow convergence unless $$N=M$$ increases)

Parameters: $$S=100$$, $$S_{\max}=260$$, $$K=100$$, $$B=60$$, $$r=0.08$$, $$\sigma=0.1$$, $$T=0.5$$. Exact value reported for rebate $$R=0.04S$$ (paid at knock-out): 5.1563.

The paper reports three grid regimes:

- Regime A: $$M = 2N$$ (time steps half asset steps)  
- Regime B: $$M = N$$  
- Regime C: $$N = 2M$$ (time steps double asset steps)

| Regime A: $$N$$ | $$M$$ | CNV | Regime B: $$N$$ | $$M$$ | CNV | Regime C: $$M$$ | $$N$$ | CNV |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 300 | 600 | 5.1556 | 300 | 300 | 5.2787 | 300 | 600 | 5.1556 |
| 325 | 650 | 5.1556 | 325 | 325 | 5.2393 | 325 | 650 | 5.1556 |
| 350 | 700 | 5.1557 | 350 | 350 | 5.2057 | 350 | 700 | 5.1557 |
| 375 | 750 | 5.1558 | 375 | 375 | 5.1832 | 375 | 750 | 5.1558 |
| 400 | 800 | 5.1558 | 400 | 400 | 5.1700 | 400 | 800 | 5.1558 |
| 425 | 850 | 5.1559 | 425 | 425 | 5.1628 | 425 | 850 | 5.1559 |
| 450 | 900 | 5.1559 | 450 | 450 | 5.1590 | 450 | 900 | 5.1559 |
| 475 | 950 | 5.1560 | 475 | 475 | 5.1572 | 475 | 950 | 5.1560 |
| 500 | 1000 | 5.1560 | 500 | 500 | 5.1563 | 500 | 1000 | 5.1560 |

Reported conclusion: convergence is slow; using equal refinement $$N=M$$ and increasing both yields convergence to the exact solution (here at $$N=M=500$$).

#### Table 4 (rebate magnitude increases option value; expiry-paid rebate cheaper)

Parameters: $$B=120$$, $$K=125$$, $$r=0.06$$, $$\sigma=0.5$$, $$T=2$$. Table reports a “zero rebate option value” (ZRDO) and “rebate terms” for rebate rates 5%, 12%, 17% of the underlying.

Interpretation consistent with the paper’s narrative: the “rebate term” columns quantify the value contribution attributable to rebates, which increases with rebate rate and decreases when rebate is paid at expiry due to discounting.

| $$S$$ | Zero rebate option value | Rebate at knock-out (5%) | (12%) | (17%) | Rebate at expiry (5%) | (12%) | (17%) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 200 | 87.3962 | 5.0691 | 12.1659 | 17.2350 | 4.4959 | 10.7902 | 15.2861 |
| 190 | 77.0044 | 5.2365 | 12.5676 | 17.8041 | 4.6443 | 11.1464 | 15.7908 |
| 180 | 66.5247 | 5.3973 | 12.9535 | 18.3507 | 4.7870 | 11.4887 | 16.2756 |
| 170 | 55.9353 | 5.5484 | 13.3161 | 18.8644 | 4.9210 | 11.8103 | 16.7313 |
| 160 | 45.2082 | 5.6860 | 13.6464 | 19.3324 | 5.0430 | 12.1033 | 17.1463 |
| 150 | 34.3070 | 5.8057 | 13.9338 | 19.7395 | 5.1492 | 12.3582 | 17.5074 |
| 140 | 23.1841 | 5.9023 | 14.1654 | 20.0677 | 5.2348 | 12.5636 | 17.7985 |
| 130 | 11.7765 | 5.9694 | 14.3266 | 20.2960 | 5.2944 | 12.7065 | 18.0009 |
| 120 | 0.0000 | 6.0000 | 14.4000 | 20.4000 | 5.3215 | 12.7717 | 18.0932 |

Reported conclusion: larger rebates increase option value; rebate-at-expiry values are lower than rebate-at-knock-out values due to discounting.

### Ablation/sensitivity axes and conclusions

- **Grid refinement ablation:** varying $$N$$ and $$M$$ (Tables 1 and 3) shows convergence to analytic price as grid is refined, with runtime trade-offs.
- **Method ablation:** comparing CN vs MCS vs AMCS (Table 2) shows CN close to analytic, and AMCS improves over MCS.
- **Rebate magnitude ablation:** increasing rebate percentage increases option value contribution (Table 4; Figures 6–7).
- **Oscillation ablation:** CN Greeks exhibit oscillations without timestep restrictions; applying the paper’s restrictions yields smoother Greeks (Figures 2 vs 4).

## 6. ASCII Architecture / Workflow Diagram(s)

### CN barrier-option PDE solver (grid + tridiagonal solve)

```text
┌─────────────────────────────────────────────────────────────────┐
│ PDE: V_t + rS V_S + 0.5 σ^2 S^2 V_SS = rV                       │
│ S∈[B,Smax], t∈[0,T]                                              │
└───────────────────────────────────────┬─────────────────────────┘
                                        ▼
┌─────────────────────────────────────────────────────────────────┐
│ Grid                                                           │
│ S_k = k·ΔS,  k=0..M      (ΔS = Smax/M)                          │
│ t_i = i·Δt,  i=0..N      (Δt = T/N)                              │
└───────────────────────────────────────┬─────────────────────────┘
                                        ▼
┌─────────────────────────────────────────────────────────────────┐
│ Conditions                                                     │
│ Terminal:  V_{N,k} = max(S_k - K, 0)                            │
│ Barrier:   V_{i,kB} = R   or   R·exp(-r·(T - t_i))              │
│ Far field: V_{i,M} = Smax - K·exp(-r·(T - t_i))                 │
└───────────────────────────────────────┬─────────────────────────┘
                                        ▼
┌─────────────────────────────────────────────────────────────────┐
│ CN step (for each i = N-1,…,0)                                  │
│ Build tridiagonal W, X so that:                                 │
│     W · V_i  =  X · V_{i+1}  (+ boundary adjustments)           │
│ Solve tridiagonal linear system (Thomas / band solver)          │
└───────────────────────────────────────┬─────────────────────────┘
                                        ▼
┌─────────────────────────────────────────────────────────────────┐
│ Outputs                                                        │
│ Price:  V(0,S0) via interpolation on {S_k}                      │
│ Greeks: finite differences on the grid                          │
│ Diagnostics: oscillations near strike/barrier; enforce Δt rule  │
└─────────────────────────────────────────────────────────────────┘
```

### Oscillation control logic (as used in the paper’s Greek study)

```text
┌──────────────────────────────┐
│ Discontinuities present?      │
│ payoff kink at K; barrier at B│
└───────────────┬──────────────┘
                ▼
┌──────────────────────────────┐
│ CN Greeks show oscillations?  │
│ (spikes near K; near B if     │
│ discrete monitoring)          │
└───────────────┬──────────────┘
                ▼
┌──────────────────────────────────────────────────┐
│ Enforce timestep restriction + positivity logic   │
│ Δt < 2 / ( r + 2(σ m)^2 )   and   σ^2 > r         │
│ (paper’s stated sufficient conditions)            │
└───────────────┬──────────────────────────────────┘
                ▼
┌──────────────────────────────┐
│ Recompute Greeks on refined/  │
│ restricted grid (smoother)    │
└──────────────────────────────┘
```

## 7. Follow-Up Works & Extensions

### A. Direct extensions that cite and build on the CN rebate-barrier baseline

Umeorah and Mba develop a feed-forward neural-network approach that reframes single-barrier Black–Scholes PDE valuation as a constrained optimization problem with a trial solution designed to satisfy boundary/terminal conditions; numerical outputs are benchmarked against Monte Carlo, Crank–Nicolson finite-difference values, and exact Black–Scholes prices, positioning CN (as in the 2019 paper) as a key numerical reference point for accuracy and speed comparisons. [Umeorah & Mba, Applied Stochastic Models in Business and Industry 2022] 

Umeorah, Mashele, Agbaeze, and Mba propose a neural-network surrogate that learns the extended Black–Scholes barrier-option pricing map and extracts Greeks, explicitly citing the 2019 CN rebate-barrier study as a PDE-based benchmark for rebate barrier options; the paper’s framing makes the 2019 CN approach a baseline “ground truth generator” for supervised learning on synthetic data and a reference for oscillation/Greek-smoothness considerations. [Umeorah et al., Axioms 2023] 

### B. Generalizations of rebate-option structures (multi-step barriers; product embedding)

Lee, Jeong, and Lee extend rebate options from constant barriers to multi-step (time-varying) barrier structures and derive closed-form formulas under Black–Scholes, explicitly referencing the 2019 CN rebate-option paper as prior work using CN discretization; the paper further embeds these generalized rebate options into equity-linked product designs (including annuity-style structures) and validates formulas numerically. [Lee et al., NAJEF 2023] 

Lee, Ha, Gaeun Lee, and Minha Lee propose approximating American option prices by optimizing explicit closed-form multi-step rebate-option pricing formulas, replacing unknown optimal exercise boundaries with optimized multi-step barriers; the connection to the 2019 paper is thematic (rebate options as first-touch payments and the barrier/first-hitting-time viewpoint), and the work advances the “rebate as an American-style digital” interpretation that the 2019 paper uses to motivate rebate terms. [Lee et al., NAJEF 2024] 

### C. Independent PDE/ML approaches for barrier options under richer dynamics

Fu and Hirsa propose an unsupervised deep-learning method to solve barrier-option PDEs under a stochastic-volatility setting (Bergomi model), adding singular terms to address the strike/barrier non-smoothness; the work is aligned with the 2019 paper’s core numerical issue (discontinuities driving stability/oscillation artifacts in sensitivities), but it replaces CN time-marching with a PDE-constrained neural approximation aimed at fast evaluation after training. [Fu & Hirsa, arXiv 2022] 

## 8. Industrial & Real-World Applications

QuantLib provides production-grade open-source implementations for barrier and rebate features in option pricing, including finite-difference engines for Black–Scholes barrier options and a dedicated finite-difference rebate engine; these implementations operationalize the same PDE discretization philosophy as the paper (solve the pricing PDE on a grid with barrier/rebate boundary conditions) and expose barrier options with an explicit rebate parameter in the instrument interface, making CN-style PDE solvers practically usable inside a broader risk and valuation framework. [GitHub: lballabio/QuantLib] 

RQuantLib exposes QuantLib’s option-pricing capabilities to R users and includes a barrier-option valuation interface with an explicit rebate parameter; this provides an applied pathway to compute barrier-option prices and Greeks (via QuantLib’s analytic implementations) and to prototype/benchmark rebate barrier option values in statistical workflows, complementing PDE-based implementations like the paper’s CN method. [GitHub: eddelbuettel/rquantlib] 

The PROJ_Option_Pricing_Matlab repository provides an open-source research-and-development codebase covering barrier options (single/double) and rebates across multiple numerical paradigms (Fourier/PROJ, PDE/finite difference, Monte Carlo, lattice); this constitutes a practical implementation ecosystem in which CN-style PDE solvers and rebate barrier option valuation can be tested and compared against alternative methods under many models. [GitHub: jkirkby3/PROJ_Option_Pricing_Matlab] 

## 9. Consolidated Reference List

[1] N. Umeorah, J. C. Mba. “Approximation of single-barrier options partial differential equations using feed-forward neural network.” *Applied Stochastic Models in Business and Industry*, 2022. DOI: 10.1002/asmb.2711. 

[2] N. Umeorah, P. Mashele, O. Agbaeze, J. C. Mba. “Barrier Options and Greeks: Modeling with Neural Networks.” *Axioms*, 2023. DOI: 10.3390/axioms12040384. 

[3] W. Fu, A. Hirsa. “Solving barrier options under stochastic volatility using deep learning.” *arXiv* preprint, 2022. arXiv:2207.00524. URL: `https://arxiv.org/abs/2207.00524`. 

[4] H. Lee, H. Jeong, G. Lee. “Valuing rebate options and equity-linked products.” *The North American Journal of Economics and Finance*, 2023, Volume 68, Article 101968. DOI: 10.1016/j.najef.2023.101968. 

[5] H. Lee, H. Ha, G. Lee, M. Lee. “Valuing American options using multi-step rebate options.” *The North American Journal of Economics and Finance*, 2024, Volume 74, Article 102227. DOI: 10.1016/j.najef.2024.102227. 

[6] QuantLib contributors. “QuantLib: The QuantLib C++ library.” GitHub repository. URL: `https://github.com/lballabio/QuantLib`. 

[7] QuantLib source browser. “fdblackscholesbarrierengine.cpp (QuantLib barrier FD engine source).” Codebrowser view of QuantLib source. URL: `https://codebrowser.dev/quantlib/quantlib/ql/pricingengines/barrier/fdblackscholesbarrierengine.cpp.html`. 

[8] QuantLib source browser. “fdblackscholesrebateengine.cpp (QuantLib rebate FD engine source).” Codebrowser view of QuantLib source. URL: `https://codebrowser.dev/quantlib/quantlib/ql/pricingengines/barrier/fdblackscholesrebateengine.cpp.html`. 

[9] D. Eddelbuettel et al. “RQuantLib: R interface to the QuantLib library.” GitHub repository. URL: `https://github.com/eddelbuettel/rquantlib`. 

[10] D. Eddelbuettel et al. “BarrierOption documentation (rebate parameter exposed in RQuantLib).” R documentation mirror for the GitHub repo. URL: `https://rdrr.io/github/eddelbuettel/rquantlib/man/BarrierOption.html`. 

[11] J. L. Kirkby. “PROJ_Option_Pricing_Matlab: Quant Option Pricing – Exotic/Vanilla (Barrier options, rebates, PDE/finite difference, Fourier/PROJ).” GitHub repository. URL: `https://github.com/jkirkby3/PROJ_Option_Pricing_Matlab`.

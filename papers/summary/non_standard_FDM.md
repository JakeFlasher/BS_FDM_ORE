## 1. Paper Identity

- **Title:** *Nonstandard Finite Difference Schemes with Application to Finance: Option Pricing*
- **Authors:** Mariyan Milev; Aldo Tagliani
- **Affiliations (as stated in the paper):**
  - Mariyan Milev — Department of Applied Mathematics, University of Venice (Venice, Italy)
  - Aldo Tagliani — Department of Computer and Management Sciences, Trento University (Trento, Italy)
- **Venue / bibliographic placement:** *Serdica Mathematical Journal*, pp. 75–88, year 2010. ([math.bas.bg](https://www.math.bas.bg/serdica/n1_10.html?utm_source=openai))  
  - **Volume numbering discrepancy (verifiable):** the Serdica table-of-contents page lists **Volume 36, Number 1 (2010)** for this article, while the EuDML record lists it as **Volume 35, Issue 1 (2010)**. ([math.bas.bg](https://www.math.bas.bg/serdica/n1_10.html?utm_source=openai))
- **Year:** 2010 ([math.bas.bg](https://www.math.bas.bg/serdica/n1_10.html?utm_source=openai))
- **DOI / arXiv ID:** none found in the Serdica/EuDML records for this article (paper predates widespread DOI assignment in some journals and is not listed with a DOI on EuDML). ([eudml.org](https://eudml.org/doc/281413?utm_source=openai))

## 2. Problem Statement & Formulation

The paper studies finite-difference valuation of option contracts whose payoff/monitoring rules induce **discontinuities in the initial data** of the associated Black–Scholes PDE at maturity and at discrete monitoring dates. The target failure modes are (i) **loss of positivity** (negative option values) and (ii) **spurious oscillations** near sharp gradients (barriers and strike), especially for widely used centered/second-order schemes such as Crank–Nicolson.

### Contract and dynamics

Under the risk-neutral measure, the underlying follows a geometric Brownian motion with constant coefficients:
$$
\frac{dS}{S} = r\,dt + \sigma\,dW_t,
$$
with constant risk-free rate $$r$$ and volatility $$\sigma$$.

The contract used as the concrete exemplar is a **discretely monitored double barrier knock-out call**. The payoff is call-like inside a corridor $$[L,U]$$ and zero outside; the option is knocked out at monitoring dates if $$S$$ lies outside the corridor.

### PDE between monitoring dates

Let $$V(S,t)$$ denote the option value, where the paper’s time variable $$t$$ is “time to expiry” on $$[0,T]$$, with the payoff posed as an “initial” condition at $$t=0$$. Between monitoring dates, $$V$$ satisfies the Black–Scholes PDE:
$$
-\frac{\partial V}{\partial t}
+ rS\frac{\partial V}{\partial S}
+ \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}
- rV
= 0.
$$

### Initial and boundary conditions; monitoring update

Initial condition at $$t=0$$:
$$
V(S,0) = \max(S-K,0)\,\mathbf{1}_{[L,U]}(S).
$$

Far-field boundary condition (as stated):
$$
V(S,t)\to 0 \quad \text{as} \quad S\to 0 \ \text{or}\ S\to\infty.
$$

Discrete monitoring dates $$0=t_0<t_1<\cdots<t_F=T$$ impose a projection (“knock-out update”) that renews discontinuities at $$S=L$$ and $$S=U$$:
$$
V(S,t_i) = V(S,t_i^-)\,\mathbf{1}_{[L,U]}(S),\qquad i=1,\ldots,F.
$$

### Positivity and (discrete) maximum-principle motivation

The paper treats nonnegativity $$V(S,t)\ge 0$$ as a financial requirement and uses a maximum-principle-style condition on the sup norm over a truncated computational interval $$[0,S_{\max}]$$:
$$
\max_{S\in[0,S_{\max}]} |V(S,t_1)|
\ge
\max_{S\in[0,S_{\max}]} |V(S,t_2)|,
\qquad
t_1\le t_2,
$$
as a diagnostic: violating the discrete analogue correlates with spurious oscillations near discontinuities.

## 3. Core Methodology

### Discretization setup (common to all schemes)

- Truncate asset domain to $$[0,S_{\max}]$$ with $$S_{\max}$$ “sufficiently large” so boundaries minimally affect interior values.
- Use a uniform grid:
  - $$S_j = j\,\Delta S,\quad j=0,1,\ldots,M,$$ with $$S_{\max}=M\,\Delta S$$.
  - $$t_n = n\,\Delta t,\quad n=0,1,\ldots,X,$$ with $$T=X\,\Delta t$$.
- Let $$V_j^n \approx V(S_j,t_n)$$ and (for interior nodes) collect $$V^n$$ into a vector over spatial indices.
- Between monitoring dates, advance the solution by solving a tridiagonal linear system each time step of the form:
  $$
  P\,V^{n+1} = N\,V^n.
  $$
- At each monitoring date, apply the knock-out projection:
  $$
  V_j \leftarrow V_j \,\mathbf{1}_{[L,U]}(S_j).
  $$

### Overall pipeline (pricing with discrete monitoring)

```
┌──────────────────────────────────────────────┐
│ Input: r, σ, K, L, U, T, {t_i}; Smax, ΔS, Δt │
└──────────────────────────────┬───────────────┘
                               ▼
┌──────────────────────────────────────────────┐
│ Mesh: S_j=jΔS (j=0..M), t_n=nΔt (n=0..X)     │
└──────────────────────────────┬───────────────┘
                               ▼
┌──────────────────────────────────────────────┐
│ Initialize at t=0: V_j^0=(S_j-K)^+ 1_[L,U](S_j)│
└──────────────────────────────┬───────────────┘
                               ▼
┌──────────────────────────────────────────────┐
│ For n=0..X-1 (piecewise between monitoring): │
│   Solve: P V^{n+1} = N V^n                    │
│   If t_{n+1} is a monitoring time:            │
│       V^{n+1} ← V^{n+1} ⊙ 1_[L,U](S)          │
└──────────────────────────────┬───────────────┘
                               ▼
┌──────────────────────────────────────────────┐
│ Output: grid values {V_j^X}; interpolate to S0│
└──────────────────────────────────────────────┘
```

### Algorithm A — Discrete double-barrier knock-out call pricing (FD + monitoring)

**Inputs:** $$r,\sigma,K,L,U,T$$; monitoring schedule $$\{t_i\}_{i=0}^F$$; truncation $$S_{\max}$$; mesh sizes $$\Delta S,\Delta t$$; choice of time-stepping scheme (Crank–Nicolson or proposed semi-implicit).

**Outputs:** approximate price $$V(S_0,T)$$ (and optionally the full grid $$\{V_j^n\}$$).

1. Set $$M = S_{\max}/\Delta S$$ and $$X=T/\Delta t$$ (integers as assumed by the paper).
2. Create $$S_j=j\Delta S$$ for $$j=0,\ldots,M$$ and $$t_n=n\Delta t$$ for $$n=0,\ldots,X$$.
3. Initialize:
   $$
   V_j^0 \leftarrow \max(S_j-K,0)\,\mathbf{1}_{[L,U]}(S_j).
   $$
4. For $$n=0,1,\ldots,X-1$$:
   1. Form tridiagonal matrices $$P$$ and $$N$$ (scheme-dependent; see Algorithms B and C below).
   2. Solve the tridiagonal linear system:
      $$
      P V^{n+1} = N V^n.
      $$
   3. If $$t_{n+1}$$ coincides with a monitoring date $$t_i$$:
      $$
      V_j^{n+1} \leftarrow V_j^{n+1}\,\mathbf{1}_{[L,U]}(S_j)\quad \forall j.
      $$
5. Return $$V(S_0,T)$$ via interpolation from $$\{V_j^X\}$$ (interpolation policy not specified in the provided text; nearest-grid evaluation is possible when $$S_0$$ lies on-grid).

**Termination criterion:** fixed number of steps $$X=T/\Delta t$$; no iterative convergence criterion is used beyond solving each linear system.

### Algorithm B — Crank–Nicolson step (as analyzed in Section 3.1)

Crank–Nicolson uses centered differences and yields:
$$
P\,V^{n+1} = N\,V^n,
$$
with $$P$$ and $$N$$ tridiagonal. Writing $$\operatorname{tridiag}\{a_j;\,b_j;\,c_j\}$$ for a tridiagonal matrix with lower diagonal $$a_j$$, main diagonal $$b_j$$, upper diagonal $$c_j$$, the paper’s coefficients (expressed in a numerically consistent form matching the analysis in Theorem 3.1) can be written as:

$$
P = \operatorname{tridiag}\left\{
\frac{r}{4}\frac{S_j}{\Delta S} - \left(\frac{\sigma S_j}{2\Delta S}\right)^2\ ;\
\frac{1}{\Delta t} + \frac{1}{2}\left(\frac{\sigma S_j}{\Delta S}\right)^2 + \frac{r}{2}\ ;\
-\frac{r}{4}\frac{S_j}{\Delta S} - \left(\frac{\sigma S_j}{2\Delta S}\right)^2
\right\},
$$

$$
N = \operatorname{tridiag}\left\{
-\frac{r}{4}\frac{S_j}{\Delta S} + \left(\frac{\sigma S_j}{2\Delta S}\right)^2\ ;\
\frac{1}{\Delta t} - \frac{1}{2}\left(\frac{\sigma S_j}{\Delta S}\right)^2 - \frac{r}{2}\ ;\
\frac{r}{4}\frac{S_j}{\Delta S} + \left(\frac{\sigma S_j}{2\Delta S}\right)^2
\right\}.
$$

**Per-step procedure:**

**Inputs:** $$V^n$$, parameters $$r,\sigma$$, mesh $$\Delta S,\Delta t$$, index set $$\{S_j\}$$.

**Output:** $$V^{n+1}$$.

1. Assemble tridiagonal $$P$$ and $$N$$ at time level $$n$$ using the formulas above.
2. Compute the right-hand side $$b \leftarrow N V^n$$.
3. Solve $$P V^{n+1} = b$$ (tridiagonal linear solve).

**Convergence / stability notes (from the paper’s analysis):**
- Positivity and a discrete maximum principle require additional time-step restrictions (Theorem 3.1).
- For parameter regimes like $$\sigma^2 < r$$ (common in the examples) and/or fine spatial grids (large $$M$$), the iteration matrix can acquire negative eigenvalues near $$-1$$, causing slowly decaying oscillatory error components.

### Algorithm C — Proposed nonstandard semi-implicit step (Section 3.2)

The proposed scheme intentionally reduces formal time accuracy to $$O(\Delta t)$$ while maintaining $$O(\Delta S^2)$$ in space and enforcing **positivity** and **non-oscillatory damping** through matrix-structure constraints.

Key nonstandard ingredient: approximate the reaction term $$-rV$$ using a **bivariate (nonlocal in space / mixed in time-level) approximation**:
$$
V(S_j,t)\ \approx\ b\left(V_{j-1}^n + V_{j+1}^n\right) + (1-2b)V_j^{n+1},
$$
with local discretization error reported as $$O(\Delta S^2,\Delta t)$$ and constant parameter $$b$$ chosen below.

Discretize:
- $$\frac{\partial V}{\partial S}$$ by an explicit centered difference (time level $$n$$).
- $$\frac{\partial^2 V}{\partial S^2}$$ by an implicit scheme (time level $$n+1$$).

The resulting semi-implicit linear system:
$$
P\,V^{n+1} = N\,V^n,
$$
where:

$$
P = \operatorname{tridiag}\left\{
-\frac{\Delta t}{2}\left(\frac{\sigma S_j}{\Delta S}\right)^2\ ;\
1 + \Delta t\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 + r(1-2b)\right]\ ;\
-\frac{\Delta t}{2}\left(\frac{\sigma S_j}{\Delta S}\right)^2
\right\},
$$

$$
N = \operatorname{tridiag}\left\{
\frac{\Delta t\,r}{2}\left(-\frac{S_j}{\Delta S}-2b\right)\ ;\
1\ ;\
\frac{\Delta t\,r}{2}\left(\frac{S_j}{\Delta S}-2b\right)
\right\}.
$$

**Choice of $$b$$ (paper’s design point):** set
$$
b = -\frac{M}{2}.
$$

This yields (as asserted in the paper) nonnegative entries in $$N$$ and makes $$P$$ an irreducible diagonally dominant **M-matrix**, implying componentwise positivity of $$P^{-1}$$ and therefore positivity of iterates when $$V^0 \ge 0$$.

> Ambiguity note (internal consistency): the provided text shows a minor inconsistency in the proof of Theorem 3.2 where $$b$$ is restated differently; the accompanying norm bounds and the stated diagonal-dominance condition $$\Delta t < 1/(rM)$$ match $$b=-M/2$$ and the explicit formula $$\|N\|_\infty = 1 + r\Delta t M$$, so $$b=-M/2$$ is the only value consistent with the stated properties.

**Per-step procedure:**

**Inputs:** $$V^n$$, $$r,\sigma$$, $$\Delta S,\Delta t$$, $$M$$.

**Output:** $$V^{n+1}$$.

1. Fix $$b=-M/2$$.
2. Assemble tridiagonal $$P$$ and $$N$$ as above.
3. Compute $$b_{\mathrm{rhs}} \leftarrow N V^n$$.
4. Solve $$P V^{n+1} = b_{\mathrm{rhs}}$$.

**Stability / monotonicity controls (paper’s statements):**
- Nonnegativity: $$N\ge 0$$ and $$P^{-1} > 0$$ imply $$V^{n+1}\ge 0$$ for $$V^n\ge 0$$.
- Norm contraction: $$\|V^{n+1}\|_\infty \le \frac{1+r\Delta t M}{1+r\Delta t(M+1)}\|V^n\|_\infty < \|V^n\|_\infty$$, yielding a discrete maximum principle and monotonicity.
- Eigenvalue positivity (Theorem 3.2): real, positive, distinct eigenvalues in $$(0,1)$$ under $$\Delta t < 1/(rM)$$.

## 4. Theoretical Results

### Theorem 3.1 (Crank–Nicolson; case $$\sigma^2 > r$$)

**Statement (restated in full with conditions):** Assume $$\sigma^2 > r$$ for the Black–Scholes parameters and consider the Crank–Nicolson matrix form $$P V^{n+1} = N V^n$$.

1. Under $$\sigma^2 > r$$, the inverse is componentwise positive:
   $$
   P^{-1} > 0,
   $$
   and the infinity norm satisfies the bound
   $$
   \|P^{-1}\|_\infty \le \frac{1}{\frac{1}{\Delta t}+\frac{r}{2}}.
   $$

2. If $$\sigma^2 > r$$ and
   $$
   \Delta t < \frac{2}{r+(\sigma M)^2},
   $$
   then:
   - (i) $$N\ge 0$$ and
     $$
     \|N\|_\infty = \frac{1}{\Delta t} - \frac{r}{2} > 0;
     $$
   - (ii) the scheme preserves **positivity** and satisfies a **discrete maximum principle** (in the paper’s sup-norm sense).

3. If $$\sigma^2 > r$$ and
   $$
   \Delta t < \frac{2}{r+2(\sigma M)^2},
   $$
   then the iteration matrix $$P^{-1}N$$ has $$M$$ distinct eigenvalues lying strictly in $$(0,1)$$:
   $$
   \lambda_i(P^{-1}N)\in(0,1)\quad\text{and are distinct}.
   $$

**Proof strategy sketch (2–4 sentences):**  
(1) The matrix $$P$$ is shown to be irreducible and row-diagonally dominant, which implies it is an **M-matrix**; standard M-matrix theory yields $$P^{-1} > 0$$ and allows bounding $$\|P^{-1}\|_\infty$$ via row dominance margins. (2) Under the stated $$\Delta t$$ restriction, all entries of $$N$$ become nonnegative, so $$V^{n+1}=P^{-1}NV^n$$ preserves positivity; the discrete maximum principle follows from bounding $$\|P^{-1}\|_\infty\|N\|_\infty \le 1$$. (3) Writing $$P=\frac{1}{\Delta t}I+C$$ and $$N=\frac{1}{\Delta t}I-C$$, the eigenvalues satisfy
$$
\lambda_i(P^{-1}N)=\frac{1-\Delta t\,\lambda_i(C)}{1+\Delta t\,\lambda_i(C)},
$$
and similarity of $$C$$ to a **Jacobi matrix** yields real, distinct eigenvalues; the time-step bound enforces $$\lambda_i(P^{-1}N)\in(0,1)$$.

**Complexity / bounds explicitly present in the paper:** norm bounds on $$\|P^{-1}\|_\infty$$ and $$\|N\|_\infty$$ as above; eigenvalue bounds via the rational transform; no explicit runtime complexity is stated beyond tridiagonal structure.

### Oscillation mechanism for $$\sigma^2 < r$$ (Crank–Nicolson; qualitative result)

The paper’s Section 3.1 also analyzes the regime $$\sigma^2 < r$$ (illustrated with $$\sigma=0.2$$, $$r=0.05$$ so $$\sigma^2=0.04<0.05$$) and argues that as spatial refinement increases (large $$M$$), the spectrum of an associated matrix $$C$$ expands so that some eigenvalues of $$P^{-1}N$$ approach $$-1$$. Decomposing the initial vector into eigenmodes then yields alternating-sign components $$\lambda_j^n v_j$$ that produce localized oscillations that decay slowly when $$\lambda_j \approx -1$$.

### Jacobi matrix definition (footnote-level definitional content)

A **Jacobi matrix** (in the sense used by the paper) is a real tridiagonal matrix $$A=\operatorname{tridiag}\{c_i,a_i,b_i\}$$ whose off-diagonal entries have the same sign, formalized as $$c_i b_{i-1} > 0$$ for the valid index range; any such matrix is diagonally similar to a symmetric Jacobi matrix, hence has real, distinct eigenvalues.

### Theorem 3.2 (semi-implicit nonstandard scheme eigenvalues)

**Statement (restated in full with conditions):** For the proposed semi-implicit scheme $$P V^{n+1} = N V^n$$ (Section 3.2) with the paper’s choice $$b=-M/2$$, if
$$
\Delta t < \frac{1}{rM},
$$
then the iteration matrix $$P^{-1}N$$ admits $$M$$ real, positive, distinct eigenvalues strictly inside the unit interval:
$$
\lambda_i(P^{-1}N)\in(0,1),\qquad i=1,\ldots,M,
$$
and all such eigenvalues are distinct.

**Proof strategy sketch (2–4 sentences):**  
The proof first uses diagonal dominance of $$N$$ under $$\Delta t < 1/(rM)$$ to show that $$N$$ is diagonally similar to a symmetric positive definite matrix, enabling a Rayleigh-quotient argument that any eigenvalue $$\lambda$$ solving $$P^{-1}Nu=\lambda u$$ must be real and positive. Distinctness is obtained by analyzing the generalized eigenvalue condition $$\det(N-\lambda P)=0$$, transforming $$N-\lambda P$$ via diagonal similarity into a symmetric **Jacobi matrix** whose eigenvalues are real and simple; the mapping from roots of these eigenvalue functions yields distinct solutions $$\lambda_1,\ldots,\lambda_M$$.

**Additional stability/monotonicity bounds stated (not part of the theorem label):**
- With $$b=-M/2$$, the paper bounds the spectral radius by a norm inequality:
  $$
  \rho(P^{-1}N) \le \|P^{-1}N\|_\infty
  \le
  \frac{1+r\Delta t M}{1+r\Delta t(M+1)} < 1 \quad \forall\,\Delta t,
  $$
  implying unconditional linear stability in this norm sense and convergence by Lax equivalence when combined with consistency.
- Local truncation error is reported as
  $$
  O(\Delta t) + O(\Delta S^2),
  $$
  with an additional accuracy constraint suggested for large $$M$$:
  $$
  \Delta t \ll \frac{1}{r(1+M)}.
  $$

## 5. Experimental Evaluation

### Experimental configuration summary (datasets, baselines, metrics, hyperparameters)

| Item | Fig. 1 / oscillation demo | Example 4.1 (Table 1 / Fig. 3) |
|---|---|---|
| Contract | Discrete double barrier knock-out call | Discrete double barrier knock-out call |
| Payoff | $$\max(S-K,0)\mathbf{1}_{[L,U]}(S)$$ | same |
| Barriers | $$L=90,\ U=110$$ | $$L=95,\ U=110$$ |
| Strike | $$K=100$$ | $$K=100$$ |
| Maturity | $$T=0.5$$ (six months) | $$T=0.5$$ (six months) |
| Rates | $$r=0.05$$ | $$r=0.05$$ |
| Volatility | $$\sigma=0.2$$ (so $$\sigma^2<r$$) | $$\sigma=0.25$$ |
| Monitoring | described as discrete; last monitoring at $$T$$ | “monitored 5 times” |
| Spatial truncation | not explicitly stated in the caption | $$S_{\max}=200$$ |
| Mesh sizes | $$\Delta S=0.01,\ \Delta t=0.001$$ (Crank–Nicolson shown oscillatory) | CN / fully implicit / Duffy: $$\Delta S=0.05,\ \Delta t=10^{-5}$$; semi-implicit: $$\Delta S=0.05,\ \Delta t=0.001$$ |
| Baselines compared | Crank–Nicolson (qualitative) | Standard fully implicit; Crank–Nicolson; Duffy exponentially fitted implicit; proposed semi-implicit; Monte Carlo |
| Evaluation metric(s) | Qualitative oscillations/negativity near barriers and strike | Pointwise option values at selected $$S_0$$; agreement with Monte Carlo (with reported standard error) |

### Key quantitative results (Table 1 reproduction)

Option prices for Example 4.1 (discrete double knock-out call; monitored 5 times; $$K=100,\sigma=0.25,T=0.5,r=0.05,L=95,U=110$$). Mesh and method settings as reported in the table caption: CN/fully implicit/Duffy use $$\Delta S=0.05,\Delta t=10^{-5}$$; semi-implicit uses $$\Delta S=0.05,\Delta t=0.001$$; $$S_{\max}=200$$.

| Underlying $$S_0$$ | Standard implicit | Crank–Nicolson | Duffy exp.-fitted implicit | Proposed semi-implicit | Monte Carlo (std. err.) |
|---:|---:|---:|---:|---:|---:|
| 95 | 0.17564 | 0.17561 | 0.17315 | 0.17398 | 0.17359 (0.00054) |
| 95.0001 | 0.17904 | 0.17963 | 0.17395 | 0.17412 | 0.17486 (0.00064) |
| 95.5 | 0.18322 | 0.18324 | 0.18109 | 0.18152 | 0.18291 (0.00066) |
| 99.5 | 0.22818 | 0.22813 | 0.22819 | 0.22902 | 0.22923 (0.00073) |
| 100 | 0.23123 | 0.23122 | 0.23137 | 0.23171 | 0.23263 (0.00036) |
| 100.5 | 0.23359 | 0.23361 | 0.23386 | 0.23246 | 0.23410 (0.00073) |
| 109.5 | 0.17582 | 0.17583 | 0.17323 | 0.17326 | 0.17426 (0.00063) |
| 109.9999 | 0.16982 | 0.16989 | 0.16656 | 0.16719 | 0.16732 (0.00062) |
| 110 | 0.16906 | 0.16912 | 0.16616 | 0.16703 | 0.16712 (0.00042) |

### Reported experimental conclusions (as stated, paraphrased)

- The Crank–Nicolson method’s nominal second-order accuracy $$O(\Delta S^2,\Delta t^2)$$ does not translate into a practical advantage under discontinuous/renewed initial data unless a very small $$\Delta t$$ is used (per the time-step restrictions in the spectral/positivity analysis).
- The proposed semi-implicit method, despite lower time accuracy $$O(\Delta S^2,\Delta t)$$, is reported to deliver values close to Monte Carlo while allowing a substantially larger $$\Delta t$$ than the Crank–Nicolson restrictions required to preserve positivity and avoid oscillations.
- The semi-implicit scheme is reported to work in both regimes $$\sigma^2<r$$ and $$\sigma^2>r$$ and to be free of the localized spikes near the strike and barriers that appear in Crank–Nicolson plots for discretely monitored barriers.

### Ablations and sensitivity axes

No explicit multi-axis ablation study is reported; the paper’s sensitivity discussion centers on:
- time-step restrictions in Theorem 3.1 vs Theorem 3.2,
- the impact of large $$M$$ on stability/oscillation behavior and on the additional accuracy constraint $$\Delta t \ll 1/(r(1+M))$$ for the semi-implicit method.

## 6. ASCII Architecture / Workflow Diagram(s)

### End-to-end workflow with monitoring updates (barrier projection)

```
┌──────────────────────────────────────────────────────────────────┐
│ Inputs                                                          │
│  r, σ, K, L, U, T; monitoring dates {t_i}; Smax, ΔS, Δt; scheme  │
└──────────────────────────────────────────────┬───────────────────┘
                                               ▼
┌──────────────────────────────────────────────────────────────────┐
│ Grid construction                                                │
│  S_j=jΔS, j=0..M (Smax=MΔS)                                      │
│  t_n=nΔt, n=0..X (T=XΔt)                                         │
└──────────────────────────────────────────────┬───────────────────┘
                                               ▼
┌──────────────────────────────────────────────────────────────────┐
│ Discontinuous initial condition (payoff)                         │
│  V_j^0 = (S_j-K)^+ · 1_[L,U](S_j)                                │
└──────────────────────────────────────────────┬───────────────────┘
                                               ▼
┌──────────────────────────────────────────────────────────────────┐
│ PDE stepping between monitoring times                            │
│  For each n: solve tridiagonal system                             │
│    P V^{n+1} = N V^n                                              │
│  (P,N = Crank–Nicolson OR proposed semi-implicit)                │
└──────────────────────────────┬───────────────────────────────────┘
                               │ if t_{n+1} ∈ {t_i}
                               ▼
┌──────────────────────────────────────────────────────────────────┐
│ Discrete knock-out update (renew discontinuities)                │
│  V^{n+1} ← V^{n+1} ⊙ 1_[L,U](S)                                  │
└──────────────────────────────┬───────────────────────────────────┘
                               ▼
┌──────────────────────────────────────────────────────────────────┐
│ Output                                                          │
│  Price curve {V_j^X}; extract/interpolate V(S0,T)                │
└──────────────────────────────────────────────────────────────────┘
```

### Single-step linear-algebra view (emphasis on positivity / spectrum)

```
┌──────────────────────────────┐
│ Known vector V^n             │
│ (after possible monitor proj)│
└───────────────┬──────────────┘
                ▼
┌──────────────────────────────┐
│ Form RHS: b = N V^n          │
│ Requirement for positivity:  │
│   N ≥ 0                      │
└───────────────┬──────────────┘
                ▼
┌──────────────────────────────┐
│ Solve tridiagonal system     │
│   P V^{n+1} = b              │
│ Requirement for positivity:  │
│   P is M-matrix ⇒ P^{-1} > 0 │
└───────────────┬──────────────┘
                ▼
┌──────────────────────────────┐
│ New vector V^{n+1}           │
│ Properties sought:           │
│  V^{n+1} ≥ 0; ρ(P^{-1}N)<1;  │
│  eigenvalues in (0,1)        │
└──────────────────────────────┘
```

## 7. Follow-Up Works & Extensions

### 7.1 Direct extensions by the same authors (positivity, damping, transform methods)

Milev and Tagliani’s *Serdica Mathematical Journal* paper on low-volatility regimes surveys and analyzes numerical diffusion introduced by nonstandard schemes for Black–Scholes, focusing on discontinuous payoffs and the trade-off between oscillation control and artificial diffusion; it positions nonstandard/modified schemes as practical tools when discontinuities make higher-order centered schemes unreliable. [Milev & Tagliani, Serdica Math. J. 2010] ([eudml.org](https://eudml.org/doc/281441?utm_source=openai))

Milev and Tagliani propose a modification of the fully implicit scheme (used in Rannacher-type damping strategies) to guarantee smooth, positivity-preserving solutions for discontinuous payoffs and low-volatility settings; the work explicitly connects spurious oscillations to positivity failure and to non-positive spectra of iteration matrices, paralleling the Serdica 2010 analysis but shifting the focus from Crank–Nicolson to the fully implicit baseline used for initial damping. [Milev & Tagliani, J. Comput. Appl. Math. 2013] ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042712004128))

Milev and Tagliani develop a mixed Laplace-transform plus finite-difference method (using Post–Widder inversion) for discretely monitored barrier options and prove an equivalence between this mixed method and certain finite-difference formulations; the method is presented as positivity-preserving and compatible with a discrete maximum principle, offering a transform-based angle on the same discontinuity-driven oscillation pathology highlighted in the Serdica 2010 paper. [Milev & Tagliani, Appl. Math. Comput. 2013] ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313007613?utm_source=openai))

Gzyl, Milev, and Tagliani apply Mellin-transform techniques and maximum-entropy inversion to price Black–Scholes options with discontinuous payoffs and time-dependent parameters, explicitly framing transform inversion as a route to bypass oscillatory artifacts that can arise in finite-difference solvers near discontinuities; the paper includes a discretely monitored barrier example and emphasizes accuracy/efficiency trade-offs relative to PDE discretizations. [Gzyl et al., Finance Research Letters 2017] ([iris.unitn.it](https://iris.unitn.it/handle/11572/190747?utm_source=openai))

### 7.2 Later NSFD / positivity-preserving schemes that explicitly cite this Serdica 2010 approach

Mehdizadeh Khalsaraei et al. review a range of finite-difference schemes for Black–Scholes and propose an “improved” nonstandard method (built around nonlocal discretizations) with sufficient conditions for positivity and claims of being free of oscillations under discontinuous initial data; their literature discussion explicitly cites the Serdica 2010 semi-implicit idea as an example of sacrificing higher-order time accuracy to reduce time-step restrictions while maintaining positivity. [Mehdizadeh Khalsaraei et al., Mathematics 2022] ([mdpi.com](https://www.mdpi.com/2227-7390/10/11/1846))

### 7.3 Independent but thematically aligned work on positivity / maximum principle for option-pricing PDEs

Chernogorova and Valkov present a positive splitting (locally one-dimensional) method for a two-dimensional Black–Scholes PDE in the Hull–White stochastic-volatility setting, proving a discrete maximum principle and hence nonnegativity preservation; the connection to the Serdica 2010 paper is methodological (maximum-principle/positivity as a design constraint) rather than contract-specific (their focus is multidimensional PDE structure rather than barrier-induced discontinuities). [Chernogorova & Valkov, arXiv 2013] ([arxiv.org](https://arxiv.org/abs/1307.0232?utm_source=openai))

Mehdizadeh Khalsaraei et al. propose a nonstandard finite-difference method for a generalized Black–Scholes equation using a nonlocal approximation of the reaction term and an implicit time step, proving positivity preservation and stability and demonstrating non-oscillatory behavior in numerical experiments; the relation is conceptual (nonlocal reaction discretization to enforce qualitative properties) and complements the Serdica 2010 choice of a bivariate/nonlocal treatment of the $$-rV$$ term. [Mehdizadeh Khalsaraei et al., Symmetry 2022] ([mdpi.com](https://www.mdpi.com/2073-8994/14/1/141))

Golbabai, Ahmadian, and Milev use meshless radial basis functions (RBF) with a Crank–Nicolson time discretization to solve an American option model with jump diffusion (a PIDE with a free boundary); the paper cites the Serdica 2010 line of work in its references and illustrates how alternative spatial discretizations (RBF rather than FD) are explored in the same broad agenda of stable, accurate derivative valuation under challenging PDE/PIDE features. [Golbabai et al., Math. Comput. Modelling 2012] ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0895717711006066))

## 8. Industrial & Real-World Applications

QuantLib provides production-grade open-source implementations of PDE/finite-difference engines for option pricing, including finite-difference barrier option engines such as `FdBlackScholesBarrierEngine`, with configurable time/space grids and damping steps (a practical mechanism often used to mitigate oscillations from nonsmooth payoffs); this directly overlaps with the paper’s practical concerns (oscillations/positivity) even though QuantLib does not publicly document implementing the specific Serdica 2010 semi-implicit discretization. [GitHub: lballabio/QuantLib] ([github.com](https://github.com/lballabio/QuantLib?utm_source=openai))

QuantLib additionally distributes widely used language bindings and binary packaging (e.g., Python wheels via `pip install QuantLib`), enabling deployment of finite-difference pricing engines in research and applied environments; public deployment scale is not quantified on the project pages, but the project explicitly positions itself for “real-life” modeling, trading, and risk management use. [QuantLib official site / downloads] ([quantlib.org](https://www.quantlib.org/?utm_source=openai))

The finmath library (`finmath-lib`) is a large Apache-licensed Java/JVM quantitative-finance library with finite-difference support (including theta-scheme infrastructure for Black–Scholes-type PDEs) and extensive model/product tooling; it represents an open-source environment where positivity/monotonicity constraints and stability considerations for PDE discretizations can be implemented and stress-tested, although it does not claim an implementation of the specific Serdica 2010 semi-implicit scheme. [GitHub: finmath/finmath-lib] ([github.com](https://github.com/finmath/finmath-lib?utm_source=openai))

## 9. Consolidated Reference List

[1] Mariyan Milev, Aldo Tagliani. “Low Volatility Options and Numerical Diffusion of Finite Difference Schemes.” *Serdica Mathematical Journal*, 2010. `https://eudml.org/doc/281441` ([eudml.org](https://eudml.org/doc/281441?utm_source=openai))

[2] Mariyan Milev, Aldo Tagliani. “Efficient implicit scheme with positivity preserving and smoothing properties.” *Journal of Computational and Applied Mathematics*, 243 (2013), 1–9. DOI: 10.1016/j.cam.2012.09.039 ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042712004128))

[3] Aldo Tagliani, Mariyan Milev. “Laplace Transform and finite difference methods for the Black–Scholes equation.” *Applied Mathematics and Computation*, 220 (2013), 649–658. DOI: 10.1016/j.amc.2013.07.011 ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313007613?utm_source=openai))

[4] H. Gzyl, M. Milev, A. Tagliani. “Discontinuous payoff option pricing by Mellin transform: A probabilistic approach.” *Finance Research Letters*, 20 (2017), 281–288. DOI: 10.1016/j.frl.2016.10.011 ([iris.unitn.it](https://iris.unitn.it/handle/11572/190747?utm_source=openai))

[5] Mohammad Mehdizadeh Khalsaraei, Ali Shokri, Higinio Ramos, Zahra Mohammadnia, Pari Khakzad. “A Positivity-Preserving Improved Nonstandard Finite Difference Method to Solve the Black-Scholes Equation.” *Mathematics* (MDPI), 10(11) (2022), 1846. DOI: 10.3390/math10111846 ([mdpi.com](https://www.mdpi.com/2227-7390/10/11/1846))

[6] Mohammad Mehdizadeh Khalsaraei, Mohammad Mehdi Rashidi, Ali Shokri, Higinio Ramos, Pari Khakzad. “A Nonstandard Finite Difference Method for a Generalized Black–Scholes Equation.” *Symmetry* (MDPI), 14(1) (2022), 141. DOI: 10.3390/sym14010141 ([mdpi.com](https://www.mdpi.com/2073-8994/14/1/141))

[7] T. Chernogorova, R. Valkov. “Positive Splitting Method for the Hull & White 2D Black-Scholes Equation.” arXiv preprint, 2013. arXiv:1307.0232 ([arxiv.org](https://arxiv.org/abs/1307.0232?utm_source=openai))

[8] Ahmad Golbabai, Davood Ahmadian, Mariyan Milev. “Radial basis functions with application to finance: American put option under jump diffusion.” *Mathematical and Computer Modelling*, 55(3–4) (2012), 1354–1362. DOI: 10.1016/j.mcm.2011.10.014 ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0895717711006066))

[9] Luigi Ballabio (maintainer). “QuantLib: The QuantLib C++ library.” GitHub repository. `https://github.com/lballabio/QuantLib` ([github.com](https://github.com/lballabio/QuantLib?utm_source=openai))

[10] QuantLib annotated source. “FdBlackScholesBarrierEngine Class Reference.” Documentation page. `https://rkapl123.github.io/QLAnnotatedSource/d8/d27/class_quant_lib_1_1_fd_black_scholes_barrier_engine.html` ([rkapl123.github.io](https://rkapl123.github.io/QLAnnotatedSource/d8/d27/class_quant_lib_1_1_fd_black_scholes_barrier_engine.html?utm_source=openai))

[11] QuantLib source browser. “fdblackscholesbarrierengine.cpp source code.” `https://codebrowser.dev/quantlib/quantlib/ql/pricingengines/barrier/fdblackscholesbarrierengine.cpp.html` ([codebrowser.dev](https://codebrowser.dev/quantlib/quantlib/ql/pricingengines/barrier/fdblackscholesbarrierengine.cpp.html?utm_source=openai))

[12] QuantLib. “Download Page (Python wheels, C# packages, source distributions).” `https://www.quantlib.org/download.shtml` ([quantlib.org](https://www.quantlib.org/download.shtml?utm_source=openai))

[13] Christian Fries et al. “finmath-lib: Mathematical Finance Library.” GitHub repository. `https://github.com/finmath/finmath-lib` ([github.com](https://github.com/finmath/finmath-lib?utm_source=openai))

---
Learn more:
1. [Volume 36, Number 1, 2010](https://www.math.bas.bg/serdica/n1_10.html?utm_source=openai)
2. [EUDML  |  Nonstandard Finite Difference Schemes with Application to Finance: Option Pricing](https://eudml.org/doc/281413?utm_source=openai)
3. [EUDML  |  Low Volatility Options and Numerical Diffusion of Finite Difference Schemes](https://eudml.org/doc/281441?utm_source=openai)
4. [Efficient implicit scheme with positivity preserving and smoothing properties - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0377042712004128)
5. [Laplace Transform and finite difference methods for the Black–Scholes equation - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0096300313007613?utm_source=openai)
6. [Discontinuous payoff option pricing by Mellin transform: A probabilistic approach](https://iris.unitn.it/handle/11572/190747?utm_source=openai)
7. [A Positivity-Preserving Improved Nonstandard Finite Difference Method to Solve the Black-Scholes Equation](https://www.mdpi.com/2227-7390/10/11/1846)
8. [Positive Splitting Method for the Hull & White 2D Black-Scholes Equation](https://arxiv.org/abs/1307.0232?utm_source=openai)
9. [A Nonstandard Finite Difference Method for a Generalized Black–Scholes Equation](https://www.mdpi.com/2073-8994/14/1/141)
10. [Radial basis functions with application to finance: American put option under jump diffusion - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0895717711006066)
11. [GitHub - lballabio/QuantLib: The QuantLib C++ library](https://github.com/lballabio/QuantLib?utm_source=openai)
12. [QuantLib, a free/open-source library for quantitative finance](https://www.quantlib.org/?utm_source=openai)
13. [GitHub - finmath/finmath-lib: Mathematical Finance Library: Algorithms and methodologies related to mathematical finance.](https://github.com/finmath/finmath-lib?utm_source=openai)
14. [QuantLib: FdBlackScholesBarrierEngine Class Reference](https://rkapl123.github.io/QLAnnotatedSource/d8/d27/class_quant_lib_1_1_fd_black_scholes_barrier_engine.html?utm_source=openai)
15. [fdblackscholesbarrierengine.cpp source code \[quantlib/ql/pricingengines/barrier/fdblackscholesbarrierengine.cpp\] - Codebrowser](https://codebrowser.dev/quantlib/quantlib/ql/pricingengines/barrier/fdblackscholesbarrierengine.cpp.html?utm_source=openai)
16. [QuantLib Download Page](https://www.quantlib.org/download.shtml?utm_source=openai)

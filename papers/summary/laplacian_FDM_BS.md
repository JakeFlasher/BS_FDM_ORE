## 1. Paper Identity

- **Title:** *Laplace Transform and finite difference methods for the Black–Scholes equation* ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313007613))  
- **Author(s):** Aldo Tagliani; Mariyan Milev (corresponding author) ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313007613))  
- **Affiliation(s):**
  - Department of Computer and Management Sciences, University of Trento, Trento, Italy ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313007613))  
  - Department of Informatics and Statistics, University of Food Technologies, Plovdiv, Bulgaria ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313007613))  
- **Venue:** *Applied Mathematics and Computation* 220 (2013) 649–658 ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313007613))  
- **Year:** 2013 ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313007613))  
- **DOI:** 10.1016/j.amc.2013.07.011 ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313007613))  
- **arXiv ID:** None reported in the source metadata. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313007613))  

## 2. Problem Statement & Formulation

Pricing a **discretely monitored double-barrier knock-out call option** under a one-factor Black–Scholes model where the option value is repeatedly “truncated” at monitoring dates by a corridor condition $$S \in [L,U]$$ creates renewed discontinuities in the solution data. Low volatility makes the Black–Scholes PDE convection-dominated, and standard higher-order time-stepping (notably Crank–Nicolson) can exhibit persistent localized spurious oscillations near discontinuities (barriers and strike), violating financial meaning (non-negativity and maximum principle).

The paper assumes the volatility and interest rate depend only on the asset level:
- Volatility: $$\sigma = \sigma(S)$$
- Interest rate: $$r = r(S)$$

Let $$V(S,t)$$ denote the option value, where $$S \in \mathbb{R}_+$$ is asset price and $$t \in [0,T]$$ is **time to expiry** (so $$t=0$$ corresponds to maturity payoff). The governing PDE is the Black–Scholes equation written in “time-to-maturity” form:
$$
-\frac{\partial V}{\partial t}(S,t)
+ r(S)\,S\,\frac{\partial V}{\partial S}(S,t)
+ \frac{1}{2}\sigma^2(S)\,S^2\,\frac{\partial^2 V}{\partial S^2}(S,t)
- r(S)\,V(S,t)
= 0.
$$

For a double-barrier knock-out call with discrete monitoring dates
$$
0 = t_0 < t_1 < \cdots < t_F = T,
$$
the payoff at maturity and the monitoring updates are posed via an indicator:
$$
\mathbf{1}_{[L,U]}(S) =
\begin{cases}
1, & S \in [L,U],\\
0, & S \notin [L,U].
\end{cases}
$$

Initial (maturity) condition:
$$
V(S,0) = \max(S-K,0)\,\mathbf{1}_{[L,U]}(S).
$$

Boundary conditions (consistent with a truncated payoff):
$$
V(S,t) \to 0 \text{ as } S \to 0,
\qquad
V(S,t) \to 0 \text{ as } S \to \infty.
$$

Discrete monitoring update (knock-out at monitoring times):
$$
V(S,t_i) = V(S,t_i^-)\,\mathbf{1}_{[L,U]}(S),
\qquad i=1,\ldots,F.
$$

The continuous solution (between monitoring instants) is smooth for any $$t>0$$, and the parabolic PDE satisfies a maximum principle stated as:
$$
\max_{S \in [0,\infty)} |V(S,t_1)|
\ge
\max_{S \in [0,\infty)} |V(S,t_2)|
\quad \text{for } t_1 \le t_2,
$$
interpreted as the option-value maximum not increasing with increasing $$t$$ (in time-to-expiry parametrization).

## 3. Core Methodology

### 3.1 Mixed Laplace-transform / finite-difference approach (“mixed method”)

The paper’s **mixed method** applies a Laplace transform in the time variable $$t$$, converts the PDE into a boundary-value ODE in $$S$$ (for each Laplace parameter value), solves that ODE with a spatial finite-difference scheme designed to preserve qualitative properties (positivity / discrete maximum principle), and then inverts the Laplace transform on the real axis via the Post–Widder formula.

#### Time Laplace transform and derivative family

Define the Laplace transform (in $$t$$) of $$V$$:
$$
U(S,\lambda)
= \mathcal{L}[V(S,\cdot)](\lambda)
= \int_0^\infty e^{-\lambda t}\,V(S,t)\,dt,
\qquad \lambda > \lambda_0.
$$

Define the $$k$$-th derivative (with respect to $$\lambda$$):
$$
U^{(k)}(S,\lambda)
:= \frac{d^k}{d\lambda^k}U(S,\lambda)
= (-1)^k \int_0^\infty t^k e^{-\lambda t}\,V(S,t)\,dt.
$$

From $$V(S,t) \ge 0$$ (financial meaning), $$U(S,\lambda)$$ is **completely monotonic** in $$\lambda$$, i.e.
$$
(-1)^k\,U^{(k)}(S,\lambda) \ge 0
\quad \text{for all integers } k \ge 0,
$$
and the numerical method is constructed to preserve this sign-alternation structure.

#### Transformed ODE recursion

Applying the Laplace transform to the PDE yields, for each integer $$k \ge 0$$, an ODE in $$S$$:
$$
-\frac{1}{2}\sigma^2(S)\,S^2\,\frac{d^2 U^{(k)}}{dS^2}(S,\lambda)
- r(S)\,S\,\frac{d U^{(k)}}{dS}(S,\lambda)
+ (r(S)+\lambda)\,U^{(k)}(S,\lambda)
=
\begin{cases}
V(S,0), & k=0,\\
-k\,U^{(k-1)}(S,\lambda), & k=1,2,\ldots
\end{cases}
$$

Boundary conditions are inherited from $$V$$:
$$
U^{(k)}(0,\lambda) = 0,
\qquad
U^{(k)}(\infty,\lambda) = 0,
\qquad k \ge 0.
$$

This is a recursion: solving for $$U^{(0)}$$ gives the base, then each $$U^{(k)}$$ is obtained from $$U^{(k-1)}$$.

### 3.2 Spatial discretization and the tridiagonal operator $$A_{\mathrm{PW}}$$

For computation, truncate $$S \in [0,\infty)$$ to $$[0,S_{\max}]$$ and enforce:
$$
U^{(k)}(0,\lambda)=0,
\qquad
U^{(k)}(S_{\max},\lambda)=0.
$$

Uniform grid:
$$
S_j = j\,\Delta S,
\qquad
j = 0,1,\ldots,M,
\qquad
S_{\max} = M\,\Delta S.
$$

Second derivative discretization (always centered):
$$
\frac{d^2 u}{dS^2}(S_j)
\approx
\frac{u_{j-1} - 2u_j + u_{j+1}}{\Delta S^2}.
$$

First derivative discretization is chosen based on convection dominance, using the paper’s rule:

- If $$\sigma^2 > r$$: use centered difference
  $$
  \frac{d u}{dS}(S_j)
  \approx
  \frac{u_{j+1}-u_{j-1}}{2\Delta S}.
  $$
- If $$\sigma^2 < r$$: use forward (upwind-type) difference
  $$
  \frac{d u}{dS}(S_j)
  \approx
  \frac{u_{j+1}-u_j}{\Delta S}.
  $$

This yields a tridiagonal linear system for each $$k$$:
$$
A_{\mathrm{PW}}(\lambda)\,U^{(0)} = V_0,
\qquad
A_{\mathrm{PW}}(\lambda)\,U^{(k)} = -k\,U^{(k-1)},
\quad k=1,2,\ldots,N,
$$
where $$V_0$$ is the vector $$\big(V(S_j,0)\big)_{j=1}^{M-1}$$ and similarly for $$U^{(k)}$$ (interior nodes).

**Case 1: $$\sigma^2 > r$$ (centered convection).** For each interior index $$j$$, define
$$
a_j = -\frac{1}{2}\left[\left(\frac{\sigma(S_j)\,S_j}{\Delta S}\right)^2 - r(S_j)\frac{S_j}{\Delta S}\right],
$$
$$
c_j = (r(S_j)+\lambda) + \left(\frac{\sigma(S_j)\,S_j}{\Delta S}\right)^2,
$$
$$
b_j = -\frac{1}{2}\left[\left(\frac{\sigma(S_j)\,S_j}{\Delta S}\right)^2 + r(S_j)\frac{S_j}{\Delta S}\right].
$$
Then $$A_{\mathrm{PW}} = \mathrm{tridiag}(a_j,c_j,b_j)$$.

**Case 2: $$\sigma^2 < r$$ (forward convection).**
$$
a_j = -\frac{1}{2}\left(\frac{\sigma(S_j)\,S_j}{\Delta S}\right)^2,
$$
$$
c_j = (r(S_j)+\lambda) + \left(\frac{\sigma(S_j)\,S_j}{\Delta S}\right)^2 + r(S_j)\frac{S_j}{\Delta S},
$$
$$
b_j = -\frac{1}{2}\left[\left(\frac{\sigma(S_j)\,S_j}{\Delta S}\right)^2 + 2\,r(S_j)\frac{S_j}{\Delta S}\right].
$$
Then $$A_{\mathrm{PW}} = \mathrm{tridiag}(a_j,c_j,b_j)$$.

With these sign patterns (negative off-diagonals, positive diagonal dominance under the stated regime choice), the paper emphasizes that $$A_{\mathrm{PW}}$$ is an **M-matrix**, implying:
$$
A_{\mathrm{PW}}^{-1} \ge 0 \quad \text{(componentwise)}.
$$
This is the discrete mechanism enforcing positivity preservation and avoiding spurious oscillations (when the M-matrix condition is maintained).

### 3.3 Laplace inversion by Post–Widder and the discrete approximation $$V_N(S,T)$$

Fix the time of interest $$T$$ (the end of a monitoring interval). The Post–Widder approximation of $$V(S,T)$$ is:
$$
V_N(S,T)
=
\frac{(-1)^N}{N!}
\left(\frac{N}{T}\right)^{N+1}
U^{(N)}\left(S,\frac{N}{T}\right),
$$
and the paper works with the discrete version at nodes $$S_j$$:
$$
V_N(S_j,T)
=
\frac{(-1)^N}{N!}
\left(\frac{N}{T}\right)^{N+1}
U^{(N)}\left(S_j,\frac{N}{T}\right).
$$

Implementation detail emphasized in Section 3: for fixed $$\lambda = \frac{N}{T}$$, $$A_{\mathrm{PW}}(\lambda)$$ is independent of the derivative index $$k$$ in the recursion, so one LU factorization of the tridiagonal matrix can be reused across all solves for $$U^{(k)}$$.

Unrolling the recursion and inserting into the Post–Widder expression yields the explicit matrix-iteration form:
$$
V_N(T)
=
\left(\frac{N}{T}\,A_{\mathrm{PW}}^{-1}\right)^{N+1} V_0,
$$
so the **iteration matrix** is:
$$
B := \frac{N}{T}\,A_{\mathrm{PW}}^{-1}.
$$

### 3.4 Monitoring-date update loop (discrete barrier enforcement)

For the full discretely monitored barrier option with dates $$t_0,\ldots,t_F$$, the computation over each interval can be treated as a repeated “propagate then truncate” operation. For interval length $$\Delta T_i := t_i - t_{i-1}$$:

1. Start from the interval’s initial condition $$V^{(i-1)}(S,0)$$ on the grid.
2. Compute $$V^{(i)}(S,\Delta T_i^-)$$ by applying the mixed method with $$T=\Delta T_i$$.
3. Enforce discrete monitoring truncation:
   $$
   V^{(i)}(S,0) \leftarrow V^{(i)}(S,\Delta T_i^-)\,\mathbf{1}_{[L,U]}(S).
   $$

The paper’s numerical illustrations focus on a single interval (interpreting “just before the first monitoring date” as the terminal computation for that subproblem).

### 3.5 Standard fully implicit finite-difference scheme (baseline and equivalence target)

Discretize time with step $$\Delta t$$ and define time levels $$t_n = n\Delta t$$. A fully implicit (backward Euler in $$t$$) discretization of the PDE yields a linear system:
$$
A_{\mathrm{FD}}\,V^{(n)} = V^{(n-1)},
\qquad n=1,2,\ldots
$$

For the same spatial grid, the paper gives (by regime) the tridiagonal form:

**If $$\sigma^2 > r$$ (centered convection):**
$$
A_{\mathrm{FD}}
=
\mathrm{tridiag}\left(
-\frac{\Delta t}{2}\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 - r\frac{S_j}{\Delta S}\right],
\;
1+\Delta t\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 + r\right],
\;
-\frac{\Delta t}{2}\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 + r\frac{S_j}{\Delta S}\right]
\right).
$$

**If $$\sigma^2 < r$$ (forward convection):**
$$
A_{\mathrm{FD}}
=
\mathrm{tridiag}\left(
-\frac{\Delta t}{2}\left(\frac{\sigma S_j}{\Delta S}\right)^2,
\;
1+\Delta t\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 + r\frac{S_j}{\Delta S} + r\right],
\;
-\frac{\Delta t}{2}\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 + 2r\frac{S_j}{\Delta S}\right]
\right).
$$

After $$N+1$$ time steps, the discrete solution at $$t=T$$ is:
$$
V^{(N+1)}(T)
=
\left(A_{\mathrm{FD}}^{-1}\right)^{N+1} V_0.
$$

### 3.6 ASCII pipeline diagram (end-to-end)

```text
┌──────────────────────────────────────────────────────────────────────────┐
│ Interval i (length ΔT = t_i - t_{i-1})                                    │
├──────────────────────────────────────────────────────────────────────────┤
│ Input: V_init(S) = V(S,0) on grid; parameters (σ(S), r(S), L,U,K, S_max)  │
└───────────────┬──────────────────────────────────────────────────────────┘
                │  Laplace in time t  (λ := N/ΔT)
                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Transformed recursion in S: for k = 0..N                                   │
│   A_PW(λ) U^(0) = V_init                                                   │
│   A_PW(λ) U^(k) = -k U^(k-1)                                               │
│ Spatial FD: centered U_SS; centered or forward U_S to keep M-matrix        │
└───────────────┬──────────────────────────────────────────────────────────┘
                │  Reuse LU(A_PW) across k
                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Post–Widder inversion at ΔT:                                               │
│   V_N(S,ΔT) = (-1)^N / N! * (N/ΔT)^(N+1) * U^(N)(S, N/ΔT)                  │
│ Equivalent closed form: V_N = ((N/ΔT) A_PW^{-1})^(N+1) V_init              │
└───────────────┬──────────────────────────────────────────────────────────┘
                │  Discrete monitoring update (knock-out)
                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ V_init(S) ← V_N(S,ΔT^-) * 1_[L,U](S)                                       │
│ Output at monitoring date t_i: V(S,t_i^-)                                  │
└──────────────────────────────────────────────────────────────────────────┘
```

## 4. Theoretical Results

### 4.1 Theorem (Post–Widder inversion formula)

**Statement (as used by the paper).** Let $$f : \mathbb{R}_+ \to \mathbb{R}$$ be continuous and of exponential order, so its Laplace transform
$$
F(\lambda) = \int_0^\infty e^{-\lambda t} f(t)\,dt
$$
converges for all $$\lambda > \lambda_0$$. Then for each continuity point $$t>0$$,
$$
f(t)
=
\lim_{k\to\infty}
\frac{(-1)^k}{k!}
\left(\frac{k}{t}\right)^{k+1}
F^{(k)}\left(\frac{k}{t}\right),
$$
where $$F^{(k)}$$ denotes the $$k$$-th derivative with respect to $$\lambda$$.

**Proof strategy sketch (paper-level).** The Post–Widder construction can be interpreted as a real-axis inversion that approximates the Dirac delta kernel by a sequence derived from Laplace-space derivatives. The growth condition ensures differentiability under the integral and legitimizes exchanging differentiation and integration. The limit extracts $$f(t)$$ at continuity points via concentration of the associated kernel around $$t$$.

### 4.2 Result (Equivalence between the mixed method and implicit time-stepping)

**Claim (Section 4).** Fix $$T$$ and choose the Post–Widder parameterization
$$
\lambda = \frac{N}{T}
\qquad \text{and} \qquad
\Delta t = \frac{T}{N}
\quad \text{so that} \quad
\lambda = \frac{1}{\Delta t}.
$$
With identical spatial discretization choices (centered vs forward convection) in both constructions, the Post–Widder mixed-method approximation equals the fully implicit finite-difference solution after $$N+1$$ time steps:
$$
V_N(S_j,T) = V^{(N+1)}(S_j,T).
$$

**Proof strategy sketch.** The spatial discretization of the transformed ODE produces the linear operator $$A_{\mathrm{PW}}(\lambda)$$, while the fully implicit time-stepping produces $$A_{\mathrm{FD}}(\Delta t)$$. The identification $$\lambda = 1/\Delta t$$ yields the scaling relation
$$
A_{\mathrm{PW}}(\lambda) = \frac{1}{\Delta t}\,A_{\mathrm{FD}}(\Delta t),
$$
and the Post–Widder recursion collapses to the same matrix-power propagation as repeated implicit Euler solves.

### 4.3 Result (Positivity preservation and discrete maximum principle via M-matrix structure)

**Claim (Sections 3–5).** Under the discretization rule that enforces $$A_{\mathrm{PW}}$$ to be an M-matrix (negative off-diagonals and appropriate diagonal dominance), the scheme preserves:

1. **Positivity:** if $$V_0 \ge 0$$, then for any computed terminal solution $$V_N(T) \ge 0$$ (componentwise).
2. **Discrete maximum principle (time-to-expiry monotonicity):** for discrete times $$t_1 < t_2$$ aligned with the scheme,
   $$
   \|V_{N_2}(t_2)\|_\infty \le \|V_{N_1}(t_1)\|_\infty.
   $$

**Proof strategy sketch.** For an M-matrix, $$A_{\mathrm{PW}}^{-1}\ge 0$$, so each recursion solve maps a nonnegative right-hand side to a nonnegative solution, matching complete monotonicity requirements $$(-1)^k U^{(k)} \ge 0$$. The iteration matrix $$B = (N/T)A_{\mathrm{PW}}^{-1}$$ satisfies an $$\ell_\infty$$ norm bound derived from diagonal dominance, giving a contraction factor less than 1 and yielding monotone decay of the sup norm across steps.

### 4.4 Result (Convergence of the mixed method)

**Claim (Section 4).** Because the mixed method is equivalent (under parameter identification) to a fully implicit finite-difference scheme, and the fully implicit scheme is consistent and unconditionally stable, the mixed-method approximation converges to the exact PDE solution as:
$$
\Delta S \to 0
\quad \text{and} \quad
\Delta t \to 0
\qquad
(\text{equivalently } M \to \infty,\, N \to \infty).
$$

**Proof strategy sketch.** The equivalence reduces convergence to that of the backward Euler time discretization coupled with the chosen spatial discretization. Standard parabolic PDE theory plus discrete stability yields convergence in the usual finite-difference sense. The Post–Widder sequence’s slow asymptotics affect efficiency but not the existence of the convergent limit when the discretization parameters shrink appropriately.

### 4.5 Analysis result (Why spurious oscillations appear if centered convection is used in convection-dominated regimes)

**Claim (Section 4, “unsuitable scheme” analysis).** If $$\sigma^2 < r$$ and one nevertheless discretizes the convection term by a centered difference (instead of the forward/upwind discretization), then in grid regions satisfying a stronger condition of the form
$$
\sigma^2 < \frac{r}{M},
$$
the subdiagonal entries of $$A_{\mathrm{PW}}$$ can become positive, so $$A_{\mathrm{PW}}$$ is no longer an M-matrix. The inverse $$A_{\mathrm{PW}}^{-1}$$ can contain negative entries, and the resulting $$V_N$$ may take negative values (financially meaningless), observed as spurious oscillations near discontinuities.

**Key spectral mechanism (Section 4).** The iteration matrix $$B = (N/T)A_{\mathrm{PW}}^{-1}$$ has eigenvalues whose real parts approach 1 as $$N \to \infty$$, and whose imaginary parts cluster near 0, so modes excited by discontinuities are not sufficiently damped. The paper rewrites the solution in eigenmodes:
$$
V_N(T)
=
\sum_{j=1}^{M} w_j \left(\frac{N}{T}\lambda_j(A_{\mathrm{PW}}^{-1})\right)^{N+1} v_j,
$$
and argues that when $$\left|\frac{N}{T}\lambda_j(A_{\mathrm{PW}}^{-1})\right| \approx 1$$, the corresponding contributions decay too slowly.

**Proof strategy sketch.** The sign change in off-diagonals breaks monotonicity, invalidating the M-matrix positivity guarantee. The eigenvalue bounds and norm bounds used show that as $$N$$ increases (which is required for convergence), the iteration operator becomes close to identity on the dominant eigenspaces, so discontinuity-driven high-frequency components persist and manifest as oscillations.

### 4.6 Numerical diffusion interpretation (forward convection discretization)

**Claim (Section 4 remarks and Section 5).** When $$\sigma^2 < r$$ and forward differencing is used to enforce an M-matrix, the method introduces numerical diffusion: the truncation error in the convection discretization yields an added diffusive contribution proportional to $$r S \Delta S$$, effectively increasing the diffusion coefficient from $$\frac{1}{2}\sigma^2 S^2$$ to
$$
\frac{1}{2}\sigma^2 S^2 + \frac{r S \Delta S}{2}.
$$
When $$\sigma^2$$ is very small, this artificial diffusion can dominate, explaining discrepancies observed in low-volatility plots.

**Proof strategy sketch.** Standard consistency analysis of one-sided differences expresses the leading truncation term as a second-derivative contribution. Interpreting this as an effective increase in the PDE’s diffusive coefficient explains “smearing” of the solution in low-volatility regimes while maintaining monotonicity and positivity.

## 5. Experimental Evaluation

### 5.1 Experimental setups (contracts, parameters, discretizations)

The paper reports numerical experiments on truncated call / discretely monitored barrier-style setups, primarily through plotted curves (Figures 1–5). No tabulated option values are provided in the body text for reproduction; parameter sets, grids, and qualitative outcomes are specified in figure captions and Section 5.

| Experiment (paper figure) | Contract / geometry | Parameters (as reported) | Methods compared | Metrics / outputs | Reported qualitative outcome |
|---|---|---|---|---|---|
| Fig. 1 | Truncated call value just before first monitoring date $$t_1 = T$$ | Upper: $$L=0, K=40, U=60, r=0.05, \sigma=0.1, T=1$$; grid $$\Delta S=0.2, \Delta t=0.001, S_{\max}=160$$ (via $$R=4$$). Lower: same but $$\sigma=0.001$$ and $$S_{\max}=80$$ (via $$R=2$$). | “Exact sol.” vs fully implicit FD | Plot of $$V(S,T)$$ vs $$S$$ | Low-volatility case shows materially larger discrepancy attributed to numerical diffusion; higher volatility case tracks “Exact sol.” more closely. |
| Fig. 2 | Same truncated call setup | $$L=0, K=40, U=60, r=0.05, \sigma=0.001, T=1$$; $$\Delta S=0.2, S_{\max}=80$$. FD: $$\Delta t=0.001$$ (so $$N=1000$$ if mapped). Mixed PW: $$N \in \{10,20,30\}$$. | Fully implicit FD vs mixed Post–Widder (“PW”) | $$V(S,T)$$ vs $$S$$ for different $$N$$ | Post–Widder sequence converges slowly; increasing $$N$$ improves agreement, consistent with known slow convergence of Post–Widder. |
| Fig. 3 | Truncated call with barriers far from strike (oscillation study) | $$L=90, K=100, U=110, r=0.05, \sigma=0.001, T=1, S_{\max}=200$$ (via $$R=2$$). Mixed method with centered convection under a regime where $$A_{\mathrm{PW}}$$ is not an M-matrix. Grid $$\Delta S=0.05$$; $$N=100,1000,10000$$. | Mixed method (different $$N$$) vs “Exact sol.” | $$V(S,T)$$ vs $$S$$; visible negativity | Spurious oscillations near discontinuities become more pronounced as $$N$$ increases, while convergence far from barriers improves. |
| Fig. 4 | Eigenvalue pattern diagnostic for $$A_{\mathrm{PW}}$$ | Illustration with parameterization $$r = \alpha M \sigma^2$$; $$\alpha \in \{0,1,1.2\}$$; $$\sigma=0.01, M=50$$; $$N/T=100$$ fixed. | Eigenvalue plots (real/complex) | Scatter of $$\mathrm{Re}(\lambda_j)$$ vs $$\mathrm{Im}(\lambda_j)$$ | Centered-convection discretization can yield complex eigenvalues depending on parameters; eigen-structure connected to oscillation persistence. |
| Fig. 5 | Extrapolation to accelerate Post–Widder | Geometry and parameters as Fig. 1 (upper). Post–Widder with $$N=1000$$ and $$2N$$; extrapolated solution $$V^{(E)}$$. | $$|V - V_N|$$ vs $$|V - V^{(E)}|$$ | Error curves vs $$S$$ | Extrapolation improves error only modestly; positivity is not generally preserved by extrapolation. |

### 5.2 Key quantitative results (what is and is not reproducible)

- Figures show comparisons labeled “Exact sol.”, “FD”, and “PW” but do not provide numerical tables of option values at specific grid points in the text excerpt, so exact numeric replication (e.g., L∞ errors) cannot be reconstructed without digitizing plots.
- Reported parameter values, grid sizes, and the qualitative behaviors (smearing from numerical diffusion; oscillations near discontinuities; very slow Post–Widder convergence; modest extrapolation gain) are fully specified in captions and Section 5 narrative.

### 5.3 Ablations / sensitivity axes explicitly varied

- Volatility magnitude: $$\sigma = 0.1$$ vs $$\sigma = 0.001$$ (Fig. 1).
- Post–Widder order (equivalently implicit time-step count): $$N \in \{10,20,30\}$$ (Fig. 2) and $$N \in \{100,1000,10000\}$$ (Fig. 3).
- Domain truncation: $$S_{\max}$$ changed via rule (5.2) with different $$R$$ (Fig. 1 and discussion in Section 5).
- Eigenvalue behavior vs parameter scaling $$r = \alpha M \sigma^2$$ (Fig. 4).
- Extrapolation using $$V_N$$ and $$V_{2N}$$ (Fig. 5).

## 6. ASCII Architecture / Workflow Diagram(s)

### 6.1 End-to-end discretely monitored barrier workflow (interval-wise propagation)

```text
┌────────────────────────────────────────────────────────────────────────────┐
│ Monitoring dates: 0=t0 < t1 < ... < tF=T                                    │
└────────────────────────────────────────────────────────────────────────────┘

For i = 1..F:

┌────────────────────────────────────────────────────────────────────────────┐
│ (1) Interval initial condition on grid: V_init(S)                           │
│     V_init(S) = payoff at i=1, else V_prev(S, t_{i-1}^-) * 1_[L,U](S)       │
└───────────────────────────────┬────────────────────────────────────────────┘
                                │  Mixed method (Laplace + FD + Post–Widder)
                                ▼
┌────────────────────────────────────────────────────────────────────────────┐
│ (2) Choose N, set λ := N/ΔT                                                 │
│ (3) Build tridiagonal A_PW(λ):                                               │
│     - centered U_SS                                                         │
│     - centered U_S if σ^2 > r; forward U_S if σ^2 < r                        │
│ (4) Solve recursion for U^(k), k=0..N with reused LU(A_PW)                   │
│ (5) Invert via Post–Widder to get V_N(S, ΔT)                                 │
└───────────────────────────────┬────────────────────────────────────────────┘
                                │  Apply monitoring (knock-out)
                                ▼
┌────────────────────────────────────────────────────────────────────────────┐
│ (6) V_prev(S, t_i^-) := V_N(S, ΔT^-)                                         │
│ (7) V_prev(S, t_i)  := V_prev(S, t_i^-) * 1_[L,U](S)                         │
└────────────────────────────────────────────────────────────────────────────┘
```

### 6.2 Equivalence map: Post–Widder order $$N$$ ↔ implicit time step $$\Delta t$$

```text
┌───────────────────────────────┐        ┌──────────────────────────────────┐
│ Mixed method parameterization │        │ Fully implicit FD time stepping  │
├───────────────────────────────┤        ├──────────────────────────────────┤
│ Choose integer N              │        │ Choose time step Δt              │
│ Set λ := N/T                  │  ⇔     │ Set N := T/Δt                    │
│ Solve A_PW(λ) recursions      │        │ Solve A_FD(Δt) V^n = V^{n-1}     │
│ V_N(T) = ((N/T)A_PW^{-1})^{N+1}V_0     │ V^{N+1}(T)=(A_FD^{-1})^{N+1}V_0  │
└───────────────────────────────┘        └──────────────────────────────────┘
with identification:  A_PW(λ) = (1/Δt) A_FD(Δt)  and  λ = 1/Δt = N/T
```

## 7. Follow-Up Works & Extensions

### 7.1 Hybrid Laplace-transform + spatial discretization schemes for option pricing

**Ma, Shi, Hon.** “Generalized Finite Integration Method with Laplace transform for European option pricing under Black–Scholes and Heston models.” *Engineering Analysis with Boundary Elements*, 2024. DOI: 10.1016/j.enganabound.2024.105751. [Ma et al., EABE 2024]. The paper applies a Laplace transform in time to remove time stepping and then solves the transformed spatial problem using a Generalized Finite Integration Method (GFIM), finally performing numerical inverse Laplace transform to recover terminal-time prices. The work positions Laplace-transform-based approaches as efficiency improvements when only terminal prices are needed, and it explicitly cites Tagliani–Milev as a prior Laplace+FD hybrid for Black–Scholes-type PDEs. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0955799724002169))

**Ma, Zhou, Cui.** “Hybrid Laplace transform and finite difference methods for pricing American options under complex models.” *Computers & Mathematics with Applications*, 2017. DOI: 10.1016/j.camwa.2017.04.018. [Ma et al., CAMWA 2017]. The work extends the Laplace-transform-then-discretize paradigm to free-boundary problems in American option pricing, proposing an iterative algorithm in Laplace space to determine the early exercise boundary and then invert to obtain prices and boundaries. The reported scope includes several underlying dynamics (CEV, jump–diffusion, regime switching, and fractional models), targeting exactly the class of situations where direct Laplace-space closed forms are unavailable. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0898122117302493))

### 7.2 Positivity/qualitative-structure preserving discretizations inspired by the “mixed method” viewpoint

**Khalsaraei, Shokri, Ramos, Mohammadnia, Khakzad.** “A Positivity-Preserving Improved Nonstandard Finite Difference Method to Solve the Black-Scholes Equation.” *Mathematics* 2022. DOI: 10.3390/math10111846. [Khalsaraei et al., Mathematics 2022]. This paper surveys multiple numerical schemes for the Black–Scholes PDE (explicitly including a Laplace-transform-based mixed method) and then proposes an improved nonstandard finite-difference approach aimed at eliminating spurious oscillations while preserving positivity. The numerical section includes truncated-call examples with barrier-style truncation, aligning with Tagliani–Milev’s emphasis that discontinuities and low volatility can break qualitative correctness unless monotone/M-matrix-like conditions are enforced. ([mdpi.com](https://www.mdpi.com/2227-7390/10/11/1846))

**Khalsaraei, Shokri, Wang, Bazm, Navidifar, Khakzad.** “Qualitatively Stable Schemes for the Black–Scholes Equation.” *Fractal and Fractional* 2023. DOI: 10.3390/fractalfract7020154. [Khalsaraei et al., Fractal Fract. 2023]. The work proposes a new scheme explicitly described as combining the Laplace transform method with a nonstandard finite-difference (NSFD) strategy, and reports qualitative properties (positivity, stability, consistency) in low-volatility regimes. The framing overlaps with the Tagliani–Milev paper’s central claim that qualitative financial constraints (non-negativity, maximum principle) are primary design targets when convection dominance and discontinuities make standard higher-order discretizations unreliable. ([doaj.org](https://doaj.org/article/d900141343ca4097a79d4e4505bdf2ff))

### 7.3 Fractional and modified Black–Scholes extensions that cite Tagliani–Milev as methodological precedent

**Sugandha, Rusyaman, Sukono, Carnia.** “Using a Mix of Finite Difference Methods and Fractional Differential Transformations to Solve Modified Black–Scholes Fractional Equations.” *Mathematics* 2024. DOI: 10.3390/math12071077. [Sugandha et al., Mathematics 2024]. This paper targets modified fractional Black–Scholes models and uses a hybrid numerical strategy (finite differences combined with a fractional differential transformation method), citing Tagliani–Milev as a precedent for mixing transform-domain techniques with finite differences in option-pricing PDEs. The connection is methodological rather than contract-specific: both works motivate hybridization as a route to stable approximations when standard discretizations are challenged by model stiffness or nonlocal operators. ([mdpi.com](https://www.mdpi.com/2227-7390/12/7/1077))

### 7.4 Variable-parameter and physics-inspired Black–Scholes generalizations citing the Laplace/FD hybrid literature

**El-Nabulsi, Anukool.** “Black–scholes equation in quantitative finance with variable parameters: a path to a generalized schrodinger equation.” *Financial Innovation* 2026 (published 14 January 2026). DOI: 10.1186/s40854-025-00877-7. [El-Nabulsi & Anukool, Financial Innovation 2026]. The article studies generalized Black–Scholes equations with time-varying parameters, using transformations and solving a diffusion equation via Laplace transform, and explicitly references Milev–Tagliani (2013) in its bibliography. The work is an extension along the “transform-based solution” axis (Laplace-domain analysis), but it shifts emphasis away from barrier discontinuities and monotone discretization toward analytical transformations and an analogy to generalized Schrödinger equations. ([link.springer.com](https://link.springer.com/article/10.1186/s40854-025-00877-7))

## 8. Industrial & Real-World Applications

**QuantLib (open-source quant finance library).** [GitHub: lballabio/QuantLib]. QuantLib is a large C++ quantitative-finance library described by its maintainers as a framework for “modeling, trading, and risk management in real-life,” with broad multi-language ecosystem support via wrappers and add-ins. The paper’s main equivalence result gives a direct interpretive bridge between the Post–Widder/Laplace-domain construction and standard fully implicit finite-difference time stepping, which is a numerics pattern that aligns with how production-grade libraries structure PDE-based pricers (even when they do not expose Laplace inversion directly). Scale: the project is large and widely visible publicly, but specific production deployment numbers are not reported on the cited pages. ([github.com](https://github.com/lballabio/QuantLib?utm_source=openai))

**Open Source Risk Engine (ORE) (open-source risk analytics framework).** [GitHub: OpenSourceRisk/Engine]. ORE provides an open-source framework for pricing and risk analysis, including contemporary risk analytics and XVAs, trade/market-data interfaces, example use cases, and comprehensive test suites; it is explicitly based on QuantLib and positioned as a benchmark/validation/training reference and as a foundation for tailored risk solutions. The Tagliani–Milev paper’s focus on positivity preservation and discrete maximum principles addresses implementation-level correctness properties that matter when embedding PDE solvers into such end-to-end analytics stacks. Scale: public sources describe sponsorship and open-source availability; deployment scale is not quantified in the repository README. ([github.com](https://github.com/OpenSourceRisk/Engine))

## 9. Consolidated Reference List

[1] Mohammad Mehdizadeh Khalsaraei, Ali Shokri, Higinio Ramos, Zahra Mohammadnia, Pari Khakzad. “A Positivity-Preserving Improved Nonstandard Finite Difference Method to Solve the Black-Scholes Equation.” *Mathematics*, 2022. DOI: 10.3390/math10111846. URL: `https://doi.org/10.3390/math10111846`. ([mdpi.com](https://www.mdpi.com/2227-7390/10/11/1846))

[2] Mohammad Mehdizadeh Khalsaraei, Ali Shokri, Yuanheng Wang, Sohrab Bazm, Giti Navidifar, Pari Khakzad. “Qualitatively Stable Schemes for the Black–Scholes Equation.” *Fractal and Fractional*, 2023. DOI: 10.3390/fractalfract7020154. URL: `https://doi.org/10.3390/fractalfract7020154`. ([doaj.org](https://doaj.org/article/d900141343ca4097a79d4e4505bdf2ff))

[3] Agus Sugandha, Endang Rusyaman, Sukono, Ema Carnia. “Using a Mix of Finite Difference Methods and Fractional Differential Transformations to Solve Modified Black–Scholes Fractional Equations.” *Mathematics*, 2024. DOI: 10.3390/math12071077. URL: `https://doi.org/10.3390/math12071077`. ([mdpi.com](https://www.mdpi.com/2227-7390/12/7/1077))

[4] Y. Ma, C.Z. Shi, Y.C. Hon. “Generalized Finite Integration Method with Laplace transform for European option pricing under Black–Scholes and Heston models.” *Engineering Analysis with Boundary Elements*, 2024. DOI: 10.1016/j.enganabound.2024.105751. URL: `https://doi.org/10.1016/j.enganabound.2024.105751`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0955799724002169))

[5] Jingtang Ma, Zhiqiang Zhou, Zhenyu Cui. “Hybrid Laplace transform and finite difference methods for pricing American options under complex models.” *Computers & Mathematics with Applications*, 2017. DOI: 10.1016/j.camwa.2017.04.018. URL: `https://doi.org/10.1016/j.camwa.2017.04.018`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0898122117302493))

[6] Rami Ahmad El-Nabulsi, Waranont Anukool. “Black–scholes equation in quantitative finance with variable parameters: a path to a generalized schrodinger equation.” *Financial Innovation*, 2026. DOI: 10.1186/s40854-025-00877-7. URL: `https://doi.org/10.1186/s40854-025-00877-7`. ([link.springer.com](https://link.springer.com/article/10.1186/s40854-025-00877-7))

[7] lballabio. “QuantLib: the free/open-source library for quantitative finance.” GitHub repository, ongoing. URL: `https://github.com/lballabio/QuantLib`. ([github.com](https://github.com/lballabio/QuantLib?utm_source=openai))

[8] OpenSourceRisk. “Open Source Risk Engine (ORE).” GitHub repository, ongoing. URL: `https://github.com/OpenSourceRisk/Engine`. ([github.com](https://github.com/OpenSourceRisk/Engine))

---
Learn more:
1. [https://www.sciencedirect.com/science/article/pii/S0096300313007613](https://www.sciencedirect.com/science/article/pii/S0096300313007613)
2. [https://www.sciencedirect.com/science/article/abs/pii/S0955799724002169](https://www.sciencedirect.com/science/article/abs/pii/S0955799724002169)
3. [https://www.sciencedirect.com/science/article/pii/S0898122117302493](https://www.sciencedirect.com/science/article/pii/S0898122117302493)
4. [https://www.mdpi.com/2227-7390/10/11/1846](https://www.mdpi.com/2227-7390/10/11/1846)
5. [https://doaj.org/article/d900141343ca4097a79d4e4505bdf2ff](https://doaj.org/article/d900141343ca4097a79d4e4505bdf2ff)
6. [https://www.mdpi.com/2227-7390/12/7/1077](https://www.mdpi.com/2227-7390/12/7/1077)
7. [https://link.springer.com/article/10.1186/s40854-025-00877-7](https://link.springer.com/article/10.1186/s40854-025-00877-7)
8. [GitHub - lballabio/QuantLib: The QuantLib C++ library](https://github.com/lballabio/QuantLib?utm_source=openai)
9. [https://github.com/OpenSourceRisk/Engine](https://github.com/OpenSourceRisk/Engine)

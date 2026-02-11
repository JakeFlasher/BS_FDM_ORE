## 1. Paper Identity

- **Title:** *Low Volatility Options and Numerical Diffusion of Finite Difference Schemes* ([eudml.org](https://eudml.org/doc/281441))  
- **Authors:** Mariyan Milev; Aldo Tagliani ([eudml.org](https://eudml.org/doc/281441))  
- **Affiliations (as stated in the paper):**
  - Mariyan Milev — Department of Applied Mathematics, University of Venice (Università Ca’ Foscari Venezia), Venice, Italy
  - Aldo Tagliani — Department of Computer and Management Sciences, Trento University, Trento, Italy
- **Venue / bibliographic record:** *Serdica Mathematical Journal*, pp. 223–236, 2010 ([eudml.org](https://eudml.org/doc/281441))  
- **Volume / issue ambiguity (must be flagged):**
  - EuDML indexes the paper as *Serdica Mathematical Journal* **35**(3), 2010, pp. 223–236. ([eudml.org](https://eudml.org/doc/281441))  
  - The University of Trento IRIS record indexes it as **36**(3), 2010, pp. 223–236. ([iris.unitn.it](https://iris.unitn.it/handle/11572/85854?utm_source=openai))  
  - The PDF header provided in the prompt states “Serdica Math. J. 36 (2010), 223–236.” (prompt source; no external citation available).
- **DOI / arXiv ID:** No DOI or arXiv identifier is provided in the EuDML bibliographic record. ([eudml.org](https://eudml.org/doc/281441))  

## 2. Problem Statement & Formulation

The paper studies finite-difference discretizations of the one-factor Black–Scholes PDE in regimes where (i) the terminal/boundary data are **discontinuous** (barrier-type clauses, truncated/digital-style payoffs, discrete monitoring updates) and (ii) volatility is **very low**, typically summarized by the regime $$\sigma^2 \ll r$$ in the numerical examples. The numerical pathology emphasized is **spurious oscillations** produced by standard centered schemes near discontinuities (prices can become negative), and the competing pathology is **artificial numerical diffusion** introduced by non-oscillatory/monotone schemes (smoothing of sharp features). The technical goal is to (a) summarize two nonstandard schemes used to avoid oscillations and preserve positivity and (b) identify the origin and size of the schemes’ artificial diffusion and how grid choices can reduce it (primarily by shrinking $$\Delta S$$, with related restrictions on $$\Delta t$$ in some cases).

### Black–Scholes PDE (paper’s baseline)
Let $$V(S,t)$$ denote the option price as a function of asset price $$S$$ and “time to expiry” $$t$$ with $$0 \le t \le T$$ (paper convention). The PDE is
$$
-\frac{\partial V}{\partial t}
+ rS \frac{\partial V}{\partial S}
+ \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}
- rV = 0.
$$

The paper also rewrites it in a generic convection–diffusion–reaction form
$$
-\frac{\partial V}{\partial t}
+ \mu(S,t)\frac{\partial V}{\partial S}
+ \sigma(S,t)\frac{\partial^2 V}{\partial S^2}
+ b(S,t)V = 0,
$$
with the Black–Scholes identification
$$
\sigma(S,t) = \frac{1}{2}\sigma^2 S^2,\quad
\mu(S,t) = rS,\quad
b(S,t) = -r.
$$
(The notation overload $$\sigma(S,t)$$ vs volatility $$\sigma$$ is inherited from the paper; $$\sigma(S,t)$$ above is the second-derivative coefficient.)

### Discontinuous payoff / monitoring formulations used in experiments

- **Truncated payoff call (Definition 4.1):** payoff at maturity (paper uses $$t=0$$ as the “initial” condition for time-marching)
  $$
  f[S(T)] =
  \begin{cases}
  S(T) - K, & S(T) \in [K,U],\\
  0, & \text{otherwise}.
  \end{cases}
  $$
  The payoff has a discontinuity at $$S=U$$ (jump from $$U-K$$ to $$0$$).

- **Discrete double barrier knock-out call (Definition 4.2 + conditions (11)–(13)):** payoff condition $$\max(S-K,0)$$ but the contract is knocked out (value set to zero) if at discrete monitoring dates the asset lies outside $$[L,U]$$. The PDE is combined with:
  $$
  V(S,0) = (S-K)^+ \, \mathbf{1}_{[L,U]}(S),
  $$
  boundary condition (as stated in the paper)
  $$
  V(S,t) \to 0 \quad \text{as} \quad S \to 0 \quad \text{or} \quad S \to \infty,
  $$
  and monitoring-date “reset” updates
  $$
  V(S,t_i) = V(S,t_i^-)\, \mathbf{1}_{[L,U]}(S),
  \quad 0=t_0 < t_1 < \cdots < t_F = T.
  $$
  Here $$\mathbf{1}_{[L,U]}(S)$$ is the indicator of the corridor.  

**Ambiguity that must be flagged:** the boundary condition $$V(S,t)\to 0$$ as $$S\to\infty$$ is not the standard far-field condition for a vanilla call, but it is stated verbatim for the discretely monitored double-barrier example; the paper does not give further justification beyond presenting it as part of the problem statement.

## 3. Core Methodology

The paper’s “core methods” are two nonstandard finite-difference schemes used to avoid oscillations and preserve positivity in the low-volatility + discontinuity regime:

1. An implicit **exponentially fitted** scheme (Duffy-style fitting factor) for convection-dominated problems.
2. A **Crank–Nicolson variant** tailored to Black–Scholes, modifying only the reaction term discretization via a six-node stencil with tuned weights.

A shared theme is enforcing a matrix structure consistent with positivity (via **M-matrix** arguments), while recognizing that monotonicity/positivity can come with extra smoothing equivalent to adding an artificial diffusion term to the PDE.

### 3.1 Common workflow elements (both schemes)

- Spatial grid in $$S$$ with spacing $$\Delta S$$ (paper also uses $$h$$) and node values $$S_j$$.
- Time grid in $$t$$ with spacing $$\Delta t$$ (paper also uses $$k$$) and times $$t_n$$.
- Backward (in calendar time) valuation implemented as forward marching in “time-to-expiry” $$t$$ from $$t=0$$ (payoff) to $$t=T$$ (present, by paper convention).
- For discrete monitoring: at each monitoring date $$t_i$$, enforce knock-out by multiplying the solution vector by $$\mathbf{1}_{[L,U]}(S_j)$$ (hard reset).

```text
┌──────────────────────────────┐
│ Contract spec + discontinuity │
│  (payoff, barriers, dates)    │
└───────────────┬──────────────┘
                │  defines
                ▼
┌──────────────────────────────┐
│ Black–Scholes PDE + IC/BC     │
│  + monitoring resets (if any) │
└───────────────┬──────────────┘
                │  choose grid (ΔS, Δt)
                ▼
┌──────────────────────────────────────────────────┐
│ Nonstandard FD time-march (pick one):             │
│  (A) Implicit exponentially fitted (Duffy)        │
│  (B) Crank–Nicolson variant (Milev–Tagliani)      │
└───────────────┬──────────────────────────────────┘
                │  at each monitoring date t_i:
                │  apply indicator 1_[L,U](S)
                ▼
┌──────────────────────────────┐
│ Output: price surface U_j^n   │
│  + inspect diffusion/oscill.  │
└──────────────────────────────┘
```

### 3.2 Implicit exponentially fitted finite-difference scheme (Section 2)

#### Discretized operator and fitting factor
Define the continuous operator
$$
\mathcal{L}V \equiv
-\frac{\partial V}{\partial t}
+ \mu(S,t)\frac{\partial V}{\partial S}
+ \sigma(S,t)\frac{\partial^2 V}{\partial S^2}
+ b(S,t)V.
$$

Define the fitted discrete operator (paper notation, time-implicit at level $$n+1$$)
$$
\mathcal{L}^{h}_{k} U^{n}_{j}
\equiv
-\frac{U^{n+1}_{j}-U^{n}_{j}}{k}
+ \mu^{n+1}_{j}\frac{U^{n+1}_{j+1}-U^{n+1}_{j-1}}{2h}
+ \rho^{n+1}_{j}\frac{\delta_x^2 U^{n+1}_{j}}{h^2}
+ b^{n+1}_{j}U^{n+1}_{j},
$$
where $$\delta_x^2 U^{n+1}_{j} = U^{n+1}_{j+1}-2U^{n+1}_{j}+U^{n+1}_{j-1}$$.

The **fitting factor** $$\rho^{n+1}_{j}$$ is chosen as
$$
\rho^{n+1}_{j}
\equiv
\frac{\mu^{n+1}_{j}h}{2}\coth\!\left(\frac{\mu^{n+1}_{j}h}{2\sigma^{n+1}_{j}}\right).
$$
The paper notes that in the centered-difference limit the fitting factor becomes identically equal to $$1$$ (interpreted as recovering the standard centered scheme behavior in the appropriate scaling).

#### Linear system form
At each time step, the scheme yields a tridiagonal linear system
$$
A U^{n+1} = U^{n},
$$
with
$$
A = \operatorname{tridiag}(a_{j,j-1},a_{j,j},a_{j,j+1}),
$$
where (consistent with the paper’s displayed tridiagonal structure)
$$
\begin{aligned}
a_{j,j-1} &= -k\left(\frac{\rho^{n}_{j}}{h^2}-\frac{\mu^{n}_{j}}{2h}\right),\\
a_{j,j}   &= 1 + k\left(\frac{2\rho^{n}_{j}}{h^2}-b^{n}_{j}\right)
          = k\left(\frac{2\rho^{n}_{j}}{h^2}-b^{n}_{j}+\frac{1}{k}\right),\\
a_{j,j+1} &= -k\left(\frac{\rho^{n}_{j}}{h^2}+\frac{\mu^{n}_{j}}{2h}\right).
\end{aligned}
$$

#### Algorithm (implicit exponentially fitted scheme)
1. **Inputs:** payoff/initial vector $$U^0_j = V(S_j,0)$$; coefficients $$\mu(S,t),\sigma(S,t),b(S,t)$$ (Black–Scholes: $$\mu=rS$$, $$\sigma=\frac{1}{2}\sigma^2 S^2$$, $$b=-r$$); grid $$\{S_j\}_{j=0}^{M}$$ with step $$h=\Delta S$$; times $$\{t_n\}_{n=0}^{N}$$ with step $$k=\Delta t$$; boundary values consistent with the contract (paper uses truncations such as $$S_{\max}$$ in experiments).
2. For $$n=0,1,\dots,N-1$$:
   1. For each interior node $$j=1,\dots,M-1$$ compute
      $$
      \rho^{n+1}_j = \frac{\mu^{n+1}_j h}{2}\coth\!\left(\frac{\mu^{n+1}_j h}{2\sigma^{n+1}_j}\right).
      $$
   2. Assemble the tridiagonal matrix $$A$$ with coefficients $$a_{j,j-1},a_{j,j},a_{j,j+1}$$ as above.
   3. Impose boundary conditions into the linear system (paper does not provide a single universal BC recipe; it depends on the option class and truncation).
   4. Solve the tridiagonal system $$A U^{n+1}=U^n$$ (one linear solve per time step).
   5. If $$t_{n+1}$$ is a discrete monitoring date, apply the knock-out update
      $$
      U^{n+1}_j \leftarrow U^{n+1}_j \,\mathbf{1}_{[L,U]}(S_j).
      $$
3. **Outputs:** discrete price surface $$U^n_j \approx V(S_j,t_n)$$ and (optionally) derived quantities (the paper discusses positivity and oscillation/diffusion behavior rather than numerical Greeks for these schemes).

### 3.3 Crank–Nicolson variant scheme (Section 3)

#### Reaction-term modification (six-node stencil)
The scheme is standard Crank–Nicolson for the derivative terms except for the Black–Scholes reaction term $$-rV$$. The mid-time value $$V(t+\Delta t/2)$$ is replaced by a nonlocal weighted average over six adjacent nodes:
$$
\begin{aligned}
V\!\left(t+\frac{\Delta t}{2}\right)
&=
\omega_1\left(U^n_{j-1}+U^n_{j+1}\right)
+\left(\frac{1}{2}-2\omega_1\right)U^n_j \\
&\quad+
\omega_2\left(U^{n+1}_{j-1}+U^{n+1}_{j+1}\right)
+\left(\frac{1}{2}-2\omega_2\right)U^{n+1}_j,
\end{aligned}
$$
with parameters $$\omega_1,\omega_2$$ chosen to enforce positivity/monotonicity properties.

#### Linear system form
The discretization yields
$$
P U^{n+1} = N U^{n},
$$
where $$P$$ and $$N$$ are tridiagonal matrices (paper gives explicit entries). The parameter choice is constrained by:
- $$P$$ is irreducibly diagonally dominant, hence an M-matrix (so $$P^{-1}>0$$),
- $$N$$ has nonnegative entries.

The paper states these are satisfied by
$$
\omega_1 = \omega_2 = -\frac{r}{16\sigma^2},
$$
together with a time-step restriction
$$
\Delta t <
\frac{1}{
r\left(\frac{1}{2}-2\omega_1\right)
+\frac{1}{2}(\sigma M)^2
},
$$
where $$M$$ denotes the number of nodes in the $$S$$ direction.

#### Algorithm (Crank–Nicolson variant)
1. **Inputs:** payoff/initial vector $$U^0_j = V(S_j,0)$$; constants/parameters $$r,\sigma$$ (as used in the paper’s examples); grid $$\{S_j\}_{j=0}^{M}$$ with step $$\Delta S$$; times $$\{t_n\}_{n=0}^{N}$$ with step $$\Delta t$$; monitoring dates $$\{t_i\}$$ (if discretely monitored); boundaries/truncation.
2. Choose weights
   $$
   \omega_1=\omega_2=-\frac{r}{16\sigma^2},
   $$
   and enforce the paper’s step constraint on $$\Delta t$$.
3. For $$n=0,1,\dots,N-1$$:
   1. Assemble tridiagonal matrices $$P$$ and $$N$$ using the paper’s coefficients (functions of $$r,\sigma,S_j,\Delta S,\Delta t,\omega_1,\omega_2$$).
   2. Impose boundary conditions in $$P$$ and $$N$$.
   3. Solve $$P U^{n+1} = N U^n$$ (one tridiagonal solve per time step).
   4. If $$t_{n+1}$$ is a monitoring time, apply the reset
      $$
      U^{n+1}_j \leftarrow U^{n+1}_j \,\mathbf{1}_{[L,U]}(S_j).
      $$
4. **Outputs:** discrete surface $$U^n_j$$; the paper’s emphasis is that solutions are positive and oscillation-free under the stated conditions, with controllable diffusion via grid refinement.

### 3.4 Numerical diffusion “modified PDE” viewpoint (paper’s diagnostic)

The paper attributes artificial diffusion to consistency terms:

- For the Duffy fitted scheme, in the convection-dominated limit $$\sigma(S,t)\to 0$$, the fitting factor tends to
  $$
  \lim_{\sigma\to 0}\rho
  =
  \begin{cases}
  \frac{\mu h}{2}, & \mu>0,\\
  -\frac{\mu h}{2}, & \mu<0,
  \end{cases}
  $$
  yielding an implicit first-order upwind scheme whose truncation error corresponds to solving
  $$
  -V_t + \mu V_S + \frac{1}{2}\mu h\, V_{SS} - bV = 0,
  $$
  i.e., adding numerical diffusion coefficient $$\frac{1}{2}\mu h$$. Under Black–Scholes drift $$\mu=rS$$ this becomes an artificial term of order $$\frac{1}{2}rS\Delta S\,V_{SS}$$.

- For the Crank–Nicolson variant, the dominant diffusion contribution in small $$\sigma$$ is attributed to the modified reaction-term discretization, giving an effective diffusion term of size
  $$
  \frac{1}{8}\left(\frac{r}{\sigma}\Delta S\right)^2 V_{SS},
  $$
  so the scheme is interpreted as solving
  $$
  -V_t + rS V_S + \frac{1}{8}\left(\frac{r}{\sigma}\Delta S\right)^2 V_{SS} - rV = 0.
  $$

## 4. Theoretical Results

No results are formally labeled “Theorem/Lemma/Corollary” in the provided text except the numbered definitions in Section 4; the claims below are the paper’s explicit mathematical guarantees (restated in full with conditions) and their immediate consequences.

### 4.1 Implicit exponentially fitted scheme (Section 2)

#### Proposition 4.1 (M-matrix structure of the fitted implicit system)
**Statement.** For the tridiagonal iteration matrix $$A$$ in
$$
A U^{n+1} = U^{n},
$$
with off-diagonal entries strictly negative $$a_{i,i+1}<0$$ and $$a_{i+1,i}<0$$ and diagonal entries strictly positive $$a_{i,i}>0$$, the paper states that $$A$$ is an irreducible diagonally dominant tridiagonal **M-matrix**. Consequently, $$A^{-1}>0$$ (all entries nonnegative).  

**Proof sketch (paper’s argument line).** The sign pattern (positive diagonal, nonpositive off-diagonals), tridiagonality, irreducibility, and diagonal dominance place $$A$$ in the standard sufficient conditions for a nonsingular M-matrix; the paper cites M-matrix results (e.g., Ortega and Windisch) to conclude $$A^{-1}>0$$.

#### Corollary 4.2 (Positivity preservation of the fitted scheme)
**Statement.** If the discrete initial condition satisfies $$U^0 \ge 0$$ componentwise, then the fitted scheme yields
$$
U^{n} = (A^{-1})^{n} U^{0} > 0
$$
for all subsequent time levels (componentwise positivity).  

**Proof sketch.** From Proposition 4.1, $$A^{-1} \ge 0$$ entrywise. Induction shows $$U^{n}=A^{-1}U^{n-1}$$ remains nonnegative and is strictly positive under irreducibility (the paper’s exposition uses repeated multiplication by $$A^{-1}$$).

#### Proposition 4.3 (Discrete maximum principle / stability in $$\|\cdot\|_\infty$$)
**Statement.** The paper states the bound
$$
\|A^{-1}\|_\infty \le \frac{1}{1-kb} = \frac{1}{1+kr} < 1
$$
(using $$b=-r$$), implying the discrete maximum principle
$$
\|U^{n+1}\|_\infty \le \|U^n\|_\infty.
$$

**Proof sketch.** The diagonal dominance of $$A$$ is used to bound the matrix inverse norm (the paper attributes this bound to Windisch). Submultiplicativity then gives $$\|U^{n+1}\|_\infty = \|A^{-1}U^n\|_\infty \le \|A^{-1}\|_\infty \|U^n\|_\infty$$.

#### Proposition 4.4 (Uniform convergence bound independent of volatility size)
**Statement.** If $$V(S,t)$$ is the analytic solution of the Black–Scholes PDE and $$U^n_j$$ the fitted discrete solution, the paper states
$$
|V(S_j,t_n) - U^n_j| \le c(h+k),
$$
where $$c$$ is independent of $$h$$, $$k$$, and $$\sigma$$ (volatility size).  

**Proof sketch.** The paper treats this as a known result for the exponentially fitted scheme: the fitting construction maintains stability and consistency uniformly in the convection-dominated parameter, giving a global error bound that does not degrade as volatility decreases.

#### Proposition 4.5 (Upwind limit and associated numerical diffusion)
**Statement.** In the limit where the PDE diffusion coefficient tends to zero, $$\sigma(S,t)\to 0$$, the fitted factor satisfies
$$
\lim_{\sigma\to 0}\rho
=
\begin{cases}
\frac{\mu h}{2}, & \mu>0,\\
-\frac{\mu h}{2}, & \mu<0,
\end{cases}
$$
so the fitted scheme becomes the first-order implicit upwind discretization. A standard consistency analysis yields that the upwind scheme solves a parabolic “modified equation”
$$
-V_t + \mu V_S + \frac{1}{2}\mu h\,V_{SS} - bV = 0,
$$
which contains an artificial diffusion term $$\frac{1}{2}\mu h\,V_{SS}$$.

**Proof sketch.** The limit follows from the asymptotic behavior of $$\coth(x)$$ as its argument grows. The modified-equation diffusion term is the standard leading truncation term of first-order upwinding for convection.

### 4.2 Crank–Nicolson variant scheme (Section 3)

#### Proposition 4.6 (Parameter choice ensuring positivity)
**Statement.** Let the Crank–Nicolson variant be written as
$$
P U^{n+1} = N U^{n}.
$$
If parameters are chosen as
$$
\omega_1 = \omega_2 = -\frac{r}{16\sigma^2}
$$
and
$$
\Delta t <
\frac{1}{
r\left(\frac{1}{2}-2\omega_1\right)
+\frac{1}{2}(\sigma M)^2
},
$$
then $$P$$ is an irreducibly diagonally dominant M-matrix (hence $$P^{-1}>0$$) and $$N$$ has nonnegative entries; therefore the update
$$
U^{n+1} = P^{-1} N U^n
$$
preserves positivity whenever $$U^0 \ge 0$$.

**Proof sketch.** The conditions are designed so that $$P$$ has the M-matrix sign/diagonal dominance structure and $$N \ge 0$$ entrywise. Positivity follows by repeated application of the nonnegative matrix $$P^{-1}N$$.

#### Proposition 4.7 (Discrete maximum principle / contractivity bound)
**Statement.** Under the same condition, the paper states
$$
\|N\|_\infty = \frac{1}{\Delta t} - \frac{r}{2},
\quad
\|P^{-1}\|_\infty \le \left(\frac{1}{\Delta t} + \frac{r}{2}\right)^{-1},
$$
and therefore
$$
\|U^{n+1}\|_\infty
\le
\frac{\frac{1}{\Delta t}-\frac{r}{2}}{\frac{1}{\Delta t}+\frac{r}{2}}
\|U^n\|_\infty
\le
\|U^n\|_\infty.
$$

**Proof sketch.** The paper again relies on matrix norm bounds for tridiagonal M-matrices (citing Windisch) and multiplies the bounds to control $$\|P^{-1}N\|_\infty$$.

#### Proposition 4.8 (Accuracy order and artificial diffusion for small volatility)
**Statement.** The paper states the scheme has discretization error $$O(\Delta S^2,\Delta t^2)$$, but when $$\sigma$$ is small the dominant artificial diffusion comes from the reaction-term discretization and has magnitude
$$
\frac{1}{8}\left(\frac{r}{\sigma}\Delta S\right)^2 V_{SS},
$$
so the scheme is interpreted as solving
$$
-V_t + rS V_S + \frac{1}{8}\left(\frac{r}{\sigma}\Delta S\right)^2 V_{SS} - rV = 0.
$$

**Proof sketch.** Second-order accuracy follows from the Crank–Nicolson backbone plus a second-order-consistent mid-time reaction approximation. The artificial diffusion term is obtained by a consistency (modified-equation) expansion, isolating the leading second-derivative contribution induced by the nonlocal reaction discretization.

## 5. Experimental Evaluation

The “experimental” section consists of controlled numerical tests (not data-driven learning experiments) on barrier-style contracts with discontinuous payoffs and very low volatility, comparing four finite-difference schemes:
- standard Crank–Nicolson,
- standard fully implicit (backward Euler in time, centered difference for $$V_S$$ as stated in Example 4.2),
- implicit exponentially fitted scheme (Duffy),
- Crank–Nicolson variant (Milev–Tagliani).

### 5.1 Setups (contracts, schemes, metrics, hyperparameters)

| Test case (paper examples) | Discontinuity source | PDE model & IC/BC | Grid / domain hyperparameters | Baselines / methods compared | Metrics / observations reported |
|---|---|---|---|---|---|
| Truncated call (Definition 4.1; Example 4.1) | Payoff discontinuity at $$S=U$$ | Black–Scholes PDE with payoff $$f[S(T)]$$ truncated to $$[K,U]$$ | Example 4.1: $$r=0.05$$, $$\sigma=0.001$$, $$T=5/12$$, $$U=70$$, $$K=50$$, $$S_{\max}=140$$, $$\Delta S=0.05$$, $$\Delta t=0.01$$ | CN vs fully implicit vs Duffy fitted vs CN-variant | Presence/absence of spurious oscillations near $$U$$; positivity (nonnegative prices); qualitative diffusion (solution smearing vs analytic solution). |
| Discrete double barrier knock-out call (Definition 4.2; Example 4.2) | Discontinuities renewed at each monitoring date via indicator reset | IC: $$V(S,0)=(S-K)^+\mathbf{1}_{[L,U]}(S)$$; updates $$V(S,t_i)=V(S,t_i^-)\mathbf{1}_{[L,U]}(S)$$; BC as stated $$V\to 0$$ as $$S\to 0,\infty$$ | Example 4.2 (text): $$K=100$$, $$\sigma=0.001$$, $$T=1$$, $$r=0.05$$, $$L=95$$, $$U=110$$; Fig. 2 caption uses $$L=90$$ and $$\Delta S=0.025$$, $$\Delta t=0.001$$ | Fully implicit (centered $$V_S$$) highlighted as failing when $$\sigma^2<r$$ | Spurious oscillations and negative values can occur when positivity not guaranteed; oscillations localized near barriers. |
| Duffy fitted scheme sensitivity to $$r$$ (Example 4.3; Fig. 3) | Numerical diffusion from upwind-like behavior in low diffusion | Same truncated call as above | $$\sigma=0.001$$, $$T=5/12$$, $$U=70$$, $$K=50$$, $$S_{\max}=140$$, $$\Delta S=0.05$$, $$\Delta t=0.01$$; compare $$r=0.01$$ vs $$r=0.5$$ | Duffy fitted | Diffusion increases with $$r$$ via term $$\frac{1}{2}rS\Delta S\,V_{SS}$$; qualitative deterioration of sharp features in plots. |
| CN-variant sensitivity to $$r$$ (Example 4.4; Fig. 4) | Numerical diffusion from reaction-term discretization | Same truncated call as above | Same grid, $$r=0.01$$ vs $$r=0.5$$, $$\sigma=0.001$$ | CN-variant | Diffusion increases strongly with larger $$r$$ and small $$\sigma$$ via term $$\frac{1}{8}\left(\frac{r}{\sigma}\Delta S\right)^2 V_{SS}$$; solution remains positive and oscillation-free in both cases. |
| Diffusion reduction by grid refinement (Example 4.5; Fig. 5) | Grid-controlled artificial diffusion | Same truncated call, high $$r$$ case | $$r=0.5$$, $$\sigma=0.001$$, $$T=5/12$$, $$U=70$$, $$K=50$$, $$S_{\max}=140$$; use $$\Delta S=0.025$$ (vs $$0.05$$), $$\Delta t=0.01$$; also notes that $$\Delta S=0.01$$, $$\Delta t=0.001$$ makes nonstandard schemes nearly indistinguishable | Duffy fitted vs CN-variant | Smaller $$\Delta S$$ diminishes numerical diffusion; accurate results require severe grid restrictions, reducing practical advantage of “uniform convergence” claims. |

**Ambiguities that must be flagged (internal consistency of Section 4):**
- Example 4.2 states lower barrier $$L=95$$, while Fig. 2 caption uses $$L=90$$; the paper excerpt does not reconcile this discrepancy.
- Fig. 2 caption text “just before first monitoring date $$t_1=T$$” suggests a degenerate monitoring schedule (first monitoring at maturity) that is unusual for “discrete monitoring”; the excerpt does not explain whether $$F=1$$ or whether this is a caption error.

### 5.2 Key reported outcomes (faithful to the paper’s narrative)

- Standard Crank–Nicolson is reported to generate spurious oscillations near payoff discontinuities for the low-volatility truncated call example, and these oscillations persist for any choice of $$\Delta S$$ and $$\Delta t$$ in the paper’s test narrative (Fig. 1).  
- The standard fully implicit scheme (with centered discretization of $$\partial V/\partial S$$) is also reported to produce oscillations/negative values when $$\sigma^2<r$$ conditions arise (Fig. 2 discussion).  
- Both nonstandard schemes (Duffy fitted and the CN-variant) are described as oscillation-free and positivity-preserving under their respective conditions, but they introduce artificial diffusion that becomes practically significant when $$r$$ is large and/or $$\sigma$$ is very small.  
- Grid refinement (smaller $$\Delta S$$, and in some remarks smaller $$\Delta t$$) is presented as the primary knob to reduce diffusion; the paper emphasizes the resulting computational cost and notes that “an optimal scheme does not exist” for extreme cases.

## 6. ASCII Architecture / Workflow Diagram(s)

```text
┌──────────────────────────────────────────────────────────────────────────────┐
│ Inputs                                                                        │
│  - Option contract: payoff f(S), barriers [L,U], monitoring dates {t_i}       │
│  - Market params: r, σ (possibly r(t,S), σ(t,S) in general form)              │
│  - Grid: S_j = j·ΔS (j=0..M),  t_n = n·Δt (n=0..N), S_max truncation          │
└───────────────────────────────┬──────────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ Initialize                                                                    │
│  U_j^0 = payoff at maturity (paper time-to-expiry convention: t=0)            │
└───────────────────────────────┬──────────────────────────────────────────────┘
                                │  for n = 0..N-1
                                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ Time step update (choose one nonstandard scheme)                              │
│                                                                              │
│  (A) Duffy implicit exponentially fitted:                                     │
│      compute ρ_j^{n+1} = (μ_j^{n+1}ΔS/2)·coth( (μ_j^{n+1}ΔS)/(2σ_j^{n+1}) )    │
│      assemble tridiagonal A(ρ,μ,b); solve  A U^{n+1} = U^n                    │
│                                                                              │
│  (B) Milev–Tagliani CN-variant:                                               │
│      choose ω1=ω2= -r/(16σ^2) + enforce Δt constraint                         │
│      assemble tridiagonal P,N; solve  P U^{n+1} = N U^n                       │
└───────────────────────────────┬──────────────────────────────────────────────┘
                                │  if t_{n+1} ∈ monitoring dates
                                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ Discrete monitoring reset (if applicable)                                     │
│  U_j^{n+1} ← U_j^{n+1} · 1_[L,U](S_j)                                         │
└───────────────────────────────┬──────────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ Outputs                                                                       │
│  - Price surface U_j^n                                                        │
│  - Diagnostics: oscillations vs diffusion                                     │
│    (modified-equation diffusion terms used for interpretation)                │
└──────────────────────────────────────────────────────────────────────────────┘
```

```text
┌─────────────────────────────────────┐
│ Artificial diffusion terms (paper)  │
└───────────────┬─────────────────────┘
                │
                ├─→ Duffy fitted (σ→0 ⇒ upwind):
                │     adds  (1/2)·μ·ΔS · V_SS
                │     Black–Scholes μ=rS ⇒ (1/2)·rSΔS·V_SS
                │
                └─→ CN-variant (small σ):
                      adds  (1/8)·( (r/σ)ΔS )^2 · V_SS
                      requires ΔS small (and Δt constrained) to be negligible
```

## 7. Follow-Up Works & Extensions

### 7.1 Positivity-preserving / oscillation-damping finite-difference refinements

**Milev & Tagliani (2013) propose a modified implicit time-marching procedure that targets the same failure mode emphasized in the 2010 paper—spurious oscillations under discontinuous payoff and low volatility—by altering an implicit scheme to ensure a smooth, oscillation-free, positivity-preserving solution, and using it only for a few initial steps within a Rannacher-style strategy.** [Milev & Tagliani, *J. Comput. Appl. Math.* 2013] The abstract explicitly frames the method as a remedy when $$\sigma^2 \ll r$$ and oscillations arise, then advocates switching to higher-order methods after smoothing the initial data. Relevance note: direct extension of the “positivity + no oscillations” design principle, now applied through an initial damping/smoothing phase rather than accepting diffusion indefinitely. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042712004128))  

**Cen & Le (2011) develop a finite-difference method for a generalized Black–Scholes equation on a piecewise uniform mesh, with stability for arbitrary volatility and interest rate and second-order spatial convergence, explicitly addressing singularities of nonsmooth payoffs (a root cause of the oscillations discussed in the source paper).** [Cen & Le, *J. Comput. Appl. Math.* 2011] The method’s stated properties (stability for arbitrary volatility/rate; handling payoff singularities) align with the source paper’s focus on convection-dominated/degenerate regimes and nonsmooth terminal data. Relevance note: alternative route to qualitative stability (mesh design + analysis) rather than exponential fitting or reaction-term nonlocality. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S037704271100029X))  

**Wang (2004) introduces a fitted finite-volume discretization for the (degenerate) Black–Scholes operator, proves stability and an error bound for the spatial discretization, and shows the resulting system matrix is an M-matrix so the discrete maximum principle holds.** [Wang, *IMA J. Numer. Anal.* 2004] This work is methodologically adjacent to the M-matrix/maximum-principle framing used in the source paper to justify positivity and non-oscillatory behavior. Relevance note: fitted discretizations for Black–Scholes that encode monotonicity via M-matrix structure, predating and thematically supporting later “qualitatively stable” schemes. ([academic.oup.com](https://academic.oup.com/imajna/article/24/4/699/687386))  

**Valkov (2012/2013) proposes a fitted finite-volume method for a generalized Black–Scholes equation transformed to a finite interval, proving unique solvability and positivity preservation for a $$\theta$$-weighted full discretization, with numerical experiments.** [Valkov, *Numerical Algorithms* 2013 / arXiv 2012] The paper explicitly highlights positivity preservation and degeneration at boundaries (a related numerical difficulty to convection-dominance/low-volatility issues). Relevance note: a rigorously analyzed positivity-preserving fitted method that generalizes the “fitting + monotonicity” approach beyond the particular fitted factor used by Duffy-style schemes. ([arxiv.org](https://arxiv.org/abs/1211.1903))  

### 7.2 Hybrid transform + PDE approaches for discontinuities and low volatility

**Tagliani & Milev (2013) combine the Laplace transform (with Post–Widder inversion) and finite differences to treat discretely monitored barrier options, prove equivalence between the mixed method and a fully implicit finite-difference scheme, and emphasize positivity preservation, discrete maximum principle, and avoidance of spurious oscillations under low volatility.** [Tagliani & Milev, *Appl. Math. Comput.* 2013] The paper’s introduction snippet explicitly discusses discontinuities renewed at monitoring dates, convection-dominated behavior for low volatility, and the need to avoid numerical diffusion/oscillations—directly overlapping the source paper’s problem statement. Relevance note: transform-based time handling + qualitative constraints provide an alternative to either (i) exponential fitting or (ii) CN reaction-term modification. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0096300313007613))  

**Gzyl, Milev & Tagliani (2017) use the Mellin transform for Black–Scholes with time-dependent parameters and discontinuous payoffs, then invert the transform using a Maximum Entropy method constrained by fractional moments, arguing that typical finite-difference drawbacks for discontinuous payoffs can be bypassed.** [Gzyl et al., *Finance Research Letters* 2017] The abstract and venue metadata explicitly position the method for discontinuous payoff pricing under Black–Scholes and include a barrier-option example with discrete monitoring. Relevance note: replaces grid-based diffusion/oscillation trade-offs with transform inversion + MaxEnt reconstruction. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S1544612316302562?utm_source=openai))  

**Khalsaraei et al. (2023) present a “qualitatively stable” approach combining Laplace transforms with a nonstandard finite-difference strategy, and explicitly cite the 2010 Milev–Tagliani numerical-diffusion paper as part of the low-volatility/positivity literature.** [Khalsaraei et al., *Fractal Fract.* 2023] The article synopsis states positivity, stability, and consistency under low volatility, and its reference list includes the 2010 source paper. Relevance note: later work that treats the same qualitative desiderata (positivity + low-volatility robustness) using NSFD design rather than fitted factors or CN variants. ([mdpi.com](https://www.mdpi.com/2504-3110/7/2/154?utm_source=openai))  

### 7.3 Discrete barrier options: alternative numerical frameworks (quadrature / recursion / integral equations)

**Hong, Lee & Li (2015) develop recombining quadrature (multinomial-tree-like) methods for discretely monitored barrier options and benchmark them against trapezoidal/Simpson integration and a “Milev–Tagliani (2010)” method, reporting faster convergence for several quadrature choices.** [Hong et al., *J. Comput. Appl. Math.* 2015] The paper is directly relevant because it targets the same discretely monitored barrier setting where discontinuities are renewed at monitoring dates, but it avoids PDE time-marching diffusion/oscillation issues by adopting quadrature-based recursion. Relevance note: an independent numerical route for discrete monitoring that explicitly compares against Milev–Tagliani methods. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042714003793?utm_source=openai))  

**Farnoosh, Fezazadeh & Sobhani (2015) price discrete double barrier options with time-dependent parameters by transforming each monitoring interval to a convenient constant-coefficient Black–Scholes/heat-equation form and computing a recursive representation efficiently, also providing Greeks.** [Farnoosh et al., *Comput. Math. Appl.* 2015] The method addresses the same discrete-monitoring structure but emphasizes analytic transforms and recursion rather than monotone FD schemes; it is compatible with repeated discontinuity resets at monitoring dates. Relevance note: transform + recursion can reduce reliance on diffusion-inducing monotone schemes when barriers are discretely monitored. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0898122115003855?utm_source=openai))  

**Ballestra & Pacelli (2011) propose a boundary element method for time-dependent (moving) double-barrier options under Black–Scholes and CEV, deriving an integral representation with unknown boundary terms solved via Volterra integral equations, and report fast/accurate performance even for non-differentiable barrier functions.** [Ballestra & Pacelli, *Appl. Math. Comput.* 2011] Although focused on moving barriers rather than low volatility per se, it directly addresses barrier-induced irregularities and provides an alternative to FD schemes that may suffer diffusion/oscillation artifacts near barriers. Relevance note: integral-equation reformulation shifts the numerical difficulty away from stabilizing convection-dominated PDE grids. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0096300311012227?utm_source=openai))  

## 8. Industrial & Real-World Applications

### QuantLib finite-difference pricing engines (barrier and vanilla options)
QuantLib provides production-oriented implementations of finite-difference engines for equity/FX-style Black–Scholes problems, including a finite-difference barrier engine that constructs an $$S$$-mesher, applies boundary/step conditions, and calls an FDM Black–Scholes solver using a solver descriptor containing a configurable number of **damping steps** (an industry-standard mechanism to reduce oscillations from nonsmooth payoffs, closely related in spirit to the source paper’s emphasis on damping oscillations near discontinuities). [GitHub: lballabio/QuantLib] The barrier engine implementation explicitly stores a damping-steps parameter and passes it into the FDM solver description used for backward time-marching. Relevance note: the paper’s oscillation-vs-diffusion trade-off is directly actionable in such engines via grid refinement and damping-step configuration for discontinuous payoffs/barrier features. ([github.com](https://github.com/lballabio/QuantLib?utm_source=openai))  

### Open Source Risk Engine (ORE): risk analytics platform built on QuantLib
The Open Source Risk Engine is an open-source platform for risk analytics and XVA that is explicitly based on QuantLib and extends it with additional models, instruments, and pricing engines. [GitHub: OpenSourceRisk/Engine] Relevance note: while ORE is not a “scheme paper” itself, it represents a real-world codebase where finite-difference engines and qualitative-stability choices (positivity, damping, grid resolution) materially affect barrier/exotic option valuation and downstream risk measures when payoffs are discontinuous or monitoring conditions introduce resets. ([github.com](https://github.com/OpenSourceRisk/Engine?utm_source=openai))  

### finmath-lib: open-source library including finite-difference methods for Black–Scholes
finmath-lib is a widely used open-source mathematical finance library that includes finite-difference methods (theta-scheme) for models including Black–Scholes and products including European options. [GitHub: finmath/finmath-lib] Relevance note: this is an accessible industrial-grade reference implementation where monotone/implicit theta schemes are used in practice, illustrating the same general stabilization principle discussed in the source paper (avoiding oscillations in difficult regimes, while managing diffusion via discretization choices). ([github.com](https://github.com/finmath/finmath-lib?utm_source=openai))  

## 9. Consolidated Reference List

[1] Mariyan Milev, Aldo Tagliani. “Efficient implicit scheme with positivity preserving and smoothing properties.” *Journal of Computational and Applied Mathematics*, 2013. DOI: `10.1016/j.cam.2012.09.039`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042712004128))  

[2] Aldo Tagliani, Mariyan Milev. “Laplace Transform and finite difference methods for the Black–Scholes equation.” *Applied Mathematics and Computation*, 2013. DOI: `10.1016/j.amc.2013.07.011`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0096300313007613))  

[3] H. Gzyl, M. Milev, A. Tagliani. “Discontinuous payoff option pricing by Mellin transform: A probabilistic approach.” *Finance Research Letters*, 2017. DOI: `10.1016/j.frl.2016.10.011`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S1544612316302562?utm_source=openai))  

[4] Yicheng Hong, Sungchul Lee, Tianguo Li. “Numerical method of pricing discretely monitored Barrier option.” *Journal of Computational and Applied Mathematics*, 2015. DOI: `10.1016/j.cam.2014.08.022`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042714003793?utm_source=openai))  

[5] R. Farnoosh, H. Fezazadeh, A. M. Sobhani. “Numerical method for discrete double barrier option pricing with time-dependent parameters.” *Computers & Mathematics with Applications*, 2015. DOI: `10.1016/j.camwa.2015.08.016`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0898122115003855?utm_source=openai))  

[6] L. Ballestra, G. Pacelli. “A boundary element method to price time-dependent double barrier options.” *Applied Mathematics and Computation*, 2011. DOI: `10.1016/j.amc.2011.09.050`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0096300311012227?utm_source=openai))  

[7] Zhongdi Cen, Anbo Le. “A robust and accurate finite difference method for a generalized Black–Scholes equation.” *Journal of Computational and Applied Mathematics*, 2011. DOI: `10.1016/j.cam.2011.01.018`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S037704271100029X))  

[8] Radoslav Valkov. “Fitted Finite Volume Method for a Generalized Black-Scholes Equation Transformed on Finite Interval.” *Numerical Algorithms*, 2013 (also arXiv:1211.1903). DOI: `10.1007/s11075-013-9701-3`; arXiv: `1211.1903`. ([arxiv.org](https://arxiv.org/abs/1211.1903))  

[9] Song Wang. “A novel fitted finite volume method for the Black–Scholes equation governing option pricing.” *IMA Journal of Numerical Analysis*, 2004. DOI: `10.1093/imanum/24.4.699`. ([academic.oup.com](https://academic.oup.com/imajna/article/24/4/699/687386))  

[10] QuantLib contributors. “QuantLib: The QuantLib C++ library.” GitHub repository. URL: `https://github.com/lballabio/QuantLib`. ([github.com](https://github.com/lballabio/QuantLib?utm_source=openai))  

[11] Open Source Risk Engine contributors. “Open Source Risk Engine (ORE).” GitHub repository. URL: `https://github.com/OpenSourceRisk/Engine`. ([github.com](https://github.com/OpenSourceRisk/Engine?utm_source=openai))  

[12] finmath contributors. “finmath-lib: Mathematical Finance Library.” GitHub repository. URL: `https://github.com/finmath/finmath-lib`. ([github.com](https://github.com/finmath/finmath-lib?utm_source=openai))

---
Learn more:
1. [EUDML  |  Low Volatility Options and Numerical Diffusion of Finite Difference Schemes](https://eudml.org/doc/281441)
2. [Low volatility options and numerical diffusion of finite difference schemes](https://iris.unitn.it/handle/11572/85854?utm_source=openai)
3. [Efficient implicit scheme with positivity preserving and smoothing properties - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0377042712004128)
4. [A robust and accurate finite difference method for a generalized Black–Scholes equation - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S037704271100029X)
5. [novel fitted finite volume method for the Black–Scholes equation governing option pricing | IMA Journal of Numerical Analysis | Oxford Academic](https://academic.oup.com/imajna/article/24/4/699/687386)
6. [\[1211.1903\] Fitted Finite Volume Method for a Generalized Black-Scholes Equation Transformed on Finite Interval](https://arxiv.org/abs/1211.1903)
7. [Laplace Transform and finite difference methods for the Black–Scholes equation - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0096300313007613)
8. [Discontinuous payoff option pricing by Mellin transform: A probabilistic approach - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1544612316302562?utm_source=openai)
9. [Qualitatively Stable Schemes for the Black–Scholes Equation | MDPI](https://www.mdpi.com/2504-3110/7/2/154?utm_source=openai)
10. [Numerical method of pricing discretely monitored Barrier option - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0377042714003793?utm_source=openai)
11. [Numerical method for discrete double barrier option pricing with time-dependent parameters - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0898122115003855?utm_source=openai)
12. [A boundary element method to price time-dependent double barrier options - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0096300311012227?utm_source=openai)
13. [GitHub - lballabio/QuantLib: The QuantLib C++ library](https://github.com/lballabio/QuantLib?utm_source=openai)
14. [GitHub - OpenSourceRisk/Engine: Open Source Risk Engine](https://github.com/OpenSourceRisk/Engine?utm_source=openai)
15. [GitHub - finmath/finmath-lib: Mathematical Finance Library: Algorithms and methodologies related to mathematical finance.](https://github.com/finmath/finmath-lib?utm_source=openai)

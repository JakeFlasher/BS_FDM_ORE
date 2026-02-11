## 1. Paper Identity

- **Title:** Efficient implicit scheme with positivity preserving and smoothing properties ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042712004128))  
- **Authors:** Mariyan Milev; Aldo Tagliani ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042712004128))  
- **Affiliations (as listed in the paper / indexed metadata):**
  - Department of Informatics and Statistics, University of Food Technologies, Plovdiv, Bulgaria (Milev) ([api.openalex.org](https://api.openalex.org/works/doi%3A10.1016/J.CAM.2012.09.039))  
  - Department of Computer and Management Sciences, University of Trento, Trento, Italy (Tagliani) ([api.openalex.org](https://api.openalex.org/works/doi%3A10.1016/J.CAM.2012.09.039))  
- **Venue:** Journal of Computational and Applied Mathematics, Volume 243, pages 1–9 (issue date shown as 1 May 2013) ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042712004128))  
- **Year:** 2013 (journal issue date); DOI encodes 2012 and the paper front-matter copyright is 2012 ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042712004128))  
- **DOI:** 10.1016/j.cam.2012.09.039 ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042712004128))  
- **Manuscript history (paper front matter):** received 20 Dec 2010; revised 25 Aug 2011.  
- **MSC / keywords (paper front matter):** 65M06, 65M12; Black–Scholes equation; finite difference schemes; fully implicit scheme; M-matrix; non-smooth initial conditions; positivity-preserving.

## 2. Problem Statement & Formulation

The paper targets numerical solution of the Black–Scholes pricing PDE for contracts whose payoff (and/or intermediate “reset” conditions) are discontinuous, a setting where standard finite-difference discretizations can generate **spurious oscillations** and even negative option values, especially in the low-volatility / high-interest regime $$\sigma^2 \ll r$$. The technical objective is an implicit time-stepping modification that (i) preserves non-negativity of the discrete solution, (ii) satisfies a discrete maximum principle (in the paper’s sense), and (iii) damps high-frequency components sufficiently fast to be usable as a smoothing phase inside **Rannacher time-stepping**.

### Risk-neutral asset dynamics

The underlying is assumed to follow a geometric Brownian motion under the risk-neutral measure:
$$
\frac{dS_t}{S_t} = r\,dt + \sigma\, dW_t.
$$

### Black–Scholes PDE (time-to-expiry variable)

Let $$t \in [0,T]$$ denote time-to-expiry and $$V(S,t)$$ the option value. The PDE is:
$$
-\frac{\partial V}{\partial t}
+ rS \frac{\partial V}{\partial S}
+ \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}
- rV
= 0,
\quad S>0,\; t \in (0,T].
$$

### Discretely monitored double-barrier knock-out call (main motivating contract)

For barriers $$L < U$$ and strike $$K$$, the payoff (at $$t=0$$ in the time-to-expiry convention) is:
$$
V(S,0) = \max(S-K,0)\,\mathbf{1}_{[L,U]}(S).
$$

The asymptotic boundary conditions are:
$$
V(S,t) \to 0 \quad \text{as} \quad S \to 0 \quad \text{or} \quad S \to \infty.
$$

Discrete monitoring dates $$0=t_0 < t_1 < \cdots < t_F = T$$ enforce a knock-out “projection”:
$$
V(S,t_i) = V(S,t_i^-)\,\mathbf{1}_{[L,U]}(S),
\quad i=1,\ldots,F.
$$

The discontinuities introduced by $$\mathbf{1}_{[L,U]}(S)$$ at each $$t_i$$ (and the kink/discontinuity in Greeks at the strike at $$t=0$$) are the numerical source of localized oscillations under centered finite differences unless additional smoothing/monotonicity structure is enforced.

## 3. Core Methodology

### 3.1 Computational grid and notation (paper’s finite-difference setup)

The semi-infinite spatial domain is truncated at $$S_{\max}$$. The computational rectangle is:
$$
(S,t) \in [0,S_{\max}] \times [0,T].
$$

Uniform mesh:
- $$S_j = j\,\Delta S$$ for $$j=0,1,\ldots,M$$, with $$S_{\max} = M\,\Delta S$$.
- $$t_n = n\,\Delta t$$ for $$n=0,1,\ldots,X$$, with $$T = X\,\Delta t$$.

Let $$V_j^n \approx V(S_j,t_n)$$.

Boundary-condition enforcement at $$S=0$$ and $$S=S_{\max}$$ is contract-dependent in the paper (e.g., a digital put uses $$V(t,0)=A e^{-rt}$$ and $$V(t,S_{\max})=0$$). The interior update is always expressed as a tridiagonal linear system for each time step.

### 3.2 Stability/quality target used throughout

The paper treats two discrete qualitative requirements as central:

- **Positivity (non-negativity):** $$V_j^n \ge 0$$ should hold for all grid points and time levels, consistent with financial interpretation.
- **Discrete maximum principle (paper’s inequality):**
  $$
  \max_{S \in [0,S_{\max}]} |V(S,t_1)| \ge \max_{S \in [0,S_{\max}]} |V(S,t_2)|
  \quad \text{for} \quad t_1 \le t_2,
  $$
  i.e., the maximum magnitude should not increase with increasing time-to-expiry $$t$$ in regimes where the continuous solution obeys such a principle.

Violation is tied (in the paper’s discussion) to spurious wiggles near sharp gradients and slow damping of discontinuity-induced errors.

### 3.3 Overall workflow (Rannacher-style smoothing with a modified implicit scheme)

The contribution is not a wholesale replacement of higher-order schemes; it is a low-order, monotonicity-oriented implicit step intended for a small number of early steps after a discontinuity is introduced (initial payoff at $$t=0$$; barrier projection at each monitoring date).

ASCII overview of the intended usage:

```text
┌─────────────────────────────────────────────────────────────────────────────┐
│ Discontinuous condition introduced                                           │
│  - at t = 0 (payoff kink / barrier truncation)                               │
│  - at monitoring date t = t_i (projection V ← V · 1_[L,U])                   │
└───────────────┬─────────────────────────────────────────────────────────────┘
                │
                ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ Smoothing phase (few time steps)                                             │
│  Modified fully implicit scheme (this paper)                                 │
│  Goal: positivity + (conditional) max principle + eigenvalues in (0,1)       │
└───────────────┬─────────────────────────────────────────────────────────────┘
                │
                ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ Main time-marching phase                                                     │
│  Switch to higher-order scheme (e.g., Crank–Nicolson family)                 │
│  Start from the smoothed V                                                   │
└─────────────────────────────────────────────────────────────────────────────┘
```

**Ambiguity to flag:** the paper states “few initial time steps” but does not fix the exact number of smoothing steps; this is left to the Rannacher strategy as usually implemented.

### 3.4 Algorithm 1 — Standard fully implicit scheme (centered convection)

This is the baseline scheme analyzed (Section 3.1.1), with the convective term $$\frac{\partial V}{\partial S}$$ discretized by a centered difference.

**Inputs**
1. Parameters $$r,\sigma$$ (constant in the main derivations).
2. Grid $$\Delta S,\Delta t$$ and truncation $$S_{\max}$$.
3. Previous-time solution vector $$V^n = (V_1^n,\ldots,V_{M-1}^n)^\top$$ (interior nodes).
4. Boundary values at $$j=0,M$$ for time levels used (contract-specific).

**Procedure (one time step $$n \to n+1$$)**
1. For each interior node index $$j$$, form the tridiagonal coefficients:
   $$
   A = \operatorname{tridiag}\{ \ell_j,\; d_j,\; u_j \},
   $$
   with
   $$
   \ell_j = -\frac{\Delta t}{2}\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 - r\frac{S_j}{\Delta S}\right],
   $$
   $$
   d_j = 1 + \Delta t\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 + r\right],
   $$
   $$
   u_j = -\frac{\Delta t}{2}\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 + r\frac{S_j}{\Delta S}\right].
   $$
2. Assemble the linear system (after incorporating boundary terms into the right-hand side as appropriate):
   $$
   A V^{n+1} = V^n.
   $$
3. Solve the tridiagonal system for $$V^{n+1}$$.

**Outputs**
- $$V^{n+1}$$ at interior nodes; with boundary values appended to represent $$V(\cdot,t_{n+1})$$ on $$[0,S_{\max}]$$.

**Termination**
- March until $$n=X$$ (i.e., $$t_X = T$$) or until the next monitoring date, where the projection update is applied.

### 3.5 Algorithm 2 — Fully implicit scheme with forward (upwind) convection

The paper notes a “viable approach” for avoiding oscillations when $$\sigma^2 < r$$: discretize $$\frac{\partial V}{\partial S}$$ using a first-order forward difference. The resulting system remains of the form:
$$
A V^{n+1} = V^n,
$$
with the tridiagonal matrix (paper’s stated coefficients):
$$
A = \operatorname{tridiag}\left\{
-\frac{\Delta t}{2}\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2\right],
\;
1+\Delta t\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 + \frac{rS_j}{\Delta S}+r\right],
\;
-\frac{\Delta t}{2}\left[\left(\frac{\sigma S_j}{\Delta S}\right)^2 + \frac{2rS_j}{\Delta S}\right]
\right\}.
$$

**Inputs / outputs / termination** mirror Algorithm 1.

**Key methodological caveat (paper):** this upwind choice is positivity-preserving and oscillation-free but introduces a numerical diffusion term that can “smear” solutions on realistic grids when $$\sigma$$ is small.

### 3.6 Algorithm 3 — Modified implicit scheme (paper’s main proposal)

Section 3.1.2 modifies the standard fully implicit method by changing how the reaction term $$-rV$$ is discretized, using a “bivariate” approximation that couples time levels $$n$$ and $$n+1$$ and neighboring spatial nodes.

#### 3.6.1 Bivariate approximation used for the reaction term

The paper states the following update relation (with local discretization error $$O(\Delta S^2,\Delta t)$$) for a constant $$a$$:
$$
V(S,t+\Delta t) = a\left(V_{j+1}^{n+1} + V_{j-1}^{n+1}\right) + (1-2a)V_j^n.
$$

Centered differences are still used for both $$\frac{\partial V}{\partial S}$$ and $$\frac{\partial^2 V}{\partial S^2}$$.

#### 3.6.2 Linear system form

One time step is written as:
$$
A V^{n+1} = \left[\frac{1}{\Delta t} - r(1-2a)\right] V^n,
$$
with
$$
A = \operatorname{tridiag}\left\{
ra + \frac{r}{2}\frac{S_j}{\Delta S} - \frac{1}{2}\left(\frac{\sigma S_j}{\Delta S}\right)^2,
\;
\frac{1}{\Delta t} + \left(\frac{\sigma S_j}{\Delta S}\right)^2,
\;
ra - \frac{r}{2}\frac{S_j}{\Delta S} - \frac{1}{2}\left(\frac{\sigma S_j}{\Delta S}\right)^2
\right\}.
$$

#### 3.6.3 Parameter choice constraints

The paper selects $$a$$ and $$\Delta t$$ to enforce M-matrix structure and a positive scaling of the right-hand side.

1. M-matrix requirement via nonpositive off-diagonals:
   $$
   ra + \frac{r}{2}\frac{S_j}{\Delta S} - \frac{1}{2}\left(\frac{\sigma S_j}{\Delta S}\right)^2 \le 0,
   $$
   yielding the sufficient constraint:
   $$
   a \le -\frac{r}{8\sigma^2}.
   $$

2. Positive right-hand scaling:
   $$
   \frac{1}{\Delta t} - r(1-2a) > 0
   \quad \Longrightarrow \quad
   \Delta t < \frac{1}{r(1-2a)}.
   $$

The paper then fixes:
$$
a = -\frac{r}{8\sigma^2},
\quad
\Delta t < \frac{1}{r\left(1 + \frac{r}{4\sigma^2}\right)}.
$$

#### 3.6.4 Step-by-step procedure

**Inputs**
1. $$r,\sigma$$; grid $$\Delta S,\Delta t$$; truncation $$S_{\max}$$.
2. Previous-time solution vector $$V^n$$ at interior nodes.
3. Boundary values as needed.
4. Parameter $$a$$ chosen as above (or satisfying the inequalities).

**Procedure (one time step $$n \to n+1$$)**
1. Set $$a = -\frac{r}{8\sigma^2}$$.
2. Check the time-step constraint:
   $$
   \Delta t < \frac{1}{r\left(1+\frac{r}{4\sigma^2}\right)}.
   $$
3. Assemble the tridiagonal matrix $$A$$ with the coefficients in Section 3.6.2.
4. Assemble the right-hand side scaling constant:
   $$
   \beta = \left[\frac{1}{\Delta t} - r(1-2a)\right].
   $$
5. Solve the tridiagonal system:
   $$
   A V^{n+1} = \beta V^n
   $$
   (with boundary contributions handled in the right-hand side).

**Outputs**
- The next-time interior vector $$V^{n+1}$$, with the paper’s proven properties (Section 4 below) under the stated constraints.

### 3.7 Algorithm 4 — Discrete monitoring + Rannacher restart (contract-level loop)

This is the contract-level PDE solve for discretely monitored barriers.

**Inputs**
1. Monitoring dates $$\{t_i\}_{i=0}^F$$ and corridor $$[L,U]$$.
2. Payoff at $$t=0$$ (time-to-expiry convention).
3. Choice of “few” smoothing steps $$m$$ after each discontinuity (not specified by the paper; must be set by the practitioner).
4. Two time-stepping schemes:
   - the modified implicit scheme (Algorithm 3) for smoothing,
   - a chosen higher-order scheme (e.g., Crank–Nicolson) for the main phase.

**Procedure**
1. Initialize $$V(S,0)$$ to the payoff (possibly discontinuous).
2. For each interval $$[t_{i-1},t_i]$$:
   1. Apply Algorithm 3 for $$m$$ time steps starting at $$t=t_{i-1}$$.
   2. Continue to $$t_i^-$$ using the chosen higher-order scheme.
   3. Apply the monitoring update (projection):
      $$
      V(S,t_i) \leftarrow V(S,t_i^-)\,\mathbf{1}_{[L,U]}(S).
      $$
3. End at $$t=T$$.

**Outputs**
- Option value surface (or slice at desired $$t$$) free of the early-time oscillations that would otherwise persist across monitoring resets.

## 4. Theoretical Results

The paper contains no formally labeled theorems/lemmas/propositions/corollaries. Section 3.1.1–3.1.2 provides a sequence of unnumbered spectral/M-matrix/maximum-principle claims; each is restated below with its stated conditions and conclusion, followed by a brief proof strategy sketch matching the paper’s reasoning.

### 4.1 Standard fully implicit (centered) scheme: regime $$\sigma^2 > r$$ (Section 3.1.1)

#### Result 4.1 (Jacobi / tridiagonal structure)
**Claim.** If $$\sigma^2 > r$$, the tridiagonal matrix $$A$$ of the standard fully implicit scheme is a Jacobi matrix (real tridiagonal with standard properties enabling spectral analysis).

**Proof sketch.** The sign pattern under $$\sigma^2 > r$$ makes the off-diagonals nonpositive while diagonal entries are positive; this matches classical Jacobi-matrix conditions invoked by the cited tridiagonal-matrix theory.

#### Result 4.2 (M-matrix, inverse positivity, discrete positivity)
**Claim.** If $$\sigma^2 > r$$, $$A$$ is a row diagonally dominant **M-matrix**, hence $$A^{-1} > 0$$ entrywise, implying $$V^{n+1} = A^{-1}V^n$$ is positive when $$V^n$$ is positive.

**Proof sketch.** M-matrix theory states that a (suitably) diagonally dominant matrix with nonpositive off-diagonals has a nonnegative inverse. The paper invokes this standard implication directly.

#### Result 4.3 ($$\ell_\infty$$ bound and conditional discrete maximum principle)
**Claim.** If $$\sigma^2 > r$$, then
$$
\|A^{-1}\|_\infty \le \frac{1}{1+r\Delta t},
$$
and therefore
$$
\|V^{n+1}\|_\infty \le \frac{1}{1+r\Delta t}\|V^n\|_\infty < \|V^n\|_\infty,
$$
so the maximum principle (in the sense of non-increasing maximum magnitude) is satisfied.

**Proof sketch.** For a row diagonally dominant matrix, standard norm bounds give the stated estimate on the inverse. Multiplying by the update relation yields the maximum-norm contraction.

#### Result 4.4 (real spectrum)
**Claim.** If $$\sigma^2 > r$$, then $$A$$ is similar to a symmetric tridiagonal matrix; hence its eigenvalues $$\lambda_i(A)$$ are real.

**Proof sketch.** A diagonal similarity transform symmetrizes the tridiagonal form when sub-/super-diagonal entries satisfy the required positivity/ratio conditions. Similarity preserves eigenvalues, so reality follows from the symmetric representative.

#### Result 4.5 (Gershgorin-based eigenvalue enclosure for $$A^{-1}$$)
**Claim.** If $$\sigma^2 > r$$, Gershgorin discs imply:
$$
\lambda_i(A^{-1}) \in
\left[
\frac{1}{1+\Delta t\left(r+2(\sigma M)^2\right)},
\;
\frac{1}{1+r\Delta t}
\right]
\subset (0,1),
$$
and the iteration is free of oscillations attributable to negative/complex eigenvalues.

**Proof sketch.** Gershgorin bounds on $$A$$ translate to interval bounds on $$A^{-1}$$ when eigenvalues are real and positive. Keeping eigenvalues in $$ (0,1) $$ avoids sign flips and oscillatory modes in repeated application.

### 4.2 Standard fully implicit (centered) scheme: regime $$\sigma^2 < r$$ (Section 3.1.1)

The paper’s main message is that the “fully implicit is always smoothing” assumption used in Rannacher-type strategies fails when $$\sigma^2 \ll r$$, due to loss of M-matrix structure and (in an extreme subcase) emergence of complex eigenvalues.

#### Result 4.6 (loss of M-matrix and negative entries of $$A^{-1}$$ for $$\frac{r}{M} < \sigma^2 < r$$)
**Claim.** If $$\frac{r}{M} < \sigma^2 < r$$, some subdiagonal entries of $$A$$ become positive, so $$A$$ is not an M-matrix; the inverse $$A^{-1}$$ exists (by the rank-one perturbation analysis referenced), but has negative entries, so $$V^{n+1}=A^{-1}V^n$$ can contain negative components even if $$V^n$$ is positive, creating spurious oscillations.

**Proof sketch.** M-matrix conditions require nonpositive off-diagonals. Once violated, inverse-positivity is not guaranteed; Sherman–Morrison is invoked to argue nonsingularity persists, and known sign-pattern characterizations for such tridiagonals imply negative inverse entries.

#### Result 4.7 (alternating sign pattern of $$A^{-1}$$ for $$\sigma^2 < \frac{r}{M}$$)
**Claim.** If $$\sigma^2 < \frac{r}{M}$$, $$A$$ has positive subdiagonal entries and can be written $$A=\operatorname{tridiag}\{a,c,-b\}$$ with $$a,c,b>0$$. The inverse $$A^{-1}$$ exists (leading principal minors are strictly positive) and exhibits an alternating sign pattern across its entries; consequently, $$V^{n+1} = A^{-1}V^n$$ is intrinsically oscillatory.

**Proof sketch.** The paper uses a recurrence for leading principal minors:
$$
A_k = c_k A_{k-1} + a_{k-1}b_{k-1}A_{k-2},
\quad
A_{-1}=0,\; A_0=1,
$$
to establish nonsingularity. It then quotes a closed-form sign pattern for $$A^{-1}$$ for this sign structure, yielding oscillation in linear combinations.

#### Result 4.8 (complex eigenvalues under $$\sigma^2 < \frac{r}{M}$$)
**Claim.** Under $$\sigma^2 < \frac{r}{M}$$, $$A$$ has $$M$$ distinct complex eigenvalues.

**Proof sketch.** A diagonal similarity transform $$A' = D^{-1}AD$$ with complex diagonal entries converts $$A$$ into a form where $$iA'$$ is symmetric tridiagonal with complex entries. The characteristic polynomial satisfies a three-term recurrence:
$$
p_k(\lambda) = (\lambda - ic_k)p_{k-1}(\lambda) - (-b_{k-1}a_{k-1})p_{k-2}(\lambda),
\quad
p_{-1}=0,\; p_0=1,
$$
with strictly positive recurrence coefficient $$(-b_{k-1}a_{k-1})$$, enabling an “orthogonal polynomial” argument that all zeros are distinct. Distinct zeros correspond to distinct eigenvalues of $$iA'$$ and hence of $$A$$.

#### Result 4.9 (eigen-expansion shows why small $$\Delta t$$ can *increase* oscillation persistence)
**Claim.** Expanding the initial vector in the eigenbasis, the iterates satisfy:
$$
V^{n+1} = (A^{-1})^n V^0 = \sum_{j=1}^M w_j \lambda_j(A^{-1})^n v_j.
$$
When $$\Delta t$$ is “not too small,” eigenvalues satisfy $$|\lambda_j(A^{-1})|<1$$ and are far from $$1$$, so oscillatory components damp quickly; when $$\Delta t$$ is very small, $$|\lambda_j(A^{-1})|\approx 1$$, so damping is slow and oscillations persist, yielding the paradoxical behavior shown in Fig. 1.

**Proof sketch.** The decomposition is standard linear-algebra for diagonalizable matrices with distinct eigenvalues. Gershgorin is used to bound the real parts of eigenvalues of $$A$$ and hence provide a lower bound on $$\Re(\lambda_j(A^{-1}))$$; the qualitative conclusion follows from how powers $$\lambda^n$$ decay depending on $$|\lambda|$$.

### 4.3 Upwind fully implicit scheme: positivity vs numerical diffusion (Section 3.1.1)

#### Result 4.10 (positivity and oscillation-freeness)
**Claim.** With a forward difference for $$\frac{\partial V}{\partial S}$$, the resulting tridiagonal matrix $$A$$ is a row diagonally dominant M-matrix, similar to a symmetric matrix, with real positive eigenvalues; therefore the scheme is positivity-preserving and avoids spurious oscillations.

**Proof sketch.** The forward differencing moves the convection contribution into coefficients that keep off-diagonals nonpositive and maintain diagonal dominance, which restores M-matrix conditions and precludes negative/complex eigenvalues.

#### Result 4.11 (modified-equation diffusion term)
**Claim.** The scheme is consistent with a modified PDE containing an artificial diffusion term:
$$
-\frac{\partial V}{\partial t}
+ rS\frac{\partial V}{\partial S}
+ \left(\frac{rS\Delta S}{2} + \frac{1}{2}\sigma^2 S^2\right)\frac{\partial^2 V}{\partial S^2}
- rV
= 0,
$$
so the additional term $$\frac{rS\Delta S}{2}\frac{\partial^2 V}{\partial S^2}$$ can smear solutions unless $$\Delta S$$ is made very small.

**Proof sketch.** Standard modified-equation analysis of first-order upwinding yields a second-derivative term proportional to the grid spacing and convection strength.

### 4.4 Modified implicit scheme: positivity, stability, spectrum (Section 3.1.2)

All statements below hold under:
$$
a = -\frac{r}{8\sigma^2},
\quad
\Delta t < \frac{1}{r\left(1+\frac{r}{4\sigma^2}\right)}.
$$

#### Result 4.12 (row diagonal dominance margin)
**Claim.** The tridiagonal matrix $$A=[a_{ij}]$$ is row diagonally dominant with margin:
$$
\delta = |a_{ii}| - \sum_{j\ne i}|a_{ij}| = \frac{1}{\Delta t} - \left(\frac{r}{2\sigma}\right)^2 > 0.
$$

**Proof sketch.** Substitute the fixed $$a$$ into the tridiagonal coefficients and compute the diagonal minus the sum of off-diagonal magnitudes; the stated inequality on $$\Delta t$$ makes the margin positive.

#### Result 4.13 (M-matrix, inverse positivity, positivity preservation)
**Claim.** Under the stated conditions, $$A$$ is an M-matrix, hence $$A^{-1}>0$$ and the update
$$
A V^{n+1} = \left[\frac{1}{\Delta t} - r(1-2a)\right]V^n
$$
yields $$V^{n+1}>0$$ whenever $$V^n \ge 0$$.

**Proof sketch.** Off-diagonals are enforced to be nonpositive, diagonal dominance provides nonsingularity, and M-matrix theory implies entrywise nonnegative inverse.

#### Result 4.14 (real spectrum of $$A$$)
**Claim.** $$A$$ is similar to a symmetric tridiagonal matrix, so $$\lambda_i(A)$$ are real.

**Proof sketch.** Same diagonal symmetrization argument as in Section 4.1, applied to this modified coefficient set.

#### Result 4.15 (conditional stability via $$\rho(A_{\mathrm{iter}})<1$$)
Define the iteration matrix:
$$
A_{\mathrm{iter}} =
\left[\frac{1}{\Delta t} - r\left(1+\frac{r}{4\sigma^2}\right)\right]A^{-1}.
$$

**Claim.** The scheme is conditionally stable with:
$$
\rho(A_{\mathrm{iter}}) \le \|A_{\mathrm{iter}}\|_\infty
\le
\frac{\frac{1}{\Delta t}-r-\left(\frac{r}{2\sigma}\right)^2}{\frac{1}{\Delta t}-\left(\frac{r}{2\sigma}\right)^2}
< 1.
$$

**Proof sketch.** Use the row-diagonal-dominance bound for $$\|A^{-1}\|_\infty$$, then bound $$\|A_{\mathrm{iter}}\|_\infty$$ and apply $$\rho(\cdot)\le \|\cdot\|_\infty$$.

#### Result 4.16 (conditional discrete maximum principle)
**Claim.** Under the same conditions,
$$
\|V^{n+1}\|_\infty = \|A_{\mathrm{iter}}V^n\|_\infty \le \|A_{\mathrm{iter}}\|_\infty \|V^n\|_\infty < \|V^n\|_\infty,
$$
so the method satisfies the discrete maximum principle (in the paper’s contraction sense) conditionally.

**Proof sketch.** Combine the stability bound of Result 4.15 with the update relation.

#### Result 4.17 (eigenvalues of $$A_{\mathrm{iter}}$$ lie in $$(0,1)$$, preventing oscillations)
**Claim.** Gershgorin yields an enclosure:
$$
\lambda_i(A_{\mathrm{iter}}) \in
\left[
\frac{\frac{1}{\Delta t}-r-\left(\frac{r}{2\sigma}\right)^2}
{\frac{1}{\Delta t}+\left(\frac{r}{2\sigma}\right)^2 + 2(\sigma M)^2},
\;
\frac{\frac{1}{\Delta t}-r-\left(\frac{r}{2\sigma}\right)^2}
{\frac{1}{\Delta t}-\left(\frac{r}{2\sigma}\right)^2}
\right]
\subset (0,1),
$$
so spurious oscillations associated with negative or complex eigenvalues are avoided.

**Proof sketch.** Gershgorin discs for a real-spectrum matrix constrain eigenvalues to a positive interval when diagonal dominance and sign structure enforce positivity of the relevant bounds.

### 4.5 Consistency, truncation error, and accuracy constraints (Section 3.1.2)

#### Result 4.18 (local truncation error order)
**Claim.** Replacing $$V_j^{n+1}$$ in the standard implicit scheme by the bivariate expression introduces a local truncation error of order:
$$
O(\Delta t,\Delta S^2),
$$
so the modified scheme is first order in time and second order in space.

**Proof sketch.** Taylor-expand about time level $$(n+1)\Delta t$$ and note that centered differences preserve $$O(\Delta S^2)$$ spatial accuracy while the time discretization remains $$O(\Delta t)$$.

#### Result 4.19 (additional low-volatility accuracy conditions)
**Claim.** When $$\sigma^2 \ll r$$, the error contributions proportional to
$$
r\Delta t\left(1+\frac{r}{4\sigma^2}\right)
\quad \text{and} \quad
\frac{1}{8}\left(\frac{r}{\sigma}\Delta S\right)^2
$$
become significant; obtaining an accurate solution requires:
$$
\Delta t\left(r+\frac{r^2}{4\sigma^2}\right)\ll 1
\quad \text{and} \quad
\frac{1}{8}\left(\frac{r}{\sigma}\Delta S\right)^2 \ll 1.
$$

**Proof sketch.** The Taylor expansion shows the effective perturbation terms in the modified-equation sense. In the regime $$\sigma^2 \ll r$$, the coefficients amplify the nominal discretization errors unless $$\Delta t$$ and $$\Delta S$$ satisfy stronger “smallness” conditions than stability alone.

### Computational complexity (inference)
Each time step in Algorithms 1–3 solves a tridiagonal linear system of dimension approximately $$M-1$$ interior nodes. Using the Thomas algorithm gives $$O(M)$$ work and $$O(M)$$ memory per time step; the paper does not state this explicitly.

## 5. Experimental Evaluation

The paper’s “experiments” are PDE test cases with analytic solutions (digital option) or benchmark qualitative behavior (truncated call / barrier truncation) rather than data-driven benchmarks. Evaluation is primarily by (i) visual detection of oscillations/negativity and (ii) comparison to analytic solutions where available.

### 5.1 Experimental setups (parameters, baselines, metrics)

| Case / figure | Contract / PDE setup | Parameters | Grid / horizon | Baselines compared | Metrics / diagnostics |
|---|---|---|---|---|---|
| Fig. 1 | Truncated call option (barrier-truncated call payoff) evolved to just before first monitoring date | $$L=90$$, $$K=100$$, $$U=110$$, $$r=0.05$$, $$\sigma=0.001$$, $$T=1$$, $$S_{\max}=200$$ | $$\Delta S=0.05$$; $$\Delta t \in \{0.01,\,0.0001\}$$ | Standard fully implicit (centered convection) vs analytic reference curve | Presence of spurious oscillations near barrier/strike; negativity |
| Fig. 2 | Same truncated call setup as Fig. 1 | Same as Fig. 1 | $$\Delta S=0.05$$; $$\Delta t \in \{0.01,\,0.0001\}$$ | Exponentially fitted scheme vs analytic | Numerical diffusion / smearing; convergence speed |
| Fig. 4 | Truncated call option, small expiry (Rannacher-style early-time smoothing regime) | $$L=90$$, $$K=100$$, $$U=110$$, $$r=0.05$$, $$\sigma=0.001$$, $$T=0.01$$, $$S_{\max}=120$$ | Exponentially fitted: $$\Delta S=0.01$$, $$\Delta t=10^{-4}$$; Modified implicit: $$\Delta S=0.01$$, $$\Delta t=10^{-6}$$ | Exponentially fitted vs modified implicit vs analytic | Agreement between schemes; oscillation absence; diffusion artifacts |
| Example 4.1 / Fig. 5 | Digital put with discontinuous payoff and known analytic solution | $$A=1$$, $$K=10$$, $$r=0.05$$, $$\sigma=0.01$$, $$T=1$$, $$S_{\max}=20$$; BCs: $$V(t,0)=A e^{-rt}$$, $$V(t,S_{\max})=0$$ | Modified implicit: $$\Delta S=0.05$$, $$\Delta t=0.001$$ | Modified implicit vs analytic; (paper notes standard implicit/CN fail when $$\sigma^2 \ll r$$) | Pointwise visual match; localized error near strike; Delta/Gamma quality; positivity |

### 5.2 Key quantitative/parametric findings reported

#### Fig. 1 (standard fully implicit, low volatility)
- Under $$\sigma=0.001$$ and $$r=0.05$$, the standard fully implicit scheme exhibits oscillations near the truncation region.
- Oscillations are shown for $$\Delta t=0.01$$ and become persistent for the smaller $$\Delta t=0.0001$$, matching the paper’s analysis that very small time steps can slow damping when complex eigenvalues appear.

#### Fig. 2 (exponentially fitted scheme)
- Exponentially fitted results are oscillation-free in the illustrated case, but show visible numerical diffusion (solution smearing) and slow convergence, especially relevant when $$\sigma^2 \ll r$$.

#### Fig. 4 (modified implicit vs exponentially fitted)
- With $$T=0.01$$, $$\sigma=0.001$$, $$r=0.05$$, and $$\Delta S=0.01$$, the modified implicit scheme (using $$\Delta t=10^{-6}$$) produces a curve practically indistinguishable from the exponentially fitted scheme (with $$\Delta t=10^{-4}$$) and the analytic reference.
- The paper interprets this as “equivalence” in practice between the modified implicit and the implicit exponentially fitted scheme for the smoothing purpose.

#### Fig. 5 (digital put; prices and Greeks)
- With $$\Delta S=0.05$$ and $$\Delta t=0.001$$, the modified implicit solution $$V_{\mathrm{imp}}$$ is visually indistinguishable from the analytic $$V_{\mathrm{ex}}$$.
- Errors in $$V$$, Delta, and Gamma are localized near $$S=K$$, consistent with the initial discontinuity.
- The paper emphasizes that standard fully implicit and Crank–Nicolson schemes can fail (negative prices and oscillations) in the low-volatility regime for such discontinuous payoffs, motivating the modified scheme as the initial smoother.

### 5.3 Ablations and sensitivity axes

No explicit ablation table is provided. Sensitivity is demonstrated along:
- **Time-step size** $$\Delta t$$ (Fig. 1: $$0.01$$ vs $$0.0001$$; Fig. 4: $$10^{-4}$$ vs $$10^{-6}$$).
- **Scheme choice** (standard implicit centered; exponentially fitted; modified implicit).
- **Volatility regime** (emphasis on $$\sigma^2 \ll r$$).

## 6. ASCII Architecture / Workflow Diagram(s)

### 6.1 End-to-end workflow for discretely monitored barriers with smoothing restarts

```text
┌──────────────────────────────────────────────────────────────────────────────┐
│ Inputs                                                                       │
│  - PDE params: r, σ; contract: K, L, U; monitoring times {t_i}; S_max         │
│  - Grid: ΔS, Δt; number of smoothing steps m (user-chosen)                   │
└──────────────────────────────────────┬───────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ Initialize at t = 0 (time-to-expiry)                                         │
│  V(S,0) = payoff (may be discontinuous / non-smooth)                         │
└──────────────────────────────────────┬───────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│ For i = 1..F monitoring intervals                                            │
│                                                                              │
│  (A) Smoothing (Rannacher-style)                                              │
│      repeat k = 1..m:                                                        │
│        Solve:  A_mod V^{n+1} = β V^n   (modified implicit, positivity)        │
│                                                                              │
│  (B) Main propagation                                                        │
│      Advance to t_i^- using higher-order scheme (e.g., Crank–Nicolson)        │
│                                                                              │
│  (C) Barrier projection at monitoring date                                   │
│      V(S,t_i) ← V(S,t_i^-) · 1_[L,U](S)                                      │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 6.2 Local stencil intuition for the modified implicit reaction-term coupling

The modified scheme couples $$V_j^n$$ to neighbors at the next time level via:
$$
V(\cdot,t+\Delta t) \approx a(V_{j-1}^{n+1}+V_{j+1}^{n+1}) + (1-2a)V_j^n.
$$

ASCII schematic (time increases downward):

```text
          t = t_{n+1}
        ┌───────────────┐
        │ V_{j-1}^{n+1}  │     │ V_{j+1}^{n+1}  │
        └───────┬───────┘     └───────┬────────┘
                │        a-coupling    │
                └──────────┬───────────┘
                           ▼
                    ┌────────────┐
                    │  V_j^n     │   (weight 1 - 2a)
                    └────────────┘
          t = t_n
```

## 7. Follow-Up Works & Extensions

### 7.1 Qualitative / positivity-oriented schemes for Black–Scholes under low volatility

**Mehdizadeh Khalsaraei, Shokri, Wang, Bazm, Navidifar, Khakzad. “Qualitatively Stable Schemes for the Black–Scholes Equation.” Fractal and Fractional, 2023 (doi:10.3390/fractalfract7020154).** The paper proposes a technique combining Laplace-transform ideas with nonstandard finite-difference strategy and explicitly discusses positivity and stability under low volatility; its reference list includes Milev & Tagliani’s 2013 “Efficient implicit scheme with positivity preserving and smoothing properties,” indicating direct scholarly linkage via citation. [Mehdizadeh Khalsaraei et al., Fractal Fract. 2023] ([mdpi.com](https://www.mdpi.com/2504-3110/7/2/154))  

**Mehdizadeh Khalsaraei, Shokri, Ramos, Mohammadnia, Khakzad. “A Positivity-Preserving Improved Nonstandard Finite Difference Method to Solve the Black-Scholes Equation.” Mathematics, 2022 (doi:10.3390/math10111846).** The work improves a pre-existing low-accuracy scheme by introducing a non-local approximation for the reaction term (their MMADE method) and states sufficient conditions for positivity preservation via a labeled theorem; this is technically close in spirit to Milev–Tagliani’s Section 3.1.2, which also modifies discretization of the reaction term $$-rV$$ to enforce M-matrix structure and positivity (the MDPI page excerpt does not, in the visible snippets, confirm that it cites the 2013 JCAM paper, so it is best categorized as an independent or potentially related development rather than a verified direct follow-up). [Mehdizadeh Khalsaraei et al., Mathematics 2022] ([mdpi.com](https://www.mdpi.com/2227-7390/10/11/1846))  

### 7.2 Verified-direct follow-ups beyond the items above

No additional direct follow-up papers (explicitly verifiable as citing the 2013 JCAM paper) were identified within the available reference set at time of writing; the OpenAlex record indicates nonzero citation activity in multiple years, but enumerating the full set of citing works requires additional verified source extraction beyond what is included here. ([api.openalex.org](https://api.openalex.org/works/doi%3A10.1016/J.CAM.2012.09.039))  

## 8. Industrial & Real-World Applications

**QuantLib (open-source quantitative finance library).** QuantLib is a widely used open-source library positioned for “modeling, trading, and risk management in real-life,” providing a production-grade codebase in which PDE finite-difference option pricing is a standard component in practice; the Milev–Tagliani scheme is directly relevant as an implementable first-step smoother for discontinuous payoffs in finite-difference engines, but no verified evidence was identified here that QuantLib includes this exact modified implicit discretization. [GitHub: lballabio/QuantLib] ([github.com](https://github.com/lballabio/QuantLib))  

**Open Source Risk Engine (ORE).** ORE is presented as a transparent framework for pricing and risk analysis, explicitly based on QuantLib and extending it with additional models, instruments, and pricing engines; this positions ORE as an industrially oriented platform where positivity-preserving and oscillation-damping discretizations (including Rannacher-style smoothing phases) are practically relevant for PDE-based pricing workflows, though no verified statement was located that ORE implements the specific Milev–Tagliani modification. [GitHub: OpenSourceRisk/Engine] ([github.com](https://github.com/OpenSourceRisk/Engine))  

**JQuantLib (Java rewrite inspired by QuantLib).** JQuantLib is a Java library intended for valuation of shares and options, and it is explicitly described as being based on QuantLib; it exemplifies real-world engineering reuse of finite-difference and related numerical machinery in production languages/environments, where monotone/positivity-preserving discretizations are operationally important even if individual research variants are not directly implemented. [GitHub: frgomes/jquantlib] ([github.com](https://github.com/frgomes/jquantlib))  

**ORE-SWIG (language bindings for ORE).** ORE-SWIG provides language bindings that allow ORE and QuantLib-side libraries to be used from environments such as Python; such bindings are a common path for deploying PDE pricers and risk analytics into research notebooks, batch risk jobs, or application stacks, making the “smoothing-only in early steps” design pattern from Milev–Tagliani a plausible integration point (no verified code-level adoption is asserted). [GitHub: OpenSourceRisk/ORE-SWIG] ([github.com](https://github.com/OpenSourceRisk/ORE-SWIG))  

## 9. Consolidated Reference List

[1] M. Mehdizadeh Khalsaraei; A. Shokri; Y. Wang; S. Bazm; G. Navidifar; P. Khakzad. “Qualitatively Stable Schemes for the Black–Scholes Equation.” *Fractal and Fractional*, 2023. doi:10.3390/fractalfract7020154. ([mdpi.com](https://www.mdpi.com/2504-3110/7/2/154))  

[2] M. Mehdizadeh Khalsaraei; A. Shokri; H. Ramos; Z. Mohammadnia; P. Khakzad. “A Positivity-Preserving Improved Nonstandard Finite Difference Method to Solve the Black-Scholes Equation.” *Mathematics*, 2022. doi:10.3390/math10111846. ([mdpi.com](https://www.mdpi.com/2227-7390/10/11/1846))  

[3] QuantLib (project maintainer: L. Ballabio et al.). “QuantLib: the free/open-source library for quantitative finance.” GitHub repository. `https://github.com/lballabio/QuantLib`. ([github.com](https://github.com/lballabio/QuantLib))  

[4] Open Source Risk. “Open Source Risk Engine (ORE).” GitHub repository. `https://github.com/OpenSourceRisk/Engine`. ([github.com](https://github.com/OpenSourceRisk/Engine))  

[5] R. Gomes (and contributors). “JQuantLib: quantitative finance library in 100% Java.” GitHub repository. `https://github.com/frgomes/jquantlib`. ([github.com](https://github.com/frgomes/jquantlib))  

[6] Open Source Risk. “ORE-SWIG: Language bindings for ORE.” GitHub repository. `https://github.com/OpenSourceRisk/ORE-SWIG`. ([github.com](https://github.com/OpenSourceRisk/ORE-SWIG))

---
Learn more:
1. [https://www.sciencedirect.com/science/article/pii/S0377042712004128](https://www.sciencedirect.com/science/article/pii/S0377042712004128)
2. [https://api.openalex.org/works/doi%3A10.1016/J.CAM.2012.09.039](https://api.openalex.org/works/doi%3A10.1016/J.CAM.2012.09.039)
3. [https://www.mdpi.com/2504-3110/7/2/154](https://www.mdpi.com/2504-3110/7/2/154)
4. [https://www.mdpi.com/2227-7390/10/11/1846](https://www.mdpi.com/2227-7390/10/11/1846)
5. [https://github.com/lballabio/QuantLib](https://github.com/lballabio/QuantLib)
6. [https://github.com/OpenSourceRisk/Engine](https://github.com/OpenSourceRisk/Engine)
7. [https://github.com/frgomes/jquantlib](https://github.com/frgomes/jquantlib)
8. [https://github.com/OpenSourceRisk/ORE-SWIG](https://github.com/OpenSourceRisk/ORE-SWIG)

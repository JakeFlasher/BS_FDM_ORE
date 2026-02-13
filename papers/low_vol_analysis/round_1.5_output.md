# Golden Reference Document

**Method:** Positivity-Preserving Finite Difference Schemes for the Black–Scholes PDE with Discontinuous Payoffs and Low Volatility

**Source paper:** M. Milev & A. Tagliani, "Low Volatility Options and Numerical Diffusion of Finite Difference Schemes," *Serdica Math. J.* **36** (2010), 223–236.

**Document status:** Corrected, completed, and self-contained. Supersedes the original paper for implementation purposes.

---

## 1. METHOD OVERVIEW

This document specifies two finite-difference schemes for the one-factor Black–Scholes PDE with constant coefficients, targeting the regime where volatility is small relative to the interest rate (\(\sigma^2 \ll r\)) and the payoff or monitoring conditions introduce discontinuities. Standard centered-difference schemes (including Crank–Nicolson) produce spurious oscillations and negative prices in this regime. **Scheme 1** is a fully implicit backward-Euler scheme in which the second-derivative coefficient is replaced by Duffy's exponentially fitted factor, ensuring unconditional positivity and a discrete maximum principle at the cost of first-order artificial diffusion \(\tfrac{1}{2}rS\,\Delta S\) in the low-volatility limit. **Scheme 2** is a Crank–Nicolson variant in which only the reaction term \(-rV\) is discretized via a six-node stencil with a tuned weight parameter, yielding second-order accuracy in both space and time with artificial diffusion \(\tfrac{1}{8}(r\Delta S/\sigma)^2\), but requiring a time-step constraint for guaranteed positivity. Both schemes produce tridiagonal linear systems solvable by the Thomas algorithm at each time step, and both support discrete barrier monitoring via a post-solve pointwise projection. The governing PDE is the Black–Scholes equation for European-exercise options (calls with truncated or barrier-modified payoffs) on a single underlying asset with constant risk-free rate \(r > 0\) and constant volatility \(\sigma > 0\), discretized on a uniform grid in the original \((S,t)\) coordinates with \(t\) denoting time-to-expiry.

---

## 2. NOTATION TABLE

| Symbol | Type | Units | Definition |
|:---|:---|:---|:---|
| \(S\) | scalar | \$ | Underlying asset price |
| \(t\) | scalar | yr | Time-to-expiry (backward variable): \(t=0\) is maturity, \(t=T\) is the present |
| \(T\) | scalar | yr | Total time to expiry |
| \(V(S,t)\) | function | \$ | Continuous option price |
| \(r\) | scalar | yr\(^{-1}\) | Risk-free interest rate, constant, positive |
| \(\sigma\) | scalar | yr\(^{-1/2}\) | Volatility parameter, constant, positive |
| \(K\) | scalar | \$ | Strike price |
| \(B_\ell\) | scalar | \$ | Lower barrier level |
| \(B_u\) | scalar | \$ | Upper barrier level |
| \(S_{\max}\) | scalar | \$ | Computational domain upper boundary; require \(S_{\max} \geq 2B_u\) |
| \(\Delta S\) | scalar | \$ | Spatial step size (uniform) |
| \(\Delta t\) | scalar | yr | Time step size (uniform) |
| \(M\) | integer | — | Number of spatial intervals: \(M = S_{\max}/\Delta S\) |
| \(N_t\) | integer | — | Number of time steps: \(N_t = T/\Delta t\) |
| \(S_j\) | scalar | \$ | Spatial grid point: \(S_j = j\,\Delta S\), \(j = 0,1,\ldots,M\) |
| \(t_n\) | scalar | yr | Temporal grid point: \(t_n = n\,\Delta t\), \(n = 0,1,\ldots,N_t\) |
| \(U_j^n\) | scalar | \$ | Discrete approximation to \(V(S_j, t_n)\) |
| \(\mathbf{U}^n\) | vector \(\in \mathbb{R}^{M-1}\) | \$ | Interior unknowns: \(\mathbf{U}^n = (U_1^n, U_2^n, \ldots, U_{M-1}^n)^T\) |
| \(\sigma_d(S)\) | function | \$\(^2\)/yr | Diffusion coefficient: \(\sigma_d(S) = \tfrac{1}{2}\sigma^2 S^2\) |
| \(\mu(S)\) | function | \$/yr | Convection coefficient: \(\mu(S) = rS\) |
| \(b\) | scalar | yr\(^{-1}\) | Reaction coefficient: \(b = -r\) |
| \(\mathrm{Pe}_j\) | scalar | — | Péclet number at node \(j\): \(\mathrm{Pe}_j = r\Delta S / (\sigma^2 S_j) = r/(\sigma^2 j)\) |
| \(\rho_j\) | scalar | \$\(^2\)/yr | Fitting factor at node \(j\) (Scheme 1 only) |
| \(\omega\) | scalar | — | Reaction-term weight (Scheme 2 only): \(\omega = -r/(16\sigma^2)\) |
| \(A\) | matrix \(\in \mathbb{R}^{(M-1)\times(M-1)}\) | — | Tridiagonal iteration matrix (Scheme 1) |
| \(P\) | matrix \(\in \mathbb{R}^{(M-1)\times(M-1)}\) | yr\(^{-1}\) | Left-hand tridiagonal matrix (Scheme 2) |
| \(N\) | matrix \(\in \mathbb{R}^{(M-1)\times(M-1)}\) | yr\(^{-1}\) | Right-hand tridiagonal matrix (Scheme 2) |
| \(F\) | integer | — | Number of monitoring intervals (discrete barriers) |
| \(\{t_i^{\mathrm{mon}}\}_{i=1}^{F}\) | set | yr | Monitoring dates: \(0 < t_1^{\mathrm{mon}} < \cdots < t_F^{\mathrm{mon}} \leq T\) |
| \(\mathbf{1}_{[a,b]}(x)\) | function | — | Indicator: 1 if \(x \in [a,b]\), 0 otherwise |

---

## 3. MATHEMATICAL SPECIFICATION

### 3(a). Governing PDE

The Black–Scholes PDE in time-to-expiry form:

$$
-\frac{\partial V}{\partial t} + rS\frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2\frac{\partial^2 V}{\partial S^2} - rV = 0, \qquad S \in (0, S_{\max}),\; t \in (0, T]
$$

Equivalently, in convection–diffusion–reaction form:

$$
-\frac{\partial V}{\partial t} + \mu(S)\frac{\partial V}{\partial S} + \sigma_d(S)\frac{\partial^2 V}{\partial S^2} + b\,V = 0
$$

with \(\mu(S) = rS\), \(\sigma_d(S) = \tfrac{1}{2}\sigma^2 S^2\), and \(b = -r\).

### 3(b). Coordinate System

No coordinate transformation is applied. All discretization is performed in the original \((S,t)\) variables on a truncated domain \([0, S_{\max}] \times [0, T]\). The variable \(t\) is time-to-expiry: the scheme marches forward from \(t = 0\) (maturity, where the payoff is known) to \(t = T\) (the present).

### 3(c). Domain and Boundary Conditions

**Spatial domain:** \([0, S_{\max}]\) with \(S_{\max} = M\,\Delta S\). The parameter \(S_{\max}\) must satisfy \(S_{\max} \geq 2B_u\) so that the far-field boundary condition is accurate. For the truncated call, interpret \(B_u\) as the upper payoff cutoff.

**Boundary conditions (all product types covered by this document):**

$$
U_0^n = 0, \qquad U_M^n = 0, \qquad \text{for all } n = 0, 1, \ldots, N_t
$$

These are homogeneous Dirichlet conditions. The condition at \(S = 0\) holds because a call on a zero-price asset is worthless. The condition at \(S = S_{\max}\) holds because \(S_{\max}\) is chosen sufficiently far above the upper barrier (or upper payoff cutoff) that the option value is negligible.

Because both boundary values are zero, the tridiagonal system for interior unknowns \(j = 1, \ldots, M{-}1\) requires **no modification** at the first or last row. For products with nonzero boundary values \(g_0^{n+1} = V(0, t_{n+1})\) and \(g_M^{n+1} = V(S_{\max}, t_{n+1})\), the right-hand side at \(j = 1\) and \(j = M{-}1\) must be adjusted by subtracting the coupling to the known boundary value; this case does not arise for the products specified below.

### 3(d). Initial / Payoff Conditions

The initial condition is set at \(t = 0\) (maturity):

**Product Type 1 — Truncated call:**

$$
U_j^0 = \begin{cases} S_j - K & \text{if } S_j \in [K,\, B_u], \\ 0 & \text{otherwise}, \end{cases} \qquad j = 0, 1, \ldots, M
$$

This payoff has a kink at \(S = K\) and a jump discontinuity at \(S = B_u\) (from \(B_u - K\) to \(0\)). There is no lower barrier; set \(B_\ell = 0\) (unused). No monitoring dates apply.

**Product Type 2 — Discrete double barrier knock-out call:**

$$
U_j^0 = \max(S_j - K,\; 0)\;\mathbf{1}_{[B_\ell,\, B_u]}(S_j), \qquad j = 0, 1, \ldots, M
$$

This payoff has discontinuities at \(S = B_\ell\) (if \(B_\ell > K\)) and \(S = B_u\). Monitoring dates \(\{t_i^{\mathrm{mon}}\}_{i=1}^F\) must be specified. The time grid \(\Delta t\) must be chosen so that every monitoring date coincides with a grid point \(t_n\) for some \(n\).

### 3(e). Spatial Discretization

**Grid construction (uniform):**

$$
S_j = j\,\Delta S, \qquad j = 0, 1, \ldots, M, \qquad M = S_{\max}/\Delta S \in \mathbb{Z}^+
$$

The parameters \(S_{\max}\) and \(\Delta S\) must be chosen so that \(M\) is a positive integer.

---

#### Scheme 1: Exponentially Fitted Fully Implicit

**Fitting factor** at each interior node \(j = 1, \ldots, M{-}1\):

$$
\mathrm{Pe}_j = \frac{r\,\Delta S}{\sigma^2 S_j} = \frac{r}{\sigma^2 j}
$$

$$
\rho_j = \sigma_d(S_j) \cdot \mathrm{Pe}_j\,\coth(\mathrm{Pe}_j) = \frac{rS_j\,\Delta S}{2}\,\coth\!\left(\frac{r\,\Delta S}{\sigma^2 S_j}\right)
$$

**Numerical guard for \(\coth\) evaluation:** If \(|\mathrm{Pe}_j| < 10^{-6}\): use \(\mathrm{Pe}\,\coth(\mathrm{Pe}) \approx 1 + \mathrm{Pe}^2/3\), giving \(\rho_j \approx \sigma_d(S_j)\bigl(1 + \mathrm{Pe}_j^2/3\bigr)\). If \(|\mathrm{Pe}_j| > 20\): use \(\coth(\mathrm{Pe}_j) \approx \mathrm{sign}(\mathrm{Pe}_j)\), giving \(\rho_j \approx rS_j\,\Delta S/2\). Otherwise: compute \(\coth\) directly as \(\coth(x) = (e^{2x}+1)/(e^{2x}-1)\).

**Tridiagonal matrix** \(A = \mathrm{tridiag}\{a_j^-;\; a_j^0;\; a_j^+\}\), with entries for \(j = 1, \ldots, M{-}1\):

$$
a_j^- = \left(-\frac{\rho_j}{\Delta S^2} + \frac{rS_j}{2\,\Delta S}\right)\Delta t
$$

$$
a_j^0 = 1 + r\,\Delta t + \frac{2\rho_j}{\Delta S^2}\,\Delta t
$$

$$
a_j^+ = -\left(\frac{\rho_j}{\Delta S^2} + \frac{rS_j}{2\,\Delta S}\right)\Delta t
$$

**Row sum** (constant across all rows):

$$
a_j^- + a_j^0 + a_j^+ = 1 + r\,\Delta t
$$

For constant \(r, \sigma\), the fitting factor \(\rho_j\) and all matrix entries depend only on the node index \(j\), not on the time level. The matrix \(A\) is assembled once and reused at every time step.

---

#### Scheme 2: Crank–Nicolson Variant

**Weight parameter:**

$$
\omega = -\frac{r}{16\sigma^2}
$$

**Left-hand matrix** \(P = \mathrm{tridiag}\{p_j^-;\; p_j^0;\; p_j^+\}\), for \(j = 1, \ldots, M{-}1\):

$$
p_j^- = -\left(\frac{\sigma S_j}{2\,\Delta S} - \frac{r}{4\sigma}\right)^{\!2}
$$

$$
p_j^0 = \frac{1}{\Delta t} + \frac{\sigma^2 S_j^2}{2\,\Delta S^2} + \frac{r}{2} + \frac{r^2}{8\sigma^2}
$$

$$
p_j^+ = -\left(\frac{\sigma S_j}{2\,\Delta S} + \frac{r}{4\sigma}\right)^{\!2}
$$

**Right-hand matrix** \(N = \mathrm{tridiag}\{n_j^-;\; n_j^0;\; n_j^+\}\), for \(j = 1, \ldots, M{-}1\):

$$
n_j^- = +\left(\frac{\sigma S_j}{2\,\Delta S} - \frac{r}{4\sigma}\right)^{\!2}
$$

$$
n_j^0 = \frac{1}{\Delta t} - \frac{\sigma^2 S_j^2}{2\,\Delta S^2} - \frac{r}{2} - \frac{r^2}{8\sigma^2}
$$

$$
n_j^+ = +\left(\frac{\sigma S_j}{2\,\Delta S} + \frac{r}{4\sigma}\right)^{\!2}
$$

**Equivalence to general-form entries.** The perfect-square forms above are obtained by substituting \(\omega = -r/(16\sigma^2)\) into the general entries:

$$
p_j^- = r\omega + \frac{rS_j}{4\,\Delta S} - \frac{\sigma^2 S_j^2}{4\,\Delta S^2}, \qquad n_j^- = -r\omega - \frac{rS_j}{4\,\Delta S} + \frac{\sigma^2 S_j^2}{4\,\Delta S^2}
$$

$$
p_j^0 = \frac{1}{\Delta t} + \frac{\sigma^2 S_j^2}{2\,\Delta S^2} + r\!\left(\frac{1}{2} - 2\omega\right), \qquad n_j^0 = \frac{1}{\Delta t} - \frac{\sigma^2 S_j^2}{2\,\Delta S^2} - r\!\left(\frac{1}{2} - 2\omega\right)
$$

$$
p_j^+ = r\omega - \frac{rS_j}{4\,\Delta S} - \frac{\sigma^2 S_j^2}{4\,\Delta S^2}, \qquad n_j^+ = -r\omega + \frac{rS_j}{4\,\Delta S} + \frac{\sigma^2 S_j^2}{4\,\Delta S^2}
$$

**Row sums** (constant across all rows):

$$
p_j^- + p_j^0 + p_j^+ = \frac{1}{\Delta t} + \frac{r}{2}, \qquad n_j^- + n_j^0 + n_j^+ = \frac{1}{\Delta t} - \frac{r}{2}
$$

For constant \(r, \sigma\), both \(P\) and \(N\) are assembled once and reused.

---

### 3(f). Temporal Discretization

Time marches forward from \(n = 0\) (maturity) to \(n = N_t\) (present). Each step advances from \(t_n\) to \(t_{n+1} = t_n + \Delta t\).

**Scheme 1 — Fully implicit update:**

$$
A\,\mathbf{U}^{n+1} = \mathbf{U}^n
$$

All spatial terms (convection, diffusion, reaction) are evaluated at the new time level \(n{+}1\). The time discretization is backward Euler, which is first-order in \(\Delta t\).

**Scheme 2 — Crank–Nicolson variant update:**

$$
P\,\mathbf{U}^{n+1} = N\,\mathbf{U}^n
$$

Convection and diffusion are averaged between levels \(n\) and \(n{+}1\) (standard Crank–Nicolson). The reaction term \(-rV\) is discretized using the six-node stencil with weight \(\omega\) at each level. The time discretization is second-order in \(\Delta t\).

---

### 3(g). Linear System Solution

Both schemes produce a tridiagonal system of dimension \((M{-}1)\):

$$
a_j^-\,x_{j-1} + a_j^0\,x_j + a_j^+\,x_{j+1} = d_j, \qquad j = 1, \ldots, M{-}1
$$

with \(x_0 = x_M = 0\) (homogeneous Dirichlet). This is solved by the **Thomas algorithm** (LU factorization for tridiagonal systems). Because both \(A\) (Scheme 1) and \(P\) (Scheme 2) are M-matrices under the stated conditions, the Thomas algorithm is stable without pivoting.

**Thomas algorithm (forward sweep and back substitution):**

For \(j = 1\): set \(\tilde{a}_1^0 = a_1^0\), \(\tilde{d}_1 = d_1\).

For \(j = 2, 3, \ldots, M{-}1\): compute the multiplier \(w = a_j^- / \tilde{a}_{j-1}^0\), then set \(\tilde{a}_j^0 = a_j^0 - w\,a_{j-1}^+\) and \(\tilde{d}_j = d_j - w\,\tilde{d}_{j-1}\).

Back substitution: \(x_{M-1} = \tilde{d}_{M-1}/\tilde{a}_{M-1}^0\). For \(j = M{-}2, M{-}3, \ldots, 1\): \(x_j = (\tilde{d}_j - a_j^+\,x_{j+1})/\tilde{a}_j^0\).

---

### 3(h). Special Treatment

**Fitting factor computation (Scheme 1):** See the numerical guard in §3(e). The fitting factor \(\rho_j\) is undefined at \(j = 0\) (where \(S_0 = 0\)), but this node is a boundary node and does not enter the interior system.

**Monitoring-date projection (Product Type 2 only):** At each monitoring date \(t_i^{\mathrm{mon}}\) that coincides with grid time \(t_n\), after computing \(\mathbf{U}^n\) via the linear solve, apply:

$$
U_j^n \;\leftarrow\; U_j^n \cdot \mathbf{1}_{[B_\ell,\, B_u]}(S_j), \qquad j = 0, 1, \ldots, M
$$

This is a pointwise post-processing step that sets to zero the solution at all nodes outside the barrier corridor. It re-introduces discontinuities at \(S = B_\ell\) and \(S = B_u\), which the subsequent time steps must resolve. No modification to the tridiagonal matrices is needed.

**Requirement on time grid:** Every monitoring date must coincide with a grid point. Choose \(\Delta t\) so that \(t_i^{\mathrm{mon}}/\Delta t\) is an integer for all \(i = 1, \ldots, F\).

---

### 3(i). Output Extraction

**Price:** After all \(N_t\) time steps, the vector \(\mathbf{U}^{N_t}\) contains the option price at the present time (\(t = T\)) for each interior grid point. The price at a specific asset value \(S^*\) is obtained by interpolation (linear or cubic) from the grid values \(\{(S_j, U_j^{N_t})\}\).

**Greeks (not specified in the paper; standard finite-difference formulas):**

$$
\Delta_j \approx \frac{U_{j+1}^{N_t} - U_{j-1}^{N_t}}{2\,\Delta S}, \qquad \Gamma_j \approx \frac{U_{j+1}^{N_t} - 2U_j^{N_t} + U_{j-1}^{N_t}}{\Delta S^2}
$$

$$
\Theta_j \approx -\frac{U_j^{N_t} - U_j^{N_t - 1}}{\Delta t}
$$

The sign convention for \(\Theta\) reflects the backward-time variable: \(\Theta = -\partial V/\partial t\) in calendar time. These formulas are not part of the original paper and are provided here as standard supplementary material for the implementer.

---

## 4. CONDITIONS AND CONSTRAINTS

### 4.1. Well-definedness (matrix invertibility)

**Scheme 1:** The matrix \(A\) is a nonsingular M-matrix for all \(r > 0\), \(\sigma > 0\), \(\Delta S > 0\), \(\Delta t > 0\). No constraint on grid parameters is needed. The Thomas algorithm completes without zero pivots.

**Scheme 2:** The matrix \(P\) is a nonsingular M-matrix unconditionally (for the stated \(\omega\) choice, any \(\Delta t > 0\), and \(r > 0\)). The system \(P\mathbf{U}^{n+1} = N\mathbf{U}^n\) is always solvable. Note: \(P\) becomes reducible (zero sub-diagonal entry) at a single row if the grid contains the critical point \(S^* = r\Delta S/(2\sigma^2)\); this falls outside the computational domain for all low-volatility parameter sets in the paper. If \(S^* \leq S_{\max}\), perturb \(\omega\) by replacing it with \(\omega - \epsilon\) for a small \(\epsilon > 0\) (e.g., \(\epsilon = 10^{-12}\)).

### 4.2. Positivity preservation

**Scheme 1:** Unconditional. If \(\mathbf{U}^0 \geq 0\) and \(\mathbf{U}^0 \neq 0\), then \(\mathbf{U}^n > 0\) (strictly) for all \(n \geq 1\). Proof: \(A^{-1} > 0\) (all entries strictly positive), so the product of a strictly positive matrix with a nonnegative nonzero vector is strictly positive.

**Scheme 2:** Guaranteed under the time-step constraint:

$$
\Delta t \;<\; \Delta t_{\max} \;\equiv\; \frac{1}{\displaystyle \frac{r}{2} + \frac{r^2}{8\sigma^2} + \frac{\sigma^2(M\!-\!1)^2}{2}}
$$

where the denominator is evaluated at the worst-case interior node \(j = M{-}1\). Under this constraint, \(N \geq 0\) entrywise and \(P^{-1} \geq 0\), so \(P^{-1}N \geq 0\) and positivity is preserved by induction.

When this constraint is violated (as in many of the paper's own examples), positivity is not formally guaranteed but is typically observed in practice because the diagonal of \(N\) becomes negative only near \(S_{\max}\), where the option value is essentially zero.

### 4.3. Discrete maximum principle (non-amplification)

**Scheme 1:** Unconditional. \(\|\mathbf{U}^{n+1}\|_\infty \leq \|A^{-1}\|_\infty\|\mathbf{U}^n\|_\infty \leq \frac{1}{1+r\Delta t}\|\mathbf{U}^n\|_\infty < \|\mathbf{U}^n\|_\infty\).

**Scheme 2:** Under the time-step constraint of §4.2:

$$
\|\mathbf{U}^{n+1}\|_\infty \;\leq\; \frac{1/\Delta t - r/2}{1/\Delta t + r/2}\;\|\mathbf{U}^n\|_\infty \;<\; \|\mathbf{U}^n\|_\infty
$$

### 4.4. Convergence order

**Scheme 1:** Global error \(|V(S_j, t_n) - U_j^n| \leq c(\Delta S + \Delta t)\), where \(c\) is independent of \(\Delta S\), \(\Delta t\), and \(\sigma\). This is first-order convergence, uniform in the volatility parameter. For fixed \(\sigma > 0\), the effective spatial order improves to \(O(\Delta S^2)\), giving \(O(\Delta S^2 + \Delta t)\) overall.

**Scheme 2:** Formal truncation error \(O(\Delta S^2 + \Delta t^2)\). However, the coefficient of the \(\Delta S^2\) spatial error term is \(r^2/(8\sigma^2)\), which is large for small \(\sigma\). For example, with \(r = 0.5\) and \(\sigma = 0.001\), this coefficient is \(31{,}250\). The formal second-order rate is achieved only when \(\Delta S\) is small enough that \(\tfrac{1}{8}(r\Delta S/\sigma)^2\) is negligible compared to the physical diffusion.

### 4.5. Artificial diffusion (diagnostic)

**Scheme 1** (low-\(\sigma\) limit): The scheme effectively adds artificial diffusion \(\tfrac{1}{2}rS\,\Delta S\,\partial^2 V/\partial S^2\), which is \(O(\Delta S)\) and dominates when \(\sigma^2 S_j \ll r\Delta S\).

**Scheme 2**: The scheme effectively adds artificial diffusion \(\tfrac{1}{8}(r\Delta S/\sigma)^2\,\partial^2 V/\partial S^2\), which is \(O(\Delta S^2)\) with coefficient \(r^2/(8\sigma^2)\), and dominates when \(\sigma^2 S \ll r\Delta S/2\).

### 4.6. Summary of parameter requirements

| Requirement | Scheme 1 | Scheme 2 |
|:---|:---|:---|
| \(r > 0\) | Required | Required |
| \(\sigma > 0\) | Required | Required |
| \(\Delta S > 0\), \(\Delta t > 0\) | Any | Any for solvability; \(\Delta t < \Delta t_{\max}\) for positivity |
| \(S_{\max} \geq 2B_u\) | Required | Required |
| \(M = S_{\max}/\Delta S \in \mathbb{Z}^+\) | Required | Required |
| Monitoring dates on grid | Required (if applicable) | Required (if applicable) |

---

## 5. MONOLITHIC PSEUDOCODE

```text
001  //=======================================================================
002  // SUBROUTINE: ThomasSolve
003  // Solves tridiagonal system with m unknowns and zero boundary values.
004  // Per §3(g).
005  //=======================================================================
006  FUNCTION ThomasSolve(
007      a_sub  : float[1..m],   // sub-diagonal; a_sub[1] unused
008      a_diag : float[1..m],   // diagonal
009      a_sup  : float[1..m],   // super-diagonal; a_sup[m] unused
010      rhs    : float[1..m],   // right-hand side
011      m      : int            // system dimension
012  ) → float[1..m]
013  BEGIN
014      DECLARE x     : float[1..m]
015      DECLARE d_mod : float[1..m]   // modified diagonal
016      DECLARE r_mod : float[1..m]   // modified RHS
017
018      // Forward sweep
019      d_mod[1] ← a_diag[1]
020      r_mod[1] ← rhs[1]
021      FOR j = 2 TO m DO
022          w ← a_sub[j] / d_mod[j-1]
023          d_mod[j] ← a_diag[j] - w * a_sup[j-1]
024          r_mod[j] ← rhs[j]   - w * r_mod[j-1]
025      END FOR
026
027      // Back substitution
028      x[m] ← r_mod[m] / d_mod[m]
029      FOR j = m-1 DOWNTO 1 DO
030          x[j] ← (r_mod[j] - a_sup[j] * x[j+1]) / d_mod[j]
031      END FOR
032
033      RETURN x
034  END FUNCTION
035
036  //=======================================================================
037  // SUBROUTINE: TridiagMatVec
038  // Computes y = T * x for tridiagonal T, with x[0] = x[m+1] = 0.
039  // Per §3(f), Scheme 2 RHS computation.
040  //=======================================================================
041  FUNCTION TridiagMatVec(
042      t_sub  : float[1..m],
043      t_diag : float[1..m],
044      t_sup  : float[1..m],
045      x      : float[1..m],
046      m      : int
047  ) → float[1..m]
048  BEGIN
049      DECLARE y : float[1..m]
050      y[1] ← t_diag[1]*x[1] + t_sup[1]*x[2]
051      FOR j = 2 TO m-1 DO
052          y[j] ← t_sub[j]*x[j-1] + t_diag[j]*x[j] + t_sup[j]*x[j+1]
053      END FOR
054      y[m] ← t_sub[m]*x[m-1] + t_diag[m]*x[m]
055      RETURN y
056  END FUNCTION
057
058  //=======================================================================
059  // SUBROUTINE: ComputeRho
060  // Computes fitting factor with numerical guards.  Per §3(e), Scheme 1.
061  //=======================================================================
062  FUNCTION ComputeRho(
063      mu_j     : float,    // = r * S_j
064      sigma_d_j: float,    // = 0.5 * sigma^2 * S_j^2
065      DeltaS   : float
066  ) → float
067  BEGIN
068      Pe ← (mu_j * DeltaS) / (2.0 * sigma_d_j)       // Péclet number
069      IF |Pe| < 1.0e-6 THEN
070          RETURN sigma_d_j * (1.0 + Pe*Pe / 3.0)      // Taylor approx
071      ELSE IF |Pe| > 20.0 THEN
072          RETURN |mu_j| * DeltaS / 2.0                 // upwind limit
073      ELSE
074          e2p ← EXP(2.0 * Pe)
075          coth_Pe ← (e2p + 1.0) / (e2p - 1.0)
076          RETURN (mu_j * DeltaS / 2.0) * coth_Pe
077      END IF
078  END FUNCTION
079
080  //=======================================================================
081  // MAIN PROCEDURE: PriceOption
082  //=======================================================================
083  PROCEDURE PriceOption(
084      // --- Market parameters ---
085      r       : float,        // risk-free rate (yr^-1), must be > 0
086      sigma   : float,        // volatility (yr^-1/2), must be > 0
087      // --- Contract parameters ---
088      K       : float,        // strike price ($)
089      B_lo    : float,        // lower barrier ($); 0 if unused
090      B_up    : float,        // upper barrier ($)
091      T       : float,        // time to expiry (yr)
092      product : {"TRUNCATED_CALL", "DISCRETE_BARRIER_CALL"},
093      // --- Monitoring dates (empty for TRUNCATED_CALL) ---
094      mon_dates : float[],    // sorted monitoring dates in time-to-expiry
095      // --- Grid parameters ---
096      S_max   : float,        // domain truncation ($), require >= 2*B_up
097      DeltaS  : float,        // spatial step ($)
098      Deltat  : float,        // time step (yr)
099      // --- Scheme choice ---
100      scheme  : {"EXP_FITTED", "CN_VARIANT"}
101  ) → float[0..M]   // price at each grid node at t = T
102  BEGIN
103      //=== GRID CONSTRUCTION  (per §3(e)) ===
104      M  ← ROUND(S_max / DeltaS)                       // int
105      Nt ← ROUND(T / Deltat)                            // int
106      m  ← M - 1                                        // interior unknowns
107
108      DECLARE S : float[0..M]
109      FOR j = 0 TO M DO
110          S[j] ← j * DeltaS                             // per §3(e)
111      END FOR
112
113      //=== PAYOFF INITIALIZATION  (per §3(d)) ===
114      DECLARE U : float[1..m]
115      FOR j = 1 TO m DO
116          IF product = "TRUNCATED_CALL" THEN
117              IF S[j] >= K AND S[j] <= B_up THEN
118                  U[j] ← S[j] - K
119              ELSE
120                  U[j] ← 0.0
121              END IF
122          ELSE IF product = "DISCRETE_BARRIER_CALL" THEN
123              call_val ← MAX(S[j] - K, 0.0)
124              IF S[j] >= B_lo AND S[j] <= B_up THEN
125                  U[j] ← call_val
126              ELSE
127                  U[j] ← 0.0
128              END IF
129          END IF
130      END FOR
131
132      //=== MATRIX ASSEMBLY  (per §3(e)) ===
133      DECLARE a_sub, a_diag, a_sup : float[1..m]
134
135      IF scheme = "EXP_FITTED" THEN
136          // Scheme 1: assemble A once (constant coefficients)
137          FOR j = 1 TO m DO
138              mu_j     ← r * S[j]
139              sigma_d_j← 0.5 * sigma*sigma * S[j]*S[j]
140              rho_j    ← ComputeRho(mu_j, sigma_d_j, DeltaS)
141
142              a_sub[j] ← (-rho_j/(DeltaS*DeltaS)
143                           + mu_j/(2.0*DeltaS)) * Deltat
144              a_diag[j]← 1.0 + r*Deltat
145                           + 2.0*rho_j*Deltat/(DeltaS*DeltaS)
146              a_sup[j] ← -(rho_j/(DeltaS*DeltaS)
147                           + mu_j/(2.0*DeltaS)) * Deltat
148          END FOR
149
150      ELSE IF scheme = "CN_VARIANT" THEN
151          // Scheme 2: assemble P and N once (constant coefficients)
152          omega ← -r / (16.0 * sigma*sigma)             // per §3(e)
153          c1    ← r / (4.0 * sigma)                     // = r/(4σ)
154
155          DECLARE p_sub, p_diag, p_sup : float[1..m]
156          DECLARE n_sub, n_diag, n_sup : float[1..m]
157
158          FOR j = 1 TO m DO
159              alpha ← sigma * S[j] / (2.0 * DeltaS)    // = σS_j/(2ΔS)
160              term_minus ← (alpha - c1)*(alpha - c1)    // perfect square
161              term_plus  ← (alpha + c1)*(alpha + c1)    // perfect square
162
163              // P matrix entries
164              p_sub[j]  ← -term_minus
165              p_diag[j] ← 1.0/Deltat + 2.0*alpha*alpha
166                           + r/2.0 + r*r/(8.0*sigma*sigma)
167              p_sup[j]  ← -term_plus
168
169              // N matrix entries
170              n_sub[j]  ← +term_minus
171              n_diag[j] ← 1.0/Deltat - 2.0*alpha*alpha
172                           - r/2.0 - r*r/(8.0*sigma*sigma)
173              n_sup[j]  ← +term_plus
174          END FOR
175
176          // Copy P into the a_sub/a_diag/a_sup arrays for Thomas solve
177          FOR j = 1 TO m DO
178              a_sub[j]  ← p_sub[j]
179              a_diag[j] ← p_diag[j]
180              a_sup[j]  ← p_sup[j]
181          END FOR
182      END IF
183
184      //=== BUILD SET OF MONITORING STEP INDICES  (per §3(h)) ===
185      DECLARE mon_steps : SET OF int ← {}
186      IF product = "DISCRETE_BARRIER_CALL" THEN
187          FOR EACH t_mon IN mon_dates DO
188              n_mon ← ROUND(t_mon / Deltat)
189              INSERT n_mon INTO mon_steps
190          END FOR
191      END IF
192
193      //=== TIME-STEPPING LOOP  (per §3(f)) ===
194      FOR n = 0 TO Nt - 1 DO
195
196          IF scheme = "EXP_FITTED" THEN
197              // Scheme 1: solve A * U_new = U_old
198              U ← ThomasSolve(a_sub, a_diag, a_sup, U, m)
199
200          ELSE IF scheme = "CN_VARIANT" THEN
201              // Scheme 2: compute RHS = N * U_old, then solve P * U_new = RHS
202              rhs ← TridiagMatVec(n_sub, n_diag, n_sup, U, m)
203              U   ← ThomasSolve(a_sub, a_diag, a_sup, rhs, m)
204          END IF
205
206          // Monitoring-date projection  (per §3(h))
207          IF (n + 1) IN mon_steps THEN
208              FOR j = 1 TO m DO
209                  IF S[j] < B_lo OR S[j] > B_up THEN
210                      U[j] ← 0.0
211                  END IF
212              END FOR
213          END IF
214
215      END FOR  // end time loop
216
217      //=== ASSEMBLE FULL OUTPUT VECTOR (per §3(i)) ===
218      DECLARE V_out : float[0..M]
219      V_out[0] ← 0.0                                    // boundary
220      FOR j = 1 TO m DO
221          V_out[j] ← U[j]
222      END FOR
223      V_out[M] ← 0.0                                    // boundary
224
225      RETURN V_out
226  END PROCEDURE
```

---

## 6. ERRATA RELATIVE TO ORIGINAL PAPER

| # | Location in paper | Original | Corrected | Reason |
|:---|:---|:---|:---|:---|
| 1 | §2, U7: tridiagonal matrix \(A\) | Coefficients carry superscript \(n\): \(\rho_j^n, \mu_j^n, b_j^n\) | Superscript should be \(n{+}1\): \(\rho_j^{n+1}, \mu_j^{n+1}, b_j^{n+1}\) | The fitted operator U6 evaluates all spatial terms at the implicit level \(n{+}1\). For the constant-coefficient Black–Scholes case the distinction is immaterial; for variable coefficients it changes the scheme. |
| 2 | §2, U11: maximum principle chain | \(\|A^{-1}U^n\|_\infty = \|A^{-1}\|_\infty\|U^n\|_\infty\) (equality) | \(\|A^{-1}U^n\|_\infty \leq \|A^{-1}\|_\infty\|U^n\|_\infty\) (inequality) | By submultiplicativity of the operator norm; equality holds only for specific vectors. The conclusion \(\|\mathbf{U}^{n+1}\|_\infty \leq \|\mathbf{U}^n\|_\infty\) is unaffected. |
| 3 | §2, U13: upwind scheme (\(\mu > 0\)) | Denominator \(2h\) in convection term: \(\mu_j^{n+1}\frac{U_{j+1}^{n+1}-U_j^{n+1}}{2h}\) | Denominator \(h\): \(\mu_j^{n+1}\frac{U_{j+1}^{n+1}-U_j^{n+1}}{h}\) | Substituting \(\rho = \mu h/2\) into the centered + artificial-diffusion combination yields \(\mu(U_{j+1}-U_j)/h\), not \(\mu(U_{j+1}-U_j)/(2h)\). The paper's own numerical-diffusion formula U15 (\(\frac{1}{2}\mu h V''\)) is consistent with denominator \(h\), confirming the typo. |
| 4 | §2, U14: upwind scheme (\(\mu < 0\)) | Denominator \(2h\): \(\mu_j^{n+1}\frac{U_j^{n+1}-U_{j-1}^{n+1}}{2h}\) | Denominator \(h\): \(\mu_j^{n+1}\frac{U_j^{n+1}-U_{j-1}^{n+1}}{h}\) | Same derivation as erratum #3 applied to the \(\mu < 0\) case. |
| 5 | §3, U24: CN variant maximum principle chain | \(\|(P^{-1}N)U^n\|_\infty = \|P^{-1}\|_\infty\|N\|_\infty\|U^n\|_\infty\) (equality) | \(\|(P^{-1}N)U^n\|_\infty \leq \|P^{-1}\|_\infty\|N\|_\infty\|U^n\|_\infty\) (inequality) | Same submultiplicativity issue as erratum #2. Conclusion unaffected. |
| 6 | §2, U10/C6: positivity claim | \((A^{-1})^n U^0 > 0\) for \(U^0 \geq 0\) | \((A^{-1})^n U^0 > 0\) requires \(U^0 \geq 0\) **and** \(U^0 \neq 0\); for \(U^0 = 0\) the result is \(= 0\), not \(> 0\) | Strict positivity of \(A^{-1}U^0\) requires at least one positive component in \(U^0\). For option pricing payoffs this is always satisfied. |
| 7 | §§2–4: grid construction, boundary conditions, \(S_{\max}\), \(M\) | Not formally specified in the paper body; \(S_{\max}\) appears only in figure captions; boundary treatment at \(j=0\) and \(j=M\) is never described; \(M\) is used without prior definition. | Uniform grid \(S_j = j\Delta S\), \(M = S_{\max}/\Delta S\); homogeneous Dirichlet BCs \(U_0^n = U_M^n = 0\); \(S_{\max} \geq 2B_u\). See §3(c) and §3(e) of this document. | Inferred from the tridiagonal structure, figure captions, and the companion paper (Tagliani & Milev, *Appl. Math. Comput.* 2013). Essential for implementation. |
| 8 | §4, Examples 4.1/4.4/4.5: time-step constraint E(9) | Examples use \(\Delta t = 0.01\); paper does not note any violation of E(9). | E(9) is violated in Examples 4.1 (\(r=0.05\): needs \(\Delta t < 0.0032\)), 4.4 (\(r=0.5\): needs \(\Delta t < 3.2 \times 10^{-5}\)), and 4.5. Only Example 4.4 with \(r=0.01\) satisfies E(9). Theoretical positivity and maximum-principle guarantees do not formally apply to the violating examples. | Arithmetic evaluation of E(9) for each parameter set. The scheme still produces acceptable results because E(9) is a sufficient condition and the violated constraint affects only nodes near \(S_{\max}\) where the solution is near zero. |
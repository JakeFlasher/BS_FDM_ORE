## 1. Paper Identity

- **Title:** *Qualitatively Stable Nonstandard Finite Difference Scheme for Numerical Solution of the Nonlinear Black–Scholes Equation*
- **Authors:** Mohammad Mehdizadeh Khalsaraei; Ali Shokri; Zahra Mohammadnia; Hamid Mohammad Sedighi
- **Affiliations:**
  - Department of Mathematics, Faculty of Science, University of Maragheh, Maragheh, Iran
  - Mechanical Engineering Department, Faculty of Engineering, Shahid Chamran University of Ahvaz, Ahvaz, Iran
- **Venue / publication model:** *Journal of Mathematics* (Hindawi; open access)
- **Year / dates:** 2021 (Received 19 Dec 2020; Accepted 26 Apr 2021; Published 11 May 2021)
- **DOI:** `10.1155/2021/6679484`
- **arXiv ID:** Not reported in the provided source.

## 2. Problem Statement & Formulation

European option prices under proportional transaction costs are modeled by a **nonlinear** Black–Scholes-type parabolic PDE in which the diffusion coefficient (effective volatility) depends on the option’s Gamma $$V_{SS}$$. The paper targets a grid-based method that (i) preserves nonnegativity of the computed option price, (ii) is monotone and stable without a Courant-type restriction, and (iii) converges (as the mesh is refined) to the viscosity solution.

### 2.1 PDE in original and transformed time

The classical (linear) Black–Scholes PDE is stated as
$$
U_t + \frac{1}{2}\sigma_0^2 S^2 U_{SS} + r S U_S - rU = 0,
$$
with constant volatility $$\sigma_0$$ and risk-free rate $$r$$.

With transaction costs and Gamma-dependent volatility, the nonlinear form is given as
$$
U_t + \frac{1}{2}\sigma^2\!\left(t,S,U_S,U_{SS}\right) S^2 U_{SS} + r S U_S - rU = 0.
$$

The paper applies the time-to-maturity change of variables
$$
\tau = T - t,
$$
and uses the notation $$V(S,\tau) = U(S, T-\tau)$$, yielding
$$
V_\tau = \frac{1}{2}\sigma^2\!\left(V_{SS}\right) S^2 V_{SS} + r S V_S - rV,
\qquad S>0,\quad 0<\tau \le T.
$$

### 2.2 Barles–Soner-type modified volatility and the function $$\Psi$$

The modified volatility is written in the paper in the Barles–Soner form
$$
\sigma^2 = \sigma_0^2\left(1+\Psi(\,\cdot\,)\right),
$$
where $$\Psi$$ is defined as the solution of the nonlinear IVP
$$
\Psi'(w)=\Psi(w)+\frac{1}{2}\frac{(w\Psi(w))^2}{w\Psi(w)-w},
\qquad w\ne 0,\qquad \Psi(0)=0.
$$

An implicit closed-form characterization is also provided in the source (piecewise in $$w>0$$ and $$w<0$$) using inverse hyperbolic/trigonometric functions.

**Ambiguity flag (verbatim-source quality):** In the provided text, the exact argument passed to $$\Psi(\cdot)$$ inside $$\sigma^2$$ is typographically corrupted (placement of factors involving $$e^{r\tau}$$, $$a^2$$, and $$S^2$$, and whether an absolute value is applied to $$V_{SS}$$). The later monotonicity proof treats expressions of the form $$\bigl(1+\Psi(Kv)\bigr)v$$ with $$K>0$$ and $$v$$ representing a discrete Gamma; this is sufficient to reproduce the scheme analysis, but the precise continuous-model argument of $$\Psi$$ cannot be reconstructed with complete certainty from the pasted text alone.

### 2.3 Spatial domain truncation and payoff/boundary conditions

The computational domain is truncated to
$$
\Omega = [0,S_{\max}] \times [0,T],
$$
with initial condition at $$\tau=0$$ (maturity payoff) and boundary conditions at $$S=0$$ and $$S=S_{\max}$$.

Define payoff functions $$L_1(S)$$, left boundary $$L_2(\tau)$$, and right boundary $$L_3(\tau)$$ as listed in the paper:

- **European call**
  $$
  L_1(S)=\max(S-K,0),\quad
  L_2(\tau)=0,\quad
  L_3(\tau)=S_{\max}-K e^{-r\tau}.
  $$

- **European put**
  $$
  L_1(S)=\max(K-S,0),\quad
  L_2(\tau)=K e^{-r\tau},\quad
  L_3(\tau)=0.
  $$

- **Butterfly spread**
  $$
  L_1(S)=\max(S-K_1,0)-2\max(S-K_2,0)+\max(S-K_3,0),\quad
  L_2(\tau)=0,\quad
  L_3(\tau)=0.
  $$

- **Cash-or-nothing (digital)**
  $$
  L_1(S)=B\,H(S-K),\quad
  L_2(\tau)=0,\quad
  L_3(\tau)=B e^{-r\tau},
  $$
  where $$H$$ is the Heaviside function and $$B>0$$ is a payoff level.

The paper states compatibility conditions
$$
L_1(0)=L_2(0),\qquad L_1(S_{\max})=L_3(0).
$$

## 3. Core Methodology

The method is an **explicit** (in the sense of a closed-form update) **nonstandard finite difference (NSFD)** time-marching scheme in which spatial derivatives are discretized in a nonstandard, “nonlocalized” way (mixing time levels inside spatial stencils) to enforce monotonicity/positivity/stability properties.

### 3.1 Grid and notation

Let
- $$S_j=jh$$ for $$j=0,1,\dots,M$$ with $$h=S_{\max}/M$$,
- $$\tau_n=n\Delta\tau$$ for $$n=0,1,\dots,N$$ with $$\Delta\tau=T/N$$,
- $$V_j^n \approx V(S_j,\tau_n)$$.

Initial and boundary data are discretized as
$$
V_j^0 = L_1(S_j),\qquad
V_0^n = L_2(\tau_n),\qquad
V_M^n = L_3(\tau_n).
$$

### 3.2 Proposed NSFD discretization and explicit update

For interior nodes $$j=1,\dots,M-1$$, the scheme is written (paper’s equation (17)) as
$$
-\frac{V_j^{n+1}-V_j^n}{\Delta\tau}
+\frac{1}{2}\sigma^2\!\left(V_{SS}^n\right)\!(j)\,S_j^2\,
\frac{V_{j-1}^n-2V_j^{n+1}+V_{j+1}^n}{h^2}
+rS_j\,\frac{V_{j+1}^n-V_j^{n+1}}{h}
-rV_j^{n+1}=0.
$$

Solving for $$V_j^{n+1}$$ yields the explicit closed-form update (paper’s equation (18)):
$$
V_j^{n+1}
=
\frac{
\left(\frac{1}{2}\sigma^2\!\left(V_{SS}^n\right)\!(j)\frac{S_j^2}{h^2}\right)V_{j-1}^n
+\left(\frac{1}{\Delta\tau}\right)V_j^n
+\left(\frac{rS_j}{h}+\frac{1}{2}\sigma^2\!\left(V_{SS}^n\right)\!(j)\frac{S_j^2}{h^2}\right)V_{j+1}^n
}{
\left(\frac{1}{\Delta\tau}\right)+\frac{rS_j}{h}+\sigma^2\!\left(V_{SS}^n\right)\!(j)\frac{S_j^2}{h^2}+r
}.
$$

Equivalently, the method is written in convex-combination form (paper’s equation (22)):
$$
V_j^{n+1} = x_j^n\,V_{j-1}^n + y_j^n\,V_j^n + z_j^n\,V_{j+1}^n,
$$
with coefficients (paper’s equations (23)–(25)):
$$
x_j^n
=
\frac{
\frac{1}{2}\sigma^2\!\left(V_{SS}^n\right)\!(j)\frac{S_j^2}{h^2}
}{
\left(\frac{1}{\Delta\tau}\right)+\frac{rS_j}{h}+\sigma^2\!\left(V_{SS}^n\right)\!(j)\frac{S_j^2}{h^2}+r
},
$$
$$
y_j^n
=
\frac{
\left(\frac{1}{\Delta\tau}\right)
}{
\left(\frac{1}{\Delta\tau}\right)+\frac{rS_j}{h}+\sigma^2\!\left(V_{SS}^n\right)\!(j)\frac{S_j^2}{h^2}+r
},
$$
$$
z_j^n
=
\frac{
\frac{rS_j}{h}+\frac{1}{2}\sigma^2\!\left(V_{SS}^n\right)\!(j)\frac{S_j^2}{h^2}
}{
\left(\frac{1}{\Delta\tau}\right)+\frac{rS_j}{h}+\sigma^2\!\left(V_{SS}^n\right)\!(j)\frac{S_j^2}{h^2}+r
}.
$$

**Ambiguity flag (explicitness vs. Gamma discretization):** The above update is explicit provided $$\sigma^2\!\left(V_{SS}^n\right)\!(j)$$ is computed from already-known quantities at time level $$n$$. In the monotonicity proof, the pasted text includes discrete second-derivative formulas that involve $$V_j^{n+1}$$ inside $$V_{SS}^n(j)$$; if taken literally, that would make $$\sigma^2$$ depend on the new unknown and would generally destroy the closed-form explicit update. The coexistence of a closed-form update (18) and a proof that appears to perturb $$V_{SS}$$ using a stencil containing $$V^{n+1}$$ suggests a typographical index error in the pasted text or an implicit definition that the paper later treats as explicit; the provided excerpt does not fully resolve this mismatch.

### 3.3 Overall pipeline (ASCII sketch)

```
┌──────────────────────────────────────────────────────────────────────┐
│ Inputs                                                               │
│  - Model: r, σ0, transaction-cost parameter a, maturity T            │
│  - Contract: payoff L1(S), boundaries L2(τ), L3(τ), truncation Smax   │
│  - Grid: M (space), N (time)                                         │
└───────────────┬──────────────────────────────────────────────────────┘
                │
                ▼
┌──────────────────────────────────────────────────────────────────────┐
│ Discretize domain Ω = [0,Smax]×[0,T]                                 │
│  S_j = j h,  h = Smax/M ;  τ_n = n Δτ,  Δτ = T/N                      │
│  Initialize V_j^0 = L1(S_j) and boundary values V_0^n, V_M^n          │
└───────────────┬──────────────────────────────────────────────────────┘
                │  for n = 0..N-1
                ▼
┌──────────────────────────────────────────────────────────────────────┐
│ For each interior node j = 1..M-1                                    │
│  1) Approximate Gamma: V_SS^n(j)                                     │
│  2) Evaluate σ^2(V_SS^n)(j) via Barles–Soner volatility               │
│  3) Form coefficients x_j^n, y_j^n, z_j^n                             │
│  4) Update: V_j^{n+1} = x_j^n V_{j-1}^n + y_j^n V_j^n + z_j^n V_{j+1}^n │
└───────────────┬──────────────────────────────────────────────────────┘
                │
                ▼
┌──────────────────────────────────────────────────────────────────────┐
│ Output grid {V_j^n}; map back to calendar time by U(S,t)=V(S,T-t)     │
└──────────────────────────────────────────────────────────────────────┘
```

### 3.4 Step-by-step procedures (as required)

#### Procedure 1 — Evaluate the Barles–Soner-type modified volatility at a grid point

**Inputs:** $$\sigma_0$$, $$r$$, transaction-cost parameter $$a$$, current time $$\tau_n$$, current node $$S_j$$, a discrete approximation $$v$$ of $$V_{SS}(S_j,\tau_n)$$.  
**Outputs:** $$\sigma^2_j{}^{\,n} \approx \sigma^2(V_{SS})(S_j,\tau_n)$$.

1. Compute the argument $$w$$ passed to $$\Psi$$ using the paper’s volatility model (source text is ambiguous about the exact formula; the scheme analysis only requires the existence of a positive scaling $$K$$ such that $$w=K v$$).
2. Compute $$\Psi(w)$$:
   - Either by solving the IVP for $$\Psi$$ numerically, or
   - By numerically inverting the paper’s implicit characterization (piecewise for $$w>0$$ and $$w<0$$).
3. Return
   $$
   \sigma^2_j{}^{\,n} = \sigma_0^2\left(1+\Psi(w)\right).
   $$

**Implementation gap flagged:** The provided source lists the IVP and an implicit formula for $$\Psi$$, but does not provide a concrete numerical subroutine (e.g., Newton tolerance, iteration cap) for evaluating $$\Psi(w)$$ in code.

#### Procedure 2 — Explicit NSFD time stepping for the nonlinear Black–Scholes PDE

**Inputs:** Grid sizes $$M,N$$; $$S_{\max}$$, $$T$$; contract data $$L_1,L_2,L_3$$; parameters $$r,\sigma_0,a$$; a routine to compute $$\sigma^2\!\left(V_{SS}^n\right)\!(j)$$.  
**Outputs:** Numerical approximation $$V_j^n$$ for all $$j,n$$.

1. Set $$h=S_{\max}/M$$ and $$\Delta\tau=T/N$$.
2. Initialize $$V_j^0=L_1(S_j)$$ for $$j=0,\dots,M$$.
3. For each time index $$n=0,\dots,N-1$$:
   1. Impose boundaries: $$V_0^n=L_2(\tau_n)$$ and $$V_M^n=L_3(\tau_n)$$.
   2. For each interior index $$j=1,\dots,M-1$$:
      1. Compute a discrete Gamma proxy $$V_{SS}^n(j)$$ (not fully specified unambiguously in the pasted text).
      2. Evaluate $$\sigma^2\!\left(V_{SS}^n\right)\!(j)$$ by Procedure 1.
      3. Compute coefficients $$x_j^n,y_j^n,z_j^n$$ as in (23)–(25).
      4. Update
         $$
         V_j^{n+1} = x_j^n V_{j-1}^n + y_j^n V_j^n + z_j^n V_{j+1}^n.
         $$
   3. Impose boundaries for the new level if needed: $$V_0^{n+1}=L_2(\tau_{n+1})$$ and $$V_M^{n+1}=L_3(\tau_{n+1})$$.
4. (Optional) Convert back to calendar time: $$U(S_j,t_m)=V(S_j,T-t_m)$$.

## 4. Theoretical Results

The paper’s Section 4 proves the Barles–Souganidis-style triad—monotonicity, stability, consistency—then invokes the standard convergence result to the viscosity solution.

### 4.1 Lemma 1 (Monotonicity)

**Claim (restated in full):** The proposed NSFD scheme is monotone in the sense required for convergence to viscosity solutions of degenerate parabolic PDEs.

**Conditions explicitly used in the proof:** For the discrete operator $$\phi_j^{n+1}$$ defined by the scheme, it suffices to verify the standard monotonicity inequalities for perturbations by any $$\xi>0$$ of the neighboring and current unknowns (the paper writes these as inequalities analogous to its (28) and (29)). The proof also uses sign information on coefficients such as $$p_j \le 0$$ and $$q_j>0$$ arising from
$$
p_j = -\frac{rS_j}{h_j},\qquad
q_j = \frac{1}{\Delta\tau_n}+\frac{rS_j}{h_j}+r.
$$

**Proof sketch (2–4 sentences):** The argument rewrites the scheme as a residual map $$\phi_j^{n+1}(\cdot)=0$$ and checks that increasing neighbor values decreases the residual while increasing the current value increases it, matching the monotonicity inequalities. The nonlinear diffusion contribution is reduced to monotonicity of the scalar map $$v \mapsto (1+\Psi(Kv))v$$ for $$K>0$$, established by differentiating this expression and using the IVP for $$\Psi$$ to show the derivative is nonnegative. Combining monotonicity of the linear terms and this increasing property of the nonlinear term yields the scheme’s monotonicity.

### 4.2 Proposition 1 (Positivity preservation)

**Claim (restated in full):** If $$V_{j-1}^n$$, $$V_j^n$$, and $$V_{j+1}^n$$ are nonnegative real numbers, then the explicit update (18) produces a nonnegative value $$V_j^{n+1}$$.

**Conditions stated:** Nonnegativity of parameters and coefficients, including
$$
S,t,T,r,h,\Delta\tau,\sigma_0 \ge 0,\qquad \sigma^2 \ge 0.
$$

**Proof sketch:** The update is written as a weighted sum
$$
V_j^{n+1} = x_j^n V_{j-1}^n + y_j^n V_j^n + z_j^n V_{j+1}^n,
$$
and the stated parameter conditions imply $$x_j^n\ge 0$$, $$y_j^n\ge 0$$, $$z_j^n\ge 0$$. Nonnegativity of the previous-step values then yields $$V_j^{n+1}\ge 0$$ by closure of $$\mathbb{R}_{\ge 0}$$ under nonnegative linear combinations.

### 4.3 Lemma 2 (Unconditional stability)

**Claim (restated in full):** The new NSFD method is unconditionally stable (the paper uses the $$\ell_\infty$$ norm).

**Stated bound form:** The paper bounds the computed solution by the maximum norm of initial and boundary data, in the form
$$
\left\|V^{n+1}\right\|_\infty \le \max\Bigl(\left\|L_1\right\|_\infty,\left\|L_2\right\|_\infty,\left\|L_3\right\|_\infty\Bigr).
$$

**Proof sketch:** The proof uses the convex-combination form $$V_j^{n+1} = x_j^n V_{j-1}^n + y_j^n V_j^n + z_j^n V_{j+1}^n$$ together with $$x_j^n,y_j^n,z_j^n\ge 0$$. It then shows (for interior indices) that $$x_j^n+y_j^n+z_j^n \le 1$$ when $$r\ge 0$$, implying a maximum-norm nonexpansiveness across time steps. Boundary indices are controlled directly by the imposed boundary conditions, and combining interior and boundary cases yields the stated uniform bound for all time levels.

### 4.4 Lemma 3 (Consistency and truncation error)

**Claim (restated in full):** The method (18) is consistent, and its local truncation error satisfies
$$
\text{LTE} = \mathcal{O}(\Delta\tau, h^2).
$$

**Proof sketch:** The paper defines a local truncation error by plugging the exact solution $$V(S_j,\tau_n)$$ into the discrete operator and then expands the temporal and spatial shifts via Taylor series in $$\Delta\tau$$ and $$h$$. The leading-order terms reproduce the continuous nonlinear Black–Scholes PDE, so they cancel. The remaining principal terms scale linearly in $$\Delta\tau$$ and quadratically in $$h$$, giving $$\mathcal{O}(\Delta\tau, h^2)$$.

### 4.5 Theorem 1 (Convergence to the viscosity solution)

**Claim (restated in full):** The discrete solution produced by method (18) converges to the viscosity solution of the nonlinear Black–Scholes equation as the mesh is refined, provided
$$
h = \max_{0\le j \le M-1} h_j \to 0,\qquad
\Delta\tau = \max_{0\le n \le N-1} \Delta\tau_n \to 0.
$$

**Proof sketch:** The theorem is an application of the standard viscosity-solution convergence result for degenerate parabolic PDEs: any numerical scheme that is monotone, stable, and consistent converges to the unique viscosity solution (under the usual comparison principle assumptions). The paper cites this result and concludes convergence by combining Lemma 1 (monotonicity), Lemma 2 (stability), and Lemma 3 (consistency).

**Complexity bounds:** No formal complexity theorem is stated. The explicit pointwise update implies a direct implementation cost of order $$\mathcal{O}(MN)$$ operations for a fixed-cost evaluation of $$\Psi$$ per grid point (this is an inference from the algorithmic form, not a paper-stated bound).

## 5. Experimental Evaluation

The numerical section is demonstrative rather than benchmark-driven: it plots option value surfaces and cross-sections for several payoffs and several transaction-cost parameter values, emphasizing positivity and absence of spurious oscillations.

### 5.1 Experimental configuration table

| Item | What the paper reports |
|---|---|
| Implementation | MATLAB (no code listing provided) |
| “Datasets” / test problems | Four option contracts/payoffs: European call, European put, butterfly spread, cash-or-nothing (digital) |
| PDE model | Nonlinear Black–Scholes in $$\tau = T-t$$ with Barles–Soner-type modified volatility depending on $$V_{SS}$$ |
| Spatial domain | $$S \in (0,S_{\max})$$ with truncation $$S_{\max}=80$$ in all shown examples |
| Temporal domain | $$\tau \in (0,T]$$ with $$T=1$$ in all shown examples |
| Baselines | None reported (no cross-method error table; no convergence-rate plots) |
| Metrics | Qualitative: positivity/non-oscillation; monotone dependence of price on transaction-cost parameter; no explicit error norms reported |
| Ablation axis | Transaction cost parameter $$a$$ varied (plots show $$a=0$$, $$a=0.02$$, $$a=0.05$$) |
| Key hyperparameters | Uniform grid steps $$h$$ and $$\Delta t$$ (equivalently $$\Delta\tau$$) per example |

### 5.2 Example-by-example setup (as stated)

- **Example 1 (European call):** Parameters include $$r=0.1$$, $$a=0.05$$, $$\sigma_0=0.2$$, $$K=40$$, $$T=1$$, $$S_{\max}=80$$, with uniform mesh $$h=4$$ and $$\Delta t=0.1$$ (time step in the plotted variable). Figures show a positive, stable surface and cross-sections at $$t=0$$ for multiple $$a$$ values.

- **Example 2 (reported as put option):** Parameters listed include $$K_1=30$$, $$K_2=40$$, $$K_3=50$$, $$r=0.1$$, $$\sigma_0=0.2$$, $$T=1$$, $$S_{\max}=80$$ with mesh $$h=2$$ and $$\Delta t=0.05$$.

  **Ambiguity flag:** The presence of $$K_1,K_2,K_3$$ is typical of a butterfly payoff, not a plain put payoff; the pasted source does not reconcile this mismatch.

- **Example 3 (butterfly spread):** Parameters listed include $$r=0.1$$, $$a=0.05$$, $$\sigma_0=0.2$$, $$K=40$$, $$T=1$$, $$S_{\max}=80$$, with mesh $$h=4$$ and $$\Delta t=0.1$$.

  **Ambiguity flag:** A butterfly payoff requires $$K_1,K_2,K_3$$; the pasted source lists a single $$K$$. The figures are labeled as butterfly spread, but the parameter listing is internally inconsistent in the provided text.

- **Example 4 (cash-or-nothing):** Parameters include $$r=0.1$$, $$\sigma_0=0.2$$, $$K=40$$, $$T=1$$, $$B=1$$, $$S_{\max}=80$$, with mesh $$h=2$$ and $$\Delta t=0.05$$.

### 5.3 Reported outcomes (qualitative and comparative)

- The plotted solutions in all four examples remain nonnegative over the grid and appear free of visible spurious oscillations for the reported step sizes.
- Cross-section comparisons at $$t=0$$ (equivalently $$\tau=T$$) show that increasing the transaction-cost parameter $$a$$ increases the computed option price (demonstrated for $$a=0$$, $$a=0.02$$, $$a=0.05$$).
- No numerical error tables against an analytical solution, no reference-solution comparisons, and no measured convergence rates are provided in the excerpt; the evaluation is primarily visual and property-based.

## 6. ASCII Architecture / Workflow Diagram(s)

### 6.1 End-to-end workflow (contract to price grid)

```
┌──────────────────────────────┐
│ Contract specification        │
│  - payoff L1(S)               │
│  - boundaries L2(τ), L3(τ)     │
│  - strikes K (or K1,K2,K3), B  │
└───────────────┬──────────────┘
                ▼
┌──────────────────────────────┐
│ Model parameters              │
│  r, σ0, transaction-cost a, T │
│  Barles–Soner Ψ(·) definition │
└───────────────┬──────────────┘
                ▼
┌──────────────────────────────┐
│ Mesh                          │
│  S_j = jh, j=0..M              │
│  τ_n = nΔτ, n=0..N              │
│  V_j^0 = L1(S_j)               │
└───────────────┬──────────────┘
                ▼   time loop n=0..N-1
┌──────────────────────────────────────────────────────────┐
│ Boundary injection                                        │
│  V_0^n = L2(τ_n),  V_M^n = L3(τ_n)                          │
└───────────────┬──────────────────────────────────────────┘
                ▼
┌──────────────────────────────────────────────────────────┐
│ Interior update (for each j=1..M-1)                        │
│  1) compute v = V_SS^n(j)                                  │
│  2) compute σ^2 = σ0^2(1 + Ψ(K v))                          │
│  3) compute x_j^n, y_j^n, z_j^n                             │
│  4) V_j^{n+1} = x_j^n V_{j-1}^n + y_j^n V_j^n + z_j^n V_{j+1}^n │
└───────────────┬──────────────────────────────────────────┘
                ▼
┌──────────────────────────────┐
│ Output                        │
│  V(S_j,τ_n) grid               │
│  U(S,t)=V(S,T-t) if desired    │
└──────────────────────────────┘
```

### 6.2 Local stencil view (single update)

```
Time level n:     V_{j-1}^n      V_j^n      V_{j+1}^n
                     │            │            │
                     └──────┬─────┴─────┬──────┘
                            │  x,y,z    │
                            ▼           │
Time level n+1:            V_j^{n+1}  (boundary nodes fixed separately)

Update:
  V_j^{n+1} = x_j^n V_{j-1}^n + y_j^n V_j^n + z_j^n V_{j+1}^n
  with x_j^n,y_j^n,z_j^n ≥ 0 and (typically) x_j^n+y_j^n+z_j^n ≤ 1
```

## 7. Follow-Up Works & Extensions

### 7.1 Direct extensions in the NSFD line (citing/continuing the approach)

**Generalized-coefficient Black–Scholes via NSFD.** The paper *A Nonstandard Finite Difference Method for a Generalized Black–Scholes Equation* proposes an NSFD scheme for a generalized Black–Scholes PDE with time- and space-dependent coefficients, focusing again on qualitative properties such as positivity and stability and reporting second-order accuracy in space; it explicitly lists the 2021 “qualitatively stable NSFD” paper in its references, positioning it as prior work in the same numerical-design program. [Khalsaraei et al., Symmetry 2022]. ([mdpi.com](https://www.mdpi.com/2073-8994/14/1/141?utm_source=openai))

**“Improved” positivity-preserving NSFD comparisons for Black–Scholes.** The paper *A Positivity-Preserving Improved Nonstandard Finite Difference Method to Solve the Black–Scholes Equation* evaluates multiple standard schemes (including variants such as $$\theta$$-methods and explicit/implicit combinations) and frames an improved NSFD-style discretization as a way to eliminate oscillations while preserving nonnegativity; it also includes the 2021 transaction-cost nonlinear NSFD paper in its reference list, indicating direct bibliographic continuity. [Khalsaraei et al., Mathematics 2022]. ([mdpi.com](https://www.mdpi.com/2227-7390/10/11/1846?utm_source=openai))

**Laplace-transform plus NSFD “qualitatively stable schemes.”** *Qualitatively Stable Schemes for the Black–Scholes Equation* constructs a scheme by combining a Laplace-transform approach with NSFD design and analyzes positivity/stability/consistency under low-volatility regimes; it explicitly states that the methodology is intended to avoid spurious oscillations typical of standard finite differences in certain parameter regimes and notes that the approach can be extended toward nonlinear Black–Scholes-type problems. [Khalsaraei et al., Fractal and Fractional 2023]. ([mdpi.com](https://www.mdpi.com/2504-3110/7/2/154?utm_source=openai))

### 7.2 Alternative discretizations for nonlinear transaction-cost Black–Scholes PDEs (same target problem class)

**High-order kernel/RBF-enhanced finite difference for nonlinear Black–Scholes.** *Hermite Finite Difference Through Kernel Approximations to Efficiently Solve Nonlinear Black-Scholes Model* develops a high-order compact discretization strategy for nonlinear Black–Scholes models under transaction-cost/feedback-type nonlinearities, aiming at higher spatial accuracy while retaining computational efficiency; its reference list includes the 2021 NSFD paper, making it a direct citing follow-up from a different numerical-approximation family. [Wang et al., Mathematics 2025]. ([mdpi.com](https://www.mdpi.com/2227-7390/13/17/2727?utm_source=openai))

### 7.3 Closely related monotone/stable/viscosity-convergent schemes for the Barles–Soner transaction-cost model (key comparators)

**Fully implicit + upwind finite differences with viscosity convergence.** *An upwind finite difference method for a nonlinear Black–Scholes equation governing European option valuation under transaction costs* (Lesmana and Wang) develops a monotone spatial discretization and a fully implicit time step, proves an M-matrix property and unconditional convergence to the viscosity solution, and uses Newton iterations to solve the nonlinear algebraic system; the 2021 NSFD paper’s main point of differentiation is avoiding such nonlinear global solves by providing a closed-form explicit update while preserving the same stability/monotonicity/consistency triad. [Lesmana & Wang, Appl. Math. Comput. 2013]. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313000180?utm_source=openai))

**Consistent explicit difference schemes with step-size stability conditions.** *Consistent stable difference schemes for nonlinear Black–Scholes equations modelling option pricing with transaction costs* analyzes explicit discretizations for the nonlinear transaction-cost PDE and derives sufficient stability conditions linking the time/space steps to the transaction-cost parameter; this contrasts with the 2021 NSFD scheme’s “unconditional” qualitative-property claims achieved via nonstandard discretization choices rather than step-size restrictions. [Company et al., ESAIM: M2AN 2009]. ([esaim-m2an.org](https://www.esaim-m2an.org/articles/m2an/abs/2009/06/m2an0861/m2an0861.html?utm_source=openai))

**Earlier finite-difference analysis for the Barles–Soner nonlinear PDE.** *A numerical method for European Option Pricing with transaction costs nonlinear equation* presents a finite difference approach for the Barles–Soner nonlinear model and studies stability/consistency with numerical illustrations; the 2021 paper positions itself as improving computational cost and avoiding oscillations while keeping monotonicity and convergence-to-viscosity-solution guarantees. [Company et al., Math. Comput. Model. 2009]. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0895717709001538?utm_source=openai))

## 8. Industrial & Real-World Applications

No verified production deployment (bank/pricing system) explicitly attributed to the specific 2021 NSFD scheme was identified in the conducted searches.

**Open-source quantitative finance libraries that provide PDE/finite-difference option-pricing infrastructure (implementation substrate).** The following projects are widely used as practical tooling for option valuation and include finite-difference or PDE-related components where a nonlinear, Gamma-dependent-volatility scheme of the type studied in the source paper could be integrated, but they do not (based on the available, high-level public pages) claim an out-of-the-box implementation of the Barles–Soner transaction-cost PDE with the 2021 NSFD stencil. [GitHub: quantlib]. ([github.com](https://github.com/quantlib?utm_source=openai))

**Java ecosystem reimplementation (industrial-style API surface).** JQuantLib is a Java rewrite inspired by QuantLib, intended as a production-friendly quantitative finance framework; it provides option valuation tooling and is positioned for integration into Java applications, making it a plausible environment for incorporating specialized finite-difference schemes when organizations standardize on JVM stacks. [GitHub: frgomes/jquantlib]. ([github.com](https://github.com/frgomes/jquantlib?utm_source=openai))

**Web/JavaScript packaging for quantitative finance primitives.** quantlib.js provides a JavaScript distribution and notebook-style usage for QuantLib-like functionality; it supports web-based experimentation and could host prototypes of finite-difference schemes for educational or lightweight deployment contexts, although no verified Barles–Soner NSFD implementation is claimed on the referenced page. [GitHub: quantlibjs/ql]. ([github.com](https://github.com/quantlibjs/ql?utm_source=openai))

## 9. Consolidated Reference List

[1] Mohammad Mehdizadeh Khalsaraei, Mohammad Mehdi Rashidi, Ali Shokri, Higinio Ramos, Pari Khakzad. “A Nonstandard Finite Difference Method for a Generalized Black–Scholes Equation.” *Symmetry*, 2022. DOI: `10.3390/sym14010141`. ([mdpi.com](https://www.mdpi.com/2073-8994/14/1/141?utm_source=openai))

[2] Mohammad Mehdizadeh Khalsaraei, Ali Shokri, Higinio Ramos, Zahra Mohammadnia, Pari Khakzad. “A Positivity-Preserving Improved Nonstandard Finite Difference Method to Solve the Black-Scholes Equation.” *Mathematics*, 2022. DOI: `10.3390/math10111846`. ([mdpi.com](https://www.mdpi.com/2227-7390/10/11/1846?utm_source=openai))

[3] Mohammad Mehdizadeh Khalsaraei, Ali Shokri, Yuanheng Wang, Sohrab Bazm, Giti Navidifar, Pari Khakzad. “Qualitatively Stable Schemes for the Black–Scholes Equation.” *Fractal and Fractional*, 2023. DOI: `10.3390/fractalfract7020154`. ([mdpi.com](https://www.mdpi.com/2504-3110/7/2/154?utm_source=openai))

[4] Shuai Wang, Jiameihui Zhu, Tao Liu. “Hermite Finite Difference Through Kernel Approximations to Efficiently Solve Nonlinear Black-Scholes Model.” *Mathematics*, 2025. DOI: `10.3390/math13172727`. ([mdpi.com](https://www.mdpi.com/2227-7390/13/17/2727?utm_source=openai))

[5] D. C. Lesmana, S. Wang. “An upwind finite difference method for a nonlinear Black–Scholes equation governing European option valuation under transaction costs.” *Applied Mathematics and Computation*, 2013. DOI: `10.1016/j.amc.2012.12.077`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0096300313000180?utm_source=openai))

[6] Rafael Company, Lucas Jódar, José-Ramón Pintos. “Consistent stable difference schemes for nonlinear Black-Scholes equations modelling option pricing with transaction costs.” *ESAIM: Mathematical Modelling and Numerical Analysis*, 2009. DOI: `10.1051/m2an/2009014`. ([esaim-m2an.org](https://www.esaim-m2an.org/articles/m2an/abs/2009/06/m2an0861/m2an0861.html?utm_source=openai))

[7] Rafael Company, Lucas Jódar, J.-R. Pintos. “A numerical method for European Option Pricing with transaction costs nonlinear equation.” *Mathematical and Computer Modelling*, 2009. DOI: `10.1016/j.mcm.2009.05.019`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0895717709001538?utm_source=openai))

[8] QuantLib (GitHub organization). “QuantLib: A free/open-source library for quantitative finance.” Repository hub: `https://github.com/quantlib`. ([github.com](https://github.com/quantlib?utm_source=openai))

[9] Richard Gomes (maintainer) and contributors. “JQuantLib: A library for Quantitative Finance written in 100% Java.” GitHub repository: `https://github.com/frgomes/jquantlib`. ([github.com](https://github.com/frgomes/jquantlib?utm_source=openai))

[10] quantlibjs contributors. “quantlib.js (ql): JavaScript packaging and notebook tooling for quantlib.js.” GitHub repository: `https://github.com/quantlibjs/ql`. ([github.com](https://github.com/quantlibjs/ql?utm_source=openai))

---
Learn more:
1. [A Nonstandard Finite Difference Method for a Generalized Black–Scholes Equation](https://www.mdpi.com/2073-8994/14/1/141?utm_source=openai)
2. [A Positivity-Preserving Improved Nonstandard Finite Difference Method to Solve the Black-Scholes Equation](https://www.mdpi.com/2227-7390/10/11/1846?utm_source=openai)
3. [Qualitatively Stable Schemes for the Black–Scholes Equation | MDPI](https://www.mdpi.com/2504-3110/7/2/154?utm_source=openai)
4. [Hermite Finite Difference Through Kernel Approximations to Efficiently Solve Nonlinear Black-Scholes Model | MDPI](https://www.mdpi.com/2227-7390/13/17/2727?utm_source=openai)
5. [An upwind finite difference method for a nonlinear Black–Scholes equation governing European option valuation under transaction costs - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0096300313000180?utm_source=openai)
6. [Consistent stable difference schemes for nonlinear Black-Scholes equations modelling option pricing with transaction costs | ESAIM: Mathematical Modelling and Numerical Analysis (ESAIM: M2AN)](https://www.esaim-m2an.org/articles/m2an/abs/2009/06/m2an0861/m2an0861.html?utm_source=openai)
7. [A numerical method for European Option Pricing with transaction costs nonlinear equation - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0895717709001538?utm_source=openai)
8. [QuantLib · GitHub](https://github.com/quantlib?utm_source=openai)
9. [GitHub - frgomes/jquantlib: JQuantLib is a library for Quantitative Finance written in 100% Java](https://github.com/frgomes/jquantlib?utm_source=openai)
10. [GitHub - quantlibjs/ql](https://github.com/quantlibjs/ql?utm_source=openai)

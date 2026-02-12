## 1. Paper Identity

- **Title:** Radial basis functions with application to finance: American put option under jump diffusion  
- **Authors:** Ahmad Golbabai; Davood Ahmadian; Mariyan Milev ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0895717711006066))  
- **Affiliations (as stated in the paper):**
  - School of Mathematics, Iran University of Science and Technology, Narmak, 16844 Tehran, Iran
  - Department of Informatics and Statistics, Faculty of Economics, University of Food Technologies – Plovdiv, bul. Maritza 26, 4002 Plovdiv, Bulgaria
- **Venue:** *Mathematical and Computer Modelling* (Elsevier), Volume 55, Issues 3–4, February 2012, pp. 1354–1362 ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0895717711006066))  
- **DOI:** 10.1016/j.mcm.2011.10.014 ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0895717711006066))  

## 2. Problem Statement & Formulation

The paper studies numerical valuation of an American put option when the underlying stock price follows a diffusion with **finite-activity jumps** (Poisson arrivals with a multiplicative jump factor). The valuation problem is posed as a **free-boundary partial integro-differential equation (PIDE)**: the continuation region satisfies a backward PIDE with a nonlocal integral term; the early-exercise feature is enforced through a moving optimal exercise boundary $$b(t)$$ (equivalently $$b(\tau)$$ after time reversal).

### 2.1 Underlying jump-diffusion dynamics (risk-neutral form used in the paper)

The underlying is modeled (paper’s notation) by the stochastic differential equation
$$
\frac{dS}{S} = r\,dt + \sigma\,dW_t + (\eta-1)\,dq,
$$
where $$r$$ is the interest rate, $$\sigma$$ is the diffusion volatility, $$W_t$$ is a Wiener process, $$q$$ is a Poisson process with intensity $$\lambda$$, and $$\eta$$ is a (positive) jump-amplitude multiplier so that a jump maps $$S \mapsto S\eta$$. The Poisson increment satisfies
$$
dq =
\begin{cases}
0, & \text{with probability } 1-\lambda\,dt,\\
1, & \text{with probability } \lambda\,dt,
\end{cases}
$$
and the paper defines the expected relative jump size
$$
k = \mathbb{E}(\eta-1).
$$

### 2.2 American put option free-boundary PIDE in the asset variable

Let $$V(S,t)$$ be the option value, strike $$K$$, maturity $$T$$. The paper writes the jump-diffusion PIDE as
$$
-\frac{\partial V}{\partial t}
=
\frac{1}{2}\sigma^2 S^2\frac{\partial^2 V}{\partial S^2}
+ (r-\lambda k)S\frac{\partial V}{\partial S}
-(r+\lambda)V
+ \lambda \int_0^\infty V(S\eta,t)\,\rho(\eta)\,d\eta,
$$
with terminal and boundary conditions
$$
V(S,T) = \max(K-S,0),
$$
$$
V(S,t)\to 0 \text{ as } S\to\infty,
\qquad
V(0,t)=K,
$$
and a free-boundary (value-matching) condition at the optimal exercise boundary $$b(t)$$:
$$
\lim_{S\to b(t)} V(S,t) = K-b(t).
$$
The text later uses an additional derivative boundary condition after front-fixing (smooth pasting in transformed variables), implying the standard smooth-pasting condition is part of the enforced boundary set, even though it is not written explicitly in the original-$$S$$ system.

### 2.3 Jump-size distribution

The jump multiplier density is specified as lognormal:
$$
\rho(\eta)\ge 0,
\qquad
\int_0^\infty \rho(\eta)\,d\eta = 1,
$$
$$
\rho(\eta) = \frac{1}{\sqrt{2\pi}\sigma_J \eta}
\exp\!\left(
-\frac{(\ln\eta-\mu)^2}{2\sigma_J^2}
\right),
$$
with log-mean parameter $$\mu$$ and log-volatility parameter $$\sigma_J$$.

### 2.4 Rescaling and front-fixing transformation (fixed boundary)

The paper rescales variables by the strike:
$$
S = K\hat{S},
\qquad
P_A(S,\tau)=\hat{P}_A(\hat{S},\tau),
\qquad
b(t)=K\hat{b}(t),
$$
and switches to time-to-maturity $$\tau=T-t$$. It then applies a **front-fixing transformation** (paper’s description) of the form
$$
y = \ln\!\left(\frac{\hat{S}}{\hat{b}(\tau)}\right),
$$
so that the free boundary $$\hat{b}(\tau)$$ maps to the fixed boundary $$y=0$$. The transformed PIDE is given (paper’s Eq. (10)) as
$$
\frac{\partial V}{\partial \tau}
-\frac{1}{b(\tau)}\frac{db(\tau)}{d\tau}\frac{\partial V}{\partial y}
=
\frac{1}{2}\sigma^2\!\left(\frac{\partial^2 V}{\partial y^2}-\frac{\partial V}{\partial y}\right)
+ (r-\lambda k)\frac{\partial V}{\partial y}
-(r+\lambda)V
+
\frac{\lambda}{b(\tau)e^{y}}
\int_0^\infty V(u)\,
\rho\!\left(\frac{u}{b(\tau)e^{y}}\right)\,du.
$$

**Ambiguity flag (notation):** in the provided paper text, the same symbol $$V$$ is used both for the transformed solution $$V(y,\tau)$$ and (inside the jump integral) as a function of the integration variable $$u$$. This is a notational compression; the implementation described later constructs a matrix acting on RBF coefficients, effectively treating the integrand as “the option value evaluated at quadrature points,” but the mapping from $$u$$ to the transformed coordinate is not written explicitly in a single consistent notation line.

### 2.5 Boundary conditions in fixed-boundary coordinates

The paper states (Eq. (12)) conditions consistent with value matching and smooth pasting at $$y=0$$ plus far-field decay:
- At expiry, for $$0<y<\infty$$, the transformed payoff becomes identically zero:
  $$
  V(y,0) = \max(1-b(0)e^y,0)=0.
  $$
  **Ambiguity flag (OCR/typo):** the text contains “as $$b(0)=10$$,” which is inconsistent with the preceding equation and with the scaling; context strongly suggests this is “$$b(0)=1.0$$,” yielding $$V(y,0)=0$$ for all $$y>0$$.
- Boundary at $$y=0$$:
  $$
  V(0,\tau)=1-b(\tau),
  \qquad
  \frac{\partial V(0,\tau)}{\partial y}=-b(\tau).
  $$
- Far field:
  $$
  \lim_{y\to\infty} V(y,\tau)=0.
  $$

They then algebraically combine the first two $$y=0$$ conditions into the fixed-boundary constraint (Eq. (13)):
$$
V(0,\tau)-\frac{\partial V(0,\tau)}{\partial y}=1,
\qquad
\lim_{y\to\infty} V(y,\tau)=0,
\qquad
V(y,0)=0 \text{ for } y>0.
$$

## 3. Core Methodology

The method combines: (i) front-fixing to eliminate the moving boundary in the spatial domain; (ii) Gauss–Laguerre quadrature (presented as Laguerre-polynomial interpolation) to discretize the semi-infinite jump integral; (iii) a global **radial basis function (RBF)** collocation method (multiquadric kernel) for spatial approximation; (iv) Crank–Nicolson time stepping; (v) a **Predictor–Corrector** iteration to address the nonlinearity introduced by $$b(\tau)$$ and $$db/d\tau$$; (vi) a projection step to enforce the early-exercise constraint.

### 3.1 Overall pipeline (paper-level view)

```text
┌──────────────────────────────────────────────┐
│ Model parameters: K, r, σ, λ, k, μ, σ_J, T    │
└───────────────────────────────┬──────────────┘
                                ▼
┌──────────────────────────────────────────────┐
│ Free-boundary American-put PIDE in (S,t)      │
│ + jump integral term + unknown exercise b(t)  │
└───────────────────────────────┬──────────────┘
                                ▼
┌──────────────────────────────────────────────┐
│ Scale by K and front-fix: y = ln(Ŝ / b̂(τ))   │
│ Fixed boundary at y=0; b(τ) encoded in BCs    │
└───────────────┬──────────────────────────────┘
                │
     ┌──────────┴──────────┐
     ▼                     ▼
┌──────────────┐   ┌───────────────────────────┐
│ Jump integral │   │ Spatial approximation via │
│ on [0,∞):     │   │ MQ-RBF collocation        │
│ Gauss-Laguerre│   │ (dense global system)     │
└──────┬───────┘   └──────────────┬────────────┘
       └──────────────┬───────────┘
                      ▼
┌──────────────────────────────────────────────┐
│ Time stepping: Crank–Nicolson +               │
│ Predictor–Corrector for nonlinear coupling    │
│ b(τ) ↔ V(0,τ) and convective boundary term    │
└───────────────────────────────┬──────────────┘
                                ▼
┌──────────────────────────────────────────────┐
│ Projection: V := max(V, intrinsic payoff)     │
│ Update b(τ) from boundary relation            │
└──────────────────────────────────────────────┘
```

### 3.2 Procedure: front-fixing transformation and fixed-boundary problem

**Procedure 1 (Front-fixing and scaling).**

1. **Inputs:** strike $$K$$; maturity $$T$$; model parameters $$r,\sigma,\lambda,k,\mu,\sigma_J$$; unknown exercise boundary $$b(t)$$.
2. Define time-to-maturity $$\tau=T-t$$.
3. Rescale asset and boundary by $$K$$:
   $$
   \hat{S}=\frac{S}{K},
   \qquad
   b(t)=K\hat{b}(t).
   $$
4. Apply the front-fixing variable transform
   $$
   y=\ln\!\left(\frac{\hat{S}}{\hat{b}(\tau)}\right),
   $$
   mapping the free boundary $$\hat{S}=\hat{b}(\tau)$$ to $$y=0$$.
5. Obtain the transformed fixed-boundary PIDE (paper’s Eq. (10)), with nonlinearity entering through the convective term
   $$
   -\frac{1}{b(\tau)}\frac{db(\tau)}{d\tau}\frac{\partial V}{\partial y}.
   $$
6. Enforce boundary conditions at $$y=0$$ via
   $$
   V(0,\tau)-V_y(0,\tau)=1,
   $$
   and far-field truncation for the semi-infinite spatial domain (paper discusses truncating the original $$S$$-domain at $$S_{\max}$$, often taking $$S_{\max}=3K$$ or $$S_{\max}=2K$$ in finite-difference practice).

**Outputs:** fixed-boundary PIDE in $$y\in[0,\infty)$$ (numerically truncated) coupled to scalar boundary function $$b(\tau)$$.

### 3.3 Procedure: approximation of the jump integral by Laguerre polynomials (Gauss–Laguerre quadrature)

The paper approximates the integral term over $$[0,\infty)$$ using Laguerre polynomials, implemented as Gauss–Laguerre quadrature with roots $$r_k$$ and weights $$w_k$$, and fixes the quadrature order to $$N_2=15$$.

**Procedure 2 (Jump integral quadrature operator assembly).**

1. **Inputs:** quadrature size $$N_2=15$$; Laguerre roots $$\{r_k\}_{k=1}^{N_2}$$ and weights $$\{w_k\}_{k=1}^{N_2}$$; current boundary scalar $$b(\tau)$$; spatial collocation points $$\{x_i\}$$ (paper uses $$x$$ and $$y$$ interchangeably here); RBF centers $$\{u_j\}$$; MQ kernel $$\varphi$$.
2. Start from the jump integral term in the transformed equation and represent the option value by an RBF expansion (see Procedure 3), so the integral becomes a weighted sum of integrals of basis functions.
3. Convert the integral to a Gauss–Laguerre form by inserting an exponential weight and compensating with an $$e^{r_k}$$ factor (paper’s Eq. (19) shows the characteristic “$$+r_k$$” term inside an exponential).
4. Define a matrix $$G(b(\tau))$$ whose entries are (paper’s Eq. (20), paraphrased structurally)
   $$
   [G(b(\tau))]_{i,j}
   =
   \sum_{k=1}^{N_2}
   w_k\;
   \frac{\varphi(r_k-u_j)}{r_k}\;
   \exp\!\bigl(-\text{kernel}(r_k,x_i,b(\tau),\mu,\sigma_J)+r_k\bigr),
   $$
   where $$\text{kernel}(\cdot)$$ corresponds to inserting the lognormal density (paper’s Eq. (8)) into the transformed integral term.
5. Use $$G(b(\tau))$$ to approximate the integral term at all collocation points as a matrix–vector product acting on the RBF coefficient vector $$\lambda(\tau)$$.

**Outputs:** a computable linear operator (matrix) approximating the jump integral term, parameterized by $$b(\tau)$$.

**Ambiguity flag (kernel expression):** in the provided text, the explicit exponent expression in the reduced-form Eq. (11) contains variable-name collisions (e.g., “$$x$$” vs. “$$y$$”), inconsistent parentheses, and terms such as “$$\ln(u-\mu)$$” that are dimensionally inconsistent with a lognormal kernel. The assembly logic in Eqs. (19)–(20) is clear (Gauss–Laguerre quadrature applied to a basis-function expansion), but the exact printed kernel should be treated as “as written in the source” rather than algebraically reliable without consulting the typeset PDF.

### 3.4 Procedure: MQ-RBF collocation for the PDE part

The paper uses a global collocation RBF approximation with a **multiquadric (MQ)** kernel:
$$
\varphi(u)=\sqrt{u^2+c^2},
$$
with **shape parameter** $$c$$ controlling accuracy/conditioning.

**Procedure 3 (RBF approximation and collocation system).**

1. **Inputs:** spatial collocation nodes $$\{x_i\}_{i=1}^N$$ on a truncated domain; partition of nodes into interior points $$1\le i\le N_\Omega$$ and boundary points $$N_\Omega+1\le i\le N$$; MQ kernel $$\varphi(\cdot)$$ with shape $$c$$.
2. Approximate the solution by an RBF expansion:
   $$
   V(u,\tau)\approx \tilde{V}(u,\tau)=\sum_{j=1}^{N}\lambda_j(\tau)\,\varphi(u-u_j).
   $$
   (The paper’s Eq. (14) uses an infinite sum in text; the implemented method is finite-dimensional with $$N$$ coefficients.)
3. Evaluate the PIDE operator at interior collocation points and enforce boundary operator conditions at boundary collocation points (paper’s Eqs. (16)–(18)).
4. Precompute interpolation and differentiation matrices built from
   $$
   \varphi_{ij}=\varphi(\lVert x_i-x_j\rVert),
   $$
   and its first and second derivatives (paper denotes these abstractly as $$L\varphi_{ij}$$ and $$L^2\varphi_{ij}$$; in the option PIDE context they correspond to spatial first/second derivatives needed for $$V_y$$ and $$V_{yy}$$, plus boundary-operator enforcement).

**Outputs:** a dense square system (global MQ-RBF collocation) coupling the coefficient vector $$\lambda(\tau)$$ to the PIDE residuals and boundary conditions.

### 3.5 Procedure: Crank–Nicolson time stepping and Predictor–Corrector nonlinear iteration

The transformed PIDE is nonlinear because $$b(\tau)$$ is unknown and enters both:
- the convective coefficient $$\frac{1}{b(\tau)}\frac{db(\tau)}{d\tau}$$, and
- the boundary condition relations tying $$b(\tau)$$ to $$V(0,\tau)$$.

The paper handles this via:
- a Crank–Nicolson (time-centered) discretization of the PDE terms, and
- a Predictor–Corrector iteration for the nonlinear coupling.

It explicitly states the boundary update:
$$
b(\tau_i)=1-V(0,\tau_i).
$$

**Procedure 4 (One time step $$\tau_i \mapsto \tau_{i+1}$$).**

1. **Inputs at time level $$\tau_i$$:** coefficient vector $$\lambda(\tau_i)$$; values and derivatives of $$V$$ implied by the RBF expansion; current boundary scalar $$b(\tau_i)=1-V(0,\tau_i)$$; time step size
   $$
   h_1=\frac{T}{N_1}.
   $$
2. Assemble the jump-integral matrices (paper denotes $$G_1(b(\tau_i))$$ and $$G_2(b(\tau_{i+1}))$$) using Procedure 2.
3. **Predictor step (paper’s Eq. (21)):** compute a provisional $$V(\cdot,\tau_{i+1})$$ using:
   - a forward-difference approximation of $$V_\tau$$,
   - a discrete approximation of the nonlinear convective coefficient using $$\frac{b(\tau_{i+1})-b(\tau_i)}{h_1 b(\tau_i)}$$ multiplying a known spatial derivative term at $$\tau_i$$,
   - Crank–Nicolson averaging of the diffusion/advection/reaction PDE terms between $$\tau_i$$ and $$\tau_{i+1}$$,
   - trapezoidal-style averaging of the jump-integral contribution (as written in Eq. (21)).
4. **Corrector iterations (paper’s Eq. (24)):** for $$\mu=1,\dots,M$$, refine the provisional solution via a corrected scheme that:
   - averages spatial derivatives at $$\tau_i$$ and the current iterate at $$\tau_{i+1}$$,
   - evaluates the jump-integral contribution using the updated iterate,
   - updates $$b^{\mu+1}(\tau_{i+1})$$ through its boundary relation with $$V^{\mu+1}(0,\tau_{i+1})$$.
   The paper fixes $$M$$ and uses $$V^{M+1}(\cdot,\tau_{i+1})\approx V(\cdot,\tau_{i+1})$$.
5. Express each iterate in terms of RBF coefficients:
   $$
   V(u_i,\tau_j)=\sum_{m=1}^{N}\lambda_m(\tau_j)\,\varphi(u_i-u_m),
   $$
   convert the predictor/corrector relations into a nonlinear algebraic system in $$\lambda(\tau_{i+1})$$ (paper’s Eqs. (25)–(26)).
6. Solve the resulting (nonlinear) collocation system, represented abstractly in the paper as
   $$
   B_1\lambda(t_{j+1})
   =
   A\!\bigl(B_1\lambda(t_{j+1}),B_2\lambda(t_{j+1}),B_3\lambda(t_{j+1}),
           C_1\lambda(t_j),C_2\lambda(t_j),C_3\lambda(t_j)\bigr),
   $$
   with matrices $$B_\ell,C_\ell$$ defined entrywise by (paper’s Eqs. (27)–(32)) using interior/boundary operator application to basis functions.

**Termination / convergence:** the corrector loop is run for a fixed, user-chosen number of iterations $$M$$; no residual-based stopping criterion is specified in the text.

### 3.6 Procedure: early-exercise (no-arbitrage) projection

After each time update, the paper enforces the American put constraint by projecting the computed continuation value onto the intrinsic payoff:
$$
V(x_i,\tau_j)=\max\!\left(V(x_i,\tau_j),\,1-b(\tau_j)\exp(x_i)\right).
$$

**Procedure 5 (Projection and boundary update).**

1. **Inputs:** computed $$V(x_i,\tau_j)$$ at collocation nodes; current $$b(\tau_j)$$.
2. Apply the pointwise projection:
   $$
   V \leftarrow \max(V,\text{intrinsic}).
   $$
3. Update boundary scalar using the boundary relation:
   $$
   b(\tau_j)=1-V(0,\tau_j).
   $$

**Outputs:** arbitrage-consistent American option value surface (numerically) and the corresponding boundary sequence $$\{b(\tau_j)\}$$.

## 4. Theoretical Results

No explicitly labeled theorem, lemma, proposition, or corollary appears in the provided paper text. No formal complexity bounds (Big-O or Big-Theta) are stated. The paper asserts that the proposed method is stable and that results agree with other numerical methods, but does not provide a stand-alone proof of stability or convergence for the full nonlinear front-fixed PIDE discretization.

## 5. Experimental Evaluation

### 5.1 Experimental setup (model parameters, numerical settings, baselines)

| Category | Item | Value / Description |
|---|---|---|
| Contract | Option type | American put |
| Contract | Strike | $$K=100$$ |
| Contract | Maturity | $$T=0.25$$ |
| Diffusion | Volatility | $$\sigma=0.15$$ |
| Rates | Interest rate | $$r=0.05$$ |
| Jumps | Jump intensity | $$\lambda=0.10$$ |
| Jumps | Lognormal jump mean | $$\mu=-0.90$$ |
| Jumps | Jump log-volatility | $$\sigma_J=0.45$$ |
| Integral quadrature | Laguerre / Gauss–Laguerre order | $$N_2=15$$ (fixed) |
| Space discretization | “Size of the space $$S$$-grid” | $$N\in\{15,20,25,30,35\}$$ (as reported) |
| Time discretization | Number of time steps | $$N_1\in\{40,50,60,70,80\}$$ (as reported) |
| Nonlinear solve | Corrector iterations | $$M$$ fixed (value not reported) |
| RBF kernel | Type | MQ multiquadric $$\varphi(u)=\sqrt{u^2+c^2}$$ |
| RBF hyperparameter | Shape parameter | $$c$$ selected empirically to minimize error |
| Baseline / reference | Reference prices | Taken from d’Halluin–Forsyth–Vetzal (paper’s ref. [9]) at $$S_0\in\{90,100,110\}$$ |
| Metric | Error | Absolute error against the reference prices |

### 5.2 Quantitative results (Table 2 reproduced)

Reference values used for error computation (from the paper’s discussion of [9]):  
- at $$S_0=90$$: 10.003822  
- at $$S_0=100$$: 3.241251  
- at $$S_0=110$$: 1.419803  

| Time steps $$N_1$$ | Space grid size $$N$$ | $$S_0=90$$ value | $$S_0=90$$ error | $$S_0=100$$ value | $$S_0=100$$ error | $$S_0=110$$ value | $$S_0=110$$ error |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 | 15 | 10.007723 | 3e-3 | 3.23832 | 2e-3 | 1.424355 | 4e-3 |
| 50 | 20 | 10.003251 | 5e-4 | 3.239984 | 1e-3 | 1.410948 | 8e-3 |
| 60 | 25 | 10.003547 | 3e-4 | 3.2406385 | 6e-4 | 1.4198438 | 4e-5 |
| 70 | 30 | 10.003968 | 1e-4 | 3.2402608 | 9e-4 | 1.4200659 | 2e-4 |
| 80 | 35 | 10.0038211 | 1e-6 | 3.2413965 | 1e-4 | 1.418803 | 3e-4 |

### 5.3 Reported ablations / sensitivity axes (Figures 1–3)

- **Shape parameter selection vs. mesh size (Fig. 1):** the paper plots “optimal” $$c$$ values versus the stock-mesh size while fixing the time mesh size (stated as $$N_1=40$$). The “optimal” $$c$$ is defined operationally: for each fixed $$S_0$$ and fixed spatial grid size, scan multiple candidate $$c$$ values and pick the one minimizing the absolute pricing error.
- **Time resolution sensitivity (Fig. 2):** error vs. time mesh size with a fixed stock mesh size (stated as $$N=15$$). The plot indicates decreasing error as time steps increase.
- **Space resolution sensitivity (Fig. 3):** error vs. stock mesh size with fixed time mesh size (stated as $$N_1=40$$). The plot indicates decreasing error as spatial resolution increases.

### 5.4 MQ vs. IMQ comparison (reported qualitative finding)

The paper reports that, under the same shape-parameter choices used for MQ in Fig. 1, the inverse multiquadric (IMQ) variant matches reference values less accurately: agreement to only the 3rd–4th decimal at $$S_0=90$$, and poor agreement for $$S_0=100$$ and $$S_0=110$$ unless $$c$$ is re-tuned. No IMQ result table is provided; the paper explicitly omits those inaccurate IMQ values.

## 6. ASCII Architecture / Workflow Diagram(s)

### 6.1 End-to-end workflow (front-fixing + quadrature + RBF + time stepping)

```text
┌──────────────────────────────────────────────────────────────────────────┐
│ Inputs: r, σ, λ, k, μ, σ_J, K, T; quadrature order N2=15; RBF shape c     │
└──────────────────────────────────────┬───────────────────────────────────┘
                                       ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Scale: Ŝ=S/K; b̂=b/K; τ=T−t                                                 │
│ Front-fix: y = ln(Ŝ / b̂(τ)); continuation domain becomes y ≥ 0            │
└──────────────────────────────────────┬───────────────────────────────────┘
                                       ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Discretize y-domain (truncate): choose collocation nodes {x_i}_{i=1}^N    │
└──────────────────────────────────────┬───────────────────────────────────┘
                                       ▼
┌───────────────────────────────┬──────────────────────────────────────────┐
│ Jump integral (nonlocal term)  │ Spatial derivatives (local PDE part)     │
│ Gauss–Laguerre nodes/weights   │ MQ-RBF collocation                        │
│ Build G(b(τ)) matrix           │ Build Φ, Φ_y, Φ_yy evaluation matrices    │
└───────────────────────────────┴──────────────────────────────────────────┘
                                       ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Time loop j = 0..N1−1                                                      │
│  • Predictor step (Crank–Nicolson form)                                    │
│  • Corrector iterations μ = 1..M                                           │
│  • Enforce boundary constraint: V(0,τ) − V_y(0,τ) = 1                      │
│  • Update b(τ)=1−V(0,τ)                                                    │
│  • Project: V := max(V, intrinsic)                                         │
└──────────────────────────────────────┬───────────────────────────────────┘
                                       ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Outputs: option values V at requested S0 (mapped into y); boundary b(τ)   │
└──────────────────────────────────────────────────────────────────────────┘
```

### 6.2 Predictor–Corrector loop structure

```text
┌──────────────────────────────────────────┐
│ Given: λ^j, V^j, b^j                     │
└───────────────────┬──────────────────────┘
                    ▼
┌──────────────────────────────────────────┐
│ Predictor: compute provisional λ^{j+1}   │
│ (uses CN averaging + jump operator G)    │
└───────────────────┬──────────────────────┘
                    ▼
┌──────────────────────────────────────────┐
│ Corrector iterations μ=1..M              │
│  • update V^{j+1,μ}                      │
│  • update b^{j+1,μ} = 1 − V^{j+1,μ}(0)   │
└───────────────────┬──────────────────────┘
                    ▼
┌──────────────────────────────────────────┐
│ Projection: V^{j+1} := max(V^{j+1},payoff)│
└───────────────────┬──────────────────────┘
                    ▼
┌──────────────────────────────────────────┐
│ Advance: (λ^j,b^j,V^j) ← (λ^{j+1},b^{j+1},V^{j+1}) │
└──────────────────────────────────────────┘
```

## 7. Follow-Up Works & Extensions

### 7.1 Direct extensions using localized RBF discretizations (conditioning and sparsity)

**Haghi–Mollapourasl–Vanmaele (2018)** propose a local **RBF-FD** discretization for European and American options under jump-diffusion (Merton and Kou) and explicitly motivate it as a way to avoid the dense/ill-conditioned linear systems of global RBF collocation; they combine RBF-FD with a three-level implicit–explicit time discretization treating the jump integral explicitly, and apply operator splitting for the American LCP formulation. Their paper’s reference list includes the 2012 Golbabai–Ahmadian–Milev article, making it a verifiable citing follow-up. [Haghi et al., *Computers & Mathematics with Applications* 2018]. ([biblio.ugent.be](https://biblio.ugent.be/publication/8593878/file/8593879.pdf))  

**Milovanović–von Sydow (2018)** develop RBF-FD schemes for option pricing problems (including American options via operator splitting) and focus on robust shape-parameter selection strategies and non-uniform node layouts for higher-dimensional problems, addressing a key practical issue emphasized but not resolved analytically in the 2012 paper (shape parameter tuning). Their results highlight how conditioning deteriorates under constant shape parameters as the node spacing shrinks, motivating $$\epsilon \propto h^{-1}$$-type rules. [Milovanović and von Sydow, *Computers & Mathematics with Applications* 2018]. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0898122117307290))  

### 7.2 RBF collocation for jump diffusions with additional market structure (regime switching)

**Bastani–Ahmadi–Damircheli (2013)** extend meshfree RBF collocation to **American options under regime-switching jump-diffusion**, resulting in a coupled system of PIDEs with free-boundary features. The abstract reports a superlinear convergence order in space and linear order in time, positioning the approach as a meshfree analogue to more traditional discretizations for regime-switching PIDEs. [Bastani et al., *Applied Numerical Mathematics* 2013]. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0168927412001961))  

### 7.3 Alternative meshless formulations for the free boundary (fixed-boundary conversions + stability analysis)

**Amani Rad–Parand–Abbasbandy (2015; preprint 2014)** introduce local weak form meshless schemes—LBIE (MLS-based) and LRPI (compactly supported RBF-based)—to compute European and American option prices, converting the American free-boundary problem to a fixed-boundary one using Richardson extrapolation and analyzing stability via matrix methods (unconditional stability claims for implicit Euler and Crank–Nicolson discretizations). This line is methodologically adjacent to Golbabai et al. via “free boundary to fixed boundary” transformations, but uses a different local weak form construction and sparse/banded linear algebra. [Amani Rad et al., *Communications in Nonlinear Science and Numerical Simulation* 2015]. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S1007570414003293))  

**Chan–Hubbert (2010)** provide a pre-2012 independent meshless RBF study for one-dimensional jump-diffusion option pricing, including the American put under both Merton and Kou jump models, using a cubic-spline RBF and emphasizing quadrature handling of the improper jump integral and computational-domain truncation for accuracy. Their reported accuracy behavior (second-order in space with first-order in time in the American case, per the abstract summary) provides a concrete comparator to the 2012 paper’s reported accuracy and shows an earlier “RBF + jump PIDE” implementation pathway. [Chan and Hubbert, arXiv 2010]. ([arxiv.org](https://arxiv.org/abs/1011.5650))  

### 7.4 Data-driven / neural approaches to jump-diffusion PIDEs (non-meshless but same PIDE class)

**Georgoulis–Papapantoleon–Smaragdakis (2024)** propose an implicit–explicit minimizing-movement time-stepping approach paired with neural network approximations per time step for European basket options under jump diffusion, with alternative quadrature strategies for the nonlocal integral operator (including sparse-grid Gauss–Hermite variants). This work targets the same governing PIDE class (jump-diffusion pricing) but replaces RBF spatial approximation with deep residual-type networks and quadrature innovations. [Georgoulis et al., arXiv 2024]. ([arxiv.org](https://arxiv.org/abs/2401.06740))  

## 8. Industrial & Real-World Applications

**QuantLib (C++)** is a large, actively maintained open-source quantitative finance library used for instrument valuation and risk management workflows; it provides a broad implementation base that industrial teams often use as a benchmark or foundation for derivative pricing engines, including American-style options priced by numerical methods (finite-difference engines exist in common QuantLib interfaces and bindings). Public deployment scale is not enumerated by the project; public adoption signals include thousands of GitHub stars and active release cadence. [GitHub: lballabio/QuantLib]. ([github.com](https://github.com/lballabio/QuantLib))  

**RQuantLib (R)** exposes parts of QuantLib to R users and includes an American option pricing interface (documentation states finite-difference engines are used for American options in the provided API). This provides a practical pathway for practitioners to run production-like pricing routines from an R environment, albeit not specific to the Golbabai et al. RBF discretization. [GitHub: eddelbuettel/rquantlib]. ([github.com](https://github.com/eddelbuettel/rquantlib))  

**QLNet (.NET / C#)** is a C# financial library derived primarily from QuantLib, enabling similar valuation/risk workflows in .NET production stacks; it is a credible “real-world” implementation substrate for American option valuation routines (again, not specific to the RBF–Laguerre scheme of the 2012 paper). [GitHub: amaggiulli/QLNet]. ([github.com](https://github.com/amaggiulli/QLNet))  

**No verified production deployment identified for the specific 2012 algorithm.** No public GitHub repository or vendor technical report was found that explicitly implements “front-fixing + Laguerre quadrature + global MQ-RBF collocation + predictor–corrector” exactly as specified by Golbabai–Ahmadian–Milev for American puts under Merton-style jump diffusion; verified public implementations in this problem class more commonly adopt finite differences, finite elements/volumes, operator splitting for the LCP, or local RBF-FD variants rather than global MQ collocation. ([biblio.ugent.be](https://biblio.ugent.be/publication/8593878/file/8593879.pdf))  

## 9. Consolidated Reference List

[1] Majid Haghi; Reza Mollapourasl; Michèle Vanmaele. “An RBF-FD method for pricing American options under jump-diffusion models.” *Computers & Mathematics with Applications*, 2018. DOI: `https://doi.org/10.1016/j.camwa.2018.08.040`. ([biblio.ugent.be](https://biblio.ugent.be/publication/8593878))  

[2] Ali Foroush Bastani; Zaniar Ahmadi; Davood Damircheli. “A radial basis collocation method for pricing American options under regime-switching jump-diffusion models.” *Applied Numerical Mathematics*, 2013. DOI: `https://doi.org/10.1016/j.apnum.2012.10.005`. ([cambridge.org](https://www.cambridge.org/core/journals/astin-bulletin-journal-of-the-iaa/article/option-pricing-in-a-jumpdiffusion-model-with-regime-switching/E60C10EBC6B096EEF8562D225E42E9FE))  

[3] Slobodan Milovanović; Lina von Sydow. “Radial Basis Function generated Finite Differences for option pricing problems.” *Computers & Mathematics with Applications*, 2018. DOI: `https://doi.org/10.1016/j.camwa.2017.11.015`. ([dblp.org](https://dblp.org/db/journals/cma/cma75))  

[4] Jamal Amani Rad; Kourosh Parand; Saeid Abbasbandy. “Local weak form meshless techniques based on the radial point interpolation (RPI) method and local boundary integral equation (LBIE) method to evaluate European and American options.” *Communications in Nonlinear Science and Numerical Simulation*, 2015. DOI: `https://doi.org/10.1016/j.cnsns.2014.07.015`. Preprint: `https://arxiv.org/abs/1412.6063`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S1007570414003293))  

[5] Ron T. L. Chan; Simon Hubbert. “A Numerical Study of Radial Basis Function Based Methods for Options Pricing under the One Dimension Jump-diffusion Model.” arXiv, 2010. `https://arxiv.org/abs/1011.5650`. ([arxiv.org](https://arxiv.org/abs/1011.5650))  

[6] Emmanuil H. Georgoulis; Antonis Papapantoleon; Costas Smaragdakis. “A deep implicit-explicit minimizing movement method for option pricing in jump-diffusion models.” arXiv, 2024. `https://arxiv.org/abs/2401.06740`. ([arxiv.org](https://arxiv.org/abs/2401.06740))  

[7] QuantLib contributors. “QuantLib: The QuantLib C++ library.” GitHub repository. `https://github.com/lballabio/QuantLib`. ([github.com](https://github.com/lballabio/QuantLib))  

[8] Dirk Eddelbuettel (and contributors). “RQuantLib: R interface to the QuantLib library.” GitHub repository. `https://github.com/eddelbuettel/rquantlib`. ([github.com](https://github.com/eddelbuettel/rquantlib))  

[9] QLNet contributors. “QLNet C# Library.” GitHub repository. `https://github.com/amaggiulli/QLNet`. ([github.com](https://github.com/amaggiulli/QLNet))

---
Learn more:
1. [https://www.sciencedirect.com/science/article/pii/S0895717711006066](https://www.sciencedirect.com/science/article/pii/S0895717711006066)
2. [https://biblio.ugent.be/publication/8593878/file/8593879.pdf](https://biblio.ugent.be/publication/8593878/file/8593879.pdf)
3. [https://www.sciencedirect.com/science/article/pii/S0898122117307290](https://www.sciencedirect.com/science/article/pii/S0898122117307290)
4. [https://www.sciencedirect.com/science/article/pii/S0168927412001961](https://www.sciencedirect.com/science/article/pii/S0168927412001961)
5. [https://www.sciencedirect.com/science/article/abs/pii/S1007570414003293](https://www.sciencedirect.com/science/article/abs/pii/S1007570414003293)
6. [https://arxiv.org/abs/1011.5650](https://arxiv.org/abs/1011.5650)
7. [https://arxiv.org/abs/2401.06740](https://arxiv.org/abs/2401.06740)
8. [https://github.com/lballabio/QuantLib](https://github.com/lballabio/QuantLib)
9. [https://github.com/eddelbuettel/rquantlib](https://github.com/eddelbuettel/rquantlib)
10. [https://github.com/amaggiulli/QLNet](https://github.com/amaggiulli/QLNet)
11. [https://biblio.ugent.be/publication/8593878](https://biblio.ugent.be/publication/8593878)
12. [https://www.cambridge.org/core/journals/astin-bulletin-journal-of-the-iaa/article/option-pricing-in-a-jumpdiffusion-model-with-regime-switching/E60C10EBC6B096EEF8562D225E42E9FE](https://www.cambridge.org/core/journals/astin-bulletin-journal-of-the-iaa/article/option-pricing-in-a-jumpdiffusion-model-with-regime-switching/E60C10EBC6B096EEF8562D225E42E9FE)
13. [https://dblp.org/db/journals/cma/cma75](https://dblp.org/db/journals/cma/cma75)

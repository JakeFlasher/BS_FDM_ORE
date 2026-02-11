## 1. Paper Identity

- **Title:** Finite difference techniques for arbitrage-free SABR ([risk.net](https://www.risk.net/journal-of-computational-finance/2465429/finite-difference-techniques-for-arbitrage-free-sabr))  
- **Authors:** Fabien Le Floc’h; Gary Kennedy ([risk.net](https://www.risk.net/journal-of-computational-finance/2465429/finite-difference-techniques-for-arbitrage-free-sabr))  
- **Affiliations (as listed on SSRN):**  
  - Fabien Le Floc’h: Calypso Technology; Independent Researcher ([papers.ssrn.com](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2402001))  
  - Gary J. Kennedy: Clarus Financial Technology ([papers.ssrn.com](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2402001))  
- **Venue / publication metadata (journal version):** Journal of Computational Finance, Volume 20, Number 3 (February 2017), pp. 51–79; first published 15 August 2016. ([risk.net](https://www.risk.net/journal-of-computational-finance/2465429/finite-difference-techniques-for-arbitrage-free-sabr))  
- **DOI (journal):** 10.21314/JCF.2016.320 ([risk.net](https://www.risk.net/journal-of-computational-finance/2465429/finite-difference-techniques-for-arbitrage-free-sabr))  
- **Preprint (SSRN):** SSRN abstract 2402001; DOI 10.2139/ssrn.2402001; posted 28 Feb 2014; last revised 13 Jan 2015; date written 8 May 2014. ([papers.ssrn.com](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2402001))  
- **Versioning ambiguity (from the provided PDF text):** the provided manuscript header states “v1.0 released October 2015” and shows a date line “October 7, 2015”; these dates do not match the SSRN “posted/revised” metadata nor the JCF “first published” date. This document treats the provided PDF content as the authoritative technical source and uses Risk/SSRN only for bibliographic fields.

## 2. Problem Statement & Formulation

Classic SABR implied-volatility approximations (notably the analytic asymptotic formula used in practice) can imply a *negative* risk-neutral density at low strikes and long maturities, creating static arbitrage in low-rate regimes. The paper targets an arbitrage-free-by-construction numerical method that remains close to the classic SABR smile while preventing negative densities, and that is fast enough to be used repeatedly across large rate-option grids (e.g., many caplet/swaption expiries). The approach discretizes a one-dimensional approximate forward (Fokker–Planck) PDE for the option-implied density and compares time-stepping schemes to eliminate Crank–Nicolson oscillations.

### 2.1 Standard (absorbing-boundary) arbitrage-free SABR density PDE (Hagan-type reduction)

Given SABR parameters $$\alpha,\beta,\rho,\nu$$ and a forward level $$f$$ at expiry $$\tau_{\mathrm{ex}}$$, the reduced one-dimensional density $$Q(T,F)$$ is defined on $$F \in [F_{\min},F_{\max}]$$ with absorbing boundary mass accumulators (left and right). The PDE is:

$$
\frac{\partial Q}{\partial T}(T,F)
=
\frac{\partial^2}{\partial F^2}\Big(M(T,F)\,Q(T,F)\Big),
$$

with boundary absorption rates (as written in the source, expressed as flux limits):

$$
\frac{\partial Q_L}{\partial T}(T)
=
\lim_{F\to F_{\min}}
\frac{\partial}{\partial F}\Big(M(T,F)\,Q(T,F)\Big),
\qquad
\frac{\partial Q_R}{\partial T}(T)
=
\lim_{F\to F_{\max}}
\frac{\partial}{\partial F}\Big(M(T,F)\,Q(T,F)\Big).
$$

The diffusion coefficient is defined by

$$
M(T,F) = \frac{1}{2} D(F)^2\,E(T,F),
\qquad
E(T,F) = \exp\big(\rho\,\nu\,\alpha\,\Gamma(F)\,T\big),
\qquad
\Gamma(F)=\frac{F^\beta-f^\beta}{F-f}.
$$

The function $$D(F)$$ is

$$
D(F)
=
\sqrt{\alpha^2 + 2\alpha\rho\nu\,y(F) + \nu^2 y(F)^2}\;F^\beta,
\qquad
y(F)=\frac{F^{1-\beta}-f^{1-\beta}}{1-\beta}.
$$

Initial condition:

$$
\lim_{T\to 0} Q(T,F) = \delta(F-f).
$$

Vanilla prices are then computed from the density and absorbed masses using a Breeden–Litzenberger integral form:

$$
V_{\mathrm{call}}(T,K)
=
\int_{K}^{F_{\max}} (F-K)\,Q(T,F)\,dF
+
(F_{\max}-K)\,Q_R(T),
$$

$$
V_{\mathrm{put}}(T,K)
=
(K-F_{\min})\,Q_L(T)
+
\int_{F_{\min}}^{K} (K-F)\,Q(T,F)\,dF.
$$

### 2.2 Nearly equivalent Dupire forward PDE (call-price formulation)

A closely related forward PDE is written directly on call prices $$V_{\mathrm{call}}(T,F)$$:

$$
\frac{\partial V_{\mathrm{call}}}{\partial T}(T,F)
=
\frac{1}{2}\,\vartheta(T,F)^2\,\frac{\partial^2 V_{\mathrm{call}}}{\partial F^2}(T,F),
\qquad
V_{\mathrm{call}}(0,F)=(f-F)^+,
$$

with a local-volatility approximation of the form $$\vartheta(T,F)=D(F)$$ in the Andreasen–Huge-style simplification, and with boundary-condition choices differing across approaches.

### 2.3 Negative rates and the free-boundary SABR backbone

To avoid an explicit shift parameter while allowing negative forwards, the paper also considers a “free-boundary” SABR backbone:

$$
dF = A\,L(F)\,dW_F,
\qquad
dA = \nu A\,dW_A,
\qquad
dW_F\,dW_A = \rho\,dt,
$$

with

$$
L(F)=|F|^\beta
$$

instead of the shifted-lognormal backbone $$L(F)=(F+b)^\beta$$ or the normal backbone $$L(F)=1$$.

## 3. Core Methodology

**Defined term (first use):** **absorbing boundary mass** refers to the cumulative probability accumulated at the left/right spatial boundaries, tracked explicitly so that the discrete scheme preserves the density’s zeroth moment (total probability) and first moment (forward) under the absorbing-boundary model.

### 3.1 End-to-end computational objective

Inputs:

- Model choice: standard arbitrage-free SABR or free-boundary SABR.
- Parameters: $$\alpha,\beta,\rho,\nu,f,\tau_{\mathrm{ex}}$$ (and discounting/numeraire details externally, since the paper works with undiscounted forward-measure-like quantities in the density integrals).
- Numerical controls: spatial resolution (number of grid points), time steps, boundary placement (either $$F_{\max}$$-based in $$F$$ or $$n_{\mathrm{sd}}$$ in transformed coordinates), and a choice of time integrator (CN, Rannacher, BDF2, Lawson–Swayne, TR-BDF2, etc.).

Outputs:

- A nonnegative (in practice, numerically stable and typically positive) discrete probability density at expiry.
- Absorbed boundary masses.
- Vanilla option prices via quadrature.
- Implied volatilities (Black or normal/Bachelier) as post-processing; the paper reports vol errors but does not specify a unique implied-vol root-finding method.

### 3.2 Pipeline (density formulation with variable transform)

```
┌────────────────────────────────────────────────────────────────────────────┐
│ Inputs: (α,β,ρ,ν,f,τ_ex), model variant (standard / free-boundary), scheme  │
└──────────────┬─────────────────────────────────────────────────────────────┘
               │
               ▼
┌────────────────────────────────────────────────────────────────────────────┐
│ Build transform z(F)=∫_f^F du / D(u); choose bounds z±=±n_sd √τ_ex          │
│ Map: z ↔ y(z) ↔ F(y) ; precompute grid midpoints and coefficients Ĉ, Γ̂     │
└──────────────┬─────────────────────────────────────────────────────────────┘
               │
               ▼
┌────────────────────────────────────────────────────────────────────────────┐
│ Spatial discretization in z: moment-preserving conservative operator L^n    │
│ + absorbing/mirror boundary constraints + boundary-mass ODE discretization  │
└──────────────┬─────────────────────────────────────────────────────────────┘
               │
               ▼
┌────────────────────────────────────────────────────────────────────────────┐
│ Time-march θ (and P_L,P_R) with L-stable 2nd-order scheme (LS or TR-BDF2)   │
│ until T = τ_ex (termination criterion)                                      │
└──────────────┬─────────────────────────────────────────────────────────────┘
               │
               ▼
┌────────────────────────────────────────────────────────────────────────────┐
│ Price calls/puts: midpoint quadrature + special cell around strike;         │
│ enforce exact put-call parity in the discrete formulas                       │
└──────────────┬─────────────────────────────────────────────────────────────┘
               │
               ▼
┌────────────────────────────────────────────────────────────────────────────┐
│ Output: prices, densities, implied vols; diagnostics: oscillations/moments  │
└────────────────────────────────────────────────────────────────────────────┘
```

### 3.3 Free-boundary SABR PDE specialization (negative rates without shift)

The free-boundary version modifies three model ingredients in the density PDE:

1. Backbone:

$$
L(F)=|F|^\beta.
$$

2. Integrated coordinate $$y(F)$$ becomes

$$
y(F)=\int_f^F \frac{du}{L(u)}
=
\frac{\operatorname{sgn}(F)|F|^{1-\beta}-\operatorname{sgn}(f)|f|^{1-\beta}}{1-\beta}.
$$

3. The “gamma-like” ratio uses absolute powers:

$$
\Gamma(F)=\frac{|F|^\beta-|f|^\beta}{F-f}.
$$

The diffusion-factor analogue is

$$
D(F)
=
\sqrt{\alpha^2 + 2\alpha\rho\nu\,y(F) + \nu^2 y(F)^2}\;|F|^\beta,
$$

with

$$
M(T,F) = \frac{1}{2}D(F)^2 E(T,F),
\qquad
E(T,F)=\exp\big(\rho\nu\alpha\Gamma(F)T\big).
$$

Boundary placement differs: $$F_{\min}$$ should be sufficiently negative, with a typical symmetric choice $$F_{\min}=-F_{\max}$$, and the absorbing behavior is retained to preserve zeroth/first moments in the discrete scheme.

### 3.4 Change of variable: transforming the Fokker–Planck PDE

**Defined term (first use):** **transformed density** is $$\theta(T,z)=Q(T,F(z))\,D(F(z))=Q(T,F(z))\,C(z)$$, designed to remove most spatial stiffness from the diffusion coefficient.

Define the Lamperti-type transform

$$
z(F)=\int_f^F \frac{du}{D(u)}.
$$

Let $$F(z)$$ denote the inverse mapping and define

$$
C(z)=D(F(z)).
$$

The transformed PDE is written in conservative form as

$$
\frac{\partial \theta}{\partial T}(T,z)
=
\frac{1}{2}\frac{\partial}{\partial z}
\left(
\frac{1}{C(z)}\frac{\partial}{\partial z}\Big(C(z)\,E(T,z)\,\theta(T,z)\Big)
\right),
$$

with absorbing boundary conditions

$$
\theta(T,z)\to 0 \text{ as } z\to z^-,
\qquad
\theta(T,z)\to 0 \text{ as } z\to z^+,
$$

where $$z^- = z(F_{\min})$$ and $$z^+=z(F_{\max})$$.

The boundary-absorbed probabilities evolve according to flux limits (signs as given in the source):

$$
\frac{\partial P_L}{\partial T}(T)
=
\lim_{z\to z^-}
\frac{1}{2}\frac{1}{C(z)}\frac{\partial}{\partial z}\Big(C(z)E(T,z)\theta(T,z)\Big),
$$

$$
\frac{\partial P_R}{\partial T}(T)
=
\lim_{z\to z^+}
-\frac{1}{2}\frac{1}{C(z)}\frac{\partial}{\partial z}\Big(C(z)E(T,z)\theta(T,z)\Big).
$$

For the standard SABR case, the paper provides explicit closed-form coordinate relations:

$$
y(z)=\frac{\alpha}{\nu}\Big(\sinh(\nu z)+\rho(\cosh(\nu z)-1)\Big),
$$

$$
F(y)=\left(f^{1-\beta}+(1-\beta)y\right)^{\frac{1}{1-\beta}}.
$$

For free-boundary SABR, $$y(z)$$ is unchanged in this coordinate, but $$F(y)$$ changes by inverting the free-boundary $$y(F)$$:

$$
\bar y = y(F=0)= -\frac{\operatorname{sgn}(f)|f|^{1-\beta}}{1-\beta},
$$

$$
F(y)=\operatorname{sgn}(y-\bar y)\left((1-\beta)|y-\bar y|\right)^{\frac{1}{1-\beta}}.
$$

### 3.5 Boundary selection in transformed coordinates

Choose

$$
z^-=-n_{\mathrm{sd}}\sqrt{\tau_{\mathrm{ex}}},\qquad z^+=+n_{\mathrm{sd}}\sqrt{\tau_{\mathrm{ex}}},
$$

with $$n_{\mathrm{sd}}$$ interpreted as “number of standard deviations” in the approximately Gaussian-like transformed coordinate. The paper reports that $$n_{\mathrm{sd}}=4$$ is highly accurate in practice for the free-boundary case, and illustrates that this transform concentrates grid resolution around the forward (located at $$z=0$$).

### 3.6 Moment-preserving spatial discretization in $$z$$

Let a uniform $$z$$-grid be defined by

$$
z_j=z^-+j h,\qquad j=0,\dots,J+1,
\qquad
h=\frac{z^+-z^-}{J+1}.
$$

Define midpoint coordinates (as in the paper’s discretization):

$$
\hat y_j = y\left(z_j-\frac{h}{2}\right),
\qquad
\hat F_j = F(\hat y_j).
$$

Define cached coefficients:

$$
\hat C_j = D(\hat F_j),
\qquad
\hat\Gamma_j=\frac{\hat F_j^\beta-f^\beta}{\hat F_j-f},
\qquad
\hat E_j(T)=\exp\big(\rho\nu\alpha \hat\Gamma_j T\big).
$$

Define a uniform time grid

$$
t_n = n\delta,
\qquad
\delta=\frac{\tau_{\mathrm{ex}}}{N},
\qquad
n=0,\dots,N.
$$

Let $$\theta_j^n=\theta(z_j,t_n)$$.

The semi-discrete form is

$$
\frac{\partial \theta}{\partial T}(z_j,t_n) = \mathcal{L}_j^n\theta(\cdot,t_n),
\qquad j=1,\dots,J,
$$

with the conservative tridiagonal operator (written directly from the paper’s coefficients):

$$
\mathcal{L}_j^n\theta
=
\frac{1}{h}\frac{\hat C_{j-1}}{\hat F_j-\hat F_{j-1}}\hat E_{j-1}(t_n)\,\theta_{j-1}
-\frac{1}{h}\left(
\frac{\hat C_j}{\hat F_{j+1}-\hat F_j}
+
\frac{\hat C_j}{\hat F_j-\hat F_{j-1}}
\right)\hat E_j(t_n)\,\theta_j
+\frac{1}{h}\frac{\hat C_{j+1}}{\hat F_{j+1}-\hat F_j}\hat E_{j+1}(t_n)\,\theta_{j+1}.
$$

Boundary constraints at fictitious endpoints are imposed in “mirror-like” form:

$$
\frac{\hat C_0}{\hat F_1-\hat F_0}\hat E_0(T)\theta_0(T)
=
-\frac{\hat C_1}{\hat F_1-\hat F_0}\hat E_1(T)\theta_1(T),
$$

$$
\frac{\hat C_{J+1}}{\hat F_{J+1}-\hat F_J}\hat E_{J+1}(T)\theta_{J+1}(T)
=
-\frac{\hat C_J}{\hat F_{J+1}-\hat F_J}\hat E_J(T)\theta_J(T).
$$

Absorbed boundary mass rates are discretized as:

$$
\frac{\partial P_L}{\partial T}(T)
=
\frac{\hat C_1}{\hat F_1-\hat F_0}\hat E_1(T)\theta_1(T),
\qquad
\frac{\partial P_R}{\partial T}(T)
=
\frac{\hat C_J}{\hat F_{J+1}-\hat F_J}\hat E_J(T)\theta_J(T).
$$

### 3.7 Pricing from the discrete transformed density (midpoint rule with parity preservation)

For strike $$K$$, define $$z^\ast=z(y(K))$$ and assume $$z^-<z^\ast<z^+$$. Let $$k$$ be the index such that

$$
z^-+(k-1)h < z^\ast \le z^-+kh.
$$

Let

$$
F_k = F\big(y(z^-+kh)\big).
$$

The paper’s discrete call and put prices are:

$$
V_{\mathrm{call}}
=
\frac{h}{4(F_k-\hat F_k)}(F_k-K)^2\,\theta_k
+
\sum_{j=k+1}^{J-1} (\hat F_j-K)\,h\,\theta_j
+
(F_{\max}-K)\,P_R,
$$

$$
V_{\mathrm{put}}
=
\frac{h}{4(\hat F_k-F_k)}(K-F_k)^2\,\theta_k
+
\sum_{j=1}^{k} (K-\hat F_j)\,h\,\theta_j
+
(K-F_{\min})\,P_L.
$$

Edge cases:

- If $$z^\ast \le z^-$$, then $$V_{\mathrm{call}}=f-K$$ and $$V_{\mathrm{put}}=0$$.
- If $$z^\ast \ge z^+$$, then $$V_{\mathrm{call}}=0$$ and $$V_{\mathrm{put}}=K-f$$.

The construction is stated to preserve put–call parity exactly at the discrete level.

### 3.8 Time-stepping schemes (all formulated to preserve moments)

Notation: spatial operator at time level $$t_n$$ is $$\mathcal{L}^n$$; boundary mass updates use the matching discrete flux expressions.

#### 3.8.1 Moment-preserving implicit Euler (IE)

**Procedure (Implicit Euler on $$\theta$$ with boundary masses):**

1. **Inputs:** $$\theta^n$$, $$P_L(t_n)$$, $$P_R(t_n)$$; coefficients defining $$\mathcal{L}^{n+1}$$ and flux terms at $$t_{n+1}$$; step size $$\delta$$.
2. **Solve** the tridiagonal linear system for interior nodes $$j=1,\dots,J$$:
   $$
   \theta_j^{n+1}-\theta_j^{n} = \delta\,\mathcal{L}_j^{n+1}\theta^{n+1}.
   $$
3. **Update** absorbed masses:
   $$
   P_L(t_{n+1})-P_L(t_n)
   =
   \delta\,
   \frac{\hat C_1}{\hat F_1-\hat F_0}\hat E_1(t_{n+1})\theta_1^{n+1},
   $$
   $$
   P_R(t_{n+1})-P_R(t_n)
   =
   \delta\,
   \frac{\hat C_J}{\hat F_{J+1}-\hat F_J}\hat E_J(t_{n+1})\theta_J^{n+1}.
   $$
4. **Outputs:** $$\theta^{n+1}$$, $$P_L(t_{n+1})$$, $$P_R(t_{n+1})$$.
5. **Termination:** repeat until $$t_N=\tau_{\mathrm{ex}}$$.

Implementation notes given by the paper for the untransformed $$Q(F)$$ discretization:

- When the left boundary is near $$0$$, a fictitious grid point $$F_0=F_{\min}-h/2$$ may be negative, which can make model functions ill-defined; the paper notes only the product $$M_0Q_0$$ is needed, and mirror-like boundary condition enforces:
  $$
  M_0Q_0 + M_1Q_1 = 0.
  $$
- The ratio $$\Gamma(F)$$ is undefined at $$F=f$$; the paper uses $$\Gamma(f)=\frac{\partial C}{\partial F}(f)$$.

#### 3.8.2 Moment-preserving Crank–Nicolson (CN)

**Procedure (Crank–Nicolson on $$\theta$$ with boundary masses):**

1. **Inputs:** $$\theta^n$$, $$P_L(t_n)$$, $$P_R(t_n)$$; operators $$\mathcal{L}^{n}$$ and $$\mathcal{L}^{n+1}$$; step size $$\delta$$.
2. **Solve** for $$\theta^{n+1}$$:
   $$
   \theta_j^{n+1}-\theta_j^{n}
   =
   \frac{\delta}{2}\left(
   \mathcal{L}_j^{n+1}\theta^{n+1}+\mathcal{L}_j^{n}\theta^{n}
   \right).
   $$
3. **Update** absorbed masses with trapezoidal averaging:
   $$
   P_L(t_{n+1})-P_L(t_n)
   =
   \frac{\delta}{2}\frac{\hat C_1}{\hat F_1-\hat F_0}
   \left(
   \hat E_1(t_{n+1})\theta_1^{n+1}+\hat E_1(t_{n})\theta_1^{n}
   \right),
   $$
   $$
   P_R(t_{n+1})-P_R(t_n)
   =
   \frac{\delta}{2}\frac{\hat C_J}{\hat F_{J+1}-\hat F_J}
   \left(
   \hat E_J(t_{n+1})\theta_J^{n+1}+\hat E_J(t_{n})\theta_J^{n}
   \right).
   $$
4. **Outputs:** $$\theta^{n+1}$$, $$P_L(t_{n+1})$$, $$P_R(t_{n+1})$$.

#### 3.8.3 Rannacher smoothing (RAN)

**Procedure (Rannacher time-stepping, moment-preserving):**

1. **Inputs:** initial $$\theta^0$$ and boundary masses; total number of full steps $$N$$ with size $$\delta$$.
2. **Perform** two half-steps of implicit Euler (size $$\delta/2$$) to damp initial-condition irregularity:
   $$
   \theta_j^{n+\frac{1}{2}}-\theta_j^{n}
   =
   \frac{\delta}{2}\,\mathcal{L}_j^{n+\frac{1}{2}}\theta^{n+\frac{1}{2}},
   $$
   with matching boundary-mass half-step updates using $$t_{n+\frac{1}{2}}$$.
3. **Switch** to Crank–Nicolson for subsequent full steps $$n=2,\dots,N-1$$, using the CN updates above.
4. **Outputs:** $$\theta^N$$, $$P_L(t_N)$$, $$P_R(t_N)$$.

#### 3.8.4 Two-step BDF2 (BDF2)

**Procedure (BDF2 with IE initialization):**

1. **Inputs:** $$\theta^0$$, $$P_L(t_0)$$, $$P_R(t_0)$$; step size $$\delta$$; operators at needed times.
2. **Initialize** $$\theta^1$$, $$P_L(t_1)$$, $$P_R(t_1)$$ via implicit Euler (Section 3.8.1).
3. **For** $$n=0,\dots,N-2$$, solve for $$\theta^{n+2}$$:
   $$
   3\theta_j^{n+2}-4\theta_j^{n+1}+\theta_j^{n}
   =
   2\delta\,\mathcal{L}_j^{n+2}\theta^{n+2}.
   $$
4. **Update** boundary masses with the matching BDF2 stencil:
   $$
   3P_L(t_{n+2})-4P_L(t_{n+1})+P_L(t_{n})
   =
   2\delta\,
   \frac{\hat C_1}{\hat F_1-\hat F_0}\hat E_1(t_{n+2})\theta_1^{n+2},
   $$
   and similarly for $$P_R$$ with index $$J$$.
5. **Outputs:** $$\theta^N$$, $$P_L(t_N)$$, $$P_R(t_N)$$.

#### 3.8.5 Implicit Richardson extrapolation (RE) applied to IE

**Procedure (global Richardson on implicit Euler):**

1. **Inputs:** initial state; target expiry $$T=\tau_{\mathrm{ex}}$$; base time step $$\delta$$.
2. **Compute** $$\bar\theta_{\delta}(z)$$, $$\bar P_{L,\delta}$$, $$\bar P_{R,\delta}$$ by implicit Euler with step $$\delta$$ to time $$T$$.
3. **Compute** $$\bar\theta_{\delta/2}(z)$$, $$\bar P_{L,\delta/2}$$, $$\bar P_{R,\delta/2}$$ by implicit Euler with step $$\delta/2$$ (twice as many steps) to time $$T$$.
4. **Extrapolate** at $$T$$:
   $$
   \theta(z)=2\bar\theta_{\delta/2}(z)-\bar\theta_{\delta}(z),
   \qquad
   P_L = 2\bar P_{L,\delta/2}-\bar P_{L,\delta},
   \qquad
   P_R = 2\bar P_{R,\delta/2}-\bar P_{R,\delta}.
   $$
5. **Outputs:** extrapolated $$\theta, P_L, P_R$$ at expiry.

#### 3.8.6 Lawson–Morris–Gourlay (LMG2, LMG3): local extrapolation based on IE

The paper describes LMG as a local (per-step) Richardson-type extrapolation that reuses the same tridiagonal matrix factorization.

**Procedure (LMG2):**

1. **At each full step** of size $$\delta$$ from $$t_n$$ to $$t_{n+1}$$, compute implicit-Euler results with step sizes $$\delta$$ and $$\delta/2$$ over the same interval.
2. **Apply** the second-order extrapolation (same form as global RE but per step):
   $$
   \theta^{n+1}=2\bar\theta_{\delta/2}^{n+1}-\bar\theta_{\delta}^{n+1},
   $$
   with the analogous linear combination for $$P_L,P_R$$.
3. **Repeat** to reach expiry.

**Procedure (LMG3):**

1. **Within each full interval** $$[t_n,t_{n+1}]$$, compute implicit Euler on substeps $$\delta/3$$ (three substeps) and on a mixed grid $$2\delta/3$$ then $$\delta/3$$.
2. **Combine** (as given):
   $$
   \theta = 4.5\,\bar\theta_{\delta/3} - 4.5\,\bar\theta_{2\delta/3} + \bar\theta_{\delta},
   $$
   with the same coefficients for $$P_L,P_R$$.
3. **Outputs:** third-order local extrapolated step values.

#### 3.8.7 Lawson–Swayne (LS)

Let

$$
b = 1-\frac{\sqrt{2}}{2}.
$$

**Procedure (Lawson–Swayne second-order scheme):**

1. **Inputs:** $$\theta^n, P_L(t_n), P_R(t_n)$$; step size $$\delta$$.
2. **Stage 1:** implicit Euler over $$b\delta$$ to obtain $$\theta^{n+b}$$ and updated boundary masses:
   $$
   \theta_j^{n+b}-\theta_j^{n} = b\delta\,\mathcal{L}_j^{n+b}\theta^{n+b},
   $$
   $$
   P_L(t_{n+b})-P_L(t_n) = b\delta\,\frac{\hat C_1}{\hat F_1-\hat F_0}\hat E_1(t_{n+b})\theta_1^{n+b},
   $$
   and similarly for $$P_R$$.
3. **Stage 2:** implicit Euler over another $$b\delta$$ to obtain $$\theta^{n+2b}$$ from $$\theta^{n+b}$$, updating masses analogously.
4. **Extrapolation to full step:** combine stage values:
   $$
   \theta_j^{n+1} = (\sqrt{2}+1)\theta_j^{n+2b}-\sqrt{2}\,\theta_j^{n+b},
   $$
   $$
   P_L(t_{n+1}) = (\sqrt{2}+1)P_L(t_{n+2b})-\sqrt{2}\,P_L(t_{n+b}),
   \qquad
   P_R(t_{n+1}) = (\sqrt{2}+1)P_R(t_{n+2b})-\sqrt{2}\,P_R(t_{n+b}).
   $$
5. **Outputs:** $$\theta^{n+1}, P_L(t_{n+1}), P_R(t_{n+1})$$.

#### 3.8.8 TR-BDF2

TR-BDF2 is described as a two-stage method: a (weighted) trapezoidal-rule stage followed by a BDF2-like stage, second order and $$L$$-stable.

Let the scheme weight be $$\alpha_{\mathrm{TR}}$$ (the paper uses the symbol $$\alpha$$; this document uses $$\alpha_{\mathrm{TR}}$$ to avoid collision with SABR’s $$\alpha$$ parameter). The paper states two choices: $$\alpha_{\mathrm{TR}}=\frac{1}{2}$$ (CN-like) or $$\alpha_{\mathrm{TR}}=2-\sqrt{2}$$ (proportional Jacobians; “optimal stability”).

**Procedure (TR-BDF2 full step):**

1. **Inputs:** $$\theta^n, P_L(t_n), P_R(t_n)$$; step size $$\delta$$.
2. **Stage 1 (TR/CN to intermediate time):**
   $$
   \theta_j^{n+\alpha_{\mathrm{TR}}}-\theta_j^n
   =
   \frac{\alpha_{\mathrm{TR}}\delta}{2}
   \left(
   \mathcal{L}_j^{n+\alpha_{\mathrm{TR}}}\theta^{n+\alpha_{\mathrm{TR}}}
   +
   \mathcal{L}_j^{n}\theta^{n}
   \right).
   $$
   The paper also provides intermediate boundary-mass update formulas; the provided text contains a symbol ambiguity in the coefficient multiplying $$\delta$$ (it appears as $$b\delta$$, plausibly an OCR error for $$\alpha_{\mathrm{TR}}\delta/2$$). This ambiguity is flagged; moment preservation requires applying the *same* trapezoidal weighting used in the interior equation to the boundary flux.
3. **Stage 2 (BDF2-type completion to $$t_{n+1}$$):**
   $$
   \theta_j^{n+1}
   =
   \frac{1}{2-\alpha_{\mathrm{TR}}}
   \left(
   \frac{1}{\alpha_{\mathrm{TR}}}\theta_j^{n+\alpha_{\mathrm{TR}}}
   -
   \frac{1-\alpha_{\mathrm{TR}}}{2\alpha_{\mathrm{TR}}}\theta_j^{n}
   +
   \delta(1-\alpha_{\mathrm{TR}})\mathcal{L}_j^{n+1}\theta^{n+1}
   \right),
   $$
   with analogous boundary-mass updates consistent with the second-stage stencil.
4. **Outputs:** $$\theta^{n+1}, P_L(t_{n+1}), P_R(t_{n+1})$$.

#### 3.8.9 Three-stage “Bathe” variant (two TR stages + BDF3-like stage)

The paper describes a three-stage extension (named “Bathe”) with stronger damping:

$$
\theta_j^{n+\frac{1}{3}}
=
\theta_j^{n}
+
\frac{\delta}{6}\left(\mathcal{L}_j^{n}\theta^{n}+\mathcal{L}_j^{n+\frac{1}{3}}\theta^{n+\frac{1}{3}}\right),
$$

$$
\theta_j^{n+\frac{2}{3}}
=
\theta_j^{n}
+
\frac{\delta}{6}\left(\mathcal{L}_j^{n+\frac{1}{3}}\theta^{n+\frac{1}{3}}+\mathcal{L}_j^{n+\frac{2}{3}}\theta^{n+\frac{2}{3}}\right),
$$

$$
\theta_j^{n+1}
=
\frac{1}{11}\left(
18\theta_j^{n+\frac{2}{3}}
-
9\theta_j^{n+\frac{1}{3}}
+
2\theta_j^{n}
+
2\delta\,\mathcal{L}_j^{n+1}\theta^{n+1}
\right).
$$

The paper does not print the corresponding $$P_L,P_R$$ formulas in the excerpted section; moment preservation requires updating boundary masses with the same linear combinations applied to the interior implicit-Euler solves.

### 3.9 Coordinate transformation for the Dupire forward PDE (non-uniform grid approach)

A uniform grid in $$z$$ induces a non-uniform grid in $$F$$:

$$
F_j = F\big(y(z_j)\big).
$$

The paper describes an additional grid shift to place the forward $$f$$ midway between two consecutive nodes:

- Find $$j_0=\left\lfloor \frac{z(y(f))-z^-}{h}\right\rfloor$$ such that $$F_{j_0}\le f < F_{j_0+1}$$.
- Define a shift $$d$$ that maps the midpoint $$\frac{1}{2}(F_{j_0}+F_{j_0+1})$$ back to $$f$$ in $$z$$ coordinates, then rebuild the grid as $$\tilde F_j = F(y(z^-+jh+d))$$ while keeping the original endpoints fixed.
- Interpolate option prices between nodes via a natural cubic spline to preserve continuity of the second derivative (and thus the implied density).

### 3.10 Performance optimization for the time-dependent exponential factor

Because

$$
E(T,F)=\exp\big(\rho\nu\alpha\Gamma(F)T\big),
$$

the paper notes that:

- $$D(F)$$ and $$\Gamma(F)$$ can be cached (time-independent).
- If $$t_{n+1}=t_n+\delta$$, then
  $$
  \exp(\rho\nu\alpha\Gamma\,t_{n+1})=\exp(\rho\nu\alpha\Gamma\,t_n)\exp(\rho\nu\alpha\Gamma\,\delta).
  $$
- Precompute $$e_j=\exp(\rho\nu\alpha\Gamma(F_j)\delta)$$ once, then update:
  $$
  E_j^{n+1} = e_j\,E_j^{n},
  \qquad E_j^{0}=1.
  $$

The paper reports that treating $$E$$ as piecewise constant over fractional stages yields small speed gains but noticeably worse accuracy for long maturities and large steps, so this approximation is not used in the reported tests.

## 4. Theoretical Results

### 4.1 Formal theorem/lemma/proposition/corollary inventory

No statements are labeled as “Theorem”, “Lemma”, “Proposition”, or “Corollary” within the provided paper text. A single external theorem is invoked by citation (Morton & Mayers, Theorem 2.2) as a sufficient non-oscillation condition for parabolic schemes.

### 4.2 Explicit guarantees and referenced stability results used by the paper (not formally numbered)

#### Claim A (moment preservation of the conservative spatial discretization)

**Statement (as used):** the discrete operator $$\mathcal{L}$$ together with absorbing/mirror boundary conditions and explicit boundary-mass tracking preserves the zeroth and first moments (total probability and forward) exactly in the discrete evolution.

**Conditions:** absorbing behavior at both boundaries; boundary mass updates use the matching discrete flux; scheme uses the conservative form (tridiagonal flux-difference structure) in $$z$$.

**Proof strategy sketch:** summing the semi-discrete equations over $$j$$ turns the tridiagonal stencil into telescoping boundary flux terms; adding the explicit boundary-mass ODE discretization cancels those fluxes to yield exact conservation of total mass. Multiplying by appropriate discrete approximations of $$F$$ (or $$\hat F_j$$ in the transformed-midpoint grid) yields the corresponding telescoping identity for the first moment; the absorbing boundary flux accounts for the mass/forward carried out of the computational domain.

#### Claim B (Crank–Nicolson oscillations explained via a Courant-number condition)

**Statement (as used):** to guarantee the absence of oscillations, the Courant number must satisfy $$\Psi \le 1$$ (cited to Morton & Mayers, Theorem 2.2), and CN oscillations in the SABR density PDE can be interpreted via the explicit-Euler component embedded in CN.

**Conditions:** parabolic diffusion equation with non-smooth initial data (delta-like), discretized on a uniform grid with large time steps.

**Proof strategy sketch:** CN can be decomposed into an averaging of explicit and implicit Euler updates; if the explicit part violates a monotonicity constraint (large $$\delta/h^2$$ relative to the diffusion scale), spurious oscillations appear near sharp initial features. The paper quantifies this using a Courant-number proxy evaluated near-the-money.

The paper’s specific Courant proxies are:

$$
\Psi_Q = \frac{M\,\delta}{h^2}
\approx
\frac{1}{2}\alpha^2 f^{2\beta}\frac{\delta}{h^2},
$$

$$
\Psi_\theta = f^\beta\frac{\delta}{h^2}.
$$

#### Claim C (relative stability of L-stable schemes on the SABR density PDE)

**Statement (as used):** TR-BDF2 and Lawson–Swayne (and related L-stable methods such as BDF2) provide strong damping of oscillations and allow accurate solutions with far fewer time steps than CN on this problem.

**Conditions:** stiff parabolic PDE with delta-like initial condition; large time steps where CN exhibits oscillations.

**Proof strategy sketch:** L-stable schemes heavily damp high-frequency components introduced by discontinuous/peaked initial data or by spatial discretization errors, while maintaining second-order accuracy in time; this reduces the need to enforce $$\Psi \le 1$$ by taking very small $$\delta$$. The paper’s numerical section provides empirical evidence via density and implied-vol error curves.

#### Complexity bounds

No explicit asymptotic runtime bound is stated. Each implicit stage results in a tridiagonal linear system per time step (or per stage), implying $$O(J)$$ work per solve with a Thomas-type solver; this is standard for 1D second-order finite differences but is not stated as a theorem in the paper.

## 5. Experimental Evaluation

### 5.1 Experimental setups (test cases, baselines, metrics, key hyperparameters)

| Test case (paper section) | Model / formulation | Parameters | Maturity | Spatial grid & boundaries | Time stepping / baselines | Metrics reported |
|---|---|---|---|---|---|---|
| CN oscillation demo (Sec. 7) | Standard arbitrage-free SABR; density PDE in $$Q(F)$$ and in transformed $$\theta(z)$$ | $$\alpha=35\%$$, $$\beta=0.25$$, $$\rho=-10\%$$, $$\nu=100\%$$, $$f=1$$ | $$\tau_{\mathrm{ex}}=1$$ | Example: 500 points in rate dimension; $$F_{\max}=5$$ for $$Q$$, and $$n_{\mathrm{sd}}=4$$ for $$\theta$$ | CN | Qualitative oscillations; Courant proxies $$\Psi_Q,\Psi_\theta$$ |
| Variable-change efficiency (Sec. 5.1; Table 1) | Standard arbitrage-free SABR; Lawson–Swayne; compare uniform $$F$$ grid vs uniform $$z$$ grid | Extreme: $$\alpha=100\%$$, $$\beta=0.30$$, $$\rho=90\%$$, $$\nu=100\%$$, $$\tau_{\mathrm{ex}}=10$$, $$f=1$$ | 10 | Uniform in $$F$$ uses varying $$F_{\max}$$; transformed uses $$n_{\mathrm{sd}}$$ | Lawson–Swayne | Call price and implied vol vs discretization size |
| Performance vs steps (Sec. 9.2.1; Fig. 5) | Density PDE in $$\theta$$ | Hagan example above | 1 | Reference: CN with 5120 points and enough steps to ensure $$\Psi_\theta<1$$ | Bathe, BDF2, CN, LMG2, LS, RAN, RE, TR-BDF2 | Max error in $$\theta$$ (proxy), implied-vol error thresholds |
| Andreasen–Huge parameter set (Sec. 9.2.2; Fig. 6) | Density PDE in $$\theta$$ | $$\alpha=0.0873$$, $$\beta=0.7$$, $$\rho=-0.48$$, $$\nu=0.47$$, $$f=0.025$$ | $$\tau_{\mathrm{ex}}=10$$ | Error measured over strikes $$[0.2f,2f]$$ | Same scheme set | Max implied-vol error over strike range |
| Dupire-vs-density comparison (Sec. 9.3; Figs. 7–8) | Density PDE vs Dupire forward PDE with matched coordinates | Andreasen–Huge parameters | 10 | Example: Lawson–Swayne with $$J=50$$, $$n_{\mathrm{sd}}=4$$ | Density PDE vs Dupire PDE | Max implied-vol difference; density oscillations at low step counts |
| Free-boundary behavior (Sec. 10; Table 2; Fig. 9–10) | Free-boundary SABR density PDE | Example: $$f=50$$ b.p., $$\beta=0.1$$, $$\alpha=0.5 f^{1-\beta}$$, $$\rho=-30\%$$, $$\nu=30\%$$ | $$\tau_{\mathrm{ex}}=3$$ | $$n_{\mathrm{sd}}=4$$ reported accurate | Lawson–Swayne; compare to Antonov et al Monte Carlo references | b.p. vols vs strike; spike at $$F=0$$; knee behavior |

### 5.2 Key quantitative results (reproduced tables)

#### Table 1 (variable change vs uniform discretization; Lawson–Swayne; extreme parameters)

Extreme parameters: $$\alpha=100\%$$, $$\beta=0.30$$, $$\rho=90\%$$, $$\nu=100\%$$, $$\tau_{\mathrm{ex}}=10$$, $$f=1$$.

| Uniform discretization of $$Q$$ in $$F$$: $$F_{\max}$$ | Points | Steps | Price | Vol |
|---:|---:|---:|---:|---:|
| 5 | 10 | 5 | 0.65010 | 87.205 |
| 50 | 100 | 10 | 0.78769 | 155.773 |
| 500 | 1000 | 20 | 0.79782 | 191.658 |
| 5000 | 10000 | 160 | 0.79835 | 196.930 |

| Discretization of $$\theta$$ in $$z$$: $$n_{\mathrm{sd}}$$ | Points | Steps | Price | Vol |
|---:|---:|---:|---:|---:|
| 3 | 10 | 5 | 0.79848 | 198.504 |
| 3 | 100 | 10 | 0.79853 | 199.148 |
| 4 | 100 | 20 | 0.79847 | 198.338 |
| 10 | 10000 | 160 | 0.79845 | 198.134 |

Reported conclusion: a uniform $$F$$ grid may require roughly $$1000\times$$ more points than the transformed-coordinate grid to reach similar accuracy in extreme cases; the practical constraint is to keep the forward away from the boundary, with this constraint being much stricter for uniform $$F$$ discretizations.

#### Table 2 (free-boundary SABR b.p.-vol accuracy vs Monte Carlo reference)

Parameters: $$f=50$$ b.p., $$\beta=0.25$$, $$\alpha=0.5 f^{1-\beta}$$, $$\rho=-30\%$$, $$\nu=30\%$$, $$\tau_{\mathrm{ex}}=3$$.

| $$K$$ | Lawson–Swayne | Reference | Difference |
|---:|---:|---:|---:|
| -0.95 | 30.72 | 30.93 | -0.21 |
| -0.80 | 29.70 | 29.95 | -0.25 |
| -0.65 | 28.69 | 28.97 | -0.28 |
| -0.50 | 27.70 | 27.99 | -0.29 |
| -0.35 | 26.74 | 27.04 | -0.30 |
| -0.20 | 25.87 | 26.15 | -0.28 |
| -0.05 | 25.21 | 25.46 | -0.25 |
| 0.10 | 25.71 | 25.85 | -0.14 |
| 0.25 | 26.65 | 26.69 | -0.04 |
| 0.40 | 27.43 | 27.39 | 0.04 |
| 0.55 | 28.06 | 27.97 | 0.09 |
| 0.70 | 28.58 | 28.45 | 0.13 |
| 0.85 | 29.03 | 28.87 | 0.16 |
| 1.00 | 29.42 | 29.25 | 0.17 |
| 1.15 | 29.79 | 29.60 | 0.19 |
| 1.30 | 30.14 | 29.94 | 0.20 |
| 1.45 | 30.48 | 30.29 | 0.19 |
| 1.60 | 30.84 | 30.63 | 0.21 |
| 1.75 | 31.20 | 30.99 | 0.21 |
| 1.90 | 31.58 | 31.37 | 0.21 |

Reported summary: maximum error vs the cited Monte Carlo references is $$0.30$$ b.p. for the PDE and $$0.20$$ b.p. for Antonov et al’s analytic approximation; the density shows a spike at $$F=0$$ that is inherent to the model rather than a propagated numerical oscillation.

#### Table B1 (reference implementation values; 500 points; 5 time-steps)

Parameters: $$\alpha=35\%$$, $$\beta=0.25$$, $$\rho=-10\%$$, $$\nu=100\%$$, $$T=1$$, $$f=1$$, $$n_{\mathrm{sd}}=4$$, $$\delta=0.2$$, $$h=0.012018637349$$.

| Scheme | ATM price | $$\theta(0)$$ | $$P_L$$ | $$P_R$$ |
|---|---:|---:|---:|---:|
| CN | 0.156227316001 | -75.391631075105 | 0.036151920718 | 0.000014227771 |
| RAN | 0.149166031026 | 0.486588975069 | 0.037035726447 | 0.000023444260 |
| BDF2 | 0.149369112184 | 0.478480554553 | 0.036571170374 | 0.000036350018 |
| RE | 0.149622595293 | 0.482424678160 | 0.036971313633 | -0.000001440333 |
| LMG2 | 0.149449019862 | 0.486727660422 | 0.037356585469 | -0.000003103630 |
| LS | 0.149701955629 | 0.482422521404 | 0.036472664324 | 0.000011244607 |
| TRBDF2 | 0.149703527234 | 0.482401023656 | 0.036469263805 | 0.000011230925 |
| Bathe | 0.149631007454 | 0.486420051293 | 0.036725562889 | 0.000008932123 |

Reported note: extrapolation schemes (RE, LMG2) can generate a very small negative accumulated probability at the upper boundary for these coarse settings, attributed to small oscillations near the boundary; values converge rapidly as time steps increase.

#### Table C1 (convergence table excerpt: ATM implied volatility, selected schemes)

Caption parameters: $$\alpha=35\%$$, $$\beta=0.25$$, $$\rho=-10\%$$, $$\nu=100\%$$, $$T=1$$, $$f=1$$, $$n_{\mathrm{sd}}=4$$.

Part 1 (explicit $$J,N$$ shown):

| Scheme | $$J$$ | $$N$$ | ATM value | Change | Ratio | Time (s) |
|---|---:|---:|---:|---:|---:|---:|
| CN | 80 | 5 | 38.92945097 | NaN | NaN | 1.1e-04 |
| CN | 160 | 10 | 37.13441397 | -1.8e+00 | NaN | 3.0e-04 |
| CN | 320 | 20 | 37.42600113 | 2.9e-01 | -6.2 | 7.8e-04 |
| CN | 640 | 40 | 37.57545976 | 1.5e-01 | 2.0 | 2.1e-03 |
| CN | 1280 | 80 | 37.65103999 | 7.6e-02 | 2.0 | 7.3e-03 |
| CN | 2560 | 160 | 37.68906248 | 3.8e-02 | 2.0 | 2.7e-02 |
| RAN | 80 | 5 | 37.59281281 | NaN | NaN | 1.4e-04 |
| RAN | 160 | 10 | 37.69385914 | 1.0e-01 | NaN | 2.9e-04 |
| RAN | 320 | 20 | 37.71890369 | 2.5e-02 | 4.0 | 6.8e-04 |
| RAN | 640 | 40 | 37.72510265 | 6.2e-03 | 4.0 | 2.2e-03 |
| RAN | 1280 | 80 | 37.72664749 | 1.5e-03 | 4.0 | 7.6e-03 |
| RAN | 2560 | 160 | 37.72703321 | 3.9e-04 | 4.0 | 2.7e-02 |
| BDF2 | 80 | 5 | 37.64496421 | NaN | NaN | 9.7e-05 |
| BDF2 | 160 | 10 | 37.70648716 | 6.2e-02 | NaN | 2.1e-04 |
| BDF2 | 320 | 20 | 37.72178602 | 1.5e-02 | 4.0 | 7.3e-04 |
| BDF2 | 640 | 40 | 37.72578924 | 4.0e-03 | 3.8 | 2.1e-03 |
| BDF2 | 1280 | 80 | 37.72681485 | 1.0e-03 | 3.9 | 7.3e-03 |
| BDF2 | 2560 | 160 | 37.72707451 | 2.6e-04 | 3.9 | 2.8e-02 |
| RE | 80 | 5 | 37.70951734 | NaN | NaN | 1.8e-04 |
| RE | 160 | 10 | 37.72310519 | 1.4e-02 | NaN | 4.2e-04 |
| RE | 320 | 20 | 37.72621751 | 3.1e-03 | 4.4 | 1.6e-03 |
| RE | 640 | 40 | 37.72693472 | 7.2e-04 | 4.3 | 5.5e-03 |
| RE | 1280 | 80 | 37.72710615 | 1.7e-04 | 4.2 | 2.0e-02 |
| RE | 2560 | 160 | 37.72714796 | 4.2e-05 | 4.1 | 7.2e-02 |

Part 2 (the provided text excerpt omits the $$J,N$$ columns for LMG2/LMG3/LS/TRBDF2; values are reproduced in order as shown, and appear aligned with the same doubling sequence used above, but that alignment should be verified against the formatted PDF table):

| Scheme | ATM value sequence | Change sequence | Ratio sequence | Time (s) sequence |
|---|---|---|---|---|
| LMG2 | 37.66503968; 37.70865944; 37.72199774; 37.72578272; 37.72680412; 37.72707055 | NaN; 4.4e-02; 1.3e-02; 3.8e-03; 1.0e-03; 2.7e-04 | NaN; NaN; 3.3; 3.5; 3.7; 3.8 | 1.4e-04; 3.9e-04; 1.3e-03; 4.9e-03; 1.8e-02; 7.1e-02 |
| LMG3 | 37.70255505; 37.72145184; 37.72587154; 37.72686398; 37.72709114; 37.72714460 | NaN; 1.9e-02; 4.4e-03; 9.9e-04; 2.3e-04; 5.3e-05 | NaN; NaN; 4.3; 4.5; 4.4; 4.2 | 1.8e-04; 5.7e-04; 2.0e-03; 7.6e-03; 2.9e-02; 1.2e-01 |
| LS | 37.72979145; 37.72762970; 37.72729564; 37.72719717; 37.72717084; 37.72716402 | NaN; -2.2e-03; -3.3e-04; -9.8e-05; -2.6e-05; -6.8e-06 | NaN; NaN; 6.5; 3.4; 3.7; 3.9 | 1.2e-04; 2.9e-04; 9.5e-04; 3.4e-03; 1.2e-02; 4.9e-02 |
| TRBDF2 | 37.73019364; 37.72772657; 37.72731963; 37.72720314; 37.72717233; 37.72716439 | NaN; -2.5e-03; -4.1e-04; -1.2e-04; -3.1e-05; -7.9e-06 | NaN; NaN; 6.1; 3.5; 3.8; 3.9 | 1.6e-04; 3.0e-04; 9.5e-04; 3.5e-03; 1.3e-02; 5.0e-02 |

### 5.3 Quantitative conclusions stated in the numerical section (time-step efficiency)

For the Hagan oscillation-parameter example (Sec. 9.2.1), the paper reports that a Black implied-volatility absolute error under $$0.1\%$$ is achieved with:

| Scheme | Time steps needed for $$<0.1\%$$ Black implied-vol abs. error (reported) |
|---|---:|
| Bathe | 3 |
| Lawson–Swayne | 6 |
| TR-BDF2 | 6 |
| LMG2 | 10 |
| BDF2 | 12 |
| Rannacher | 12 |

For the Andreasen–Huge parameter set (Sec. 9.2.2), the paper reports that a Black implied-volatility accuracy better than $$0.1\%$$ is achieved with:

| Scheme | Time steps needed for $$<0.1\%$$ max implied-vol error on $$[0.2f,2f]$$ (reported) |
|---|---:|
| Lawson–Swayne | 3 |
| Bathe | 3 |
| TR-BDF2 | 4 |
| LMG2 | 5 |
| Rannacher | 12 |
| BDF2 | 12 |

Dupire-forward-PDE vs density-PDE comparison (Sec. 9.3): with Lawson–Swayne, 5 time steps and $$J=50$$, $$n_{\mathrm{sd}}=4$$, the implied-vol difference between the two formulations is reported as always under $$0.04\%$$; with only 2 steps, the density formulation is reported to be more accurate and more stable, with the Dupire approach showing a large density oscillation that disappears with 3+ steps.

## 6. ASCII Architecture / Workflow Diagram(s)

### 6.1 Full workflow (transformed density + absorbed boundary masses + pricing)

```
┌──────────────────────────────────────────────────────────────────────────┐
│ Inputs                                                                    │
│  • SABR params: α,β,ρ,ν                                                   │
│  • Forward: f, expiry: τ_ex                                               │
│  • Variant: standard / free-boundary (L(F)=F^β vs |F|^β)                  │
│  • Numerics: n_sd, J (space), N (time), scheme (LS / TR-BDF2 / …)         │
└───────────────┬──────────────────────────────────────────────────────────┘
                │
                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Coordinate & coefficient build                                             │
│  z(F)=∫_f^F du/D(u) ; choose z±=±n_sd √τ_ex                                │
│  Uniform grid z_j ; midpoints z_j-h/2                                      │
│  Map: z → y(z) → F(y) ; compute F̂_j, Ĉ_j, Γ̂_j                            │
│  Time factor: Ê_j(t)=exp(ρν α Γ̂_j t) (updated recursively in time)        │
└───────────────┬──────────────────────────────────────────────────────────┘
                │
                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Semi-discrete conservative operator                                        │
│  θ_T(z_j,t_n) = L_j^n θ(·,t_n)  (tridiagonal flux-difference form)         │
│  Mirror/absorbing endpoint constraints for θ_0, θ_{J+1}                   │
│  Boundary masses: P_L'(t), P_R'(t) from discrete boundary fluxes           │
└───────────────┬──────────────────────────────────────────────────────────┘
                │
                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Time integration to T=τ_ex                                                 │
│  Use L-stable scheme (e.g., Lawson–Swayne or TR-BDF2)                      │
│  Solve tridiagonal systems per stage; update P_L,P_R consistently          │
└───────────────┬──────────────────────────────────────────────────────────┘
                │
                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Pricing                                                                    │
│  For each strike K:                                                        │
│   • locate cell containing z(y(K))                                         │
│   • midpoint quadrature + special near-strike term                         │
│   • add absorbed-mass payoff terms (F_max−K)P_R or (K−F_min)P_L            │
└───────────────┬──────────────────────────────────────────────────────────┘
                │
                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Outputs                                                                    │
│  • Arbitrage-free discrete density θ (typically nonnegative)               │
│  • Boundary absorbed probabilities P_L, P_R                                │
│  • Vanilla prices; implied vols (Black / Bachelier as appropriate)         │
└──────────────────────────────────────────────────────────────────────────┘
```

### 6.2 Lawson–Swayne single-step structure (two IE stages + extrapolation)

```
t_n ──IE(bδ)──► t_{n+b} ──IE(bδ)──► t_{n+2b} ──extrapolate──► t_{n+1}
  θ^n            θ^{n+b}              θ^{n+2b}                  θ^{n+1}
  P^n            P^{n+b}              P^{n+2b}                  P^{n+1}

θ^{n+1} = (√2+1) θ^{n+2b} − √2 θ^{n+b}   (same coefficients for P_L,P_R)
```

## 7. Follow-Up Works & Extensions

### 7.1 Arbitrage-free density construction via stochastic collocation (alternative to PDE time-stepping)

Oosterlee and Grzelak propose an arbitrage-free density implied by the Hagan-formula smile using stochastic collocation: a small set of collocation points on the implied survival function is projected onto an arbitrage-free variable (Gaussian) so that the resulting density is arbitrage free while implied vols stay close to the original input smile; the method is positioned as fast and straightforward to implement via one-dimensional interpolation and a linear-system inversion. [Oosterlee & Grzelak, Journal of Computational Finance 2017] ([risk.net](https://www.risk.net/journal-computational-finance/2463775/arbitrage-arbitrage-free-implied-volatilities))

### 7.2 Model-free arbitrage-free implied-volatility interpolation (extensions of the collocation theme)

Le Floc’h and Oosterlee develop “model-free stochastic collocation” methods for constructing arbitrage-free implied volatilities without relying on a parametric stochastic-vol model: Part I formalizes collocation-based interpolation with arbitrage constraints in a peer-reviewed setting and positions it as a surface-construction tool whose outputs are consistent with an arbitrage-free density. [Le Floc’h & Oosterlee, Decisions in Economics and Finance 2019] ([link.springer.com](https://link.springer.com/article/10.1007/s10203-019-00238-x))

Le Floc’h and Oosterlee extend the collocation approach in Part II by collocating on monotonic splines (including B-spline parameterizations) and emphasizing calibration stability and regularization; the paper explicitly frames the spline-based collocation as a richer arbitrage-free representation capable of fitting complex implied distributions (including multimodal shapes) while respecting moment constraints. [Le Floc’h & Oosterlee, Risks 2019] ([mdpi.com](https://www.mdpi.com/2227-9091/7/1/30))

### 7.3 Negative-rate extensions building on free-boundary SABR ideas

Xiong, Deng, and Wang extend the “free-boundary SABR” idea (using $$|F_t|^\beta$$ in the rate dynamics to cross zero) to a Libor Market Model setting, defining an FB-SABR-LMM intended to recover implied-volatility surfaces and price instruments under negative rates without a shift; the paper positions FB-SABR as the key modeling ingredient for negative-rate support. [Xiong et al., Quantitative Finance and Economics 2020] ([doaj.org](https://doaj.org/article/70f8af9cad514975b2adb40068020259))

### 7.4 Potentially related: asymptotic “arbitrage-free” properties of SABR-type formulas

Fukasawa proves asymptotic arbitrage-free properties (in an asymptotic sense) for several implied-volatility approximations including the SABR formula and rough SABR formula; this is conceptually related to the paper’s motivation (preventing static arbitrage induced by approximations) but does not directly extend the finite-difference density-PDE methodology. [Fukasawa, arXiv 2022] ([arxiv.org](https://arxiv.org/abs/2201.02752))

## 8. Industrial & Real-World Applications

Open-source reference implementation of the paper’s arbitrage-free SABR PDE methodology is available as a Julia package that accompanies the paper and implements the transformed-density finite-difference solver (including Lawson–Swayne time stepping and a free-boundary SABR variant); the repository explicitly states it is not production-ready and is meant to illustrate the techniques. [GitHub: fabienlefloch/ArbitrageFreeSABR.jl] ([github.com](https://github.com/fabienlefloch/ArbitrageFreeSABR.jl))

No verified production deployments (e.g., publicly documented bank/vendor systems explicitly stating use of the Lawson–Swayne or TR-BDF2 arbitrage-free SABR density-PDE approach) were identified at time of writing under the constraint that application claims must be backed by peer-reviewed papers, arXiv papers, or verifiable public GitHub repositories.

## 9. Consolidated Reference List

[1] Cornelis W. Oosterlee, Lech A. Grzelak. “From arbitrage to arbitrage-free implied volatilities.” Journal of Computational Finance, 2017. doi: 10.21314/JCF.2016.316. `https://www.risk.net/journal-computational-finance/2463775/arbitrage-arbitrage-free-implied-volatilities`

[2] Fabien Le Floc’h, Cornelis W. Oosterlee. “Model-free stochastic collocation for an arbitrage-free implied volatility: Part I.” Decisions in Economics and Finance, 2019. doi: 10.1007/s10203-019-00238-x. `https://link.springer.com/article/10.1007/s10203-019-00238-x`

[3] Fabien Le Floc’h, Cornelis W. Oosterlee. “Model-Free Stochastic Collocation for an Arbitrage-Free Implied Volatility, Part II.” Risks, 2019. doi: 10.3390/risks7010030. `https://www.mdpi.com/2227-9091/7/1/30`

[4] Jie Xiong, Geng Deng, Xindong Wang. “Extension of SABR Libor Market Model to handle negative interest rates.” Quantitative Finance and Economics, 2020. doi: 10.3934/qfe.2020007. `https://doaj.org/article/70f8af9cad514975b2adb40068020259`

[5] Masaaki Fukasawa. “On asymptotically arbitrage-free approximations of the implied volatility.” arXiv, 2022. arXiv:2201.02752. `https://arxiv.org/abs/2201.02752`

[6] GitHub repository. “fabienlefloch/ArbitrageFreeSABR.jl.” Julia reference implementation accompanying the paper. `https://github.com/fabienlefloch/ArbitrageFreeSABR.jl`

---
Learn more:
1. [https://www.risk.net/journal-of-computational-finance/2465429/finite-difference-techniques-for-arbitrage-free-sabr](https://www.risk.net/journal-of-computational-finance/2465429/finite-difference-techniques-for-arbitrage-free-sabr)
2. [https://papers.ssrn.com/sol3/papers.cfm?abstract\_id=2402001](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2402001)
3. [https://www.risk.net/journal-computational-finance/2463775/arbitrage-arbitrage-free-implied-volatilities](https://www.risk.net/journal-computational-finance/2463775/arbitrage-arbitrage-free-implied-volatilities)
4. [https://link.springer.com/article/10.1007/s10203-019-00238-x](https://link.springer.com/article/10.1007/s10203-019-00238-x)
5. [https://www.mdpi.com/2227-9091/7/1/30](https://www.mdpi.com/2227-9091/7/1/30)
6. [https://doaj.org/article/70f8af9cad514975b2adb40068020259](https://doaj.org/article/70f8af9cad514975b2adb40068020259)
7. [https://arxiv.org/abs/2201.02752](https://arxiv.org/abs/2201.02752)
8. [https://github.com/fabienlefloch/ArbitrageFreeSABR.jl](https://github.com/fabienlefloch/ArbitrageFreeSABR.jl)

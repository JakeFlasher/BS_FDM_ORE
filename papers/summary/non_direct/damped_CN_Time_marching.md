## 1. Paper Identity

- **Title:** *The Damped Crank-Nicolson Time-Marching Scheme for the Adaptive Solution of the Black-Scholes Equation* ([katalog.ub.uni-heidelberg.de](https://katalog.ub.uni-heidelberg.de/titel/68618822))  
- **Authors:** Christian Goll; Rolf Rannacher; Winnifried Wollner ([katalog.ub.uni-heidelberg.de](https://katalog.ub.uni-heidelberg.de/titel/68618822))  
- **Affiliations (as stated in the paper text provided):**  
  - Institut für Angewandte Mathematik, Universität Heidelberg, Germany (Goll; Rannacher)  
  - Department of Mathematics, University of Hamburg, Germany (Wollner)  
- **Venue / publication metadata:** *The Journal of Computational Finance*, Volume 18, Issue 4, pp. 1–37, published 2015 (library catalog records “30 April 2015”). ([katalog.ub.uni-heidelberg.de](https://katalog.ub.uni-heidelberg.de/titel/68618822))  
- **DOI:** `10.21314/JCF.2015.301` ([katalog.ub.uni-heidelberg.de](https://katalog.ub.uni-heidelberg.de/titel/68618822))  
- **Version note (date discrepancy):** the manuscript text provided is dated “January 3, 2013” (preprint); the peer-reviewed journal version is 2015. ([katalog.ub.uni-heidelberg.de](https://katalog.ub.uni-heidelberg.de/titel/68618822))  

## 2. Problem Statement & Formulation

The paper targets goal-oriented, adaptive space-time discretization for parabolic PDEs with irregular data, using the multivariate Black–Scholes equation (after time reversal) as the model. The numerical objective is an accurate approximation of a linear quantity of interest $$J(u(T))$$ such as a point value $$u(T,x_0)$$ or a sensitivity (“Greek”), under nonsmooth payoffs and nonsmooth dual terminal data that trigger weak damping in standard second-order schemes.

### 2.1 Continuous PDE model (unbounded domain; time-reversed)

Let $$d \ge 1$$, maturity $$T>0$$, strike $$K>0$$, interest rate $$r>0$$, volatilities $$\sigma_i>0$$, and correlation matrix $$\rho=(\rho_{ij})$$ symmetric positive definite with $$\rho_{ii}=1$$ and $$-1 \le \rho_{ij} \le 1$$. Define the symmetric positive definite matrix
$$
\Xi := \bigl(\sigma_i \rho_{ij} \sigma_j\bigr)_{i,j=1}^d.
$$
For $$x \in \mathbb{R}_+^d$$ and $$t \in (0,T]$$, the time-reversed Black–Scholes PDE is
$$
\partial_t u
-\frac{1}{2}\sum_{i,j=1}^d \Xi_{ij} x_i x_j \, \partial_{x_i}\partial_{x_j} u
- r \sum_{i=1}^d x_i \partial_{x_i} u
+ r u
= 0
\quad \text{in } (0,T] \times \mathbb{R}_+^d,
$$
with initial condition (payoff) at $$t=0$$
$$
u(0,x) = u_0(x).
$$
For the basket call/put payoffs used in the text:
$$
u_0(x) =
\begin{cases}
\max\Bigl(\sum_{i=1}^d x_i - K,\, 0\Bigr), & \text{call}, \\
\max\Bigl(K - \sum_{i=1}^d x_i,\, 0\Bigr), & \text{put}.
\end{cases}
$$

### 2.2 Domain truncation to a bounded polygonal domain

A bounded domain is introduced by choosing (in the paper’s notation) upper bounds per coordinate and setting
$$
\Omega := (0,\bar{x}_1)\times \cdots \times (0,\bar{x}_d) \subset \mathbb{R}_+^d,
$$
with Dirichlet boundary subset $$\Gamma_D := \partial \Omega \cap \mathbb{R}_+^d$$ and boundary data compatible with asymptotics:
$$
g(t,x) =
\begin{cases}
\sum_{i=1}^d x_i - K e^{-rt}, & \text{call},\\
0, & \text{put}.
\end{cases}
$$
The bounded-domain problem is
$$
\begin{aligned}
\partial_t u
-\frac{1}{2}\sum_{i,j=1}^d \Xi_{ij} x_i x_j \, \partial_{x_i}\partial_{x_j} u
- r \sum_{i=1}^d x_i \partial_{x_i} u
+ r u
&= 0 && \text{in } (0,T]\times \Omega,\\
u|_{\Gamma_D} &= g && \text{on } [0,T]\times \Gamma_D,\\
u(0,\cdot) &= u_0 && \text{in } \Omega.
\end{aligned}
$$

**Ambiguity flag (notation in Theorem 1 as provided):** the theorem statement in the provided text reuses symbols that appear to denote both the spatial coordinate $$x_i$$ and the truncation boundary; the restatement in Section 4 uses $$\bar{x}_i$$ for the truncation boundary to disambiguate.

### 2.3 Variational setting and weak formulation

Define the weighted Sobolev-like space
$$
V := \Bigl\{ v \in L^2(\Omega)\ \big|\ x_i \partial_{x_i} v \in L^2(\Omega)\ \text{for } i=1,\dots,d \Bigr\},
$$
with scalar product
$$
(u,v)_V := (u,v)_{L^2(\Omega)} + \sum_{i=1}^d (x_i \partial_{x_i} u,\ x_i \partial_{x_i} v)_{L^2(\Omega)}.
$$
Let
$$
V_0 := \overline{C_0^\infty(\Omega)}^{\|\cdot\|_V} = \{ v \in V\ |\ v|_{\Gamma_D}=0\}.
$$
With $$H:=L^2(\Omega)$$ and $$Y \in \{V,V_0\}$$, the Gelfand triple is $$Y \hookrightarrow H \hookrightarrow Y'$$. Define the Bochner space
$$
W(0,T)(Y) := \Bigl\{ w \in L^2(0,T;Y)\ \big|\ \partial_t w \in L^2(0,T;Y') \Bigr\}.
$$

Define coefficient mappings (time-independent coefficients)
$$
A(x) := \Bigl(\frac{1}{2}\Xi_{ij} x_i x_j\Bigr)_{i,j=1}^d,
\qquad
\beta(x) := \Bigl( \bigl(\Xi_{ii} + \frac{1}{2}\sum_{\substack{j=1\\ j\ne i}}^d \Xi_{ij} - r \bigr) x_i \Bigr)_{i=1}^d.
$$
Define the bilinear form for $$v,w \in V$$:
$$
a(v,w) := (A\nabla v,\nabla w)_{L^2(\Omega)} + (\beta \cdot \nabla v,\ w)_{L^2(\Omega)} + r(v,w)_{L^2(\Omega)}.
$$
Using $$((\cdot,\cdot)) := (\cdot,\cdot)_{L^2((0,T)\times \Omega)}$$, the weak form is:

Given $$u_0 \in L^2(\Omega)$$ and an extension $$\tilde{g} \in W(0,T)(V)$$ of the Dirichlet data, find
$$
u \in \tilde{g} + W(0,T)(V_0)
$$
such that
$$
((\partial_t u,\phi)) + \int_0^T a(u(t),\phi(t))\,dt + (u(0),\phi(0))_{L^2(\Omega)} = (u_0,\phi(0))_{L^2(\Omega)}
\quad \forall \phi \in W(0,T)(V_0).
$$

## 3. Core Methodology

**Defined term (first use):** **Dual Weighted Residual (DWR) method** denotes goal-oriented a posteriori error representation/estimation based on primal residuals weighted by dual (adjoint) solution sensitivities.

**Defined term (first use):** **damped Crank–Nicolson** denotes a second-order Crank–Nicolson time stepping in which selected Crank–Nicolson steps are replaced by implicit Euler substeps (“Rannacher damping”) to smooth irregular initial/terminal data.

### 3.1 End-to-end pipeline (space-time FE + damped CN + DWR)

```
┌──────────────────────────────┐
│  Inputs                       │
│  Ω, T, r, σ, ρ, payoff u0      │
│  boundary g(t,x), QoI J(·)     │
│  initial meshes Th, Tk         │
│  damping set JL ⊂ {1,…,M}      │
└───────────────┬───────────────┘
                │
                ▼
┌──────────────────────────────┐
│ Build damped time partition   │
│ Tk → T̂k via JL (split steps) │
│ define J0 (Euler), J1 (CN)    │
└───────────────┬───────────────┘
                │
                ▼
┌──────────────────────────────┐
│ Solve PRIMAL u_kh             │
│ space FE (conforming) +       │
│ time: CN on J1, Euler on J0   │
└───────────────┬───────────────┘
                │ uses same space-time form
                ▼
┌──────────────────────────────┐
│ Solve DUAL z_kh (adjoint)     │
│ backward time marching with   │
│ damping-consistent scheme     │
└───────────────┬───────────────┘
                │
                ▼
┌──────────────────────────────┐
│ Compute DWR estimators        │
│ η ≈ J(u) − J(u_kh) = ηk + ηh  │
│ + local indicators Σk, Σh     │
└───────────────┬───────────────┘
                │
                ▼
┌──────────────────────────────┐
│ Adaptive refinement (Alg. 1)  │
│ balance ηh and ηk via κ       │
│ refine Th and/or Tk           │
└───────────────┬───────────────┘
                │ stop if TOL or caps hit
                ▼
┌──────────────────────────────┐
│ Output                         │
│ u_kh, z_kh, η, refined meshes  │
└──────────────────────────────┘
```

### 3.2 Time discretization: Crank–Nicolson and damping (Section 2.2.1)

Let time points satisfy $$0=t_0 < t_1 < \cdots < t_M=T$$, define intervals $$I_m := (t_{m-1},t_m]$$ and step sizes $$k_m := t_m - t_{m-1}$$.

#### Procedure 1 — Standard Crank–Nicolson time marching (semi-discrete in time)

**Inputs:** initial data $$u_0 \in V$$; boundary extension values $$\tilde{g}_k^m \in V$$ at times $$t_m$$; bilinear form $$a(\cdot,\cdot)$$; time grid $$\{t_m\}_{m=0}^M$$.

**Outputs:** values $$u_k^m \in \tilde{g}_k^m + V_0$$ for $$m=0,\dots,M$$.

1. Initialize $$u_k^0$$ by $$L^2(\Omega)$$ projection:
   $$
   (u_k^0 - u_0,\psi)_{L^2(\Omega)} = 0 \quad \forall \psi \in L^2(\Omega).
   $$
2. For each $$m=1,\dots,M$$, solve for $$u_k^m \in \tilde{g}_k^m + V_0$$:
   $$
   (u_k^m - u_k^{m-1},\psi)_{L^2(\Omega)} + \frac{1}{2}k_m\, a(u_k^m + u_k^{m-1},\psi) = 0
   \quad \forall \psi \in V_0.
   $$
3. Terminate after $$m=M$$; the quantity of interest is evaluated at $$u_k^M \approx u(T)$$.

**Irregular-data limitation stated in the paper:** standard CN requires $$u_0 \in V$$ for well-posedness in this semi-discrete form and fails to achieve optimal temporal convergence for low-regularity initial data; damping is introduced to recover expected accuracy and suppress spurious oscillations.

#### Procedure 2 — Damped Crank–Nicolson (Rannacher-type start-up/generalized damping set)

**Inputs:** initial data $$u_0$$ (allowing lower regularity than $$V$$ in the continuous model); damping index set $$J_L = \{l_1,\dots,l_L\} \subset \{1,\dots,M\}$$; base time grid $$\{t_m\}_{m=0}^M$$.

**Outputs:** time-marching values $$u_k^m$$ on an augmented index set $$J$$ containing half-step indices.

1. For each $$l \in J_L$$, split interval $$I_l$$ by introducing midpoint time
   $$
   t_{l-1/2} := t_l - \frac{1}{2}k_l.
   $$
2. Define index sets (paper’s definition):
   $$
   J_0 := \{l,\ l-1/2\ |\ l \in J_L\},\qquad
   J := J_0 \cup \{1,2,\dots,M\},\qquad
   J_1 := J \setminus J_0.
   $$
3. Define the refined partition $$\hat{\mathcal{T}}_k := \{I_m\ |\ m \in J\}$$ by
   $$
   I_m :=
   \begin{cases}
   (t_{m-1},t_m], & m \in J_1,\\
   (t_{m-1/2},t_m], & m \in J_0,
   \end{cases}
   \qquad
   k_m := |I_m|.
   $$
4. Compute $$u_k^0$$ by
   $$
   (u_k^0,\psi)_{L^2(\Omega)} = (u_0,\psi)_{L^2(\Omega)} \quad \forall \psi \in L^2(\Omega).
   $$
5. For each time index $$m \in J$$, compute $$u_k^m \in \tilde{g}_k^m + V_0$$ by selecting the step type:
   - **Implicit Euler (damped) step** for $$m \in J_0$$:
     $$
     (u_k^m,\psi)_{L^2(\Omega)} + k_m\, a(u_k^m,\psi) = (u_k^{m-1/2},\psi)_{L^2(\Omega)}
     \quad \forall \psi \in V_0.
     $$
   - **Crank–Nicolson step** for $$m \in J_1$$:
     $$
     (u_k^m,\psi)_{L^2(\Omega)} + \frac{1}{2}k_m\, a(u_k^m + u_k^{m-1},\psi) = (u_k^{m-1},\psi)_{L^2(\Omega)}
     \quad \forall \psi \in V_0.
     $$
6. Terminate at the largest index corresponding to $$t=T$$.

**Damping-step count guideline (paper remark):** if initial data belong to $$H^m$$ with $$m \le 2$$, at least $$2-m$$ damping steps are required to recover optimal rates; Dirac-type data in adjoint problems motivate more damping than the primal.

### 3.3 Spatial discretization (Section 2.3)

For each time point $$t_m$$ with $$m \in J \cup \{0\}$$, define a shape-regular mesh $$\mathcal{T}_h^m$$ of $$\Omega$$ into cells $$K$$ (simplices or quads). Define polynomial degree $$s \ge 1$$ and FE spaces
$$
V_h^{s,m} := \{ v \in C^0(\Omega)\ |\ v|_K \in Q_s(K)\ \forall K \in \mathcal{T}_h^m\} \subset V,
\qquad
V_{h,0}^{s,m} := V_h^{s,m} \cap V_0.
$$

#### Procedure 3 — Fully discrete primal damped CN (equations (2.9))

**Inputs:** $$u_0$$, boundary continuation $$\tilde{g}_k^m$$, time partition $$\hat{\mathcal{T}}_k$$ with sets $$J_0,J_1$$, spatial meshes $$\mathcal{T}_h^m$$, FE degree $$s$$.

**Outputs:** $$u_{kh}^m \in \tilde{g}_k^m + V_{h,0}^{s,m}$$ for $$m \in J \cup \{0\}$$.

1. Initialize at $$m=0$$ using the discrete space:
   $$
   (u_{kh}^0,\psi_h)_{L^2(\Omega)} = (u_0 - \tilde{g}_k^0,\psi_h)_{L^2(\Omega)}
   \quad \forall \psi_h \in V_{h,0}^{s,0}.
   $$
2. For each $$m \in J_0$$ (implicit Euler step), solve
   $$
   (u_{kh}^m,\psi_h)_{L^2(\Omega)} + k_m\, a(u_{kh}^m,\psi_h) = (u_{kh}^{m-1/2},\psi_h)_{L^2(\Omega)}
   \quad \forall \psi_h \in V_{h,0}^{s,m}.
   $$
3. For each $$m \in J_1$$ (Crank–Nicolson step), solve
   $$
   (u_{kh}^m,\psi_h)_{L^2(\Omega)} + \frac{1}{2}k_m\, a(u_{kh}^m + u_{kh}^{m-1},\psi_h) = (u_{kh}^{m-1},\psi_h)_{L^2(\Omega)}
   \quad \forall \psi_h \in V_{h,0}^{s,m}.
   $$
4. Terminate at $$t=T$$.

### 3.4 Galerkin-in-time reformulations (Section 3)

The paper makes adjoint consistency explicit by embedding time marching into a (mixed) space-time Galerkin formulation.

#### 3.4.1 cG(r) method for parabolic problems (Section 3.1)

For $$r \ge 1$$ and $$Y \in \{V,V_0\}$$, define trial and test spaces on the original partition $$\mathcal{T}_k=\{I_m\}_{m=1}^M$$:
$$
X_k^r(Y) := \{ \phi_k \in C^0(I;L^2(\Omega))\ |\ \phi_k|_{I_m} \in P_r(I_m;Y)\},
$$
$$
X_{e,k}^{r-1}(Y) := \{ \phi_k : I \to Y\ |\ \phi_k|_{I_m} \in P_{r-1}(I_m;Y),\ \phi_k(0) \in L^2(\Omega)\}.
$$
With boundary interpolation $$g_k$$ and extension $$\tilde{g}_k \in X_k^r(V)$$, the cG(r) semi-discretization reads: given $$u_0 \in V$$, find $$u_k \in \tilde{g}_k + X_k^r(V_0)$$ such that
$$
((\partial_t u_k,\phi_k)) + \int_0^T a(u_k(t),\phi_k(t))\,dt + (u_k(0)-u_0,\phi_{k}^{0,-})_{L^2(\Omega)} = 0
\quad \forall \phi_k \in X_{e,k}^{r-1}(V_0).
$$
For $$r=1$$, this reproduces Crank–Nicolson.

#### 3.4.2 Damped CN as a mixed Galerkin method (Section 3.2)

Damping is modeled as switching from cG(1) on undamped intervals to dG(0) on damped half-steps by allowing controlled discontinuities.

Define time points where discontinuities may occur:
$$
J_{dc} := \{l-1/2\ |\ l \in J_0\},
\qquad
J_c := (J \cup \{0\}) \setminus J_{dc}.
$$
Define the damped trial space for $$r \ge 1$$:
$$
\hat{X}_k^r(Y) :=
\Bigl\{ \phi_k \in L^2(I;Y)\ \Big|\ 
\phi_k|_{I_m} \in P_r(I_m;Y)\ (m \in J_1),\ 
\phi_k|_{I_l} \in P_{r-1}(I_l;Y)\ (l \in J_0),\
[\phi_k]_n = 0\ \forall n \in J_c
\Bigr\}.
$$
With this, the primal damped scheme is posed as: find $$u_k \in \tilde{g}_k + \hat{X}_k^1(V_0)$$ such that
$$
\sum_{m \in J} ((\partial_t u_k,\phi_k))_m
+ \int_0^T a(u_k(t),\phi_k(t))\,dt
+ \sum_{l \in J_{dc}} ([u_k]_l,\phi_k^{l,+})_{L^2(\Omega)}
+ (u_k(0)-u_0,\phi_k^{0,-})_{L^2(\Omega)}
= 0
$$
for all $$\phi_k \in X_{e,k}^{0}(V_0)$$, where $$((\cdot,\cdot))_m$$ denotes the $$L^2(I_m\times \Omega)$$ inner product.

### 3.5 Dual (adjoint) problem and adjoint-consistent damping (Section 3.4)

Let $$J \in V'$$ be a linear functional representing the quantity of interest. The continuous dual problem is: find $$z \in W(0,T)(V_0)$$ such that
$$
-((\phi,\partial_t z)) + \int_0^T a(\phi(t),z(t))\,dt + (\phi(T),z(T))_{L^2(\Omega)} = J(\phi(T))
\quad \forall \phi \in W(0,T)(V_0).
$$

The discrete dual must (i) be usable as a test function in the discrete primal space-time forms and (ii) reflect damping so that “dual of the discrete” matches “discrete of the dual” in the Galerkin sense. The paper derives the dual time marching scheme as the formal dual of the mixed Galerkin-in-time primal scheme.

#### Procedure 4 — Dual time marching scheme (semi-discrete form; equations (3.11))

**Inputs:** same time index sets $$J,J_0,J_1,J_c,J_{dc}$$; step sizes $$k_l$$; bilinear form $$a(\cdot,\cdot)$$; linear functional $$J(\cdot)$$.

**Outputs:** backward-time sequence $$\{z_k^l\}_{l \in J \cup \{0\}} \subset V_0$$.

Let $$\psi \in V_0$$ be an arbitrary test function. Compute unknowns in descending time order using the case distinctions:

1. **Terminal step $$l=M$$**
   - If $$M \in J_0$$:
     $$
     (\psi,z_k^M)_{L^2(\Omega)} + k_M a(\psi,z_k^M) = J(\psi).
     $$
   - If $$M \in J_1$$:
     $$
     (\psi,z_k^M)_{L^2(\Omega)} + \frac{1}{2}k_M a(\psi,z_k^M) = J(\psi).
     $$
2. **Intermediate steps $$0 < l < M$$**
   - If $$l \in J_{dc}$$ and $$l \in J_1$$:
     $$
     (\psi,z_k^l)_{L^2(\Omega)} + \frac{1}{2}k_l a(\psi,z_k^l) = (\psi,z_k^{l+1/2})_{L^2(\Omega)}.
     $$
   - If $$l \in J_{dc}$$ and $$l \in J_0$$:
     $$
     (\psi,z_k^l)_{L^2(\Omega)} + k_l a(\psi,z_k^l) = (\psi,z_k^{l+1/2})_{L^2(\Omega)}.
     $$
   - If $$l \in J_c$$ and $$l \in J_0$$:
     $$
     (\psi,z_k^l)_{L^2(\Omega)} + k_l a(\psi,z_k^l) + \frac{1}{2}k_{l+1} a(\psi,z_k^{l+1})
     = (\psi,z_k^{l+1})_{L^2(\Omega)}.
     $$
   - If $$l \in J_c$$ and $$l \in J_1$$:
     $$
     (\psi,z_k^l)_{L^2(\Omega)} + \frac{1}{2}k_l a(\psi,z_k^l)
     = (\psi,z_k^{l+1})_{L^2(\Omega)} - \frac{1}{2}k_{l+1} a(\psi,z_k^{l+1}).
     $$
3. **Initial step $$l=0$$**
   - If $$0 \in J_{dc}$$:
     $$
     (\psi,z_k^0)_{L^2(\Omega)} = (\psi,z_k^1)_{L^2(\Omega)}.
     $$
   - If $$0 \in J_c$$:
     $$
     (\psi,z_k^0)_{L^2(\Omega)} = (\psi,z_k^1)_{L^2(\Omega)} - \frac{1}{2}k_1 a(\psi,z_k^1).
     $$

The fully discrete dual scheme is identical with $$\psi_h,z_{kh}^l \in V_{h,0}^{s,l}$$.

**Shifted-grid observation (paper remark):** for constant step size and no damping, the dual time marching coincides with CN on time intervals shifted by $$k/2$$.

### 3.6 DWR a posteriori error estimation (Section 4)

The discrete goal error is decomposed as
$$
J(u(T)) - J(u_{kh}(T))
= \bigl(J(u(T)) - J(u_k(T))\bigr) + \bigl(J(u_k(T)) - J(u_{kh}(T))\bigr)
\approx \eta_k + \eta_h,
$$
where $$\eta_k$$ estimates time discretization error and $$\eta_h$$ estimates spatial discretization error.

The paper restricts to continuous linear $$J:V \to \mathbb{R}$$ for exposition, and notes that point evaluations and derivatives (e.g., Delta) are not continuous on $$V$$; this is handled by using local averages or modified function spaces without changing indicator structure.

#### 3.6.1 Global estimator construction via weight reconstruction (Section 4.1)

Define
$$
B(v,w) := ((\partial_t v,w)) + \int_0^T a(v(t),w(t))\,dt + (v(0),w(0))_{L^2(\Omega)},
$$
and the discrete space-time bilinear form
$$
\tilde{B}_e(v_k,w_k)
:=
\sum_{m \in J} ((\partial_t v_k,w_k))_m
+ \int_0^T a(v_k(t),w_k(t))\,dt
+ \sum_{l \in J_{dc}} ([v_k]_l,w_k^{l,+})_{L^2(\Omega)}
+ (v_k(0),w_k(0))_{L^2(\Omega)}.
$$
Primal and dual solutions satisfy (continuous and discrete) variational identities in these forms; the DWR framework yields an error identity (Theorem 2, Section 4).

Practical evaluation introduces (i) substitution of semi-discrete unknowns by fully discrete ones inside residuals, (ii) neglect of a Dirichlet-data oscillation term $$B_k$$, and (iii) reconstruction operators to approximate unknown weights.

Define temporal reconstructions:
$$
\Pi_k^{(u)} := \hat{i}_{2k}^{(2)} - \mathrm{id},
\qquad
\Pi_k^{(z)} := \hat{i}_k^{(1)} - \mathrm{id},
$$
mapping low-order discrete solutions into higher-order-in-time approximations used as weight surrogates. Define spatial reconstructions (patch-wise quadratic enrichment):
$$
\Pi_h^{(u)} := i_{2h}^{(2)} - \mathrm{id},
\qquad
\Pi_h^{(z)} := i_{2h}^{(2)} - \mathrm{id}.
$$
The computable estimators are
$$
\eta_k := \frac{1}{2}\Bigl( \rho(u_{kh})(\Pi_k^{(z)} z_{kh}) + \rho'(z_{kh})(\Pi_k^{(u)} u_{kh}) \Bigr),
\qquad
\eta_h := \frac{1}{2}\Bigl( \rho(u_{kh})(\Pi_h^{(z)} z_{kh}) + \rho'(z_{kh})(\Pi_h^{(u)} u_{kh}) \Bigr),
$$
with residuals defined from $$\tilde{B}_e$$ and $$J(\cdot)$$ as in the paper.

**Special handling at the first damped half step (paper Remark 4.2):** $$\hat{i}_k^{(1)}$$ is modified near $$t_{1/2}$$ when damping is applied at the first primal step, because the dual scheme enforces $$z_{kh}^{1/2}=z_{kh}^0$$, making naive interpolation-error surrogates degenerate.

#### 3.6.2 Localization to time intervals and spatial elements (Section 4.2)

Time localization uses local residuals $$\rho_m(u_{kh},\cdot)$$ and $$\rho_m'(z_{kh},\cdot)$$ on each $$I_m$$ with case distinctions for damped/undamped intervals and for $$m=0$$. The global time estimator splits as
$$
\eta_k = \frac{1}{2}\sum_{m \in J} \hat{\eta}_k^m,
\qquad
\hat{\eta}_k^m := \hat{\eta}_{p,k}^m + \hat{\eta}_{d,k}^m,
$$
and is mapped back to the original unsplit partition $$\mathcal{T}_k$$ by summing the two half-step indicators on each damped base interval.

Space localization avoids oscillatory cell residuals by integration by parts. For each time mesh $$\mathcal{T}_h^l$$ and cell $$K \in \mathcal{T}_h^l$$, the bilinear form is rewritten using the strong operators
$$
A_{BS} v := -\nabla \cdot (A\nabla v) - \beta \cdot \nabla v + r v,
\qquad
A_{BS}' w := -\nabla \cdot (A\nabla w) - \nabla \cdot (\beta w) + r w,
$$
and edge flux jumps $$[n \cdot (A\nabla \cdot)]$$. Elementwise indicators $$\eta_{h,K}$$ are formed by summing over time and then used to drive mesh refinement.

### 3.7 Space-time refinement cycle (Section 4.3, Algorithm 1)

The adaptive strategy balances temporal and spatial contributions by enforcing
$$
\frac{1}{\kappa} \le \frac{\eta_h}{\eta_k} \le \kappa
$$
for a user parameter $$\kappa \ge 1$$ (paper uses $$\kappa=4$$ in experiments).

#### Procedure 5 — Space-time adaptive loop (Algorithm 1, paraphrased with all decision logic)

**Inputs:** caps $$M_{\max},N_{\max}\in \mathbb{N}$$, balance parameter $$\kappa \ge 1$$, tolerance $$\mathrm{TOL}>0$$, initial meshes $$\mathcal{T}_{k}^{0},\mathcal{T}_{h}^{0}$$, initial damped set $$J_L^0$$.

**Outputs:** adapted meshes $$\mathcal{T}_k,\mathcal{T}_h$$, primal/dual solutions, estimator values.

1. Set iteration counter $$n:=0$$.
2. Construct refined time partition $$\hat{\mathcal{T}}_{k^n}$$ from $$\mathcal{T}_{k^n}$$ and damping set $$J_L^n$$ (split each damped base interval into two Euler substeps).
3. Solve primal and dual discrete problems to obtain $$u_{k^n h^n}$$ and $$z_{k^n h^n}$$.
4. Compute local indicator sets $$\Sigma_{k^n}$$ and $$\Sigma_{h^n}$$ and global estimates $$\eta_{k^n}$$ and $$\eta_{h^n}$$.
5. Stop and return if any termination condition holds:
   - $$M_n \ge M_{\max}$$, where $$M_n := \mathrm{card}(\mathcal{T}_{k^n})$$;
   - $$N_n \ge N_{\max}$$, where $$N_n := \mathrm{card}(\mathcal{T}_{h^n})$$;
   - $$|\eta_{k^n}| + |\eta_{h^n}| \le \mathrm{TOL}$$.
6. If $$|\eta_{k^n}| > \kappa |\eta_{h^n}|$$, refine the time mesh $$\mathcal{T}_{k^n}$$ guided by $$\Sigma_{k^n}$$.
7. Else if $$|\eta_{h^n}| > \kappa |\eta_{k^n}|$$, refine the spatial mesh $$\mathcal{T}_{h^n}$$ guided by $$\Sigma_{h^n}$$.
8. Else (balanced), refine both time and space meshes guided by both indicator sets.
9. Update the set of damped intervals $$J_L^{n+1}$$ (paper leaves the detailed selection rule problem-dependent and ties it to data regularity and the current partition).
10. Increment $$n := n+1$$ and repeat.

**Ambiguity flag (selection of $$J_L^{n+1}$$):** Algorithm 1 includes “Compute $$J_L^{n+1}$$” but the provided text does not specify an explicit automatic rule; numerical sections instead parameterize damping by a fixed number of damped intervals at the beginning and/or end (notation $$dk\text{-}m_p\text{-}m_d$$).

## 4. Theoretical Results

### Theorem 1 (domain truncation error; paper Section 2.1)

**Statement (restated in full, with disambiguated boundary notation):**  
Let $$u$$ be a regular solution of the unbounded-domain Black–Scholes problem on $$\mathbb{R}_+^d$$, and let $$u_{co}$$ be a regular solution of the truncated-domain problem on $$\Omega=(0,\bar{x}_1)\times\cdots\times(0,\bar{x}_d)$$ with Dirichlet boundary data $$g$$ on $$\Gamma_D$$. For any point $$(t,x) \in [0,T]\times \Omega$$ such that, for each $$i=1,\dots,d$$,
$$
\ln\Bigl(\frac{\bar{x}_i}{x_i}\Bigr) \ge (\sigma_i^2 + 2r)t \ge 0,
$$
the truncation error satisfies
$$
|u(t,x) - u_{co}(t,x)|
\le
\|u-g\|_{L^\infty((t,T)\times \Gamma_D)}
\cdot
\sum_{i=1}^d
\exp\Biggl(
-\frac{\ln(\bar{x}_i/x_i)\bigl(\ln(\bar{x}_i/x_i) - (\sigma_i^2+2r)t\bigr)^2}{2\sigma_i^2 t}
\Biggr).
$$

**Ambiguity flag:** the provided theorem text contains typographic/OCR artifacts (notably symbol reuse and a variable mismatch $$s$$ vs $$x$$); the restatement above preserves the inequality structure and exponential-decay dependence on distance-to-boundary in log-coordinates but treats the truncation boundary as $$\bar{x}_i$$.

**Proof strategy sketch (as implied by cited source and stated form):**  
The bound has the maximum-principle flavor: the interior error is controlled by the boundary mismatch $$\|u-g\|_{L^\infty}$$ multiplied by exponentially small factors tied to the probability of reaching the artificial boundary under the associated diffusion. The log-distance terms are consistent with transforming the multiplicative diffusion in asset coordinates into additive diffusion in log-coordinates and bounding barrier-hitting contributions. The condition $$\ln(\bar{x}_i/x_i) \ge (\sigma_i^2+2r)t$$ enforces that the point is sufficiently far from the truncation boundary relative to the time horizon.

### Theorem 2 (DWR error representations for time and space discretization; paper Section 4)

**Statement (restated in full):**  
Let $$(u,z) \in \tilde{g}+W(0,T)(V_0)\times W(0,T)(V_0)$$ be stationary points of the continuous primal/dual Lagrangian, let $$(u_k,z_k)\in \tilde{g}_k+\hat{X}_k^1(V_0)\times X_{e,k}^0(V_0)$$ be stationary points of the semi-discrete Lagrangian, and let $$(u_{kh},z_{kh}) \in \tilde{g}_k+\hat{X}_{kh,0}^{1,1}\times X_{e,kh,0}^{0,1}$$ be stationary points of the fully discrete Lagrangian (all in the sense of the paper’s stationarity relations).

Then the temporal discretization error satisfies
$$
J(u(T)) - J(u_k(T))
=
\frac{1}{2}\Bigl(
\rho(u_k)(z - \psi_k) + \rho'(z_k)(u - \varphi_k)
\Bigr)
+ B_k,
$$
and the spatial discretization error satisfies
$$
J(u_k(T)) - J(u_{kh}(T))
=
\frac{1}{2}\Bigl(
\rho(u_{kh})(z_k - \psi_{kh}) + \rho'(z_{kh})(u_k - \varphi_{kh})
\Bigr).
$$
The functions $$(\varphi_k,\psi_k) \in \tilde{g}_k+\hat{X}_k^1(V_0)\times X_{e,k}^0(V_0)$$ and $$(\varphi_{kh},\psi_{kh}) \in \tilde{g}_k+\hat{X}_{kh,0}^{1,1}\times X_{e,kh,0}^{0,1}$$ are arbitrary. The oscillation term $$B_k$$ is
$$
B_k
=
\inf_{\varphi \in \tilde{g}+W(0,T)(V_0)}
\Bigl(
-\frac{1}{2}\tilde{B}_e(\varphi - \tilde{g}_k)(z)
\Bigr),
$$
measuring the effect of approximating time-dependent Dirichlet data.

**Proof strategy sketch (2–4 sentences, paper-consistent):**  
The argument uses the DWR framework: define a primal–dual Lagrangian whose stationary points encode the primal and adjoint equations, then compare continuous and discrete stationary points. Taylor expansion around stationary points and Galerkin orthogonality yield an identity expressing goal error as residuals applied to dual/primal weight errors plus higher-order remainder; choosing arbitrary interpolants $$\varphi_k,\psi_k$$ inserts computable approximation errors. The Dirichlet oscillation term appears because the discrete space enforces an approximation of time-dependent boundary data and breaks exact consistency unless corrected.

**Complexity/guarantee notes explicitly stated in the paper text:**
- Existence/uniqueness of the weak solution follows from continuity and a Gårding inequality for $$a(\cdot,\cdot)$$ on $$V_0$$ (constants $$C,c,\lambda>0$$ with $$a(v,v)\ge c|v|_V^2 - \lambda \|v\|_{L^2(\Omega)}^2$$), giving well-posedness of the continuous problem.
- The refinement-loop stopping criteria and the balance condition $$\frac{1}{\kappa} \le \frac{\eta_h}{\eta_k} \le \kappa$$ provide an algorithmic convergence/termination rule, not a formal mathematical convergence proof.
- The paper’s work-model estimate compares uniform vs adaptive work under an optimal linear solver assumption:
  $$
  \omega_{glob} = N_{glob} M_{glob},
  \qquad
  \omega_{adap} = 6 N_{adap} M_{adap},
  $$
  where the factor $$6$$ accounts for primal solve, dual solve, indicator evaluation, and extra adaptivity iterations.

## 5. Experimental Evaluation

### 5.1 Experimental configuration summary (paper Section 5)

| Item | Value(s) in the paper |
|---|---|
| Model PDE | Time-reversed Black–Scholes in $$d=1$$ and $$d=2$$ on bounded $$\Omega$$ with Dirichlet boundary data |
| Spatial discretization | Conforming FE, degree $$s=1$$ (Q1/P1 implied by examples), with patch-wise quadratic reconstructions for DWR weights |
| Time discretization | Damped CN: implicit Euler on selected half-steps $$J_0$$, CN on $$J_1$$ |
| Damping strategies compared | $$dk\text{-}0\text{-}0$$, $$dk\text{-}1\text{-}0$$, $$dk\text{-}1\text{-}1$$, $$dk\text{-}1\text{-}2$$ (notation: damp first $$m_p$$ and last $$m_d$$ base intervals) |
| Quantities of interest $$J$$ | Option value at $$t=T$$ and $$x=x_0$$; Delta $$\partial_x u(T,x_0)$$ |
| Error metric | True goal error $$|J(e)|$$ with $$e := u - u_{kh}$$ (exact in 1D, reference in 2D); effectivity index $$I_{eff} := \eta / J(e(T))$$ |
| Adaptivity balancing parameter | $$\kappa = 4$$ |
| Implementation | deal.II (C++), UMFPACK direct solver; Gaussian quadrature, order 7 used for projecting 2D put payoff onto FE space |

### 5.2 Test case parameters (paper Table 1)

| Test | $$x_0$$ | $$K$$ | $$T$$ | $$\sigma$$ | Domain upper bound $$\bar{x}$$ | $$r$$ | Reference $$u(T,x_0)$$ | Reference Delta $$\partial_x u(T,x_0)$$ |
|---|---:|---:|---:|---|---|---:|---:|---:|
| 1D Call | $$100$$ | $$100$$ | $$1$$ | $$0.2$$ | $$200$$ | $$\log(1.1)$$ | $$\approx 12.9927$$ | $$\approx 0.7179$$ |
| 2D Put | $$(25,25)$$ | $$25$$ | $$1$$ | $$(1/2,\ 3/10)$$ | $$(100,100)$$ | $$0.05$$ | $$\approx 2.2692$$ (reference mesh) | not reported |

### 5.3 Key quantitative results (paper Tables 2–12)

#### Table 2(a): invariance of $$\eta_h$$ under temporal refinement (1D call, price functional)

| $$M$$ | $$\eta_h$$ (dk-1-0) | $$\eta_h$$ (dk-1-1) |
|---:|---:|---:|
| 4 | $$1.68\cdot 10^{-3}$$ | $$1.88\cdot 10^{-3}$$ |
| 8 | $$1.79\cdot 10^{-3}$$ | $$1.85\cdot 10^{-3}$$ |
| 16 | $$1.87\cdot 10^{-3}$$ | $$1.84\cdot 10^{-3}$$ |
| 32 | $$1.84\cdot 10^{-3}$$ | $$1.84\cdot 10^{-3}$$ |
| 64 | $$1.84\cdot 10^{-3}$$ | $$1.84\cdot 10^{-3}$$ |
| 128 | $$1.84\cdot 10^{-3}$$ | $$1.84\cdot 10^{-3}$$ |

Setup: $$N=129$$.

#### Table 2(b): invariance of $$\eta_k$$ under spatial refinement (1D call, price functional)

| $$N$$ | $$\eta_k$$ (dk-1-0) | $$\eta_k$$ (dk-1-1) |
|---:|---:|---:|
| 33 | $$1.02\cdot 10^{-4}$$ | $$3.20\cdot 10^{-4}$$ |
| 65 | $$1.02\cdot 10^{-4}$$ | $$3.18\cdot 10^{-4}$$ |
| 129 | $$1.01\cdot 10^{-4}$$ | $$3.17\cdot 10^{-4}$$ |
| 257 | $$-2.46\cdot 10^{-4}$$ | $$3.18\cdot 10^{-4}$$ |
| 513 | $$-1.59\cdot 10^{-5}$$ | $$3.14\cdot 10^{-4}$$ |
| 1025 | $$-1.36\cdot 10^{-4}$$ | $$3.17\cdot 10^{-4}$$ |

Setup: $$M=50$$.

**Interpretation stated in the paper:** $$\eta_k$$ is stable under spatial refinement only when damping is applied to the dual (dk-1-1), illustrating the need for dual damping for reliable temporal error estimation.

#### Table 3(a): invariance of $$\eta_h$$ under temporal refinement (2D put, price functional)

| $$M$$ | $$\eta_h$$ (dk-1-0) | $$\eta_h$$ (dk-1-1) |
|---:|---:|---:|
| 8 | $$7.55\cdot 10^{-4}$$ | $$7.74\cdot 10^{-4}$$ |
| 16 | $$7.81\cdot 10^{-4}$$ | $$7.72\cdot 10^{-4}$$ |
| 32 | $$7.72\cdot 10^{-4}$$ | $$7.71\cdot 10^{-4}$$ |
| 64 | $$7.71\cdot 10^{-4}$$ | $$7.71\cdot 10^{-4}$$ |
| 128 | $$7.71\cdot 10^{-4}$$ | $$7.71\cdot 10^{-4}$$ |

Setup: $$N=16641$$.

#### Table 3(b): invariance of $$\eta_k$$ under spatial refinement (2D put, price functional)

| $$N$$ | $$\eta_k$$ (dk-1-0) | $$\eta_k$$ (dk-1-1) |
|---:|---:|---:|
| 81 | $$2.12\cdot 10^{-4}$$ | $$1.09\cdot 10^{-3}$$ |
| 289 | $$2.56\cdot 10^{-4}$$ | $$9.77\cdot 10^{-4}$$ |
| 1089 | $$2.48\cdot 10^{-4}$$ | $$9.44\cdot 10^{-4}$$ |
| 4225 | $$-4.76\cdot 10^{-4}$$ | $$9.43\cdot 10^{-4}$$ |
| 16641 | $$-8.49\cdot 10^{-4}$$ | $$9.38\cdot 10^{-4}$$ |

Setup: $$M=16$$.

#### Table 4: temporal estimator effectivity for dominant temporal error (1D call, price functional)

Setup note in paper: $$\eta_h \approx 1.1\cdot 10^{-7}$$ (spatial error negligible).

| $$M$$ | $$|J(e)|$$ (dk-1-0) | $$\eta_k$$ (dk-1-0) | $$I_{eff}$$ (dk-1-0) | $$|J(e)|$$ (dk-1-1) | $$\eta_k$$ (dk-1-1) | $$I_{eff}$$ (dk-1-1) |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | $$4.17\cdot 10^{-3}$$ | $$-3.67\cdot 10^{-3}$$ | $$-0.88$$ | $$1.35\cdot 10^{-2}$$ | $$1.23\cdot 10^{-2}$$ | $$0.91$$ |
| 16 | $$1.13\cdot 10^{-3}$$ | $$-1.02\cdot 10^{-3}$$ | $$-0.90$$ | $$3.36\cdot 10^{-3}$$ | $$3.07\cdot 10^{-3}$$ | $$0.92$$ |
| 32 | $$2.94\cdot 10^{-4}$$ | $$-2.69\cdot 10^{-4}$$ | $$-0.91$$ | $$8.38\cdot 10^{-4}$$ | $$7.73\cdot 10^{-4}$$ | $$0.92$$ |
| 64 | $$7.51\cdot 10^{-5}$$ | $$-6.94\cdot 10^{-5}$$ | $$-0.92$$ | $$2.09\cdot 10^{-4}$$ | $$1.94\cdot 10^{-4}$$ | $$0.93$$ |
| 128 | $$1.90\cdot 10^{-5}$$ | $$-1.77\cdot 10^{-5}$$ | $$-0.92$$ | $$5.24\cdot 10^{-5}$$ | $$4.86\cdot 10^{-5}$$ | $$0.93$$ |
| 256 | $$4.87\cdot 10^{-6}$$ | $$-4.55\cdot 10^{-6}$$ | $$-0.91$$ | $$1.32\cdot 10^{-5}$$ | $$1.22\cdot 10^{-5}$$ | $$0.93$$ |

**Conclusion stated in paper:** dk-1-0 produces sign errors and poor estimator quality when temporal error dominates; dk-1-1 maintains stable, near-unity effectivity.

#### Table 5: spatial estimator effectivity for dominant spatial error (1D call, dk-1-1)

Setup note in paper: $$\eta_k \approx 3.2\cdot 10^{-6}$$.

| $$N$$ | $$|J(e)|$$ | $$\eta_h$$ | $$I_{eff}$$ |
|---:|---:|---:|---:|
| 9 | $$5.42\cdot 10^{-1}$$ | $$6.79\cdot 10^{-1}$$ | $$1.25$$ |
| 17 | $$1.20\cdot 10^{-1}$$ | $$1.57\cdot 10^{-1}$$ | $$1.30$$ |
| 33 | $$2.95\cdot 10^{-2}$$ | $$3.08\cdot 10^{-2}$$ | $$1.04$$ |
| 65 | $$7.35\cdot 10^{-3}$$ | $$7.42\cdot 10^{-3}$$ | $$1.01$$ |
| 129 | $$1.84\cdot 10^{-3}$$ | $$1.84\cdot 10^{-3}$$ | $$1.00$$ |
| 257 | $$4.62\cdot 10^{-4}$$ | $$4.59\cdot 10^{-4}$$ | $$1.00$$ |
| 513 | $$1.18\cdot 10^{-4}$$ | $$1.15\cdot 10^{-4}$$ | $$1.00$$ |

#### Table 6: spatial estimator effectivity for dominant spatial error (2D put, dk-1-1)

Setup note in paper: $$\eta_k \approx 2.7\cdot 10^{-6}$$.

| $$N$$ | $$|J(e)|$$ | $$\eta_h$$ | $$I_{eff}$$ |
|---:|---:|---:|---:|
| 81 | $$2.62\cdot 10^{-1}$$ | $$1.58\cdot 10^{-1}$$ | $$0.60$$ |
| 289 | $$7.08\cdot 10^{-2}$$ | $$4.83\cdot 10^{-2}$$ | $$0.68$$ |
| 1089 | $$1.80\cdot 10^{-2}$$ | $$1.23\cdot 10^{-2}$$ | $$0.68$$ |
| 4225 | $$4.51\cdot 10^{-3}$$ | $$3.08\cdot 10^{-3}$$ | $$0.68$$ |
| 16641 | $$1.13\cdot 10^{-3}$$ | $$7.71\cdot 10^{-4}$$ | $$0.69$$ |
| 66049 | $$2.82\cdot 10^{-4}$$ | $$1.93\cdot 10^{-4}$$ | $$0.69$$ |

**Paper diagnosis:** effectivity below $$1$$ is attributed to insufficient recovery of the singular dual solution; using only the dual residual (a known DWR practice) improves effectivity (paper reports rising to $$0.98$$).

#### Table 7: adaptive temporal refinement comparison under different damping (1D call, fixed fine space mesh)

Setup: uniform spatial mesh with $$16385$$ cells; spatial error about $$10^{-7}$$.

| Method | $$M$$ | $$|J(e)|$$ | $$I_{eff}$$ |
|---|---:|---:|---:|
| dk-0-0 | 4 | $$4.8\cdot 10^{-1}$$ | $$1230.79$$ |
| dk-0-0 | 8 | $$2.5\cdot 10^{-1}$$ | $$1202.86$$ |
| dk-0-0 | 16 | $$1.2\cdot 10^{-1}$$ | $$1189.18$$ |
| dk-0-0 | 32 | $$6.2\cdot 10^{-2}$$ | $$1179.15$$ |
| dk-0-0 | 64 | $$3.1\cdot 10^{-2}$$ | $$1167.05$$ |
| dk-0-0 | 128 | $$1.6\cdot 10^{-2}$$ | $$1148.03$$ |
| dk-0-0 | 256 | $$7.6\cdot 10^{-3}$$ | $$1117.78$$ |
| dk-0-0 | 512 | $$3.6\cdot 10^{-3}$$ | $$1079.07$$ |
| dk-1-0 | 4 | $$1.4\cdot 10^{-2}$$ | $$-0.89$$ |
| dk-1-0 | 8 | $$4.2\cdot 10^{-3}$$ | $$-0.88$$ |
| dk-1-0 | 10 | $$1.5\cdot 10^{-3}$$ | $$5.33$$ |
| dk-1-0 | 20 | $$2.9\cdot 10^{-4}$$ | $$7.32$$ |
| dk-1-0 | 22 | $$7.4\cdot 10^{-4}$$ | $$4.31$$ |
| dk-1-0 | 24 | $$8.7\cdot 10^{-4}$$ | $$4.17$$ |
| dk-1-0 | 48 | $$2.0\cdot 10^{-4}$$ | $$4.75$$ |
| dk-1-0 | 50 | $$2.2\cdot 10^{-4}$$ | $$6.14$$ |
| dk-1-1 | 4 | $$5.6\cdot 10^{-2}$$ | $$0.93$$ |
| dk-1-1 | 8 | $$1.4\cdot 10^{-2}$$ | $$0.91$$ |
| dk-1-1 | 12 | $$2.0\cdot 10^{-3}$$ | $$0.92$$ |
| dk-1-1 | 20 | $$3.3\cdot 10^{-4}$$ | $$0.91$$ |
| dk-1-1 | 24 | $$4.5\cdot 10^{-4}$$ | $$0.94$$ |
| dk-1-1 | 40 | $$1.4\cdot 10^{-4}$$ | $$0.96$$ |
| dk-1-1 | 72 | $$3.7\cdot 10^{-5}$$ | $$0.98$$ |
| dk-1-1 | 76 | $$4.9\cdot 10^{-5}$$ | $$0.94$$ |

**Paper conclusion:** adaptive temporal refinement is effective only when the dual is damped consistently (dk-1-1).

#### Table 8: adaptive spatial refinement (2D put, dk-1-1, fixed fine time mesh)

Setup: equidistant time mesh $$M=300$$ with $$\eta_k \approx 3\cdot 10^{-6}$$.

| $$N$$ | $$|J(e)|$$ | $$\eta_h$$ | $$I_{eff}$$ |
|---:|---:|---:|---:|
| 81 | $$2.62\cdot 10^{-1}$$ | $$1.58\cdot 10^{-1}$$ | $$0.60$$ |
| 137 | $$7.13\cdot 10^{-2}$$ | $$4.59\cdot 10^{-2}$$ | $$0.64$$ |
| 207 | $$2.22\cdot 10^{-2}$$ | $$1.10\cdot 10^{-2}$$ | $$0.50$$ |
| 523 | $$6.07\cdot 10^{-3}$$ | $$3.99\cdot 10^{-3}$$ | $$0.66$$ |
| 1157 | $$2.06\cdot 10^{-3}$$ | $$1.41\cdot 10^{-3}$$ | $$0.69$$ |
| 2503 | $$7.60\cdot 10^{-4}$$ | $$5.99\cdot 10^{-4}$$ | $$0.79$$ |
| 4991 | $$3.16\cdot 10^{-4}$$ | $$2.63\cdot 10^{-4}$$ | $$0.84$$ |
| 10255 | $$1.41\cdot 10^{-4}$$ | $$1.23\cdot 10^{-4}$$ | $$0.89$$ |

**Paper conclusion:** spatial adaptivity resolves dual singularities better than uniform refinement, improving effectivity.

#### Table 9: simultaneous space-time adaptivity (1D call, dk-1-1, price functional)

| Strategy | $$N$$ | $$M$$ | $$|J(e)|$$ | $$I_{eff}$$ |
|---|---:|---:|---:|---:|
| global uniform | 65 | 32 | $$8.19\cdot 10^{-3}$$ | $$1.00$$ |
| global uniform | 129 | 64 | $$2.04\cdot 10^{-3}$$ | $$0.99$$ |
| global uniform | 257 | 128 | $$5.11\cdot 10^{-4}$$ | $$0.99$$ |
| global uniform | 513 | 256 | $$1.28\cdot 10^{-4}$$ | $$0.99$$ |
| adaptive | 28 | 12 | $$1.12\cdot 10^{-2}$$ | $$0.94$$ |
| adaptive | 48 | 12 | $$4.36\cdot 10^{-3}$$ | $$0.96$$ |
| adaptive | 79 | 20 | $$1.05\cdot 10^{-3}$$ | $$0.98$$ |
| adaptive | 113 | 24 | $$1.21\cdot 10^{-4}$$ | $$0.83$$ |

#### Table 10: simultaneous space-time adaptivity (2D put, dk-1-1, price functional)

| Strategy | $$N$$ | $$M$$ | $$|J(e)|$$ | $$I_{eff}$$ |
|---|---:|---:|---:|---:|
| global uniform | 1089 | 16 | $$1.90\cdot 10^{-2}$$ | $$0.70$$ |
| global uniform | 4225 | 32 | $$4.77\cdot 10^{-3}$$ | $$0.70$$ |
| global uniform | 16641 | 64 | $$1.19\cdot 10^{-3}$$ | $$0.70$$ |
| global uniform | 66049 | 128 | $$2.95\cdot 10^{-4}$$ | $$0.70$$ |
| global uniform | 263169 | 256 | $$7.15\cdot 10^{-5}$$ | $$0.73$$ |
| adaptive | 207 | 8 | $$2.64\cdot 10^{-2}$$ | $$0.56$$ |
| adaptive | 523 | 12 | $$6.57\cdot 10^{-3}$$ | $$0.68$$ |
| adaptive | 1157 | 12 | $$2.56\cdot 10^{-3}$$ | $$0.73$$ |
| adaptive | 2515 | 20 | $$8.14\cdot 10^{-4}$$ | $$0.80$$ |
| adaptive | 4971 | 20 | $$3.71\cdot 10^{-4}$$ | $$0.84$$ |

#### Table 11: Greek Delta estimation (1D call, dk-1-2, goal $$J(u)=\partial_x u(T,x_0)$$)

| Strategy | $$N$$ | $$M$$ | $$|J(e)|$$ | $$I_{eff}$$ |
|---|---:|---:|---:|---:|
| global uniform | 9 | 4 | $$4.66\cdot 10^{-2}$$ | $$0.52$$ |
| global uniform | 17 | 8 | $$1.36\cdot 10^{-2}$$ | $$1.03$$ |
| global uniform | 33 | 16 | $$3.68\cdot 10^{-3}$$ | $$1.20$$ |
| global uniform | 65 | 32 | $$9.41\cdot 10^{-4}$$ | $$1.26$$ |
| global uniform | 129 | 64 | $$2.37\cdot 10^{-4}$$ | $$1.27$$ |
| global uniform | 257 | 128 | $$5.93\cdot 10^{-5}$$ | $$1.28$$ |
| global uniform | 513 | 256 | $$1.48\cdot 10^{-5}$$ | $$1.28$$ |
| global uniform | 1025 | 512 | $$3.71\cdot 10^{-6}$$ | $$1.28$$ |
| global uniform | 2049 | 1024 | $$9.26\cdot 10^{-7}$$ | $$1.28$$ |
| adaptive | 14 | 4 | $$4.47\cdot 10^{-3}$$ | $$1.40$$ |
| adaptive | 21 | 8 | $$1.59\cdot 10^{-3}$$ | $$1.51$$ |
| adaptive | 28 | 12 | $$4.84\cdot 10^{-6}$$ | $$34.05$$ |
| adaptive | 33 | 20 | $$1.04\cdot 10^{-4}$$ | $$-0.60$$ |
| adaptive | 39 | 24 | $$1.62\cdot 10^{-4}$$ | $$0.50$$ |
| adaptive | 50 | 40 | $$9.50\cdot 10^{-5}$$ | $$0.84$$ |
| adaptive | 60 | 40 | $$7.01\cdot 10^{-5}$$ | $$1.28$$ |
| adaptive | 68 | 40 | $$4.62\cdot 10^{-5}$$ | $$1.10$$ |
| adaptive | 97 | 40 | $$1.34\cdot 10^{-5}$$ | $$0.98$$ |

**Paper note:** higher damping near $$t=T$$ is used for Delta due to higher singularity of the dual terminal data.

#### Table 12: work model comparison (paper’s complexity discussion)

Work model (paper): $$\omega_{glob} = N_{glob}M_{glob}$$ and $$\omega_{adap} = 6N_{adap}M_{adap}$$.

| Example | global $$N$$ | global $$M$$ | global $$|J(e)|$$ | adaptive $$N$$ | adaptive $$M$$ | adaptive $$|J(e)|$$ | $$\omega_{glob}/\omega_{adap}$$ |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1D call (price) | 513 | 256 | $$1.28\cdot 10^{-4}$$ | 113 | 24 | $$1.21\cdot 10^{-4}$$ | $$\approx 8$$ |
| 2D put (price) | 66049 | 128 | $$2.95\cdot 10^{-4}$$ | 4971 | 20 | $$3.71\cdot 10^{-4}$$ | $$\approx 28$$ |
| 1D call (Delta) | 513 | 256 | $$1.48\cdot 10^{-5}$$ | 97 | 40 | $$1.34\cdot 10^{-5}$$ | $$\approx 6$$ |

## 6. ASCII Architecture / Workflow Diagram(s)

### 6.1 Space-time DWR workflow with damping-aware dual

```
┌──────────────────────────────────────────────────────────────────────────────┐
│  Continuous problem (Section 2.1)                                            │
│  Black–Scholes on Ω, time-reversed, weak form in W(0,T)(V0)                  │
└───────────────────────────────┬──────────────────────────────────────────────┘
                                │ choose QoI J(·)
                                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│  Discretize time (Section 2.2)                                               │
│  Tk = {Im}  → choose JL (damped base steps)                                  │
│  build T̂k with J0 (Euler half-steps) and J1 (CN steps)                      │
└───────────────────────────────┬──────────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│  Discretize space (Section 2.3, 3.3)                                         │
│  meshes Thm, FE spaces Vh^{s,m} ⊂ V, impose ΓD via continuation g̃k           │
└───────────────────────────────┬──────────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│  Solve primal u_kh (Section 2.3)                                             │
│  forward in time on T̂k:                                                     │
│   - m ∈ J0: implicit Euler                                                   │
│   - m ∈ J1: Crank–Nicolson                                                   │
└───────────────────────────────┬──────────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│  Solve dual z_kh (Section 3.4)                                               │
│  backward scheme derived as formal dual of the mixed Galerkin-in-time form   │
│  with damping-consistent case distinctions (3.11)                            │
└───────────────────────────────┬──────────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│  DWR estimator (Section 4)                                                   │
│  η ≈ (1/2)[ρ(u_kh)(Π z_kh) + ρ'(z_kh)(Π u_kh)]                               │
│  split into ηk (time) + ηh (space) and localize to Σk, Σh                    │
└───────────────────────────────┬──────────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│  Adaptivity loop (Algorithm 1)                                               │
│  refine Tk and/or Th based on κ-balance, stop on TOL or caps                 │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 6.2 Dual vs primal damping locations (dk-$$m_p$$-$$m_d$$ notation)

```
Time axis (primal forward):   t=0 (payoff u0) ───────────────→ t=T (QoI at “today”)

Base steps:                     1      2      …      M-md   …    M-1    M
Damping for primal u:          [Euler] [Euler] …              (CN)  (CN)  (CN)
Damping for dual z:             (CN)   (CN)   …              [Euler][Euler]… (near t=T)
                               ^ dual is marched backward, so terminal data J acts here
```

## 7. Follow-Up Works & Extensions

### 7.1 Space-time DWR extensions with open-source implementations

**Partition-of-unity localized space-time DWR estimators with reproducible code** [Thiele & Wick, J Sci Comput 2024].  
Thiele and Wick develop space-time goal-oriented a posteriori estimators for parabolic problems with partition-of-unity (PU) localization, compare “split” vs “joint” estimators, and provide open-source implementations and reproducibility artifacts. The paper explicitly cites Goll–Rannacher–Wollner’s damped CN/DWR work as a relevant space-time DWR reference and positions its own contributions as an extension to broader parabolic settings and software realizations. ([link.springer.com](https://link.springer.com/article/10.1007/s10915-024-02485-6))

**Variational discretization with measure controls citing the damped CN/DWR Black–Scholes work** [Herberg et al., MCRF 2020].  
Herberg, Hinze, and Schumacher study a parabolic optimal control problem with space-time measure controls and propose a variational discretization in which discrete controls become Dirac measures located at space-time grid points. The work cites the damped Crank–Nicolson/DWR Black–Scholes paper among its references and provides convergence results (strong in $$L^q$$ for states, weak-* in measures for controls) for their discrete solutions. ([aimsciences.org](https://www.aimsciences.org/article/doi/10.3934/mcrf.2020018))

### 7.2 Alternative strategies for nonsmooth initial/terminal data in Crank–Nicolson-type schemes (finance-relevant)

**Time-variable transformation to restore convergence for Dirac-type initial data without Rannacher steps** [Reisinger & Whitley, arXiv 2012].  
Reisinger and Whitley analyze Crank–Nicolson for the heat equation with Dirac delta initial data, prove divergence in original variables under mesh-ratio constraints, and show that a square-root time change yields convergence with the mesh ratio controlling the order; their experiments indicate quadratic convergence for option price/Greeks without Rannacher start-up in some regimes. The work is a methodologically independent response to the same irregular-data pathology that motivates damping in the source paper. ([arxiv.org](https://arxiv.org/abs/1210.5487))

**Smoothing Crank–Nicolson for barrier options with repeated time-discontinuities** [Wade et al., JCAM 2007].  
Wade et al. address oscillations from nonsmooth payoffs and barrier-triggered time discontinuities, describe how standard Luskin–Rannacher smoothing must be modified when discontinuities occur repeatedly in time, and note that the number of smoothing steps can exceed two in practice. The paper is independent of the DWR/adjoint-consistency focus but reinforces the empirical “damping amount depends on irregularity” message used in the source paper for both primal and dual problems. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042706002469))

### 7.3 Finite element option-pricing models using Rannacher-type damping (related discretization pattern)

**Finite elements for transaction-cost option pricing with Rannacher time stepping** [Wei et al., arXiv 2020].  
Wei, Erlangga, and Zhumakhanova apply finite elements to Leland’s option pricing model (transaction costs) and combine spatial $$P_1/P_2$$ elements with a Crank–Nicolson-type time discretization implemented via a Rannacher approach, reporting favorable agreement with finite-difference results and emphasizing stability control via space–time mesh-size ratios. The setting is different (nonlinear model modification) but uses the same FE + Rannacher-damping pattern for nonsmooth payoffs. ([arxiv.org](https://arxiv.org/abs/2010.13541))

## 8. Industrial & Real-World Applications

**Open-source space-time DWR implementations for parabolic PDEs** [GitHub: jpthiele/pu-dwr-diffusion; GitHub: jpthiele/pu-dwr-combustion].  
The J. Sci. Comput. 2024 work provides open-source code bases implementing PU-localized space-time DWR estimators and adaptive algorithms for diffusion and nonlinear combustion benchmarks, with accompanying reproducibility artifacts on Zenodo. The practical role aligns with the source paper’s emphasis on adjoint-based error estimation and adaptive mesh/time-step control for time-dependent PDEs with potentially singular adjoint data. ([link.springer.com](https://link.springer.com/article/10.1007/s10915-024-02485-6))

**Adaptive space-time PDE simulation in computational biology using damped Crank–Nicolson and DWR** [Krämer et al., PLOS Comput Biol 2015].  
A PLOS Computational Biology study reports a problem-specific simulation code (built on deal.II) for 3D cytokine signaling between T cells, discretized in time by a damped Crank–Nicolson method and controlled by adaptive space and time grids via the DWR method. This constitutes a non-finance deployment of the same numerical design pattern (damped CN + DWR-driven adaptivity) motivated in the source paper for irregular-data parabolic PDEs. ([journals.plos.org](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004206))

**Open-source quantitative-finance software framework (contextual PDE pricing ecosystem)** [GitHub: lballabio/QuantLib].  
QuantLib is a large open-source C++ library aimed at providing a comprehensive framework for quantitative finance (modeling, trading, risk management), commonly used as an engineering substrate in which PDE-based pricing components (including Crank–Nicolson-type evolvers in some language bindings/ports) can be implemented. The source paper’s techniques target a specialized FEM+DWR adaptive pipeline rather than a general-purpose library architecture; no direct verified incorporation of the paper’s adjoint-consistent dual damping is identified here. ([github.com](https://github.com/lballabio/QuantLib))

## 9. Consolidated Reference List

[1] J. P. Thiele, T. Wick. “Numerical Modeling and Open-Source Implementation of Variational Partition-of-Unity Localizations of Space-Time Dual-Weighted Residual Estimators for Parabolic Problems.” *Journal of Scientific Computing*, 2024. DOI: `10.1007/s10915-024-02485-6`. URL: `https://doi.org/10.1007/s10915-024-02485-6`. ([link.springer.com](https://link.springer.com/article/10.1007/s10915-024-02485-6))

[2] GitHub: jpthiele/pu-dwr-diffusion. “PU-DWR diffusion code (space-time goal-oriented adaptivity).” URL: `https://github.com/jpthiele/pu-dwr-diffusion`. ([link.springer.com](https://link.springer.com/article/10.1007/s10915-024-02485-6))

[3] GitHub: jpthiele/pu-dwr-combustion. “PU-DWR combustion code (space-time goal-oriented adaptivity).” URL: `https://github.com/jpthiele/pu-dwr-combustion`. ([link.springer.com](https://link.springer.com/article/10.1007/s10915-024-02485-6))

[4] E. Herberg, M. Hinze, H. Schumacher. “Maximal discrete sparsity in parabolic optimal control with measures.” *Mathematical Control and Related Fields*, 10(4), 735–759, 2020. DOI: `10.3934/mcrf.2020018`. URL: `https://doi.org/10.3934/mcrf.2020018`. ([aimsciences.org](https://www.aimsciences.org/article/doi/10.3934/mcrf.2020018))

[5] C. Reisinger, A. Whitley. “The impact of a natural time change on the convergence of the Crank-Nicolson scheme.” arXiv, 2012. arXiv:1210.5487. URL: `https://arxiv.org/abs/1210.5487`. ([arxiv.org](https://arxiv.org/abs/1210.5487))

[6] B. Wade, A. Khaliq, M. Yousuf, J. Vigo-Aguiar, R. Deininger. “On smoothing of the Crank–Nicolson scheme and higher order schemes for pricing barrier options.” *Journal of Computational and Applied Mathematics*, 204(1), 144–158, 2007. DOI: `10.1016/j.cam.2006.04.034`. URL: `https://doi.org/10.1016/j.cam.2006.04.034`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042706002469))

[7] D. Wei, Y. A. Erlangga, G. Zhumakhanova. “A Finite Element Approach to the Numerical Solutions of Leland’s Model.” arXiv, 2020. arXiv:2010.13541. URL: `https://arxiv.org/abs/2010.13541`. ([arxiv.org](https://arxiv.org/abs/2010.13541))

[8] (Authors as listed on the article page) “Three-Dimensional Gradients of Cytokine Signaling between T Cells.” *PLOS Computational Biology*, 2015. DOI: `10.1371/journal.pcbi.1004206`. URL: `https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004206`. ([journals.plos.org](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004206))

[9] GitHub: lballabio/QuantLib. “QuantLib: the free/open-source library for quantitative finance.” URL: `https://github.com/lballabio/QuantLib`. ([github.com](https://github.com/lballabio/QuantLib))

[10] quantlib.js documentation. “Class CrankNicolson (finite difference evolver).” URL: `https://quantlib.js.org/docs/classes/_ql_methods_finitedifferences_cranknicolson_.cranknicolson.html`. ([quantlib.js.org](https://quantlib.js.org/docs/classes/_ql_methods_finitedifferences_cranknicolson_.cranknicolson.html))

---
Learn more:
1. [https://katalog.ub.uni-heidelberg.de/titel/68618822](https://katalog.ub.uni-heidelberg.de/titel/68618822)
2. [https://link.springer.com/article/10.1007/s10915-024-02485-6](https://link.springer.com/article/10.1007/s10915-024-02485-6)
3. [https://www.aimsciences.org/article/doi/10.3934/mcrf.2020018](https://www.aimsciences.org/article/doi/10.3934/mcrf.2020018)
4. [https://arxiv.org/abs/1210.5487](https://arxiv.org/abs/1210.5487)
5. [https://www.sciencedirect.com/science/article/pii/S0377042706002469](https://www.sciencedirect.com/science/article/pii/S0377042706002469)
6. [https://arxiv.org/abs/2010.13541](https://arxiv.org/abs/2010.13541)
7. [https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004206](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004206)
8. [https://github.com/lballabio/QuantLib](https://github.com/lballabio/QuantLib)
9. [https://quantlib.js.org/docs/classes/\_ql\_methods\_finitedifferences\_cranknicolson\_.cranknicolson.html](https://quantlib.js.org/docs/classes/_ql_methods_finitedifferences_cranknicolson_.cranknicolson.html)

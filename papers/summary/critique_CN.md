## 1. Paper Identity

- **Title:** *A Critique of the Crank Nicolson Scheme: Strengths and Weaknesses for Financial Instrument Pricing*. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
- **Author(s):** Daniel J. Duffy. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
- **Affiliation(s):** Datasim Component Technology BV (Datasim; Amsterdam-based trainer and software developer, per author bio). ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
- **Venue:** *Wilmott* (Wilmott magazine; bibliographic record lists “Wilmott, Vol. 2004, Issue 4”). ([osti.gov](https://www.osti.gov/biblio/2502172))  
- **Year / date:** July 2004. ([osti.gov](https://www.osti.gov/biblio/2502172))  
- **Pages:** 68–76. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
- **DOI:** `10.1002/wilm.42820040417`. ([osti.gov](https://www.osti.gov/biblio/2502172))  
- **arXiv ID:** None reported in the paper/metadata available in the accessed sources. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

## 2. Problem Statement & Formulation

The paper analyzes finite-difference discretizations for option-pricing PDEs, focusing on why the widely used Crank–Nicolson (CN) time discretization can generate non-physical/spurious oscillations in convection-diffusion regimes and around non-smooth payoff/condition data, and it proposes more robust alternatives (notably exponentially fitted schemes and the Keller box scheme) that preserve stability/monotonicity properties more uniformly, including in small-volatility (convection-dominated) limits. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

The canonical one-factor Black–Scholes PDE is stated as: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
-\frac{\partial V}{\partial t}
+\frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2}
+r S \frac{\partial V}{\partial S}
-rV
= 0.
$$

For generality, the paper recasts the discussion in terms of a generic parabolic initial–boundary value problem on a finite spatial interval $$\,(A,B)\,$$ and time horizon $$\,(0,T)\,$$: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\begin{aligned}
Lu &\equiv -\frac{\partial u}{\partial t}
+ \sigma(x,t)\frac{\partial^2 u}{\partial x^2}
+ \mu(x,t)\frac{\partial u}{\partial x}
+ b(x,t)u
= f(x,t), \quad (x,t)\in (A,B)\times(0,T),\\
u(x,0) &= \varphi(x), \quad x\in (A,B),\\
u(A,t) &= g_0(t),\quad u(B,t)=g_1(t), \quad t\in (0,T).
\end{aligned}
$$

The technical objective underlying the critique is not just stability in the von Neumann sense, but avoiding oscillations and preserving maximum-principle/monotonicity behavior when $$\sigma$$ is small (singular perturbation / boundary-layer behavior) or when data are only piecewise smooth (typical payoffs and certain boundary conditions). ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

## 3. Core Methodology

### 3.1 Discretization primitives and notation

Uniform tensor-product grids are used repeatedly. For parabolic problems, the paper writes a uniform mesh on $$[A,B]\times[0,T]$$ as: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
A=x_0 < x_1 < \cdots < x_J = B,\quad h=x_j-x_{j-1}\ \text{constant},
$$
$$
0=t_0 < t_1 < \cdots < t_N = T,\quad k=\frac{T}{N}\ \text{constant}.
$$

CN-style averaging in time is defined (for the one-factor Black–Scholes discretization) by: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
V_j^{n+\frac{1}{2}} \equiv \frac{1}{2}\left(V_j^{n+1}+V_j^{n}\right).
$$

The paper uses operator notation $$D_+D_-$$ and $$D_0$$ in its fitted-scheme presentation. On a uniform grid, these correspond to the standard centered second derivative and centered first derivative:
$$
D_+D_- U_j \equiv \frac{U_{j+1}-2U_j+U_{j-1}}{h^2},\quad
D_0 U_j \equiv \frac{U_{j+1}-U_{j-1}}{2h}.
$$
This matches the explicit centered-difference formulas written earlier for $$u_{xx}$$ and $$u_x$$. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

### 3.2 Fully implicit time stepping (implicit Euler in time; centered in space)

**Named method:** fully implicit finite-difference scheme for the one-factor Black–Scholes equation.

**Discrete equation (as presented):** at mesh index $$j$$ and time level $$n\to n+1$$, the scheme replaces the time derivative by a one-sided difference and uses three-point centered differences in $$S$$, yielding (paper notation): ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
-\frac{V_j^{n+1}-V_j^{n}}{k}
+ r\,j\Delta S\left(\frac{V_{j+1}^{n+1}-V_{j-1}^{n+1}}{2\Delta S}\right)
+\frac{1}{2}\sigma^2 j^2 \Delta S^2
\left(\frac{V_{j+1}^{n+1}-2V_j^{n+1}+V_{j-1}^{n+1}}{(\Delta S)^2}\right)
= rV_j^{n+1}.
$$

**Tridiagonal linear system form:** the paper rewrites the fully implicit scheme into: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
a_j^{n+1}V_{j-1}^{n+1}+b_j^{n+1}V_j^{n+1}+c_j^{n+1}V_{j+1}^{n+1}=F_j^{n+1},
$$
with coefficients stated as: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
a_j^{n+1}=\left(\frac{1}{2}\sigma^2 j^2 k-\frac{k r j}{2}\right),\quad
c_j^{n+1}=\left(\frac{1}{2}\sigma^2 j^2 k+\frac{k r j}{2}\right),\quad
F_j^{n+1}=-V_j^{n}.
$$
The PDF extraction shows a line-break ambiguity in $$b_j^{n+1}$$; it is presented as a negative quantity of the form $$-(1+\sigma^2 j^2 k + r\cdot(\cdot))$$, and the intended dimensional form is consistent with a term proportional to $$rk$$. Ambiguity is flagged rather than reconstructed beyond what is visible. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

**Algorithm (step-by-step procedure).**

1. **Inputs:** parameters $$r,\sigma$$; spatial grid $$\{S_j\}$$ with step $$\Delta S$$; time step $$k$$; boundary conditions and a terminal/initial condition depending on time orientation of the pricing PDE.
2. **Initialize** $$V_j^{0}$$ from the payoff (or appropriate initial condition for the chosen time direction) and impose boundary values at $$j=0$$ and $$j=J$$ for all time levels.
3. **For** $$n=0,1,\ldots,N-1$$:  
   3.1. Form the tridiagonal coefficients $$a_j^{n+1}, b_j^{n+1}, c_j^{n+1}$$ and right-hand side $$F_j^{n+1}$$ for interior indices $$j=1,\ldots,J-1$$. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
   3.2. Solve the tridiagonal linear system for $$\{V_j^{n+1}\}_{j=1}^{J-1}$$ using a direct tridiagonal solver (the paper mentions LU decomposition for the tridiagonal system). ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
   3.3. Enforce boundary conditions at time level $$n+1$$.
4. **Outputs:** grid function $$V_j^{n}$$ at all time levels, particularly the desired price at the valuation time and spatial location; Greeks via post-processing finite differences if needed.

**Properties claimed:** unconditional stability (no restriction coupling $$k$$ and $$\Delta S$$), and “no spurious oscillations” for this fully implicit baseline; drawback is first-order accuracy in $$k$$, mitigated via extrapolation to recover second order. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

### 3.3 Crank–Nicolson time discretization (“what is CN, really?”) and failure modes

**Named method:** Crank–Nicolson (CN), framed as “centred differencing in space combined with averaging in time” for convection–diffusion-type equations. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

**Definition of half-step state:** ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
V_j^{n+\frac{1}{2}} = \frac{1}{2}\left(V_j^{n+1}+V_j^{n}\right).
$$

**CN scheme for the one-factor Black–Scholes discretization (as presented):** ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\begin{aligned}
-\frac{V_j^{n+1}-V_j^{n}}{k}
&+ r\,j\Delta S\left(\frac{V_{j+1}^{n+\frac{1}{2}}-V_{j-1}^{n+\frac{1}{2}}}{2\Delta S}\right)
+\frac{1}{2}\sigma^2 j^2 \Delta S^2
\left(\frac{V_{j+1}^{n+\frac{1}{2}}-2V_j^{n+\frac{1}{2}}+V_{j-1}^{n+\frac{1}{2}}}{(\Delta S)^2}\right)
\\
&= \tau V_j^{n+\frac{1}{2}}.
\end{aligned}
$$
The symbol on the right-hand side appears as $$\tau$$ in the extracted PDF text; context suggests it is the coefficient of the reaction term (for Black–Scholes, that coefficient is $$r$$). The output treats this as an OCR/encoding ambiguity rather than redefining the model. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

**CN scheme for the generic parabolic problem (6):** the paper applies centered differences in space at $$x_j$$ and CN-style averaging in time, yielding: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\begin{aligned}
-\frac{u_j^{n+1}-u_j^n}{k}
&+\sigma_j^{n+\frac{1}{2}}\frac{u_{j+1}^{n+\frac{1}{2}}-2u_j^{n+\frac{1}{2}}+u_{j-1}^{n+\frac{1}{2}}}{h^2}
+\mu_j^{n+\frac{1}{2}}\frac{u_{j+1}^{n+\frac{1}{2}}-u_{j-1}^{n+\frac{1}{2}}}{2h}
\\
&+ b_j^{n+\frac{1}{2}}u_j^{n+\frac{1}{2}}
= f_j^{n+\frac{1}{2}}.
\end{aligned}
$$
It is then rewritten into a tridiagonal linear system at each time step. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

**Algorithm (step-by-step procedure) for CN time stepping (generic form).**

1. **Inputs:** coefficients $$\sigma(x,t),\mu(x,t),b(x,t),f(x,t)$$; grid sizes $$h,k$$; boundary values $$g_0(t_n),g_1(t_n)$$; initial condition $$\varphi(x_j)$$.
2. **Initialize:** set $$u_j^0=\varphi(x_j)$$ for $$j=1,\ldots,J-1$$; set boundary nodes $$u_0^n=g_0(t_n)$$ and $$u_J^n=g_1(t_n)$$ for all $$n$$.
3. **For** $$n=0,1,\ldots,N-1$$:  
   3.1. Form the tridiagonal system implied by the CN discretization at the half-time level $$n+\tfrac{1}{2}$$ (the paper denotes the resulting standard form as equation (8)). ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
   3.2. Solve the tridiagonal system for $$u^{n+1}$$.
4. **Outputs:** discrete solution $$u_j^n$$ at mesh points; derived quantities (option sensitivities) by finite differences.

**Failure/weakness mechanisms catalogued (Section 3):**

- **Convection-dominated regime:** the Black–Scholes PDE is treated as a convection–diffusion equation; centered differences can oscillate when diffusion is small or convection is large unless the mesh is sufficiently fine (a “critical” spatial mesh size is referenced). ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
- **Loss of monotonicity in the linear system:** for tridiagonal discretizations, the paper stresses that if the coefficient (denoted $$a_j^n$$ in the standard form) is not positive, the computed solution can oscillate or become non-physical; this is highlighted as a key diagnostic for breakdown. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
- **Time-varying/decaying volatility:** a concrete problematic pattern is given as $$\sigma(t)=\sigma_0 e^{-\alpha(T-t)}$$, where diffusion becomes small as maturity approaches, exacerbating convection dominance. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
- **Degenerate limit $$\sigma\to 0$$:** formally setting volatility to zero turns the PDE into a first-order hyperbolic equation; the paper describes the resulting centered-difference CN-type discretization as “weakly stable,” with poor dissipation of initial errors and sensitivity to rounding errors, motivating one-sided (upwind) discretizations. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
- **Boundary-condition discretization effects:** Neumann/Robin conditions require boundary-derivative approximations; first-order boundary discretization or nonuniform meshes can spoil the expected second-order behavior of CN in the interior. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
- **Non-smooth payoff (initial data) and compatibility:** typical payoff functions are not smooth (the call payoff has a kink at the strike), causing early-time oscillations and degraded delta/gamma accuracy; incompatibility at corners $$\varphi(A)=g_0(0)$$ and $$\varphi(B)=g_1(0)$$ is explicitly singled out as a source of local inaccuracy. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
- **Stability-proof mismatch:** von Neumann analysis is criticized as formally applicable only to constant-coefficient linear initial-value problems, whereas Black–Scholes settings often violate these conditions; maximum-principle/matrix-theoretic approaches are advocated. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

### 3.4 Exponentially fitted finite-difference schemes

#### 3.4.1 Fitting-factor derivation (steady 1D convection–diffusion)

The “fitting factor” is introduced by forcing the discrete scheme to reproduce the exact solution at mesh points for the constant-coefficient boundary-value problem: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\sigma \frac{d^2u}{dx^2} + \mu \frac{du}{dx}=0,\quad x\in (A,B),\qquad
u(A)=\beta_0,\quad u(B)=\beta_1,
$$
with $$\sigma>0$$ and $$\mu>0$$.

A fitted discrete scheme is written as: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\sigma \rho\, D_+D_- U_j + \mu\, D_0 U_j = 0,\quad j=1,\ldots,J-1;\qquad
U_0=\beta_0,\ U_J=\beta_1.
$$

The fitting factor is chosen so that solutions coincide at the mesh points, yielding: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\rho = \frac{\mu h}{2\sigma}\coth\left(\frac{\mu h}{2\sigma}\right),
\quad
\coth x = \frac{e^x + e^{-x}}{e^x - e^{-x}} = \frac{e^{2x}+1}{e^{2x}-1}.
$$

#### 3.4.2 Variable-coefficient fitted scheme for a stationary problem

The variable-coefficient problem is stated as: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\sigma(x)\frac{d^2u}{dx^2} + \mu(x)\frac{du}{dx} + b(x)u = f(x),\quad
u(A)=\beta_0,\ u(B)=\beta_1,
$$
under conditions
$$
\sigma(x)\ge 0,\quad \mu(x)\ge \alpha > 0,\quad b(x)\le 0,\quad x\in (A,B).
$$

The fitted discrete scheme is: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\rho_j^h\, D_+D_- U_j + \mu_j D_0 U_j + b_j U_j = f_j,\quad j=1,\ldots,J-1,
\qquad
U_0=\beta_0,\ U_J=\beta_1,
$$
with the fitted coefficient
$$
\rho_j^h = \frac{\mu_j h}{2}\coth\left(\frac{\mu_j h}{2\sigma_j}\right),
\quad
\sigma_j=\sigma(x_j),\ \mu_j=\mu(x_j),\ b_j=b(x_j),\ f_j=f(x_j).
$$
The paper emphasizes that this choice enforces uniform stability/monotonicity properties that standard centered schemes can lose as $$\sigma\to 0$$. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

#### 3.4.3 Fully discrete exponentially fitted scheme for the parabolic problem (6)

A fitted discrete operator is defined (implicit in time, fitted in space) as: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
L_h^k U_j^n
\equiv
-\frac{U_j^{n+1}-U_j^n}{k}
+ \rho_j^{n+1} D_+D_- U_j^{n+1}
+ \mu_j^{n+1} D_0 U_j^{n+1}
+ b_j^{n+1} U_j^{n+1},
$$
with
$$
\rho_j^{n+1} \equiv \frac{\mu_j^{n+1} h}{2}\coth\left(\frac{\mu_j^{n+1} h}{2\sigma_j^{n+1}}\right),
\quad
\phi_j^{n+1}=\phi(x_j,t_{n+1}).
$$

The full scheme is then: find $$\{U_j^n\}$$ such that: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\begin{aligned}
L_h^k U_j^n &= f_j^{n+1},\quad j=1,\ldots,J-1,\ \ n=0,\ldots,N-1,\\
U_0^n &= g_0(t_n),\quad U_J^n=g_1(t_n),\quad n=0,\ldots,N,\\
U_j^0 &= \varphi(x_j),\quad j=1,\ldots,J-1.
\end{aligned}
$$

**Algorithm: exponentially fitted implicit solver for (6) (procedure).**

1. **Inputs:** $$\sigma(x,t),\mu(x,t),b(x,t),f(x,t)$$; grid parameters $$J,N,h,k$$; boundary data $$g_0,g_1$$; initial data $$\varphi$$.
2. **Mesh:** construct $$x_j=A+jh$$ and $$t_n=nk$$; evaluate initial line $$U_j^0=\varphi(x_j)$$ and boundary lines $$U_0^n=g_0(t_n),\ U_J^n=g_1(t_n)$$. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
3. **Time stepping:** for each $$n=0,\ldots,N-1$$:  
   3.1. Evaluate coefficients at $$t_{n+1}$$: $$\sigma_j^{n+1},\mu_j^{n+1},b_j^{n+1},f_j^{n+1}$$. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
   3.2. Compute fitted factors $$\rho_j^{n+1} = \frac{\mu_j^{n+1} h}{2}\coth\!\left(\frac{\mu_j^{n+1} h}{2\sigma_j^{n+1}}\right)$$ for interior nodes. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
   3.3. Assemble the tridiagonal linear system corresponding to $$L_h^k U_j^n = f_j^{n+1}$$ for $$j=1,\ldots,J-1$$ (implicit unknowns at level $$n+1$$). ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
   3.4. Solve for $$U^{n+1}$$ (direct tridiagonal solve implied by the stencil structure).
4. **Outputs:** discrete solution $$U_j^n$$ over the grid; option price and sensitivities can be computed from $$U$$ using appropriate discrete derivatives (the paper’s stated practical emphasis is improved price/delta behavior and reduced oscillations). ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

### 3.5 Limiting (“graceful degradation”) behavior for small volatility / small drift

The paper explicitly checks two extreme limits of the fitted parabolic scheme (16): pure convection $$\sigma\to 0$$ and pure diffusion $$\mu\to 0$$, using asymptotics of the hyperbolic cotangent. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

- **Limit $$\sigma\to 0$$ (pure convection):** the fitting coefficient satisfies: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
  $$
  \lim_{\sigma\to 0}\frac{\mu h}{2}\coth\left(\frac{\mu h}{2\sigma}\right)
  =
  \begin{cases}
  +\frac{\mu h}{2}, & \mu>0,\\
  -\frac{\mu h}{2}, & \mu<0.
  \end{cases}
  $$
  Substitution into (16) yields implicit upwind schemes:
  $$
  \text{if }\mu>0:\quad
  -\frac{U_j^{n+1}-U_j^n}{k}
  + \mu_j^{n+1}\frac{U_{j+1}^{n+1}-U_j^{n+1}}{h}
  + b_j^{n+1}U_j^{n+1}
  = f_j^{n+1},
  $$
  $$
  \text{if }\mu<0:\quad
  -\frac{U_j^{n+1}-U_j^n}{k}
  + \mu_j^{n+1}\frac{U_j^{n+1}-U_{j-1}^{n+1}}{h}
  + b_j^{n+1}U_j^{n+1}
  = f_j^{n+1}.
  $$
  The paper cites these as standard stable/convergent discretizations for the degenerate first-order limit. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

- **Limit $$\mu\to 0$$ (pure diffusion):** using $$\lim_{x\to 0}x\coth x=1$$, the fitted term reduces to the standard diffusion discretization:
  $$
  -\frac{U_j^{n+1}-U_j^n}{k}
  + \sigma_j^{n+1}D_+D_-U_j^{n+1}
  + b_j^{n+1}U_j^{n+1}
  = f_j^{n+1}.
  $$
  This is described as a standard approximation for pure diffusion problems. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

### 3.6 Keller “box” scheme (first-order system + box averaging)

Section 8.1 sketches a distinct alternative to CN: the Keller box scheme, motivated via a self-adjoint parabolic PDE in divergence form with boundary conditions involving derivatives. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

**Model PDE (as presented for motivation):** ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\begin{aligned}
\frac{\partial u}{\partial t}
&=
\frac{\partial}{\partial x}\left(a\frac{\partial u}{\partial x}\right)
+ cu + S,\quad 0<x<1,\ t>0,\\
u(x,0)&=g(x),\quad 0<x<1,\\
\alpha_0 u(0,t) + \alpha_1 a(0,t)u_x(0,t) &= g_0(t),\\
\beta_0 u(1,t) + \beta_1 a(1,t)u_x(1,t) &= g_1(t).
\end{aligned}
$$

**First-order reformulation:** introduce $$v$$ so that boundary conditions no longer involve derivatives: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
a u_x = v,\qquad
v_x = u_t - cu - S,
$$
with boundary conditions
$$
\alpha_0 u(0,t)+\alpha_1 v(0,t)=g_0(t),\qquad
\beta_0 u(1,t)+\beta_1 v(1,t)=g_1(t),
$$
and initial condition $$u(x,0)=g(x)$$.

**Box-averaging notation:** midpoints and averages are defined (paper notation): ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
x_{j\pm\frac{1}{2}}=\frac{1}{2}(x_j+x_{j\pm 1}),\quad
t_{n\pm\frac{1}{2}}=\frac{1}{2}(t_n+t_{n\pm 1}),
$$
$$
\phi_{j\pm\frac{1}{2}}^n=\frac{1}{2}\left(\phi_j^n+\phi_{j\pm 1}^n\right),\quad
\phi_j^{n\pm\frac{1}{2}}=\frac{1}{2}\left(\phi_j^n+\phi_j^{n\pm 1}\right).
$$
One-sided divided differences are: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
D_x^- \phi_j^n = \frac{\phi_j^n-\phi_{j-1}^n}{h_j},\qquad
D_t^- \phi_j^n = \frac{\phi_j^n-\phi_j^{n-1}}{k_n}.
$$

**Keller box scheme (as written):** solve for $$u$$ and $$v$$ simultaneously:
$$
a_{j-\frac{1}{2}}^n\, D_x^- u_j^n = v_{j-\frac{1}{2}}^n,
$$
$$
D_x^- v_j^{n-\frac{1}{2}}
=
D_t^- u_{j-\frac{1}{2}}^n
- c_{j-\frac{1}{2}}^{n-\frac{1}{2}}u_{j-\frac{1}{2}}^{n-\frac{1}{2}}
- S_{j-\frac{1}{2}}^{n-\frac{1}{2}},
\qquad
1\le j\le J,\ 1\le n\le N.
$$
Boundary conditions are: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
\alpha_0 u_0^n + \alpha_1 v_0^n = g_0^n,\qquad
\beta_0 u_J^n + \beta_1 v_J^n = g_1^n,\qquad
1\le n\le N.
$$

**Initialization for piecewise smooth data:** the scheme defines an approximate initial condition for $$v$$ via the derivative of $$g$$: ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  
$$
v_{j-\frac{1}{2}}^0 = a_{j-\frac{1}{2}}^0 \frac{dg}{dx}\bigg|_{x_{j-\frac{1}{2}}},\qquad 1\le j\le J.
$$
For piecewise smooth boundary conditions, mid-time boundary enforcement is proposed:
$$
\alpha_0 u_0^{n-\frac{1}{2}} + \alpha_1 v_0^{n-\frac{1}{2}} = g_0^{n-\frac{1}{2}},\qquad
\beta_0 u_J^{n-\frac{1}{2}} + \beta_1 v_J^{n-\frac{1}{2}} = g_1^{n-\frac{1}{2}},
$$
with the explicit caveat that mesh points are assumed to coincide with discontinuity times $$t=t_n$$. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

**Properties claimed:** unconditional stability; second-order accuracy for $$u$$ and $$u_x$$ (interpreted for Black–Scholes as price and delta); applicability of Richardson extrapolation with two orders of improvement per extrapolation (including on nonuniform meshes); robustness for piecewise smooth data (payoffs). ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

### 3.7 ASCII overview diagram (end-to-end pipeline)

```
┌──────────────────────────────────────────┐
│ PDE coefficients (σ, μ, b, f) + IC/BC     │
│ (payoff ϕ, boundaries g0,g1)              │
└───────────────────────┬──────────────────┘
                        │ choose mesh: {xj},{tn}, h,k
                        ▼
┌──────────────────────────────────────────┐
│ (Option A) CN / centered FD              │
│ (Option B) Fitted implicit FD            │
│ (Option C) Keller box (u,v) system       │
└───────────────────────┬──────────────────┘
                        │ build linear system per time step
                        ▼
┌──────────────────────────────────────────┐
│ Solve per level (typically tridiagonal)  │
│ + enforce discrete boundary conditions    │
└───────────────────────┬──────────────────┘
                        │ post-process
                        ▼
┌──────────────────────────────────────────┐
│ Price grid + Greeks (Δ, Γ) near strike    │
│ with oscillation control as primary goal  │
└──────────────────────────────────────────┘
```
([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

## 4. Theoretical Results

### 4.1 Stationary fitted scheme (12): uniform stability, monotonicity, convergence (unnumbered “fundamental results”)

**Claim 1 (uniform stability bound).** Under assumptions $$\sigma(x)\ge 0$$, $$\mu(x)\ge \alpha>0$$, $$b(x)\le 0$$ on $$x\in(A,B)$$, the paper states that the solution of scheme (12) is **uniformly stable**, i.e.:
$$
|U_j|
\le
|\beta_0| + |\beta_1| + \frac{1}{\alpha}\max_{k=1,\ldots,J}|f_k|,
\qquad j=1,\ldots,J-1.
$$
**Proof sketch:** the result is attributed to the fitted-scheme literature (Il’in; Duffy) and follows from a discrete maximum principle argument: lower-bounding convection by $$\alpha$$ and enforcing the M-matrix sign structure yields an $$\ell_\infty$$ bound controlled by boundary magnitudes plus a source-term contribution. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

**Claim 2 (monotonicity / positivity preservation).** The matrix representation $$AU=F$$ for scheme (12) is said to be **monotone**, producing positive discrete solutions from positive inputs, with tridiagonal coefficients:
$$
a_{j,j-1}=\frac{\rho_j^h}{h^2}-\frac{\mu_j}{2h} > 0,\qquad
a_{j,j}=-\frac{2\rho_j^h}{h^2}+b_j < 0,\qquad
a_{j,j+1}=\frac{\rho_j^h}{h^2}+\frac{\mu_j}{2h} > 0.
$$
**Proof sketch:** positivity of off-diagonal entries and negativity of the diagonal entry are the structural conditions that place $$A$$ in the class of M-matrices used in monotone finite-difference theory; the fitted coefficient $$\rho_j^h$$ is constructed precisely to maintain these inequalities independent of the (potentially small) diffusion parameter. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

**Claim 3 (first-order uniform convergence in space).** Let $$u$$ be the exact solution of (11) and $$U$$ the fitted discrete solution of (12). The paper states:
$$
|u(x_j)-U_j| \le Mh,
$$
where $$M$$ is independent of $$h$$ and $$\sigma$$.  
**Proof sketch:** the fitted discretization is designed to match boundary-layer structure induced by small diffusion, yielding a uniform-in-$$\sigma$$ truncation and stability bound; combining stability with consistency gives a uniform error estimate (the statement is explicitly attributed to Il’in). ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

### 4.2 Lemma 1 (discrete maximum principle for the fitted parabolic scheme)

**Lemma 1 (statement).** If a discrete function $$w_j^n$$ satisfies:
- $$L_h^k w_j^n \le 0$$ in the interior of the mesh, and  
- $$w_j^n \ge 0$$ on the boundary (denoted $$\partial\Omega$$ in the PDF text),
then:
$$
w_j^n \ge 0,\quad \forall j=0,\ldots,J,\ \forall n=0,\ldots,N.
$$
**Proof sketch (paper’s structure).** The inequality is rewritten as a vector inequality $$A^n W^{n+1}\ge W^n$$ where $$W^n=(w_1^n,\ldots,w_{J-1}^n)^T$$ and $$A^n$$ is a tridiagonal matrix whose off-diagonal elements are non-positive, whose diagonal elements are strictly positive, and which is irreducibly diagonally dominant. The cited matrix-theory result (Varga) implies $$A^n$$ is nonsingular and $$\left(A^n\right)^{-1}\ge 0$$, giving the discrete maximum principle. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

### 4.3 Lemma 2 (sup-norm bound via barrier function)

**Lemma 2 (statement).** Let $$\{U_j^n\}$$ solve scheme (16) and assume:
$$
\max |U_j^n| \le m \ \text{on the boundary for all }j,n,\qquad
\max |f_j^n| \le N \ \text{in }D\text{ for all }j,n.
$$
Then:
$$
\max_j |U_j^n| \le -\frac{N}{\beta} + m \quad \text{in }\overline{D}.
$$
**Ambiguity flag:** the constant $$\beta$$ is used but not explicitly defined in the displayed lemma text; it is consistent with a negative bound on the reaction coefficient $$b(x,t)$$ (or a related coercivity constant) needed to build the barrier term $$-N/\beta$$. The statement is reproduced as written. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

**Proof sketch:** define a barrier function $$w_j^n=-\frac{N}{\beta}+m \pm U_j^n$$, verify nonnegativity on the boundary, and show $$L_h^k w_j^n \le 0$$ in the interior; Lemma 1 then yields $$w_j^n\ge 0$$ everywhere, which implies the stated bound. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

### 4.4 Uniform (in $$\sigma$$) convergence estimate for the parabolic fitted scheme

**Claim (error estimate (18)).** Let $$u(x,t)$$ solve (6) and $$U_j^n$$ solve (16). Then:
$$
|u(x_j,t_n)-U_j^n|\le M(h+k),
$$
where $$M$$ is independent of $$h$$, $$k$$, and $$\sigma$$.  
**Proof sketch:** the paper’s argument is that stability is proven via the discrete maximum principle (not von Neumann analysis), yielding bounds independent of small diffusion; combining stability with consistency yields an error estimate where the constant does not deteriorate as $$\sigma\to 0$$, unlike classical centered-in-space + CN-in-time discretizations. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

**Complexity bounds:** no explicit asymptotic computational complexity bounds are stated beyond the repeated solution of tridiagonal linear systems at each time level. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

## 5. Experimental Evaluation

### 5.1 Test problems, baselines, metrics, and hyperparameters (as stated)

| Category | Items explicitly mentioned in the paper | Notes / gaps (not filled in) |
|---|---|---|
| “Datasets” / instruments | Double barrier call; single barrier call; time-dependent volatility cases (example: linear function of time is mentioned in the option list; decaying exponential $$\sigma(t)=\sigma_0 e^{-\alpha(T-t)}$$ is discussed as problematic for CN); asymmetric plain-vanilla power call; asymmetric capped power call. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) | No strike/tenor/spot parameter sets tabulated; no grid-resolution table. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) |
| Boundary/initial-condition features | Dirichlet boundaries on finite interval for barriers; non-smooth payoff at strike (European call kink); compatibility conditions at corners $$\varphi(A)=g_0(0)$$ and $$\varphi(B)=g_1(0)$$. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) | No explicit discrete boundary formulas beyond discussion and box-scheme tactics. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) |
| Baselines / comparators | Haug (1998) formulas; Topper (1998); Monte Carlo; comparisons stated to be “favourable”. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) | No explicit numeric error table, confidence intervals, or convergence plots provided in the article text. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) |
| Metrics | Option price; delta; gamma is discussed as problematic under CN; “good approximations to price and delta” claimed for the fitted method. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) | No explicit norm definitions (e.g., $$\ell_\infty$$ error) given for experiments. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) |
| Hyperparameters | Uniform mesh sizes $$h$$ and $$k$$ are defined for fitted scheme derivations; tridiagonal solves per time step. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) | Concrete choices of $$J$$, $$N$$, domain truncation $$[A,B]$$, and smoothing/damping parameters are not tabulated. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) |
| Ablations | None presented as formal ablation studies; section-by-section critique isolates mechanisms (boundary conditions, initial discontinuities, mesh nonuniformity, volatility degeneracy). ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) | No controlled ablation axes with measured deltas are reported. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)) |

### 5.2 Quantitative results

No tables/figures with explicit numerical values (prices, deltas, gammas, convergence rates) are provided in the article text; performance is reported qualitatively (e.g., “compare favourably”) with reference to external sources and a working paper. ([ma.imperial.ac.uk](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf))  

## 6. ASCII Architecture / Workflow Diagram(s)

### 6.1 Exponentially fitted implicit solver for (6) / Black–Scholes-type PDEs (scheme (16))

```
┌──────────────────────────────────────────────────────────┐
│ Inputs                                                   │
│  σ(x,t), μ(x,t), b(x,t), f(x,t)                           │
│  IC: Uj^0 = ϕ(xj)                                         │
│  BC: U0^n = g0(tn), UJ^n = g1(tn)                          │
└───────────────────────────────┬──────────────────────────┘
                                │ uniform mesh: xj, tn (h,k)
                                ▼
┌──────────────────────────────────────────────────────────┐
│ For each time step n→n+1                                  │
│  Compute ρj^{n+1} = (μj^{n+1} h / 2) coth( μj^{n+1} h /     │
│                                           (2 σj^{n+1}) )    │
└───────────────────────────────┬──────────────────────────┘
                                │ assemble fitted operator L_h^k
                                ▼
┌──────────────────────────────────────────────────────────┐
│ Tridiagonal linear system for {Uj^{n+1}}_{j=1}^{J-1}       │
│  -(Uj^{n+1}-Uj^n)/k + ρj^{n+1} D+ D- Uj^{n+1}              │
│  + μj^{n+1} D0 Uj^{n+1} + b_j^{n+1} Uj^{n+1} = f_j^{n+1}    │
└───────────────────────────────┬──────────────────────────┘
                                │ solve (direct tridiagonal)
                                ▼
┌──────────────────────────────────────────────────────────┐
│ Outputs                                                   │
│  Discrete surface Uj^n                                     │
│  Greeks via finite differences (Δ, Γ) with oscillation      │
│  control emphasized near strike and in small-σ regimes      │
└──────────────────────────────────────────────────────────┘
```

### 6.2 Keller box scheme workflow (first-order system + box averaging)

```
┌──────────────────────────────────────────────────────────┐
│ Start: self-adjoint parabolic PDE + derivative BCs         │
│  ut = (a ux)x + c u + S                                    │
│  α0 u(0,t)+α1 a(0,t) ux(0,t)=g0(t)                          │
│  β0 u(1,t)+β1 a(1,t) ux(1,t)=g1(t)                          │
└───────────────────────────────┬──────────────────────────┘
                                │ introduce v = a ux
                                ▼
┌──────────────────────────────────────────────────────────┐
│ First-order system                                         │
│  a ux = v                                                  │
│  vx = ut - c u - S                                         │
│  α0 u(0,t)+α1 v(0,t)=g0(t), β0 u(1,t)+β1 v(1,t)=g1(t)       │
└───────────────────────────────┬──────────────────────────┘
                                │ box averages + one-sided diffs
                                ▼
┌──────────────────────────────────────────────────────────┐
│ Solve simultaneously for u and v per time level             │
│  a_{j-1/2}^n D_x^- u_j^n = v_{j-1/2}^n                      │
│  D_x^- v_j^{n-1/2} = D_t^- u_{j-1/2}^n - c u - S            │
└───────────────────────────────┬──────────────────────────┘
                                ▼
┌──────────────────────────────────────────────────────────┐
│ Outputs: u (price) and v≈a ux (delta-like) with claimed     │
│ second-order accuracy and reduced oscillations              │
└──────────────────────────────────────────────────────────┘
```

## 7. Follow-Up Works & Extensions

### 7.1 CN oscillation control via smoothing and time transformations

Wade, Khaliq, Yousuf, Vigo-Aguiar, and Deininger develop a smoothing strategy for CN tailored to barrier options, arguing that standard start-up smoothing (e.g., Luskin–Rannacher-style Euler start) can be insufficient when discontinuities are repeatedly injected at discrete barrier-monitoring times; they propose reapplying damping at each discontinuity time and extend the idea to higher-order Padé-based schemes. This directly operationalizes the paper’s critique that CN oscillations are triggered by nonsmooth data and that hedging Greeks (delta/gamma) are particularly sensitive. [Wade et al., JCAM 2007]. ([pure.kfupm.edu.sa](https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/))  

Reisinger and Whitley analyze the convergence pathology of CN with singular/nonsmooth initial data and show that a square-root time change can restore convergence properties; their numerical evidence targets European and American options and explicitly reports improved convergence for price, delta, and gamma without relying on Rannacher start-up steps. This extends the paper’s theme that CN’s “stability” does not imply accuracy for Greeks near kinks and that modifications to time discretization can be decisive. [Reisinger & Whitley, IMA JNA 2014]. ([academic.oup.com](https://academic.oup.com/imajna/article/34/3/1156/717033))  

Ma and Zhou study “moving mesh” implicit finite differences for Asian-option PDEs with moving boundaries and build the time stepping around Rannacher’s idea (a few backward-Euler steps followed by CN), proving second-order convergence in both time and space under their graded-time/moving-space construction and validating with numerical examples. This is a concrete follow-up in the same design space as the paper’s critique: non-smooth initial/final data plus convection/diffusion imbalances motivate damping/robustness modifications beyond vanilla CN. [Ma & Zhou, J. Comput. Math. 2016]. ([global-sci.org](https://global-sci.org/index.php/JCM/article/view/12230))  

### 7.2 Monotone/fitted discretizations and maximum-principle analysis in option-pricing PDEs

Wang introduces a fitted finite-volume spatial discretization for the degenerate Black–Scholes PDE and proves stability and an M-matrix (discrete maximum principle) property for the resulting system, with numerical experiments demonstrating effectiveness. This aligns closely with the paper’s “exponential fitting” narrative: preserving monotonicity/DMP behavior is the mechanism used to avoid spurious oscillations in convection-dominated/degenerate regimes. [Wang, IMA JNA 2004]. ([academic.oup.com](https://academic.oup.com/imajna/article/24/4/699/687386))  

Pooley, Forsyth, and Vetzal analyze convergence for option-pricing PDEs under uncertain volatility, emphasizing that non-monotone discretizations (the paper gives standard CN as an example) can converge to incorrect solutions or become unstable, and they provide both convergence theory for an iterative implicit approach and numerical examples. This provides a finance-specific theoretical backing for the paper’s warning that “stable” schemes can still be unreliable when monotonicity/viscosity-solution consistency is violated. [Pooley et al., IMA JNA 2003]. ([academic.oup.com](https://academic.oup.com/imajna/article-abstract/23/2/241/684490))  

Zhang and Wang extend fitted finite-volume spatial discretization to partial integro-differential equations arising under jump-diffusion dynamics, using CN for time discretization and combining an iterative solver with FFT acceleration for the integral term; they also address American-style constraints via a penalty term and highlight an M-matrix property for the discretized system matrix. This connects to the paper’s broader claim that robust discretizations in space (monotone/M-matrix) can be paired with time stepping choices, and that oscillation control is intertwined with preserving positivity/monotonicity under more complex dynamics. [Zhang & Wang, Appl. Math. Comput. 2008]. ([research-repository.uwa.edu.au](https://research-repository.uwa.edu.au/en/publications/pricing-options-under-jump-diffusion-processes-with-fitted-finite))  

Zhang and Wang propose a fitted finite-volume-based scheme for the nonlinear uncertain-volatility PDE, proving consistency, stability, and monotonicity (hence convergence to the viscosity solution) and presenting numerical experiments. This is a direct continuation of the monotonicity/maximum-principle line that the paper advocates as more appropriate than von Neumann analysis for variable-coefficient/nonlinear finance PDEs. [Zhang & Wang, Appl. Numer. Math. 2009]. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0168927409000051))  

Ramírez-Espinoza and Ehrhardt compare conservative/finite-volume approaches for convection-dominated Black–Scholes-type problems, explicitly including exponentially fitted methods among the examined schemes and reporting experiments across different linear/nonlinear pricing PDEs. This work functions as a later synthesis that situates exponential fitting (as promoted in the paper) among other high-resolution/conservative discretizations for the same convection-dominance failure mode. [Ramírez-Espinoza & Ehrhardt, AAMM 2013]. ([global-sci.com](https://global-sci.com/index.php/aamm/article/view/4690))  

### 7.3 Keller box scheme usage in option pricing (domain-specific adaptations)

Mardianto, Putra, Pratama, and Putri apply the Keller box method to the Black–Scholes model for a European put, presenting a numerical study and checking stability via von Neumann analysis. This reflects direct downstream adoption of the paper’s “Keller scheme” alternative in a finance PDE setting, even though the paper’s own emphasis is on the box scheme’s unconditional stability and second-order accuracy for both the solution and its spatial derivative. [Mardianto et al., J. Mat. UNAND 2024]. ([doaj.org](https://doaj.org/article/8994136a0b9e4c8bbb2998eb9159db45))  

## 8. Industrial & Real-World Applications

QuantLib is a widely used open-source C++ quantitative finance library with an explicit finite-differences framework that includes CN-style schemes (e.g., a Crank–Nicolson evolver in the finite-differences module) as well as alternative time-stepping methods (e.g., explicit/implicit Euler and other splitting schemes, depending on the interface). This directly matches the paper’s practical setting (production-grade PDE pricers where CN is common but requires care), and it provides a concrete, verifiable codebase where the paper’s recommended alternatives (implicit steps, alternative schemes, monotone discretizations) can be implemented and compared. [GitHub: lballabio/QuantLib]. ([github.com](https://github.com/lballabio/QuantLib))  

The Open Source Risk Engine (ORE) is a large open-source risk analytics system (C++), intended for real-life risk calculations across portfolios; its public repository and activity indicate industrial-scale engineering around quantitative models and numerical methods. It is relevant as an example of production-oriented infrastructure into which robust PDE solvers (including CN variants or fitted/monotone schemes) can be integrated when pricing or risk models reduce to PDE/PIDE forms. [GitHub: OpenSourceRisk/Engine]. ([github.com](https://github.com/OpenSourceRisk/Engine))  

TF Quant Finance (google/tf-quant-finance) provides open-source implementations of quantitative finance components including PDE solvers and example notebooks (the repository explicitly lists “ODE & PDE solvers” as a mid-level tier and includes an “American Option pricing under the Black-Scholes model” example), offering an industrial-grade, GPU-friendly environment where issues raised in the paper (CN damping, payoff nonsmoothness, stability vs. accuracy for Greeks) can be explored at scale. The repository is explicitly marked archived/unmaintained, which constrains “production deployment” claims to the existence and reproducibility of the code rather than ongoing supported use. [GitHub: google/tf-quant-finance]. ([github.com](https://github.com/google/tf-quant-finance))  

## 9. Consolidated Reference List

[1] B. A. Wade, A. Q. M. Khaliq, M. Yousuf, J. Vigo-Aguiar, R. Deininger. “On smoothing of the Crank–Nicolson scheme and higher order schemes for pricing barrier options.” *Journal of Computational and Applied Mathematics*, 2007. DOI: `https://doi.org/10.1016/j.cam.2006.04.034`. ([pure.kfupm.edu.sa](https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/))  

[2] C. Reisinger, A. Whitley. “The impact of a natural time change on the convergence of the Crank–Nicolson scheme.” *IMA Journal of Numerical Analysis*, 2014. DOI: `https://doi.org/10.1093/imanum/drt029`. arXiv: `https://arxiv.org/abs/1210.5487`. ([academic.oup.com](https://academic.oup.com/imajna/article/34/3/1156/717033))  

[3] Jingtang Ma, Zhiqiang Zhou. “Convergence Rates of Moving Mesh Rannacher Methods for PDEs of Asian Options Pricing.” *Journal of Computational Mathematics*, Vol. 34 No. 3 (2016). DOI: `https://doi.org/10.4208/jcm.1601-m2014-0217`. ([global-sci.org](https://global-sci.org/index.php/JCM/article/view/12230))  

[4] D. M. Pooley, P. A. Forsyth, K. R. Vetzal. “Numerical convergence properties of option pricing PDEs with uncertain volatility.” *IMA Journal of Numerical Analysis*, 2003. DOI: `https://doi.org/10.1093/imanum/23.2.241`. ([academic.oup.com](https://academic.oup.com/imajna/article-abstract/23/2/241/684490))  

[5] Song Wang. “A novel fitted finite volume method for the Black–Scholes equation governing option pricing.” *IMA Journal of Numerical Analysis*, 2004. DOI: `https://doi.org/10.1093/imanum/24.4.699`. ([academic.oup.com](https://academic.oup.com/imajna/article/24/4/699/687386))  

[6] K. Zhang, Song Wang. “Pricing options under jump diffusion processes with fitted finite volume method.” *Applied Mathematics and Computation*, 2008. DOI: `https://doi.org/10.1016/j.amc.2007.12.043`. ([research-repository.uwa.edu.au](https://research-repository.uwa.edu.au/en/publications/pricing-options-under-jump-diffusion-processes-with-fitted-finite))  

[7] K. Zhang, Song Wang. “A computational scheme for uncertain volatility model in option pricing.” *Applied Numerical Mathematics*, 2009. DOI: `https://doi.org/10.1016/j.apnum.2009.01.004`. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0168927409000051))  

[8] Germán I. Ramírez-Espinoza, Matthias Ehrhardt. “Conservative and Finite Volume Methods for the Convection-Dominated Pricing Problem.” *Advances in Applied Mathematics and Mechanics*, Vol. 5 No. 6 (2013). DOI: `https://doi.org/10.4208/aamm.12-m1216`. ([global-sci.com](https://global-sci.com/index.php/aamm/article/view/4690))  

[9] Lutfi Mardianto, Gusrian Putra, Benediktus Ivan Pratama, Endah R. M. Putri. “Numerical Solution of European Put Option for Black-Scholes Model Using Keller Box Method.” *Jurnal Matematika UNAND*, 2024. DOI: `https://doi.org/10.25077/jmua.13.3.188-197.2024`. ([doaj.org](https://doaj.org/article/8994136a0b9e4c8bbb2998eb9159db45))  

[10] QuantLib contributors. “QuantLib (C++ library).” GitHub repository: `https://github.com/lballabio/QuantLib`. ([github.com](https://github.com/lballabio/QuantLib))  

[11] The QuantLib contributors. “QuantLib: a free/open-source library for quantitative finance.” Zenodo software record (example release record). DOI: `https://doi.org/10.5281/zenodo.10986753`. ([zenodo.org](https://zenodo.org/records/10986753))  

[12] Open Source Risk Engine contributors. “Open Source Risk Engine (ORE).” GitHub repository: `https://github.com/OpenSourceRisk/Engine`. ([github.com](https://github.com/OpenSourceRisk/Engine))  

[13] Google. “TF Quant Finance (tf-quant-finance).” GitHub repository: `https://github.com/google/tf-quant-finance`. ([github.com](https://github.com/google/tf-quant-finance))

---
Learn more:
1. [https://www.ma.imperial.ac.uk/~ajacquie/IC\_Num\_Methods/IC\_Num\_Methods\_Docs/Literature/DuffyCN.pdf](https://www.ma.imperial.ac.uk/~ajacquie/IC_Num_Methods/IC_Num_Methods_Docs/Literature/DuffyCN.pdf)
2. [https://www.osti.gov/biblio/2502172](https://www.osti.gov/biblio/2502172)
3. [https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/](https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/)
4. [https://academic.oup.com/imajna/article/34/3/1156/717033](https://academic.oup.com/imajna/article/34/3/1156/717033)
5. [https://global-sci.org/index.php/JCM/article/view/12230](https://global-sci.org/index.php/JCM/article/view/12230)
6. [https://academic.oup.com/imajna/article/24/4/699/687386](https://academic.oup.com/imajna/article/24/4/699/687386)
7. [https://academic.oup.com/imajna/article-abstract/23/2/241/684490](https://academic.oup.com/imajna/article-abstract/23/2/241/684490)
8. [https://research-repository.uwa.edu.au/en/publications/pricing-options-under-jump-diffusion-processes-with-fitted-finite](https://research-repository.uwa.edu.au/en/publications/pricing-options-under-jump-diffusion-processes-with-fitted-finite)
9. [https://www.sciencedirect.com/science/article/pii/S0168927409000051](https://www.sciencedirect.com/science/article/pii/S0168927409000051)
10. [https://global-sci.com/index.php/aamm/article/view/4690](https://global-sci.com/index.php/aamm/article/view/4690)
11. [https://doaj.org/article/8994136a0b9e4c8bbb2998eb9159db45](https://doaj.org/article/8994136a0b9e4c8bbb2998eb9159db45)
12. [https://github.com/lballabio/QuantLib](https://github.com/lballabio/QuantLib)
13. [https://github.com/OpenSourceRisk/Engine](https://github.com/OpenSourceRisk/Engine)
14. [https://github.com/google/tf-quant-finance](https://github.com/google/tf-quant-finance)
15. [https://zenodo.org/records/10986753](https://zenodo.org/records/10986753)

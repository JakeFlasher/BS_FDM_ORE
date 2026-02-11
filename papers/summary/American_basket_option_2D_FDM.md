## 1. Paper Identity

- **Title:** *American-type basket option pricing: a simple two-dimensional partial differential equation*  
- **Authors:** Hamza Hanbali; Daniel Linders  
- **Affiliations:**  
  - Hamza Hanbali: Department of Accountancy, Finance and Insurance, KU Leuven, Naamsestraat 69, Leuven 3000, Belgium  
  - Daniel Linders: Department of Mathematics, University of Illinois, Urbana, IL, USA  
- **Venue:** *Quantitative Finance*, Volume 19, Issue 10, pages 1689–1704  
- **Year / publication dates:** 2019 (published online 12 April 2019; received 13 April 2018; accepted 22 February 2019)  
- **DOI:** 10.1080/14697688.2019.1588987  

## 2. Problem Statement & Formulation

The paper targets efficient pricing (and extraction of hedging ratios and exercise regions) for **American-type** derivatives written on a **basket** of $$n$$ correlated stocks in the multivariate Black–Scholes setting. The central computational obstacle is the curse of dimensionality in PDE-based valuation: the classical PDE for a basket derivative is $$n$$-dimensional with mixed second derivatives, making stable finite-difference valuation impractical when $$n$$ is moderate-to-large.

### Multivariate Black–Scholes market (as specified in Section 2.1)

There are $$n$$ non-dividend-paying stocks with prices $$S_i(t)$$, $$i \in \{1,\dots,n\}$$, whose real-world dynamics (paper’s notation) are
$$
S_i(t) = S_i(0) + \int_0^t \mu_i S_i(s)\,ds + \int_0^t \sigma_i S_i(s)\,dB_i(s),
\quad t>0,
$$
with correlated Brownian motion increments satisfying
$$
\mathbb{E}\!\left[dB_i(t)\,dB_j(t)\right] = \rho_{i,j}\,dt.
$$

The bank account $$D(t)$$ follows
$$
\frac{dD(t)}{D(t)} = r\,dt,
$$
with constant deterministic rate $$r$$.

The basket price is a weighted sum
$$
S(t) = \sum_{i=1}^n w_i S_i(t),
$$
with fixed positive weights $$w_i>0$$.

### Basket derivative pricing problem

Let maturity be $$T$$ and payoff function be $$H$$ applied to the basket. For a European claim, payoff is $$H(S(T))$$. Denote the arbitrage-free price at time $$t$$ by
$$
V(t,S_1(t),\dots,S_n(t)).
$$

In the multivariate Black–Scholes framework, the paper states that the European basket derivative price satisfies the $$n$$-dimensional PDE
$$
\frac{\partial V}{\partial t}
+ \frac{1}{2}\sum_{i=1}^n\sum_{j=1}^n \sigma_i\sigma_j\rho_{i,j} w_i w_j S_i S_j
\frac{\partial^2 V}{\partial S_i \partial S_j}
+ r \sum_{i=1}^n S_i \frac{\partial V}{\partial S_i}
- rV
= 0,
$$
together with terminal and boundary conditions determined by $$H$$ and contract features.

For an American-type claim, the price can be expressed as an optimal stopping problem (standard risk-neutral formulation, consistent with the paper’s early-exercise treatment in finite differences):
$$
V(t,S_1(t),\dots,S_n(t))
=
\sup_{\tau \in \mathcal{T}_{t,T}}
\mathbb{E}^Q\!\left[
e^{-r(\tau-t)} H(S(\tau))
\,\middle|\,
\mathcal{F}_t
\right],
$$
where $$\mathcal{T}_{t,T}$$ is the set of stopping times valued in $$[t,T]$$ and $$Q$$ is the risk-neutral measure.

### Paper’s dimensionality-reduction formulation via comonotonicity (Sections 2.3–2.4)

For **convex** payoff $$H$$, the paper uses comonotonicity to obtain bounds:
- a **comonotonic upper bound** price $$V^c(t,S^c(t))$$ obtained by forcing perfect positive dependence (all assets driven by one Brownian motion with unchanged marginal volatilities),
- a **comonotonic lower bound** price $$V^l(t,S^l(t))$$ obtained by forcing comonotonicity but reducing marginal volatilities via factors $$\nu_i \in [0,1]$$.

The paper then approximates the “real” (correlated) basket option by solving only **two one-factor (two-dimensional) problems** (upper and lower), and optionally combining them via a variance-matching weight (Section 5).

## 3. Core Methodology

### 3.1 Comonotonic bound construction (Sections 2.3–2.4)

#### Comonotonic upper-bound market (Section 2.3)

Define an **artificial comonotonic market** with stock prices $$S_i^c(t)$$ driven by a single Brownian motion $$B(t)$$:
$$
S_i^c(t) = S_i(0) + \int_0^t \mu_i S_i^c(s)\,ds + \int_0^t \sigma_i S_i^c(s)\,dB(s).
$$

The comonotonic basket is
$$
S^c(t) = \sum_{i=1}^n w_i S_i^c(t).
$$

Key structural property used by the paper: for convex $$H$$, the comonotonic basket derivative price $$V^c(t,S^c(t))$$ is an **upper bound** for the original basket derivative price.

#### Comonotonic lower-bound market (Section 2.4)

Define another artificial comonotonic market with reduced volatilities:
$$
S_i^l(t) = S_i(0) + \int_0^t \mu_i S_i^l(s)\,ds + \int_0^t \nu_i \sigma_i S_i^l(s)\,dB(s),
\quad 0 \leq \nu_i \leq 1.
$$
The basket is $$S^l(t)=\sum_{i=1}^n w_i S_i^l(t)$$ and the derivative price is $$V^l(t,S^l(t))$$.

The paper cites a sufficient parameterization for $$\nu_i$$ (Lemma 2.1) that yields a **lower bound** $$V^l \leq V$$ for convex payoffs. For numerical work in Section 4.2 onward, a specific choice of $$\nu_i$$ is used (equation (26)) that plugs weights and initial spots into the lower-bound construction.

### 3.2 One-factor representation and state variable choice (Section 3.1–3.2)

Because all stocks are driven by one Brownian motion in the comonotonic market, the basket can be written as a deterministic increasing function of that driver. The paper defines functions
$$
f_i(x) = w_i S_i(0)\exp\!\left( (\mu_i - \tfrac{1}{2}\sigma_i^2)t + \sigma_i x \right),
\quad
f(x) = \sum_{i=1}^n f_i(x),
$$
and uses strict monotonicity of $$f$$ to argue that the filtration generated by $$S^c(t)$$ matches the filtration of the single Brownian factor, enabling reduction to a 2D PDE in $$ (t, B_\lambda) $$ (Section 3.2).

### 3.3 Risk-neutral drift handling via a shifted Brownian motion (Section 3.1–3.2)

The paper assumes no-arbitrage in the comonotonic market via a common market price of risk $$\lambda$$ satisfying
$$
\lambda = \frac{\mu_i - r}{\sigma_i},
\quad i=1,\dots,n.
$$

Define the shifted process $$B_\lambda(t)$$ by
$$
dB_\lambda(t) = dB(t) + \lambda\,dt.
$$

Under the risk-neutral measure $$Q$$, $$B_\lambda$$ becomes a standard Brownian motion; the paper proceeds under the objective measure but treats $$B_\lambda$$ as a Brownian motion with drift $$\lambda$$ under $$P$$, while using risk-neutral arguments for pricing.

Under this representation, the comonotonic basket can be expressed as
$$
S^c(t) = \sum_{i=1}^n w_i S_i(0)
\exp\!\left( (r - \tfrac{1}{2}\sigma_i^2)t + \sigma_i B_\lambda(t) \right).
$$

### 3.4 European comonotonic PDE and replicating hedge (Section 3.2.1)

Define the option price as a function of $$ (t, B_\lambda) $$:
$$
V_\lambda^c(t,B_\lambda) =
V^c\!\left(t,\;\sum_{i=1}^n w_i S_i(0)\exp\!\left( (r - \tfrac{1}{2}\sigma_i^2)t + \sigma_i B_\lambda \right)\right).
$$

The paper derives a **constant-coefficient** backward PDE (Theorem 3.1):
$$
\frac{\partial V_\lambda^c}{\partial t}
+\frac{1}{2}\frac{\partial^2 V_\lambda^c}{\partial B_\lambda^2}
- r V_\lambda^c
=0,
$$
with terminal condition (payoff at maturity)
$$
V_\lambda^c(T,B_\lambda) =
H\!\left(\sum_{i=1}^n w_i S_i(0)\exp\!\left( (r - \tfrac{1}{2}\sigma_i^2)T + \sigma_i B_\lambda \right)\right).
$$

The replicating portfolio in the comonotonic market holds:
- $$\pi_1(t)$$ units of the comonotonic basket $$S^c(t)$$
- $$\pi_2(t)$$ units of the bank account $$D(t)$$

with (paper’s statement)
$$
\pi_1(t)=\frac{1}{\sigma^c(t)}\frac{\partial V^c}{\partial B_\lambda},
\quad
\pi_2(t)=\frac{1}{D(t)}\left(V^c(t,B_\lambda)-\pi_1(t)S^c(t)\right),
$$
where $$\sigma^c(t)=\sum_{i=1}^n \sigma_i w_i S_i^c(t)$$ (from Lemma 3.1).  
Notation ambiguity: The theorem labels the PDE unknown as $$V_\lambda^c$$ but writes $$\partial V^c/\partial B_\lambda$$ in the hedge; the intended object is the derivative of the price function with respect to the Brownian state variable.

### 3.5 Path-dependent extension (continuously sampled Asian-type basket options, Section 3.2.2)

Define the path functional
$$
I(t)=\int_0^t g(u,S^c)\,du,
$$
where $$g$$ specifies the averaging rule. The price depends on $$ (t,I,B_\lambda) $$ and satisfies
$$
\frac{\partial V_\lambda^c}{\partial t}
+ g(t,S^c)\frac{\partial V_\lambda^c}{\partial I}
+\frac{1}{2}\frac{\partial^2 V_\lambda^c}{\partial B_\lambda^2}
- rV_\lambda^c
=0.
$$
The paper does not provide a closed form for this PDE, but notes its structural similarity to the 1D Asian-option PDE.

### 3.6 Closed-form integral representation for the European comonotonic price (Section 3.3)

The paper provides an explicit solution to the European comonotonic PDE as a 1D Gaussian integral (Theorem 3.2), enabling deterministic quadrature if desired. This representation is also used as conceptual support for the finite-difference PDE solver.

### 3.7 Explicit finite difference scheme for the comonotonic PDE (Section 4.1)

The explicit scheme solves the constant-coefficient PDE in variables $$ (t,B_\lambda) $$.

- Time grid (backward):  
  $$t_k = T - k\delta t,\quad k=0,1,\dots,N.$$
- Space grid for $$B_\lambda$$:  
  $$b_j = (j-I)\delta B,\quad j=0,1,\dots,J,$$
  with chosen integers $$I,J$$.

Let $$V_j^k$$ approximate $$V_\lambda^c(t_k,b_j)$$. The paper uses finite differences:
$$
\frac{\partial V_\lambda^c}{\partial t}(t_k,b_j)
\approx
\frac{V_j^k - V_j^{k+1}}{\delta t},
\quad
\frac{\partial^2 V_\lambda^c}{\partial B_\lambda^2}(t_k,b_j)
\approx
\frac{V_{j+1}^k - 2V_j^k + V_{j-1}^k}{\delta B^2}.
$$

Terminal condition at $$t_0=T$$:
$$
V_j^0
=
H\!\left(\sum_{i=1}^n w_i S_i(0)\exp\!\left( (r - \tfrac{1}{2}\sigma_i^2)T + \sigma_i b_j \right)\right).
$$

Update rule for interior nodes (paper’s equation (24)):
$$
V_j^{k+1}
=
\frac{1}{2}\frac{\delta t}{\delta B^2}V_{j-1}^k
+
\left(1-r\delta t - \frac{\delta t}{\delta B^2}\right)V_j^k
+
\frac{1}{2}\frac{\delta t}{\delta B^2}V_{j+1}^k,
\quad j=1,\dots,J-1.
$$
Boundary values $$V_0^k$$ and $$V_J^k$$ are contract-dependent (the paper references standard option-PDE boundary condition design).

Stability condition (Theorem 4.1):
$$
\delta t \leq \frac{2\delta B^2}{r\delta B^2 + 2}.
$$

### 3.8 American early exercise via nodewise maximization (Section 4.3)

For an American put, the scheme introduces early exercise at each time step:

1. Compute continuation value $$\tilde V_j^{k+1}$$ by the European-like explicit update.
2. Compute intrinsic value at the node via the basket mapping from $$b_j$$ to $$S_j^{k+1}$$ (in the lower-bound market, volatilities become $$\nu_i\sigma_i$$):
   $$
   S_j^{k+1}
   =
   \sum_{i=1}^n w_i S_i(0)
   \exp\!\left( (r - \tfrac{1}{2}(\nu_i\sigma_i)^2)t_{k+1} + \nu_i\sigma_i b_j \right).
   $$
3. Set the American value by
   $$
   V_j^{k+1} = \max\!\left(\tilde V_j^{k+1},\; (K - S_j^{k+1})_+ \right).
   $$

The exercise boundary $$S^\ast(t_k)$$ is extracted from the grid as the basket level below which immediate exercise dominates continuation.

### 3.9 Delta extraction in the comonotonic lower-bound grid (Section 4.2)

The paper computes deltas with respect to the basket spot $$S^l$$ from the grid by applying the chain rule. With nodewise approximation
$$
\Delta_j^k \approx
\frac{V_{j+1}^k - V_{j-1}^k}{2\delta B \sum_{i=1}^n \nu_i\sigma_i w_i (S_i)_j^k},
$$
where
$$
(S_i)_j^k = S_i(0)\exp\!\left( (r - \tfrac{1}{2}(\nu_i\sigma_i)^2)t_k + \nu_i\sigma_i b_j \right).
$$

### 3.10 Approximation for correlated baskets via upper/lower bound combination (Section 5)

For convex payoff (paper develops it for basket puts), define a combined approximation for a European basket option:
$$
\bar V(t,S_1,\dots,S_n;K)
=
z\,V^l(t,S^l;K) + (1-z)\,V^c(t,S^c;K),
$$
with weight
$$
z
=
\frac{\mathrm{Var}[S^c(T)\mid \mathcal{F}_t] - \mathrm{Var}[S(T)\mid \mathcal{F}_t]}
{\mathrm{Var}[S^c(T)\mid \mathcal{F}_t] - \mathrm{Var}[S^l(T)\mid \mathcal{F}_t]}.
$$

The conditional variances under $$Q$$ are (paper’s formulas):
$$
\mathrm{Var}[S(T)\mid \mathcal{F}_t]
=
\sum_{i=1}^n\sum_{j=1}^n
w_i w_j S_i(t)S_j(t)e^{2r(T-t)}
\left(e^{\rho_{i,j}\sigma_i\sigma_j(T-t)}-1\right),
$$
$$
\mathrm{Var}[S^c(T)\mid \mathcal{F}_t]
=
\sum_{i=1}^n\sum_{j=1}^n
w_i w_j S_i(t)S_j(t)e^{2r(T-t)}
\left(e^{\sigma_i\sigma_j(T-t)}-1\right),
$$
$$
\mathrm{Var}[S^l(T)\mid \mathcal{F}_t]
=
\sum_{i=1}^n\sum_{j=1}^n
w_i w_j S_i(t)S_j(t)e^{2r(T-t)}
\left(e^{\nu_i\nu_j\sigma_i\sigma_j(T-t)}-1\right).
$$

The paper states that choosing $$z$$ as above enforces the moment identity
$$
\int_0^{+\infty} \bar V\,dK = \int_0^{+\infty} V\,dK
$$
for European options. The paper then applies the same bound-computation machinery (finite differences with early exercise for each bound) in the American setting and compares against LSM.

### 3.11 End-to-end pipeline (paper’s architecture) — ASCII overview

```
┌──────────────────────────────────────────────────────────────────────────┐
│ Inputs                                                                    │
│  - Spots: S_i(0), vols: σ_i, weights: w_i, corr: ρ_{i,j}, rate: r, T, K   │
│  - Payoff H(·) (convex in paper’s main tests: put)                         │
└───────────────────────────────┬──────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Comonotonic bounds (Section 2)                                            │
│  Upper bound: comonotonic market with vols σ_i                             │
│  Lower bound: comonotonic market with adjusted vols ν_i σ_i                │
└───────────────────────────────┬──────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ 2D PDE in (t, B_λ) with constant coefficients (Section 3)                  │
│  ∂_t V + (1/2) ∂^2_{B_λ} V - rV = 0                                        │
│  Terminal condition encodes basket mapping and payoff                      │
└───────────────────────────────┬──────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Finite difference solver (Section 4)                                      │
│  - Explicit update (24), stable if (25)                                    │
│  - American: V := max(continuation, intrinsic) per time step               │
└───────────────────────────────┬──────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Outputs                                                                   │
│  - Prices: V^l, V^c, and approximation \bar V (Section 5)                  │
│  - Deltas from grid (Section 4.2)                                         │
│  - Exercise boundary S*(t) from max-step (Section 4.3)                     │
└──────────────────────────────────────────────────────────────────────────┘
```

### 3.12 Step-by-step procedures (all named algorithms)

#### Algorithm 1 — Lower-bound volatility factors $$\nu_i$$ used in the paper’s numerics (equation (26))

**Inputs:** $$n$$, weights $$w_i$$, initial spots $$S_i(0)$$, volatilities $$\sigma_i$$, correlation matrix $$\rho_{i,j}$$.  
**Outputs:** $$\nu_i$$ for $$i=1,\dots,n$$ (used as effective volatilities $$\nu_i\sigma_i$$ in the comonotonic lower-bound market).

1. Compute the denominator:
   $$
   D
   =
   \sqrt{
   \sum_{j=1}^n\sum_{k=1}^n
   w_j w_k S_j(0) S_k(0)\rho_{j,k}\sigma_j\sigma_k
   }.
   $$
2. For each $$i \in \{1,\dots,n\}$$, compute
   $$
   \nu_i
   =
   \frac{\sum_{j=1}^n w_j S_j(0)\rho_{i,j}\sigma_j}{D}.
   $$
3. Use the resulting effective marginal volatilities $$\nu_i\sigma_i$$ in the comonotonic SDE (lower-bound market).

#### Algorithm 2 — Explicit finite difference solver for the comonotonic European PDE (Section 4.1)

**Inputs:** payoff $$H$$, $$r$$, $$T$$, stock parameters $$\{S_i(0),\sigma_i,w_i\}_{i=1}^n$$, grid parameters $$\delta t,\delta B,N,J,I$$, boundary conditions $$V_0^k,V_J^k$$.  
**Outputs:** grid values $$V_j^k \approx V_\lambda^c(t_k,b_j)$$, in particular price at time $$t_N=0$$ and the relevant $$b_{j_0}$$ corresponding to $$B_\lambda(0)=0$$.

1. Construct backward time grid $$t_k=T-k\delta t$$, $$k=0,\dots,N$$.
2. Construct space grid $$b_j=(j-I)\delta B$$, $$j=0,\dots,J$$.
3. Initialize terminal vector at maturity:
   $$
   V_j^0 = H\!\left(\sum_{i=1}^n w_i S_i(0)\exp\!\left((r-\tfrac{1}{2}\sigma_i^2)T+\sigma_i b_j\right)\right).
   $$
4. For $$k=0,\dots,N-1$$:
   - Impose boundary conditions $$V_0^k,V_J^k$$.
   - For each interior node $$j=1,\dots,J-1$$, update using
     $$
     V_j^{k+1}
     =
     \frac{1}{2}\frac{\delta t}{\delta B^2}V_{j-1}^k
     +
     \left(1-r\delta t-\frac{\delta t}{\delta B^2}\right)V_j^k
     +
     \frac{1}{2}\frac{\delta t}{\delta B^2}V_{j+1}^k.
     $$
5. Return $$V^{N}$$ (time $$0$$) and interpolate/extract the node corresponding to $$B_\lambda(0)=0$$ (typically $$b_I=0$$).

**Termination / convergence control:** choose $$\delta t,\delta B$$ such that stability condition (25) holds; increase $$N,J$$ until observed convergence (paper reports using 15,000 time steps per year in American computations).

#### Algorithm 3 — American early exercise embedding in the comonotonic finite difference method (Section 4.3)

**Inputs:** As Algorithm 2, plus strike $$K$$, American payoff (paper uses put intrinsic value), and effective volatilities (use $$\sigma_i$$ for the upper bound and $$\nu_i\sigma_i$$ for the lower bound).  
**Outputs:** American grid prices $$V_j^k$$ and an exercise region estimate $$S^\ast(t_k)$$.

1. Initialize at maturity:
   $$
   V_j^0 = H\!\left(S(T;b_j)\right)
   $$
   with the appropriate basket mapping $$S(T;b_j)$$ (upper or lower bound parameterization).
2. For each time step $$k=0,\dots,N-1$$:
   1. Compute continuation values $$\tilde V_j^{k+1}$$ using the explicit update rule (24).
   2. Compute nodewise basket levels $$S_j^{k+1}$$ from $$b_j$$ using the comonotonic basket mapping at time $$t_{k+1}$$.
   3. Overwrite with early exercise:
      $$
      V_j^{k+1} = \max\!\left(\tilde V_j^{k+1},\; (K-S_j^{k+1})_+\right).
      $$
3. Extract the exercise boundary $$S^\ast(t_k)$$ as the grid-implied threshold where $$V_j^k$$ switches from continuation to intrinsic.

#### Algorithm 4 — Approximate correlated-basket price via bound combination (Section 5)

**Inputs:** original basket parameters including correlations $$\rho_{i,j}$$, payoff (convex), and maturity/strike.  
**Outputs:** approximation $$\bar V$$.

1. Compute lower-bound factors $$\nu_i$$ (Algorithm 1).
2. Compute upper-bound price $$V^c$$ using Algorithm 3 with effective volatilities $$\sigma_i$$ and comonotonic dependence.
3. Compute lower-bound price $$V^l$$ using Algorithm 3 with effective volatilities $$\nu_i\sigma_i$$ and comonotonic dependence.
4. Compute $$z$$ from conditional variances (28).
5. Output
   $$
   \bar V = z V^l + (1-z)V^c.
   $$

## 4. Theoretical Results

### Lemma 2.1 (Section 2.4): sufficient conditions for a comonotonic lower bound

**Statement (full claim).**  
Assume the factors $$\nu_i$$ in the comonotonic lower-bound stock model
$$
dS_i^l(t) = \mu_i S_i^l(t)\,dt + \nu_i\sigma_i S_i^l(t)\,dB(t)
$$
are defined by
$$
\nu_i
=
\frac{\sum_{j=1}^n \beta_j\rho_{i,j}\sigma_j}
{\sqrt{\sum_{j=1}^n\sum_{k=1}^n \beta_j\beta_k\rho_{j,k}\sigma_j\sigma_k}},
\quad i=1,\dots,n,
\quad \beta_j>0.
$$
Then there exists a choice of $$\beta_j>0$$ such that $$0 \leq \nu_i \leq 1$$ for all $$i$$ and the corresponding comonotonic lower-bound derivative price satisfies
$$
V^l(t,S^l) \leq V(t,S_1,\dots,S_n)
$$
in the setting of Section 2 (paper’s context is convex payoff).

**Proof strategy sketch.**  
The lemma is a convex-order/comonotonicity bound: by reducing marginal volatilities appropriately and imposing comonotonic dependence, the basket distribution is made “less variable” in a sense compatible with convex ordering, so expectations of convex payoffs decrease. The parameters $$\beta_j$$ are tuning degrees of freedom ensuring the resulting scaled volatilities are admissible $$([0,1])$$ while preserving required ordering relations. The paper defers the detailed constructive choice of $$\beta_j$$ and full proof to prior work (Kaas et al. (2000), Deelstra et al. (2004), Linders and Stassen (2016)).

### Lemma 3.1 (Section 3.1): Ito/SIE dynamics of the comonotonic basket

**Statement (full claim).**  
In the comonotonic market of (4), the basket $$S^c(t)=\sum_{i=1}^n w_i S_i^c(t)$$ satisfies
$$
S^c(t) = S(0) + \int_0^t \mu^c(s)\,ds + \int_0^t \sigma^c(s)\,dB(s),
$$
where
$$
\mu^c(s) = \sum_{i=1}^n \mu_i w_i S_i^c(s),
\quad
\sigma^c(s) = \sum_{i=1}^n \sigma_i w_i S_i^c(s).
$$
Equivalently,
$$
dS^c(t) = \mu^c(t)\,dt + \sigma^c(t)\,dB(t).
$$

**Proof strategy sketch.**  
Linearity of stochastic integrals and finite sums yields the result directly: substitute $$S^c(t)=\sum_i w_i S_i^c(t)$$ into the SIEs for $$S_i^c$$ and interchange summation with integration.

### Theorem 3.1 (Section 3.2.1; Appendix 1): PDE for the European comonotonic price and unique replicating strategy

**Statement (full claim).**  
For a European-type derivative with payoff
$$
V_\lambda^c(T,B_\lambda)
=
H\!\left(
\sum_{i=1}^n w_i S_i(0)\exp\!\left((r-\tfrac{1}{2}\sigma_i^2)T+\sigma_i B_\lambda\right)
\right),
$$
the price function $$V_\lambda^c(t,B_\lambda)$$ satisfies the backward PDE
$$
\frac{\partial V_\lambda^c}{\partial t}
+\frac{1}{2}\frac{\partial^2 V_\lambda^c}{\partial B_\lambda^2}
-rV_\lambda^c
=0.
$$
Moreover, there exists a unique self-financing replicating strategy $$\pi(t)=(\pi_1(t),\pi_2(t))$$ in the comonotonic basket and bank account, with
$$
\pi_1(t)=\frac{1}{\sigma^c(t)}\frac{\partial V_\lambda^c}{\partial B_\lambda}(t,B_\lambda),
\quad
\pi_2(t)=\frac{1}{D(t)}\left(V_\lambda^c(t,B_\lambda)-\pi_1(t)S^c(t)\right).
$$

**Proof strategy sketch.**  
The appendix follows the classical Black–Scholes replication argument adapted to the state variable $$B_\lambda$$. Apply Ito’s formula to $$V_\lambda^c(t,B_\lambda)$$ and write the portfolio value process for a self-financing strategy in $$S^c$$ and $$D$$. Choose $$\pi_1$$ to match the diffusion term (by equating the coefficient of $$dB_\lambda$$), which forces the hedging error’s stochastic integral to be zero and yields uniqueness. Enforcing zero drift for the hedging error yields the PDE.

### Theorem 3.2 (Section 3.3; Appendix 2): explicit integral solution to the comonotonic PDE

**Statement (full claim).**  
The European comonotonic price admits the explicit representation
$$
V_\lambda^c(t,B_\lambda)
=
e^{-r(T-t)}
\int_{-\infty}^{+\infty}
H\!\left(
\sum_{i=1}^n w_i S_i^c(t)\exp\!\left((r-\tfrac{1}{2}\sigma_i^2)(T-t)+\sigma_i y\right)
\right)\,
\phi_{T-t}(y)\,dy,
$$
where
$$
S_i^c(t)=S_i(0)\exp\!\left((r-\tfrac{1}{2}\sigma_i^2)t+\sigma_i B_\lambda(t)\right),
$$
and $$\phi_{T-t}$$ is the normal density with mean $$0$$ and variance $$T-t$$:
$$
\phi_{T-t}(y)=\frac{\exp\!\left(-\frac{y^2}{2(T-t)}\right)}{\sqrt{2\pi (T-t)}}.
$$

**Proof strategy sketch.**  
Appendix 2 transforms the PDE into the 1D heat equation via (i) discounting the value function and (ii) time reversal/affine rescaling of time. The heat equation solution is written as convolution with the Gaussian fundamental solution, matching the terminal condition at $$t=T$$. Substituting back yields a discounted Gaussian expectation over future Brownian increments.

### Theorem 4.1 (Section 4.1; Appendix 3): stability of the explicit finite difference scheme

**Statement (full claim).**  
The explicit finite difference scheme
$$
V_j^{k+1}
=
\frac{1}{2}\frac{\delta t}{\delta B^2}V_{j-1}^k
+
\left(1-r\delta t-\frac{\delta t}{\delta B^2}\right)V_j^k
+
\frac{1}{2}\frac{\delta t}{\delta B^2}V_{j+1}^k
$$
is stable provided
$$
\delta t \leq \frac{2\delta B^2}{r\delta B^2+2}.
$$

**Proof strategy sketch.**  
Appendix 3 performs a von Neumann/Fourier stability analysis: introduce a perturbed solution with oscillatory error $$\epsilon_j^k=\nu^k e^{ij\omega}$$ and derive the amplification factor $$\nu$$ as a function of $$\delta t,\delta B,r,\omega$$. Stability requires $$|\nu|\leq 1$$ for all frequencies; bounding the amplification factor yields the condition.

### Lemma A.1 (Appendix 4, A.2–A.3): CDF inversion characterization for $$S^c(T)$$

**Statement (full claim).**  
For any $$x>0$$, there exists $$F_{S^c}(x;T,t)\in (0,1)$$ satisfying
$$
\sum_{i=1}^n
w_i S_i^c(t)
\exp\!\left((r-\tfrac{1}{2}\sigma_i^2)(T-t)
+\sigma_i\sqrt{T-t}\,\Phi^{-1}(F_{S^c}(x;T,t))\right)
= x,
$$
where $$\Phi$$ is the standard normal CDF. Moreover,
$$
F_{S^c}(x;T,t)=Q\!\left(S^c(T)\leq x \mid \mathcal{F}_t\right),
$$
the time-$$t$$ risk-neutral CDF of the comonotonic basket evaluated at $$x$$.

**Proof strategy sketch.**  
Under $$Q$$, the Brownian increment $$B_\lambda(T)-B_\lambda(t)$$ is Gaussian, so it can be represented as $$\sqrt{T-t}\,\Phi^{-1}(U)$$ with uniform $$U$$. The comonotonic basket is a strictly increasing function of that increment, so its CDF is invertible. Additivity of inverse CDFs for comonotonic sums yields an explicit expression for $$F_{S^c}^{-1}$$; composing with $$F_{S^c}$$ gives the stated identity.

### Theorem A.1 (Appendix 4, A.2–A.4): risk-neutral density for $$S^c(T)$$ and pricing integral over basket states

**Statement (full claim).**  
The European comonotonic price can be written as
$$
V^c(t,S^c)
=
e^{-r(T-t)}\int_0^{+\infty} H(s)\,f_{S^c}(s;T,t)\,ds,
$$
where $$f_{S^c}(s;T,t)$$ is the time-$$t$$ risk-neutral density of $$S^c(T)$$ given by
$$
f_{S^c}(s;T,t)
=
\frac{
\phi_{T-t}\!\left(\sqrt{T-t}\,\Phi^{-1}(F_{S^c}(s;T,t))\right)
}{
\sum_{i=1}^n
\sigma_i w_i S_i^c(t)
\exp\!\left((r-\tfrac{1}{2}\sigma_i^2)(T-t)
+\sigma_i\sqrt{T-t}\,\Phi^{-1}(F_{S^c}(s;T,t))\right)
}.
$$

**Proof strategy sketch.**  
Start from the Gaussian-integral representation in Theorem 3.2 and change variables from the Brownian increment coordinate to the future basket level $$s$$. Lemma A.1 provides the inverse mapping and its derivative, yielding the Jacobian required to express the integral as an expectation over $$s$$. The appendix verifies density properties (positivity and unit integral) and shows that integrating the density up to $$x$$ recovers $$F_{S^c}(x;T,t)$$.

### Additional result (Appendix 4, A.1): alternative PDE in coordinates $$ (t,S^c) $$

The paper notes the comonotonic PDE can also be rewritten as
$$
\frac{\partial V^c}{\partial t}
+\frac{1}{2}(\sigma^c(t))^2\frac{\partial^2 V^c}{\partial (S^c)^2}
+rS^c\frac{\partial V^c}{\partial S^c}
-rV^c
=0,
$$
with final condition $$V^c(T,S^c)=H(S^c(T))$$ and boundary conditions derived from limits $$S^c \to 0$$ and $$S^c \to +\infty$$. The paper reports that solving the constant-coefficient PDE in $$ (t,B_\lambda) $$ is faster in practice than solving this time-dependent-coefficient form.

## 5. Experimental Evaluation

### 5.1 Experimental setups, baselines, metrics, and hyperparameters (Sections 4.2–5)

| Item | Setting(s) used in the paper | Notes / parameters |
|---|---|---|
| Underlyings | Synthetic baskets (examples with $$n=4$$ and $$n=8$$; also varying $$n$$ up to 20 in timing test) | No historical dataset; all tests are parameterized Black–Scholes baskets |
| Payoffs | Basket puts (convex payoff) | Early exercise only for American-type |
| Baselines | (i) Monte Carlo simulation (European “reference” in Section 4.2); (ii) LSM (Longstaff–Schwartz) for American prices | LSM used as standard benchmark |
| PDE method | Explicit finite differences on the 2D comonotonic PDE in $$ (t,B_\lambda) $$ | Upper/lower comonotonic bounds; lower bound used as close approximation in Section 4; combined bound approximation in Section 5 |
| FD time discretization | 15,000 time steps per year | Paper reports convergence tests indicating sufficiency |
| FD stability | Choose $$\delta B$$ to satisfy stability condition $$\delta t \leq \frac{2\delta B^2}{r\delta B^2+2}$$ | Explicit scheme stability |
| LSM exercise grid | 50 exercise times per year | Bermudan-style discretization of American exercise for simulation |
| LSM regression basis | Quadratic polynomials in normalized prices | Regression model: $$\beta_0 + \sum_i \beta_{1,i}\frac{S_i(t)}{S_i(0)} + \sum_i \beta_{2,i}\left(\frac{S_i(t)}{S_i(0)}\right)^2$$ |
| LSM simulation size | 100,000 paths; reported LSM price is average of 10 runs | Intended to reduce Monte Carlo variance |
| Hardware / software | R on a Dell machine with CPU “Intel(R) Core(TM) i7-4910MQ CPU 2.90GHz” | Random normals via `mvrnorm()` (MASS); regression via `lm()` |

### 5.2 European 4-asset basket put: parameters (Table 1) and bound behavior (Section 4.2)

**Table 1 (paper).** Four-asset basket put option parameters:

| Parameter | Stock 1 | Stock 2 | Stock 3 | Stock 4 | Basket-level |
|---|---:|---:|---:|---:|---|
| Volatility $$\sigma_i$$ | 0.25 | 0.30 | 0.45 | 0.80 | Pairwise correlation = 0.3 |
| Spot $$S_i(0)$$ | 100 | 100 | 100 | 100 | Maturity $$T=1$$ |
| Weight $$w_i$$ | 0.25 | 0.25 | 0.25 | 0.25 | Risk-free rate $$r=0.01$$ |

Section 4.2 computes:
- $$V^{\mathrm{Sim}}$$ via Monte Carlo with 1,000,000 simulations (treated as reference),
- comonotonic upper bound $$V^c$$ and lower bound $$V^l$$ via the finite difference scheme, with lower-bound factors $$\nu_i$$ chosen as
  $$
  \nu_i
  =
  \frac{\sum_{j=1}^n w_j S_j(0)\rho_{i,j}\sigma_j}
  {\sqrt{\sum_{j=1}^n\sum_{k=1}^n w_j w_k S_j(0)S_k(0)\rho_{j,k}\sigma_j\sigma_k}},
  \quad i=1,\dots,n.
  $$

**Reported qualitative findings (Figure 1):**
- For each strike $$K \in \{95,96,\dots,115\}$$, the put’s convex payoff implies bounds $$V^l \leq V^{\mathrm{Sim}} \leq V^c$$.
- The lower bound is reported to be a close approximation $$V^l \approx V^{\mathrm{Sim}}$$ for these parameters.
- The upper bound is reported to be less accurate because it removes dependence information by forcing correlations to 1.
- Computation time: Monte Carlo about 4 seconds; bound computation about 0.4 seconds (Section 4.2).

### 5.3 Price and delta surfaces from the lower-bound FD grid (Section 4.2)

For the four-asset case with $$K=100$$, the paper uses the full FD grid to produce:
- a price surface across spot values and times to maturity (Figure 2, left),
- the corresponding delta surface (Figure 2, right),
with deltas computed from grid values via the chain-rule approximation in Section 4.2.

### 5.4 American lower-bound FD scheme and exercise region (Section 4.3)

The paper implements early exercise by applying the nodewise max operator
$$
V_j^{k+1}=\max\!\left(\tilde V_j^{k+1},(K-S_j^{k+1})_+\right)
$$
at each time step (explicit scheme continuation value $$\tilde V$$ plus intrinsic payoff). Using 15,000 time steps per year and a stability-compliant $$\delta B$$ grid, the paper extracts an exercise boundary $$S^\ast(t)$$ (Figure 3) for the four-asset American put.

### 5.5 American 8-asset basket put: FD vs LSM comparison (Table 2; Section 5)

Setup described in Section 5:
- Basket: 8 equally weighted stocks, initial $$S_i(0)=40$$ for all.
- Volatilities: stocks 2–8 have volatilities 60%, 10%, 90%, 30%, 70%, 80%, 20%; stock 1 volatility is either 30% or 90%.
- Pairwise correlation uniform: $$\rho \in \{0.3,0.8\}$$.
- Risk-free rate: 1%.
- Maturities: 6 months, 1 year, 2 years.
- Strikes: 35, 40, 45.
- LSM: 50 exercise times per year; 100,000 paths; quadratic regression in normalized stock prices; average of 10 runs.
- FD method: comonotonic-PDE finite difference machinery (upper/lower bound computations and approximation as described).

**Table 2 (paper): American-type basket put option prices and computation times.**

| Maturity | Strike | $$\sigma_1$$ | $$\rho$$ | FD price | LSM price | Rel. diff. (%) | Comp. time FD (s) | Comp. time LSM (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 months | 35 | 0.3 | 0.3 | 1.279 | 1.282 | 0.219 | 9.36 | 79.81 |
| 6 months | 35 | 0.3 | 0.8 | 2.408 | 2.413 | 0.228 | 9.17 | 96.44 |
| 6 months | 35 | 0.9 | 0.3 | 1.700 | 1.703 | 0.181 | 9.17 | 86.71 |
| 6 months | 35 | 0.9 | 0.8 | 3.053 | 3.064 | 0.369 | 9.17 | 108.04 |
| 6 months | 40 | 0.3 | 0.3 | 3.461 | 3.468 | 0.183 | 9.36 | 164.88 |
| 6 months | 40 | 0.3 | 0.8 | 4.866 | 4.879 | 0.256 | 10.33 | 175.81 |
| 6 months | 40 | 0.9 | 0.3 | 3.989 | 4.008 | 0.453 | 10.39 | 163.91 |
| 6 months | 40 | 0.9 | 0.8 | 5.618 | 5.620 | 0.023 | 9.57 | 157.50 |
| 6 months | 45 | 0.3 | 0.3 | 6.742 | 6.739 | 0.037 | 9.33 | 216.22 |
| 6 months | 45 | 0.3 | 0.8 | 8.109 | 8.118 | 0.104 | 9.49 | 201.32 |
| 6 months | 45 | 0.9 | 0.3 | 7.236 | 7.245 | 0.132 | 9.33 | 203.58 |
| 6 months | 45 | 0.9 | 0.8 | 8.860 | 8.881 | 0.232 | 8.98 | 196.07 |
| 1 year | 35 | 0.3 | 0.3 | 2.386 | 2.374 | 0.486 | 18.78 | 231.91 |
| 1 year | 35 | 0.3 | 0.8 | 4.015 | 4.019 | 0.100 | 18.34 | 270.99 |
| 1 year | 35 | 0.9 | 0.3 | 3.033 | 3.032 | 0.030 | 18.09 | 246.30 |
| 1 year | 35 | 0.9 | 0.8 | 4.945 | 4.946 | 0.022 | 18.22 | 280.50 |
| 1 year | 40 | 0.3 | 0.3 | 4.848 | 4.837 | 0.237 | 17.97 | 1601.22 |
| 1 year | 40 | 0.3 | 0.8 | 6.743 | 6.732 | 0.159 | 84.94 | 1557.22 |
| 1 year | 40 | 0.9 | 0.3 | 5.596 | 5.596 | 0.010 | 85.59 | 1515.75 |
| 1 year | 40 | 0.9 | 0.8 | 7.784 | 7.787 | 0.035 | 81.05 | 1541.76 |
| 1 year | 45 | 0.3 | 0.3 | 8.093 | 8.079 | 0.169 | 84.80 | 1977.86 |
| 1 year | 45 | 0.3 | 0.8 | 10.027 | 10.016 | 0.109 | 82.61 | 1281.59 |
| 1 year | 45 | 0.9 | 0.3 | 8.093 | 8.842 | 0.077 | 18.24 | 474.82 |
| 1 year | 45 | 0.9 | 0.8 | 11.098 | 11.085 | 0.120 | 17.85 | 455.61 |
| 2 years | 35 | 0.3 | 0.3 | 4.035 | 3.974 | 1.521 | 37.36 | 729.34 |
| 2 years | 35 | 0.3 | 0.8 | 6.196 | 6.183 | 0.216 | 37.30 | 2374.89 |
| 2 years | 35 | 0.9 | 0.3 | 4.978 | 4.954 | 0.492 | 167.00 | 3028.90 |
| 2 years | 35 | 0.9 | 0.8 | 7.481 | 7.477 | 0.065 | 163.18 | 3454.50 |
| 2 years | 40 | 0.3 | 0.3 | 6.775 | 6.704 | 1.041 | 169.78 | 4023.14 |
| 2 years | 40 | 0.3 | 0.8 | 9.204 | 9.176 | 0.299 | 165.61 | 4148.22 |
| 2 years | 40 | 0.9 | 0.3 | 7.822 | 7.814 | 0.101 | 167.06 | 4055.97 |
| 2 years | 40 | 0.9 | 0.8 | 10.614 | 10.594 | 0.187 | 163.75 | 4245.87 |
| 2 years | 45 | 0.3 | 0.3 | 10.065 | 10.011 | 0.538 | 170.83 | 4990.59 |
| 2 years | 45 | 0.3 | 0.8 | 12.602 | 12.576 | 0.209 | 165.97 | 4862.00 |
| 2 years | 45 | 0.9 | 0.3 | 11.136 | 11.130 | 0.049 | 167.44 | 4911.06 |
| 2 years | 45 | 0.9 | 0.8 | 14.075 | 14.049 | 0.188 | 164.21 | 4836.68 |

**Ambiguity flag (faithfulness vs transcription).**  
The two rows for maturity “1 year”, strike “45”, $$\sigma_1=0.9$$ show FD/LSM computation times (18.24, 474.82) and (17.85, 455.61) that are materially smaller than nearby 1-year rows in the same table; the values above are reproduced exactly as present in the provided paper text, but these two rows are candidates for typographical or transcription misalignment.

**Reported conclusions from Table 2 (Section 5):**
- Prices from FD and LSM are close, with relative differences generally below about 1.5% in the table.
- Computation times favor the FD method strongly; the paper reports speedups up to about 30 times versus LSM in tested settings.
- LSM time increases with maturity (more exercise opportunities) and strike (more in-the-money paths and regressions), while the comonotonic-PDE machinery remains comparatively fast.

### 5.6 Sensitivity study: relative differences vs maturity/strike/volatility/correlation (Figure 4; Section 5)

Basket: 4 equally weighted stocks, initial prices 40. Base case: maturity 1 year, strike 40, pairwise correlation 0.3, volatilities $$\sigma_1=40\%$$, $$\sigma_2=60\%$$, $$\sigma_3=10\%$$, $$\sigma_4=90\%$$.

Reported observations:
- Relative differences are “relatively small” overall but increase for high maturities (paper suggests compounded approximation error and potential bias in LSM).
- Relative differences can become large for low strikes (paper attributes to fewer in-the-money paths and weaker regressions in LSM).
- Relative differences increase for high volatilities (same mechanism: regressions degrade when ITM sample size changes).
- Differences are often very small for moderate correlation values, but extreme correlations close to 1 reveal major LSM instability (next subsection).

### 5.7 LSM deficiency near the comonotonic limit and an adjusted LSM fix (Tables 3–5; Section 5)

**Table 3 (paper): high pairwise correlations (4-stock equally weighted basket; initial prices 40).**

| Pairwise correlation $$\rho$$ | FD price | LSM price | Comonotonic LSM |
|---:|---:|---:|---:|
| 0.95 | 7.427 | 7.433 | – |
| 0.96 | 7.452 | 7.444 | – |
| 0.97 | 7.476 | 7.476 | – |
| 0.98 | 7.500 | 7.511 | – |
| 0.99 | 7.524 | 8.670 | – |
| 0.992 | 7.529 | 10.250 | – |
| 0.994 | 7.534 | 10.234 | – |
| 0.996 | 7.538 | 14.702 | – |
| 0.998 | 7.543 | 17.212 | – |
| 1.000 | 7.548 | 259.274 | 7.550 |

Reported interpretation (Section 5):
- FD prices remain smooth as $$\rho \to 1$$.
- Standard LSM becomes unstable as correlations approach 1 and explodes in the comonotonic case.
- Proposed fix: in the comonotonic case, regress continuation value on only one asset (the first), avoiding overfitting from redundant state variables; this yields a comonotonic-LSM price aligned with FD.

**Table 4 (paper): regression coefficients and p-values at exercise time 25 for $$\rho \in \{0.5,0.95,1\}$$.**

| Coefficient | $$\rho=0.5$$ | $$\rho=0.95$$ | $$\rho=1$$ |
|---|---|---|---|
| $$\beta_0$$ | 56.36 (p < 2e-16) | 47.85 (p = 4e-05) | 9383.1 (p = 0.508) |
| $$\beta_{1,1}$$ | -12.17 (p < 2e-16) | -14.86 (p = 4e-09) | 101543.8 (p = 0.521) |
| $$\beta_{1,2}$$ | -11.43 (p < 2e-16) | -12.82 (p < 2e-16) | -268656.9 (p = 0.530) |
| $$\beta_{1,3}$$ | -46.10 (p = 7e-09) | -22.96 (p = 0.369) | -24649.0 (p = 0.511) |
| $$\beta_{1,4}$$ | -9.83 (p < 2e-16) | -13.78 (p < 2e-16) | -351436.8 (p = 0.547) |
| $$\beta_{2,1}$$ | 2.44 (p = 4e-10) | 4.79 (p = 0.0009) | 460461.3 (p = 0.541) |
| $$\beta_{2,2}$$ | 2.42 (p < 2e-16) | 3.90 (p = 2e-05) | 28192.3 (p = 0.569) |
| $$\beta_{2,3}$$ | 19.54 (p = 1.5e-06) | 7.44 (p = 0.575) | Not available |
| $$\beta_{2,4}$$ | 2.27 (p < 2e-16) | 7.74 (p < 2e-16) | -909.1 (p = 0.624) |

**Table 5 (paper): adjusted comonotonic LSM regression at exercise time 25 (only stock 1 regressors).**

| Coefficient | $$\rho=1$$ |
|---|---|
| $$\beta_0$$ | 42.54 (p < 2e-16) |
| $$\beta_{1,1}$$ | -48.52 (p < 2e-16) |
| $$\beta_{2,1}$$ | 11.78 (p < 2e-16) |

### 5.8 Computation time vs basket size (Figure 5; Section 5)

Setting: maturity 1 year, strike 40, all spots 40, all volatilities 30%, pairwise correlation 0.3. Reported conclusion: FD time grows mildly with basket size (some overhead in mapping the $$B_\lambda$$ grid through all assets), while LSM time increases steeply due to simulating higher-dimensional Brownian drivers and regressions; paper reports speedups up to about 30 times for small baskets and about 10 times for larger baskets.

## 6. ASCII Architecture / Workflow Diagram(s)

### Diagram 1 — Full workflow: correlated basket to comonotonic-PDE approximation

```
┌──────────────────────────────────────────────────────────────────────────┐
│ Correlated basket model (n assets)                                         │
│  Inputs: {S_i(0), σ_i, w_i}_{i=1..n}, correlation matrix ρ, r, T, payoff H │
└───────────────────────────────┬──────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Comonotonic bound construction                                            │
│  Upper: set dependence to comonotonic, use vols σ_i                        │
│  Lower: compute ν_i (eq. 26) and use vols ν_i σ_i                          │
└───────────────────────────────┬──────────────────────────────────────────┘
                                │ (two 1-factor problems)
                                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ State reduction to (t, B_λ)                                                │
│  Basket mapping:                                                          │
│   S(t;B_λ)= Σ_i w_i S_i(0) exp((r-0.5 σ_i^2)t + σ_i B_λ)                    │
│  PDE: ∂_t V + 0.5 ∂^2_{B_λ}V - rV = 0                                      │
└───────────────────────────────┬──────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Explicit finite differences on (t,B_λ) grid                                 │
│  - Update: V^{k+1}_j from V^k_{j-1},V^k_j,V^k_{j+1} (eq. 24)                │
│  - Stability constraint on δt,δB (Thm 4.1)                                  │
│  - American: V := max(continuation, intrinsic) per time step                │
└───────────────────────────────┬──────────────────────────────────────────┘
                                │
                                ▼
┌──────────────────────────────────────────────────────────────────────────┐
│ Outputs                                                                    │
│  - Prices: V^l, V^c, optional combination \bar V via z (Section 5)          │
│  - Greeks: delta from grid + chain rule                                     │
│  - Exercise region: S*(t) extracted from max operator                        │
└──────────────────────────────────────────────────────────────────────────┘
```

### Diagram 2 — One FD time step for the American put (continuation + early exercise)

```
Time step: t_k  →  t_{k+1}   (backward in calendar time)

┌──────────────────────────────┐
│ Given grid slice V^k_j        │
│ at time t_k                  │
└───────────────┬──────────────┘
                │  explicit PDE update (eq. 24)
                ▼
┌──────────────────────────────┐
│ Continuation:  V~^{k+1}_j     │
│ from neighbors j-1,j,j+1      │
└───────────────┬──────────────┘
                │  compute intrinsic payoff
                ▼
┌──────────────────────────────┐
│ Intrinsic: (K - S^{k+1}_j)_+  │
│ where S^{k+1}_j = Σ_i w_i S_i │
│   exp((r-0.5(ν_iσ_i)^2)t_{k+1}│
│       + ν_iσ_i b_j)           │
└───────────────┬──────────────┘
                │  early exercise projection
                ▼
┌──────────────────────────────┐
│ American value: V^{k+1}_j     │
│ = max(V~^{k+1}_j, intrinsic)  │
└──────────────────────────────┘
```

## 7. Follow-Up Works & Extensions

### 7.1 PDE/PDCP-based extensions and comparisons for American basket options

**[in’t Hout & Snoeijer, Mathematics 2021]** Karel J. in’t Hout and Jacob Snoeijer, *Numerical Valuation of American Basket Options via Partial Differential Complementarity Problems*. The paper explicitly studies the comonotonic approach “considered by Hanbali and Linders (2019)” alongside PCA-based dimension reduction for American basket options formulated as PDCPs, where each approximation reduces valuation to a limited set of low-dimensional PDCP solves. It reports (via numerical experiments) that the PCA-based and comonotonic approximations “always yield approximations … that lie close to each other” for both European and American basket put options, and it develops an efficient discretization (finite differences with ADI time-stepping and an Ikonen–Toivanen technique for complementarity handling) with near second-order convergence observed empirically. This work functions as a direct methodological follow-up by repositioning the 2019 comonotonic reduction inside a PDCP framework and benchmarking it against PCA reduction. ([mdpi.com](https://www.mdpi.com/2227-7390/9/13/1498))

### 7.2 Machine-learning alternatives for high-dimensional optimal stopping (contextual to the paper’s LSM comparison)

**[Herrera et al., arXiv 2021 (rev. 2023)]** Calypso Herrera, Florian Krach, Pierre Ruyssen, and Josef Teichmann, *Optimal Stopping via Randomized Neural Networks*. The paper proposes randomized neural network constructions (random hidden layers, trained output layer) to approximate continuation values with linear regression-like training, targeting high-dimensional optimal stopping problems including American option pricing under multiple models. Relative to Hanbali–Linders (2019), this line of work tackles the same “high-dimensional American option pricing” bottleneck but replaces PDE reduction by data-driven regression architectures designed to remain computationally tractable as state dimension increases. The relevance to the source paper is direct at the level of the shared bottleneck (continuation value approximation for American-style exercise) and the diagnostic focus on speed/Greeks, but the modeling and numerical approach is different (random-feature neural regressors rather than comonotonic-PDE reduction). ([arxiv.org](https://arxiv.org/abs/2104.13669))

### 7.3 Dynamic-programming / low-dimensional numerical design frameworks citing comonotonic/PDE dimension-reduction

**[Ben-Abdellatif et al., Risks 2024]** Malek Ben-Abdellatif, Hatem Ben-Ameur, Rim Chérif, and Bruno Rémillard, *Dynamic Programming for Designing and Valuing Two-Dimensional Financial Derivatives*. The paper develops a dynamic programming framework that couples finite elements, parallel computing, and local polynomial approximations in backward recursion, emphasizing two-dimensional derivative valuation and illustrating with American options. It cites Hanbali and Linders (2019) as an example of dimension reduction for option pricing problems and positions its own contribution as a modular two-component evaluation pipeline (underlying dynamics vs derivative structure). This is a follow-on in the sense that it re-frames “low-dimensional surrogate evaluation” as a design principle and explicitly connects to the dimension-reduction theme that motivates the 2019 comonotonic PDE approach. ([mdpi.com](https://www.mdpi.com/2227-9091/12/12/183))

### 7.4 Analytical/free-boundary approaches for multidimensional American options (citing the 2019 paper)

**[Agliardi & Agliardi, IJFS 2023]** Elettra Agliardi and Rossella Agliardi, *Pricing Multidimensional American Options*. The paper provides an explicit-form solution approach for optimal stopping problems involving a multidimensional geometric Brownian motion, using a free-boundary value method and fundamental-solution techniques, with emphasis on perpetual American-style contracts. It cites Hanbali and Linders (2019) among the referenced American-option valuation literature and addresses the same broad challenge category (multidimensional early-exercise valuation), but focuses on analytical structure for perpetual options rather than numerical dimension reduction for finite-maturity basket options. ([mdpi.com](https://www.mdpi.com/2227-7072/11/1/51))

## 8. Industrial & Real-World Applications

### 8.1 Open-source production-grade libraries implementing American basket option valuation (industry-standard baseline methods)

**[GitHub: lballabio/QuantLib]** QuantLib is a widely used open-source C++ quantitative finance library (about 6.8k GitHub stars at time of writing) that includes Monte Carlo engines for early-exercise products and explicitly exposes an “American basket-option engine” factory (e.g., `MakeMCAmericanBasketEngine`) built on a Longstaff–Schwartz-type Monte Carlo engine. This constitutes a real-world implementation path for the baseline methodology used for comparison in Hanbali–Linders (2019), namely LSM for multifactor American-style derivatives, including basket structures. No verified evidence was found (in public repositories) that QuantLib implements the specific comonotonic two-dimensional PDE approximation of Hanbali–Linders; the practical overlap is that both target feasible American basket valuation with Greeks/exercise handling, but via different numerical paradigms. ([github.com](https://github.com/lballabio/QuantLib))

### 8.2 Open-source risk engines embedding pricing libraries (deployment-oriented tooling)

**[GitHub: OpenSourceRisk/Engine]** Open Source Risk Engine (ORE) is an open-source framework for pricing and risk analysis (about 671 GitHub stars at time of writing) that is explicitly “based on QuantLib” and extends it with additional simulation models, instruments, and pricing engines, with interfaces (API/XML) and application launchers (Excel, Python, Jupyter). This provides an example of how American-style valuation and risk analytics are embedded into end-to-end risk systems, aligning with Hanbali–Linders’ emphasis on fast computation and the extraction of exercise/Greek information in practical workflows. No verified evidence was found (in the public ORE repository) of a built-in implementation of the specific comonotonic-PDE method; ORE is relevant as an industrial-grade integration layer where such a method could be deployed if implemented. ([github.com](https://github.com/OpenSourceRisk/Engine))

### 8.3 Open-source academic/industry-adjacent libraries implementing regression Monte Carlo and bounds for early exercise

**[GitHub: finmath/finmath-lib]** The finmath-lib Java library (about 551 GitHub stars at time of writing) includes Monte Carlo simulation for multi-asset Black–Scholes-type models and explicitly lists “American Monte-Carlo: Estimation of conditional expectations in a Monte-Carlo framework” plus Bermudan option lower/upper bound methods (regression continuation estimation and dual methods). This is industrially relevant as an implementation base for the LSM-type regression machinery that Hanbali–Linders compare against, and as an engineering reference for producing Greeks and bounds in early-exercise problems. No verified evidence was found that finmath-lib implements the comonotonic volatility-adjustment PDE approximation proposed by Hanbali–Linders. ([github.com](https://github.com/finmath/finmath-lib))

### 8.4 Open-source differentiable/PDE/MC toolkits including American option examples

**[GitHub: google/tf-quant-finance]** TF Quant Finance is an archived open-source TensorFlow-based quant library (about 5.2k GitHub stars at time of writing) that advertises tiers including PDE solvers, Ito process tooling, and provides an “American Option pricing under the Black–Scholes model” tutorial notebook. This is relevant as a practical implementation substrate for alternative fast solvers (automatic differentiation, GPU acceleration) that address similar scaling pressures as Hanbali–Linders’ method, but it does not constitute a verified implementation of the comonotonic-PDE approximation for basket options. ([github.com](https://github.com/google/tf-quant-finance))

## 9. Consolidated Reference List

[1] Karel J. in’t Hout; Jacob Snoeijer. “Numerical Valuation of American Basket Options via Partial Differential Complementarity Problems.” *Mathematics*, 2021, 9(13):1498. DOI: 10.3390/math9131498. `https://doi.org/10.3390/math9131498`

[2] Karel in’t Hout; Jacob Snoeijer. “Numerical valuation of American basket options via partial differential complementarity problems.” arXiv, 2021. arXiv:2106.01200. `https://arxiv.org/abs/2106.01200`

[3] Calypso Herrera; Florian Krach; Pierre Ruyssen; Josef Teichmann. “Optimal Stopping via Randomized Neural Networks.” arXiv, 2021 (revised Dec 2023). arXiv:2104.13669. `https://arxiv.org/abs/2104.13669`

[4] Malek Ben-Abdellatif; Hatem Ben-Ameur; Rim Chérif; Bruno Rémillard. “Dynamic Programming for Designing and Valuing Two-Dimensional Financial Derivatives.” *Risks*, 2024, 12(12):183. DOI: 10.3390/risks12120183. `https://doi.org/10.3390/risks12120183`

[5] Elettra Agliardi; Rossella Agliardi. “Pricing Multidimensional American Options.” *International Journal of Financial Studies*, 2023, 11(1):51. DOI: 10.3390/ijfs11010051. `https://doi.org/10.3390/ijfs11010051`

[6] GitHub repository: lballabio/QuantLib. “The QuantLib C++ library.” `https://github.com/lballabio/QuantLib`

[7] Klaus Spanderen; StatPro Italia srl (contributors in QuantLib source). “mclongstaffschwartzengine.cpp” (QuantLib test-suite source file including Longstaff–Schwartz engine headers). Debian Sources mirror of QuantLib source. `https://sources.debian.org/src/quantlib/1.9.1-1/test-suite/mclongstaffschwartzengine.cpp`

[8] QuantLib annotated source documentation. “mcamericanbasketengine.hpp Source File” (documents `MakeMCAmericanBasketEngine` and its Longstaff–Schwartz Monte Carlo basis). `https://rkapl123.github.io/QLAnnotatedSource/d5/d33/mcamericanbasketengine_8hpp_source.html`

[9] GitHub repository: OpenSourceRisk/Engine. “Open Source Risk Engine (ORE).” `https://github.com/OpenSourceRisk/Engine`

[10] GitHub repository: finmath/finmath-lib. “Mathematical Finance Library: Algorithms and methodologies related to mathematical finance.” `https://github.com/finmath/finmath-lib`

[11] GitHub repository: google/tf-quant-finance. “TF Quant Finance: TensorFlow based Quant Finance Library (ARCHIVED).” `https://github.com/google/tf-quant-finance`

---
Learn more:
1. [https://www.mdpi.com/2227-7390/9/13/1498](https://www.mdpi.com/2227-7390/9/13/1498)
2. [https://arxiv.org/abs/2104.13669](https://arxiv.org/abs/2104.13669)
3. [https://www.mdpi.com/2227-9091/12/12/183](https://www.mdpi.com/2227-9091/12/12/183)
4. [https://www.mdpi.com/2227-7072/11/1/51](https://www.mdpi.com/2227-7072/11/1/51)
5. [https://github.com/lballabio/QuantLib](https://github.com/lballabio/QuantLib)
6. [https://github.com/OpenSourceRisk/Engine](https://github.com/OpenSourceRisk/Engine)
7. [https://github.com/finmath/finmath-lib](https://github.com/finmath/finmath-lib)
8. [https://github.com/google/tf-quant-finance](https://github.com/google/tf-quant-finance)

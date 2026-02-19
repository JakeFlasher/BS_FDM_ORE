
# From Black-Scholes to Schrödinger: A Tutorial for Computer Scientists

---

## §0 — Global Roadmap

The entire argument fits into four transformations. Trace them left-to-right and top-to-bottom:

```
┌──────────────────────┐   x = ln S        ┌────────────────────────┐
│   Black-Scholes PDE  │   τ = T − t        │  BS in log-price /     │
│                      │ ─────────────────►  │  residual-time coords  │
│  V_t + ½σ²S²V_SS    │   (§3 Step 1)       │                        │
│  + rSV_S − rV = 0   │                     │  U_τ = ½σ²U_xx         │
└──────────────────────┘                     │   + (r−½σ²)U_x − rU   │
                                             └───────────┬────────────┘
                                                         │
                                          U = e^{−αx−βτ}Ψ  (§3 Step 2)
                                                         │
                                                         ▼
                          ┌──────────────────────────────────────────────┐
                          │          Heat (Diffusion) Equation           │
                          │         Ψ_τ  =  (σ²/2) Ψ_xx                 │
                          └──────────────────┬───────────────────────────┘
                                             │
                                      Wick rotation:  τ ──► it   (§4, §6)
                                             │
                                             ▼
                          ┌──────────────────────────────────────────────┐
                          │    Free-particle Schrödinger Equation        │
                          │    iℏ ψ_t  =  −(ℏ²/2m) ψ_xx                │
                          │    with   σ² ↔ ℏ/m                          │
                          └─────────────────────────────────────────────┘
```

The Schrödinger equation tells you how a quantum particle travels from now into the future; the Black-Scholes equation tells you how a future payoff collapses back to a present-day price. The two are the same equation — one in real time, the other in imaginary time.

---

## §1 — Prerequisites

### 1.1 Options and the No-Arbitrage Principle

A **European call option** on a stock with price \(S\) gives its holder the right — but not the obligation — to buy the stock at a pre-agreed **strike price** \(K\) on a future **expiry date** \(T\). At expiry the holder either exercises (if \(S_T > K\)) or walks away, so the payoff is \(\max(S_T - K, 0)\). A **European put** is the mirror image: the right to sell, with payoff \(\max(K - S_T, 0)\).

The **no-arbitrage principle** is the foundational axiom of mathematical finance. It states that no trading strategy should yield a riskless profit from zero initial capital. A direct corollary is the **law of one price**: two portfolios that guarantee the same payoff at time \(T\) must have the same price at every earlier time. This single principle, combined with a model for how \(S\) evolves, suffices to derive the Black-Scholes PDE.

### 1.2 Geometric Brownian Motion and Itô's Lemma

We model the stock price as a **Geometric Brownian Motion (GBM)**:

$$
dS = \mu\, S\, dt + \sigma\, S\, dW,
$$

where \(\mu\) is the expected return, \(\sigma > 0\) is the **volatility**, and \(W(t)\) is a standard Wiener process — a continuous-time random walk whose increments \(dW\) are independent Gaussian variables with mean zero and variance \(dt\). GBM guarantees \(S > 0\) and says that log-returns \(d(\ln S)\) are normally distributed, which matches empirical market data to a first approximation.

To differentiate a function of a stochastic variable we need:

**Lemma 1 (Itô's Lemma).** *Let \(f(S,t)\) be twice continuously differentiable and let \(S\) follow \(dS = \mu S\,dt + \sigma S\,dW\). Then*

$$
df = \Bigl(\frac{\partial f}{\partial t} + \mu S\frac{\partial f}{\partial S} + \tfrac{1}{2}\sigma^2 S^2 \frac{\partial^2 f}{\partial S^2}\Bigr)dt + \sigma S\frac{\partial f}{\partial S}\,dW.
$$

*Proof sketch.* Taylor-expand \(f(S+dS,\,t+dt)\) to second order in \(dS\):

$$
df = f_t\,dt + f_S\,dS + \tfrac{1}{2}f_{SS}(dS)^2 + \cdots
$$

The stochastic calculus rule that distinguishes Itô from ordinary calculus is \((dW)^2 = dt\) (in mean-square). Because \(dS = \mu S\,dt + \sigma S\,dW\), we compute \((dS)^2 = \sigma^2 S^2 (dW)^2 + \text{higher-order terms} = \sigma^2 S^2\,dt\). Substituting and dropping terms of order \((dt)^{3/2}\) and higher, we arrive at the formula above. \(\square\)

### 1.3 Classification of Second-Order Linear PDEs

A general second-order linear PDE in two variables \((x_1, x_2)\) has the form \(A\,u_{x_1 x_1} + 2B\,u_{x_1 x_2} + C\,u_{x_2 x_2} + \text{lower-order} = 0.\) Its character is determined by the discriminant \(\Delta = B^2 - AC\). The PDE is **elliptic** if \(\Delta < 0\) (exemplified by the Laplace equation \(\nabla^2 u = 0\)), **parabolic** if \(\Delta = 0\) (exemplified by the heat equation \(u_t = D\,u_{xx}\)), and **hyperbolic** if \(\Delta > 0\) (exemplified by the wave equation \(u_{tt} = c^2 u_{xx}\)). Parabolic equations describe diffusive, irreversible processes — exactly the character one expects for option prices that spread out over time.

### 1.4 The Free-Particle Schrödinger Equation

In quantum mechanics (QM), the state of a particle is encoded in a complex-valued **wave function** \(\psi(x,t)\). The probability of finding the particle near position \(x\) is \(|\psi(x,t)|^2\). The time-evolution of \(\psi\) is governed by the **Schrödinger equation**; for a free particle (no external potential) in one dimension it reads

$$
i\hbar\,\frac{\partial \psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2 \psi}{\partial x^2},
$$

where \(\hbar\) is the reduced Planck constant and \(m\) is the particle's mass. Dividing both sides by \(i\hbar\):

$$
\frac{\partial \psi}{\partial t} = \frac{i\hbar}{2m}\frac{\partial^2 \psi}{\partial x^2}.
$$

Notice the structure: it is a first-order equation in time, second-order in space — parabolic. The crucial novelty compared with the heat equation is the imaginary unit \(i\) multiplying the diffusion coefficient. That single factor is the entire gap between quantum mechanics and finance, as we shall see.

---

## §2 — Derivation of the Black-Scholes PDE

**Why this step?** We want a deterministic equation for \(V(S,t)\), the option price, even though \(S\) is random. The trick is **delta hedging**: construct a portfolio in which the randomness of the stock exactly cancels the randomness of the option, producing a riskless instrument. No-arbitrage then forces that instrument to earn the risk-free rate \(r\).

**Proposition 1 (Black-Scholes PDE).** *Under GBM with constant \(\sigma\) and risk-free rate \(r\), the price \(V(S,t)\) of any European derivative satisfies*

$$
\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2\frac{\partial^2 V}{\partial S^2} + r\,S\frac{\partial V}{\partial S} - r\,V = 0.
$$

*Proof.* Construct the portfolio \(\Pi = V - \Delta\,S\), where \(\Delta\) is a number of shares to be chosen momentarily. Its instantaneous change is

$$
d\Pi = dV - \Delta\,dS.
$$

Apply Itô's lemma (Lemma 1) to \(V(S,t)\):

$$
dV = \Bigl(V_t + \mu S V_S + \tfrac{1}{2}\sigma^2 S^2 V_{SS}\Bigr)dt + \sigma S V_S\,dW.
$$

Substitute and group terms:

$$
d\Pi = \Bigl(V_t + \mu S V_S + \tfrac{1}{2}\sigma^2 S^2 V_{SS} - \Delta\,\mu S\Bigr)dt + \sigma S(V_S - \Delta)\,dW.
$$

Choose \(\Delta = V_S\) (the "delta" of the option). The entire stochastic \(dW\) term vanishes, leaving a purely deterministic increment:

$$
d\Pi = \Bigl(V_t + \tfrac{1}{2}\sigma^2 S^2 V_{SS}\Bigr)dt.
$$

A riskless portfolio must earn the risk-free rate, so \(d\Pi = r\,\Pi\,dt = r(V - V_S\,S)\,dt\). Equating:

$$
V_t + \tfrac{1}{2}\sigma^2 S^2 V_{SS} = r\,V - r\,S\,V_S.
$$

Rearranging yields the Black-Scholes PDE. Note that \(\mu\) has dropped out completely: the expected return on the stock is irrelevant to option pricing, which is the mathematical content of **risk-neutral valuation**. \(\square\)

**Sanity check.** Set \(\sigma = 0\). The PDE becomes \(V_t + rSV_S - rV = 0\), which is the first-order transport equation whose solution is \(V(S,t) = e^{-r(T-t)}\,\text{payoff}(Se^{r(T-t)})\) — precisely the discounted payoff along the deterministic growth path \(S \to Se^{r(T-t)}\). This is correct: with zero volatility, the stock price is deterministic and the option can be priced by simple discounting.

---

## §3 — Core Derivation: Black-Scholes → Heat Equation

### Step 1: Log-Price Substitution

**Why this step?** The BS PDE has variable coefficients \(S^2\) and \(S\) because price changes are *multiplicative*: a stock moving from 100 to 50 is the same proportional loss as 10 to 5. In log-coordinates the equation becomes *translation-invariant*. We also reverse the time arrow because option value is determined by a future boundary condition (the payoff at expiry \(T\)).

Define

$$
x = \ln S, \qquad \tau = T - t, \qquad U(x,\tau) = V(S,t).
$$

**Lemma 2 (Change of variables).** *Under the substitution above, the Black-Scholes PDE becomes*

$$
\frac{\partial U}{\partial \tau} = \frac{\sigma^2}{2}\frac{\partial^2 U}{\partial x^2} + \Bigl(r - \frac{\sigma^2}{2}\Bigr)\frac{\partial U}{\partial x} - r\,U.
$$

*Proof.* Since \(S = e^x\), we have \(\partial x/\partial S = 1/S\). Apply the chain rule to each derivative in the BS PDE.

**Time derivative.** Because \(\tau = T - t\) with \(T\) constant, \(d\tau = -dt\), so

$$
\frac{\partial V}{\partial t} = \frac{\partial U}{\partial \tau}\cdot\frac{\partial \tau}{\partial t} = -\frac{\partial U}{\partial \tau}.
$$

**First spatial derivative.** With \(x = \ln S\):

$$
\frac{\partial V}{\partial S} = \frac{\partial U}{\partial x}\cdot\frac{\partial x}{\partial S} = \frac{1}{S}\frac{\partial U}{\partial x}.
$$

**Second spatial derivative.** Differentiate the above once more in \(S\):

$$
\frac{\partial^2 V}{\partial S^2} = \frac{\partial}{\partial S}\!\Bigl(\frac{1}{S}\frac{\partial U}{\partial x}\Bigr) = -\frac{1}{S^2}\frac{\partial U}{\partial x} + \frac{1}{S}\cdot\frac{1}{S}\frac{\partial^2 U}{\partial x^2} = \frac{1}{S^2}\Bigl(\frac{\partial^2 U}{\partial x^2} - \frac{\partial U}{\partial x}\Bigr).
$$

Now substitute into the BS PDE \(V_t + \tfrac{1}{2}\sigma^2 S^2 V_{SS} + rSV_S - rV = 0\):

$$
-\frac{\partial U}{\partial \tau} + \frac{\sigma^2}{2}\,S^2\cdot\frac{1}{S^2}\Bigl(U_{xx} - U_x\Bigr) + r\,S\cdot\frac{1}{S}\,U_x - r\,U = 0.
$$

The factors of \(S\) and \(S^2\) cancel identically, leaving

$$
-U_\tau + \frac{\sigma^2}{2}\bigl(U_{xx} - U_x\bigr) + r\,U_x - r\,U = 0.
$$

Collect the \(U_x\) terms: \(-\sigma^2/2 + r\). Rearrange to isolate \(U_\tau\):

$$
\frac{\partial U}{\partial \tau} = \frac{\sigma^2}{2}\,U_{xx} + \Bigl(r - \frac{\sigma^2}{2}\Bigr)U_x - r\,U. \qquad\square
$$

This is a constant-coefficient parabolic PDE — a drift–diffusion equation with a decay term. The log-price substitution has removed all variable coefficients.

**Sanity check (dimensional).** Every term has dimensions of \([\text{currency}]/[\text{time}]\). The coefficient \(\sigma^2/2\) has units \(\text{time}^{-1}\), and \(U_{xx}\) has units of \([\text{currency}]\) because \(x\) is dimensionless (a logarithm). Likewise \(r\) has units \(\text{time}^{-1}\). All terms are consistent. ✓

### Step 2: Gauge (Eigenfunction) Transform

**Why this step?** The equation still contains a first-order derivative (drift) and a zeroth-order term (decay). Both are physically meaningful — the drift encodes the risk-neutral expected return, and the decay encodes discounting — but mathematically they clutter the PDE. We perform an exponential change of dependent variable, much like "completing the square," to absorb these terms into a prefactor and leave behind the pure heat equation.

Write

$$
U(x,\tau) = e^{-\alpha x - \beta\tau}\,\Psi(x,\tau),
$$

where constants \(\alpha,\beta\) are to be determined so that the equation for \(\Psi\) contains only \(\Psi_\tau\) and \(\Psi_{xx}\).

**Proposition 2 (Gauge reduction to the heat equation).** *With*

$$
\alpha = \frac{r}{\sigma^2} - \frac{1}{2}, \qquad \beta = r + \frac{(r - \sigma^2/2)^2}{2\sigma^2} = \frac{(r + \sigma^2/2)^2}{2\sigma^2},
$$

*the function \(\Psi\) satisfies the heat equation*

$$
\frac{\partial \Psi}{\partial \tau} = \frac{\sigma^2}{2}\,\frac{\partial^2 \Psi}{\partial x^2}.
$$

*Proof.* Compute the needed derivatives of \(U = e^{-\alpha x - \beta\tau}\Psi\) (we suppress arguments for clarity).

**Time:**

$$
U_\tau = e^{-\alpha x - \beta\tau}\bigl(\Psi_\tau - \beta\,\Psi\bigr).
$$

**First space derivative:**

$$
U_x = e^{-\alpha x - \beta\tau}\bigl(\Psi_x - \alpha\,\Psi\bigr).
$$

**Second space derivative:** Differentiate \(U_x\) once more:

$$
U_{xx} = e^{-\alpha x - \beta\tau}\bigl(\Psi_{xx} - 2\alpha\,\Psi_x + \alpha^2\Psi\bigr).
$$

Define \(\nu := r - \sigma^2/2\) for brevity. Substitute all three expressions into the log-price BS equation \(U_\tau = \tfrac{\sigma^2}{2}U_{xx} + \nu\,U_x - r\,U\) and divide both sides by the common factor \(e^{-\alpha x - \beta\tau}\):

$$
\Psi_\tau - \beta\Psi = \frac{\sigma^2}{2}\bigl(\Psi_{xx} - 2\alpha\Psi_x + \alpha^2\Psi\bigr) + \nu\bigl(\Psi_x - \alpha\Psi\bigr) - r\Psi.
$$

Expand and sort by derivative order:

$$
\Psi_\tau = \frac{\sigma^2}{2}\Psi_{xx} + \bigl(\underbrace{-\sigma^2\alpha + \nu}_{\text{coeff. of }\Psi_x}\bigr)\Psi_x + \bigl(\underbrace{\tfrac{\sigma^2}{2}\alpha^2 - \nu\alpha - r + \beta}_{\text{coeff. of }\Psi}\bigr)\Psi.
$$

**Eliminating the first-order term.** Set the coefficient of \(\Psi_x\) to zero:

$$
-\sigma^2\alpha + \nu = 0 \;\;\Longrightarrow\;\; \alpha = \frac{\nu}{\sigma^2} = \frac{r - \sigma^2/2}{\sigma^2} = \frac{r}{\sigma^2} - \frac{1}{2}.
$$

**Eliminating the zeroth-order term.** Set the coefficient of \(\Psi\) to zero. Note that, with the value of \(\alpha\) just found, \(\nu = \sigma^2 \alpha\), so \(\nu\alpha = \sigma^2\alpha^2\). The zeroth-order coefficient becomes

$$
\frac{\sigma^2\alpha^2}{2} - \sigma^2\alpha^2 - r + \beta = -\frac{\sigma^2\alpha^2}{2} - r + \beta.
$$

Setting this to zero:

$$
\beta = r + \frac{\sigma^2\alpha^2}{2} = r + \frac{(r - \sigma^2/2)^2}{2\sigma^2}.
$$

To see the compact form, expand \(\alpha + 1 = r/\sigma^2 + 1/2 = (r + \sigma^2/2)/\sigma^2\), whence \(\sigma^2(\alpha+1)^2/2 = (r+\sigma^2/2)^2/(2\sigma^2)\). Adding and simplifying shows this equals \(r + \sigma^2\alpha^2/2\), confirming the identity

$$
\beta = \frac{(r + \sigma^2/2)^2}{2\sigma^2}.
$$

With both unwanted terms eliminated, we are left with

$$
\frac{\partial\Psi}{\partial\tau} = \frac{\sigma^2}{2}\,\frac{\partial^2\Psi}{\partial x^2}. \qquad\square
$$

**Sanity check (special case).** Let \(r = 0\). Then \(\alpha = -1/2\), \(\beta = \sigma^2/8\), and the gauge factor becomes \(e^{x/2 - \sigma^2\tau/8}\). The original PDE with \(r=0\) is a pure diffusion with drift \(-\sigma^2/2\), and the exponential tilt \(e^{x/2}\) is exactly the Girsanov-type change of measure that removes that drift. ✓

---

## §4 — Structural Comparison: Heat Equation vs. Schrödinger Equation

Place the two equations side by side:

```
  ┌──────────────────────────────────────┬──────────────────────────────────────┐
  │          HEAT EQUATION               │      SCHRÖDINGER EQUATION            │
  │        (Black-Scholes core)          │       (free particle, 1-D)           │
  ├──────────────────────────────────────┼──────────────────────────────────────┤
  │                                      │                                      │
  │   ∂Ψ         σ²  ∂²Ψ                │        ∂ψ        iℏ  ∂²ψ            │
  │  ──── =  ─── ─────                  │   ──── =  ─── ─────                 │
  │   ∂τ         2  ∂x²                 │    ∂t         2m ∂x²                │
  │                                      │                                      │
  ├──────────────────────────────────────┼──────────────────────────────────────┤
  │  "time" variable:  τ (real)          │  time variable:  t (real)            │
  │  diffusion coeff:  σ²/2  (real >0)  │  diffusion coeff:  iℏ/2m  (imag.)   │
  │  solution Ψ:  real-valued            │  solution ψ:  complex-valued         │
  │  evolution:  dissipative             │  evolution:  unitary                  │
  └──────────────────────────────────────┴──────────────────────────────────────┘
```

Performing a **Wick rotation** \(t \to -i\tau\) on the Schrödinger equation converts \(i\hbar/(2m)\) into the real quantity \(\hbar/(2m)\). Matching diffusion coefficients then yields the **parameter correspondence**

$$
\frac{\sigma^2}{2} = \frac{\hbar}{2m} \qquad\Longleftrightarrow\qquad \sigma^2 = \frac{\hbar}{m}.
$$

The formal derivation is carried out in Proposition 3 below.

**Proposition 3 (Wick rotation).** *Under the substitution \(t = -i\tau\), the free-particle Schrödinger equation \(i\hbar\,\psi_t = -(\hbar^2/2m)\psi_{xx}\) becomes the heat equation \(\psi_\tau = (\hbar/2m)\psi_{xx}\).*

*Proof.* From \(t = -i\tau\) we get \(\partial/\partial t = (\partial\tau/\partial t)\,\partial/\partial\tau = i\,\partial/\partial\tau\). Substituting:

$$
i\hbar\cdot i\,\psi_\tau = -\frac{\hbar^2}{2m}\psi_{xx} \;\;\Longrightarrow\;\; -\hbar\,\psi_\tau = -\frac{\hbar^2}{2m}\psi_{xx} \;\;\Longrightarrow\;\; \psi_\tau = \frac{\hbar}{2m}\psi_{xx}. \qquad\square
$$

---

## §5 — Physical Intuition: \(\sigma^2 \leftrightarrow \hbar/m\)

The correspondence \(\sigma^2 = \hbar/m\) admits a vivid physical reading. The constant \(\hbar\) quantifies the intrinsic randomness of the quantum world — it sets the scale below which deterministic trajectories cease to exist. Volatility \(\sigma\) plays the same role in markets: it measures the irreducible randomness of log-returns. It is therefore natural to set \(\hbar = 1\) (as physicists routinely do), giving

$$
m = \frac{1}{\sigma^2}.
$$

In quantum mechanics, the mass \(m\) of a particle measures its resistance to changes in velocity. A heavy particle's wave-packet spreads slowly; a light particle's wave-packet spreads quickly. Now translate this into finance: a **low-volatility** asset (a large-cap index, say) behaves like a **heavy** particle — its probability distribution broadens slowly over time, concentrated near its expected value. A **high-volatility** asset (a meme coin or penny stock) is a **light** particle — the faintest breath of news sends its probability density sprawling across a vast range of prices. The ASCII sketch below illustrates one unit of time's worth of probability-density evolution for two hypothetical assets:

```
  LOW-σ  (heavy particle,  σ = 0.1,  m = 100)        HIGH-σ  (light particle,  σ = 0.8,  m ≈ 1.6)

  P(x)                                                P(x)
   │      ▓▓                                           │
   │     ▓▓▓▓                                          │       ░░░░░░░░░░░░░░░░░░░░░
   │    ▓▓▓▓▓▓                                         │     ░░░░░░░░░░░░░░░░░░░░░░░░░
   │   ▓▓▓▓▓▓▓▓                                        │   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
   │  ▓▓▓▓▓▓▓▓▓▓                                       │  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
   └──────────────────► x                               └──────────────────────────────────► x
      (narrow Gaussian)                                     (wide, flat Gaussian)
```

Both distributions are Gaussians with variance \(\sigma^2 \tau\). For the heavy particle the variance is small (tight peak); for the light particle it is large (spread out). This directly mirrors the quantum-mechanical uncertainty principle: for a given "momentum," a lighter particle has greater positional uncertainty.

---

## §6 — The Role of the Imaginary Unit \(i\): Probability vs. Probability Amplitude

The single factor of \(i\) separating the heat equation from the Schrödinger equation has profound consequences for the character of solutions.

**Why quantum mechanics needs complex numbers.** In QM, the wave function \(\psi\) is a **probability amplitude**, not a probability. Probabilities are obtained by taking \(|\psi|^2\). Because \(\psi\) is complex, two contributions to \(\psi\) at the same point can cancel — this is **interference**, responsible for phenomena from electron diffraction to quantum computing. The Schrödinger time-evolution operator \(e^{-iHt/\hbar}\) is **unitary**: it preserves the norm \(\int|\psi|^2 dx = 1\). Information is never lost; probability is reshuffled among positions but never created or destroyed.

**Why finance does not need complex numbers.** In the risk-neutral pricing framework the option price is *directly* a probability-weighted expectation — no squaring is required. The heat-equation propagator \(e^{-H\tau}\) (where \(H = -\sigma^2\partial_{xx}/2\)) is **not** unitary; it is a contraction semigroup. It smooths and damps: sharp features in the initial payoff get blurred as \(\tau\) grows. This reflects the economic reality that the further into the future the payoff lies, the more uncertainty there is about its value, and the more it is discounted. In physical language, the heat equation is **dissipative** while the Schrödinger equation is **conservative**.

The table below summarises the distinction:

```
  ┌─────────────────────────┬────────────────────────────────┐
  │  SCHRÖDINGER  (QM)      │  HEAT / BLACK-SCHOLES          │
  ├─────────────────────────┼────────────────────────────────┤
  │  ψ  ∈  ℂ               │  Ψ  ∈  ℝ                      │
  │  |ψ|² = probability    │  Ψ  itself ∝ probability       │
  │  evolution is unitary   │  evolution is dissipative       │
  │  interference possible  │  no interference                │
  │  norm ∫|ψ|²=1 conserved│  total "mass" can decay         │
  └─────────────────────────┴────────────────────────────────┘
```

---

## §7 — Path Integrals, Wick Rotation, and Option Pricing

### 7.1 Feynman's Path Integral

Richard Feynman reformulated quantum mechanics by replacing the Schrödinger equation with a summation over histories. The **quantum-mechanical propagator** — the Green's function that evolves a wave function from \((x_i, t_i)\) to \((x_f, t_f)\) — is written as

$$
K(x_f,t_f;\,x_i,t_i) = \int \mathcal{D}[x(t)]\;\exp\!\Bigl(\frac{i}{\hbar}S[x(t)]\Bigr),
$$

where \(S[x(t)] = \int_{t_i}^{t_f}\frac{1}{2}m\dot{x}^2\,dt\) is the classical **action** for a free particle and the integral is taken over *all* continuous paths \(x(t)\) connecting \(x_i\) to \(x_f\). (Defining the path-integral measure \(\mathcal{D}[x]\) rigorously requires careful limiting procedures; in practice one discretises time into \(N\) steps and takes \(N\to\infty\). We note the issue but do not belabor it here.)

Paths near the classical trajectory (which extremises \(S\)) contribute coherently; paths far from it oscillate wildly and cancel — this is the **stationary-phase** mechanism behind classical mechanics emerging from quantum mechanics.

### 7.2 Wick Rotation: From Oscillation to Damping

Apply the Wick rotation \(t = -i\tau\). As derived in the proof of Proposition 3, the time derivative transforms as \(\partial_t = i\,\partial_\tau\). For the action one finds, by tracking the substitution \(dt = -i\,d\tau\) and \(\dot{x} = dx/dt = i\,dx/d\tau\):

$$
\frac{i}{\hbar}S = \frac{i}{\hbar}\int\frac{m}{2}\dot{x}^2\,dt = \frac{i}{\hbar}\int\frac{m}{2}(i\dot{x}_E)^2(-i\,d\tau) = -\frac{1}{\hbar}\int\frac{m}{2}\dot{x}_E^2\,d\tau = -\frac{S_E}{\hbar},
$$

where \(\dot{x}_E = dx/d\tau\) and \(S_E = \int\frac{m}{2}\dot{x}_E^2\,d\tau\) is the **Euclidean action**. The oscillatory weight \(e^{iS/\hbar}\) has become a real, decaying weight \(e^{-S_E/\hbar}\). The path integral is now a mathematically well-defined Wiener integral, and extreme paths are suppressed exponentially rather than cancelled by interference.

Inserting the financial identifications \(\hbar = 1\), \(m = 1/\sigma^2\):

$$
e^{-S_E/\hbar} = \exp\!\Bigl(-\int_0^\tau \frac{1}{2\sigma^2}\dot{x}_E^2\,d\tau'\Bigr).
$$

This is precisely the **Wiener measure** on Brownian-motion paths in log-price space, weighted by volatility \(\sigma\).

### 7.3 The Propagator as Risk-Neutral Transition Density

The Green's function of the heat equation \(\Psi_\tau = \frac{\sigma^2}{2}\Psi_{xx}\), subject to the initial condition \(\Psi(x,0) = \delta(x - x_0)\), is the Gaussian

$$
G(x,\tau;\,x_0) = \frac{1}{\sqrt{2\pi\sigma^2\tau}}\;\exp\!\Bigl(-\frac{(x - x_0)^2}{2\sigma^2\tau}\Bigr).
$$

**Lemma 3.** *The function \(G\) defined above satisfies \(G_\tau = \frac{\sigma^2}{2}G_{xx}\).*

*Proof sketch.* Direct computation gives \(G_\tau = G\bigl[(x-x_0)^2/(2\sigma^2\tau^2) - 1/(2\tau)\bigr]\) and \(G_{xx} = G\bigl[(x-x_0)^2/(\sigma^4\tau^2) - 1/(\sigma^2\tau)\bigr]\). Multiplying \(G_{xx}\) by \(\sigma^2/2\) reproduces \(G_\tau\). \(\square\)

This Gaussian is the **transition density of Brownian motion** in log-price space: it gives the probability that the log-price moves from \(x_0\) to \(x\) in time \(\tau\). The general solution is

$$
\Psi(x,\tau) = \int_{-\infty}^{\infty}G(x,\tau;\,x_0)\;\Psi(x_0,0)\,dx_0,
$$

where \(\Psi(x_0,0)\) is the gauge-transformed payoff at expiry. Unwinding the gauge and log-price substitutions yields the celebrated **risk-neutral pricing formula** familiar from finance textbooks:

$$
V(S,t) = e^{-r(T-t)}\,\mathbb{E}^{Q}\!\bigl[\text{payoff}(S_T)\;\big|\;S_t = S\bigr],
$$

where the expectation is taken under the risk-neutral measure \(Q\), under which the log-price drifts at rate \(r - \sigma^2/2\) (not \(\mu - \sigma^2/2\)). The factor \(e^{-r(T-t)}\) is part of the gauge prefactor \(e^{-\alpha x - \beta\tau}\), and the Gaussian convolution is the risk-neutral expected payoff. The Feynman path-integral viewpoint thus gives the same answer: the option's present value is a sum over all possible log-price paths, each weighted by the exponential of minus its Euclidean action, divided by \(\sigma^2\). Extreme price paths (enormous rallies, devastating crashes) carry exponentially small weight; typical paths dominate, just as near-classical trajectories dominate the quantum propagator.

---

## §8 — Summary: Concept-by-Concept Correspondence

```
┌──────────────────────────┬────────────────────────────┬───────────────────────────────┐
│        CONCEPT           │   QUANTUM MECHANICS        │   BLACK-SCHOLES / FINANCE     │
├──────────────────────────┼────────────────────────────┼───────────────────────────────┤
│ Fundamental equation     │ iℏ ψ_t = −(ℏ²/2m)ψ_xx    │ Ψ_τ = (σ²/2)Ψ_xx             │
│ Time variable            │ t  (real, forward)         │ τ = T−t  (real, backward)     │
│ Relation of times        │ t = −iτ  (Wick rotation)  │ τ = it   (Wick rotation)      │
│ "Position" variable      │ x  (spatial coordinate)    │ x = ln S  (log-price)         │
│ State function           │ ψ(x,t)  ∈ ℂ              │ Ψ(x,τ)  ∈ ℝ                  │
│ Probability              │ |ψ|²                      │ Ψ directly (after gauge)      │
│ Diffusion coefficient    │ iℏ / 2m                    │ σ² / 2                        │
│ Intrinsic randomness     │ ℏ  (Planck's constant)    │ σ  (volatility)               │
│ Inertia                  │ m  (particle mass)         │ 1/σ²  ("market mass")         │
│ Heavy particle ↔         │ Slow spreading             │ Blue-chip / large-cap index   │
│ Light particle ↔         │ Fast spreading             │ Meme coin / penny stock       │
│ External potential V(x)  │ Shifts energy levels       │ Risk-free rate r (discount)   │
│ Drift / vector potential │ Constant force field        │ Risk-neutral drift r − σ²/2  │
│ Path-integral weight     │ e^{iS/ℏ}  (oscillatory)   │ e^{−S_E/σ²}  (decaying)      │
│ Propagator / Green's fn  │ K(x_f,t; x_i,0)           │ G(x,τ; x₀,0)  (Gaussian)     │
│ Principle                │ Least action (stationary φ)│ Risk-neutral pricing          │
│ Superposition character  │ Interference (complex ψ)   │ No interference (real Ψ)      │
│ Norm conservation        │ Unitarity: ∫|ψ|²=1        │ Dissipation: value decays     │
│ Interpretation           │ "Where does particle go?"  │ "What is the option worth?"   │
└──────────────────────────┴────────────────────────────┴───────────────────────────────┘
```

To summarise in one sentence: the Black-Scholes model *is* a quantum-mechanical free particle living in imaginary time, where Planck's constant becomes the volatility and the particle's mass becomes the inverse of variance. The Schrödinger equation asks "given a particle here and now, where might it be later?"; the Black-Scholes equation asks "given a payoff there and later, what is it worth here and now?" The mathematics is identical; only the direction of the question — and one factor of \(i\) — differs.
## Appendices and Worked Examples

What follows extends the main tutorial with material that solidifies the concepts: explicit payoff diagrams promised in §1, a fully worked numerical example tracing the entire transformation chain, a deeper proof of the Green's function, and pointers for CS students interested in computational applications.

---

## Appendix A — Option Payoff Diagrams

The prerequisite discussion in §1 mentioned payoffs verbally. The ASCII diagrams below make them concrete.

**European Call** with strike \(K = 100\):

```
  Payoff
    │
 60 ┤                                          ╱
    │                                        ╱
 40 ┤                                      ╱
    │                                    ╱
 20 ┤                                  ╱
    │                                ╱
  0 ┤━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╱─────────────
    └───────┬───────┬──────┬──────┬──────┬──────► S
           40      60     80    100    120    160
                                  K
         payoff = max(S − K, 0)
```

**European Put** with strike \(K = 100\):

```
  Payoff
    │
 80 ┤╲
    │  ╲
 60 ┤    ╲
    │      ╲
 40 ┤        ╲
    │          ╲
 20 ┤            ╲
    │              ╲
  0 ┤────────────────╲━━━━━━━━━━━━━━━━━━━━━━━━━━
    └───────┬───────┬──────┬──────┬──────┬──────► S
           20      40     60    100    120    160
                                 K
         payoff = max(K − S, 0)
```

The kink at \(S = K\) is the source of all the mathematical richness in option pricing. At expiry the payoff is a piecewise-linear, non-smooth function. The role of the Black-Scholes PDE (equivalently, the heat equation after our transformations) is to "smooth out" this kink as we move backward in time from expiry to the present — exactly the way the heat equation smooths a discontinuous initial temperature profile. This is not merely an analogy; it is the same mathematical mechanism.

---

## Appendix B — Full Worked Example: European Call with Concrete Numbers

We price a European call option with parameters \(K = 100\), \(S_0 = 100\), \(T = 1\) year, \(r = 0.05\), \(\sigma = 0.20\). We trace every transformation explicitly.

### B.1 The Initial (Terminal) Condition

At expiry (\(t = T\), i.e.\ \(\tau = 0\)), the call payoff is

$$
V(S, T) = \max(S - K, 0).
$$

### B.2 Step 1 — Log-Price Substitution

Set \(x = \ln S\) and \(\tau = T - t\). At \(\tau = 0\):

$$
U(x, 0) = V(e^x, T) = \max(e^x - 100,\; 0).
$$

The kink is now at \(x_0 = \ln 100 \approx 4.6052\). For \(x < x_0\), the payoff is zero; for \(x > x_0\), it is \(e^x - 100\).

### B.3 Step 2 — Gauge Transform

Compute the constants with our parameter values:

$$
\alpha = \frac{r}{\sigma^2} - \frac{1}{2} = \frac{0.05}{0.04} - 0.5 = 1.25 - 0.5 = 0.75.
$$

$$
\beta = \frac{(r + \sigma^2/2)^2}{2\sigma^2} = \frac{(0.05 + 0.02)^2}{2 \times 0.04} = \frac{(0.07)^2}{0.08} = \frac{0.0049}{0.08} = 0.06125.
$$

The gauge-transformed initial condition at \(\tau = 0\) is therefore

$$
\Psi(x, 0) = e^{\alpha x + \beta \cdot 0}\,U(x, 0) = e^{0.75\,x}\,\max(e^x - 100,\; 0).
$$

For \(x > \ln 100\), this is \(e^{0.75x}(e^x - 100) = e^{1.75x} - 100\,e^{0.75x}\).

### B.4 Heat-Equation Propagation

With diffusion coefficient \(D = \sigma^2/2 = 0.02\), the solution at time \(\tau\) is

$$
\Psi(x, \tau) = \int_{-\infty}^{\infty} \frac{1}{\sqrt{2\pi \cdot 0.04\,\tau}}\; \exp\!\Bigl(-\frac{(x - x')^2}{0.08\,\tau}\Bigr)\;\Psi(x', 0)\;dx'.
$$

At \(\tau = 1\) (i.e.\ the present, since \(T = 1\)) with \(x = \ln 100\):

$$
\Psi(\ln 100,\; 1) = \int_{\ln 100}^{\infty} \frac{1}{\sqrt{0.04\pi \cdot 2}}\; \exp\!\Bigl(-\frac{(\ln 100 - x')^2}{0.08}\Bigr)\;\bigl(e^{1.75x'} - 100\,e^{0.75x'}\bigr)\;dx'.
$$

These are Gaussian integrals of the form \(\int e^{ax' - b(x')^2}dx'\), which evaluate to expressions involving the standard normal CDF \(\Phi\). After carrying through the algebra (which is precisely how the famous Black-Scholes formula is derived), one obtains

$$
V(S_0, 0) = S_0\,\Phi(d_1) - K\,e^{-rT}\,\Phi(d_2),
$$

with

$$
d_1 = \frac{\ln(S_0/K) + (r + \sigma^2/2)\,T}{\sigma\sqrt{T}}, \qquad d_2 = d_1 - \sigma\sqrt{T}.
$$

### B.5 Numerical Evaluation

Plug in numbers:

$$
d_1 = \frac{\ln(100/100) + (0.05 + 0.02)\cdot 1}{0.20 \cdot 1} = \frac{0 + 0.07}{0.20} = 0.35.
$$

$$
d_2 = 0.35 - 0.20 = 0.15.
$$

From standard normal tables: \(\Phi(0.35) \approx 0.6368\) and \(\Phi(0.15) \approx 0.5596\). Therefore

$$
V = 100 \times 0.6368 - 100 \times e^{-0.05} \times 0.5596 = 63.68 - 95.12 \times 0.5596 \approx 63.68 - 53.23 = 10.45.
$$

The at-the-money call with these parameters is worth approximately \(\$10.45\).

**Sanity check.** A rough lower bound is the "intrinsic value" discounted: \(\max(S_0 - Ke^{-rT}, 0) = 100 - 95.12 = 4.88\). The option must be worth more than this (because it also has "time value" from volatility), and indeed \(10.45 > 4.88\). An upper bound is \(S_0 = 100\) (no option can be worth more than the stock), and \(10.45 < 100\). Both bounds hold. ✓

---

## Appendix C — Proof that the Gaussian is the Fundamental Solution

In §7.3 we stated Lemma 3 without full detail. Here we supply the complete verification.

**Lemma 3 (restated).** *The function*

$$
G(x, \tau;\, x_0) = \frac{1}{\sqrt{2\pi\sigma^2\tau}}\;\exp\!\Bigl(-\frac{(x - x_0)^2}{2\sigma^2\tau}\Bigr)
$$

*satisfies \(G_\tau = \frac{\sigma^2}{2}\,G_{xx}\) for all \(\tau > 0\), and \(G(\cdot, \tau;\, x_0) \to \delta(\cdot - x_0)\) as \(\tau \to 0^+\).*

*Proof.* Denote \(\xi = x - x_0\) and \(D = \sigma^2/2\). Then \(G = (4\pi D\tau)^{-1/2}\exp(-\xi^2/(4D\tau))\).

**Compute \(G_\tau\).** Write \(G = (4\pi D\tau)^{-1/2} \exp\bigl(-\xi^2/(4D\tau)\bigr)\). The logarithmic derivative in \(\tau\) is

$$
\frac{\partial}{\partial\tau}\ln G = -\frac{1}{2\tau} + \frac{\xi^2}{4D\tau^2}.
$$

Therefore

$$
G_\tau = G\Bigl(-\frac{1}{2\tau} + \frac{\xi^2}{4D\tau^2}\Bigr). \tag{C.1}
$$

**Compute \(G_x\) and \(G_{xx}\).** Since \(\xi = x - x_0\), we have \(\partial\xi/\partial x = 1\):

$$
G_x = G\cdot\Bigl(-\frac{\xi}{2D\tau}\Bigr).
$$

Differentiate once more using the product rule:

$$
G_{xx} = G_x\cdot\Bigl(-\frac{\xi}{2D\tau}\Bigr) + G\cdot\Bigl(-\frac{1}{2D\tau}\Bigr) = G\Bigl(\frac{\xi^2}{4D^2\tau^2} - \frac{1}{2D\tau}\Bigr). \tag{C.2}
$$

**Verify the PDE.** Multiply (C.2) by \(D\):

$$
D\,G_{xx} = G\Bigl(\frac{\xi^2}{4D\tau^2} - \frac{1}{2\tau}\Bigr).
$$

This is identical to (C.1). Hence \(G_\tau = D\,G_{xx} = \frac{\sigma^2}{2}\,G_{xx}\). \(\square\)

The delta-function initial condition follows from the standard fact that a normalised Gaussian with variance \(\sigma^2\tau \to 0\) concentrates all its unit mass at the origin. (Formally, for any smooth test function \(\phi\), one shows \(\int G\,\phi\,dx \to \phi(x_0)\) by the substitution \(u = \xi/\sqrt{2\sigma^2\tau}\) and dominated convergence.)

---

## Appendix D — Physical Meaning of the Drift and Decay Terms

In §3, Step 2, we stripped away the drift and decay to reach the pure heat equation. Here we interpret those removed terms more carefully, since they carry financial meaning that the gauge transform conceals.

Return to the log-price BS equation:

$$
U_\tau = \frac{\sigma^2}{2}\,U_{xx} + \underbrace{\Bigl(r - \frac{\sigma^2}{2}\Bigr)}_{\text{drift } \nu}\,U_x \;-\; \underbrace{r\,U}_{\text{decay}}.
$$

**The drift term** \(\nu\,U_x\). In a reference frame co-moving with the drift, the center of the log-price distribution slides at speed \(\nu = r - \sigma^2/2\). This is the well-known **risk-neutral drift**: under the risk-neutral measure \(Q\), the expected log-return per unit time is not the real-world \(\mu - \sigma^2/2\) but rather \(r - \sigma^2/2\). The subtraction of \(\sigma^2/2\) is the **Itô correction** — the difference between the expectation of \(\ln S\) and the log of the expectation of \(S\). In the quantum analogy, a constant drift corresponds to a particle in a uniform force field (constant vector potential), causing an overall translation of the wave packet without changing its shape.

**The decay term** \(-rU\). This multiplicative damping shrinks the amplitude of \(U\) at rate \(r\). Financially, this is the **time value of money**: a dollar received in the future is worth less than a dollar today, discounted at rate \(r\). In the Schrödinger analogy, this corresponds to a constant scalar potential \(V_0 = r\), which shifts all energy eigenvalues by \(r\) without altering the spatial profile of eigenstates. The gauge transform \(U = e^{-\alpha x - \beta\tau}\Psi\) absorbs both the drift (via the \(e^{-\alpha x}\) spatial tilt) and the decay (via the \(e^{-\beta\tau}\) temporal factor), leaving \(\Psi\) to evolve under pure diffusion.

**Proposition 4 (Drift removal by exponential tilt).** *If \(U\) satisfies \(U_\tau = D\,U_{xx} + \nu\,U_x\), then \(\tilde{U}(x,\tau) := e^{\nu x/(2D) + \nu^2\tau/(4D)}\,U(x,\tau)\) satisfies \(\tilde{U}_\tau = D\,\tilde{U}_{xx}\).*

*Proof sketch.* This is a special case of Proposition 2 with the decay term absent. Set \(\alpha = -\nu/(2D) = -\nu/\sigma^2\) and \(\beta = -\nu^2/(4D) = -\nu^2/(2\sigma^2)\) (note the sign difference from our earlier convention because here we write \(\tilde{U} = e^{+\alpha' x + \beta'\tau}\,U\)). The coefficient-matching proceeds identically: the exponential tilt generates exactly the cross-terms needed to cancel \(\nu\,U_x\). \(\square\)

This exponential-tilt technique appears throughout applied mathematics. In machine learning, it is closely related to the **exponential family** change-of-measure trick; in large-deviations theory, it is the **Esscher transform**; in physics, it is a **gauge transformation**. All describe the same algebraic mechanism: an exponential change of dependent variable to shift the "natural" parameter of a distribution.

---

## Appendix E — Summary for the Computationally Minded

For CS graduate students, the deepest lesson of this correspondence may be computational rather than conceptual. Here we draw explicit connections to algorithms and numerical methods.

**Finite-difference PDE solvers.** The heat equation \(\Psi_\tau = \frac{\sigma^2}{2}\Psi_{xx}\) is the prototypical PDE for which finite-difference schemes (explicit Euler, implicit Euler, Crank-Nicolson) are taught. Pricing options numerically amounts to running such a scheme on a uniform grid in \((x, \tau)\)-space, then undoing the gauge and log-price transforms. The Crank-Nicolson scheme, second-order accurate in both \(x\) and \(\tau\), is the workhorse of production option-pricing code.

**Monte Carlo and path integrals.** The Feynman–Kac theorem formalises the link between the PDE and an expectation over random paths. In practice, one simulates \(N\) independent Brownian paths of \(x(\tau)\), evaluates the payoff at the end of each path, multiplies by the gauge factor, and averages. This is Monte Carlo option pricing — and it is precisely a discretised path integral with the Wiener measure. Variance-reduction techniques (antithetic variates, importance sampling) are the numerical analyst's version of choosing a clever gauge.

**Spectral methods and quantum computing.** The heat-equation propagator \(e^{D\tau\,\partial_{xx}}\) diagonalises in Fourier space: each mode \(e^{ikx}\) simply decays as \(e^{-Dk^2\tau}\). This is exploited by FFT-based option pricers. The Schrödinger propagator, by contrast, has purely imaginary eigenvalues \(e^{-iDk^2 t}\) — modes rotate rather than decay. Quantum algorithms for option pricing (e.g.\ amplitude estimation) exploit this unitary structure directly. The BS–Schrödinger dictionary is not just a curiosity; it is the mathematical backbone of quantum-finance algorithms.

The table below maps each computational method to its role in the two theories:

```
┌─────────────────────────┬────────────────────────────┬─────────────────────────────┐
│  NUMERICAL METHOD       │  FINANCE APPLICATION       │  QM / PHYSICS APPLICATION   │
├─────────────────────────┼────────────────────────────┼─────────────────────────────┤
│  Finite differences     │  PDE option pricing        │  TDSE on a lattice          │
│  (Crank-Nicolson)       │  (log-price grid)          │  (position-space grid)      │
│                         │                            │                             │
│  Monte Carlo sampling   │  Path-integral pricing     │  Lattice QCD,               │
│  of Brownian paths      │  (Feynman-Kac)             │  quantum Monte Carlo        │
│                         │                            │                             │
│  FFT / spectral methods │  Fourier option pricing    │  Momentum-space propagation │
│                         │  (Carr-Madan)              │  (plane-wave expansion)     │
│                         │                            │                             │
│  Quantum amplitude      │  Quadratic speedup for     │  Native: time evolution     │
│  estimation             │  MC option pricing         │  of quantum states          │
└─────────────────────────┴────────────────────────────┴─────────────────────────────┘
```

---

## Appendix F — Annotated Further Reading

For readers who wish to go deeper, the following directions are natural next steps. Paul Wilmott's *Paul Wilmott on Quantitative Finance* develops the PDE approach to derivatives pricing with great care and remains the standard practitioner reference. For the physics side, Feynman and Hibbs's *Quantum Mechanics and Path Integrals* is the original source on path integrals and is written with remarkable clarity. The connection between the two is developed rigorously in Jean-Philippe Bouchaud and Marc Potters's *Theory of Financial Risk and Derivative Pricing*, which builds derivatives theory from the statistical-physics perspective. Finally, for quantum computing applications to finance, a good entry point is the 2019 paper by Stamatopoulos et al.\ on option pricing using quantum amplitude estimation.

---

## Final Remark

The conceptual arc of this tutorial can be compressed into three sentences. The Black-Scholes PDE, once stripped of its variable coefficients by the log-price substitution and of its drift and decay by a gauge transform, is the one-dimensional heat equation with diffusion coefficient \(\sigma^2/2\). The heat equation is the Wick-rotated (imaginary-time) Schrödinger equation, with the correspondence \(\sigma^2 \leftrightarrow \hbar/m\). Everything else — delta hedging, risk-neutral pricing, Feynman path integrals, option Greeks, quantum propagation — is a corollary of this single structural identity, viewed from one side or the other of the \(i\).
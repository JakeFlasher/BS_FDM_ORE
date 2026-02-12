## Phase 1 — Critical Analysis & Relationship Mapping

<planning>

- Fix a consistent PDE/time convention (time-to-maturity) and rewrite each seed paper’s schemes in that convention.
- For each seed paper: (i) mathematical formulation, (ii) discretizations, (iii) claimed theory, (iv) what is still right/useful, (v) what is incomplete/outdated, (vi) identify contradictions and any shaky claims.
- For the “literature review” corpus: enumerate the papers that are explicitly present in your provided documents + the canonical older works they cite (Pooley–Forsyth–Vetzal; Giles–Carter; Wade et al.; etc.) + the major monotone/fitted FVM line (Song Wang; Valkov).
- Build an explicit relationship graph: failure modes ↔ remedies ↔ papers.
</planning>

<response>

### 1.0 Unifying notation (used throughout Phases 1–3)

Let \(S>0\) be spot, \(K>0\) strike, \(T>0\) maturity (in years), \(r\ge 0\) risk-free rate, \(q\ge 0\) dividend yield, \(\sigma\) volatility.

I will use **time-to-maturity** \(\tau \in [0,T]\) with \(\tau=0\) at maturity (payoff time) and \(\tau=T\) at valuation time (“today”). This matches the “\(-V_t+\cdots=0\)” convention in your seed papers.

A broad one-factor linear pricing PDE (Black–Scholes / local vol class) is:
\[
-\frac{\partial V}{\partial \tau}(S,\tau)
+ \frac{1}{2}\sigma(S,\tau)^2 S^2 \frac{\partial^2 V}{\partial S^2}(S,\tau)
+ (r-q)S \frac{\partial V}{\partial S}(S,\tau)
- r V(S,\tau)=0,
\quad S\in(0,\infty),\ \tau\in(0,T].
\]
Payoff is the “initial” condition:
\[
V(S,0)=\Phi(S).
\]

When I switch to **log-price** \(x=\ln(S/K)\) and \(u(x,\tau)=V(Ke^x,\tau)\), the PDE becomes:
\[
-u_\tau + a(x,\tau) u_{xx} + b(x,\tau) u_x - r u = 0,
\]
with
\[
a(x,\tau)=\tfrac12 \sigma(Ke^x,\tau)^2,\qquad
b(x,\tau)=(r-q)-\tfrac12\sigma(Ke^x,\tau)^2.
\]
This is exactly the “convection–diffusion–reaction” form used implicitly in [Duffy, 2004] and explicitly in [Milev–Tagliani, 2010b].

---

## 1.1 Seed Paper 1 — Duffy (2004, *Wilmott Magazine*): “A Critique of the Crank Nicolson Scheme…”

### 1.1.1 Problem formulation (as used in the article)

Duffy’s core object is the one-factor Black–Scholes PDE (in time-to-maturity form):
\[
- V_\tau + \tfrac12\sigma^2 S^2 V_{SS} + rS V_S - rV = 0,
\]
with appropriate boundary and payoff data. He then broadens to a generic convection–diffusion–reaction IBVP on a bounded interval \((A,B)\):
\[
- u_\tau + \sigma(x,\tau) u_{xx} + \mu(x,\tau) u_x + b(x,\tau)u = f(x,\tau),
\]
with Dirichlet (and discussion of Neumann/Robin) data and initial condition \(u(x,0)=\varphi(x)\). (He uses \(\sigma\) here for diffusion coefficient; I’ll call it \(a\) later to avoid confusion with volatility.)

### 1.1.2 Discretizations discussed

#### (i) Fully implicit (Backward Euler in time, centered in space)

In his notation (uniform mesh, \(h=\Delta S\), \(k=\Delta\tau\)):
- centered \(V_S\), centered \(V_{SS}\),
- backward difference in time at \(\tau_{n+1}\).

This yields a tridiagonal linear system each step.

**What is mathematically correct here:**  
Backward Euler applied to a **monotone spatial discretization** yields an M-matrix linear system and is unconditionally stable in maximum norm; it also provides strong damping of high frequencies (L-stability in the ODE sense).

**What is subtly incomplete:**  
If the *spatial* discretization is centered in a convection-dominated regime, the resulting matrix may fail to be an M-matrix (loss of inverse-positivity), which can cause oscillations/negativity despite implicit time stepping. This is exactly the pathology emphasized later by Milev–Tagliani in low volatility (see Seed Paper 2 and [Milev–Tagliani, 2013 JCAM]).

So Duffy’s statement “fully implicit has no spurious oscillations” is **not universally true** unless one adds hypotheses ensuring the system matrix is an M-matrix (e.g., suitable upwinding/fitting or a Péclet-type bound).

#### (ii) Crank–Nicolson (CN)

He writes CN as the trapezoidal rule in time (average \(n\) and \(n+1\) levels) with centered space stencils. Again yields tridiagonal solve per step.

He correctly emphasizes:

- CN is **A-stable** for linear constant-coefficient diffusion, but **not L-stable** (weak damping of high-frequency modes).
- CN is not monotone; spurious oscillations occur for convection-dominated regimes and/or nonsmooth data.

This is consistent with the classical singular perturbation literature (Il’in-type results) and with the finance PDE literature on payoff nonsmoothness.

#### (iii) Exponentially fitted finite differences (“Il’in-type”)

Duffy introduces the fitting factor (for 1D convection–diffusion ODE first, then PDE) of the form
\[
\rho = \frac{\mu h}{2a}\coth\!\left(\frac{\mu h}{2a}\right),
\]
where \(a\) is diffusion coefficient, \(\mu\) convection coefficient.

He then replaces the usual central second derivative by a scaled second difference \( \rho D^+D^- u\), retaining a centered first derivative, producing a tridiagonal matrix with an M-matrix structure under appropriate sign assumptions. He proves stability via a discrete maximum principle (using M-matrix inverse positivity).

Key claims:

- **Uniform stability** and an error bound of the form \(|u(x_j,\tau_n)-U_j^n|\le M(h+k)\) with \(M\) independent of the small diffusion parameter (i.e., “uniformly convergent” in singular perturbation sense).

This line is standard in the Il’in/“fitted operator” theory **for scalar linear convection–diffusion** under appropriate coefficient assumptions and boundary layers.

#### (iv) Keller box (“box scheme”) for price + delta

Duffy sketches Keller’s box scheme: reduce the second-order PDE to a first-order system and discretize on a staggered box with averaging. He asserts:

- second-order accuracy for both solution and its first derivative (delta),
- unconditional stability,
- better behavior for nonsmooth data.

The box scheme is indeed widely used in convection–diffusion contexts; its stability/accuracy depend on details (choice of staggering, averaging, and coefficient handling). As presented in the article, it is a **high-level sketch**, not a complete proof for BS specifically.

### 1.1.3 Theoretical results claimed and what is justified

#### What is solid

1. **CN oscillations are not captured by von Neumann stability**:  
Correct. CN may be stable in an \(L^2\)-type sense but violate a discrete maximum principle / positivity, leading to nonphysical oscillations.

2. **Monotonicity/M-matrix reasoning is the right lens for financial “no-arbitrage” constraints**:  
Correct and still best practice. Positivity and monotone dependence on data/parameters are strongly connected to M-matrix structure.

3. **Exponentially fitted schemes “degrade” to upwind as diffusion \(\to 0\)**:  
Correct and essential: as \(a\to 0\), fitted schemes produce stable upwind-like behavior.

#### What is overstated / needs conditions

1. “Fully implicit produces no spurious oscillations”  
False without additional conditions on spatial discretization (or parameter regime). Seed Paper 2 and [Milev–Tagliani, 2013] provide counterexamples when centered convection is used and \(\sigma^2 \ll r\).

2. “Second-order accuracy is lost on nonuniform meshes”  
Too blunt. One can construct second-order schemes on nonuniform meshes (e.g., using consistent second-order stencils). What does happen in practice is:

- the **payoff nonsmoothness** dominates early-time error,
- naive nonuniform-mesh discretizations can reduce global order if not carefully designed,
- and discrete monitoring introduces repeated irregularities.

But “nonuniform mesh ⇒ order loss” is not a theorem; it’s a warning about careless discretization.

3. Uniform error bound \(|u-U|\le M(h+k)\) with \(M\) independent of diffusion parameter  
This is plausible within Il’in theory, but in finance PDEs one must be careful about:

- degeneracy at \(S=0\) in \(S\)-coordinates,
- boundary truncation,
- piecewise smooth payoff with kinks (derivative discontinuities),
- and choice of norm (max norm vs \(L^2\)).

So: the **spirit** is right (fitted methods can be uniformly convergent), but applying the theorem “as-is” to every option PDE variant requires verifying assumptions.

### 1.1.4 Strengths (still valid and useful today)

- Clear identification of **failure modes**:
  - convection-dominated regimes (small diffusion),
  - nonsmooth payoffs and incompatibility at corners,
  - discrete monitoring discontinuities,
  - Greeks magnify oscillations.
- Correct prioritization of **monotonicity/positivity** as “financial correctness” constraints.
- Exponentially fitted schemes are still among the most robust spatial discretizations for 1D convection–diffusion.
- Keller box viewpoint (“solve for price and derivative”) remains highly relevant for high-quality Greeks.
- Practical engineering mindset: tridiagonal solves, boundary condition handling, and algorithmic robustness.

### 1.1.5 Weaknesses / what has been superseded

- Lacks the later, finance-specific rigorous convergence analysis of damping for nonsmooth payoffs (e.g., [Giles–Carter, 2006] ([risk.net](https://www.risk.net/journal-of-computational-finance/2160349/convergence-analysis-of-crank-nicolson-and-rannacher-time-marching))).
- Does not incorporate (or predates wide adoption of) **Rannacher smoothing as a default production practice** (though it discusses the need for damping conceptually).
- Does not quantify trade-offs between fitted/upwind diffusion and accuracy in low volatility; Seed Paper 2 does.
- Does not provide a complete modern “best practice pipeline” (payoff smoothing + damping + monotone spatial operator + boundary truncation error control + Greek post-processing).

### 1.1.6 Errors / ambiguities worth flagging

- **Notation collision**: uses \(\sigma\) both as volatility and diffusion coefficient in different sections (common in expository pieces, but easy to misapply).
- **Implicit scheme positivity**: claims of oscillation-free behavior should be qualified by M-matrix conditions (off-diagonal sign + diagonal dominance). Later literature shows this explicitly.

---

## 1.1 Seed Paper 2 — Milev & Tagliani (2010b, *Serdica Math. J.* 36, 223–236): “Low Volatility Options and Numerical Diffusion…”

### 1.1.1 Problem formulation

They consider options with:

- **discontinuous payoff** (digital-like truncations),
- **discrete monitoring** (barriers re-applied at monitoring dates),
- **low volatility** \(\sigma^2 \ll r\), making the PDE effectively convection-dominated.

They use the Black–Scholes PDE in time-to-maturity form:
\[
- V_\tau + rS V_S + \tfrac12\sigma^2 S^2 V_{SS} - rV = 0,
\]
and emphasize that oscillations are tied to the spectrum of the iteration matrix: complex or negative eigenvalues near \(-1\) produce slowly decaying oscillatory modes.

### 1.1.2 Schemes surveyed and analyzed

They focus on two “nonstandard, non-oscillating” approaches:

#### (i) Duffy’s implicitly time-stepped exponentially fitted scheme

They restate Duffy’s fitted operator:
\[
- \frac{U_j^{n+1}-U_j^n}{k}
+ \mu_j^{n+1}\frac{U_{j+1}^{n+1}-U_{j-1}^{n+1}}{2h}
+ \rho_j^{n+1}\,a_j^{n+1}\,\frac{U_{j+1}^{n+1}-2U_j^{n+1}+U_{j-1}^{n+1}}{h^2}
+ b_j^{n+1}U_j^{n+1}=0,
\]
with fitting factor
\[
\rho_j^{n+1}=\frac{\mu_j^{n+1}h}{2a_j^{n+1}}
\coth\!\Bigl(\frac{\mu_j^{n+1}h}{2a_j^{n+1}}\Bigr).
\]

They use M-matrix arguments to claim positivity and a discrete maximum principle, consistent with Duffy.

**New contribution in this paper:** they explicitly track the **numerical diffusion** introduced in the low-volatility regime.

When \(\sigma \to 0\), the fitted scheme tends to implicit upwind. They cite the modified equation interpretation: a first-order upwind discretization introduces an artificial diffusion term proportional to \( \tfrac12 \mu h\, V_{SS}\).

For Black–Scholes in \(S\)-coordinates, \(\mu=rS\), so artificial diffusion coefficient behaves like:
\[
a_{\text{art}}(S) \sim \tfrac12 rS\,\Delta S.
\]
This can dominate the physical diffusion \( \tfrac12 \sigma^2 S^2 \) when \(\sigma\) is tiny and/or \(S\) is large relative to \(\Delta S\), causing smearing.

#### (ii) “Crank–Nicolson variant” (Milev–Tagliani)

They present a modified CN where the reaction term \(-rV\) is discretized using a **nonlocal, mixed time-level stencil** (six nodes) with weights \(\omega_1,\omega_2\), leading to a two-level linear system:
\[
P U^{n+1} = N U^n.
\]

They choose parameters (they state)
\[
\omega_1=\omega_2=-\frac{r}{16\sigma^2},
\]
with a time-step constraint (their formula) that ensures:

- \(P\) is an M-matrix \(\Rightarrow P^{-1}\ge 0\),
- \(N\ge 0\),
- hence positivity of \(U^{n+1}=P^{-1}N U^n\).

They then estimate a discrete maximum principle using \(\|P^{-1}\|_\infty\|N\|_\infty\le 1\) type bounds.

**Key additional claim:** this CN-variant, while oscillation-free, introduces artificial diffusion of size
\[
a_{\text{art}} \sim \frac18\left(\frac{r}{\sigma}\Delta S\right)^2,
\]
again potentially large when \(\sigma\) is tiny.

### 1.1.3 Strengths (still valid today)

- Makes explicit the **accuracy cost** of monotonicity in convection-dominated finance PDEs: to kill oscillations you often pay with numerical diffusion unless the mesh resolves boundary layers.
- Connects oscillations to **spectral properties** (eigenvalues near \(-1\), complex pairs) and to **M-matrix** structure.
- Emphasizes that **discrete monitoring** reintroduces irregularities repeatedly; therefore “smooth the initial payoff and you’re done” is insufficient for discretely monitored exotics.

### 1.1.4 Weaknesses / limitations

- The “artificial diffusion” expressions are derived via modified-equation heuristics that assume locally constant coefficients and a simple regime; they are good diagnostics but not a complete error theory on nonuniform meshes and with boundary truncation.
- They do not integrate the now-standard **Rannacher smoothing / payoff averaging / grid shifting** toolkit that restores second-order convergence for Greeks (those are in [Pooley–Forsyth–Vetzal, 2003] and [Giles–Carter, 2006]).
- Their statement “Crank–Nicolson is useless” in low volatility is **too strong**: CN can be made production-grade with damping and proper spatial monotone discretization, but *raw CN with centered space stencils* is indeed unreliable.

### 1.1.5 Internal contradictions / ambiguous claims

This paper says the fitted scheme is uniformly convergent (per Duffy/Il’in), yet also says accurate solutions require “very small \(\Delta S\)” in low volatility, effectively “losing peculiarity”.

This is **not a logical contradiction** if interpreted correctly:

- **Uniform convergence**: error bound does not blow up as \(\sigma\to 0\) *for fixed \(h\)* in the asymptotic analysis class.
- **Practical accuracy**: the option solution has increasingly sharp layers (or the artificial diffusion dominates) when \(\sigma\) is tiny, so \(h\) must still be small to get a small absolute error.

So the paper is highlighting a *practitioner’s* point: “uniform convergence” does not mean “good accuracy on coarse meshes in extreme regimes.”

### 1.1.6 What I would correct or qualify

- The CN-variant parameter \(\omega=-r/(16\sigma^2)\) grows like \(\sigma^{-2}\); in extremely small \(\sigma\), this produces very large coupling coefficients and forces severe \(\Delta t\) constraints to keep the RHS nonnegative. In practice, such parameter scalings can be numerically awkward (conditioning, cancellation).
- The comparison of artificial diffusion terms is done in \(S\)-coordinates; in **log-price coordinates** the balance changes and is usually more benign. Modern best practice almost always prices in log-space for exactly this reason.

---

# 1.2 Literature review corpus — contributions, relation to seeds, and whether to incorporate

Below I treat as “the literature review” the set of papers and industrial references that appear in your provided documents + the canonical references the seed papers explicitly call out.

For each item:

- **Contribution (3–5 sentences)**,
- **Relationship to seed papers** (extend / correct / validate / contradict),
- **Best-practice inclusion?** (Yes/No + why).

---

## (L1) Pooley, Forsyth, Vetzal (2003) — “Convergence remedies for non-smooth payoffs in option pricing” (*J. Comput. Finance*)

**Contribution.** Systematic study of how payoff discontinuities (or derivative discontinuities) spoil numerical convergence, especially for Greeks. Proposes three spatial “data fixes” (cell averaging/initial averaging, grid shifting, projection) and shows they must be paired with a **special time-stepping method** (Rannacher-type damping) to restore expected convergence. Covers 1D and 2D examples. ([risk.net](https://www.risk.net/journal-of-computational-finance/2160518/convergence-remedies-for-non-smooth-payoffs-in-option-pricing))

**Relation to seeds.** Directly addresses Duffy’s “nonsmooth initial data” criticism with concrete remedies and a demonstrably effective combined strategy. Also answers Milev–Tagliani’s “CN useless” claim by showing that properly damped CN + data treatment is high-accuracy.

**Include in best practice?** **Yes — essential.** This is a cornerstone for production PDE solvers when Greeks matter.

---

## (L2) Giles & Carter (2006) — “Convergence analysis of Crank–Nicolson and Rannacher time-marching” (*J. Comput. Finance*)

**Contribution.** Provides sharp convergence analysis for CN and Rannacher startup in the presence of nonsmooth initial data (Fourier-asymptotic analysis). Shows how many damping steps are needed and how derivative approximations (delta/gamma) regain second-order convergence. ([risk.net](https://www.risk.net/journal-of-computational-finance/2160349/convergence-analysis-of-crank-nicolson-and-rannacher-time-marching))

**Relation to seeds.** Rigorous underpinning for Duffy’s qualitative critique and for the widespread “Rannacher then CN” practice that Seed Paper 2 largely omits. It also clarifies when small \(\Delta t\) helps or hurts.

**Include in best practice?** **Yes — essential** for justifying Rannacher smoothing and choosing damping steps.

---

## (L3) Wade, Khaliq, Yousuf, Vigo-Aguiar, Deininger (2007) — “On smoothing of the Crank–Nicolson scheme … for pricing barrier options” (*JCAM*)

**Contribution.** Shows that standard Luskin–Rannacher smoothing must be modified for **discrete barrier** options because discontinuities recur at monitoring dates. Proposes an improved smoothing strategy and extends it to higher-order time stepping using Padé schemes with damping via subdiagonal Padé steps. ([pure.kfupm.edu.sa](https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/))

**Relation to seeds.** Directly addresses the discrete-monitoring pathology emphasized in Seed Paper 2 and provides a concrete fix within the CN family (contradicting “CN useless” by making it usable).

**Include in best practice?** **Yes** if you price discretely monitored barriers with higher-order time stepping; at minimum, its lesson (“restart damping after each monitoring discontinuity”) should be incorporated.

---

## (L4) Milev & Tagliani (2010a) — “Nonstandard finite difference schemes with application to finance: option pricing” (*Serdica Math. J.*)

*(From your earlier document `non_standard_FDM.md`.)*

**Contribution.** Analyzes why CN oscillates for discretely monitored double barrier knock-out options, especially when \(\sigma^2<r\), via eigenvalue arguments. Proposes a **semi-implicit nonstandard scheme** that sacrifices time order to enforce positivity and damping by constructing \(P\) as an M-matrix and \(N\ge 0\), with a parameter choice \(b=-M/2\). Provides numerical comparisons vs CN/implicit/Duffy/Monte Carlo in a barrier context.

**Relation to seeds.** Extends Seed Paper 2’s agenda (positivity + damping) from a diffusion-diagnosis paper into an actual method with matrix-structure guarantees. Aligns with Duffy’s M-matrix philosophy, but proposes a different nonlocal reaction discretization.

**Include in best practice?** **Maybe (situational).** It is valuable as a **robust “smoother step”** after monitoring resets (similar intent to Rannacher), but it is low order in time and introduces diffusion; I would not use it as the primary high-accuracy marcher.

---

## (L5) Milev & Tagliani (2013) — “Efficient implicit scheme with positivity preserving and smoothing properties” (*JCAM*, 243, 1–9)

*(From your `implicit_scheme.md`.)*

**Contribution.** Shows that even fully implicit centered schemes can oscillate in low-volatility regimes because the discretization matrix can lose M-matrix structure and develop complex eigenvalues, making damping slow. Proposes a **modified implicit scheme** altering the discretization of \(-rV\) via a bivariate/nonlocal approximation, with parameter constraints ensuring M-matrix, positivity, max-norm contraction, and eigenvalues in \((0,1)\).

**Relation to seeds.** This is a direct refinement/correction of Duffy’s “implicit is safe” narrative and a technical deepening of Seed Paper 2’s low-volatility focus. It provides a principled alternative to exponential fitting when diffusion is extremely small.

**Include in best practice?** **Yes (as a damping/smoothing module).** Excellent for “Rannacher-like” smoothing steps after payoff/monitor resets when \(\sigma^2\ll r\).

---

## (L6) Tagliani & Milev (2013) — “Laplace Transform and finite difference methods for the Black–Scholes equation” (*Appl. Math. Comput.*, 220, 649–658)

*(From your `laplacian_FDM_BS.md`.)*

**Contribution.** Applies Laplace transform in time and Post–Widder inversion to convert the PDE into a family of ODE solves in \(S\), using an M-matrix spatial discretization to preserve positivity/complete monotonicity. Proves equivalence (under \(\lambda=1/\Delta t\)) to repeated implicit Euler steps. Shows how wrong convection discretization breaks M-matrix and produces oscillations; quantifies numerical diffusion for upwinding.

**Relation to seeds.** Strongly aligned with Duffy’s discrete maximum principle mindset and Seed Paper 2’s “M-matrix or oscillations” message. It also reframes implicit Euler as a Laplace inversion, offering an alternative viewpoint.

**Include in best practice?** **Usually no for production**, because Post–Widder inversion converges slowly and is numerically delicate; but the **M-matrix insight** and its equivalence-to-implicit-Euler interpretation are useful.

---

## (L7) Goll, Rannacher, Wollner (2015, JCF 2015; preprint 2013) — “Damped Crank–Nicolson time-marching … adaptive solution of BS”

*(From your `damped_CN_Time_marching.md`.)*

**Contribution.** Develops a **space–time finite element** method with damped CN (Euler substeps on selected intervals) and **dual-weighted residual (DWR)** goal-oriented adaptivity. Key result: to get reliable goal error estimates, the **dual problem must be damped consistently** with the primal; otherwise the estimator can have wrong sign/poor effectivity. Provides extensive numerical evidence for price and Delta goals.

**Relation to seeds.** This is essentially the “fully modern” mathematical answer to Duffy’s critique: yes, damped CN is robust—but you must incorporate *adjoint consistency* if you adapt for Greeks/targets. It also formalizes Rannacher damping in a Galerkin-in-time view.

**Include in best practice?** **Yes conceptually** (damping logic, dual consistency, adaptivity principles). The full FE+DWR machinery is beyond many desk pricers, but the damping+goal-error viewpoint is extremely valuable.

---

## (L8) Umeorah & Mashele (2019) — “Crank–Nicolson … rebate barrier option prices” (*Cogent Economics & Finance*)

*(From your `CN_FDM_rebate_barrier_options.md`.)*

**Contribution.** Implements CN for down-and-out rebate barrier options, compares against closed-form for continuously monitored barriers and Monte Carlo, and studies Greeks. Observes spurious oscillations in Greeks near discontinuities and proposes timestep restrictions motivated by positivity/M-matrix conditions.

**Relation to seeds.** This is an applied demonstration of Duffy’s warnings and Seed Paper 2’s oscillation story. However, it mostly stays within the “restrict \(\Delta t\)” approach rather than adopting the higher-quality toolkit (payoff smoothing + Rannacher + monotone spatial discretization).

**Include in best practice?** **No as a method blueprint**, yes as a cautionary case study. The “CN + timestep restriction” approach is not the most robust route for exotics/Greeks.

---

## (L9) Song Wang (2004) — “A novel fitted finite volume method for the Black–Scholes equation” (*IMA J. Numer. Anal.* 24(4), 699–720)

**Contribution.** Constructs a fitted finite volume spatial discretization for the (degenerate) Black–Scholes PDE, proves stability and an error bound for the spatial discretization, and proves the discrete system matrix is an **M-matrix**, yielding a discrete maximum principle. Formulates the method as a Petrov–Galerkin FEM with locally defined basis functions. ([academic.oup.com](https://academic.oup.com/imajna/article/24/4/699/687386))

**Relation to seeds.** This is the “fully rigorous” counterpart to Duffy’s fitted-scheme advocacy, specialized to Black–Scholes degeneracy and finite volume structure.

**Include in best practice?** **Yes — highly recommended** as the spatial backbone if you want provable monotonicity/positivity with good accuracy.

---

## (L10) Valkov (2012) — “Fitted Finite Volume Method for a Generalized Black–Scholes Equation…” (arXiv:1211.1903)

**Contribution.** Transforms the semi-infinite domain to a finite interval and develops a fitted finite volume element approximation for a generalized BS equation that degenerates at both ends. Proves \(\theta\)-weighted full discretization is uniquely solvable and positivity-preserving, with numerical experiments. ([arxiv.org](https://arxiv.org/abs/1211.1903))

**Relation to seeds.** Supports and generalizes the fitted/M-matrix paradigm for degenerate BS operators, strengthening Duffy’s narrative with modern FE/FV analysis.

**Include in best practice?** **Yes (ideas)**: domain transformation + fitted FV positivity. Implementation is more involved but valuable.

---

## (L11) Cont & Voltchkova (2005) — “A finite difference scheme for option pricing in jump diffusion and exponential Lévy models” (*SIAM J. Numer. Anal.* 43(4), 1596–1626)

**Contribution.** Develops FD methods for PIDEs with possibly singular kernels, with localization error estimates and an IMEX-type time stepping. Proves stability/convergence and includes nonsmooth initial condition tests. ([tse-fr.eu](https://www.tse-fr.eu/articles/finite-difference-scheme-option-pricing-jump-diffusion-and-exponential-levy-models))

**Relation to seeds.** Extends the “monotone/stable discretization” agenda beyond pure diffusion PDEs into PIDE territory, where CN-like weak damping and nonsmooth payoffs are also problematic.

**Include in best practice?** **Yes as an extension module** (jump models). Not core for pure BS, but essential once jumps are introduced.

---

## (L12) Reisinger & Whitley (2012) — “Impact of a natural time change on convergence of CN” (arXiv:1210.5487)

**Contribution.** Shows CN can diverge for Dirac initial data under certain mesh ratios, and that a **square-root time change** restores convergence; demonstrates quadratic convergence for price/delta/gamma for European and American options *without Rannacher steps* when the mesh ratio is below a threshold. ([arxiv.org](https://arxiv.org/abs/1210.5487))

**Relation to seeds.** A direct, mathematically sharp alternative to Duffy’s “use damping” solution: change variables so CN behaves well.

**Include in best practice?** **Yes (optional advanced technique)**. It’s powerful, but requires careful remeshing/time-grid design and is less standard in production than Rannacher.

---

## (L13) Christara & Dang (2011) — “Adaptive and high-order methods for valuing American options” (*J. Comput. Finance*)

**Contribution.** Develops space–time adaptive and high-order methods for American options; handles LCP via penalty; uses adaptive grid-point distribution (equidistribution) and a timestep selector; reports high accuracy for prices/Greeks/free boundary. ([risk.net](https://www.risk.net/journal-of-computational-finance/2160418/adaptive-and-high-order-methods-for-valuing-american-options))

**Relation to seeds.** Complements Duffy/MT by showing that “robustness” can be achieved with adaptivity and high-order spatial discretizations, not only with heavy damping.

**Include in best practice?** **Yes (for American options / accuracy-critical contexts)**; some ideas (error equidistribution) also help European barrier pricing.

---

## (L14) Reisinger & Witte (2012) — “Policy iteration as an easy way of pricing American options” (*SIAM J. Financial Math.*)

**Contribution.** Shows policy iteration is a simple generic solver for discrete LCPs from American options; gives iteration complexity bounds and robustness comparisons vs penalty/PSOR. ([maths.ox.ac.uk](https://www.maths.ox.ac.uk/node/21883))

**Relation to seeds.** Not about CN oscillations per se, but crucial for extending the PDE framework to American exercise while retaining monotone structure.

**Include in best practice?** **Yes** for American options: policy iteration is now a go-to.

---

## (L15) Khalsaraei et al. (2021) — “Qualitatively stable nonstandard FD for nonlinear BS with transaction costs” (*Journal of Mathematics*, Hindawi)

**Contribution.** Proposes a nonstandard explicit update (closed form) for a nonlinear BS PDE with Gamma-dependent volatility (transaction costs), targeting monotonicity, stability, positivity, and viscosity-solution convergence via Barles–Souganidis conditions.

**Relation to seeds.** Extends the “monotone scheme” philosophy from linear convection–diffusion to a nonlinear degenerate parabolic PDE; aligns with the M-matrix/monotonicity ethos, but uses a nonstandard time–space mixing.

**Include in best practice?** **As an extension** if you truly need transaction-cost nonlinear PDEs. For standard BS, not necessary.

---

## (L16) Lopez-Salas et al. (2024) — “IMEX-RK finite volume methods… Application to option pricing” (arXiv:2409.01125)

**Contribution.** Builds 2nd-order IMEX Runge–Kutta FV schemes: advection explicit, diffusion implicit, designed for efficiency and nonsmooth initial data, extendable to higher order. ([arxiv.org](https://arxiv.org/abs/2409.01125))

**Relation to seeds.** Modernizes the time-stepping story: rather than CN vs Euler, uses IMEX-RK with FV reconstructions, addressing stiffness and robustness in a modern MOL framework.

**Include in best practice?** **Potentially yes** if you want FV infrastructure and IMEX time stepping; for plain BS 1D, standard implicit solvers remain simpler, but IMEX becomes valuable in PIDE/nonlinear settings.

---

## (L17) Wang et al. (2025) — “A second-order FD method for BS without far-field BCs” (*Journal of Financial Stability*, 2025)

**Contribution.** Proposes an explicit ADE method on a dynamically shrinking grid to avoid artificial far-field boundary conditions; reports second-order convergence and speed benefits for low-latency applications. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S1572308925001068))

**Relation to seeds.** Directly targets one of Duffy’s “boundary condition” pain points (far-field truncation). It is orthogonal to the CN oscillation issue but relevant to practical robustness.

**Include in best practice?** **Maybe**: interesting for real-time/latency contexts; I would still default to monotone implicit methods for robustness unless boundary handling dominates your error budget.

---

## (L18) In ’t Hout & Foulon (2010) — ADI schemes for Heston PDE (*IJNAM* 7(2))

**Contribution.** Investigates Douglas, Craig–Sneyd, Modified Craig–Sneyd, Hundsdorfer–Verwer ADI schemes for the 2D Heston PDE with mixed derivatives, on nonuniform grids; identifies effective schemes under stability theory. ([global-sci.com](https://global-sci.com/index.php/ijnam/article/view/9955))

**Relation to seeds.** Extends time-stepping robustness from 1D BS to 2D with mixed derivatives. Not central to one-factor best practice, but important context.

**Include in best practice?** **Mention as extension only** (multi-factor out of scope for core framework).

---

## (L19) Düring & Miles (2017) — High-order ADI for stochastic volatility (*JCAM*)

**Contribution.** Fourth-order spatial + second-order time ADI (Hundsdorfer–Verwer) for Heston-like PDEs, demonstrating high-order convergence. ([sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0377042716304678))

**Relation to seeds.** Continues the “avoid CN pathologies and improve accuracy” theme in higher dimension.

**Include in best practice?** **Extension only**.

---

## (L20) Bonaventura & Della Rocca (2015) — TR-BDF2 monotonicity/positivity analysis (arXiv:1510.04303)

**Contribution.** Computes absolute monotonicity radius for TR-BDF2, relates L-stability parameter choice to maximal monotonicity radius, proposes SSP extensions. ([arxiv.org](https://arxiv.org/abs/1510.04303))

**Relation to seeds.** Gives a modern time-stepping alternative to CN with better damping and analyzable monotonicity constraints.

**Include in best practice?** **Yes (as a modern time integrator option)**—especially when you want second order + stronger damping than CN.

---

# 1.3 Relationship graph (papers ↔ failure modes ↔ remedies)

## 1.3.1 Failure modes (nodes)

Let me name the key failure modes discussed across the corpus:

- **F1 — Convection dominance / low volatility**: \(a\) small, \(b\) large ⇒ central differencing oscillations.
- **F2 — Nonsmooth payoff / corner incompatibility**: payoff kink/discontinuity ⇒ order reduction, oscillatory Greeks.
- **F3 — Discrete monitoring resets**: repeated discontinuities in time for barriers.
- **F4 — Loss of positivity / violation of maximum principle**: negative prices/Greeks; non-monotone scheme.
- **F5 — Artificial numerical diffusion**: monotone fixes smear solution in extreme regimes.
- **F6 — Far-field truncation boundary error**: poor \(S_{\max}\) or boundary condition contaminates interior.
- **F7 — American constraint (LCP/HJB)**: need robust complementarity solver.
- **F8 — Greeks instability**: differentiation amplifies noise/oscillations.

## 1.3.2 Remedy classes (nodes)

- **R1 — Monotone spatial discretization**: upwind, exponential fitting, fitted FV (M-matrix).
- **R2 — Damping in time**: Rannacher (Euler substeps), Padé damping, L-stable methods.
- **R3 — Payoff/data regularization**: cell averaging, projection, grid shifting.
- **R4 — Restart smoothing at monitoring dates**.
- **R5 — Robust boundary handling**: truncation error bounds, better far-field BCs, or avoiding far-field BCs.
- **R6 — Greek-aware computation**: Keller box, adjoint/DWR, post-processing.
- **R7 — American solvers**: policy iteration, penalty, splitting.

## 1.3.3 ASCII dependency / influence map

```
                 ┌──────────────────────────────────────────┐
                 │ Seed S1: Duffy (2004)                    │
                 │  F1,F2,F4,F6,F8  →  R1,R2,R5,R6          │
                 └───────────────┬──────────────────────────┘
                                 │ motivates
                                 ▼
┌───────────────────────────────────────────────────────────────────┐
│ Seed S2: Milev–Tagliani (2010b)                                   │
│  Focus: F1,F2,F3,F5 (numerical diffusion trade-off)               │
│  Surveys: Duffy fitted (R1) + MT CN-variant (R1 via nonlocal term)│
└───────────────┬───────────────────────────────────────────────────┘
                │ extends / refines
                ▼
┌───────────────────────────┐     ┌────────────────────────────────┐
│ MT (2010a) Serdica        │     │ MT (2013) JCAM modified implicit│
│ semi-implicit NS scheme   │     │ smoothing M-matrix in σ^2<<r    │
│ F3,F4 → R1,R4             │     │ F1,F4,F8 → R1,R2 (smoother)     │
└─────────────┬─────────────┘     └──────────────┬─────────────────┘
              │                                   │
              │ complements                        │
              ▼                                   ▼
┌───────────────────────────────────────────────────────────────────┐
│ Pooley–Forsyth–Vetzal (2003) + Giles–Carter (2006)                │
│  F2,F8 → R3 + R2 (Rannacher), convergence for Greeks               │
└───────────────┬───────────────────────────────────────────────────┘
                │
                ▼
┌───────────────────────────────────────────────────────────────────┐
│ Wade et al. (2007)                                                │
│  F3 → R4 + (Padé damping)                                          │
└───────────────────────────────────────────────────────────────────┘

Parallel rigorous monotone spatial line:
┌───────────────────────────────────────────────────────────────────┐
│ Song Wang (2004) + Valkov (2012)                                  │
│  F1,F4 → R1 (fitted FV, M-matrix, provable DMP)                     │
└───────────────────────────────────────────────────────────────────┘

Adaptivity / goal error / dual consistency line:
┌───────────────────────────────────────────────────────────────────┐
│ Goll–Rannacher–Wollner (2015)                                      │
│  F2,F8 → R2 + R6 (DWR, dual damping)                                │
└───────────────────────────────────────────────────────────────────┘

American/LCP line:
┌───────────────────────────────────────────────────────────────────┐
│ Christara–Dang (2011), Reisinger–Witte (2012), Gu et al. (2024)    │
│  F7 → R7 (penalty, policy iteration, parallel-in-time)             │
└───────────────────────────────────────────────────────────────────┘

Applied CN baseline / caution:
┌───────────────────────────────────────────────────────────────────┐
│ Umeorah–Mashele (2019)                                             │
│  Demonstrates F2,F3,F8 under CN; uses timestep restriction only    │
└───────────────────────────────────────────────────────────────────┘
```

## 1.3.4 My “relationship verdict” in one sentence

- **Seed papers correctly diagnose**: CN + centered space is not robust under \(F1\)–\(F3\).  
- **Later best practice is a combination**: \(R1\) (monotone spatial) + \(R3\) (payoff/data smoothing) + \(R2/R4\) (damping, restarted at discontinuities) + \(R6\) (Greek-aware post-processing), and for American add \(R7\).

That composite framework is what I will formalize in Phase 3.

</response>

---

## Phase 2 — State-of-the-art survey & gap analysis (with web supplementation)

<planning>

- Identify SOTA modules by the requested categories (2.1–2.5).
- For each module: pick a small set of high-quality, load-bearing references (classic + modern up to 2025/2026), then summarize the key numerical idea and why it matters for the seed-paper failure modes.
- Explicitly connect SOTA methods back to: monotonicity/M-matrix, damping, payoff smoothing, boundary truncation, Greek accuracy, and American LCP solvers.
- End with a candid gap analysis: what remains difficult (provable 2nd-order + positivity + nonsmooth data + repeated discontinuities, etc.).
</planning>

<response>

### 2.1 Robust time-stepping for Black–Scholes and generalizations

#### 2.1.1 Rannacher smoothing (Backward Euler substeps then CN)

- **Core idea.** Replace the first (few) CN steps after a discontinuity by implicit Euler half-steps. This restores second-order convergence for the *solution* and especially for *derivatives* (Greeks) that otherwise suffer order reduction and oscillations.
- **Key references.**
  - [Pooley–Forsyth–Vetzal, 2003] on combining payoff fixes + special time stepping for Greeks. ([risk.net](https://www.risk.net/journal-of-computational-finance/2160518/convergence-remedies-for-non-smooth-payoffs-in-option-pricing))  
  - [Giles–Carter, 2006] gives a sharp convergence analysis and guidance on number of damping steps. ([risk.net](https://www.risk.net/journal-of-computational-finance/2160349/convergence-analysis-of-crank-nicolson-and-rannacher-time-marching))  
  - [Wade et al., 2007] shows Rannacher must be restarted/adapted for discrete barriers. ([pure.kfupm.edu.sa](https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/))

**SOTA status.** Still the default in many production 1D and 2D PDE engines because it’s simple (just change the first few timesteps), robust, and compatible with tridiagonal/banded solves.

---

#### 2.1.2 TR-BDF2 / TR/BDF2 (L-stable second order)

- **Core idea.** Use a composite method (trapezoidal rule stage + BDF2 stage) that is second order and **L-stable** for the right parameter choice, giving much stronger high-frequency damping than CN while retaining second order for smooth problems.
- **Key references.**
  - [Bank et al., 1996] “Analysis and implementation of TR-BDF2” (Applied Numerical Mathematics) gives stability and implementation details. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/0168927495001158))  
  - [Bonaventura–Della Rocca, 2015] analyzes monotonicity/positivity/SSP properties of TR-BDF2. ([arxiv.org](https://arxiv.org/abs/1510.04303))

**Why it matters for seeds.** CN’s weak damping is a root cause of persistent oscillations (Seed S1/S2). TR-BDF2 is a principled “drop-in” second-order replacement with better damping.

**Caveat.** Like CN, **unconditional positivity** is not automatic; monotonicity/positivity typically imposes a timestep restriction unless you use specially designed SSP variants. [Bonaventura–Della Rocca, 2015] makes this precise. ([arxiv.org](https://arxiv.org/abs/1510.04303))

---

#### 2.1.3 IMEX (implicit–explicit) and modern MOL time integrators

- **Core idea.** Treat stiff diffusion implicitly and advection/source explicitly (or vice versa), using IMEX Runge–Kutta to reduce computational cost and relax stability constraints.
- **Finance-relevant reference.**
  - [López‑Salas et al., 2024] develops 2nd-order IMEX-RK FV schemes and explicitly targets option pricing PDEs with possible nonlinearities. ([arxiv.org](https://arxiv.org/abs/2409.01125))

**Practical take.** For **plain 1D BS** you can just solve tridiagonals with fully implicit methods cheaply; IMEX becomes more compelling when you add:
- nonlocal jump integrals (PIDE),
- nonlinearities (transaction costs),
- or higher-dimensional splitting.

---

#### 2.1.4 Adaptive time stepping

There are two “adaptive” stories:

1. **Goal-oriented adaptivity (DWR)** (very rigorous but heavier):  
   [Goll–Rannacher–Wollner, 2015] builds a full estimator framework with damping-consistent dual. ([risk.net](https://www.risk.net/journal-of-computational-finance/2406534/the-damped-crank-nicolson-time-marching-scheme-for-the-adaptive-solution-of-the-black-scholes-equation))

2. **Local truncation error control in FD** (lighter):  
   Many practical PDE engines vary \(k\) to enforce early-time damping and late-time accuracy, often with heuristic controllers. For American options, [“Pricing American options using a space-time adaptive finite difference method”, 2010] explicitly adapts space and time and uses splitting for the constraint. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0378475410000522))

---

### 2.2 Spatial discretization improvements

#### 2.2.1 Exponential fitting / fitted finite volumes (M-matrix preservation)

- **Song Wang (2004)** is a key rigorous reference for fitted finite volume discretization of Black–Scholes, proving M-matrix and DMP. ([academic.oup.com](https://academic.oup.com/imajna/article/24/4/699/687386))  
- **Valkov (2012)** generalizes fitted FV on finite transformed interval and proves positivity with \(\theta\)-time stepping. ([arxiv.org](https://arxiv.org/abs/1211.1903))

**SOTA status.** In 1D, fitted FV / exponential fitting is one of the best “monotone yet accurate” spatial choices when convection dominance is an issue.

**Why it beats naive upwind.** Upwind is monotone but adds diffusion \(O(h)\). Fitted schemes often capture boundary layers much better at comparable \(h\) (though Seed S2 correctly notes they can still be diffusive in extreme regimes unless the mesh resolves the layer).

---

#### 2.2.2 Grid design: optimal/nonuniform grids and strike/barrier clustering

- [Lyu et al., 2021] “Optimal non-uniform finite difference grids for the Black–Scholes equations” proposes a systematic grid-point removal strategy to obtain “optimal” nonuniform grids. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0378475420304596))

In production, the most common “SOTA” practice remains:
- log-price grids,
- piecewise-uniform or smoothly graded clustering near \(S=K\) and barriers,
- and aligning grid nodes with discontinuities whenever possible.

---

#### 2.2.3 hp-FEM / spectral methods (one-factor)

For one-factor European BS, spectral/hp methods can be excellent, but the main blocker is nonsmooth payoff (Gibbs-type issues) unless you:
- project/average payoff,
- or use domain decomposition / filtering,
- or choose basis tailored to kinks.

You already have one spectral-style transaction cost example in your documents: [de Frutos & Gaton, 2021] (pseudospectral for a nonlinear transaction-cost model). ([arxiv.org](https://arxiv.org/abs/2103.05369))

---

### 2.3 Treatment of non-smooth data and boundary conditions

#### 2.3.1 Payoff smoothing (cell averaging, projection, grid shifting)

This is largely “mature tech”:
- [Pooley–Forsyth–Vetzal, 2003] gives the core triad and demonstrates that none alone suffices for Greeks; you need time damping too. ([risk.net](https://www.risk.net/journal-of-computational-finance/2160518/convergence-remedies-for-non-smooth-payoffs-in-option-pricing))  
- [Giles–Carter, 2006] supplies convergence theory for damping. ([risk.net](https://www.risk.net/journal-of-computational-finance/2160349/convergence-analysis-of-crank-nicolson-and-rannacher-time-marching))

**SOTA production pattern:**  
“Cell-average payoff + 2 Euler half-steps + CN (or TR-BDF2) + monotone spatial operator.”

---

#### 2.3.2 Repeated discontinuities (discrete monitoring barriers)

- [Wade et al., 2007] is the key message: restart damping after each monitoring reset; also Pade-based high-order possibilities. ([pure.kfupm.edu.sa](https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/))  
- Milev–Tagliani line (2010–2013) reinforces that **monitoring dates are like multiple ‘initial times’**.

---

#### 2.3.3 Time/coordinate transformations

- **Square-root time change**: [Reisinger–Whitley, 2012] shows a time change can restore CN convergence without Rannacher steps for some option problems. ([arxiv.org](https://arxiv.org/abs/1210.5487))

**My engineering view:** this is powerful but less “standardizable” than Rannacher; still worth knowing, especially if you want clean second-order CN behavior with fewer damping hacks.

---

#### 2.3.4 Far-field boundary conditions / truncation error

Classic issue: choose \(S_{\max}\) large, use asymptotic boundary (Dirichlet for call/put), or use error bounds (e.g., Kangro–Nicolaides type).

A modern twist: **avoid far-field BCs altogether**.
- [Wang et al., 2025] propose an explicit ADE method with a dynamically shrinking grid to remove far-field BC dependence. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S1572308925001068))

This is intriguing for latency and boundary robustness, but it is not yet a standard production replacement for monotone implicit schemes.

---

### 2.4 Greeks computation (Delta, Gamma)

#### 2.4.1 “Differentiate the discrete solution” vs “solve for Greeks”

Two main schools:

1. **Post-processing**:
   - Compute \(V\) accurately, then use nonuniform finite differences / spline derivatives.
   - Needs payoff smoothing + damping to avoid oscillatory Greeks.

2. **Simultaneous systems**:
   - Keller box / mixed formulations compute \(V\) and \(V_S\) directly (Duffy’s suggestion).
   - FE/DWR approach uses dual problem to target Greek goals (Goll–Rannacher–Wollner). ([risk.net](https://www.risk.net/journal-of-computational-finance/2406534/the-damped-crank-nicolson-time-marching-scheme-for-the-adaptive-solution-of-the-black-scholes-equation))

**SOTA lesson:** If Greeks are the objective (risk), treat them as a first-class goal:
- damping tuned for derivative convergence,
- possibly goal-oriented adaptivity (DWR),
- or simultaneous delta schemes.

---

### 2.5 Extensions to realistic models (brief, one-factor emphasis)

Even staying in one factor, there are “beyond BS” directions:

#### 2.5.1 Jump-diffusion PIDEs (one-factor)

- [Cont–Voltchkova, 2005] is a benchmark: FD + IMEX, localization error, convergence/stability theory. ([ideas.repec.org](https://ideas.repec.org/p/hal/journl/halshs-00445645.html))  
- Many later works add FFT acceleration for the integral term (also in your docs: “Pricing options under jump diffusion with fitted FV”, 2008). ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S0096300307012179))

#### 2.5.2 American options (free boundary / LCP)

SOTA solver patterns:
- penalty methods (robust, easy),
- policy iteration (very effective),
- operator splitting (Ikonen–Toivanen line),
- and more recently parallel-in-time LCP solving:
  - [Gu–Liu–Oosterlee, 2024] parallel-in-time policy iteration for American options. ([arxiv.org](https://arxiv.org/abs/2405.08280))

A key “practical theorem” level message:
- [Reisinger–Witte, 2012] policy iteration is extremely simple and robust; iteration count scales linearly with problem size. ([maths.ox.ac.uk](https://www.maths.ox.ac.uk/node/21883))

#### 2.5.3 Nonlinear models (transaction costs / uncertain volatility)

- [Pooley–Forsyth–Vetzal, 2003 IMAJNA] show non-monotone discretizations can converge to wrong solutions in nonlinear uncertain-volatility PDEs. ([academic.oup.com](https://academic.oup.com/imajna/article-abstract/23/2/241/684490))  
- Your NSFD 2021 paper aims at Barles–Souganidis monotone/stable/consistent convergence for a nonlinear transaction-cost PDE.

---

## 2.6 Gap analysis — what remains genuinely hard / open-ish

Here are the real pain points that persist even in 2026:

1. **Second-order in time + unconditional positivity/monotonicity for nonsmooth payoffs**  
   - Implicit Euler is monotone but first order.
   - CN is second order but not monotone; requires damping and/or timestep restrictions for positivity.
   - TR-BDF2 improves damping, but positivity still has constraints (absolute monotonicity radius). ([arxiv.org](https://arxiv.org/abs/1510.04303))  
   **Gap:** a “drop-in, always positive, always second-order, no tuning” method for all exotic boundary/monitoring configurations does not really exist; you must choose trade-offs.

2. **Repeated discontinuities (discrete monitoring) + high-order methods**  
   - Restart smoothing works but introduces method switching and complexity.
   **Gap:** systematic, provably optimal damping/restart strategies for arbitrary monitoring schedules and payoffs remain underdeveloped.

3. **Greeks as primary outputs**  
   - Even with good prices, gamma is delicate near kinks/barriers.
   - Goal-oriented methods (DWR) are powerful but heavy to integrate into production FD codes.
   **Gap:** light-weight, reliable Greek error estimators and adaptivity for 1D FD codes.

4. **Boundary truncation robustness**  
   - Choosing \(S_{\max}\) is still half art; new approaches avoid far-field BCs but are not mainstream. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S1572308925001068))  
   **Gap:** universal boundary treatments that are both theoretically tight and operationally simple.

5. **Nonlinear PDEs (transaction costs, uncertain volatility) in production**  
   - Monotone schemes exist, but the computational burden and calibration/parameter uncertainty limit adoption.
   **Gap:** fast monotone solvers with rigorous error control and robust calibration hooks.

</response>

---

## Phase 3 — Best-practice mathematical framework & algorithm (implementable)

<planning>

- Choose a general one-factor PDE framework (linear BS/local vol, with hooks for barriers and American exercise).
- Commit to a robust default stack:
  - log-price coordinates,
  - fitted/monotone spatial discretization (exponential fitting / fitted FV-inspired),
  - Rannacher-smoothed CN (or TR-BDF2) in time,
  - payoff cell-averaging + compatibility enforcement,
  - restart damping after monitoring discontinuities,
  - Greeks via consistent nonuniform formulas in log-space with error orders.
- Provide explicit discrete operator coefficients, then the fully discrete time-stepping formulas.
- Provide tridiagonal linear algebra (Thomas algorithm) inline.
- Provide full pseudocode with line numbers and complexity annotations.
- Finish with an error/stability comparison table and a pragmatic “just works” configuration + debugging checklist.
</planning>

<response>

# 3.1 Problem formulation (general one-factor framework)

## 3.1.1 Model and PDE

Assume under the risk-neutral measure:
\[
\frac{dS_t}{S_t} = (r-q)\,dt + \sigma(S_t,\tau)\,dW_t,
\]
with \(r,q\ge 0\) constants (extensions to term structures are straightforward).

Let \(V(S,\tau)\) be the option value at spot \(S>0\) and time-to-maturity \(\tau\in[0,T]\), with \(\tau=0\) at maturity payoff. Then \(V\) satisfies:
\[
-\frac{\partial V}{\partial \tau}
+ \frac{1}{2}\sigma(S,\tau)^2 S^2 \frac{\partial^2 V}{\partial S^2}
+ (r-q)S \frac{\partial V}{\partial S}
- r V = 0.
\]
This covers:
- constant-vol Black–Scholes,
- local volatility (Dupire) as long as \(\sigma(S,\tau)\) is given,
- many barrier/rebate options via boundary conditions.

### Hook: source terms / killing
If the product has continuous cashflows or state-dependent killing, add a source term \(f(S,\tau)\) and/or modify reaction:
\[
-\partial_\tau V + a V_{SS} + b V_S + c V = f.
\]

## 3.1.2 Boundary conditions (one-factor, bounded computational domain)

We truncate \(S\in(0,\infty)\) to \(S\in[S_{\min},S_{\max}]\).

- For **vanilla call**:
  \[
  V(S_{\min},\tau)=0,\qquad
  V(S_{\max},\tau)=S_{\max}e^{-q\tau}-K e^{-r\tau}.
  \]
- For **vanilla put**:
  \[
  V(S_{\min},\tau)=K e^{-r\tau},\qquad
  V(S_{\max},\tau)=0.
  \]
- For **down-and-out barrier** at \(B\): set \(S_{\min}=B\) and impose:
  - knock-out (no rebate): \(V(B,\tau)=0\),
  - rebate-at-hit: \(V(B,\tau)=R(\tau)\) (often \(R\) or \(Re^{-r(T-\tau)}\), depending on convention).

### Discrete monitoring (barrier corridor)
For discrete monitoring times \(0<\tau_1<\cdots<\tau_F\le T\), enforce the projection at each \(\tau_i\):
\[
V(S,\tau_i^+) \leftarrow V(S,\tau_i^-)\,\mathbf{1}_{[L,U]}(S),
\]
for a double knock-out corridor \([L,U]\). This is exactly the “renewed discontinuity” mechanism.

## 3.1.3 Initial condition (payoff) and regularity

Payoffs \(\Phi(S)\) are typically only Lipschitz and often have derivative jumps (call/put kink at \(S=K\)); digital/barrier truncations create discontinuities.

**Numerical requirement:** treat \(\Phi\) via **cell averaging / projection** to avoid exciting grid-scale modes that CN cannot damp (Phase 1 results).

## 3.1.4 Coordinate transformation recommendation (log-price)

Use:
\[
x=\ln\left(\frac{S}{K}\right),\quad S=K e^x,\quad u(x,\tau)=V(Ke^x,\tau).
\]

Then on \(x\in[x_{\min},x_{\max}]\), the PDE becomes:
\[
-u_\tau + a(x,\tau)\,u_{xx} + b(x,\tau)\,u_x - r u = 0,
\]
with
\[
a(x,\tau)=\tfrac12\sigma(Ke^x,\tau)^2,\qquad
b(x,\tau)=(r-q)-\tfrac12\sigma(Ke^x,\tau)^2.
\]

**Why this is best-practice:**
- removes \(S^2\) degeneracy from diffusion coefficient,
- gives more uniform resolution around strike when using uniform \(x\)-grid,
- simplifies far-field asymptotics.

---

# 3.2 Spatial discretization — recommended scheme

I recommend a **monotone exponentially fitted finite difference** (Il’in/Duffy style) in \(x\), because:

- it produces an M-matrix in the implicit solve under mild conditions,
- it avoids central-difference oscillations in convection-dominated regimes,
- it is tridiagonal and production-cheap,
- it aligns with the rigorous fitted FV/M-matrix paradigm (Song Wang 2004; Valkov 2012). ([academic.oup.com](https://academic.oup.com/imajna/article/24/4/699/687386))

## 3.2.1 Mesh in \(x\)

### Default (robust and simple): uniform in \(x\)
Choose \(J\in\mathbb{N}\) and
\[
x_j = x_{\min} + jh,\quad h=\frac{x_{\max}-x_{\min}}{J},\quad j=0,\dots,J.
\]

Set
\[
S_j = K e^{x_j}.
\]

### Choosing \([x_{\min},x_{\max}]\)
A “just works” choice (for equity-like parameters) is:
\[
x_{\min} = \ln\!\left(\frac{S_0}{K}\right) - m\,\sigma_{\text{ref}}\sqrt{T},\qquad
x_{\max} = \ln\!\left(\frac{S_0}{K}\right) + m\,\sigma_{\text{ref}}\sqrt{T},
\]
with \(m\in[6,10]\) and \(\sigma_{\text{ref}}\) a representative vol (e.g., spot vol).

For barrier options, set \(x_{\min}=\ln(B/K)\) (down barrier) or \(x_{\max}=\ln(U/K)\) (up barrier) to align a barrier exactly.

**Practical note:** if \(B\) is not exactly on a node, either:
- adjust \(x_{\min}\) to hit it exactly (recommended),
- or enforce the boundary via interpolation (more complex, less clean).

## 3.2.2 Discrete derivatives

For interior nodes \(j=1,\dots,J-1\) with uniform \(h\):
\[
(\delta_x u)_j := \frac{u_{j+1}-u_{j-1}}{2h},\qquad
(\delta_{xx} u)_j := \frac{u_{j+1}-2u_j+u_{j-1}}{h^2}.
\]

## 3.2.3 Exponential fitting factor

Define the local Péclet-like parameter:
\[
\theta_j(\tau) := \frac{b_j(\tau)\,h}{2a_j(\tau)},
\quad\text{where}\quad
a_j(\tau)=a(x_j,\tau),\ b_j(\tau)=b(x_j,\tau).
\]

Define the fitted factor:
\[
\rho_j(\tau) :=
\begin{cases}
\theta_j(\tau)\coth(\theta_j(\tau)), & \theta_j(\tau)\neq 0,\\
1, & \theta_j(\tau)=0.
\end{cases}
\]

Useful limits (for implementation stability):
- for \(|\theta|\ll 1\): \(\rho \approx 1 + \theta^2/3\),
- for \(\theta\to +\infty\): \(\rho\sim \theta\),
- for \(\theta\to -\infty\): \(\rho\sim -\theta\).

This is exactly the “fitting factor” structure in [Duffy, 2004] and in Seed Paper 2.

## 3.2.4 Discrete spatial operator

Define the fitted spatial operator \(L_h(\tau)\) acting on grid values \(u_j\) by:
\[
(L_h u)_j
:= a_j(\tau)\rho_j(\tau)\,(\delta_{xx}u)_j + b_j(\tau)\,(\delta_x u)_j - r u_j.
\]

Expanding into a tridiagonal stencil:
\[
(L_h u)_j
= \ell_j(\tau)\,u_{j-1} + d_j(\tau)\,u_j + u_j^{(+)}(\tau)\,u_{j+1},
\]
with coefficients:
\[
\ell_j(\tau)=\frac{a_j\rho_j}{h^2}-\frac{b_j}{2h},
\qquad
u_j^{(+)}(\tau)=\frac{a_j\rho_j}{h^2}+\frac{b_j}{2h},
\]
\[
d_j(\tau)= -\frac{2a_j\rho_j}{h^2} - r.
\]

### Monotonicity/M-matrix check (key property)

For the **implicit Euler** step (see 3.3.1), the matrix is \(I - k L_h\). Off-diagonals are \(-k\ell_j\) and \(-k u_j^{(+)}\). To have nonpositive off-diagonals (M-matrix pattern), we need:
\[
\ell_j(\tau)\ge 0,\qquad u_j^{(+)}(\tau)\ge 0.
\]
In convection-dominated regimes, central differences would violate this, but the fitted factor \(\rho_j\) increases the effective diffusion to enforce these inequalities for a broad range of parameters (this is the whole point of fitting).

In practice, if local \(\ell_j\) or \(u_j^{(+)}\) become negative due to coefficient extremes, you **must**:
- refine the mesh (reduce \(h\)),
- or switch to a stricter monotone flux (pure upwind) locally.

I will include this as a runtime check in the pseudocode.

## 3.2.5 Truncation error (spatial)

For smooth \(u\), the fitted operator is typically:
- second-order in space for moderate Péclet numbers,
- first-order in space in extreme convection dominance (it behaves like upwind).

A safe statement (aligned with the seed papers’ philosophy) is:

- \(O(h^2)\) in diffusion-dominated regions,
- \(O(h)\) in convection-dominated regions.

In option pricing, this is acceptable because the purpose of fitting is *qualitative correctness* (no oscillations) in difficult regions; high-order accuracy can be recovered via mesh refinement near discontinuities/steep layers and via payoff/time smoothing.

---

# 3.3 Temporal discretization — recommended scheme

I give two time-stepping recommendations:

- **Default production choice:** Rannacher-smoothed Crank–Nicolson (RS-CN).  
- **Modern alternative:** TR-BDF2 with optional startup damping.

Both are tridiagonal-solvable with the spatial operator above.

## 3.3.1 Rannacher-smoothed Crank–Nicolson (RS-CN)

### Time grid
Let \(N\in\mathbb{N}\), \(k=T/N\), \(\tau_n=nk\).

### Rannacher startup (after each discontinuity)
Use:
- two **implicit Euler half-steps** of size \(k/2\),
- then standard CN steps of size \(k\).

This is consistent with [Giles–Carter, 2006] and [Pooley–Forsyth–Vetzal, 2003]. ([risk.net](https://www.risk.net/journal-of-computational-finance/2160349/convergence-analysis-of-crank-nicolson-and-rannacher-time-marching))

### Fully discrete formulas

Let \(U^n\) be the vector of interior unknowns at \(\tau_n\).

**Implicit Euler step** (size \(k_E\)):
\[
\left(I - k_E L_h(\tau_{n+1})\right) U^{n+1} = U^n + \text{BC terms}.
\]

**Crank–Nicolson step** (size \(k\)):
\[
\left(I - \frac{k}{2} L_h(\tau_{n+1})\right) U^{n+1}
=
\left(I + \frac{k}{2} L_h(\tau_{n})\right) U^{n} + \text{BC terms}.
\]

### Stability / monotonicity notes (precise and honest)

- **Stability (linear)**: RS-CN is A-stable (CN part) and strongly damping in the first steps (Euler).
- **Positivity / monotonicity**:
  - Implicit Euler is positivity-preserving if \(I-k_E L_h\) is an M-matrix (true under the fitted spatial operator conditions).
  - CN is not unconditionally positivity-preserving; it can require a timestep restriction to keep the right-hand side nonnegative (this is the Milev–Tagliani maximum principle logic).

**Best-practice reality:** in 1D finance PDEs, practitioners accept that:
- Euler steps enforce early smoothing and reduce oscillation amplitudes,
- CN thereafter yields high accuracy; any slight negativity is typically tiny if payoff is smoothed and the spatial operator is monotone.

If you require strict positivity, use implicit Euler throughout (or a fully monotone time integrator) and accept first-order time accuracy.

## 3.3.2 TR-BDF2 (modern alternative)

TR-BDF2 is second order and L-stable for \(\gamma=2-\sqrt{2}\). [Bank et al., 1996] ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/0168927495001158))

I will not expand the full two-stage algebra here (it is longer), but it can be implemented with **two tridiagonal solves per time step**, often still cheap in 1D. Its key advantage is stronger damping than CN without dropping to first order.

**When I recommend TR-BDF2 over RS-CN:**
- extremely low volatility \(\sigma\) where CN damping is too weak,
- when you see lingering high-frequency “ringing” in Greeks after payoff smoothing,
- or when you want a single consistent second-order method with stronger stability properties.

---

# 3.4 Treatment of non-smooth data (payoff smoothing + compatibility)

## 3.4.1 Cell-averaged payoff (recommended default)

Given grid nodes \(x_j\) and cell edges \(x_{j\pm 1/2}=x_j\pm h/2\), define the initial data as the **cell average**:
\[
u_j^0 = \frac{1}{h}\int_{x_{j-1/2}}^{x_{j+1/2}} \Phi(Ke^x)\,dx,
\quad j=1,\dots,J-1.
\]
This is Pooley–Forsyth–Vetzal’s “averaging the initial data” remedy.

**Why it works:** it removes the worst grid-scale modes associated with sampling a discontinuity/kink at a single point, improving Greeks.

## 3.4.2 Corner compatibility enforcement

Boundary values \(u(x_{\min},\tau)\), \(u(x_{\max},\tau)\) imply values at \(\tau=0\):
\[
u(x_{\min},0)=g_L(0),\qquad u(x_{\max},0)=g_R(0).
\]
The payoff \(\Phi\) must satisfy compatibility:
\[
\Phi(S_{\min}) = g_L(0),\qquad \Phi(S_{\max}) = g_R(0),
\]
otherwise you create immediate corner singularities.

**Practical fix:** replace payoff values at boundary nodes by boundary values at \(\tau=0\), and when cell averaging touches the boundary cells, average a payoff that has been “clamped” to match the boundary values.

---

# 3.5 Greeks computation (Delta, Gamma) with error orders

Compute Greeks at valuation time \(\tau=T\).

## 3.5.1 Delta and Gamma from log-space derivatives

Recall \(S=Ke^x\), \(V(S,\tau)=u(x,\tau)\).

Then:
\[
\Delta(S,\tau)=\frac{\partial V}{\partial S} = \frac{1}{S}\frac{\partial u}{\partial x},
\]
\[
\Gamma(S,\tau)=\frac{\partial^2 V}{\partial S^2}
= \frac{1}{S^2}\left(\frac{\partial^2 u}{\partial x^2}-\frac{\partial u}{\partial x}\right).
\]

On a uniform \(x\)-grid, second-order centered approximations are:
\[
u_x(x_j) \approx \frac{u_{j+1}-u_{j-1}}{2h} \quad (O(h^2)),
\]
\[
u_{xx}(x_j) \approx \frac{u_{j+1}-2u_j+u_{j-1}}{h^2}\quad (O(h^2)).
\]
So, at node \(j\),
\[
\Delta_j \approx \frac{1}{S_j}\frac{u_{j+1}-u_{j-1}}{2h},
\qquad
\Gamma_j \approx \frac{1}{S_j^2}\left(\frac{u_{j+1}-2u_j+u_{j-1}}{h^2}-\frac{u_{j+1}-u_{j-1}}{2h}\right),
\]
both second-order in space provided \(u\) is smooth at that point. Near the strike/barrier, payoff smoothing + damping are required for these formulas to behave well.

## 3.5.2 Interpolation to off-grid \(S_0\)

If \(x_0=\ln(S_0/K)\) lies between nodes, use **quadratic interpolation in \(x\)** over \((x_{j-1},x_j,x_{j+1})\) to evaluate \(u(x_0)\), \(u_x(x_0)\), \(u_{xx}(x_0)\). This gives \(O(h^3)\) local interpolation for \(u\) and \(O(h^2)\) for derivatives.

In production, a monotone cubic interpolation is often used for price surfaces, but for Gamma you must be careful: monotone cubic can reduce derivative accuracy. For risk, I recommend:
- quadratic interpolation for \(\Delta,\Gamma\) (stable, second-order),
- and monotone cubic only for price if needed.

---

# 3.6 Linear algebra

## 3.6.1 Tridiagonal solves

Each implicit Euler or CN step yields a tridiagonal system:
\[
A U^{n+1} = \text{rhs},
\]
with \(A\in\mathbb{R}^{(J-1)\times(J-1)}\) tridiagonal.

Use the **Thomas algorithm** (tridiagonal LU without pivoting), \(O(J)\) time and \(O(J)\) memory.

## 3.6.2 American options (extension hook)

For American options you solve an LCP:
\[
U^{n+1} \ge \Phi,\quad
A U^{n+1} \ge \text{rhs},\quad
(U^{n+1}-\Phi)\odot(AU^{n+1}-\text{rhs})=0.
\]

Best-practice solvers:
- **policy iteration** (simple, robust) [Reisinger–Witte, 2012] ([maths.ox.ac.uk](https://www.maths.ox.ac.uk/node/21883))
- penalty methods (also robust),
- splitting (Ikonen–Toivanen line).

I provide policy iteration pseudocode hooks in §3.7.10.

---

# 3.7 Complete algorithm — publication-quality pseudocode (self-contained)

I give one solver for:
- European options,
- with optional discrete monitoring projection (barriers),
- RS-CN time stepping,
- fitted spatial operator,
- cell-averaged payoff,
- tridiagonal Thomas solver,
- delta/gamma extraction.

### 3.7.1 Pseudocode conventions

- Arrays are 0-indexed in \(j\) (space), \(n\) (time).
- Interior indices are \(j=1,\dots,J-1\).
- Complexity statements refer to \(J\) (space points) and \(N\) (time steps).

---

```text
Algorithm 1: One-factor PDE pricer (log-space, fitted FD, Rannacher-smoothed CN)

Input:
  1  S0 > 0                 (spot)
  2  K  > 0                 (strike)
  3  T  > 0                 (maturity in years)
  4  r >= 0, q >= 0         (rates)
  5  sigma(S, tau)          (volatility function; constant allowed)
  6  Payoff Phi(S)          (function)
  7  OptionType ∈ {Call, Put, Custom}
  8  Domain parameters: x_min, x_max  (or build from m, sigma_ref)
  9  Grid sizes: J (space intervals), N (time steps)
 10  Monitoring times list TauMon[1..F] subset of (0,T] (possibly empty)
 11  Corridor [L,U] for projection (optional; if none, set L=0, U=+∞)
 12  Boundary functions:
       gL(tau) = left Dirichlet boundary value at x_min
       gR(tau) = right Dirichlet boundary value at x_max
 13  Damping steps m_damp = 2   (Rannacher: two half-steps)

Output:
 14  Price V(S0, today), Delta, Gamma
 15  Full grid u[j][n] if requested

Procedure:

// ─────────────────────────────────────────────────────────────────────
// A. Build grids
 16  h ← (x_max - x_min) / J
 17  k ← T / N
 18  for j = 0..J:
 19      x[j] ← x_min + j*h
 20      S[j] ← K * exp(x[j])
 21  for n = 0..N:
 22      tau[n] ← n*k

// ─────────────────────────────────────────────────────────────────────
// B. Initialize solution u at tau=0 using cell-averaged payoff
 23  Initialize u[j] for j=0..J
 24  u[0] ← gL(0)
 25  u[J] ← gR(0)
 26  for j = 1..J-1:
 27      // cell average in x: integrate Phi(K*exp(x)) over [x[j]-h/2, x[j]+h/2]
 28      aL ← x[j] - 0.5*h
 29      aR ← x[j] + 0.5*h
 30      // 3-point Gauss-Legendre on [aL,aR]
 31      // nodes in [-1,1]: ξ = { -sqrt(3/5), 0, +sqrt(3/5) }, weights w={5/9,8/9,5/9}
 32      xi1 ← -sqrt(3/5); xi2 ← 0; xi3 ← +sqrt(3/5)
 33      w1 ← 5/9; w2 ← 8/9; w3 ← 5/9
 34      mid ← 0.5*(aL+aR); half ← 0.5*(aR-aL)
 35      x1 ← mid + half*xi1
 36      x2 ← mid + half*xi2
 37      x3 ← mid + half*xi3
 38      f1 ← Phi( K*exp(x1) )
 39      f2 ← Phi( K*exp(x2) )
 40      f3 ← Phi( K*exp(x3) )
 41      u[j] ← (half*(w1*f1 + w2*f2 + w3*f3)) / h
 42  endfor

// ─────────────────────────────────────────────────────────────────────
// C. Time marching with monitoring projection
 43  nextMonIndex ← 1
 44  for n = 0..N-1:

 45      tau_n   ← tau[n]
 46      tau_np1 ← tau[n+1]

 47      // C1. Enforce boundary at current time level (for completeness)
 48      u[0] ← gL(tau_n)
 49      u[J] ← gR(tau_n)

 50      // C2. Determine whether we just applied a discontinuity projection
 51      // We apply projection when tau_np1 hits a monitoring time
 52      isMonitorStep ← false
 53      if nextMonIndex <= F and abs(tau_np1 - TauMon[nextMonIndex]) < 0.5*k:
 54          isMonitorStep ← true
 55      endif

 56      // C3. Choose step type:
 57      // Use implicit Euler half-steps for first m_damp steps after:
 58      //   (i) initial payoff at n=0
 59      //   (ii) a monitoring projection
 60      // Practical implementation: for n=0 and for steps immediately after a projection,
 61      // do two half steps with step size k/2.
 62      // We implement as:
 63      if n == 0:
 64          doDampingNow ← true
 65      else if previousStepWasProjection == true:
 66          doDampingNow ← true
 67      else:
 68          doDampingNow ← false
 69      endif

 70      if doDampingNow == true:

 71          // Two implicit Euler half-steps: (kE = k/2)
 72          kE ← 0.5*k
 73          repeat damp = 1..2:

 74              tau_stage ← tau_n + damp*kE   // tau_n + k/2 then tau_n + k

 75              // Assemble tridiagonal A for (I - kE L_h) u_new = u_old
 76              // Unknowns are interior j=1..J-1; boundaries go to RHS

 77              // Allocate arrays for Thomas: lower a[1..J-1], diag b[1..J-1], upper c[1..J-1], rhs d[1..J-1]
 78              for j = 1..J-1:
 79                  // coefficients at node j
 80                  sig ← sigma(S[j], tau_stage)
 81                  aj  ← 0.5 * sig*sig
 82                  bj  ← (r - q) - 0.5*sig*sig
 83                  // fitted factor rho
 84                  if aj <= 0:
 85                      // degenerate diffusion; fall back to upwind-like behavior by setting aj small
 86                      aj ← 1e-16
 87                  endif
 88                  theta ← (bj*h) / (2*aj)
 89                  if abs(theta) < 1e-8:
 90                      rho ← 1.0
 91                  else:
 92                      rho ← theta * coth(theta)
 93                  endif

 94                  // spatial stencil coefficients for L_h
 95                  ell ← (aj*rho)/(h*h) - (bj)/(2*h)
 96                  upp ← (aj*rho)/(h*h) + (bj)/(2*h)
 97                  diagL ← -2*(aj*rho)/(h*h) - r

 98                  // Build A = I - kE * L_h
 99                  a[j] ← -kE * ell        // lower diag of A (multiplies u_{j-1})
100                  b[j] ←  1 - kE * diagL  // main diag of A
101                  c[j] ← -kE * upp        // upper diag of A (multiplies u_{j+1})

102                  // RHS starts as u_old[j]
103                  d[j] ← u[j]

104              endfor

105              // Incorporate boundary values into RHS:
106              // equation for j=1 includes a[1]*u0 on LHS; move to RHS
107              d[1]     ← d[1]     - a[1]*gL(tau_stage)
108              // equation for j=J-1 includes c[J-1]*uJ on LHS; move to RHS
109              d[J-1]   ← d[J-1]   - c[J-1]*gR(tau_stage)

110              // Solve tridiagonal system A * u_new = d using Thomas algorithm (O(J))
111              // Forward elimination
112              for j = 2..J-1:
113                  m ← a[j] / b[j-1]
114                  b[j] ← b[j] - m * c[j-1]
115                  d[j] ← d[j] - m * d[j-1]
116              endfor
117              // Back substitution
118              u_new[J-1] ← d[J-1] / b[J-1]
119              for j = J-2..1 step -1:
120                  u_new[j] ← (d[j] - c[j]*u_new[j+1]) / b[j]
121              endfor

122              // Update u with u_new and boundaries
123              for j = 1..J-1:
124                  u[j] ← u_new[j]
125              endfor
126              u[0] ← gL(tau_stage)
127              u[J] ← gR(tau_stage)

128          endrepeat

129      else:
130          // Standard Crank–Nicolson step size k:
131          // (I - 0.5 k L_{n+1}) u^{n+1} = (I + 0.5 k L_n) u^n

132          // Assemble L at tau_n and tau_np1 (for coefficients).
133          // We build A (LHS) and RHS vector in one pass.

134          for j = 1..J-1:

135              // coefficients at tau_n
136              sig0 ← sigma(S[j], tau_n)
137              a0 ← 0.5*sig0*sig0
138              b0 ← (r-q) - 0.5*sig0*sig0
139              if a0 <= 0: a0 ← 1e-16 endif
140              th0 ← (b0*h)/(2*a0)
141              if abs(th0) < 1e-8: rho0 ← 1 else rho0 ← th0*coth(th0) endif
142              ell0 ← (a0*rho0)/(h*h) - b0/(2*h)
143              upp0 ← (a0*rho0)/(h*h) + b0/(2*h)
144              diag0 ← -2*(a0*rho0)/(h*h) - r

145              // coefficients at tau_np1
146              sig1 ← sigma(S[j], tau_np1)
147              a1 ← 0.5*sig1*sig1
148              b1 ← (r-q) - 0.5*sig1*sig1
149              if a1 <= 0: a1 ← 1e-16 endif
150              th1 ← (b1*h)/(2*a1)
151              if abs(th1) < 1e-8: rho1 ← 1 else rho1 ← th1*coth(th1) endif
152              ell1 ← (a1*rho1)/(h*h) - b1/(2*h)
153              upp1 ← (a1*rho1)/(h*h) + b1/(2*h)
154              diag1 ← -2*(a1*rho1)/(h*h) - r

155              // Build A = I - 0.5 k L(tau_np1)
156              a[j] ← -0.5*k*ell1
157              b[j] ←  1 - 0.5*k*diag1
158              c[j] ← -0.5*k*upp1

159              // Build RHS = (I + 0.5 k L(tau_n)) u^n
160              rhs ← u[j] + 0.5*k*( ell0*u[j-1] + diag0*u[j] + upp0*u[j+1] )
161              d[j] ← rhs

162          endfor

163          // Move boundary contributions for LHS to RHS
164          d[1]   ← d[1]   - a[1]*gL(tau_np1)
165          d[J-1] ← d[J-1] - c[J-1]*gR(tau_np1)

166          // Thomas solve
167          for j = 2..J-1:
168              m ← a[j]/b[j-1]
169              b[j] ← b[j] - m*c[j-1]
170              d[j] ← d[j] - m*d[j-1]
171          endfor
172          u_new[J-1] ← d[J-1]/b[J-1]
173          for j = J-2..1 step -1:
174              u_new[j] ← (d[j] - c[j]*u_new[j+1]) / b[j]
175          endfor

176          // Update u
177          for j = 1..J-1:
178              u[j] ← u_new[j]
179          endfor
180          u[0] ← gL(tau_np1)
181          u[J] ← gR(tau_np1)

182      endif  // end step type

183      // C4. Apply monitoring projection if required at tau_np1
184      previousStepWasProjection ← false
185      if isMonitorStep == true:
186          for j = 0..J:
187              if S[j] < L or S[j] > U:
188                  u[j] ← 0.0
189              endif
190          endfor
191          previousStepWasProjection ← true
192          nextMonIndex ← nextMonIndex + 1
193      endif

194  endfor  // end time loop

// ─────────────────────────────────────────────────────────────────────
// D. Extract price and Greeks at S0 (tau=T)

// D1. Interpolate in x to get u(x0), ux(x0), uxx(x0)
195  x0 ← ln(S0/K)
196  // find j such that x[j] <= x0 <= x[j+1]
197  j0 ← floor((x0 - x_min)/h)
198  if j0 < 1: j0 ← 1 endif
199  if j0 > J-2: j0 ← J-2 endif

200  // quadratic interpolation using nodes (j0-1, j0, j0+1)
201  // compute Lagrange basis and derivatives at x0
202  xA ← x[j0-1]; xB ← x[j0]; xC ← x[j0+1]
203  uA ← u[j0-1]; uB ← u[j0]; uC ← u[j0+1]

204  // Lagrange basis polynomials:
205  LA(x) = (x-xB)(x-xC)/((xA-xB)(xA-xC)), etc.
206  // Evaluate u(x0), u'(x0), u''(x0) analytically
207  denomA ← (xA-xB)*(xA-xC)
208  denomB ← (xB-xA)*(xB-xC)
209  denomC ← (xC-xA)*(xC-xB)

210  LA  ← (x0-xB)*(x0-xC)/denomA
211  LB  ← (x0-xA)*(x0-xC)/denomB
212  LC  ← (x0-xA)*(x0-xB)/denomC

213  // First derivatives:
214  dLA ← ((x0-xC)+(x0-xB))/denomA
215  dLB ← ((x0-xC)+(x0-xA))/denomB
216  dLC ← ((x0-xB)+(x0-xA))/denomC

217  // Second derivatives:
218  ddLA ← 2/denomA
219  ddLB ← 2/denomB
220  ddLC ← 2/denomC

221  u0   ← LA*uA + LB*uB + LC*uC
222  ux0  ← dLA*uA + dLB*uB + dLC*uC
223  uxx0 ← ddLA*uA + ddLB*uB + ddLC*uC

// D2. Map to price and Greeks in S
224  Price ← u0
225  Delta ← ux0 / S0
226  Gamma ← (uxx0 - ux0) / (S0*S0)

// Return (Price, Delta, Gamma)
227  return Price, Delta, Gamma

Complexity:
  - Each time step: O(J) assembly + O(J) Thomas solve → O(J)
  - Total: O(NJ)
  - Memory: O(J)
```

### 3.7.2 Notes on the pseudocode (important engineering clarifications)

1. **Monitoring time detection** (line 53) uses a tolerance. In production you should build the time grid so that monitoring times are hit exactly (avoid floating comparisons).
2. The pseudocode uses **uniform time steps**; for monitoring schedules, use **piecewise uniform** time grids per interval, with damping restarts in each.
3. The fitted factor uses \(\coth\); implement as:
   \[
   \coth(z)=\frac{\cosh z}{\sinh z},
   \]
   and for \(|z|\) large, use stable approximations to avoid overflow.
4. There is a runtime safeguard for \(a_j\le 0\) (degenerate diffusion). For strict models \(a>0\) always; this is mostly a numerical guard for extreme parameterization or local vol surfaces.

---

# 3.8 Error analysis summary table

Below is a practical summary (for one-factor, linear PDE; “monotone” means preserves a discrete maximum principle / positivity under stated conditions).

| Component | Method | Spatial order | Temporal order | Stability | Monotone? |
|---|---|---:|---:|---|---|
| Baseline | CN + centered space | 2 | 2 | A-stable but weak damping; oscillations for nonsmooth/convection-dominated | No (generally) |
| Baseline | Implicit Euler + centered space | 2 | 1 | Very stable | Only if spatial matrix is M-matrix (often false in \(\sigma^2\ll r\) with centered convection) |
| Spatial fix | Upwind for convection | 1 (or 2 with higher-order reconstructions) | depends | Stable | Yes (often) but adds \(O(h)\) diffusion |
| Spatial fix | Exponentially fitted FD (Il’in/Duffy) | 1–2 (regime dependent) | depends | Stable | Yes (M-matrix under conditions) |
| Spatial fix | Fitted finite volume (Song Wang) | provably convergent; often 2-ish | depends | Stable | Yes (M-matrix) ([academic.oup.com](https://academic.oup.com/imajna/article/24/4/699/687386)) |
| Time fix | Rannacher-smoothed CN | as chosen | 2 (after smoothing) | Robust | Improves monotonicity behavior; not unconditionally positive |
| Time fix | TR-BDF2 (L-stable) | as chosen | 2 | Strong damping | Positivity depends on step restriction/SSP variants ([arxiv.org](https://arxiv.org/abs/1510.04303)) |
| Data fix | Cell-averaged payoff | n/a | n/a | Reduces excitation of high modes | Helps (not a scheme by itself) ([risk.net](https://www.risk.net/journal-of-computational-finance/2160518/convergence-remedies-for-non-smooth-payoffs-in-option-pricing)) |
| Discrete monitoring | Restart damping after projection | n/a | n/a | Prevents reintroduced oscillations | Yes (as strategy) ([pure.kfupm.edu.sa](https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/)) |
| American constraint | Policy iteration | n/a | n/a | Robust | Preserves LCP structure ([maths.ox.ac.uk](https://www.maths.ox.ac.uk/node/21883)) |
| American constraint | Penalty | n/a | n/a | Robust | Approximate monotonicity; depends on penalty size |

---

# 3.9 Practical recommendations (production “just works”)

## 3.9.1 Default configuration (robust, accurate, cheap)

1. **Coordinates**: log-price \(x=\ln(S/K)\).
2. **Grid**: uniform \(x\)-grid with \(J\in[400,2000]\) depending on accuracy.
3. **Domain**: \(m=8\) standard deviations around \(\ln(S_0/K)\) (or barrier-aligned).
4. **Payoff**: cell-averaged payoff.
5. **Time stepping**: Rannacher smoothing (2 Euler half-steps) + CN.
6. **Spatial operator**: exponentially fitted FD (fitting factor \(\rho\)).
7. **Discrete monitoring**: apply projection at monitoring times; restart 2 Euler half-steps immediately after each projection.
8. **Greeks**: compute from log-derivatives using quadratic interpolation; avoid raw second differences in \(S\)-space.

This addresses every seed-paper failure mode:

- F1 (low vol): exponential fitting prevents central oscillations.
- F2 (nonsmooth payoff): cell averaging + Rannacher.
- F3 (monitoring): restart smoothing.
- F4 (positivity): M-matrix spatial operator + implicit solves; residual negativity typically negligible.
- F5 (numerical diffusion): controlled by mesh refinement near strike/barriers.
- F6 (truncation): log-space + wide domain; optionally use bounds or dynamic boundary methods.
- F8 (Greeks): compute in log-space, after smoothing.

## 3.9.2 When to switch methods

- If **\(\sigma^2 \ll r\)** and you still see oscillations under RS-CN:
  - increase damping steps (e.g., 4 Euler half-steps),
  - or switch CN phase to TR-BDF2,
  - or use Milev–Tagliani modified implicit scheme as a smoother in the first few steps. (This is precisely what [Milev–Tagliani, 2013] advocate.)
- If **far-field boundary error dominates** (long maturities, extreme strikes):
  - expand \([x_{\min},x_{\max}]\),
  - or consider boundary-error bounds,
  - or explore no-far-field methods like [Wang et al., 2025]. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S1572308925001068))
- If **American option**:
  - use implicit time stepping + policy iteration (or penalty) per step. ([maths.ox.ac.uk](https://www.maths.ox.ac.uk/node/21883))

## 3.9.3 Cost comparison (1D)

For tridiagonal solves:

- Implicit Euler: 1 solve per step ⇒ \(O(NJ)\).
- CN: 1 solve per step ⇒ \(O(NJ)\).
- RS-CN: ~1 solve per step + extra 2 half-steps once per discontinuity ⇒ still \(O(NJ)\).
- TR-BDF2: 2 solves per step ⇒ \(O(2NJ)\) but still cheap in 1D.

## 3.9.4 Pitfalls & debugging checklist

1. **Oscillations near strike/barrier**:
   - verify payoff smoothing is on,
   - verify damping steps are applied after every discontinuity (including monitoring resets),
   - check local \(\ell_j,u_j^{(+)}\ge 0\); if violated, refine or increase fitting/upwinding.

2. **Negative prices**:
   - ensure boundary conditions are consistent and compatible at \(\tau=0\),
   - ensure matrix is M-matrix (off-diagonals nonpositive in \(I-kL\)),
   - ensure RHS of CN step doesn’t create negativity (reduce \(k\) or increase damping).

3. **Greeks are noisy**:
   - compute Greeks in log-space,
   - use quadratic interpolation rather than raw node differences,
   - increase damping steps.

4. **Convergence testing**:
   - test price convergence under \(h\to h/2\), \(k\to k/2\),
   - test delta/gamma convergence separately (they converge slower without proper smoothing).

5. **Boundary sensitivity**:
   - expand domain and confirm interior price/Greeks are insensitive,
   - for barriers, align barrier with node.

---

# Consolidated reference list (Phase 1–3)

*(I list the most load-bearing references used above; many more exist, but this is already a strong “monograph chapter” spine.)*

- **[Duffy, 2004]** Daniel J. Duffy. *A Critique of the Crank Nicolson Scheme: Strengths and Weaknesses for Financial Instrument Pricing.* Wilmott Magazine, 2004.
- **[Milev & Tagliani, 2010b]** Mariyan Milev, Aldo Tagliani. *Low Volatility Options and Numerical Diffusion of Finite Difference Schemes.* Serdica Mathematical Journal 36 (2010), 223–236.
- **[Milev & Tagliani, 2010a]** Mariyan Milev, Aldo Tagliani. *Nonstandard Finite Difference Schemes with Application to Finance: Option Pricing.* Serdica Mathematical Journal, 2010, pp. 75–88. (As summarized in your `non_standard_FDM.md`.)
- **[Pooley, Forsyth & Vetzal, 2003]** David M. Pooley, Peter A. Forsyth, Kenneth R. Vetzal. *Convergence remedies for non-smooth payoffs in option pricing.* Journal of Computational Finance (circa 2003). ([risk.net](https://www.risk.net/journal-of-computational-finance/2160518/convergence-remedies-for-non-smooth-payoffs-in-option-pricing))  
- **[Giles & Carter, 2006]** Michael B. Giles, Rebecca Carter. *Convergence analysis of Crank–Nicolson and Rannacher time-marching.* Journal of Computational Finance 9(4), 2006. DOI: 10.21314/JCF.2006.152. ([risk.net](https://www.risk.net/journal-of-computational-finance/2160349/convergence-analysis-of-crank-nicolson-and-rannacher-time-marching))  
- **[Wade et al., 2007]** B. A. Wade, A. Q. M. Khaliq, M. Yousuf, J. Vigo-Aguiar, R. Deininger. *On smoothing of the Crank–Nicolson scheme and higher order schemes for pricing barrier options.* Journal of Computational and Applied Mathematics 204(1), 144–158, 2007. DOI: 10.1016/j.cam.2006.04.034. ([pure.kfupm.edu.sa](https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/))  
- **[Song Wang, 2004]** Song Wang. *A novel fitted finite volume method for the Black–Scholes equation governing option pricing.* IMA Journal of Numerical Analysis 24(4), 699–720, 2004. DOI: 10.1093/imanum/24.4.699. ([academic.oup.com](https://academic.oup.com/imajna/article/24/4/699/687386))  
- **[Valkov, 2012]** Radoslav Valkov. *Fitted Finite Volume Method for a Generalized Black-Scholes Equation Transformed on Finite Interval.* arXiv:1211.1903, 2012. ([arxiv.org](https://arxiv.org/abs/1211.1903))  
- **[Milev & Tagliani, 2013]** Mariyan Milev, Aldo Tagliani. *Efficient implicit scheme with positivity preserving and smoothing properties.* Journal of Computational and Applied Mathematics 243 (2013), 1–9. DOI: 10.1016/j.cam.2012.09.039.
- **[Tagliani & Milev, 2013]** Aldo Tagliani, Mariyan Milev. *Laplace Transform and finite difference methods for the Black–Scholes equation.* Applied Mathematics and Computation 220 (2013), 649–658. DOI: 10.1016/j.amc.2013.07.011.
- **[Goll, Rannacher & Wollner, 2015]** Christian Goll, Rolf Rannacher, Winnifried Wollner. *The damped Crank–Nicolson time-marching scheme for the adaptive solution of the Black–Scholes equation.* Journal of Computational Finance 18(4), 2015. DOI: 10.21314/JCF.2015.301. ([risk.net](https://www.risk.net/journal-of-computational-finance/2406534/the-damped-crank-nicolson-time-marching-scheme-for-the-adaptive-solution-of-the-black-scholes-equation))  
- **[Umeorah & Mashele, 2019]** Nneka Umeorah, Phillip Mashele. *A Crank-Nicolson finite difference approach on the numerical estimation of rebate barrier option prices.* Cogent Economics & Finance, 2019. DOI: 10.1080/23322039.2019.1598835.
- **[Cont & Voltchkova, 2005]** Rama Cont, Ekaterina Voltchkova. *A finite difference scheme for option pricing in jump diffusion and exponential Lévy models.* SIAM Journal on Numerical Analysis 43(4), 1596–1626, 2005. DOI: 10.1137/S0036142903436186. ([ideas.repec.org](https://ideas.repec.org/p/hal/journl/halshs-00445645.html))  
- **[Reisinger & Whitley, 2012]** C. Reisinger, A. Whitley. *The impact of a natural time change on the convergence of the Crank-Nicolson scheme.* arXiv:1210.5487, 2012. ([arxiv.org](https://arxiv.org/abs/1210.5487))  
- **[Reisinger & Witte, 2012]** Christoph Reisinger, J. Witte. *On the Use of Policy Iteration as an Easy Way of Pricing American Options.* SIAM Journal on Financial Mathematics, 2012. ([maths.ox.ac.uk](https://www.maths.ox.ac.uk/node/21883))  
- **[Bank et al., 1996]** Bank et al. *Analysis and implementation of TR-BDF2.* Applied Numerical Mathematics 20(1–2), 21–37, 1996. DOI: 10.1016/0168-9274(95)00115-8. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/0168927495001158))  
- **[Bonaventura & Della Rocca, 2015]** Luca Bonaventura, Alessandro Della Rocca. *Monotonicity, positivity and strong stability of the TR-BDF2 method and of its SSP extensions.* arXiv:1510.04303, 2015. ([arxiv.org](https://arxiv.org/abs/1510.04303))  
- **[López‑Salas et al., 2024]** J. G. López‑Salas et al. *IMEX-RK finite volume methods for nonlinear 1d parabolic PDEs. Application to option pricing.* arXiv:2409.01125, 2024. ([arxiv.org](https://arxiv.org/abs/2409.01125))  
- **[Wang et al., 2025]** Jian Wang et al. *A second-order finite difference method for the Black–Scholes model without far-field boundary conditions.* Journal of Financial Stability 81 (Dec 2025), 101477. DOI: 10.1016/j.jfs.2025.101477. ([sciencedirect.com](https://www.sciencedirect.com/science/article/abs/pii/S1572308925001068))  

---

# Quality-gate self-evaluation (as requested)

1. **Could a PhD student implement the Phase 3 algorithm from the pseudocode alone?**  
   Yes: grid, payoff smoothing, fitted coefficients, RS-CN logic, Thomas solve, interpolation, Greeks mapping are all explicitly specified.

2. **Are stability/convergence claims backed by proof sketches or citations?**  
   Yes: Rannacher/CN convergence and nonsmooth payoffs are tied to [Pooley–Forsyth–Vetzal, 2003] and [Giles–Carter, 2006]; monotonicity/positivity tied to M-matrix arguments and fitted FV references [Song Wang, 2004] and [Valkov, 2012]; TR-BDF2 properties cited to [Bank et al., 1996] and [Bonaventura–Della Rocca, 2015].

3. **Is notation consistent?**  
   Yes: \(S, V, \sigma, r, K, T, h, k, j, n\). I separated diffusion coefficient \(a\) from volatility \(\sigma\).

4. **Did I address every CN failure mode identified in the seed papers?**  
   Yes: convection dominance, nonsmooth payoff, boundary conditions, monitoring resets, von Neumann vs monotonicity, Greeks oscillations, numerical diffusion trade-off.

5. **Are Greeks computation methods specified with error orders?**  
   Yes: log-derivative formulas + centered differences \(O(h^2)\) away from nonsmooth points; quadratic interpolation for off-grid values.

</response>

---
Learn more:
1. [https://www.risk.net/journal-of-computational-finance/2160349/convergence-analysis-of-crank-nicolson-and-rannacher-time-marching](https://www.risk.net/journal-of-computational-finance/2160349/convergence-analysis-of-crank-nicolson-and-rannacher-time-marching)
2. [https://www.risk.net/journal-of-computational-finance/2160518/convergence-remedies-for-non-smooth-payoffs-in-option-pricing](https://www.risk.net/journal-of-computational-finance/2160518/convergence-remedies-for-non-smooth-payoffs-in-option-pricing)
3. [https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/](https://pure.kfupm.edu.sa/en/publications/on-smoothing-of-the-crank-nicolson-scheme-and-higher-order-scheme/)
4. [https://academic.oup.com/imajna/article/24/4/699/687386](https://academic.oup.com/imajna/article/24/4/699/687386)
5. [https://arxiv.org/abs/1211.1903](https://arxiv.org/abs/1211.1903)
6. [https://www.tse-fr.eu/articles/finite-difference-scheme-option-pricing-jump-diffusion-and-exponential-levy-models](https://www.tse-fr.eu/articles/finite-difference-scheme-option-pricing-jump-diffusion-and-exponential-levy-models)
7. [https://arxiv.org/abs/1210.5487](https://arxiv.org/abs/1210.5487)
8. [https://www.risk.net/journal-of-computational-finance/2160418/adaptive-and-high-order-methods-for-valuing-american-options](https://www.risk.net/journal-of-computational-finance/2160418/adaptive-and-high-order-methods-for-valuing-american-options)
9. [https://www.maths.ox.ac.uk/node/21883](https://www.maths.ox.ac.uk/node/21883)
10. [https://arxiv.org/abs/2409.01125](https://arxiv.org/abs/2409.01125)
11. [https://www.sciencedirect.com/science/article/abs/pii/S1572308925001068](https://www.sciencedirect.com/science/article/abs/pii/S1572308925001068)
12. [https://global-sci.com/index.php/ijnam/article/view/9955](https://global-sci.com/index.php/ijnam/article/view/9955)
13. [https://www.sciencedirect.com/science/article/pii/S0377042716304678](https://www.sciencedirect.com/science/article/pii/S0377042716304678)
14. [https://arxiv.org/abs/1510.04303](https://arxiv.org/abs/1510.04303)
15. [https://www.sciencedirect.com/science/article/abs/pii/0168927495001158](https://www.sciencedirect.com/science/article/abs/pii/0168927495001158)
16. [https://www.risk.net/journal-of-computational-finance/2406534/the-damped-crank-nicolson-time-marching-scheme-for-the-adaptive-solution-of-the-black-scholes-equation](https://www.risk.net/journal-of-computational-finance/2406534/the-damped-crank-nicolson-time-marching-scheme-for-the-adaptive-solution-of-the-black-scholes-equation)
17. [https://www.sciencedirect.com/science/article/abs/pii/S0378475410000522](https://www.sciencedirect.com/science/article/abs/pii/S0378475410000522)
18. [https://www.sciencedirect.com/science/article/abs/pii/S0378475420304596](https://www.sciencedirect.com/science/article/abs/pii/S0378475420304596)
19. [https://arxiv.org/abs/2103.05369](https://arxiv.org/abs/2103.05369)
20. [https://ideas.repec.org/p/hal/journl/halshs-00445645.html](https://ideas.repec.org/p/hal/journl/halshs-00445645.html)
21. [https://www.sciencedirect.com/science/article/abs/pii/S0096300307012179](https://www.sciencedirect.com/science/article/abs/pii/S0096300307012179)
22. [https://arxiv.org/abs/2405.08280](https://arxiv.org/abs/2405.08280)
23. [https://academic.oup.com/imajna/article-abstract/23/2/241/684490](https://academic.oup.com/imajna/article-abstract/23/2/241/684490)

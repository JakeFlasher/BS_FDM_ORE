

<planning>
The 3–5 most critical logical chains that must hold for the method to work:

1. **Exp-fitted core chain**: E(1) → E(2)/U4 → E(3) (fitting factor) → U6 (fitted operator) → corrected U7 (tridiagonal matrix A) → U8 (sign conditions) → C3 (M-matrix) → C4 (A⁻¹ > 0) → E(5) (norm bound) → C6 (positivity) + C7 (max principle). This is the foundational chain for Scheme 1.

2. **CN variant core chain**: E(1) → E(8) (six-node reaction term) → U17–U19 (P, N matrices) → G1.5 (completing-the-square for ω) → E(9) (parameter choice + Δt constraint) → C16 (P is M-matrix, N ≥ 0) → C17 (positivity) + C20 (max principle). This is the foundational chain for Scheme 2.

3. **Low-σ degradation chain**: E(3) → U12 (σ→0 limit of ρ) → corrected upwind schemes (denominator h, not 2h) → G1.2 (Taylor expansion) → U15/C12/C14 (numerical diffusion ½rSΔS) → E(7) (modified equation). This justifies the paper's central concern about artificial diffusion.

4. **CN variant diffusion chain**: E(8) + E(9) → G1.4 (Taylor expansion of modified reaction term) → C22 (artificial diffusion ⅛(rΔS/σ)²) → E(10) (modified equation). This quantifies the CN variant's diffusion.

5. **Practical validity chain**: E(9) constraint evaluated at examples → G7.3 (most examples violate E(9)) → explanation of why the scheme still works → guidance on when to enforce E(9) vs. accept relaxed guarantees.

Verification strategy: trace each chain step-by-step, performing explicit algebra at every transition. Flag any step where the transition depends on an external result or an unverified assumption.
</planning>

---

# Round 1.4 — End-to-End Verification

---

## 1. CORRECTED EQUATION CHAIN

This chain is **authoritative**: it supersedes the paper where Round 1.3 corrections were made. A reader should be able to implement both schemes from this chain alone.

### Notation conventions

Throughout: \(h \equiv \Delta S\) (spatial step), \(k \equiv \Delta t\) (time step), \(\sigma\) denotes volatility (scalar), and \(\sigma_d(S) = \tfrac{1}{2}\sigma^2 S^2\) denotes the diffusion coefficient when disambiguation is needed. The variable \(t\) is **time-to-expiry** (backward): \(t=0\) is maturity, and the scheme marches forward from \(t=0\) to \(t=T\) (the present).

---

### CE-1. Black–Scholes PDE

$$-\frac{\partial V}{\partial t} + rS\frac{\partial V}{\partial S} + \frac{1}{2}\sigma^2 S^2\frac{\partial^2 V}{\partial S^2} - rV = 0, \qquad S \in (0, S_{\max}),\; t \in (0, T]$$

**Source:** E(1). Unchanged. Standard BS PDE under the risk-neutral measure with constant \(r > 0\), \(\sigma > 0\).

---

### CE-2. Generalized form and coefficient identification

$$-\frac{\partial V}{\partial t} + \mu(S)\frac{\partial V}{\partial S} + \sigma_d(S)\frac{\partial^2 V}{\partial S^2} + b\,V = 0$$

with

$$\mu(S) = rS, \qquad \sigma_d(S) = \tfrac{1}{2}\sigma^2 S^2, \qquad b = -r$$

**Source:** E(2) and U4. Unchanged. The subscript \(d\) disambiguates the diffusion coefficient from the volatility parameter.

---

### CE-3. Computational grid

**Spatial** (uniform):

$$S_j = j\,\Delta S, \qquad j = 0, 1, \ldots, M, \qquad S_{\max} = M\,\Delta S$$

**Temporal** (uniform):

$$t_n = n\,\Delta t, \qquad n = 0, 1, \ldots, N_t, \qquad T = N_t\,\Delta t$$

**Source:** G2.2/G2.3/G6.3 resolution. Not explicitly in the paper; inferred from tridiagonal structure and figure captions.

**Requirements:**

- \(S_{\max} \gg U\) (upper barrier) so that \(V(S_{\max},t) \approx 0\). Rule of thumb: \(S_{\max} \geq 2U\).
- The interior unknowns are at nodes \(j = 1, 2, \ldots, M-1\). Nodes \(j=0\) and \(j=M\) are boundary nodes with prescribed values.

---

### CE-4. Boundary conditions

For the truncated call (Definition 4.1) and discrete double barrier knock-out call (Definition 4.2):

$$U_0^n = 0, \qquad U_M^n = 0, \qquad \text{for all } n$$

**Source:** G2.4/G6.1 resolution. Follows from E(12) and the option contracts (payoff is zero at \(S=0\) and well above the upper barrier). Homogeneous Dirichlet conditions leave the interior tridiagonal system unmodified.

For options with non-zero boundary values (e.g., digital puts with \(V(0,t) = Ae^{-rt}\)): modify the right-hand side at \(j=1\) and \(j=M-1\) by subtracting the coupling to the known boundary value.

---

### CE-5. Payoff initialization (\(t = 0\), i.e., maturity)

**Truncated call:**

$$U_j^0 = \begin{cases} S_j - K & \text{if } S_j \in [K, U] \\ 0 & \text{otherwise} \end{cases}$$

**Discrete double barrier knock-out call:**

$$U_j^0 = \max(S_j - K,\, 0)\;\mathbf{1}_{[L,U]}(S_j)$$

**Source:** Definitions 4.1, 4.2 and E(11). Unchanged.

---

### CE-6. Monitoring update (discrete barriers only)

At each monitoring date \(t_i\) (\(i = 1, \ldots, F\), with \(0 = t_0 < t_1 < \cdots < t_F = T\)), after computing the solution at \(t = t_i\):

$$U_j \;\leftarrow\; U_j \cdot \mathbf{1}_{[L,U]}(S_j), \qquad j = 0, 1, \ldots, M$$

**Source:** E(13). Unchanged. This is a **post-processing** step between time solves that re-introduces discontinuities at the barrier levels.

---

## Scheme 1: Exponentially Fitted Implicit Scheme

### CE-7. Fitting factor

$$\rho_j = \frac{\mu_j\,\Delta S}{2}\,\coth\!\left(\frac{\mu_j\,\Delta S}{2\,\sigma_d^{(j)}}\right)$$

where \(\mu_j = rS_j\) and \(\sigma_d^{(j)} = \tfrac{1}{2}\sigma^2 S_j^2\).

**Péclet number:** Define \(\text{Pe}_j = \frac{\mu_j\,\Delta S}{2\,\sigma_d^{(j)}} = \frac{r\,\Delta S}{\sigma^2 S_j} = \frac{r}{\sigma^2 j}\).

Then \(\rho_j = \sigma_d^{(j)} \cdot \text{Pe}_j\,\coth(\text{Pe}_j)\).

**Source:** E(3). Unchanged (notation adapted to use \(\Delta S\) and constant-coefficient identifications).

**Numerical guard:** For \(|\text{Pe}_j| < \epsilon\) (e.g., \(\epsilon = 10^{-6}\)): use \(\text{Pe}\,\coth(\text{Pe}) \approx 1 + \text{Pe}^2/3\), so \(\rho_j \approx \sigma_d^{(j)}\). For \(|\text{Pe}_j| > 20\): \(\coth(\text{Pe}_j) \approx \text{sign}(\text{Pe}_j)\), so \(\rho_j \approx |\mu_j|\Delta S/2\). Otherwise: compute directly.

---

### CE-8. Fitted operator (set to zero)

$$-\frac{U_j^{n+1} - U_j^n}{\Delta t} + \mu_j\frac{U_{j+1}^{n+1} - U_{j-1}^{n+1}}{2\Delta S} + \rho_j\frac{U_{j+1}^{n+1} - 2U_j^{n+1} + U_{j-1}^{n+1}}{\Delta S^2} + b\,U_j^{n+1} = 0$$

**Source:** U6 with \(L_k^h U_j^n = 0\). All spatial terms evaluated at level \(n{+}1\) (fully implicit).

**Justification:** This is a backward-Euler discretization in time, centered in space, with \(\rho_j\) replacing the true diffusion coefficient \(\sigma_d^{(j)}\).

---

### CE-9. Tridiagonal system

$$A\,U^{n+1} = U^n$$

where \(U^n = (U_1^n, U_2^n, \ldots, U_{M-1}^n)^T\) (interior nodes only, since \(U_0 = U_M = 0\)), and

$$\boxed{A = \operatorname{tridiag}\!\left\{\underbrace{\left(-\frac{\rho_j}{\Delta S^2} + \frac{\mu_j}{2\Delta S}\right)\!\Delta t}_{\displaystyle a_j^-};\;\;\underbrace{\left(\frac{2\rho_j}{\Delta S^2} - b + \frac{1}{\Delta t}\right)\!\Delta t}_{\displaystyle a_j^0};\;\;\underbrace{-\left(\frac{\rho_j}{\Delta S^2} + \frac{\mu_j}{2\Delta S}\right)\!\Delta t}_{\displaystyle a_j^+}\right\}}$$

**⚠ CORRECTION (G4.1):** The paper's U7 writes superscript \(n\) on \(\rho_j^n, \mu_j^n, b_j^n\). The fitted operator CE-8 evaluates all coefficients at level \(n{+}1\). For the constant-coefficient BS case (\(r, \sigma\) constant), both are identical. For variable coefficients, **use level \(n{+}1\)**.

**Original (U7):** \(\rho_j^n, \mu_j^n, b_j^n\) in matrix entries.
**Corrected:** \(\rho_j^{n+1}, \mu_j^{n+1}, b_j^{n+1}\) (or simply \(\rho_j, \mu_j, b\) for constant coefficients).

**Derivation from CE-8:** Multiply CE-8 by \(\Delta t\), move all \(U^{n+1}\) terms to the left and \(U^n\) to the right, negate to make the diagonal positive:

Coefficient of \(U_{j-1}^{n+1}\): convection contributes \(-\mu_j\Delta t/(2\Delta S)\), diffusion contributes \(+\rho_j\Delta t/\Delta S^2\). Total: \(a_j^- = (-\rho_j/\Delta S^2 + \mu_j/(2\Delta S))\Delta t\).

Coefficient of \(U_j^{n+1}\): time contributes \(+1\), diffusion contributes \(+2\rho_j\Delta t/\Delta S^2\), reaction contributes \(-b\Delta t = r\Delta t\). Total: \(a_j^0 = (2\rho_j/\Delta S^2 - b + 1/\Delta t)\Delta t = 1 + r\Delta t + 2\rho_j\Delta t/\Delta S^2\).

Coefficient of \(U_{j+1}^{n+1}\): \(a_j^+ = -(\rho_j/\Delta S^2 + \mu_j/(2\Delta S))\Delta t\).

**Row sum:** \(a_j^- + a_j^0 + a_j^+ = (1/\Delta t - b)\Delta t = 1 - b\Delta t = 1 + r\Delta t\), independent of \(j\).

---

### CE-10. Sign conditions and M-matrix structure

For all interior nodes \(j = 1, \ldots, M-1\) with \(\mu_j = rS_j > 0\) and \(\sigma_d^{(j)} > 0\):

**Super-diagonal:** \(a_j^+ = -(\rho_j/\Delta S^2 + \mu_j/(2\Delta S))\Delta t < 0\). ✓ (Both terms in parenthesis are positive.)

**Sub-diagonal:** \(a_j^- = (-\rho_j/\Delta S^2 + \mu_j/(2\Delta S))\Delta t\). Since \(\coth(x) > 1\) for \(x > 0\), we have \(\rho_j > \mu_j\Delta S/2\), hence \(\rho_j/\Delta S^2 > \mu_j/(2\Delta S)\), so \(a_j^- < 0\). ✓

**Diagonal:** \(a_j^0 = 1 + r\Delta t + 2\rho_j\Delta t/\Delta S^2 > 0\). ✓

**Off-diagonals nonzero** (since \(\rho_j > 0\) and \(\mu_j > 0\)): \(A\) is irreducible.

**Strict diagonal dominance:** Row sum \(= 1 + r\Delta t > 0\), and diagonal \(= \text{row sum} + |a_j^-| + |a_j^+| > |a_j^-| + |a_j^+|\).

**Conclusion** (Ortega [5], §§6.2.3, 6.2.17): \(A\) is a nonsingular, irreducibly diagonally dominant M-matrix, so \(A^{-1} > 0\) (all entries strictly positive).

---

### CE-11. Norm bound (Windisch [9])

$$\|A^{-1}\|_\infty \;\leq\; \frac{1}{1 + r\Delta t} \;<\; 1$$

**Derivation:** For an M-matrix with constant diagonal dominance excess \(d = 1 + r\Delta t\) (the row sum), Windisch's theorem gives \(\|A^{-1}\|_\infty \leq 1/d\).

---

### CE-12. Positivity and maximum principle

**Positivity:** For \(U^0 \geq 0\) with at least one positive component: \(U^n = (A^{-1})^n U^0 > 0\) for all \(n \geq 1\).

**⚠ CORRECTION (G4.5):** The paper writes \((A^{-1})^n U^0 > 0\) for \(U^0 \geq 0\). Strict positivity requires \(U^0 \neq 0\). For \(U^0 = 0\): \(U^n = 0\).

**Maximum principle:**

$$\|U^{n+1}\|_\infty = \|A^{-1}U^n\|_\infty \;\leq\; \|A^{-1}\|_\infty\,\|U^n\|_\infty \;\leq\; \|U^n\|_\infty$$

**⚠ CORRECTION (G4.2):** The paper writes "=" for the second relation. The correct relation is "\(\leq\)" by submultiplicativity. The conclusion is unaffected.

---

### CE-13. Convergence bound

$$|V(S_j, t_n) - U_j^n| \;\leq\; c\,(h + k)$$

where \(c\) is independent of \(h, k, \sigma\).

**Source:** E(6). Stated without proof; justified by the general theory for exponentially fitted schemes (Duffy [2], Ch. 19; Miller–O'Riordan–Shishkin framework). The \(O(h)\) spatial order (rather than \(O(h^2)\)) is the worst case as \(\sigma \to 0\); for fixed \(\sigma > 0\), the scheme achieves \(O(h^2 + k)\).

**Status:** Partially verified (proof strategy established in G3.3/G7.1; complete proof requires external reference).

---

### CE-14. Low-volatility limit (\(\sigma \to 0\))

As \(\sigma \to 0\), \(\text{Pe}_j \to \infty\), and \(\coth(\text{Pe}_j) \to 1\), giving:

$$\rho_j \;\to\; \frac{|\mu_j|\,\Delta S}{2} = \frac{rS_j\,\Delta S}{2}$$

Substituting into CE-8 (and combining convection with artificial diffusion terms, as derived in G4.3):

**⚠ CORRECTION (G4.3):** The upwind scheme has denominator \(\Delta S\), not \(2\Delta S\) as in the paper's U13/U14.

For \(\mu_j > 0\) (always true for BS with \(r > 0\)):

$$\boxed{-\frac{U_j^{n+1} - U_j^n}{\Delta t} + \mu_j\frac{U_{j+1}^{n+1} - U_j^{n+1}}{\Delta S} + b\,U_j^{n+1} = 0}$$

**Original (U13):** denominator \(2h\). **Corrected:** denominator \(h = \Delta S\).

**Derivation (G4.3):** Substitute \(\rho = \mu\Delta S/2\) into the centered convection + fitted diffusion terms of CE-8:

\[\frac{\mu}{2\Delta S}(U_{j+1} - U_{j-1}) + \frac{\mu\Delta S/2}{\Delta S^2}(U_{j+1} - 2U_j + U_{j-1})\]
\[= \frac{\mu}{2\Delta S}(U_{j+1} - U_{j-1}) + \frac{\mu}{2\Delta S}(U_{j+1} - 2U_j + U_{j-1})\]
\[= \frac{\mu}{2\Delta S}[2U_{j+1} - 2U_j] = \frac{\mu}{\Delta S}(U_{j+1} - U_j)\]

---

### CE-15. Numerical diffusion of the exp-fitted scheme

The forward difference \((U_{j+1} - U_j)/\Delta S\) has truncation error \(\frac{\Delta S}{2}V'' + O(\Delta S^2)\). Therefore the upwind scheme (and by extension the exp-fitted scheme for large Pe) effectively solves:

$$-\frac{\partial V}{\partial t} + \mu(S)\frac{\partial V}{\partial S} + \frac{1}{2}\mu(S)\,\Delta S\,\frac{\partial^2 V}{\partial S^2} + b\,V = 0$$

Substituting \(\mu = rS\):

$$\text{Artificial diffusion (exp-fitted):}\quad \frac{1}{2}rS\,\Delta S\,\frac{\partial^2 V}{\partial S^2}$$

Valid when \(\text{Pe}_j \gg 1\), i.e., \(\sigma^2 S_j \ll r\Delta S\).

**Source:** U15, C12, C14, E(7). Confirmed by Taylor expansion (G1.2).

---

## Scheme 2: Crank–Nicolson Variant

### CE-16. Discretization structure

$$P\,U^{n+1} = N\,U^n$$

where \(P\) and \(N\) are \((M{-}1)\times(M{-}1)\) tridiagonal matrices. Standard Crank–Nicolson averaging is used for convection and diffusion; the reaction term \(-rV\) is discretized via the six-node approximation E(8).

---

### CE-17. Matrix \(P\) (left-hand side)

$$P = \operatorname{tridiag}\!\left\{p_j^-;\; p_j^0;\; p_j^+\right\}$$

with

$$p_j^- = r\omega_2 + \frac{rS_j}{4\Delta S} - \left(\frac{\sigma S_j}{2\Delta S}\right)^{\!2}$$

$$p_j^0 = \frac{1}{\Delta t} + \frac{1}{2}\!\left(\frac{\sigma S_j}{\Delta S}\right)^{\!2} + r\!\left(\frac{1}{2} - 2\omega_2\right)$$

$$p_j^+ = r\omega_2 - \frac{rS_j}{4\Delta S} - \left(\frac{\sigma S_j}{2\Delta S}\right)^{\!2}$$

**Source:** U18. Unchanged. **Derivation:** G1.3 (full algebra verified, all 6 coefficient components independently checked).

---

### CE-18. Matrix \(N\) (right-hand side)

$$N = \operatorname{tridiag}\!\left\{n_j^-;\; n_j^0;\; n_j^+\right\}$$

with

$$n_j^- = -r\omega_1 - \frac{rS_j}{4\Delta S} + \left(\frac{\sigma S_j}{2\Delta S}\right)^{\!2}$$

$$n_j^0 = \frac{1}{\Delta t} - \frac{1}{2}\!\left(\frac{\sigma S_j}{\Delta S}\right)^{\!2} - r\!\left(\frac{1}{2} - 2\omega_1\right)$$

$$n_j^+ = -r\omega_1 + \frac{rS_j}{4\Delta S} + \left(\frac{\sigma S_j}{2\Delta S}\right)^{\!2}$$

**Source:** U19. Unchanged. **Derivation:** G1.3.

**Row sums** (verified: all \(S_j\)-dependent terms cancel pairwise):

$$\text{Row sum of } P = \frac{1}{\Delta t} + \frac{r}{2}, \qquad \text{Row sum of } N = \frac{1}{\Delta t} - \frac{r}{2}$$

Both are constant across rows (independent of \(j\)).

---

### CE-19. Parameter choice

$$\boxed{\omega_1 = \omega_2 = -\frac{r}{16\sigma^2}}$$

**Derivation (G1.5):** Complete the square in \(S_j\) for the sub-diagonal of \(P\):

$$p_j^- = -\frac{\sigma^2}{4\Delta S^2}\!\left(S_j - \frac{r\Delta S}{2\sigma^2}\right)^{\!2} + \frac{r^2}{16\sigma^2} + r\omega_2$$

The maximum over all \(S_j \geq 0\) is \(\frac{r^2}{16\sigma^2} + r\omega_2\). Requiring \(p_j^- \leq 0\) for all \(j\) gives \(\omega_2 \leq -\frac{r}{16\sigma^2}\). The tightest choice (minimizing artificial diffusion) is equality. By symmetry of the argument applied to \(N\)'s sub-diagonal (requiring \(\geq 0\)): \(\omega_1 = -\frac{r}{16\sigma^2}\).

**With this choice, the off-diagonals simplify to perfect squares:**

$$p_j^- = -\left(\frac{\sigma S_j}{2\Delta S} - \frac{r}{4\sigma}\right)^{\!2} \leq 0, \qquad n_j^- = +\left(\frac{\sigma S_j}{2\Delta S} - \frac{r}{4\sigma}\right)^{\!2} \geq 0$$

---

### CE-20. Time-step constraint

$$\boxed{\Delta t \;<\; \frac{1}{r\!\left(\frac{1}{2} + \frac{r}{8\sigma^2}\right) + \frac{1}{2}(\sigma M)^2}}$$

**Derivation:** Requires the diagonal of \(N\) to be non-negative at the worst-case node \(j = M{-}1\):

$$n_j^0 = \frac{1}{\Delta t} - \frac{\sigma^2 S_j^2}{2\Delta S^2} - r\!\left(\frac{1}{2} + \frac{r}{8\sigma^2}\right) \;\geq\; 0$$

Evaluating at \(S_{M-1} \approx M\Delta S\) and solving for \(\Delta t\) gives CE-20.

**Source:** E(9). Unchanged.

**Practical note (G7.3):** For the paper's low-volatility examples, CE-20 is extremely restrictive (e.g., \(\Delta t < 3.2 \times 10^{-5}\) for \(r=0.5\), \(\sigma=0.001\), \(M=2800\)). The paper's examples use \(\Delta t = 0.01\), violating CE-20. The scheme still produces non-oscillatory results because: (a) the diagonal of \(N\) becomes negative only near \(S_{\max}\), where the solution is essentially zero; (b) \(P\) remains an M-matrix regardless of \(\Delta t\); (c) CE-20 is a sufficient condition, not a sharp stability boundary.

---

### CE-21. M-matrix and positivity (under CE-19 and CE-20)

**\(P\) is an M-matrix:** Off-diagonals \(\leq 0\) (CE-19); diagonal \(> 0\) (unconditional); row sum \(= 1/\Delta t + r/2 > 0\) gives strict diagonal dominance; off-diagonals nonzero for almost all \(j\) (see G7.4 for the edge case). Hence \(P^{-1} \geq 0\).

**\(N \geq 0\):** Sub-diagonal is a perfect square \(\geq 0\); super-diagonal \(> 0\) for \(S_j > 0\); diagonal \(\geq 0\) under CE-20.

**Positivity:** \(U^{n+1} = P^{-1}NU^n = (P^{-1}N)^n U^0 \geq 0\) for \(U^0 \geq 0\).

---

### CE-22. Maximum principle (under CE-19 and CE-20)

$$\|U^{n+1}\|_\infty \;\leq\; \|P^{-1}\|_\infty\,\|N\|_\infty\,\|U^n\|_\infty \;\leq\; \frac{1/\Delta t - r/2}{1/\Delta t + r/2}\,\|U^n\|_\infty \;\leq\; \|U^n\|_\infty$$

using \(\|N\|_\infty = 1/\Delta t - r/2\) (row sum, since \(N \geq 0\)) and \(\|P^{-1}\|_\infty \leq (1/\Delta t + r/2)^{-1}\) (Windisch).

**⚠ CORRECTION (G4.4):** The paper's U24 writes "=" for the first relation. Correct: "\(\leq\)" by submultiplicativity.

---

### CE-23. Numerical diffusion of CN variant

The six-node approximation E(8) introduces an additional spatial error from the Taylor expansion (G1.4):

$$\omega_1(U_{j-1}^n + U_{j+1}^n) + (\tfrac{1}{2}-2\omega_1)U_j^n = \tfrac{1}{2}V + \omega_1\Delta S^2\,V'' + O(\Delta S^4)$$

Summing the \(n\) and \(n{+}1\) contributions and multiplying by \(-r\):

$$\text{Artificial diffusion (CN variant):}\quad \frac{1}{8}\!\left(\frac{r\,\Delta S}{\sigma}\right)^{\!2}\frac{\partial^2 V}{\partial S^2}$$

When this dominates the physical diffusion \(\frac{1}{2}\sigma^2 S^2\,\partial^2 V/\partial S^2\) (i.e., when \(\sigma^2 S \ll r\Delta S/2\)), the scheme effectively solves:

$$-\frac{\partial V}{\partial t} + rS\frac{\partial V}{\partial S} + \frac{1}{8}\!\left(\frac{r\,\Delta S}{\sigma}\right)^{\!2}\frac{\partial^2 V}{\partial S^2} - rV = 0$$

**Source:** C22, U26, E(10). Confirmed by G1.4 derivation.

**Discretization error:** Formally \(O(\Delta S^2 + \Delta t^2)\), but the \(\Delta S^2\) coefficient is \(r^2/(8\sigma^2)\), which is enormous for small \(\sigma\) (e.g., 31,250 for \(r=0.5\), \(\sigma=0.001\)). Accurate results require \(\Delta S\) small enough that \(\frac{1}{8}(r\Delta S/\sigma)^2\) is negligible.

---

### CE-24. Time-stepping algorithm (both schemes)

**Exp-fitted scheme:** At each time step \(n \to n{+}1\):
1. Set RHS \(= U^n\) (vector of interior values).
2. Solve \(AU^{n+1} = \text{RHS}\) by the Thomas algorithm.

**CN variant:** At each time step \(n \to n{+}1\):
1. Compute RHS \(= NU^n\) (tridiagonal matrix–vector product).
2. Solve \(PU^{n+1} = \text{RHS}\) by the Thomas algorithm.

For **constant coefficients**, both \(A\) and \((P, N)\) are assembled once and reused.

At **monitoring dates** (discrete barriers): apply CE-6 after the solve.

---

### CE-25. Comparison of artificial diffusion

| Scheme | Artificial diffusion coefficient | Depends on | Order in \(\Delta S\) |
|:---|:---|:---|:---|
| Exp-fitted (low \(\sigma\)) | \(\frac{1}{2}rS\,\Delta S\) | \(r\), \(S\), \(\Delta S\) | \(O(\Delta S)\) |
| CN variant | \(\frac{1}{8}(r\Delta S/\sigma)^2\) | \(r\), \(\sigma\), \(\Delta S\) | \(O(\Delta S^2)\) |

The CN variant's diffusion decreases faster under grid refinement (\(\Delta S^2\) vs. \(\Delta S\)), but has a larger coefficient for small \(\sigma\). For a given \(\Delta S\), the crossover in magnitude occurs at \(S = r\Delta S/(4\sigma^2)\), which for the paper's parameters lies far outside the computational domain.

---

## 2. DIMENSIONAL CONSISTENCY

Units: \(V\) in \(\$\); \(S\) in \(\$\); \(t\) in years; \(r\) in yr\(^{-1}\); \(\sigma\) (volatility) in yr\(^{-1/2}\); \(\sigma^2\) in yr\(^{-1}\).

| Eq ID | LHS dimension | RHS dimension | Status |
|:---|:---|:---|:---|
| CE-1 / E(1) | Each term: \(\$/\)yr | \(\$/\)yr | ✅ Consistent |
| CE-2 / E(2) | Each term: \(\$/\)yr | \(\$/\)yr | ✅ (\([\mu]=\$/\text{yr}\), \([\sigma_d]=\$^2/\text{yr}\), \([b]=\text{yr}^{-1}\)) |
| CE-7 / E(3) | \([\rho]=\$^2/\text{yr}\) | \([\mu\Delta S/2]=\$^2/\text{yr}\); \(\coth(\cdot)\) dimensionless | ✅ |
| CE-8 / U6 | Each term: \(\$/\)yr | \(\$/\)yr | ✅ |
| CE-9 / U7 | \([A]\) dimensionless | all entries: \([\text{yr}^{-1}] \cdot [\text{yr}]=1\) | ✅ (\(AU^{n+1}=U^n\): both sides in \(\$\)) |
| CE-11 / E(5) | \(\|A^{-1}\|_\infty\): dimensionless | \(1/(1+r\Delta t)\): dimensionless | ✅ |
| CE-13 / E(6) | \(|V - U|: \$\) | \(c(h+k)\): conventional — \(c\) absorbs derivative bounds and units | ✅ (standard convention) |
| CE-15 / E(7) | Each term: \(\$/\)yr | \(\$/\)yr; \([\tfrac{1}{2}\mu\Delta S]=\$^2/\text{yr}\) | ✅ |
| CE-17 / U18 | \([P]\): yr\(^{-1}\) per entry | Sub: \([r\omega]+[rS/\Delta S]+[\sigma^2 S^2/\Delta S^2]\): all yr\(^{-1}\) | ✅ (\(PU=\$/\text{yr}\)) |
| CE-18 / U19 | \([N]\): yr\(^{-1}\) per entry | Same structure as \(P\) | ✅ |
| CE-19 / E(9) \(\omega\) | dimensionless | \([r/\sigma^2]=(\text{yr}^{-1})/(\text{yr}^{-1})=1\) | ✅ |
| CE-20 / E(9) \(\Delta t\) | years | \(1/[\text{yr}^{-1}]\) = years | ✅ |
| CE-23 / U26 | \(\$/\)yr (as PDE term) | \([(r\Delta S/\sigma)^2]=\$^2/\text{yr}\); \(\times[\partial^2 V/\partial S^2]=\$/\$^2\): total \(\$/\)yr | ✅ |
| CE-23 / E(10) | Each term: \(\$/\)yr | \(\$/\)yr | ✅ |
| U34 (\(\tau_d\)) | years | \([\Delta S^2/(\sigma S)^2]=\$^2/((\text{yr}^{-1})\$^2)=\text{yr}\) | ✅ |

**No dimensional inconsistencies found.** The convergence bound E(6) uses the standard convention where \(c\) absorbs dimension-carrying constants (solution derivative bounds).

---

## 3. LIMITING CASE VERIFICATION

### (a) \(\sigma \to 0\): degradation to upwind

**Fitting factor limit:**
\(\text{Pe}_j = r/(\sigma^2 j) \to \infty\). Then \(\coth(\text{Pe}_j) \to 1\), so

$$\rho_j \;\to\; \frac{rS_j\,\Delta S}{2}$$

**Matrix \(A\) in the limit:** Substituting \(\rho_j = rS_j\Delta S/2\) into CE-9:

$$a_j^- = \left(-\frac{rS_j\Delta S/2}{\Delta S^2} + \frac{rS_j}{2\Delta S}\right)\Delta t = \left(-\frac{rS_j}{2\Delta S} + \frac{rS_j}{2\Delta S}\right)\Delta t = 0$$

$$a_j^0 = 1 + r\Delta t + \frac{2 \cdot rS_j\Delta S/2}{\Delta S^2}\Delta t = 1 + r\Delta t + \frac{rS_j}{\Delta S}\Delta t$$

$$a_j^+ = -\left(\frac{rS_j}{2\Delta S} + \frac{rS_j}{2\Delta S}\right)\Delta t = -\frac{rS_j}{\Delta S}\Delta t$$

The system at row \(j\) becomes:

$$(1 + r\Delta t + rS_j\Delta t/\Delta S)\,U_j^{n+1} - (rS_j\Delta t/\Delta S)\,U_{j+1}^{n+1} = U_j^n$$

Rearranging:

$$-\frac{U_j^{n+1} - U_j^n}{\Delta t} + rS_j\frac{U_{j+1}^{n+1} - U_j^{n+1}}{\Delta S} - rU_j^{n+1} = 0$$

This is the **first-order implicit forward-difference upwind scheme** with denominator \(\Delta S\). ✅ Confirms CE-14 and the G4.3 correction.

**Verification of the CN variant as \(\sigma \to 0\):** \(\omega_1 = \omega_2 = -r/(16\sigma^2) \to -\infty\). The \(\Delta t\) constraint CE-20 becomes \(\Delta t < 8\sigma^2/r^2 \to 0\). So the CN variant is not directly usable in the pure \(\sigma = 0\) limit — it requires increasingly small \(\Delta t\). This is consistent with the paper's finding that the CN variant introduces artificial diffusion proportional to \((r/\sigma)^2\).

---

### (b) \(r \to 0\) (no dividend yield \(q\) in this paper)

The paper does not include a dividend yield. For \(r \to 0\): \(\mu_j = rS_j \to 0\), \(b = -r \to 0\).

**Exp-fitted scheme:** \(\text{Pe}_j = r/(\sigma^2 j) \to 0\). Using \(x\coth(x) \to 1\) as \(x \to 0\):

$$\rho_j \;\to\; \sigma_d^{(j)} = \tfrac{1}{2}\sigma^2 S_j^2$$

Matrix \(A\) becomes the standard fully implicit centered scheme for pure diffusion:

$$a_j^- = -\frac{\sigma^2 S_j^2}{2\Delta S^2}\Delta t, \qquad a_j^0 = 1 + \frac{\sigma^2 S_j^2}{\Delta S^2}\Delta t, \qquad a_j^+ = -\frac{\sigma^2 S_j^2}{2\Delta S^2}\Delta t$$

✅ The convection and reaction terms vanish; pure diffusion with coefficient \(\frac{1}{2}\sigma^2 S^2\) is recovered.

**CN variant:** \(\omega \to 0\). Matrices simplify to the standard Crank–Nicolson for \(-\partial_t V + \frac{1}{2}\sigma^2 S^2\partial_{SS}V = 0\):

$$p_j^- = n_j^- = \frac{\sigma^2 S_j^2}{4\Delta S^2}, \qquad p_j^0 = \frac{1}{\Delta t} + \frac{\sigma^2 S_j^2}{2\Delta S^2}, \qquad n_j^0 = \frac{1}{\Delta t} - \frac{\sigma^2 S_j^2}{2\Delta S^2}$$

Wait — actually the signs differ between \(P\) and \(N\). Let me recalculate:

With \(r = 0\), \(\omega = 0\):

\(p_j^- = 0 + 0 - (\sigma S_j/(2\Delta S))^2 = -\alpha_j\). \(n_j^- = 0 - 0 + (\sigma S_j/(2\Delta S))^2 = +\alpha_j\).

\(p_j^0 = 1/\Delta t + 2\alpha_j + 0 = 1/\Delta t + 2\alpha_j\). \(n_j^0 = 1/\Delta t - 2\alpha_j\).

\(p_j^+ = 0 - 0 - \alpha_j = -\alpha_j\). \(n_j^+ = 0 + 0 + \alpha_j = +\alpha_j\).

where \(\alpha_j = (\sigma S_j/(2\Delta S))^2\).

✅ This is the standard CN scheme for the diffusion equation \(-\partial_t V + \frac{1}{2}\sigma^2 S^2\partial_{SS}V = 0\).

---

### (c) \(\Delta S \to 0\) with \(\Delta t\) fixed: consistency

**Exp-fitted scheme:** As \(\Delta S \to 0\), \(\text{Pe}_j = r\Delta S/(\sigma^2 S_j) \to 0\). Using the expansion \(\text{Pe}\,\coth(\text{Pe}) = 1 + \text{Pe}^2/3 + O(\text{Pe}^4)\):

$$\rho_j = \sigma_d^{(j)}\!\left(1 + \frac{\text{Pe}_j^2}{3} + O(\text{Pe}_j^4)\right) = \sigma_d^{(j)} + \frac{\mu_j^2\Delta S^2}{12\sigma_d^{(j)}} + O(\Delta S^4)$$

Taylor-expanding the spatial differences in CE-8 about \(S_j\) (using \(V_{j\pm 1} = V \pm V'\Delta S + \frac{1}{2}V''\Delta S^2 \pm \frac{1}{6}V'''\Delta S^3 + \cdots\)):

$$\mu_j\frac{V_{j+1}-V_{j-1}}{2\Delta S} = \mu_j\!\left(V' + \frac{\Delta S^2}{6}V''' + O(\Delta S^4)\right)$$

$$\rho_j\frac{V_{j+1}-2V_j+V_{j-1}}{\Delta S^2} = \left(\sigma_d^{(j)} + \frac{\mu_j^2\Delta S^2}{12\sigma_d^{(j)}}\right)\!\left(V'' + \frac{\Delta S^2}{12}V'''' + O(\Delta S^4)\right)$$

The leading terms give \(\mu V' + \sigma_d V'' = \) the spatial part of the PDE. The residual is:

$$\text{Spatial truncation error} = O(\Delta S^2)$$

Combined with the backward-Euler time truncation error \(O(\Delta t)\):

$$\text{Total truncation error} = O(\Delta S^2 + \Delta t) \;\to\; O(\Delta t) \text{ as } \Delta S \to 0$$

✅ Consistency with CE-1 is recovered. (For the uniform-in-\(\sigma\) worst case, the effective spatial order degrades to \(O(\Delta S)\) as analyzed in CE-13.)

**CN variant:** The standard CN components give \(O(\Delta S^2 + \Delta t^2)\) truncation error. The modification from \(\omega\) adds \(O(\Delta S^2)\) with coefficient \(r^2/(8\sigma^2)\). As \(\Delta S \to 0\):

$$\text{Total truncation error} \to O(\Delta t^2)$$

✅ Consistency recovered.

---

### (d) Constant coefficients: textbook form

For constant \(r, \sigma\), the matrix entries in CE-9 (exp-fitted) and CE-17/CE-18 (CN variant) vary across rows only through \(S_j = j\Delta S\) and \(\rho_j\). The parameters \(r, \sigma, \omega\) are constant, so the matrices are assembled once. The time-level correction (G4.1) is immaterial since \(\mu_j^n = \mu_j^{n+1} = rS_j\) for all \(n\).

The exp-fitted matrix CE-9 matches the formulation in Duffy [2], Ch. 19. The CN variant CE-17/CE-18 matches the formulation in Milev–Tagliani [4].

✅ Both reduce to their known textbook forms.

---

## 4. PROPERTY VERIFICATION

| # | Property | Corrected conditions | Grade | Notes |
|:---|:---|:---|:---|:---|
| 1 | **M-matrix: \(A\) (exp-fitted)** | \(r > 0\), \(\sigma > 0\), \(S_j > 0\) for interior nodes | **Fully verified** | Sign conditions from \(\coth(x)>1\); diagonal dominance from row sum \(1+r\Delta t > 0\); irreducibility from nonzero off-diagonals |
| 2 | **\(A^{-1} > 0\)** | M-matrix property (row 1) | **Fully verified** | Ortega [5] §6.2.17 |
| 3 | **Norm bound \(\|A^{-1}\|_\infty \leq 1/(1+r\Delta t)\)** | M-matrix with constant row sum | **Fully verified** | Windisch [9] |
| 4 | **Positivity of exp-fitted solution** | \(U^0 \geq 0\), \(U^0 \neq 0\) | **Fully verified** | Induction; correction G4.5 (strict positivity needs \(U^0 \neq 0\)) |
| 5 | **Max principle (exp-fitted)** | Unconditional (\(r > 0\)) | **Fully verified** | Correction G4.2 (\(\leq\) not \(=\) in norm chain); conclusion unaffected |
| 6 | **Convergence \(O(h+k)\) uniform in \(\sigma\)** | Fitted scheme with E(3) | **Verified with caveats** | Proof strategy in G3.3; complete proof requires Duffy [2] Ch. 19 or Miller–O'Riordan–Shishkin |
| 7 | **Oscillation-free (exp-fitted)** | M-matrix property | **Fully verified** | Consequence of \(A^{-1} > 0\) and max principle |
| 8 | **M-matrix: \(P\) (CN variant)** | CE-19 (\(\omega\) choice); critical \(S^* = r\Delta S/(2\sigma^2)\) outside grid | **Verified with caveats** | If \(S^*\) falls on the grid, \(P\) becomes reducible at that row. For all paper examples, \(S^* \gg S_{\max}\). Implementer safeguard: check \(S^* > S_{\max}\); if not, perturb \(\omega\) slightly |
| 9 | **\(N \geq 0\) (CN variant)** | CE-19 + CE-20 | **Fully verified** | Sub-diagonal is a perfect square \(\geq 0\); super \(> 0\); diagonal \(\geq 0\) under CE-20 |
| 10 | **Positivity (CN variant)** | CE-19 + CE-20 + \(U^0 \geq 0\) | **Fully verified** | \(P^{-1} \geq 0\) and \(N \geq 0\) implies \(P^{-1}N \geq 0\) |
| 11 | **Max principle (CN variant)** | CE-19 + CE-20 | **Fully verified** | Correction G4.4 (\(\leq\) not \(=\)); requires \(N \geq 0\) for \(\|N\|_\infty = \text{row sum}\) |
| 12 | **Accuracy \(O(\Delta S^2 + \Delta t^2)\) (CN variant)** | Standard CN + modification | **Verified with caveats** | Formal order correct; coefficient of \(\Delta S^2\) term is \(r^2/(8\sigma^2)\), which is large for small \(\sigma\). "Second order" may require impractically small \(\Delta S\). |
| 13 | **Artificial diffusion \(\frac{1}{2}rS\Delta S\) (exp-fitted)** | \(\text{Pe}_j \gg 1\), i.e., \(\sigma^2 S_j \ll r\Delta S\) | **Fully verified** | Asymptotic expansion of \(\coth\); Taylor expansion of upwind difference |
| 14 | **Artificial diffusion \(\frac{1}{8}(r\Delta S/\sigma)^2\) (CN variant)** | Low \(\sigma\), i.e., \(\sigma^2 S \ll r\Delta S/2\) | **Fully verified** | Taylor expansion of six-node reaction term (G1.4) |

---

## 5. IMPLEMENTATION READINESS CHECKLIST

| Component | Fully specified? | Remaining ambiguity / action required |
|:---|:---|:---|
| **Coordinate transform** | ✅ Yes (none used) | Paper works in original \((S,t)\) coordinates. No log-transform. |
| **Grid construction** | ✅ Yes (after G2.2/G6.3) | Uniform grids \(S_j = j\Delta S\), \(t_n = n\Delta t\). \(M = S_{\max}/\Delta S\), \(N_t = T/\Delta t\). Must choose \(S_{\max} \geq 2U\). |
| **Payoff initialization** | ✅ Yes | CE-5 specifies both truncated call and discrete barrier payoffs. |
| **Operator assembly (exp-fitted)** | ✅ Yes | CE-7 (fitting factor with numerical guards) → CE-9 (tridiagonal matrix). Assemble once for constant coefficients. Must handle \(\coth\) overflow for large Péclet. |
| **Operator assembly (CN variant)** | ✅ Yes | CE-19 (\(\omega\)) → CE-17/CE-18 (\(P, N\)). Assemble once for constant coefficients. Check CE-20 constraint; if violated, acknowledge loss of formal guarantees. |
| **Time stepping** | ✅ Yes | CE-24: Thomas algorithm for tridiagonal solve. Both \(A\) and \(P\) are M-matrices, so Thomas is stable without pivoting. |
| **Boundary conditions** | ✅ Yes (after G2.4/G6.1) | Homogeneous Dirichlet \(U_0 = U_M = 0\) for all paper examples. System is \((M{-}1)\times(M{-}1)\) for interior unknowns. For non-zero BCs, modify RHS at \(j=1\) and \(j=M{-}1\). |
| **Monitoring update** | ✅ Yes | CE-6: pointwise projection \(U_j \leftarrow U_j \cdot \mathbf{1}_{[L,U]}(S_j)\) at each monitoring date. Post-processing step; no system modification needed. |
| **Greeks extraction** | ❌ **Not specified** | Paper does not discuss Greeks. Standard approach: \(\Delta \approx (U_{j+1}-U_{j-1})/(2\Delta S)\), \(\Gamma \approx (U_{j+1}-2U_j+U_{j-1})/\Delta S^2\), \(\Theta\) from the PDE or \(-(U_j^{n+1}-U_j^n)/\Delta t\). Must be implemented by the user. |
| **Parameter constraints** | ✅ Yes (after G7.3) | **Exp-fitted:** No \(\Delta S, \Delta t\) restriction for stability (unconditional). For accuracy: \(\Delta S\) small enough that \(\frac{1}{2}rS_{\max}\Delta S\) is negligible vs. desired accuracy. **CN variant:** CE-19 for \(\omega\); CE-20 for \(\Delta t\) (sufficient condition). If CE-20 is relaxed, monitor for negative values. |
| **Variable coefficients** | ⚠️ Partially specified (G6.2) | For \(r(S,t)\), \(\sigma(S,t)\): exp-fitted scheme generalizes naturally (evaluate \(\mu_j^{n+1}, \sigma_d^{(j),n+1}\) at each node/time). CN variant requires node-dependent \(\omega_j = -r(S_j)/(16\sigma^2(S_j))\) or a conservative global choice \(\omega = -\max r/(16\min\sigma^2)\). Matrices must be reassembled at each time step. |

**Summary:** All components needed for a constant-coefficient implementation of both schemes are fully specified after the Round 1.3 resolutions. The only unspecified component is Greeks extraction (not discussed in the paper). Variable-coefficient extensions are partially specified with clear guidance.

# Round 1.3 — Gap Resolution

---

## Implementation Blockers

### GAP G2.2 + G2.3 + G6.3: Grid construction unstated

**RESOLUTION TYPE:** derived + sourced

**RESOLUTION:**

The paper never formally specifies the computational grid. From the tridiagonal matrix entries (U7, U18, U19, U33), which depend on \(S_j\), and from figure captions (e.g., \(S_{\max}=140\), \(\Delta S=0.05\)), the following construction is inferred and confirmed by the companion paper (Note-3, §3 "Discretization setup"):

A **uniform spatial grid** is used:

$$S_j = j\,\Delta S, \qquad j = 0, 1, \ldots, M, \qquad S_{\max} = M\,\Delta S$$

A **uniform temporal grid** is used:

$$t_n = n\,\Delta t, \qquad n = 0, 1, \ldots, X, \qquad T = X\,\Delta t$$

The number of spatial nodes is therefore:

$$M = \frac{S_{\max}}{\Delta S}$$

For Example 4.1: \(S_{\max} = 140\), \(\Delta S = 0.05\), so \(M = 2800\). For Example 4.2: \(S_{\max}\) is not stated in the text but can be inferred from the figure axis; a typical choice for a double barrier option with \(U=110\) would be \(S_{\max} = 200\) (matching Note-3's Example 4.1).

The choice of \(S_{\max}\) must satisfy \(S_{\max} \gg U\) (the upper barrier) so that the far-field boundary condition \(V(S_{\max},t) \approx 0\) is valid. A common rule of thumb (from Note-2, §5.1) is \(S_{\max} = R \cdot U\) with \(R \geq 2\).

**VERIFICATION:**

- Consistent with U7 entries at \(S_j = j\Delta S\), U18/U19 entries involving \(S_j/\Delta S = j\), and U33.
- Confirmed by Note-3 §3: identical grid specification for the same authors' companion paper.
- Example 4.1 figure caption gives \(S_{\max}=140\), \(\Delta S=0.05\), yielding integer \(M = 2800\). ✓

**STATUS:** Resolved.

---

### GAP G2.4 + G6.1: Boundary treatment unspecified

**RESOLUTION TYPE:** derived + sourced

**RESOLUTION:**

The tridiagonal systems (E(4), U17, U32) are stated for interior nodes \(j=1,\ldots,M-1\) only. Boundary conditions at \(j=0\) and \(j=M\) must be specified to close the system.

**Step 1: Identify the boundary conditions from the PDE problem.**

From E(12): \(V(S,t) \to 0\) as \(S \to 0\) or \(S \to \infty\).

For the truncated call (Definition 4.1) and the discrete barrier knock-out call (Definition 4.2), both contracts satisfy:

$$V(0, t) = 0, \qquad V(S_{\max}, t) = 0 \quad \text{for all } t \in [0,T]$$

These are **homogeneous Dirichlet conditions**. The first holds because a call option on an asset with \(S=0\) is worthless; the second holds because \(S_{\max}\) is chosen far enough above the upper barrier \(U\) that the truncated/knocked-out payoff is zero.

**Step 2: Incorporate into the tridiagonal system.**

For any of the three schemes (exponentially fitted, CN variant, standard implicit), the equation at interior node \(j\) couples nodes \(j-1\), \(j\), and \(j+1\). At the boundary nodes:

*At \(j=1\):* The sub-diagonal term involves \(U_0^{n+1}\). With \(U_0^{n+1} = 0\), this term vanishes. The first row of the system becomes:

$$a_{1,1}\,U_1^{n+1} + a_{1,2}\,U_2^{n+1} = (\text{RHS})_1 - a_{1,0} \cdot 0 = (\text{RHS})_1$$

*At \(j=M-1\):* The super-diagonal term involves \(U_M^{n+1}\). With \(U_M^{n+1} = 0\), the last row becomes:

$$a_{M-1,M-2}\,U_{M-2}^{n+1} + a_{M-1,M-1}\,U_{M-1}^{n+1} = (\text{RHS})_{M-1} - a_{M-1,M} \cdot 0 = (\text{RHS})_{M-1}$$

Since the boundary values are zero for all the paper's examples, the tridiagonal system requires **no modification** — the interior-node formulation as written is complete.

**Step 3: Barrier monitoring update (discrete monitoring).**

At each monitoring date \(t_i\), the projection E(13) is applied pointwise:

$$V_j^{n+1} \leftarrow V_j^{n+1} \cdot \mathbf{1}_{[L,U]}(S_j), \qquad j = 0, 1, \ldots, M$$

This sets \(V_j = 0\) for all nodes outside \([L, U]\), re-introducing discontinuities at \(S = L\) and \(S = U\). No modification to the tridiagonal structure is needed — this is a post-processing step after each solve.

**For non-zero boundary conditions** (e.g., digital put options as in Note-1, Example 4.1, where \(V(0,t) = Ae^{-rt}\)):

The boundary contribution modifies the right-hand side. At \(j=1\), the RHS becomes:

$$(\text{RHS})_1 \leftarrow (\text{RHS})_1 - a_{1,0} \cdot g_0^{n+1}$$

where \(g_0^{n+1} = V(0, t_{n+1})\) is the known boundary value. Similarly at \(j=M-1\) with \(g_M^{n+1}\).

**VERIFICATION:**

- Consistent with E(12) and the option contract definitions (Definitions 4.1, 4.2).
- Confirmed by Note-3 §3.4 and Note-1 §3.1, which use identical Dirichlet conditions.
- The homogeneous case \(g_0 = g_M = 0\) leaves the tridiagonal system unmodified, explaining why the paper never mentions boundary treatment explicitly.

**STATUS:** Resolved.

---

### GAP G4.1: Time-level indices in U7 — \(n\) vs. \(n{+}1\)

**RESOLUTION TYPE:** corrected

**RESOLUTION:**

The fitted operator U6 evaluates all spatial coefficients at time level \(n{+}1\):

$$L_k^h U_j^n \equiv -\frac{U_j^{n+1} - U_j^n}{k} + \mu_j^{n+1}\frac{U_{j+1}^{n+1} - U_{j-1}^{n+1}}{2h} + \rho_j^{n+1}\frac{\delta_x^2 U_j^{n+1}}{h^2} + b_j^{n+1} U_j^{n+1}$$

Setting \(L_k^h U_j^n = 0\) and rearranging (this also fills **G1.1**):

$$\frac{U_j^{n+1} - U_j^n}{k} = \mu_j^{n+1}\frac{U_{j+1}^{n+1} - U_{j-1}^{n+1}}{2h} + \rho_j^{n+1}\frac{U_{j+1}^{n+1} - 2U_j^{n+1} + U_{j-1}^{n+1}}{h^2} + b_j^{n+1} U_j^{n+1}$$

Multiply by \(k\) and move all \(U^{n+1}\) terms to the left:

**Coefficient of \(U_{j-1}^{n+1}\):**

From convection: \(-k\mu_j^{n+1} \cdot \frac{-1}{2h} = +\frac{k\mu_j^{n+1}}{2h}\)

From diffusion: \(-k\rho_j^{n+1} \cdot \frac{1}{h^2} = -\frac{k\rho_j^{n+1}}{h^2}\)

Total: \(\displaystyle\left(-\frac{\rho_j^{n+1}}{h^2} + \frac{\mu_j^{n+1}}{2h}\right)k\)

**Coefficient of \(U_j^{n+1}\):**

From time: \(+1\)

From diffusion: \(-k\rho_j^{n+1} \cdot \frac{-2}{h^2} = \frac{2k\rho_j^{n+1}}{h^2}\)

From reaction: \(-kb_j^{n+1}\)

Total: \(\displaystyle\left(\frac{2\rho_j^{n+1}}{h^2} - b_j^{n+1} + \frac{1}{k}\right)k\)

**Coefficient of \(U_{j+1}^{n+1}\):**

From convection: \(-\frac{k\mu_j^{n+1}}{2h}\)

From diffusion: \(-\frac{k\rho_j^{n+1}}{h^2}\)

Total: \(\displaystyle -\left(\frac{\rho_j^{n+1}}{h^2} + \frac{\mu_j^{n+1}}{2h}\right)k\)

The **corrected** tridiagonal matrix is:

$$A = \operatorname{tridiag}\!\left\{\left(-\frac{\rho_j^{n+1}}{h^2} + \frac{\mu_j^{n+1}}{2h}\right)k;\;\left(\frac{2\rho_j^{n+1}}{h^2} - b_j^{n+1} + \frac{1}{k}\right)k;\;-\left(\frac{\rho_j^{n+1}}{h^2} + \frac{\mu_j^{n+1}}{2h}\right)k\right\}$$

The only change from U7 is replacing superscript \(n\) with \(n{+}1\) on \(\rho_j\), \(\mu_j\), and \(b_j\). For the constant-coefficient Black–Scholes case (\(r, \sigma\) constant), \(\mu_j^n = \mu_j^{n+1} = rS_j\), so U7 is numerically correct despite the notational error.

**VERIFICATION:**

- Each coefficient is independently derived from U6 by explicit algebra. ✓
- Row sum: \(\left(-\frac{\rho}{h^2}+\frac{\mu}{2h}\right)k + \left(\frac{2\rho}{h^2}-b+\frac{1}{k}\right)k + \left(-\frac{\rho}{h^2}-\frac{\mu}{2h}\right)k = (1/k - b)k = 1 - kb\). With \(b = -r\): row sum \(= 1+kr\), matching E(5). ✓
- For constant coefficients: \(n+1 \equiv n\), so the corrected matrix coincides with U7. ✓

**STATUS:** Resolved. **Downstream impact:** For variable-coefficient implementations, the corrected superscripts must be used. All subsequent results (C2–C10, E(5), E(6)) are unaffected for constant coefficients.

---

### GAP G7.3: E(9) time-step constraint violated in examples

**RESOLUTION TYPE:** derived

**RESOLUTION:**

E(9) requires \(\Delta t < \frac{1}{r(\frac{1}{2} - 2\omega_1) + \frac{1}{2}(\sigma M)^2}\) with \(\omega_1 = -\frac{r}{16\sigma^2}\).

Substituting \(\frac{1}{2} - 2\omega_1 = \frac{1}{2} + \frac{r}{8\sigma^2}\):

$$\Delta t < \frac{1}{r\!\left(\frac{1}{2} + \frac{r}{8\sigma^2}\right) + \frac{1}{2}(\sigma M)^2}$$

**Evaluation for each example** (using \(\sigma = 0.001\), \(S_{\max}=140\), \(\Delta S = 0.05\), so \(M = 2800\)):

| Example | \(r\) | \(\frac{r}{8\sigma^2}\) | Denominator | \(\Delta t_{\max}\) | Used \(\Delta t\) | Satisfied? |
|:---|:---|:---|:---|:---|:---|:---|
| 4.1 (\(r=0.05\)) | 0.05 | 6250 | 316.4 | 0.00316 | 0.01 | **No** |
| 4.4 (\(r=0.01\)) | 0.01 | 1250 | 16.4 | 0.0609 | 0.01 | **Yes** |
| 4.4 (\(r=0.5\)) | 0.5 | 62500 | 31254 | \(3.2 \times 10^{-5}\) | 0.01 | **No** |
| 4.5 (\(r=0.5\)) | 0.5 | 62500 | 31254 | \(3.2 \times 10^{-5}\) | 0.01 | **No** |

Only the \(r=0.01\) case satisfies E(9).

**Why the scheme still works when E(9) is violated:**

The \(\Delta t\) constraint in E(9) comes from requiring the **diagonal** of \(N\) to be non-negative at the largest interior node \(j = M-1\):

$$\frac{1}{\Delta t} - \frac{1}{2}\!\left(\frac{\sigma S_{M-1}}{\Delta S}\right)^{\!2} - r\!\left(\frac{1}{2} - 2\omega_1\right) \geq 0$$

When violated, some diagonal entries of \(N\) become negative for nodes near \(S_{\max}\). However, two structural properties survive:

**(a) Row sums of \(N\) and \(P\) are independent of \(j\).** By direct computation (verified in G1.3's derivation below), all \(S_j\)-dependent terms cancel in the row sums:

$$\text{Row sum of } N = \frac{1}{\Delta t} - \frac{r}{2}, \qquad \text{Row sum of } P = \frac{1}{\Delta t} + \frac{r}{2}$$

Both are positive for any \(\Delta t < 2/r\) (always satisfied in practice).

**(b) The max-principle contraction** \(\|P^{-1}\|_\infty \|N\|_\infty \leq \frac{1/\Delta t - r/2}{1/\Delta t + r/2} < 1\) **holds only when \(N \geq 0\)** (entrywise), because \(\|N\|_\infty = \frac{1}{\Delta t} - \frac{r}{2}\) requires all entries to be non-negative. When some diagonal entries are negative, \(\|N\|_\infty > \frac{1}{\Delta t} - \frac{r}{2}\), and the simple bound breaks.

**(c) Practical mitigation:** In the paper's examples, the option value near \(S_{\max}\) is essentially zero (far from the payoff region), so the negative diagonal entries multiply negligible values and don't propagate significant errors into the interior.

**Implementer guidance:**

For guaranteed positivity and max-principle: enforce E(9). For practical use with moderate accuracy requirements: monitor \(\min_j V_j^n\) at each time step; if negative values appear, reduce \(\Delta t\) or use Rannacher-style smoothing (as in Note-1, §3.3).

**VERIFICATION:**

- Arithmetic for each example checked against E(9). ✓
- Row-sum independence of \(j\) verified by explicit algebra (see G1.3 resolution). ✓
- Practical explanation consistent with the paper's non-oscillatory results in Figs. 3–5. ✓

**STATUS:** Resolved. **Downstream impact:** C17 (positivity) and C20 (max principle) formally hold only when E(9) is satisfied. The paper's examples 4.1, 4.4 (\(r=0.5\)), and 4.5 lack these guarantees but work due to the zero-boundary effect near \(S_{\max}\).

---

### GAP G4.3: Denominators \(2h\) in U13/U14 should be \(h\)

**RESOLUTION TYPE:** corrected

**RESOLUTION:**

Substitute the limit \(\rho = \frac{\mu h}{2}\) (from U12, \(\mu > 0\)) into U6 with \(L_k^h U_j^n = 0\):

$$-\frac{U_j^{n+1} - U_j^n}{k} + \mu \frac{U_{j+1}^{n+1} - U_{j-1}^{n+1}}{2h} + \frac{\mu h}{2}\cdot\frac{U_{j+1}^{n+1} - 2U_j^{n+1} + U_{j-1}^{n+1}}{h^2} + b\,U_j^{n+1} = 0$$

**Step 1:** Simplify the artificial diffusion term:

$$\frac{\mu h}{2}\cdot\frac{U_{j+1}^{n+1} - 2U_j^{n+1} + U_{j-1}^{n+1}}{h^2} = \frac{\mu}{2h}\left(U_{j+1}^{n+1} - 2U_j^{n+1} + U_{j-1}^{n+1}\right)$$

**Step 2:** Combine with the convection term:

$$\frac{\mu}{2h}\left(U_{j+1}^{n+1} - U_{j-1}^{n+1}\right) + \frac{\mu}{2h}\left(U_{j+1}^{n+1} - 2U_j^{n+1} + U_{j-1}^{n+1}\right)$$

$$= \frac{\mu}{2h}\left[\left(U_{j+1} - U_{j-1}\right) + \left(U_{j+1} - 2U_j + U_{j-1}\right)\right]$$

$$= \frac{\mu}{2h}\left[2U_{j+1}^{n+1} - 2U_j^{n+1}\right] = \frac{\mu}{h}\left(U_{j+1}^{n+1} - U_j^{n+1}\right)$$

**Step 3:** The corrected upwind scheme (\(\mu > 0\)):

$$\boxed{-\frac{U_j^{n+1} - U_j^n}{k} + \mu_j^{n+1}\frac{U_{j+1}^{n+1} - U_j^{n+1}}{h} + b_j^{n+1}U_j^{n+1} = 0}$$

The denominator is \(h\), **not** \(2h\) as written in U13.

**Step 4:** Similarly for \(\mu < 0\), substitute \(\rho = -\mu h/2 = |\mu|h/2\):

$$\frac{\mu}{2h}(U_{j+1} - U_{j-1}) + \frac{|\mu|}{2h}(U_{j+1} - 2U_j + U_{j-1})$$

$$= \frac{\mu}{2h}(U_{j+1} - U_{j-1}) - \frac{\mu}{2h}(U_{j+1} - 2U_j + U_{j-1})$$

$$= \frac{\mu}{2h}\left[(U_{j+1}-U_{j-1}) - (U_{j+1}-2U_j+U_{j-1})\right] = \frac{\mu}{2h}(2U_j - 2U_{j-1}) = \frac{\mu}{h}(U_j - U_{j-1})$$

The corrected upwind scheme (\(\mu < 0\)):

$$\boxed{-\frac{U_j^{n+1} - U_j^n}{k} + \mu_j^{n+1}\frac{U_j^{n+1} - U_{j-1}^{n+1}}{h} + b_j^{n+1}U_j^{n+1} = 0}$$

Again, denominator \(h\), not \(2h\).

**VERIFICATION:**

- The numerical diffusion formula U15 states \(\frac{1}{2}\mu h \frac{\partial^2 V}{\partial S^2}\). A forward difference \(\frac{U_{j+1}-U_j}{h}\) has truncation error \(\frac{h}{2}V''\), giving diffusion \(\frac{1}{2}\mu h V''\). This is consistent with denominator \(h\). With denominator \(2h\), the truncation error would be \(\frac{h}{4}V''\), giving \(\frac{1}{4}\mu h V''\), inconsistent with U15. ✓
- Note-1 (§3.5, Algorithm 2) gives the upwind fully implicit matrix with terms \(\frac{rS_j}{\Delta S}\) (not \(\frac{rS_j}{2\Delta S}\)) in the off-diagonals, confirming denominator \(\Delta S = h\). ✓

**STATUS:** Resolved. **Downstream impact:** The numerical diffusion C12/U15/E(7)/C14 are all consistent with the corrected denominator \(h\). Only the displayed equations U13 and U14 contain the typo; the rest of the paper uses the correct expressions.

---

### GAP G3.3 + G7.1: Convergence bound E(6) unproved

**RESOLUTION TYPE:** sourced + derived (partial)

**RESOLUTION:**

E(6) states: \(|V(S_j, t_n) - U_j^n| \leq c(h+k)\) with \(c\) independent of \(h, k, \sigma\).

**Proof strategy (standard for exponentially fitted schemes on singularly perturbed problems):**

**Step 1: Uniform stability.** From E(5): \(\|A^{-1}\|_\infty \leq \frac{1}{1+kr} < 1\), independent of \(\sigma\). This provides unconditional \(\ell_\infty\) stability.

**Step 2: Local truncation error.** The fitted operator's truncation error must be bounded uniformly in \(\sigma\). For the centered difference scheme (\(\rho = \sigma\)), the spatial error is \(O(h^2)\). For the upwind limit (\(\rho = |\mu|h/2\)), the spatial error is \(O(h)\) — as shown in the G4.3 resolution, the forward difference has truncation error \(\frac{h}{2}V''\). For intermediate \(\sigma\), the fitting factor smoothly interpolates, and a Taylor expansion of the truncation error gives:

$$|\tau_j^n| \leq C_1 k + C_2 h$$

where \(C_1, C_2\) depend on the solution's derivatives but not on \(\sigma\), because the fitting factor absorbs the \(\sigma\)-dependence of the convection-diffusion balance.

**Step 3: Convergence via discrete stability.** The Lax-equivalence-type argument: if the scheme is stable (\(\|A^{-1}\|_\infty \leq 1\)) and the local truncation error is \(O(h+k)\), then the global error satisfies:

$$\|e^n\|_\infty \leq \frac{1}{r}\max_{m \leq n}\|\tau^m\|_\infty \leq \frac{C}{r}(h+k)$$

where \(\frac{1}{r}\) comes from the stability bound \(\|A^{-1}\|_\infty \leq \frac{1}{1+kr}\) accumulated over \(n = T/k\) steps.

**Source:** The theory of uniform convergence for exponentially fitted schemes is developed in Duffy [2] (the paper's reference), specifically Chapter 19 ("Exponentially Fitted Schemes"), which establishes first-order uniform convergence for the fitting factor E(3). The general framework is from Miller, O'Riordan, and Shishkin, *Fitted Numerical Methods for Singular Perturbation Problems* (World Scientific, 1996/2012), Chapters 8–9.

**Key limitation:** The bound \(O(h+k)\) (first order) applies uniformly. For fixed \(\sigma > 0\), the scheme is actually \(O(h^2 + k)\) since the fitting factor \(\rho \to \sigma\) and the centered difference is recovered. The \(O(h+k)\) bound is the worst-case (over all \(\sigma\)) result.

**VERIFICATION:**

- Consistency with E(5): the stability constant \(1/(1+kr)\) is independent of \(\sigma\). ✓
- Limiting cases: as \(\sigma \to \infty\), \(\rho \to \sigma\), spatial error \(\to O(h^2)\); as \(\sigma \to 0\), \(\rho \to |\mu|h/2\), spatial error \(\to O(h)\). The bound \(O(h+k)\) covers both. ✓
- Consistent with Note-2 §4.4 (convergence claim for the mixed method/implicit scheme). ✓

**STATUS:** Partially resolved. The proof strategy is established and the key steps are verified, but a complete line-by-line proof of the uniform truncation error bound requires access to Duffy [2] Chapter 19 or the Miller–O'Riordan–Shishkin framework, which is not included in the supplementary materials.

---

### GAP G6.2: Variable-coefficient case unaddressed

**RESOLUTION TYPE:** generalized

**RESOLUTION:**

The paper states \(r = r(t,S)\) and \(\sigma = \sigma(t,S)\) in §1 but derives everything for constant coefficients. Extension to variable coefficients requires:

**Exponentially fitted scheme (§2):** The formulation is already general. The fitting factor E(3) evaluates \(\mu_j^{n+1}\) and \(\sigma_j^{n+1}\) at each node and time level, so it naturally handles \(\sigma(S,t)\) and \(r(S,t)\). The M-matrix structure is verified node-by-node: at each interior node \(j\), the sign conditions on U7 (corrected per G4.1) hold provided \(\sigma_j^{n+1} > 0\), which is assured whenever \(\sigma(S_j, t_{n+1}) > 0\). The stability bound E(5) generalizes to:

$$\|A^{-1}\|_\infty \leq \frac{1}{1 + k\min_j |b_j^{n+1}|} = \frac{1}{1 + k\min_j r(S_j, t_{n+1})}$$

assuming \(b = -r < 0\) everywhere.

**CN variant (§3):** The parameter choice E(9) requires modification. With variable \(r(S)\) and \(\sigma(S)\), the M-matrix condition on the sub-diagonal of \(P\) becomes node-dependent:

$$r(S_j)\omega_2(S_j) + \frac{r(S_j)S_j}{4\Delta S} - \left(\frac{\sigma(S_j)S_j}{2\Delta S}\right)^2 \leq 0$$

**Option 1 (node-dependent \(\omega\)):** Set

$$\omega_2(S_j) = -\frac{r(S_j)}{16\sigma^2(S_j)}, \qquad j = 1, \ldots, M-1$$

This preserves the completing-the-square structure at each node but makes the matrices node-dependent (which they already are through \(S_j\)).

**Option 2 (conservative global \(\omega\)):** Set

$$\omega_2 = -\frac{\max_j r(S_j)}{16 \min_j \sigma^2(S_j)}$$

This is a worst-case choice ensuring the M-matrix condition holds globally but introduces more numerical diffusion than necessary.

The \(\Delta t\) constraint from the diagonal of \(N\) becomes:

$$\Delta t < \frac{1}{\max_j\left[r(S_j)\!\left(\frac{1}{2} - 2\omega_1(S_j)\right)\right] + \frac{1}{2}\max_j\!\left(\frac{\sigma(S_j)S_j}{\Delta S}\right)^2}$$

**VERIFICATION:**

- For constant coefficients, both options reduce to E(9). ✓
- The completing-the-square identity used in G1.5 applies identically at each node with local \(r(S_j)\) and \(\sigma(S_j)\). ✓
- Consistent with Note-1 §3.6 Algorithm 3, where node-dependent coefficients are used. ✓

**STATUS:** Resolved. **Downstream impact:** No equations change for constant-coefficient implementations. For variable-coefficient implementations, the modified E(9) and node-dependent \(\omega\) must be used.

---

## Missing Steps (G1 Series)

### GAP G1.1: Rearrangement from U6 to U7

**RESOLUTION TYPE:** derived

**RESOLUTION:** Fully resolved in the G4.1 resolution above. The step-by-step algebra from \(L_k^h U_j^n = 0\) to the tridiagonal form \(AU^{n+1} = U^n\) is shown there.

**STATUS:** Resolved.

---

### GAP G1.2: Numerical diffusion of upwind scheme

**RESOLUTION TYPE:** derived

**RESOLUTION:**

The corrected upwind scheme (from G4.3, \(\mu > 0\)) is:

$$-\frac{U_j^{n+1} - U_j^n}{k} + \mu\frac{U_{j+1}^{n+1} - U_j^{n+1}}{h} + bU_j^{n+1} = 0$$

Taylor-expand the forward difference about \(S_j\):

$$\frac{U_{j+1} - U_j}{h} = V'(S_j) + \frac{h}{2}V''(S_j) + \frac{h^2}{6}V'''(S_j) + O(h^3)$$

Substituting:

$$-\frac{\partial V}{\partial t} + \mu\left[\frac{\partial V}{\partial S} + \frac{h}{2}\frac{\partial^2 V}{\partial S^2} + O(h^2)\right] + bV = 0$$

$$\implies -\frac{\partial V}{\partial t} + \mu\frac{\partial V}{\partial S} + \frac{1}{2}\mu h\frac{\partial^2 V}{\partial S^2} + bV = O(h^2)$$

The leading-order modified equation is E(7):

$$-\frac{\partial V}{\partial t} + \mu(S,t)\frac{\partial V}{\partial S} + \frac{1}{2}\mu(S,t)h\frac{\partial^2 V}{\partial S^2} + b(S,t)V = 0$$

The numerical diffusion term is \(\frac{1}{2}\mu(S,t)h\frac{\partial^2 V}{\partial S^2}\).

Specializing to B-S via U4 (\(\mu = rS\), \(h = \Delta S\)):

$$\text{Numerical diffusion} = \frac{1}{2}rS\Delta S\frac{\partial^2 V}{\partial S^2}$$

This confirms C14 and resolves U2.

**VERIFICATION:** The Taylor expansion is standard. Dimensional check: \([\mu][h] = (\text{length/time})(\text{length}) = \text{length}^2/\text{time}\), matching the diffusion coefficient dimension. ✓

**STATUS:** Resolved.

---

### GAP G1.3: Derivation of P and N matrices for CN variant

**RESOLUTION TYPE:** derived

**RESOLUTION:**

Starting from E(1), discretize using Crank–Nicolson averaging for convection and diffusion, and E(8) for the reaction term. The full discrete equation at node \(j\) is:

$$-\frac{U_j^{n+1} - U_j^n}{\Delta t} + \frac{rS_j}{4\Delta S}\!\left(U_{j+1}^{n+1} - U_{j-1}^{n+1} + U_{j+1}^n - U_{j-1}^n\right) + \frac{\sigma^2 S_j^2}{4\Delta S^2}\!\left(\delta^2 U_j^{n+1} + \delta^2 U_j^n\right)$$

$$-\,r\!\left[\omega_2\!\left(U_{j-1}^{n+1}+U_{j+1}^{n+1}\right) + \!\left(\tfrac{1}{2}-2\omega_2\right)U_j^{n+1} + \omega_1\!\left(U_{j-1}^n+U_{j+1}^n\right) + \!\left(\tfrac{1}{2}-2\omega_1\right)U_j^n\right] = 0$$

where \(\delta^2 U_j = U_{j+1} - 2U_j + U_{j-1}\).

Collecting all \(n{+}1\) coefficients and negating (convention: \(PU^{n+1} = NU^n\) with positive diagonal of \(P\)):

**Matrix \(P\) (negated \(n{+}1\) coefficients):**

| Entry | Convection | Diffusion | Reaction | Time | Total |
|:---|:---|:---|:---|:---|:---|
| Sub (\(U_{j-1}^{n+1}\)) | \(+\frac{rS_j}{4\Delta S}\) | \(-\frac{\sigma^2 S_j^2}{4\Delta S^2}\) | \(+r\omega_2\) | 0 | \(r\omega_2 + \frac{rS_j}{4\Delta S} - \left(\frac{\sigma S_j}{2\Delta S}\right)^2\) |
| Diag (\(U_j^{n+1}\)) | 0 | \(+\frac{\sigma^2 S_j^2}{2\Delta S^2}\) | \(+r(\frac{1}{2}-2\omega_2)\) | \(+\frac{1}{\Delta t}\) | \(\frac{1}{\Delta t} + \frac{1}{2}\!\left(\frac{\sigma S_j}{\Delta S}\right)^2 + r\!\left(\frac{1}{2}-2\omega_2\right)\) |
| Super (\(U_{j+1}^{n+1}\)) | \(-\frac{rS_j}{4\Delta S}\) | \(-\frac{\sigma^2 S_j^2}{4\Delta S^2}\) | \(+r\omega_2\) | 0 | \(r\omega_2 - \frac{rS_j}{4\Delta S} - \left(\frac{\sigma S_j}{2\Delta S}\right)^2\) |

**Matrix \(N\) (un-negated \(n\) coefficients):**

| Entry | Convection | Diffusion | Reaction | Time | Total |
|:---|:---|:---|:---|:---|:---|
| Sub (\(U_{j-1}^n\)) | \(-\frac{rS_j}{4\Delta S}\) | \(+\frac{\sigma^2 S_j^2}{4\Delta S^2}\) | \(-r\omega_1\) | 0 | \(-r\omega_1 - \frac{rS_j}{4\Delta S} + \left(\frac{\sigma S_j}{2\Delta S}\right)^2\) |
| Diag (\(U_j^n\)) | 0 | \(-\frac{\sigma^2 S_j^2}{2\Delta S^2}\) | \(-r(\frac{1}{2}-2\omega_1)\) | \(+\frac{1}{\Delta t}\) | \(\frac{1}{\Delta t} - \frac{1}{2}\!\left(\frac{\sigma S_j}{\Delta S}\right)^2 - r\!\left(\frac{1}{2}-2\omega_1\right)\) |
| Super (\(U_{j+1}^n\)) | \(+\frac{rS_j}{4\Delta S}\) | \(+\frac{\sigma^2 S_j^2}{4\Delta S^2}\) | \(-r\omega_1\) | 0 | \(-r\omega_1 + \frac{rS_j}{4\Delta S} + \left(\frac{\sigma S_j}{2\Delta S}\right)^2\) |

All six entries match U18 and U19 exactly.

**Row sum verification:**

For \(N\): All \(S_j\)-dependent terms (convection ±, diffusion ±, reaction \(\omega_1\) terms) cancel pairwise, leaving \(\frac{1}{\Delta t} - \frac{r}{2}\). This confirms C18.

For \(P\): Similarly, all \(S_j\)-dependent terms cancel, leaving \(\frac{1}{\Delta t} + \frac{r}{2}\). This confirms C19 via Windisch's bound.

**VERIFICATION:**

- Each of the 12 coefficient components (4 sources × 3 positions for each matrix) independently verified. ✓
- Row sums confirmed: \(\frac{1}{\Delta t} \pm \frac{r}{2}\), independent of \(j\). ✓
- Matches U18/U19 term-by-term. ✓

**STATUS:** Resolved.

---

### GAP G1.4: Artificial diffusion of CN variant

**RESOLUTION TYPE:** derived

**RESOLUTION:**

The modification E(8) replaces the standard CN reaction-term discretization \(-r \cdot \frac{1}{2}(U_j^n + U_j^{n+1})\) with the six-node approximation. The deviation introduces spatial error.

**Step 1: Taylor-expand the level-\(n\) part of E(8).**

$$\omega_1(U_{j-1}^n + U_{j+1}^n) + \left(\tfrac{1}{2}-2\omega_1\right)U_j^n$$

Using \(U_{j \pm 1}^n = V \pm V'\Delta S + \frac{1}{2}V''\Delta S^2 \pm \frac{1}{6}V'''\Delta S^3 + O(\Delta S^4)\):

$$U_{j-1}^n + U_{j+1}^n = 2V + V''\Delta S^2 + O(\Delta S^4)$$

Substituting:

$$\omega_1(2V + V''\Delta S^2) + \left(\tfrac{1}{2}-2\omega_1\right)V + O(\Delta S^4) = \tfrac{1}{2}V + \omega_1 V''\Delta S^2 + O(\Delta S^4)$$

**Step 2: Similarly for the level-\(n{+}1\) part:** \(\frac{1}{2}V + \omega_2 V''\Delta S^2 + O(\Delta S^4)\)

**Step 3: Total reaction-term approximation:**

$$\frac{1}{2}V + \omega_1 V''\Delta S^2 + \frac{1}{2}V + \omega_2 V''\Delta S^2 + O(\Delta S^4) = V + (\omega_1+\omega_2)\Delta S^2 V'' + O(\Delta S^4)$$

The standard CN reaction term gives just \(V\). The modification adds:

$$-r \cdot (\omega_1+\omega_2)\Delta S^2\frac{\partial^2 V}{\partial S^2}$$

**Step 4: Substitute \(\omega_1 = \omega_2 = -\frac{r}{16\sigma^2}\):**

$$-r \cdot \left(-\frac{r}{8\sigma^2}\right)\Delta S^2\frac{\partial^2 V}{\partial S^2} = \frac{r^2\Delta S^2}{8\sigma^2}\frac{\partial^2 V}{\partial S^2} = \frac{1}{8}\!\left(\frac{r}{\sigma}\Delta S\right)^{\!2}\frac{\partial^2 V}{\partial S^2}$$

This confirms U26/C22.

**VERIFICATION:**

- Limiting case: as \(\sigma \to \infty\), \(\omega \to 0\) and the artificial diffusion vanishes (standard CN recovered). ✓
- Dimensional check: \([r/\sigma]^2[\Delta S]^2 = \text{dimensionless}\) times \([V/S^2] = \text{money}/\text{price}^2\). ✓
- Consistent with E(10) (the modified equation for the CN variant). ✓

**STATUS:** Resolved.

---

### GAP G1.5: Derivation of \(\omega_1 = \omega_2 = -r/(16\sigma^2)\)

**RESOLUTION TYPE:** derived

**RESOLUTION:**

**Objective:** Choose \(\omega_2\) so that the sub-diagonal of \(P\) is \(\leq 0\) for all interior nodes.

The sub-diagonal (from G1.3) is:

$$P_{\text{sub},j} = r\omega_2 + \frac{rS_j}{4\Delta S} - \frac{\sigma^2 S_j^2}{4\Delta S^2}$$

**Step 1: Complete the square in \(S_j\).**

$$P_{\text{sub},j} = -\frac{\sigma^2}{4\Delta S^2}\!\left(S_j^2 - \frac{r\Delta S}{\sigma^2}S_j\right) + r\omega_2$$

$$= -\frac{\sigma^2}{4\Delta S^2}\!\left(S_j - \frac{r\Delta S}{2\sigma^2}\right)^{\!2} + \frac{\sigma^2}{4\Delta S^2}\cdot\frac{r^2\Delta S^2}{4\sigma^4} + r\omega_2$$

$$= -\frac{\sigma^2}{4\Delta S^2}\!\left(S_j - \frac{r\Delta S}{2\sigma^2}\right)^{\!2} + \frac{r^2}{16\sigma^2} + r\omega_2$$

**Step 2: The maximum of \(P_{\text{sub},j}\) over all \(S_j \geq 0\)** is achieved at \(S_j = \frac{r\Delta S}{2\sigma^2}\), where the squared term vanishes:

$$\max_{S_j} P_{\text{sub},j} = \frac{r^2}{16\sigma^2} + r\omega_2$$

**Step 3: Require \(\max P_{\text{sub},j} \leq 0\):**

$$\frac{r^2}{16\sigma^2} + r\omega_2 \leq 0 \implies \omega_2 \leq -\frac{r}{16\sigma^2}$$

**Step 4: The tightest choice (preserving accuracy):**

$$\boxed{\omega_2 = -\frac{r}{16\sigma^2}}$$

By symmetry of the argument (applied to the sub-diagonal of \(N\) requiring \(\geq 0\)), the same bound applies to \(\omega_1\), giving \(\omega_1 = -\frac{r}{16\sigma^2}\).

**Step 5: \(\Delta t\) constraint.** The diagonal of \(N\) must be \(\geq 0\). Its minimum occurs at \(j = M-1\) (largest \(S_j\)):

$$\frac{1}{\Delta t} - \frac{(\sigma M)^2}{2} - r\!\left(\frac{1}{2} + \frac{r}{8\sigma^2}\right) \geq 0$$

$$\implies \Delta t < \frac{1}{r\!\left(\frac{1}{2} - 2\omega_1\right) + \frac{1}{2}(\sigma M)^2}$$

This is E(9).

**VERIFICATION:**

- At the critical \(S_j = r\Delta S/(2\sigma^2)\), \(P_{\text{sub},j} = 0\) exactly. ✓
- For Example 4.1 (\(r=0.05\), \(\sigma=0.001\), \(\Delta S = 0.05\)): critical \(S_j = 0.05 \times 0.05/(2 \times 10^{-6}) = 1250 > S_{\max} = 140\), so the sub-diagonal is strictly negative on the grid. ✓
- Consistent with Note-1 §3.6 which uses \(a = -r/(8\sigma^2)\) for a slightly different scheme (the factor of 2 difference arises from the different discretization structure). ✓

**STATUS:** Resolved.

---

### GAP G1.6: Bridge from exp-fitted to numerical diffusion for small finite \(\sigma\)

**RESOLUTION TYPE:** derived

**RESOLUTION:**

The fitted scheme uses \(\rho_j\) (E(3)) in place of \(\sigma_j = \frac{1}{2}\sigma^2 S_j^2\) for the diffusion coefficient. The effective numerical diffusion is \((\rho_j - \sigma_j)\frac{\partial^2 V}{\partial S^2}\).

**Step 1: Asymptotic expansion of \(\coth\) for large argument.**

Define the Péclet number \(\text{Pe}_j = \frac{\mu_j h}{2\sigma_j} = \frac{rS_j\Delta S}{\sigma^2 S_j^2} = \frac{r\Delta S}{\sigma^2 S_j}\).

For \(\text{Pe}_j \gg 1\): \(\coth(\text{Pe}_j) = 1 + \frac{2}{e^{2\text{Pe}_j}-1} = 1 + O(e^{-2\text{Pe}_j})\).

**Step 2: Compute \(\rho_j - \sigma_j\).**

$$\rho_j = \frac{\mu_j h}{2}\coth(\text{Pe}_j) = \sigma_j \cdot \text{Pe}_j \cdot \coth(\text{Pe}_j) \approx \sigma_j \cdot \text{Pe}_j = \frac{\mu_j h}{2}$$

$$\rho_j - \sigma_j \approx \frac{\mu_j h}{2} - \sigma_j = \frac{rS_j \Delta S}{2} - \frac{\sigma^2 S_j^2}{2} = \frac{S_j}{2}(r\Delta S - \sigma^2 S_j)$$

**Step 3: When \(\sigma^2 S_j \ll r\Delta S\):**

$$\rho_j - \sigma_j \approx \frac{1}{2}rS_j\Delta S$$

This confirms C14: the artificial diffusion is \(\frac{1}{2}rS\Delta S\frac{\partial^2 V}{\partial S^2}\).

**Quantitative condition (resolving G2.5):** The approximation requires \(\text{Pe}_j \gg 1\), i.e.:

$$\sigma^2 S_j \ll r\Delta S$$

This is stronger than just "\(\sigma \to 0\)" — it depends on the grid point \(S_j\) and the spatial step \(\Delta S\).

**VERIFICATION:**

- Limiting case \(\sigma \to 0\): \(\text{Pe}_j \to \infty\), \(\rho_j \to \frac{\mu_j h}{2}\), exact match with U12. ✓
- Limiting case \(\text{Pe}_j \to 0\) (\(\sigma \to \infty\)): \(\coth(x) \approx 1/x + x/3\), so \(\rho_j \approx \sigma_j + \frac{\mu_j^2 h^2}{12\sigma_j}\), and numerical diffusion \(\approx \frac{\mu_j^2 h^2}{12\sigma_j}\frac{\partial^2 V}{\partial S^2}\), which vanishes relative to \(\sigma_j\). ✓

**STATUS:** Resolved.

---

## Unstated Assumptions (G2 Series)

### GAP G2.1: \(h = \Delta S\) and \(k = \Delta t\) equivalence

**RESOLUTION TYPE:** disambiguated

**RESOLUTION:** Section 2 uses \((h, k)\) for the spatial and temporal steps; Sections 3–4 use \((\Delta S, \Delta t)\). These are identical: \(h \equiv \Delta S\) and \(k \equiv \Delta t\). The equivalence is confirmed by the fact that both refer to the same uniform grid applied to the same PDE, and the matrix structures in §2 and §§3–4 are consistent under this identification. The companion paper Note-2 §3.6 explicitly performs this identification: "\(\lambda = 1/\Delta t\)" and uses both notations interchangeably.

**STATUS:** Resolved.

---

### GAP G2.5: Quantitative condition for exp-fitted diffusion

**RESOLUTION TYPE:** derived

**RESOLUTION:** Resolved in G1.6 above. The required condition is \(\text{Pe}_j = \frac{r\Delta S}{\sigma^2 S_j} \gg 1\), i.e., \(\sigma^2 S_j \ll r\Delta S\).

**STATUS:** Resolved.

---

### GAP G2.6: E(10) validity condition

**RESOLUTION TYPE:** derived

**RESOLUTION:** E(10) drops the physical diffusion \(\frac{1}{2}\sigma^2 S^2\frac{\partial^2 V}{\partial S^2}\) and retains only the numerical diffusion \(\frac{1}{8}(r\Delta S/\sigma)^2\frac{\partial^2 V}{\partial S^2}\). The full modified equation is:

$$-\frac{\partial V}{\partial t} + rS\frac{\partial V}{\partial S} + \left[\frac{1}{2}\sigma^2 S^2 + \frac{1}{8}\!\left(\frac{r\Delta S}{\sigma}\right)^{\!2}\right]\frac{\partial^2 V}{\partial S^2} - rV = 0$$

E(10) is valid when:

$$\frac{1}{8}\!\left(\frac{r\Delta S}{\sigma}\right)^{\!2} \gg \frac{1}{2}\sigma^2 S^2 \quad\implies\quad \frac{r^2\Delta S^2}{4\sigma^4 S^2} \gg 1 \quad\implies\quad \sigma^2 S \ll \frac{r\Delta S}{2}$$

This is essentially the same large-Péclet condition as G2.5, confirming the paper's qualification "when \(\sigma S \to 0\)."

**STATUS:** Resolved.

---

## External Dependencies (G3 Series)

### GAP G3.1: M-matrix theory

**RESOLUTION TYPE:** sourced

**RESOLUTION:**

The results used from Ortega [5] are:

**Definition:** A real square matrix \(A\) is a **Z-matrix** if \(a_{ij} \leq 0\) for all \(i \neq j\).

**Theorem (Ortega §6.2.3, §6.2.17):** Let \(A\) be a Z-matrix. If \(A\) is irreducibly diagonally dominant (i.e., \(|a_{ii}| \geq \sum_{j \neq i}|a_{ij}|\) for all \(i\), with strict inequality for at least one \(i\), and the associated directed graph is strongly connected), then \(A\) is a nonsingular M-matrix, and \(A^{-1} > 0\) (strictly positive entries).

For a tridiagonal Z-matrix, irreducibility holds if and only if all sub-diagonal and super-diagonal entries are nonzero. The sign conditions U8 (\(a_{i,i+1} < 0\), \(a_{i+1,i} < 0\), \(a_{i,i} > 0\)) combined with strict diagonal dominance (row sum \(= 1 + kr > 0\)) give an irreducibly diagonally dominant Z-matrix, hence a nonsingular M-matrix with \(A^{-1} > 0\).

**STATUS:** Resolved.

---

### GAP G3.2: Windisch bound

**RESOLUTION TYPE:** sourced + derived

**RESOLUTION:**

**Theorem (Windisch [9]):** For a nonsingular M-matrix \(A\) with \(a_{ij} \leq 0\) for \(i \neq j\):

$$\|A^{-1}\|_\infty \leq \frac{1}{\min_i\!\left(a_{ii} - \sum_{j \neq i}|a_{ij}|\right)}$$

**Application to the exp-fitted scheme:** The diagonal dominance margin at row \(i\) is the row sum (since off-diagonals are negative):

$$a_{ii} - |a_{i,i-1}| - |a_{i,i+1}| = a_{ii} + a_{i,i-1} + a_{i,i+1} = 1 - kb = 1 + kr$$

This is constant across rows. Therefore \(\|A^{-1}\|_\infty \leq \frac{1}{1+kr}\).

**Application to the CN variant:** The diagonal dominance margin of \(P\) is:

$$P_{jj} + P_{j,j-1} + P_{j,j+1} = \frac{1}{\Delta t} + \frac{r}{2}$$

(constant across rows, verified in G1.3). Therefore \(\|P^{-1}\|_\infty \leq \left(\frac{1}{\Delta t} + \frac{r}{2}\right)^{-1}\).

**STATUS:** Resolved.

---

### GAP G3.4: Characteristic diffusion time

**RESOLUTION TYPE:** sourced

**RESOLUTION:** From Tavella–Randall [8], the characteristic diffusion time is \(\tau_d = \frac{\Delta S^2}{(\sigma S)^2}\). This represents the timescale over which diffusion acts across one spatial cell. When \(\Delta t \gg \tau_d\), the time step exceeds the natural diffusion timescale, and the scheme may not adequately resolve diffusive dynamics near discontinuities. For the paper's low-volatility examples (\(\sigma = 0.001\), \(S \approx 100\)): \(\tau_d = 0.05^2/(0.001 \times 100)^2 = 0.25\). With \(\Delta t = 0.01\): \(\Delta t/\tau_d = 0.04 \ll 1\), so the diffusion criterion is satisfied; the issue is convection dominance, not diffusion-timescale mismatch.

**STATUS:** Resolved.

---

### GAP G3.5: \(L^2\) vs \(L^\infty\) convergence

**RESOLUTION TYPE:** sourced

**RESOLUTION:** Giles and Carter [3] establish that for the Black–Scholes equation with non-smooth initial data (e.g., the kink at the strike of a vanilla call), the standard Crank–Nicolson scheme converges in \(L^2\) but fails to converge in \(L^\infty\) (supremum norm). This means that while the integrated squared error decreases with grid refinement, pointwise oscillations near the non-smoothness persist. For barrier options with jump discontinuities (even sharper than kinks), the \(L^\infty\) non-convergence is more severe, which is why the paper's C25 claim (oscillations "for every choice of steps") is consistent with this result, even though the paper only demonstrates it for one parameter set.

**STATUS:** Resolved.

---

## Suspected Typos (G4 Series, remaining)

### GAP G4.2: Submultiplicativity in U11

**RESOLUTION TYPE:** corrected

**RESOLUTION:**

U11 writes: \(\|A^{-1}U^n\|_\infty = \|A^{-1}\|_\infty\|U^n\|_\infty\).

The correct relationship is:

$$\|A^{-1}U^n\|_\infty \leq \|A^{-1}\|_\infty\|U^n\|_\infty$$

The second "\(=\)" should be "\(\leq\)" by the submultiplicativity property of the operator norm. The conclusion \(\|U^{n+1}\|_\infty \leq \|U^n\|_\infty\) is unaffected since the chain of inequalities still holds.

**STATUS:** Resolved. No downstream impact.

---

### GAP G4.4: Submultiplicativity in U24

**RESOLUTION TYPE:** corrected

**RESOLUTION:** Same issue as G4.2. U24 writes \(\|(P^{-1}N)U^n\|_\infty = \|P^{-1}\|_\infty\|N\|_\infty\|U^n\|_\infty\). The "\(=\)" should be "\(\leq\)." Conclusion C20 is unaffected.

**STATUS:** Resolved. No downstream impact.

---

### GAP G4.5: Strict vs. non-strict positivity

**RESOLUTION TYPE:** corrected

**RESOLUTION:**

The paper assumes \(U^0 \geq 0\) (A10) and concludes \((A^{-1})^n U^0 > 0\) (U10).

**Corrected statement:** Since \(A^{-1} > 0\) (all entries strictly positive, from C4), for any vector \(U^0 \geq 0\):

- If \(U^0\) has at least one positive component, then \(A^{-1}U^0 > 0\) (strictly), because each component of \(A^{-1}U^0\) is a sum of strictly positive terms times non-negative values, at least one of which is positive.
- If \(U^0 = 0\), then \(A^{-1}U^0 = 0\).

The correct conclusion is:
$$(A^{-1})^n U^0 \geq 0 \text{ for all } U^0 \geq 0; \quad (A^{-1})^n U^0 > 0 \text{ if } U^0 \geq 0 \text{ and } U^0 \neq 0$$

For option pricing, \(U^0\) is the payoff vector, which is nonzero (at least one node has positive payoff), so strict positivity holds after one time step.

**STATUS:** Resolved. No downstream impact on the paper's conclusions.

---

## Ambiguities (G5 Series)

### GAP G5.1: "Identically equal to 1"

**RESOLUTION TYPE:** disambiguated

**RESOLUTION:** The statement "This factor is identically equal to 1 in the centred difference scheme" means \(\rho_j / \sigma_j = 1\), i.e., \(\rho_j = \sigma_j = \sigma(S_j, t)\), the original diffusion coefficient. It does **not** mean \(\rho = 1\) numerically. The centered difference scheme corresponds to using the unmodified diffusion coefficient, which is the \(\rho/\sigma = 1\) case of the fitted scheme. This can be verified from E(3): for small Péclet (\(\mu h/(2\sigma) \to 0\)), \(\coth(x) \approx 1/x + x/3\), so \(\rho = \frac{\mu h}{2}\cdot\frac{2\sigma}{\mu h}(1 + O((\mu h/\sigma)^2)) = \sigma(1 + O(\text{Pe}^2))\), confirming \(\rho \to \sigma\).

**STATUS:** Resolved.

---

### GAP G5.2: E(10) omits physical diffusion

**RESOLUTION TYPE:** disambiguated

**RESOLUTION:** Fully resolved in G2.6 above. The full modified equation includes both physical and numerical diffusion; E(10) is the limiting form when numerical diffusion dominates, valid when \(\sigma^2 S \ll r\Delta S/2\).

**STATUS:** Resolved.

---

### GAP G5.3: Symbol \(\sigma\) overloading

**RESOLUTION TYPE:** disambiguated

**RESOLUTION:** In E(2), \(\sigma(S,t)\) denotes the general diffusion coefficient. In E(1) and §§3–4, \(\sigma\) denotes the volatility parameter. The connection is U4: \(\sigma(S,t) = \frac{1}{2}\sigma^2 S^2\). For implementation: in §2 formulas (U6, E(3), U7), use the general coefficient \(\sigma(S,t)\). In §§3–4 formulas (U18, U19, E(9)), use the volatility parameter \(\sigma\) directly (with \(\sigma^2 S_j^2/(4\Delta S^2)\) appearing as \((\sigma S_j/(2\Delta S))^2\)). The context is always unambiguous from the section.

**STATUS:** Resolved.

---

### GAP G5.4: "For every choice" universal quantifier

**RESOLUTION TYPE:** disambiguated

**RESOLUTION:** C25 should be interpreted as: "For every choice of \(\Delta S\) and \(\Delta t\), *when the standard CN scheme is applied directly to discontinuous initial data without smoothing*, spurious oscillations occur in the \(L^\infty\) sense." This is consistent with Giles–Carter [3] (G3.5): CN converges in \(L^2\) but not \(L^\infty\) for non-smooth data. The "every" qualifier does not apply if smoothing (e.g., Rannacher steps) is used first.

**STATUS:** Resolved.

---

### GAP G5.5: Time convention

**RESOLUTION TYPE:** disambiguated

**RESOLUTION:** The sign \(-\partial V/\partial t\) in E(1) and the "initial condition" being the payoff at \(t=0\) confirm that \(t\) is **time-to-expiry** (backward variable): \(t=0\) corresponds to maturity (payoff known), and the scheme marches forward in \(t\) toward the present (\(t=T\)). This is consistent with the monitoring update E(13) where \(0 = t_0 < t_1 < \cdots < t_F = T\), and with the boundary condition E(12). The convention is standard in the finite-difference option pricing literature (cf. Duffy [2], Tavella–Randall [8]).

**STATUS:** Resolved.

---

## Incompleteness (G6 Series, remaining)

### GAP G6.4: Behavior when E(9) is violated

**RESOLUTION TYPE:** derived

**RESOLUTION:** Largely resolved in the G7.3 resolution above. Summary of findings:

When E(9) is violated, the theoretical guarantees (positivity via \(N \geq 0\), max principle via \(\|P^{-1}\|_\infty\|N\|_\infty < 1\)) formally break down. However:

**(a)** The matrix \(P\) remains an M-matrix regardless of \(\Delta t\) (the off-diagonal conditions depend only on \(\omega_2\), not on \(\Delta t\)).

**(b)** The row sum of \(N\) is always \(\frac{1}{\Delta t} - \frac{r}{2} > 0\) (for \(\Delta t < 2/r\)), so \(N\) is not "catastrophically" negative.

**(c)** The diagonal of \(N\) becomes negative only at nodes with large \(S_j\) (near \(S_{\max}\)), where the solution is already near zero for the paper's barrier option examples.

**(d)** The scheme remains stable (in the sense that \(\rho(P^{-1}N) < 1\)) for a wider range of \(\Delta t\) than E(9) permits, because E(9) is a sufficient condition derived from norm inequalities, not a sharp spectral condition.

E(9) is best understood as a **sufficient condition for unconditional positivity preservation**, not as a stability boundary.

**STATUS:** Resolved.

---

## Proof Gaps (G7 Series, remaining)

### GAP G7.2: Discretization error \(O(\Delta S^2, \Delta t^2)\) of CN variant

**RESOLUTION TYPE:** derived

**RESOLUTION:**

The CN variant differs from standard CN only in the reaction-term discretization (E(8)). The standard CN scheme has truncation error \(O(\Delta S^2 + \Delta t^2)\).

**Step 1: Convection and diffusion terms.** These use standard CN averaging (half at level \(n\), half at \(n+1\)) with centered spatial differences, giving \(O(\Delta S^2 + \Delta t^2)\) truncation error (standard result).

**Step 2: Time derivative.** The backward difference \((U_j^{n+1} - U_j^n)/\Delta t\) approximates \(\partial V/\partial t\) at the midpoint \(t + \Delta t/2\) with error \(O(\Delta t^2)\) (standard result for midpoint evaluation).

**Step 3: Reaction term.** From G1.4, the six-node approximation gives:

$$V(t+\Delta t/2) = V + (\omega_1+\omega_2)\Delta S^2 V'' + O(\Delta S^4) + O(\Delta t^2)$$

The \(O(\Delta t^2)\) comes from the time-averaging (half at \(n\), half at \(n+1\) evaluates at the midpoint to second order). The spatial error is \(O(\Delta S^2)\) with coefficient \(|\omega_1+\omega_2| = \frac{r}{8\sigma^2}\).

**Step 4: Overall truncation error:**

$$\tau_j^n = O(\Delta t^2) + C\Delta S^2, \qquad C = O\!\left(\frac{r^2}{\sigma^2}\max|V''|\right)$$

The formal order is \(O(\Delta S^2 + \Delta t^2)\), but the coefficient \(C\) of the \(\Delta S^2\) term grows as \(r^2/(8\sigma^2)\), which is large for small \(\sigma\). Thus the order is technically correct but practically misleading: for \(\sigma = 0.001\) and \(r = 0.5\), the coefficient is \(0.25/(8 \times 10^{-6}) = 31250\), requiring very small \(\Delta S\) for acceptable accuracy.

**VERIFICATION:**

- Standard CN components: well-established \(O(\Delta S^2 + \Delta t^2)\). ✓
- Modification contributes \(O(\Delta S^2)\) with large coefficient, consistent with C22 and E(10). ✓
- For \(\sigma = O(1)\): \(|\omega| = O(r/\sigma^2) = O(1)\), and the modification has \(O(\Delta S^2)\) with moderate coefficient — standard second-order behavior. ✓

**STATUS:** Resolved. **Downstream impact:** C24 (accuracy constraints) follows directly from the \(\sigma\)-dependent coefficient.

---

### GAP G7.4: Irreducibility of \(P\) at critical \(S_j\)

**RESOLUTION TYPE:** derived

**RESOLUTION:**

From G1.5, the sub-diagonal of \(P\) is:

$$P_{\text{sub},j} = -\left(\frac{\sigma S_j}{2\Delta S} - \frac{r}{4\sigma}\right)^{\!2} \leq 0$$

This equals zero when \(S_j = S^* := \frac{r\Delta S}{2\sigma^2}\).

**Analysis of when \(S^*\) falls on the grid:**

\(S^* = j^*\Delta S\) requires \(j^* = \frac{r}{2\sigma^2}\), which is an integer only for specific parameter combinations.

| Parameters | \(S^*\) | \(S_{\max}\) | On grid? |
|:---|:---|:---|:---|
| \(r=0.05, \sigma=0.001, \Delta S=0.05\) | 1250 | 140 | **No** (outside domain) |
| \(r=0.5, \sigma=0.001, \Delta S=0.05\) | 12500 | 140 | **No** (outside domain) |
| \(r=0.5, \sigma=0.5, \Delta S=1.0\) | 1.0 | — | **Possibly** (depends on \(S_{\max}\)) |

For the paper's low-volatility examples, \(S^* \gg S_{\max}\), so the sub-diagonal is strictly negative on the entire grid, and \(P\) is irreducible.

**When \(S^*\) does fall on the grid:** \(P\) becomes reducible at row \(j^*\). However, \(P\) remains nonsingular (diagonal is strictly positive, and the determinant can be shown positive by the leading-principal-minor recurrence). The conclusion \(P^{-1} > 0\) (strict) would weaken to \(P^{-1} \geq 0\), but for positivity preservation (\(V^{n+1} \geq 0\)), \(P^{-1} \geq 0\) suffices.

**Implementer safeguard:** Check whether \(S^* \leq S_{\max}\). If so, either perturb \(\omega_2\) slightly (e.g., \(\omega_2 = -r/(16\sigma^2) - \epsilon\) for small \(\epsilon > 0\)) to ensure strict negativity, or verify that the numerical solution remains acceptable.

**STATUS:** Resolved.

---

## Summary Table

| Gap ID | Status | Downstream Impact |
|:---|:---|:---|
| **G1.1** | Resolved | Fills step U6 → U7/E(4); confirms matrix structure |
| **G1.2** | Resolved | Confirms U15, C12, E(7), C14 |
| **G1.3** | Resolved | Confirms all entries of U18, U19; confirms C18, C19 |
| **G1.4** | Resolved | Confirms U26, C22, E(10) |
| **G1.5** | Resolved | Confirms E(9) parameter choice |
| **G1.6** | Resolved | Confirms C14 for finite \(\sigma\); fills bridge to G2.5 |
| **G2.1** | Resolved | Notation only; no equation changes |
| **G2.2** | Resolved | Enables grid construction for implementation |
| **G2.3** | Resolved | Defines \(M = S_{\max}/\Delta S\) |
| **G2.4** | Resolved | Specifies Dirichlet BCs; tridiagonal system unchanged for zero BCs |
| **G2.5** | Resolved | Quantitative condition \(\sigma^2 S \ll r\Delta S\) |
| **G2.6** | Resolved | Condition \(\sigma^2 S \ll r\Delta S/2\) for E(10) validity |
| **G3.1** | Resolved | Restates M-matrix theorem supporting C3, C4, C16 |
| **G3.2** | Resolved | Restates Windisch bound supporting E(5), C19 |
| **G3.3** | Partially resolved | Proof strategy established; full proof requires Duffy [2] Ch. 19 |
| **G3.4** | Resolved | Defines \(\tau_d\); confirms not the binding constraint in examples |
| **G3.5** | Resolved | \(L^2\) convergence without \(L^\infty\) for non-smooth data |
| **G4.1** | Resolved | Corrects \(n \to n{+}1\) in U7; matters only for variable coefficients |
| **G4.2** | Resolved | Cosmetic: \(= \to \leq\) in U11; conclusion unaffected |
| **G4.3** | Resolved | Corrects \(2h \to h\) in U13/U14; consistent with U15, C12 |
| **G4.4** | Resolved | Cosmetic: \(= \to \leq\) in U24; conclusion unaffected |
| **G4.5** | Resolved | Clarifies \(> 0\) requires \(U^0 \neq 0\); conclusion unaffected |
| **G5.1** | Resolved | \(\rho/\sigma = 1\), not \(\rho = 1\) |
| **G5.2** | Resolved | Full modified equation given; E(10) is limiting form |
| **G5.3** | Resolved | Context disambiguates; §2 = general coeff, §§3–4 = volatility |
| **G5.4** | Resolved | "Every" applies to direct CN without smoothing |
| **G5.5** | Resolved | \(t\) = time-to-expiry (backward variable) |
| **G6.1** | Resolved | Homogeneous Dirichlet BCs; no system modification needed |
| **G6.2** | Resolved | Node-dependent \(\omega_j\) or conservative global choice |
| **G6.3** | Resolved | \(S_j = j\Delta S\), \(M = S_{\max}/\Delta S\) |
| **G6.4** | Resolved | E(9) is sufficient, not necessary; scheme empirically works beyond it |
| **G7.1** | Partially resolved | Proof strategy given; complete proof needs Duffy [2] reference |
| **G7.2** | Resolved | \(O(\Delta S^2 + \Delta t^2)\) correct; coefficient \(\propto r^2/\sigma^2\) |
| **G7.3** | Resolved | Examples violate E(9); scheme works due to zero-boundary effect |
| **G7.4** | Resolved | Critical \(S^*\) outside grid for all examples; safeguard given |
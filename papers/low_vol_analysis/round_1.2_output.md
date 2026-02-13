
# Round 1.2 — Systematic Mathematical Audit

## 1. GAP REGISTRY

| Gap ID | Location (§/Eq/Claim) | Category | Description | Severity |
|:---|:---|:---|:---|:---|
| G1.1 | §2, U6/E(3) → E(4)/U7 | G1 — Missing Step | Algebraic rearrangement of the fitted operator set to zero into the tridiagonal form \(A\,U^{n+1}=U^n\) is omitted. I have verified the result is correct by performing the rearrangement explicitly. | Minor |
| G1.2 | §2, U13/U14 → U15/C12 | G1 — Missing Step | Taylor expansion establishing that the upwind finite difference introduces numerical diffusion \(\tfrac{1}{2}\mu\,h\,\partial^2 V/\partial S^2\) is never shown; paper states "through a standard analysis of consistency." I verified: expanding \((U_{j+1}-U_j)/h = V' + \tfrac{h}{2}V'' + O(h^2)\) gives the claimed term. | Minor |
| G1.3 | §3, E(1)+E(8) → U17/U18/U19 | G1 — Missing Step | Full derivation of the tridiagonal matrices \(P\) and \(N\) from the Crank–Nicolson discretization with six-node reaction term E(8) is omitted. I have verified all six matrix entries by explicit algebra (see audit notes). | Minor |
| G1.4 | §3, E(8)/E(9) → U26/C22 | G1 — Missing Step | Taylor expansion showing that the modified reaction-term discretization introduces artificial diffusion \(\tfrac{1}{8}(r\Delta S/\sigma)^2 \partial^2V/\partial S^2\) is not shown. I verified: the modification adds \(\omega(\Delta S^2)\partial^2V/\partial S^2\) at each time level, totalling \((\omega_1+\omega_2)r\Delta S^2\partial^2V/\partial S^2 = \tfrac{r^2\Delta S^2}{8\sigma^2}\partial^2V/\partial S^2\). | Minor |
| G1.5 | §3, U20a/U20b → E(9) | G1 — Missing Step | Derivation of \(\omega_1=\omega_2=-r/(16\sigma^2)\) from the M-matrix / nonnegativity criteria is omitted. I verified: sub-diagonal of \(P\) equals \(-(\sigma S_j/(2\Delta S) - r/(4\sigma))^2\), whose maximum over \(S_j\) is zero precisely at \(\omega_2=-r/(16\sigma^2)\). | Minor |
| G1.6 | §2, C14 claim vs. derivation | G1 — Missing Step | Paper derives numerical diffusion only for the upwind limit (\(\sigma\to 0\)) but applies it to the exp-fitted scheme for small finite \(\sigma\). The bridge requires showing \(\rho - \sigma(S,t) \approx \tfrac{1}{2}\mu h\) for large Péclet number, which involves an asymptotic expansion of \(\coth\) never performed. The result is standard and correct. | Minor |
| G2.1 | §2 vs. §§3–4 | G2 — Unstated Assumption | \(h=\Delta S\) and \(k=\Delta t\) are never stated to be equivalent. Section 2 uses \((h,k)\), Sections 3–4 use \((\Delta S,\Delta t)\). | Minor |
| G2.2 | Throughout §§2–4 | G2 — Unstated Assumption | A uniform spatial grid \(S_j = j\,\Delta S\), \(j=0,\dots,M\), is never formally stated. It is inferred from the tridiagonal structure and the appearance of \(S_j\) in matrix entries. | Significant |
| G2.3 | Throughout §§2–4 | G2 — Unstated Assumption | \(S_{\max}=M\Delta S\) is never defined in the body text. \(M\) appears first in E(9) without definition; \(S_{\max}\) appears only in figure captions. | Significant |
| G2.4 | Throughout §§2–4 | G2 — Unstated Assumption | Boundary treatment at \(S=0\) and \(S=S_{\max}\) for the tridiagonal systems is never specified. The matrices U7, U18/U19, U33 are written for interior nodes only. Implementation requires knowing how the first and last rows are modified. | Critical |
| G2.5 | §2, C14 | G2 — Unstated Assumption | The formula \(\tfrac{1}{2}rS\Delta S\) for the exp-fitted diffusion requires \(\sigma^2 S \ll r\Delta S\) (large Péclet number), not merely \(\sigma\to 0\). This quantitative condition is never stated. | Minor |
| G2.6 | §3, E(10) | G2 — Unstated Assumption | E(10) drops the original diffusion \(\tfrac{1}{2}\sigma^2 S^2\partial^2V/\partial S^2\) and retains only the numerical diffusion. This requires \(\sigma^2 S^2 \ll \tfrac{1}{4}(r\Delta S/\sigma)^2\), i.e., \(\sigma^4 S^2 \ll r^2\Delta S^2/4\). The paper says "when \(\sigma S\to 0\)" without quantifying. | Minor |
| G3.1 | §2, C3/C4 | G3 — External Dep. | M-matrix theory: \(A\) irreducible diag. dominant \(\Rightarrow\) nonsingular M-matrix \(\Rightarrow A^{-1}>0\). Relies on Ortega [5] §§6.2.3, 6.2.17 without restating the theorem. | Minor |
| G3.2 | §2 E(5), §3 C19 | G3 — External Dep. | Windisch [9] bound: for diag. dominant M-matrix, \(\|A^{-1}\|_\infty \le 1/\min_i(\text{diag dominance excess})\). Used twice, never restated. | Minor |
| G3.3 | §2, E(6)/C8 | G3 — External Dep. | Convergence bound \(\lvert V(S_j,t_n)-U_j^n\rvert \le c(h+k)\) with \(c\) independent of \(h,k,\sigma\) is stated without proof and without citing a specific theorem. This is the paper's central "uniform convergence" claim for the exp-fitted scheme. | Significant |
| G3.4 | §4, U34/C34 | G3 — External Dep. | Characteristic diffusion time \(\tau_d\) and oscillation criterion from Tavella–Randall [8]. Not restated. | Minor |
| G3.5 | §4, C35 | G3 — External Dep. | \(L^2\) convergence but not \(L^\infty\) convergence for non-smooth data from Giles–Carter [3]. Not restated. | Minor |
| G4.1 | §2, U7 | G4 — Suspected Typo | Matrix \(A\) writes coefficients with superscript \(n\) (e.g. \(\rho_j^n, \mu_j^n, b_j^n\)), but the fitted operator U6 evaluates all spatial terms at level \(n{+}1\) (e.g. \(\mu_j^{n+1}, \rho_j^{n+1}\)). For constant-coefficient B-S these coincide, but for variable coefficients the discrepancy matters. Should be \(n{+}1\). | Significant |
| G4.2 | §2, U11 | G4 — Suspected Typo | Chain writes \(\|A^{-1}U^n\|_\infty = \|A^{-1}\|_\infty\|U^n\|_\infty\). By sub-multiplicativity of operator norms, the second relation should be \(\le\), not \(=\). The conclusion \(\|U^{n+1}\|_\infty \le \|U^n\|_\infty\) is unaffected. | Minor |
| G4.3 | §2, U13/U14 | G4 — Suspected Typo | Both upwind schemes show denominator \(2h\) in the convection term. Derivation from U12 substituted into U6 yields denominator \(h\): combining centered \(\mu(U_{j+1}-U_{j-1})/(2h)\) with artificial diffusion \((\mu h/2)(U_{j+1}-2U_j+U_{j-1})/h^2\) gives \(\mu(U_{j+1}-U_j)/h\), not \(\mu(U_{j+1}-U_j)/(2h)\). Numerical diffusion formula U15 is consistent with \(h\), confirming the typo. | Significant |
| G4.4 | §3, U24 | G4 — Suspected Typo | Same submultiplicativity issue as G4.2: \(\|(P^{-1}N)U^n\|_\infty = \|P^{-1}\|_\infty\|N\|_\infty\|U^n\|_\infty\) should use \(\le\). Conclusion unaffected. | Minor |
| G4.5 | §2, U10/C6 | G4 — Suspected Typo | Conclusion states \((A^{-1})^n U^0 > 0\) (strict) but assumption A10 gives \(U^0\ge 0\) (non-strict). For \(U^0\ge 0\) with \(U^0\neq 0\), \(A^{-1}>0\) gives \(A^{-1}U^0 > 0\). For \(U^0 = 0\), result is \(=0\) not \(>0\). Should conclude \(\ge 0\) in general, \(> 0\) if \(U^0\) has at least one positive component. | Minor |
| G5.1 | §2, p.226 top (C1) | G5 — Ambiguity | "This factor is identically equal to 1 in the centred difference scheme." Means the ratio \(\rho/\sigma(S,t)=1\), i.e. \(\rho=\sigma(S,t)\). Could be misread as \(\rho=1\). | Minor |
| G5.2 | §3, E(10) | G5 — Ambiguity | E(10) is presented as the modified equation solved by the CN variant, but it omits the original physical diffusion \(\tfrac{1}{2}\sigma^2 S^2\partial^2V/\partial S^2\). The full modified equation is \(-\partial_t V + rS\partial_S V + [\tfrac{1}{2}\sigma^2 S^2 + \tfrac{1}{8}(r\Delta S/\sigma)^2]\partial_S^2 V - rV = 0\). E(10) is valid only in the limit where numerical diffusion dominates physical diffusion. | Minor |
| G5.3 | §§1–2 | G5 — Ambiguity | Symbol \(\sigma\) is overloaded: denotes both volatility parameter (scalar, §§1,3,4) and general diffusion coefficient \(\sigma(S,t)\) in E(2). Context disambiguates, but the identification U4 (\(\sigma(S,t)=\tfrac{1}{2}\sigma^2 S^2\)) means the same symbol appears on both sides of the equation. | Minor |
| G5.4 | §4, p.229 (C25) | G5 — Ambiguity | "The Crank-Nicolson scheme produces undesired spurious oscillations in the numerical solution for every choice of steps \(\Delta S\) and \(\Delta t\)." This is demonstrated empirically for one parameter set. It is likely true in the \(L^\infty\) sense (per Giles–Carter [3], C35), but the universal quantifier "every" is not proved. | Minor |
| G5.5 | §1, before E(1) | G5 — Ambiguity | "\(t\) is the time to expiry \(T\)" with \(0\le t\le T\). Could mean \(t\) is current calendar time with expiry at \(T\), or \(t\) is time remaining to expiry. The sign of \(-\partial V/\partial t\) and the initial condition at \(t=0\) being the payoff confirm \(t\) is time-to-expiry (backward variable). The paper never states this convention explicitly. | Minor |
| G6.1 | Throughout §§2–4 | G6 — Incompleteness | Boundary condition treatment for the finite difference schemes is never specified. The tridiagonal matrices are written for interior nodes only. For implementation, the first and last rows of the linear system require modification to incorporate Dirichlet conditions \(V(0,t)=0\), \(V(S_{\max},t)\approx 0\), or barrier conditions. | Critical |
| G6.2 | §1 vs. §§2–4 | G6 — Incompleteness | §1 states \(r=r(t,S)\) and \(\sigma=\sigma(t,S)\) (variable coefficients), but all subsequent analysis assumes constant \(r,\sigma\). The CN variant matrices U18/U19 and parameter choice E(9) are derived specifically for constant coefficients. No guidance for the variable-coefficient case is given. | Significant |
| G6.3 | Throughout | G6 — Incompleteness | Grid construction: \(S_{\max}\), relationship \(S_{\max}=M\Delta S\), choice of \(M\), and the uniform-grid assumption are not formally specified. Appears only implicitly through figure captions. | Significant |
| G6.4 | §3/§4 | G6 — Incompleteness | Behavior of the CN variant scheme when E(9)'s \(\Delta t\) constraint is violated is never discussed, although the paper's own examples violate it (see G7.3). In practice the scheme still produces non-oscillating solutions, suggesting E(9) is a sufficient condition whose violation doesn't necessarily cause failure. No analysis of this regime is provided. | Significant |
| G7.1 | §2, E(6)/C8 | G7 — Proof Gap | The convergence bound \(\lvert V(S_j,t_n)-U_j^n\rvert \le c(h+k)\) with \(c\) independent of \(\sigma\) is the paper's key "uniform convergence" claim. No proof is given, no specific theorem is cited. This is the foundational result justifying the exp-fitted approach. | Significant |
| G7.2 | §3, C21/U25 | G7 — Proof Gap | Discretization error \(O(\Delta S^2,\Delta t^2)\) for the CN variant is stated without derivation. While the standard CN achieves this order and the modification is \(O(\Delta S^2)\), the constant in the modified term is \(r^2/(8\sigma^2)\) which can dominate for small \(\sigma\). The order claim is technically correct but practically misleading without noting the \(\sigma\)-dependent constant. | Significant |
| G7.3 | §4, Examples 4.1, 4.4, 4.5 | G7 — Proof Gap | E(9)'s time step constraint is violated in most numerical examples. Specifically: Example 4.1 (\(r=0.05\)): E(9) requires \(\Delta t < 0.0032\), but \(\Delta t=0.01\) is used. Example 4.4 (\(r=0.5\)): E(9) requires \(\Delta t < 3.2\times 10^{-5}\), but \(\Delta t=0.01\) is used. Only Example 4.4 (\(r=0.01\)) satisfies E(9): requires \(\Delta t < 0.061\). The theoretical guarantees (positivity, max principle) do not apply to these examples, yet the paper does not acknowledge this. | Significant |
| G7.4 | §3, C16/U20a | G7 — Proof Gap | With \(\omega_2=-r/(16\sigma^2)\), the sub-diagonal of \(P\) equals \(-(\sigma S_j/(2\Delta S)-r/(4\sigma))^2 \le 0\) and is zero when \(S_j = r\Delta S/(2\sigma^2)\). If this value coincides with a grid point, \(P\) is reducible (zero sub-diagonal entry). The paper claims \(P\) is "irreducibly diagonally dominant." For the paper's parameters, \(r\Delta S/(2\sigma^2)\) falls far outside \([0,S_{\max}]\), so this is not triggered in the examples. | Minor |

---

## 2. VERIFICATION LOG

### §1 Foundations

| From → To | Status | Gap ID | Notes |
|:---|:---|:---|:---|
| A1 → E(1) | VERIFIED | — | Standard B-S PDE under risk-neutral measure. Sign convention consistent with \(t\) = time-to-expiry (G5.5 ambiguity noted). |
| E(1) → E(2) via U4 | VERIFIED | — | Substituting \(\sigma(S,t)=\tfrac12\sigma^2S^2\), \(\mu=rS\), \(b=-r\) into E(2) recovers E(1). |
| E(1) → U1 | VERIFIED | — | Generality statement; \(r,\sigma\) treated as constant hereafter (G6.2). |

### §2 Exponentially Fitted Scheme

| From → To | Status | Gap ID | Notes |
|:---|:---|:---|:---|
| E(2) → U5 | VERIFIED | — | Operator \(L\) is the left-hand side of E(2); definition only. |
| U5 → U6 | VERIFIED | — | Standard FD replacements: backward-Euler in time, centered in space, with fitting factor \(\rho\) replacing \(\sigma(S,t)\). |
| — → E(3) | VERIFIED | — | Duffy's definition of fitting factor. Limit \(\rho\to\sigma\) as Péclet → 0 confirmed via \(\coth(x)\approx 1/x\) for small \(x\). |
| U6 + E(3) → E(4) + U7 | GAP | G1.1, G4.1 | Rearrangement verified by explicit algebra (correct). Time-level superscript \(n\) in U7 should be \(n{+}1\) per U6. |
| U7 → U8 | VERIFIED | — | For \(\mu>0\): \(\rho>\mu h/2\) since \(\coth(x)>1\) for \(x>0\), giving sub-diagonal \(<0\). Super-diagonal \(<0\) trivially. Diagonal \(>0\) from \(\rho>0\), \(-b=r>0\), \(1/k>0\). Also verified for \(\mu<0\) and \(\mu=0\). |
| U8 → C2 | VERIFIED | — | Sign conditions hold for all \(\mu\) (not just \(\mu>0\)) given \(\sigma>0\). |
| C2 → C3 | VERIFIED | G3.1 | Irreducibility from nonzero off-diagonals; strict diagonal dominance because \(-b+1/k = r+1/k > 0\) uniformly. References Ortega [5]. |
| C3 → C4 | VERIFIED | G3.1 | Standard M-matrix result. |
| C3 + A7 + A12 → E(5)/C5 | VERIFIED | G3.2 | Diagonal dominance excess = row sum = \(1-kb = 1+kr\). Windisch bound gives \(\|A^{-1}\|_\infty \le 1/(1+kr) < 1\). |
| C4 + A10 → U10 → C6 | VERIFIED | G4.5 | Induction correct. Minor: conclusion \(>0\) requires \(U^0\) to have ≥ 1 positive component, not just \(U^0\ge 0\). |
| E(5) + C5 → U11 → C7 | VERIFIED | G4.2 | Second "=" should be "≤". Conclusion \(\|U^{n+1}\|_\infty \le \|U^n\|_\infty\) is correct. |
| — → E(6)/C8 | GAP | G3.3, G7.1 | Stated without proof or specific reference. Central uniform-convergence claim. |
| C3 + C5 → C9, C10 | VERIFIED | — | Uniform stability from E(5); oscillation-freedom from M-matrix positivity. Qualitative consequences of prior results. |
| E(3) → U12 | VERIFIED | — | \(\lim_{\sigma\to 0}\coth(\mu h/(2\sigma))=\pm 1\). Correctly computed for both signs of \(\mu\). |
| U12 + U6 → U13, U14 | GAP | G4.3 | Substitution gives denominator \(h\), not \(2h\). The algebraic combination of centered + artificial diffusion yields \(\mu(U_{j+1}-U_j)/h\) for \(\mu>0\). |
| U13/U14 → U15 → C12 | GAP | G1.2 | Taylor expansion omitted but is standard and I verified: \((U_{j+1}-U_j)/h = V'+\tfrac{h}{2}V''+O(h^2)\), giving numerical diffusion \(\tfrac{1}{2}\mu h V''\). |
| C12 → E(7) | VERIFIED | — | Adding numerical diffusion \(\tfrac{1}{2}\mu h\partial^2V/\partial S^2\) to the hyperbolic equation U16. |
| C12 + U4 → C14 (= U2) | VERIFIED | G2.5, G1.6 | Substituting \(\mu=rS\), \(h=\Delta S\). Valid approximation requires \(\sigma^2 S \ll r\Delta S\). |
| C12 → C13 | VERIFIED | — | Qualitative: diffusion \(\tfrac{1}{2}rS\Delta S\) grows with \(r\) and \(\Delta S\). |

### §3 Crank–Nicolson Variant

| From → To | Status | Gap ID | Notes |
|:---|:---|:---|:---|
| E(1) → E(8) | VERIFIED | — | Definition of six-node reaction-term approximation. Design choice. |
| E(1) + E(8) → U17/U18/U19 | GAP | G1.3 | All six matrix entries verified by explicit algebra. Omission is a missing step only. |
| U18/U19 + U20a/U20b → E(9) | GAP | G1.5 | Sub-diagonal of \(P\) is \(-(\sigma S_j/(2\Delta S)-r/(4\sigma))^2\), maximized at zero when \(\omega_2=-r/(16\sigma^2)\). \(\Delta t\) bound from requiring \(N\) diagonal \(\ge 0\) at \(S_j=S_{\max}=M\Delta S\). Verified. |
| E(9) → C16 | VERIFIED | G3.1, G3.2, G7.4 | \(P\) is M-matrix: diagonal dominance excess \(=1/\Delta t+r/2>0\) (unconditional). \(N\ge 0\): requires \(\Delta t\) bound. Irreducibility of \(P\) holds except at one specific \(S_j\) (G7.4). |
| C16 + A10 → U21 → C17 | VERIFIED | — | \(P^{-1}\ge 0\) and \(N\ge 0\) imply \(P^{-1}N\ge 0\), and \((P^{-1}N)^n U^0 \ge 0\). |
| U19 + E(9) → U22/C18 | VERIFIED | — | Row sum of \(N\): all \(S_j\)-dependent terms cancel; result \(=1/\Delta t - r/2\) for every row. Explicitly computed and confirmed. |
| U18 + E(9) → U23/C19 | VERIFIED | G3.2 | Row sum of \(P = 1/\Delta t + r/2\) (constant across rows). Windisch bound gives \(\|P^{-1}\|_\infty \le (1/\Delta t + r/2)^{-1}\). |
| C18 + C19 → U24 → C20 | VERIFIED | G4.4 | Same submultiplicativity issue as G4.2. Product \(\frac{1/\Delta t - r/2}{1/\Delta t + r/2} < 1\) since \(r>0\). |
| U17–U19 → C21/U25 | GAP | G7.2 | \(O(\Delta S^2,\Delta t^2)\) stated without derivation. Order is correct (modification is \(O(\Delta S^2)\)), but the constant includes \(r^2/(8\sigma^2)\) which is large for small \(\sigma\). |
| E(8) + E(9) → U26/C22 | GAP | G1.4 | Taylor expansion of modified reaction term: deviation from standard CN is \(-r(\omega_1+\omega_2)\Delta S^2 \partial^2V/\partial S^2 = \tfrac{r^2\Delta S^2}{8\sigma^2}\partial^2V/\partial S^2\). Verified. |
| C22 → E(10)/C23 | VERIFIED | G5.2, G2.6 | Valid when numerical diffusion \(\gg\) physical diffusion. Paper qualifies with "When \(\sigma S\to 0\)". |
| E(9) + C22 → U27/U28/C24 | VERIFIED | — | \(\tfrac{1}{8}(r\Delta S/\sigma)^2\) small requires small \(\Delta S\). From E(9): \(\Delta t \sim 8(\sigma/r)^2\). Algebra confirmed. |
| C14 + C22 → U29 | VERIFIED | — | Comparison of \(\tfrac{1}{2}rS\Delta S\) (Duffy) vs. \(\tfrac{1}{8}(r\Delta S/\sigma)^2\) (CN variant). |

### §4 Numerical Results

| From → To | Status | Gap ID | Notes |
|:---|:---|:---|:---|
| E(1) + centered FD → U32/U33 | VERIFIED | — | Standard fully implicit centered scheme. All three tridiagonal entries verified by explicit algebra. |
| U33 → C27 | VERIFIED | — | Sub-diagonal \(\ge 0\) when \(\sigma^2 S_j < r\Delta S\). For uniform grid, condition reduces to \(\sigma^2 < r\) at \(j=1\). Cancellation of \(\Delta S\) confirmed. |
| — → C25, C26 | VERIFIED (empirical) | G5.4 | Oscillations demonstrated in Figs. 1, 2. "For every choice" not proved; consistent with Giles–Carter \(L^\infty\) non-convergence (C35). |
| C14 → C28 | VERIFIED (empirical) | — | Higher \(r\) increases \(\tfrac{1}{2}rS\Delta S\). Figs. 3 confirm. |
| C22 → C29, C30 | VERIFIED (empirical) | — | Higher \(r\) increases \(\tfrac{1}{8}(r\Delta S/\sigma)^2\). Figs. 4 confirm. CN variant remains non-oscillating. |
| E(9) vs. examples → C31, C32 | GAP | G7.3 | Examples use \(\Delta t=0.01\) which violates E(9) for \(r=0.05\) (needs \(\Delta t<0.0032\)) and \(r=0.5\) (needs \(\Delta t<3.2\times10^{-5}\)). Only \(r=0.01\) case satisfies E(9). |
| — → C33 | VERIFIED | — | Qualitative conclusion: no universally optimal scheme. |
| — → U34/C34 | VERIFIED | G3.4 | Definition from Tavella–Randall [8]. |
| — → C35 | VERIFIED | G3.5 | External result from Giles–Carter [3]. |

---

## 3. PROOF AUDIT

### Proof of C6 (Positivity of numerical solution)

| Step | Statement | Status | Notes |
|:---|:---|:---|:---|
| 1 | From E(4): \(U^{n+1}=A^{-1}U^n\) | VERIFIED | Direct consequence of \(AU^{n+1}=U^n\) with \(A\) nonsingular. |
| 2 | \(A^{-1}>0\) (componentwise) | VERIFIED | From C4, which follows from C3 (M-matrix) via Ortega [5]. |
| 3 | By induction: \(U^n=(A^{-1})^nU^0>0\) | VERIFIED with caveat | Requires \(U^0>0\) (strictly). Paper assumes \(U^0\ge 0\) (A10). Strictly, should conclude \(\ge 0\). See G4.5. |
| **Verdict** | **Sound** | — | Fixable minor issue (strict vs. non-strict inequality). |

### Proof of C7 (Discrete maximum principle, exp-fitted scheme)

| Step | Statement | Status | Notes |
|:---|:---|:---|:---|
| 1 | \(\|U^{n+1}\|_\infty = \|A^{-1}U^n\|_\infty\) | VERIFIED | From E(4). |
| 2 | \(\|A^{-1}U^n\|_\infty \le \|A^{-1}\|_\infty\|U^n\|_\infty\) | VERIFIED | Sub-multiplicativity. Paper writes "=" (G4.2). |
| 3 | \(\|A^{-1}\|_\infty \le 1/(1+kr)<1\) | VERIFIED | From E(5)/C5. |
| 4 | Conclusion: \(\|U^{n+1}\|_\infty \le \|U^n\|_\infty\) | VERIFIED | Follows from steps 1–3. |
| **Verdict** | **Sound** | — | Cosmetic typo in step 2 (G4.2). |

### Proof of C11 (Limit σ→0 gives upwind)

| Step | Statement | Status | Notes |
|:---|:---|:---|:---|
| 1 | \(\lim_{\sigma\to 0}\rho=\|\mu\|h/2\) | VERIFIED | From \(\coth(x)\to\pm 1\) as \(x\to\pm\infty\). |
| 2 | Substitution into U6 yields U13 (\(\mu>0\)) and U14 (\(\mu<0\)) | VERIFIED with caveat | Algebra yields denominator \(h\), not \(2h\) as written (G4.3). |
| **Verdict** | **Sound** | — | Result correct; displayed equations have typo in denominator. |

### Proof of C17 (CN variant positivity)

| Step | Statement | Status | Notes |
|:---|:---|:---|:---|
| 1 | Under E(9), \(P\) is a (nonsingular) M-matrix with \(P^{-1}\ge 0\) | VERIFIED | Diagonal dominance excess \(=1/\Delta t+r/2>0\) unconditionally. Off-diagonals \(\le 0\) from E(9). |
| 2 | Under E(9), \(N\ge 0\) | VERIFIED | Sub-diagonal is perfect square \(\ge 0\); super-diagonal always \(>0\); diagonal \(\ge 0\) from \(\Delta t\) bound. |
| 3 | \(U^{n+1}=P^{-1}NU^n\), and \(P^{-1}\ge 0\), \(N\ge 0\), \(U^0\ge 0\) | VERIFIED | Product of nonneg matrix and nonneg vector is nonneg. |
| 4 | Conclusion: \(U^{n+1}\ge 0\) for all \(n\) | VERIFIED | By induction. |
| **Verdict** | **Sound** | — | Note: paper claims \(P^{-1}>0\) (strict, requiring irreducibility), which can fail at one grid point (G7.4). For positivity of solution, \(P^{-1}\ge 0\) suffices. |

### Proof of C20 (Discrete maximum principle, CN variant)

| Step | Statement | Status | Notes |
|:---|:---|:---|:---|
| 1 | \(\|N\|_\infty = 1/\Delta t - r/2\) | VERIFIED | Row sum of \(N\) computed: all \(S_j\)-dependent terms cancel exactly. |
| 2 | \(\|P^{-1}\|_\infty \le (1/\Delta t+r/2)^{-1}\) | VERIFIED | Row sum of \(P\) is \(1/\Delta t + r/2\) (constant); Windisch bound. |
| 3 | \(\|U^{n+1}\|_\infty \le \|P^{-1}\|_\infty\|N\|_\infty\|U^n\|_\infty\) | VERIFIED | Sub-multiplicativity. Paper writes "=" (G4.4). |
| 4 | Bound \(\le \frac{1/\Delta t-r/2}{1/\Delta t+r/2}\|U^n\|_\infty\) | VERIFIED | Ratio \(<1\) since \(r>0\). |
| 5 | Conclusion: \(\|U^{n+1}\|_\infty \le \|U^n\|_\infty\) | VERIFIED | |
| **Verdict** | **Sound** | — | Cosmetic typo in step 3 (G4.4). Valid only under E(9). |

### Proof of C24 (Accuracy constraints for CN variant)

| Step | Statement | Status | Notes |
|:---|:---|:---|:---|
| 1 | Numerical diffusion is \(\tfrac{1}{8}(r\Delta S/\sigma)^2\) | VERIFIED | From C22/G1.4 (Taylor expansion I performed). |
| 2 | For accuracy, need this term to be small | VERIFIED | Qualitative. |
| 3 | From E(9), \(\Delta t \sim 8(\sigma/r)^2\) | VERIFIED | Dominant term in denominator is \(r^2/(8\sigma^2)\) for small \(\sigma\). |
| **Verdict** | **Sound** | — | |

---

## 4. IMPLEMENTATION BLOCKERS

Ordered by dependency (resolve upstream first):

| Priority | Gap ID | Description | Reason it blocks implementation |
|:---|:---|:---|:---|
| 1 | G2.2 + G2.3 + G6.3 | Grid construction: uniform grid \(S_j=j\Delta S\), relationship \(M=S_{\max}/\Delta S\), choice of \(S_{\max}\) | Cannot construct the computational domain or allocate arrays without knowing the grid. All matrix entries depend on \(S_j\). |
| 2 | G2.4 + G6.1 | Boundary treatment at \(S=0\) and \(S=S_{\max}\) | Cannot form the complete linear system. First/last rows of tridiagonal matrices require known boundary values; the modification is unspecified. For the double-barrier case, barrier conditions replace boundary conditions at monitoring dates — the projection E(13) is defined but the FD-level incorporation is not described. |
| 3 | G4.1 | Time-level indices in U7: \(n\) vs. \(n{+}1\) | For variable-coefficient implementations (\(\sigma(S,t)\), \(r(S,t)\)), using the wrong time level for coefficients gives a different scheme. Must be clarified before implementing the general case (though harmless for constant coefficients). |
| 4 | G7.3 | E(9) time-step constraint violated in examples | Implementer must decide: (a) enforce E(9), which for small \(\sigma\) and moderate \(r\) requires extremely small \(\Delta t\) (e.g. \(3.2\times10^{-5}\) for \(r=0.5\), \(\sigma=0.001\)); or (b) relax E(9) and accept that theoretical positivity/max-principle guarantees do not hold. The paper provides no guidance for case (b). |
| 5 | G4.3 | Denominators \(2h\) → \(h\) in U13/U14 | Blocks correct implementation of the pure-upwind limit (\(\sigma=0\)). Does not affect the exponentially fitted scheme (which uses E(3) directly), but would cause a factor-of-2 error if someone implements the upwind case from U13/U14 as written. |
| 6 | G3.3 + G7.1 | Convergence bound E(6) unproved | Does not block coding, but blocks verification: without a reference for the error bound, there is no way to validate convergence rates or set grid-refinement targets. |
| 7 | G6.2 | Variable-coefficient case unaddressed | Blocks implementation for stochastic volatility or local volatility models. The CN variant's parameter choice E(9) with \(\omega=-r/(16\sigma^2)\) depends on constant \(\sigma\); extension to \(\sigma(S,t)\) requires re-deriving the M-matrix conditions node-by-node. |
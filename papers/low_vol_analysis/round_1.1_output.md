
# Complete Mathematical Extraction

## Seed Paper
**Title:** Low Volatility Options and Numerical Diffusion of Finite Difference Schemes
**Authors:** Mariyan Milev; Aldo Tagliani
**Venue:** Serdica Mathematical Journal, **36** (2010), 223–236

---

## 1. NOTATION REGISTRY

| Symbol | Meaning | First appearance (§/Eq) | Type | Notes |
|:---|:---|:---|:---|:---|
| \(V(S,t)\) | Option price | §1, before Eq(1) | scalar function | |
| \(S\) | Underlying asset price | §1, before Eq(1) | scalar | |
| \(t\) | Time to expiry | §1, before Eq(1) | scalar | \(0 \le t \le T\) |
| \(T\) | Expiry (maturity) | §1, before Eq(1) | scalar | |
| \(r\) | Risk-free interest rate | §1, Eq(1) | scalar/function | Stated \(r = r(t,S)\) in §1; treated as constant in §§2–4 |
| \(\sigma\) (volatility) | Volatility parameter | §1, Eq(1) | scalar/function | Stated \(\sigma = \sigma(t,S)\) in §1; treated as constant in §§3–4. **OVERLOADED** — see next row |
| \(\sigma(S,t)\) (diffusion coeff.) | General diffusion coefficient in Eq(2) | §2, Eq(2) | scalar function | **OVERLOADED** with volatility \(\sigma\). Identified with \(\tfrac{1}{2}\sigma^2 S^2\) in U4 |
| \(\mu(S,t)\) | General convection/drift coefficient | §2, Eq(2) | scalar function | Identified with \(rS\) in U4 |
| \(b(S,t)\) | General reaction coefficient | §2, Eq(2) | scalar function | Identified with \(-r\) in U4 |
| \(L\) (operator) | Continuous differential operator | §2, after Eq(2) | differential operator | **OVERLOADED** with lower barrier \(L\) in §4 |
| \(L_k^h\) | Fitted finite-difference operator | §2, after \(L\) def | difference operator | |
| \(U_j^n\) | Discrete approximation to \(V(S_j,t_n)\) | §2, def of \(L_k^h\) | scalar | **OVERLOADED** — symbol \(U\) also used for upper barrier in §4 |
| \(h\) | Fixed spatial step | §2, def of \(L_k^h\) | scalar | Equivalent to \(\Delta S\) used in §§3–4; equivalence not stated |
| \(k\) | Fixed time step | §2, def of \(L_k^h\) | scalar | Equivalent to \(\Delta t\) used in §§3–4; equivalence not stated |
| \(\rho_j^{n+1}\) | Fitting factor | §2, Eq(3) | scalar | |
| \(\delta_x^2\) | Second central difference operator | §2, def of \(L_k^h\) | difference operator | **UNDEFINED** — standard meaning \(\delta_x^2 U_j = U_{j+1}-2U_j+U_{j-1}\) assumed |
| \(\mu_j^{n+1}\) | \(\mu\) evaluated at \((S_j,t_{n+1})\) | §2, def of \(L_k^h\) | scalar | |
| \(\sigma_j^{n+1}\) | \(\sigma(S,t)\) evaluated at \((S_j,t_{n+1})\) | §2, Eq(3) | scalar | This is the general diffusion coefficient, not volatility |
| \(b_j^{n+1}\) | \(b\) evaluated at \((S_j,t_{n+1})\) | §2, def of \(L_k^h\) | scalar | |
| \(A\) (§2 matrix) | Tridiagonal iteration matrix (exp.-fitted) | §2, Eq(4) | matrix | **OVERLOADED** — reused for standard implicit matrix in §4 |
| \(a_{i,j}\) | Entries of matrix \(A\) | §2, Eq(4) | scalar | |
| \(U^n\) | Solution vector at time level \(n\) | §2, Eq(4) | vector | |
| \(U^0\) | Initial solution vector | §2, positivity chain | vector | |
| \(c\) | Constant in convergence bound | §2, Eq(6) | scalar | Independent of \(h,k,\sigma\) |
| \(\Delta S\) | Spatial step | §3/p.225 top | scalar | Equivalent to \(h\) |
| \(\Delta t\) | Time step | §3/p.225 top | scalar | Equivalent to \(k\) |
| \(S_j\) | Spatial grid point | §3, matrices \(P,N\) | scalar | \(S_j = j\Delta S\) implied |
| \(\omega_1, \omega_2\) | Parameters in CN variant reaction-term discretization | §3, Eq(8) | scalar | |
| \(P\) | Left-hand-side tridiagonal matrix (CN variant) | §3, after Eq(8) | matrix | |
| \(N\) (matrix) | Right-hand-side tridiagonal matrix (CN variant) | §3, after Eq(8) | matrix | |
| \(M\) | Number of nodes in \(S\)-direction | §3, Eq(9) | integer | **UNDEFINED** prior to first use; no grid construction stated |
| \(K\) | Strike price | §4, Def 4.1 | scalar | |
| \(U\) (barrier) | Upper barrier | §4, Def 4.1 | scalar | **OVERLOADED** with discrete solution \(U_j^n\) |
| \(L\) (barrier) | Lower barrier | §4, Def 4.2 | scalar | **OVERLOADED** with differential operator \(L\) |
| \(f[S(T)]\) | Payoff function (truncated call) | §4, Def 4.1 | scalar function | |
| \(\mathbf{1}_{[L,U]}(S)\) | Indicator function of interval \([L,U]\) | §4, Eq(11)/Eq(13) | scalar function | Defined in text after Eq(13) |
| \(t_i\) | Monitoring dates | §4, Eq(13) | scalar | \(i=1,\ldots,F\) |
| \(F\) | Total number of monitoring intervals | §4, Eq(13) | integer | |
| \(S_{\max}\) | Domain truncation bound | §4, Fig 1 caption | scalar | **UNDEFINED** in main text; appears only in figure captions |
| \(\lambda_i(A^{-1})\) | Eigenvalues of \(A^{-1}\) | §4, p.232 | scalar (possibly complex) | |
| \(\tau_d\) | Characteristic diffusion time | §4, p.235 | scalar | Referenced to Tavella–Randall [8] |

### Overloaded symbols summary

| Symbol | Meaning 1 | Meaning 2 | Disambiguation |
|:---|:---|:---|:---|
| \(\sigma\) / \(\sigma(S,t)\) | Volatility parameter (§1, §3–4) | General diffusion coefficient in Eq(2) | Context: Eq(2) uses \(\sigma(S,t)\) for diffusion coeff.; elsewhere \(\sigma\) is volatility |
| \(L\) | Differential operator (§2) | Lower barrier (§4, Def 4.2) | Context-dependent |
| \(U\) | Discrete solution \(U_j^n\) (§2–3) | Upper barrier (§4, Defs 4.1–4.2) | Context-dependent |
| \(A\) | Exp.-fitted tridiag. matrix (§2) | Standard fully implicit tridiag. matrix (§4, p.232) | Section-dependent; different matrices |
| \(h,k\) vs \(\Delta S, \Delta t\) | Spatial/time steps in §2 | Spatial/time steps in §3–4 | Implicitly equivalent; never stated |

---

## 2. EQUATION CHAIN

### Numbered equations

| Eq ID | Expression | Depends on | Role |
|:---|:---|:---|:---|
| E(1) | \(-\dfrac{\partial V}{\partial t} + rS\dfrac{\partial V}{\partial S} + \dfrac{1}{2}\sigma^2 S^2\dfrac{\partial^2 V}{\partial S^2} - rV = 0\) | — | **Definition** (Black–Scholes PDE) |
| E(2) | \(-\dfrac{\partial V}{\partial t} + \mu(S,t)\dfrac{\partial V}{\partial S} + \sigma(S,t)\dfrac{\partial^2 V}{\partial S^2} + b(S,t)\,V = 0\) | E(1), U4 | **Definition** (generalized form of E(1)) |
| E(3) | \(\rho_j^{n+1} \equiv \dfrac{\mu_j^{n+1}\,h}{2}\coth\dfrac{\mu_j^{n+1}\,h}{2\,\sigma_j^{n+1}}\) | E(2), U6 | **Definition** (Duffy fitting factor) |
| E(4) | \(A\,U^{n+1} = U^n\) | E(3), U5, U6, U7 | **Discretization** (exp.-fitted finite difference equation) |
| E(5) | \(\|A^{-1}\|_\infty \;\le\; \dfrac{1}{1-kb} \;=\; \dfrac{1}{1+kr} \;<\; 1\) | E(4), U7, U8, U9 | **Derived result** (norm bound; references Windisch [9]) |
| E(6) | \(\lvert V(S_j,t_n) - U_j^n\rvert \;\le\; c(h+k)\) | E(4), E(3) | **Stated without proof** (convergence; \(c\) independent of \(h,k,\sigma\)) |
| E(7) | \(-\dfrac{\partial V}{\partial t} + \mu(S,t)\dfrac{\partial V}{\partial S} + \dfrac{1}{2}\mu(S,t)\,h\,\dfrac{\partial^2 V}{\partial S^2} - b(S,t)\,V = 0\) | U10, U13, U14, U15 | **Derived result** (modified equation solved by upwind scheme) |
| E(8) | \(V\!\left(t+\tfrac{\Delta t}{2}\right) = \omega_1(U_{j-1}^n+U_{j+1}^n)+\left(\tfrac{1}{2}-2\omega_1\right)U_j^n + \omega_2(U_{j-1}^{n+1}+U_{j+1}^{n+1})+\left(\tfrac{1}{2}-2\omega_2\right)U_j^{n+1}\) | E(1) | **Definition** (six-node reaction-term approximation for CN variant) |
| E(9) | \(\omega_1 = \omega_2 = -\dfrac{r}{16\sigma^2}\) and \(\Delta t < \dfrac{1}{r\!\left(\frac{1}{2}-2\omega_1\right)+\frac{1}{2}(\sigma M)^2}\) | E(8), U17–U20 | **Derived result** (parameter choice for positivity/M-matrix) |
| E(10) | \(-\dfrac{\partial V}{\partial t} + rS\dfrac{\partial V}{\partial S} + \dfrac{1}{8}\!\left(\dfrac{r}{\sigma}\,\Delta S\right)^{\!2}\dfrac{\partial^2 V}{\partial S^2} - rV = 0\) | E(8), E(9), U26 | **Derived result** (modified equation solved by CN variant when \(\sigma S\to 0\)) |
| E(11) | \(V(S,0) = (S-K)^+\,\mathbf{1}_{[L,U]}(S)\) | — | **Initial condition** (discrete double barrier knock-out call) |
| E(12) | \(V(S,t) \to 0 \;\text{as}\; S\to 0 \;\text{or}\; S\to\infty\) | — | **Boundary condition** |
| E(13) | \(V(S,t_i) = V(S,t_i^-)\,\mathbf{1}_{[L,U]}(S),\quad 0=t_0<t_1<\cdots<t_F=T\) | E(11) | **Boundary/monitoring condition** (knock-out projection at monitoring dates) |

### Unnumbered equations (in order of appearance)

| Eq ID | Expression | Location | Depends on | Role |
|:---|:---|:---|:---|:---|
| U1 | \(r = r(t,S)\) and \(\sigma = \sigma(t,S)\) | §1, after E(1) | E(1) | **Stated generality** (coefficient dependencies) |
| U2 | Exp.-fitted artificial diffusion: \(\tfrac{1}{2}rS\,\Delta S\,\dfrac{\partial^2 V}{\partial S^2}\) | p.225, top | — | **Stated without proof** (forward reference to §2 derivation) |
| U3 | CN variant artificial diffusion: \(\tfrac{1}{8}\!\left(\dfrac{r\,\Delta S}{\sigma}\right)^{\!2}\dfrac{\partial^2 V}{\partial S^2}\) | p.225, top | — | **Stated without proof** (forward reference to §3 derivation) |
| U4 | \(\sigma(S,t)=\tfrac{1}{2}\sigma^2 S^2,\;\mu(S,t)=rS,\;b(S,t)=-r\) | §2, after E(2) | E(1), E(2) | **Definition** (coefficient identification) |
| U5 | \(LV \equiv -\dfrac{\partial V}{\partial t}+\mu(S,t)\dfrac{\partial V}{\partial S}+\sigma(S,t)\dfrac{\partial^2 V}{\partial S^2}+b(S,t)\,V\) | §2 | E(2) | **Definition** (continuous operator) |
| U6 | \(L_k^h U_j^n \equiv -\dfrac{U_j^{n+1}-U_j^n}{k}+\mu_j^{n+1}\dfrac{U_{j+1}^{n+1}-U_{j-1}^{n+1}}{2h}+\rho_j^{n+1}\dfrac{\delta_x^2 U_j^{n+1}}{h^2}+b_j^{n+1}U_j^{n+1}\) | §2 | U5, E(3) | **Definition** (fitted difference operator) |
| U7 | \(A = [a_{i,j}]=\operatorname{tridiag}\!\left\{\!\left(-\dfrac{\rho_j^n}{h^2}+\dfrac{\mu_j^n}{2h}\right)\!k;\;\left(\dfrac{2\rho_j^n}{h^2}-b_j^n+\dfrac{1}{k}\right)\!k;\;-\!\left(\dfrac{\rho_j^n}{h^2}+\dfrac{\mu_j^n}{2h}\right)\!k\right\}\) | §2, p.226 | U6, E(3), E(4) | **Discretization** (iteration matrix). **Flag:** time-level index is \(n\) in matrix but \(n{+}1\) in U6 |
| U8 | \(a_{i,i+1}<0,\;a_{i+1,i}<0,\;a_{i,i}>0\) | §2, p.226 | U7 | **Derived result** (sign conditions on \(A\)) |
| U9 | \(A\) is an irreducible diagonally dominant tridiagonal M-matrix \(\Rightarrow A^{-1}>0\) | §2, p.226 | U7, U8 | **Derived result** (references Ortega [5]) |
| U10 | \(U^n = A^{-1}U^{n-1}=(A^{-1})(A^{-1}U^{n-2})=\cdots=(A^{-1})^n U^0 > 0\) | §2, p.226 | U9, E(4) | **Derived result** (positivity by induction) |
| U11 | \(\|U^{n+1}\|_\infty = \|A^{-1}U^n\|_\infty = \|A^{-1}\|_\infty\|U^n\|_\infty \le 1\cdot\|U^n\|_\infty \le \|U^n\|_\infty\) | §2, p.226 | E(5), U9 | **Derived result** (discrete max principle). **Flag:** second "\(=\)" should be "\(\le\)" |
| U12 | \(\displaystyle\lim_{\sigma\to 0}\rho = \lim_{\sigma\to 0}\frac{\mu h}{2}\coth\frac{\mu h}{2\sigma}=\begin{cases}\frac{\mu h}{2}&\mu>0\\[4pt]-\frac{\mu h}{2}&\mu<0\end{cases}\) | §2, p.227 | E(3) | **Derived result** (limit of fitting factor) |
| U13 | \(-\dfrac{U_j^{n+1}-U_j^n}{k}+\mu_j^{n+1}\dfrac{U_{j+1}^{n+1}-U_j^{n+1}}{2h}+b_j^{n+1}U_j^{n+1}=0,\;\mu>0\) | §2, p.227 | U12, U6 | **Derived result** (upwind scheme, \(\mu>0\)). **Flag:** denominator \(2h\) is inconsistent with derivation giving \(h\); likely typographical error |
| U14 | \(-\dfrac{U_j^{n+1}-U_j^n}{k}+\mu_j^{n+1}\dfrac{U_j^{n+1}-U_{j-1}^{n+1}}{2h}+b_j^{n+1}U_j^{n+1}=0,\;\mu<0\) | §2, p.227 | U12, U6 | **Derived result** (upwind scheme, \(\mu<0\)). **Same flag** as U13 |
| U15 | Numerical diffusion: \(\tfrac{1}{2}\mu(S,t)\,h\,\dfrac{\partial^2 V}{\partial S^2}\) | §2, p.227 | U13, U14 | **Derived result** (consistency analysis of upwind) |
| U16 | Hyperbolic target: \(-\dfrac{\partial V}{\partial t}+\mu(S,t)\dfrac{\partial V}{\partial S}-b(S,t)\,V=0\) | §2, p.227 | E(2) | **Definition** (pure convection limit of E(2) as \(\sigma\to 0\)) |
| U17 | \(P\,U^{n+1}=N\,U^n\) | §3, p.228 | E(1), E(8) | **Discretization** (CN variant linear system) |
| U18 | \(P=\operatorname{tridiag}\!\left\{r\omega_2+\dfrac{r}{4}\dfrac{S_j}{\Delta S}-\left(\dfrac{\sigma S_j}{2\,\Delta S}\right)^{\!2};\;\dfrac{1}{\Delta t}+\dfrac{1}{2}\!\left(\dfrac{\sigma S_j}{\Delta S}\right)^{\!2}+r\!\left(\dfrac{1}{2}-2\omega_2\right);\;r\omega_2-\dfrac{r}{4}\dfrac{S_j}{\Delta S}-\left(\dfrac{\sigma S_j}{2\,\Delta S}\right)^{\!2}\right\}\) | §3, p.228 | E(8), U17 | **Discretization** (left-hand matrix of CN variant) |
| U19 | \(N=\operatorname{tridiag}\!\left\{-r\omega_1-\dfrac{r}{4}\dfrac{S_j}{\Delta S}+\left(\dfrac{\sigma S_j}{2\,\Delta S}\right)^{\!2};\;\dfrac{1}{\Delta t}-\dfrac{1}{2}\!\left(\dfrac{\sigma S_j}{\Delta S}\right)^{\!2}-r\!\left(\dfrac{1}{2}-2\omega_1\right);\;-r\omega_1+\dfrac{r}{4}\dfrac{S_j}{\Delta S}+\left(\dfrac{\sigma S_j}{2\,\Delta S}\right)^{\!2}\right\}\) | §3, p.228 | E(8), U17 | **Discretization** (right-hand matrix of CN variant) |
| U20a | \(P\) is irreducibly diagonally dominant M-matrix \(\Rightarrow P^{-1}>0\) | §3, p.228 | U18, E(9) | **Design criterion** |
| U20b | \(N\) has nonnegative entries | §3, p.228 | U19, E(9) | **Design criterion** |
| U21 | \(U^{n+1}=P^{-1}N\,U^n=(P^{-1}N)^n U^0\) is positive since \(U^0\ge 0\) | §3, p.228 | U20a, U20b, U17 | **Derived result** (positivity preservation) |
| U22 | \(\|N\|_\infty = \dfrac{1}{\Delta t}-\dfrac{r}{2}\) | §3, p.228 | U19, E(9) | **Derived result** |
| U23 | \(\|P^{-1}\|_\infty \le \left(\dfrac{1}{\Delta t}+\dfrac{r}{2}\right)^{\!-1}\) | §3, p.228 | U18, E(9) | **Derived result** (references Windisch [9]) |
| U24 | \(\|U^{n+1}\|_\infty = \|(P^{-1}N)U^n\|_\infty = \|P^{-1}\|_\infty\|N\|_\infty\|U^n\|_\infty \le \dfrac{\frac{1}{\Delta t}-\frac{r}{2}}{\frac{1}{\Delta t}+\frac{r}{2}}\|U^n\|_\infty \le \|U^n\|_\infty\) | §3, p.228 | U22, U23 | **Derived result** (discrete max principle, CN variant). **Flag:** "\(=\)" in chain should be "\(\le\)" |
| U25 | Discretization error: \(O(\Delta S^2,\Delta t^2)\) | §3, p.228 | U17, U18, U19 | **Stated without proof** (accuracy order of CN variant) |
| U26 | Artificial diffusion from \(-rV\) discretization: \(\dfrac{1}{8}\!\left(\dfrac{r}{\sigma}\,\Delta S\right)^{\!2}\dfrac{\partial^2 V}{\partial S^2}\) | §3, p.229 top | E(8), E(9) | **Derived result** (stated as "standard analysis of consistency") |
| U27 | Requirement: \(\tfrac{1}{8}\!\left(\tfrac{r}{\sigma}\Delta S\right)^2\) must become insignificant | §3, p.229 | U26, E(10) | **Accuracy constraint** |
| U28 | From E(9): small \(\Delta t \sim 8\!\left(\dfrac{\sigma}{r}\right)^{\!2}\) | §3, p.229 | E(9), U27 | **Derived result** (time-step implied by accuracy constraint) |
| U29 | Comparison: \(\tfrac{1}{2}rS\,\Delta S\,\dfrac{\partial^2 V}{\partial S^2}\) vs \(\tfrac{1}{8}\!\left(\tfrac{r}{\sigma}\Delta S\right)^{\!2}\dfrac{\partial^2 V}{\partial S^2}\) | §3, p.229 | U2, U3 | **Derived result** (relative magnitudes of artificial diffusion) |
| U30 | \(f[S(T)]=\begin{cases}S(T)-K&S(T)\in[K,U]\\0&\text{otherwise}\end{cases}\) | §4, Def 4.1 | — | **Definition** (truncated call payoff) |
| U31 | \(\mathbf{1}_{[L,U]}(S)=\begin{cases}1&S\in[L,U]\\0&S\notin[L,U]\end{cases}\) | §4, after E(13) | — | **Definition** (indicator function) |
| U32 | \(A\,V^{n+1}=V^n\) | §4, p.232 | E(1) | **Discretization** (standard fully implicit, centered) |
| U33 | \(A=\operatorname{tridiag}\!\left\{-\dfrac{\Delta t}{2}\!\left[\!\left(\dfrac{\sigma S_j}{\Delta S}\right)^{\!2}\!-r\dfrac{S_j}{\Delta S}\right];\;1+\Delta t\!\left[\!\left(\dfrac{\sigma S_j}{\Delta S}\right)^{\!2}\!+r\right];\;-\dfrac{\Delta t}{2}\!\left[\!\left(\dfrac{\sigma S_j}{\Delta S}\right)^{\!2}\!+r\dfrac{S_j}{\Delta S}\right]\right\}\) | §4, p.232 | U32, E(1) | **Discretization** (standard fully implicit tridiagonal matrix with centered differences) |
| U34 | \(\tau_d = \dfrac{\Delta S^2}{(\sigma S)^2}\) | §4, p.235 | — | **Definition** (characteristic diffusion time; from Tavella–Randall [8]) |

---

## 3. THEOREM/CLAIM INVENTORY

| Claim ID | Statement | Proof provided? | Conditions required |
|:---|:---|:---|:---|
| C1 | The fitting factor \(\rho\) defined by E(3) is identically equal to 1 in the centred difference scheme (i.e., when \(\rho_j^{n+1}=\sigma_j^{n+1}\)). | No (stated as known fact) | — |
| C2 | With fitting factor E(3), the matrix \(A\) (U7) satisfies \(a_{i,i+1}<0\), \(a_{i+1,i}<0\), \(a_{i,i}>0\). | Partial (sign conditions argued from \(\coth(x)>1\) for \(x>0\)) | \(\mu_j > 0\) (holds for B-S since \(\mu=rS>0\) on interior grid) |
| C3 | \(A\) is an irreducible diagonally dominant tridiagonal M-matrix. | No (references Ortega [5], §6.2.3 p.104 and §6.2.17 p.110) | C2 (sign conditions) |
| C4 | \(A^{-1}>0\) (componentwise). | External, from [5] | C3 |
| C5 | \(\|A^{-1}\|_\infty \le \frac{1}{1-kb}=\frac{1}{1+kr}<1\). | External, from Windisch [9] | C3 and \(b=-r<0\) (i.e., \(r>0\)) |
| C6 | The numerical solution is positive: \((A^{-1})^n U^0>0\) for \(U^0>0\). | Yes (induction using C4) | C4, \(U^0>0\) |
| C7 | The scheme satisfies the discrete maximum principle: \(\|U^{n+1}\|_\infty \le \|U^n\|_\infty\). | Yes (norm bound argument using C5) | C5 |
| C8 | Convergence: \(\lvert V(S_j,t_n)-U_j^n\rvert\le c(h+k)\) where \(c\) is independent of \(h,k,\sigma\). | No (stated as known result) | Fitted scheme with E(3), stability and consistency of scheme. **Gap:** no proof or explicit reference given for this specific bound |
| C9 | The fitted scheme is uniformly stable for all values of \(h,k,\sigma\). | No (stated as consequence of preceding analysis) | C3, C5 |
| C10 | The fitted scheme is oscillation-free. | No (stated as consequence of M-matrix structure) | C3, C4 |
| C11 | In the limit \(\sigma\to 0\), the fitted scheme reduces to first-order implicit upwind schemes (U13, U14). | Yes (limit computation U12 substituted into U6) | E(3) |
| C12 | The upwind scheme introduces numerical diffusion \(\frac{1}{2}\mu(S,t)\,h\,\frac{\partial^2 V}{\partial S^2}\), so that it solves E(7) instead of U16. | Partially (stated as "through a standard analysis of consistency"; no Taylor expansion shown) | U13/U14 |
| C13 | The numerical diffusion of the upwind/exp.-fitted scheme is significant whenever \(r\) is large or \(h\) is not small enough. | No (qualitative observation) | C12 |
| C14 | Specializing C12 to B-S (\(\mu=rS\), \(h=\Delta S\)): artificial diffusion is \(\frac{1}{2}rS\,\Delta S\,\frac{\partial^2 V}{\partial S^2}\). | Follows from C12 and U4 | U4, C12 |
| C15 | The CN variant scheme (§3) is spurious-oscillation-free. | No (references Milev–Tagliani [4]) | E(9) |
| C16 | Under condition E(9), the matrices satisfy: \(P\) is irred. diag. dominant M-matrix with \(P^{-1}>0\); \(N\ge 0\). | No (references [4]; M-matrix theory invoked) | E(9) |
| C17 | The CN variant scheme is positivity-preserving under E(9). | Yes (combine \(P^{-1}>0\), \(N\ge 0\), and \(U^0\ge 0\)) | C16, \(U^0\ge 0\) |
| C18 | \(\|N\|_\infty = \frac{1}{\Delta t}-\frac{r}{2}\). | No (stated without derivation) | E(9) |
| C19 | \(\|P^{-1}\|_\infty \le \left(\frac{1}{\Delta t}+\frac{r}{2}\right)^{-1}\). | External, from Windisch [9] | C16 |
| C20 | The CN variant satisfies the discrete maximum principle: \(\|U^{n+1}\|_\infty\le\|U^n\|_\infty\). | Yes (norm product argument) | C18, C19 |
| C21 | The CN variant has discretization error \(O(\Delta S^2,\Delta t^2)\). | No (stated without derivation) | Standard Crank–Nicolson accuracy + E(8) modification |
| C22 | When \(\sigma\) is small, the artificial diffusion of the CN variant from the \(-rV\) discretization amounts to \(\frac{1}{8}\left(\frac{r}{\sigma}\Delta S\right)^2\frac{\partial^2 V}{\partial S^2}\). | Partially ("from a standard analysis of consistency"; no Taylor expansion shown) | E(8), E(9) |
| C23 | When \(\sigma S\to 0\), the CN variant effectively solves E(10). | Follows from C22 (dominant diffusion term) | C22 |
| C24 | An accurate solution from the CN variant requires \(\frac{1}{8}\left(\frac{r}{\sigma}\Delta S\right)^2\) to be insignificant, implying very small \(\Delta S\) and from E(9) small \(\Delta t\sim 8(\sigma/r)^2\). | Yes (algebraic consequence of C22 and E(9)) | C22, E(9) |
| C25 | Standard Crank–Nicolson produces spurious oscillations for every choice of \(\Delta S\) and \(\Delta t\) when \(\sigma^2\ll r\) with discontinuous payoff. | Empirical (demonstrated in Fig 1) | \(\sigma^2\ll r\), discontinuous initial data |
| C26 | Standard fully implicit (centered) also produces spurious oscillations when \(\sigma^2\ll r\). | Empirical (demonstrated in Fig 2) | Same as C25 |
| C27 | For the standard fully implicit centered scheme (U33): if \(\sigma^2>r\) is violated, positivity is not guaranteed; some \(\lambda_i(A^{-1})\) may become complex. | No (stated without proof; references spectral analysis) | \(\sigma^2 < r\) |
| C28 | The numerical diffusion of the exp.-fitted scheme depends on \(r\) and \(S\); higher \(r\) gives more diffusion. | Empirical (demonstrated in Fig 3) | C14 |
| C29 | The numerical diffusion of the CN variant depends on \(r\), \(\sigma\), and \(S\); higher \(r\) gives more diffusion. | Empirical (demonstrated in Fig 4) | C22 |
| C30 | For low volatility, the CN variant is oscillation-free and positive in both \(r=0.01\) and \(r=0.5\) cases. | Empirical (Fig 4) | E(9) |
| C31 | Numerical diffusion of both schemes can be diminished by choosing smaller \(\Delta S\). | Empirical (Fig 5; \(\Delta S=0.025\) vs \(\Delta S=0.05\)) | C14, C22 |
| C32 | With \(\Delta S=0.01\), \(\Delta t=0.001\), the solutions of Duffy's scheme and the CN variant are "practically indistinguishable and much more accurate." | Empirical (stated p.235) | C14, C22 |
| C33 | An optimal finite difference scheme does not exist because the numerical diffusion depends on different parameters for each scheme. | Qualitative conclusion (§§4–5) | — |
| C34 | Whenever \(\Delta t\gg\tau_d\), an oscillating behavior may arise. | External, from Tavella–Randall [8] | U34 |
| C35 | Giles–Carter (2006) demonstrated that for non-smooth initial data in B-S, convergence holds in \(L^2\) norm but not in supremum norm. | External, from [3] | Non-smooth payoff |

---

## 4. ASSUMPTION INVENTORY

| Assumption ID | Statement | Explicit / Implicit | Used in (Eq/Claim refs) |
|:---|:---|:---|:---|
| A1 | The underlying follows a geometric Brownian motion under the risk-neutral measure, yielding the B-S PDE E(1). | **Explicit** (§1, "the price \(V(S,t)\) of the option satisfies (1)") | E(1), all subsequent |
| A2 | \(r\) and \(\sigma\) are stated as functions \(r(t,S)\), \(\sigma(t,S)\) in §1, but are treated as **constants** throughout §§2–4 (all discretizations and parameter choices use constant \(r,\sigma\)). | **Implicit** (generality stated but never used; all formulas use constant \(r,\sigma\)) | E(9), U4, U18, U19, U33, all numerical examples |
| A3 | The computational domain is truncated to \([0, S_{\max}]\). | **Implicit** (\(S_{\max}\) appears in figure captions but no truncation is formally introduced in the text; \(M\) nodes are referenced without defining \(S_{\max}=M\Delta S\)). | U7, U18, U19, U33, E(9), all examples |
| A4 | A uniform spatial grid \(S_j = j\,\Delta S\), \(j=0,\ldots,M\), is used. | **Implicit** (never formally stated; inferred from the tridiagonal structure and use of \(S_j\) in matrix entries) | U7, U18, U19, U33 |
| A5 | A uniform time grid \(t_n = n\,\Delta t\) is used. | **Implicit** (never formally stated; inferred from the difference operators) | U6, U7, U17 |
| A6 | Boundary conditions at \(S=0\) and \(S=S_{\max}\) are compatible with the tridiagonal system (interior-node formulation). | **Implicit** (no boundary treatment is described; the tridiagonal systems are written for interior nodes only) | E(4), U17, U32 |
| A7 | \(r > 0\). | **Implicit** (used in E(5) to get \(\frac{1}{1+kr}<1\), and in sign conditions; never stated) | E(5), C5, C7, C20 |
| A8 | \(\sigma > 0\) (volatility is positive). | **Implicit** (the fitting factor E(3) requires \(\sigma_j^{n+1} > 0\) in the denominator; the limit \(\sigma\to 0\) is treated separately) | E(3), U12 |
| A9 | \(\mu_j > 0\) on all interior grid points, i.e., \(rS_j > 0\) for \(j\ge 1\). | **Implicit** (required for sign condition \(a_{i+1,i}<0\) in U8, since \(\coth(x)>1\) for \(x>0\) needs \(\mu>0\)) | C2, U8 |
| A10 | \(U^0 \ge 0\) (nonnegative initial data). | **Explicit** (stated in context of positivity proofs: "since \(U^0\ge 0\)") | C6, C17 |
| A11 | The fitting factor \(\rho\) is identically 1 for the centered difference scheme. | **Explicit** (stated p.226, top) | C1 |
| A12 | For Eq(5)/(C5): \(b = -r\) (specialization to Black–Scholes). | **Explicit** (used in simplification \(1-kb = 1+kr\)) | E(5) |
| A13 | For E(9): \(\omega_1 = \omega_2\) (equal weighting on old/new time levels). | **Explicit** (stated in E(9)) | E(9), U18, U19 |
| A14 | For CN variant analysis: \(r\) and \(\sigma\) are constant (the matrix entries in U18/U19 use constant \(r,\sigma\)). | **Implicit** (the general coefficients from E(2) are specialized) | U18, U19, E(9) |
| A15 | \(S_{\max}\) is "sufficiently large" so that boundary effects at \(S_{\max}\) are negligible. | **Implicit** (standard in FD option pricing but never stated; only \(V\to 0\) as \(S\to\infty\) is given) | A3, E(12) |
| A16 | The exponentially fitted scheme's convergence bound E(6) holds with constant \(c\) independent of \(h,k,\sigma\). | **Stated without proof** (no reference to a specific proof; described as "the following result") | C8, E(6) |
| A17 | The "standard analysis of consistency" (modified equation analysis via Taylor expansion) applies to derive U15, C12, C22, U26. | **Implicit** (the paper invokes this analysis twice but never performs it) | C12, C22 |
| A18 | The discrete maximum principle (C7, C20) is stated in the sup-norm sense \(\|U^{n+1}\|_\infty \le \|U^n\|_\infty\). The paper equates this with "the scheme satisfies the discrete maximum principle." | **Explicit** | C7, C20 |
| A19 | For Section 4 examples: discontinuities in the initial/monitoring conditions are the primary source of spurious oscillations (combined with low volatility). | **Explicit** (stated in §§1, 4) | C25, C26 |
| A20 | The standard second central difference operator is used for \(\partial^2 V/\partial S^2\) in all schemes (exponentially fitted, CN variant, standard implicit). | **Implicit** (stated only through the discrete operators; never explicitly declared as a design choice) | U6, U7, U18, U19, U33 |
| A21 | All schemes use a fully implicit treatment of the spatial terms (except CN variant which uses Crank–Nicolson averaging for diffusion/convection). | **Implicit** (the time-level structure is visible in the operators but the implicit/CN characterization is implicit in the matrix structure) | U6, U17 |

---

## 5. LOGICAL DEPENDENCY GRAPH

```
────────────────────────────────────
FOUNDATIONS
────────────────────────────────────

A1 (risk-neutral GBM)
 └─→ E(1) [Black–Scholes PDE]
      ├─→ U1 [r(t,S), σ(t,S) dependencies]
      └─→ E(2) [generalized form]
           └─→ U4 [coefficient identification: σ(S,t)=½σ²S², μ=rS, b=−r]
                └─→ [connects E(2) back to E(1)]

────────────────────────────────────
SECTION 2: EXPONENTIALLY FITTED SCHEME
────────────────────────────────────

E(2) ──→ U5 [operator L definition]
          └─→ U6 [fitted operator L_k^h]
               │
               ├─→ E(3) [fitting factor ρ]
               │    ├── A8 (σ > 0)
               │    └── A9 (μ > 0)
               │
               └─→ E(4) [FD equation A U^{n+1} = U^n]
                    │
                    ├─→ U7 [tridiagonal matrix A]
                    │    │   ├── A3 (domain truncation)
                    │    │   ├── A4 (uniform spatial grid)
                    │    │   └── A5 (uniform time grid)
                    │    │
                    │    └─→ U8 [sign conditions on A]
                    │         │   └── A9 (μ > 0)
                    │         │       └── E(3) [ρ > |μh/2| when μ > 0]
                    │         │
                    │         └─→ C2 [sign conditions hold]
                    │
                    ├── C2 ──→ C3 [A is M-matrix]
                    │           │   └── External: Ortega [5]
                    │           │
                    │           └─→ C4 [A⁻¹ > 0]
                    │
                    ├── A7 (r > 0), A12 (b = −r)
                    │    └── C4 ──→ E(5) [‖A⁻¹‖∞ ≤ 1/(1+kr) < 1]
                    │               │   └── External: Windisch [9]
                    │               │
                    │               ├─→ C5 [norm bound]
                    │               │
                    │               └─→ U11 [discrete max principle]
                    │                    └─→ C7 [max principle holds]
                    │
                    ├── C4, A10 (U⁰ ≥ 0)
                    │    └─→ U10 [positivity induction chain]
                    │         └─→ C6 [solution is positive]
                    │
                    └── C3, C5 ──→ C9 [uniformly stable]
                                   C10 [oscillation-free]

E(3) ──→ U12 [limit σ → 0]
          ├─→ U13 [upwind scheme μ > 0]
          └─→ U14 [upwind scheme μ < 0]
               │
               └─→ C11 [fitted → upwind in limit]

U13, U14 ──→ U15 [numerical diffusion ½μh ∂²V/∂S²]
              │    └── A17 (consistency analysis)
              │
              └─→ C12 [upwind introduces diffusion]
                   │
                   ├─→ U16 [hyperbolic target equation]
                   ├─→ E(7) [modified parabolic equation]
                   │
                   └── U4 ──→ C14 [B-S specialization: ½rSΔS ∂²V/∂S²]
                              └─→ U2 [forward ref resolved]

C12, C14 ──→ C13 [diffusion significant for large r or large h]

E(6) ──→ C8 [convergence] ← A16 (stated without proof)

────────────────────────────────────
SECTION 3: CRANK–NICOLSON VARIANT
────────────────────────────────────

E(1) ──→ E(8) [six-node reaction-term discretization]
          │    └── A14 (constant r, σ)
          │
          └─→ U17 [P U^{n+1} = N U^n]
               ├─→ U18 [matrix P]
               └─→ U19 [matrix N]
                    │
                    ├── U20a, U20b [design criteria: P is M-matrix, N ≥ 0]
                    │
                    └─→ E(9) [ω₁ = ω₂ = −r/(16σ²); Δt bound]
                         │
                         ├─→ C16 [M-matrix / nonnegativity verified]
                         │    │   └── External: [4], [5], [9]
                         │    │
                         │    ├─→ C17 [positivity-preserving]
                         │    │    └── A10 (U⁰ ≥ 0)
                         │    │         └─→ U21
                         │    │
                         │    ├─→ C18 [‖N‖∞ = 1/Δt − r/2]
                         │    ├─→ C19 [‖P⁻¹‖∞ bound]
                         │    │    └── External: Windisch [9]
                         │    │
                         │    └── C18, C19 ──→ U24 [max principle inequality]
                         │                      └─→ C20 [discrete max principle]
                         │
                         └─→ C15 [oscillation-free] ← External: [4]

U17, U18, U19 ──→ C21 [O(ΔS², Δt²) error] ← U25

E(8), E(9), A17 ──→ U26 [artificial diffusion ⅛(r/σ ΔS)² ∂²V/∂S²]
                     │
                     └─→ C22 [CN variant artificial diffusion]
                          │
                          ├─→ E(10) [modified PDE for CN variant]
                          │    └─→ C23 [effectively solves E(10) when σS → 0]
                          │
                          ├─→ U27, U28 [accuracy constraints]
                          │    └─→ C24 [needs very small ΔS, Δt ~ 8(σ/r)²]
                          │
                          └─→ U3 [forward ref resolved]

C14, C22 ──→ U29 [comparison of diffusion terms]

────────────────────────────────────
SECTION 4: NUMERICAL RESULTS
────────────────────────────────────

E(1), A3, A4, A5 ──→ U32, U33 [standard fully implicit centered scheme]
                      │
                      └─→ C27 [σ² > r violated ⟹ loss of positivity,
                               possible complex eigenvalues]

E(1), E(11), E(12), E(13) ──→ [problem specifications for examples]

U30 [truncated call payoff] ← Definition 4.1
U31 [indicator function] ← Definition 4.2, after E(13)

C25, C26 ──→ [empirical: CN and fully implicit oscillate
              when σ² ≪ r, discontinuous payoff]
              └── A19 (discontinuities + low vol cause oscillations)

C14 ──→ C28 [empirical: higher r ⟹ more Duffy diffusion]
C22 ──→ C29 [empirical: higher r ⟹ more CN variant diffusion]
C22, E(9) ──→ C30 [CN variant oscillation-free for both r values]

C14, C22 ──→ C31 [smaller ΔS diminishes diffusion]
              └─→ C32 [ΔS=0.01, Δt=0.001 ⟹ schemes indistinguishable]

C14, C22, C33 ──→ [no optimal scheme exists]

U34 [τ_d] ← External: Tavella–Randall [8]
      └─→ C34 [Δt ≫ τ_d ⟹ oscillation risk]

External: Giles–Carter [3] ──→ C35 [L² convergence but not sup-norm
                                     for non-smooth data]
```

### Key gaps and unresolved items

| Item | Nature | Location |
|:---|:---|:---|
| E(6) / C8 | Convergence bound stated without proof or explicit reference; constant \(c\) independence from \(\sigma\) is the key claim enabling "uniform convergence" | §2, p.226 |
| U13/U14 denominator | Paper writes \(2h\); mathematical derivation from U12+U6 yields \(h\). Likely typographical error; the numerical diffusion formula U15 is consistent with \(h\) | §2, p.227 |
| U7 time-level index | Matrix \(A\) uses superscript \(n\) on coefficients; fitted operator U6 evaluates at \(n{+}1\). Inconsistency not addressed | §2, p.226 |
| U11/U24 equality chains | Both use "\(=\)" where "\(\le\)" is mathematically required in the submultiplicativity step | §2 p.226, §3 p.228 |
| C12 / C22 | Artificial diffusion terms derived via "standard analysis of consistency" — no Taylor expansion is shown | §2 p.227, §3 p.229 |
| C21 / U25 | \(O(\Delta S^2, \Delta t^2)\) accuracy of CN variant stated without derivation | §3, p.228 |
| A3 / A15 | Domain truncation \(S_{\max}\) and grid construction \(M, \Delta S\) never formally introduced in the body text | Throughout |
| A6 | Boundary condition treatment in the tridiagonal systems is never specified | Throughout |

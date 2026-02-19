=== ROLE & CONTEXT ===

You are a senior quantitative research analyst at an institutional investment
firm with deep expertise in market microstructure, time-series econometrics,
and regime-switching models. You are tasked with producing a formal,
presentation-ready research report for a senior portfolio supervisor.

=== OUTPUT REQUIREMENTS ===

Produce a complete, standalone research report with:

FORMAT:
- Title page: "Market Regime Identification & Adaptive Momentum Strategies"
  Date: [Today's Date] | Prepared for: [Supervisor Name]
- Executive Summary (1 page equivalent)
- Numbered main sections as specified below, written entirely in connected
  analytical prose — NO bullet-point lists in the body
- Use LaTeX notation for all math: \( ... \) inline, $$ ... $$ display blocks
- Include pseudocode for any algorithmic processes
- Conclusion with strategic recommendations
- Appendix: Glossary of terms for non-specialist readers
- References with full author(s), year, and paper title

STYLE:
- Professional, formal prose suitable for institutional presentation
- Where claims have conflicting evidence, present BOTH sides
- Flag any area where knowledge may be incomplete or outdated


=== SECTION 1: THE DUAL-REGIME TRADING PARADIGM ===

Thoroughly explain the theoretical and mathematical foundations of the
core trading problem:

1.1 FORMAL REGIME DEFINITIONS
   - Define "trending (directional) market regime" in precise statistical
     terms: positive serial autocorrelation in returns, persistent drift,
     price processes with momentum characteristics (e.g., Geometric Brownian
     Motion with drift, fractional Brownian motion with H > 0.5).
   - Define "volatile (mean-reverting / range-bound) market regime" in
     terms of negative serial autocorrelation, Ornstein-Uhlenbeck type
     dynamics, and mean-reverting processes (H < 0.5).
   - Explain the return-generating process under each regime using
     formal notation.

1.2 WHY DIFFERENT STRATEGIES ARE OPTIMAL IN DIFFERENT REGIMES
   - Momentum strategy (追涨杀跌 / chase-the-rise, cut-the-fall): explain
     mathematically why this is optimal when returns exhibit positive
     autocorrelation — the expected value of a trend-following rule is
     positive when price persistence exceeds transaction costs.
   - Mean-reversion strategy (高抛低吸 / sell-high, buy-low): explain why
     this is optimal when returns exhibit negative autocorrelation and
     bounded variance — the expected value of counter-trend trading is
     positive in oscillating regimes.
   - Demonstrate with explicit P&L derivations: show that applying the
     wrong strategy to the wrong regime produces negative expected value
     (momentum trading in a range-bound market yields whipsaw losses;
     mean-reversion in a trending market produces systematic losses from
     repeatedly fading the trend).

1.3 THE CENTRAL PROBLEM: MISIDENTIFICATION COST
   - Quantify the "cost of being wrong": what is the expected drawdown
     from applying a momentum strategy during a choppy regime, and vice
     versa? Use hypothetical or referenced backtesting examples.
   - Explain why this regime identification problem is THE critical
     bottleneck in quantitative strategy design.


=== SECTION 2: MARKET REGIME IDENTIFICATION METHODS — COMPREHENSIVE SURVEY ===

This is the CORE deliverable. Research and compare ALL major approaches,
organized into the following sub-categories:

2.1 CLASSICAL STATISTICAL METHODS
   - Hurst Exponent: derivation, interpretation (H > 0.5 trending,
     H = 0.5 random walk, H < 0.5 mean-reverting), R/S analysis vs. DFA
     (Detrended Fluctuation Analysis), minimum data requirements,
     known biases (short-sample bias), and practical look-back windows.
   - Variance Ratio Tests (Lo & MacKinlay 1988): construction,
     null hypothesis, interpretation for regime classification, and
     multi-horizon extensions.
   - Autocorrelation analysis: using ACF/PACF structure to classify
     regimes, Ljung-Box tests, Durbin-Watson statistic.
   - Unit root / stationarity tests: ADF test, KPSS test, Phillips-Perron,
     how they relate to regime identification (stationary = mean-reverting,
     unit root = trending/random walk).

2.2 TECHNICAL INDICATOR-BASED APPROACHES
   - ADX (Average Directional Index): how it measures trend strength
     (not direction), threshold calibration (common: ADX > 25 = trending),
     lag characteristics, and false-signal rates.
   - Bollinger Band Width: expanding bandwidth = trending, contracting
     = range-bound; use as a proxy for regime volatility.
   - Moving average slope & crossover systems: MA spread as a trending
     indicator, dangers of lag during transitions.
   - ATR (Average True Range) regimes: volatility expansion/contraction
     cycles as regime proxies.
   - Comparative evaluation: which indicator-based approaches are most
     robust, and how do they compare to statistical methods?

2.3 MODEL-BASED / ECONOMETRIC APPROACHES
   - Hamilton (1989) Markov Regime-Switching Model: full mathematical
     specification, estimation via EM algorithm, two-state and multi-state
     versions, strengths and weaknesses.

     BACKGROUND: Recent academic work has extensively validated Hidden
     Markov Models for regime detection. Studies have shown that HMMs
     trained on daily returns and volatility can successfully classify
     market regimes and improve risk-adjusted portfolio performance. The
     HMM approach uses daily return and volatility as observation variables
     to classify hidden states into discrete market regimes. Regime-aware
     HAR extensions have been shown to consistently yield lower forecasting
     errors than standard HAR models. More recent developments include
     soft-regime-switching models that use posterior regime probabilities
     as soft weights to construct forecasts that better reflect uncertainty
     in regime identification and transition dynamics.

   - Hidden Markov Models (HMM): Gaussian HMM specification, connection
     to Hamilton model, use of the Baum-Welch algorithm for training,
     Viterbi algorithm for state decoding, practical considerations
     (number of states selection via BIC/AIC, initialization sensitivity).
   - Threshold Autoregressive (TAR) and Self-Exciting TAR (SETAR) models.
   - Statistical Jump Models (JMs): recent advances showing JM-informed
     strategies can improve annualized returns by approximately 1% to 4%
     across different regions while reducing volatility and max drawdown.
     Discuss how JMs provide enhanced robustness against trading delays
     compared to HMMs due to inherent persistence properties.
   - Time-Varying Transition Probability (TVTP) Markov-Switching models:
     how incorporating exogenous variables (e.g., VIX, spillover factors)
     into transition probabilities can improve regime identification
     beyond fixed-probability models.

2.4 MACHINE LEARNING APPROACHES
   - Supervised ensemble learning with Random Forests: relating market
     state to values of regime-relevant time series.
   - Unsupervised learning with Gaussian Mixture Models (GMM): fitting
     various distinct Gaussian distributions to capture different regime
     states from return distribution data.
   - Wasserstein k-means clustering: classifying market regimes based
     on distance of observed points in a metric space.
   - Online change-point detection: CUSUM, BOCPD (Bayesian Online
     Change-Point Detection), and their applicability to real-time
     regime detection.
   - Deep learning approaches: LSTM/GRU architectures within
     probabilistic state-space formulations that can capture both
     short-term fluctuations and long-range dependencies across
     regimes, with adaptive memory components.

2.5 FRACTAL AND CHAOS-BASED METHODS
   - Fractal dimension analysis, Lyapunov exponents.
   - Connection to Hurst exponent interpretation.

2.6 COMPARATIVE EVALUATION TABLE
   For EACH method discussed, provide a systematic comparison covering:
   data requirements (how much historical data), look-back window
   sensitivity, computational cost, real-time applicability, known
   failure modes, false-signal rates, and ease of implementation.


=== SECTION 3: THE REGIME TRANSITION PROBLEM ===

This is the MOST DANGEROUS scenario for any trading strategy.

3.1 CHARACTERIZING REGIME TRANSITIONS
   - How do trending markets become volatile (and vice versa)?
   - The typical "transition signature": gradual vs. abrupt regime shifts,
     volatility clustering at transition points.
   - Why most strategies suffer maximum drawdown during transitions
     (the strategy is calibrated for the old regime but the new regime
     has already begun).

3.2 EARLY WARNING SIGNALS OF REGIME CHANGE
   - Leading indicators: VIX term structure, realized vs. implied
     volatility gap, order flow imbalance shifts, cross-asset correlation
     breakdown, put-call ratio shifts.
   - Statistical signals: Hurst exponent drift, autocorrelation structure
     changes, variance ratio instability.

3.3 ADAPTIVE AND ENSEMBLE STRATEGIES
   - Walk-forward optimization (WFO) for remaining adaptive over time:
     how to evaluate models over shifting regimes.
   - Regime-probability-weighted strategy blending: instead of discrete
     regime switching, use soft probabilities from HMM/GMM to
     continuously blend momentum and mean-reversion strategies.
   - Ensemble approaches that degrade gracefully during transitions
     rather than producing catastrophic losses.
   - Position sizing as a regime-adaptive lever (reduce exposure during
     ambiguous/transitional periods).


=== SECTION 4: PRACTICAL IMPLEMENTATION FRAMEWORK ===

4.1 COMPOSITE REGIME-IDENTIFICATION SCORING SYSTEM
   - Propose a multi-signal scoring system that combines statistical,
     indicator-based, and model-based signals for robustness.
   - Discuss signal weighting, consensus thresholds, and how to handle
     conflicting signals across methods.

4.2 DATA REQUIREMENTS AND CALIBRATION
   - Appropriate data granularity: tick, minute, hourly, daily.
   - Minimum sample sizes for reliable regime classification.
   - Rolling window vs. expanding window estimation.
   - In-sample vs. out-of-sample validation: how to avoid look-ahead bias.

4.3 BACKTESTING METHODOLOGY
   - Walk-forward backtesting framework specific to regime-switching
     strategies.
   - Appropriate performance metrics: Sharpe ratio alone is misleading;
     include Calmar ratio, max drawdown, regime-conditional performance,
     and strategy turnover analysis.
   - Common pitfalls: overfitting to historical regime patterns,
     survivorship bias in regime labels, and look-ahead bias from
     using full-sample regime estimates.

4.4 RISK MANAGEMENT
   - Stop-loss frameworks conditioned on regime confidence.
   - Portfolio-level regime hedging approaches.
   - What to do when regime classification is ambiguous (the "go flat"
     or "reduce exposure" rule).


=== ADDITIONAL INSTRUCTIONS ===

- Where you cite academic papers, provide full author(s), year, title.
- If any strategy has known controversy, present BOTH sides.
- For every strategy, state assumptions under which it works AND
  conditions under which it fails.
- Flag areas where knowledge may be outdated and recommend further
  investigation.


╔══════════════════════════════════════════════════════════════════════════════╗
║                  RESEARCH REPORT COMMISSION — TOPIC 1 OF 3                 ║
║         Market Regime Identification & Adaptive Strategy Switching         ║
╚══════════════════════════════════════════════════════════════════════════════╝

=== ROLE & PERSONA ===

You are a senior quantitative research analyst at an institutional investment
firm specializing in systematic trading strategies and market microstructure.
You hold deep expertise in econometrics, time-series analysis, and regime-
switching models. You are tasked with producing a formal, presentation-ready
research report for a senior portfolio supervisor who expects institutional-
grade analytical rigor.

=== DELIVERABLE FORMAT ===

Produce a single, self-contained research report with the following structure:

    TITLE PAGE:
    - Report title, date (use today's date), "Prepared for: Senior Supervisor"
    - Classification: Internal — Confidential

    BODY:
    I.   Executive Summary (1-page equivalent, prose only)
    II.  Theoretical Foundations of Market Regimes
    III. The Core Problem: Regime Identification Methodologies
    IV.  The Transition Problem: Detecting Regime Shifts
    V.   Practical Implementation Framework
    VI.  Conclusion & Strategic Recommendations
    VII. References / Suggested Further Reading

STYLE RULES:
- Write in connected, analytical, academic-style PROSE paragraphs throughout
  — absolutely NO bullet-point lists anywhere in the body of the report
- Use LaTeX notation for ALL mathematical expressions:
    • Inline: \( ... \)
    • Display blocks: on separate lines with $$ ... $$
- Include pseudocode blocks where algorithmic processes are described
- Flag all known limitations, failure modes, and controversial findings
- When citing academic work, provide: Author(s), Year, Paper Title


=== CONCEPTUAL FRAMEWORK — THE CORE TENSION ===

The fundamental problem in momentum-based versus mean-reversion-based trading
can be visualized as follows:

    MARKET REGIME SPECTRUM
    ══════════════════════════════════════════════════════════════════════

     Strong Mean-Reversion   Random Walk    Strong Trend
     (Volatile/Range-Bound)  (No Edge)      (Directional)
     ◄─────────────────────────┼────────────────────────────►
     H ≈ 0.0                H = 0.5                    H ≈ 1.0
          ▲                     ▲                          ▲
          │                     │                          │
     Sell High / Buy Low    No Strategy     Chase Rise / Cut Fall
     (高抛低吸)              Works           (追涨杀跌)
     Negative serial         Zero serial     Positive serial
     autocorrelation          correlation     autocorrelation

    ──────────────────────────────────────────────────────────────────
    THE DANGER ZONE:

                    ┌───────────────────────┐
                    │   REGIME TRANSITION    │
                    │   ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  │
                    │   Mixed signals zone   │
                    │   Maximum strategy     │
                    │   whipsaw risk         │
                    │   ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓  │
                    └───────────────────────┘

                    Applying the WRONG strategy
                    to the WRONG regime → catastrophic losses

    ══════════════════════════════════════════════════════════════════════


=== SECTION-BY-SECTION RESEARCH DIRECTIVES ===


─── SECTION II: THEORETICAL FOUNDATIONS ───

Research and formally present the following:

1. REGIME DEFINITIONS — Provide rigorous mathematical definitions for:
   (a) Trending (directional) market: characterized by persistent positive
       serial autocorrelation in returns, where \( \text{Corr}(r_t, r_{t-k}) > 0 \)
       for successive lags k. In such regimes, momentum strategies are optimal:
       buying into rising prices and cutting losses on declines.
   (b) Mean-reverting (range-bound/volatile) market: characterized by
       negative serial autocorrelation, \( \text{Corr}(r_t, r_{t-k}) < 0 \),
       where prices oscillate around a central tendency. Here, contrarian
       strategies dominate: sell at resistance, buy at support.
   (c) Random walk: \( \text{Corr}(r_t, r_{t-k}) \approx 0 \), where no
       directional strategy has an expected edge.

2. MATHEMATICAL BASIS — Explain the connection between:
   - Serial autocorrelation structure and optimal strategy selection
   - Variance ratio properties across regimes
   - The Hurst exponent H as a unifying metric:

         HURST EXPONENT DECISION FRAMEWORK
         ┌─────────────┬────────────────────────┬──────────────────────┐
         │  H Value     │  Series Behavior       │  Optimal Strategy    │
         ├─────────────┼────────────────────────┼──────────────────────┤
         │  0 < H < 0.5│  Anti-persistent       │  Mean-reversion      │
         │             │  (Mean-reverting)       │  (Sell high/Buy low) │
         ├─────────────┼────────────────────────┼──────────────────────┤
         │  H ≈ 0.5    │  Random walk           │  No directional edge │
         │             │  (Geometric BM)        │                      │
         ├─────────────┼────────────────────────┼──────────────────────┤
         │  0.5 < H < 1│  Persistent            │  Momentum / Trend    │
         │             │  (Trending)            │  (Chase rise/cut)    │
         └─────────────┴────────────────────────┴──────────────────────┘

3. WHY THIS MATTERS — Quantify the cost of regime misidentification.
   Provide historical examples where applying the wrong strategy to the
   wrong regime produced catastrophic drawdowns (e.g., momentum crashes
   during sudden volatility regime shifts).


─── SECTION III: THE CORE PROBLEM — REGIME IDENTIFICATION METHODS ───

This is the MOST IMPORTANT section. The speaker emphasized that "the crucial
problem lies in identifying whether the market structure is trending-based or
volatile-based" and that "good amount of historical data" is required.

Research, compare, and critically evaluate ALL major approaches:

CATEGORY A — STATISTICAL / ECONOMETRIC METHODS

(A1) Hurst Exponent (H)
     - BACKGROUND: The Hurst exponent is a measure of long-term memory of
       time series, originally from hydrology (Hurst 1951), adapted to
       finance. It uses the variance of log prices to assess diffusive
       behavior relative to geometric Brownian motion. The relationship
       is: ⟨|log(t+τ) − log(t)|²⟩ ~ τ^(2H)
     - Calculation methods: Rescaled Range (R/S) analysis, Detrended
       Fluctuation Analysis (DFA), Generalized Hurst Exponent (GHE)
     - KNOWN ISSUE: H is highly sensitive to the choice of lag parameter
       and lookback window. A series can appear mean-reverting at short
       lags and trending at long lags simultaneously. Discuss this
       multi-scale dependency in detail.
     - KNOWN ISSUE: Different estimation methods (R/S vs. DFA vs.
       variance-based) can yield contradictory conclusions on the SAME
       data. Address this and recommend reconciliation approaches.

(A2) Variance Ratio Tests (Lo & MacKinlay 1988)
     - Test whether the variance of k-period returns equals k times the
       variance of 1-period returns (which holds only under random walk)
     - Discuss overlapping vs. non-overlapping estimators and debiasing

(A3) Augmented Dickey-Fuller (ADF) & Phillips-Perron Tests
     - Unit root testing as a proxy for regime identification
     - KPSS test as a complementary approach (tests stationarity directly)
     - Discuss their suitability as real-time vs. ex-post classification

(A4) Autocorrelation Function (ACF) Analysis
     - Direct measurement of serial correlation at multiple lags
     - Rolling-window ACF as a real-time regime proxy

CATEGORY B — MODEL-BASED / PROBABILISTIC METHODS

(B1) Markov Regime-Switching Models (Hamilton 1989)
     - BACKGROUND: Since Hamilton (1989) and Kim & Nelson (2017), the
       assumption of stationarity in market data has been challenged.
       MRS models allow the return-generating process to switch between
       discrete states, each with its own mean and variance parameters.
     - Transition probability matrix estimation and interpretation
     - Limitations: requires pre-specification of number of regimes;
       real-time state inference suffers from look-ahead bias in
       smoothed probabilities vs. filtered probabilities

(B2) Hidden Markov Models (HMM)
     - BACKGROUND: HMMs make two fundamental assumptions: (1) all
       observations depend solely on the current hidden state and are
       conditionally independent of other variables; (2) transition
       probabilities depend only on the current state (Markov property).
       HMMs have been demonstrated to improve portfolio performance by
       yielding higher returns and lower maximum drawdown.
     - Discuss observation variables: daily returns + volatility (e.g.,
       MSE from 10-day sliding window) as a standard configuration
     - Baum-Welch (EM) algorithm for parameter estimation
     - Viterbi algorithm for optimal state sequence decoding
     - Out-of-sample regime filtering for live trading

(B3) Threshold Autoregressive (TAR) / Self-Exciting TAR (SETAR) Models
     - Regime determined by an observable threshold variable
     - Compare advantages vs. latent-state HMM approach

CATEGORY C — TECHNICAL INDICATOR-BASED METHODS

(C1) ADX (Average Directional Index)
     - ADX > 25 as trend indicator; ADX < 20 as range-bound
(C2) Bollinger Band Width as volatility-regime proxy
(C3) Moving Average slope and crossover regime classification
(C4) ATR (Average True Range) regime filtering

For each: discuss simplicity vs. sophistication tradeoff and lag behavior

CATEGORY D — MACHINE LEARNING APPROACHES

(D1) Gaussian Mixture Models (GMM) for unsupervised regime clustering
     - BACKGROUND: GMMs fit distinct Gaussian distributions to capture
       different parts of the return distribution. Entirely data-driven
       with no need for labeled training data. Implementable via
       scikit-learn.
(D2) Random Forest supervised ensemble learning for regime classification
(D3) Online change-point detection: CUSUM, BOCPD (Bayesian Online
     Change Point Detection)
(D4) Wasserstein k-means clustering for regime classification

CATEGORY E — FRACTAL / CHAOS METHODS

(E1) Fractal dimension analysis
(E2) Lyapunov exponents
(E3) Detrended Fluctuation Analysis (DFA)

For EACH method across all categories, provide a structured assessment:

         COMPARATIVE ASSESSMENT TEMPLATE (apply to every method)
         ┌──────────────────────────────────────────────────────┐
         │  Method Name: ___________                           │
         │                                                      │
         │  Historical Data Requirement:                        │
         │    Minimum sample size: ___ bars / observations      │
         │    Optimal lookback window: ___                      │
         │    Granularity sensitivity: tick / minute / daily     │
         │                                                      │
         │  Computational Cost:     □ Low  □ Medium  □ High     │
         │  Real-Time Applicability: □ Yes  □ Partial  □ No     │
         │  Lag / Delay:            ___ bars typical             │
         │  False-Signal Rate:      □ Low  □ Medium  □ High     │
         │  Known Failure Modes:    ___________                 │
         │  Implementation Complexity: □ Low  □ Medium  □ High  │
         │  Academic Support:       □ Strong  □ Moderate □ Weak │
         └──────────────────────────────────────────────────────┘


─── SECTION IV: THE TRANSITION PROBLEM ───

The most dangerous scenario: when the market transitions between regimes.

         REGIME TRANSITION TIMELINE (schematic)
         ════════════════════════════════════════════════════════

         Time ──────────────────────────────────────────────────►

          TRENDING REGIME          TRANSITION         VOLATILE REGIME
         ┌──────────────────┐  ┌──────────────┐  ┌──────────────────┐
         │ Momentum works   │  │  ▓ DANGER ▓   │  │ Mean-reversion   │
         │ H > 0.5          │──│  Both fail    │──│ works            │
         │ Strong ADX       │  │  Whipsaw zone │  │ H < 0.5          │
         │ Clear MA slope   │  │  ADX ≈ 20-25  │  │ Flat MA          │
         └──────────────────┘  │  H ≈ 0.5      │  └──────────────────┘
                               │  Max drawdown  │
                               │  risk          │
                               └──────────────┘

         ════════════════════════════════════════════════════════

Research the following:

1. EARLY WARNING SIGNALS of regime change:
   - Declining ADX from trending territory
   - Hurst exponent crossing through 0.5
   - Volatility expansion / compression signals
   - Breakdown in autocorrelation structure
   - Changes in volume patterns preceding regime shifts

2. ADAPTIVE / ENSEMBLE STRATEGIES for graceful degradation:
   - Strategy blending based on regime probability
   - Position sizing as a function of regime confidence
   - Stop-loss tightening during ambiguous regime periods
   - Ensemble models combining multiple regime indicators

3. HISTORICAL CASE STUDIES of regime transition failures:
   - 2007-2008 financial crisis (trending → volatile)
   - COVID-19 crash of March 2020
   - Flash crashes and micro-regime transitions


─── SECTION V: PRACTICAL IMPLEMENTATION FRAMEWORK ───

Propose a COMPOSITE REGIME-IDENTIFICATION SCORING SYSTEM:

         COMPOSITE REGIME SCORE ARCHITECTURE
         ════════════════════════════════════════════════════════

         Input Signals                  Aggregation         Output
         ┌─────────────────┐
         │ Hurst Exponent  │──┐
         ├─────────────────┤  │
         │ Variance Ratio  │──┤      ┌──────────────┐     ┌──────────┐
         ├─────────────────┤  ├─────►│  Weighted     │────►│ Regime   │
         │ ADX Level       │──┤      │  Ensemble     │     │ Score    │
         ├─────────────────┤  ├─────►│  Scoring      │────►│ [-1, +1]│
         │ MA Slope         │──┤     │  Engine       │     │          │
         ├─────────────────┤  │      └──────────────┘     │ -1 = MR  │
         │ HMM State Prob  │──┤                           │  0 = RW  │
         ├─────────────────┤  │                           │ +1 = TR  │
         │ Autocorrelation │──┘                           └──────────┘
         └─────────────────┘
                                                           │
                                                           ▼
                                                   ┌──────────────┐
                                                   │ Strategy      │
                                                   │ Allocation    │
                                                   │ Engine        │
                                                   │               │
                                                   │ Momentum wt:  │
                                                   │ max(0, Score) │
                                                   │ MR weight:    │
                                                   │ max(0,-Score) │
                                                   └──────────────┘

         ════════════════════════════════════════════════════════

Address:
- Recommended data granularity (tick / minute / hourly / daily) for
  different holding periods
- Minimum historical sample sizes for reliable classification
- Parameter calibration approaches (walk-forward optimization)
- Backtesting methodology: walk-forward out-of-sample testing to
  avoid look-ahead bias
- Confidence thresholds: when the composite score is too ambiguous
  to act upon, recommend a "flat" or reduced-exposure posture


=== ADDITIONAL INSTRUCTIONS ===

- Where academic papers or known empirical results are cited, provide the
  full Author(s), Year, and Paper Title so they can be located
- If any claim or method has known controversy or conflicting evidence,
  present BOTH sides objectively
- For every strategy discussed, explicitly state the assumptions under
  which it works AND the conditions under which it fails
- Flag any area where your knowledge may be incomplete or outdated and
  recommend further investigation
- This report should connect forward to Topics 2 and 3: briefly note that
  regime identification (this topic) influences factor timing (Topic 2) and
  that volume patterns analyzed here connect to order flow analysis (Topic 3)
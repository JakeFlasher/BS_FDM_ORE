## 1) `ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp`

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
    Copyright (C) 2026 QuantLib contributors
*/

/*! \file fdmfittedblackscholesop.hpp
    \brief Exponentially fitted Black-Scholes linear operator (log-space)
*/

#ifndef quantlib_fdm_fitted_black_scholes_op_hpp
#define quantlib_fdm_fitted_black_scholes_op_hpp

#include <ql/processes/blackscholesprocess.hpp>
#include <ql/utilities/null.hpp>
#include <ql/methods/finitedifferences/utilities/fdmquantohelper.hpp>
#include <ql/methods/finitedifferences/operators/firstderivativeop.hpp>
#include <ql/methods/finitedifferences/operators/triplebandlinearop.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp>

#include <vector>

namespace QuantLib {

    //! Exponentially fitted linear operator for the log-space Black-Scholes PDE
    /*! This operator implements Duffy's exponentially fitted diffusion
        coefficient in QuantLib's log-space convention \(x = \ln(S)\).

        The fitted diffusion is computed per node using the local Péclet number
        and the fitting factor
        \f[
            \rho = Pe \coth(Pe) = \frac{Pe}{\tanh(Pe)}.
        \f]

        \note When local volatility is enabled, the fitting factor is recomputed
              per node and per time step in setTime().
    */
    class FdmFittedBlackScholesOp : public FdmLinearOpComposite {
      public:
        FdmFittedBlackScholesOp(
            const ext::shared_ptr<FdmMesher>& mesher,
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
            Real strike,
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>(),
            Size direction = 0,
            ext::shared_ptr<FdmQuantoHelper> quantoHelper =
                ext::shared_ptr<FdmQuantoHelper>());

        Size size() const override;
        void setTime(Time t1, Time t2) override;

        Array apply(const Array& r) const override;
        Array apply_mixed(const Array& r) const override;
        Array apply_direction(Size direction, const Array& r) const override;
        Array solve_splitting(Size direction, const Array& r, Real s) const override;
        Array preconditioner(const Array& r, Real s) const override;

        std::vector<SparseMatrix> toMatrixDecomp() const override;

        //! True if all interior-node off-diagonals of the current mapT_ are non-negative.
        bool mMatrixSatisfied() const;

        //! Number of interior nodes where lower or upper off-diagonal is negative.
        Size mMatrixViolationCount() const;

      private:
        static Real fittingFactor(Real Pe);

        const ext::shared_ptr<FdmMesher> mesher_;
        const ext::shared_ptr<YieldTermStructure> rTS_, qTS_;
        const ext::shared_ptr<BlackVolTermStructure> volTS_;
        const ext::shared_ptr<LocalVolTermStructure> localVol_;
        const Array x_;
        const FirstDerivativeOp  dxMap_;
        const TripleBandLinearOp dxxMap_;
        TripleBandLinearOp mapT_;
        const Real strike_;
        const Real illegalLocalVolOverwrite_;
        const Size direction_;
        const ext::shared_ptr<FdmQuantoHelper> quantoHelper_;
    };
}

#endif
```

---

## 2) `ql/methods/finitedifferences/operators/fdmfittedblackscholesop.cpp`

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
    Copyright (C) 2026 QuantLib contributors
*/

#include <ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp>

#include <ql/errors.hpp>
#include <ql/math/array.hpp>
#include <ql/math/functional.hpp>
#include <ql/mathconstants.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/operators/modtriplebandlinearop.hpp>
#include <ql/methods/finitedifferences/operators/secondderivativeop.hpp>
#include <ql/utilities/null.hpp>

#include <algorithm>
#include <cmath>
#include <utility>

namespace QuantLib {

    namespace {
        Real interiorSpacing(const ext::shared_ptr<FdmMesher>& mesher,
                             const FdmLinearOpIterator& iter,
                             Size direction) {
            const Real hm = mesher->dminus(iter, direction);
            const Real hp = mesher->dplus(iter, direction);

            if (hm == Null<Real>() || hp == Null<Real>())
                return Null<Real>();

            return 0.5*(hm + hp);
        }

        bool isInteriorNode(const ext::shared_ptr<FdmMesher>& mesher,
                            const FdmLinearOpIterator& iter,
                            Size direction) {
            const Size co  = iter.coordinates()[direction];
            const Size dim = mesher->layout()->dim()[direction];
            return (co != 0U && co != dim - 1U);
        }
    }

    FdmFittedBlackScholesOp::FdmFittedBlackScholesOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& bsProcess,
        Real strike,
        bool localVol,
        Real illegalLocalVolOverwrite,
        Size direction,
        ext::shared_ptr<FdmQuantoHelper> quantoHelper)
    : mesher_(mesher),
      rTS_(bsProcess->riskFreeRate().currentLink()),
      qTS_(bsProcess->dividendYield().currentLink()),
      volTS_(bsProcess->blackVolatility().currentLink()),
      localVol_((localVol) ? bsProcess->localVolatility().currentLink()
                           : ext::shared_ptr<LocalVolTermStructure>()),
      x_((localVol) ? Array(Exp(mesher->locations(direction))) : Array()),
      dxMap_(FirstDerivativeOp(direction, mesher)),
      dxxMap_(SecondDerivativeOp(direction, mesher)),
      mapT_(direction, mesher),
      strike_(strike),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      direction_(direction),
      quantoHelper_(std::move(quantoHelper)) {}

    Size FdmFittedBlackScholesOp::size() const { return 1U; }

    Real FdmFittedBlackScholesOp::fittingFactor(Real Pe) {
        const Real absPe = std::fabs(Pe);

        // Taylor: Pe*coth(Pe) = 1 + Pe^2/3 - Pe^4/45 + ...
        if (absPe < 1e-8)
            return 1.0 + (Pe*Pe)/3.0;

        // tanh(x) saturates to +/-1 for moderate |x| in IEEE-754 double,
        // so avoid returning Pe/tanh(Pe) when tanh() is exactly +/-1.
        if (absPe > 300.0)
            return absPe;

        return Pe/std::tanh(Pe);
    }

    void FdmFittedBlackScholesOp::setTime(Time t1, Time t2) {
        const Rate r = rTS_->forwardRate(t1, t2, Continuous).rate();
        const Rate q = qTS_->forwardRate(t1, t2, Continuous).rate();

        const Size n = mesher_->layout()->size();
        Array convection(n), fittedDiffusion(n);

        // Guard against division by (near) zero when computing Pe = b*h/(2*a).
        // This also provides a numerically stable way to realize the a->0 upwind limit.
        const Real minA = 1e-20;

        if (localVol_ != nullptr) {

            Array variance(n);
            const Time tMid = 0.5*(t1 + t2);

            for (const auto& iter : *mesher_->layout()) {
                const Size i = iter.index();

                if (illegalLocalVolOverwrite_ < 0.0) {
                    variance[i] = squared(localVol_->localVol(tMid, x_[i], true));
                } else {
                    try {
                        variance[i] = squared(localVol_->localVol(tMid, x_[i], true));
                    } catch (Error&) {
                        variance[i] = squared(illegalLocalVolOverwrite_);
                    }
                }
            }

            Array quantoAdj;
            Real quantoAdjScalar = 0.0;
            bool quantoAdjIsScalar = false;

            if (quantoHelper_ != nullptr) {
                quantoAdj = quantoHelper_->quantoAdjustment(Sqrt(variance), t1, t2);
                QL_REQUIRE(!quantoAdj.empty(), "empty quanto adjustment");

                if (quantoAdj.size() == 1U) {
                    quantoAdjIsScalar = true;
                    quantoAdjScalar = quantoAdj[0];
                } else {
                    QL_REQUIRE(quantoAdj.size() == n,
                               "inconsistent quanto adjustment size ("
                                   << quantoAdj.size() << ", expected " << n << ")");
                }
            }

            for (const auto& iter : *mesher_->layout()) {
                const Size i = iter.index();

                const Real aRaw = 0.5*variance[i];               // physical a = σ^2/2
                const Real a = std::max(aRaw, minA);              // guarded diffusion for Pe and a*rho
                Real b = (r - q) - aRaw;                          // log-space drift coefficient
                if (quantoHelper_ != nullptr)
                    b -= (quantoAdjIsScalar ? quantoAdjScalar : quantoAdj[i]);

                convection[i] = b;

                const Real h = interiorSpacing(mesher_, iter, direction_);
                if (h == Null<Real>()) {
                    fittedDiffusion[i] = a; // boundary node: skip fitting
                    continue;
                }

                const Real Pe = b*h/(2.0*a);
                Real rho = fittingFactor(Pe);

                // Enforce the theoretical lower bound rho >= 1 (protect against round-off).
                rho = std::max<Real>(1.0, rho);

                fittedDiffusion[i] = a*rho;
            }

        } else {

            const Real v = volTS_->blackForwardVariance(t1, t2, strike_)/(t2 - t1);
            const Real aRaw = 0.5*v;
            const Real a = std::max(aRaw, minA);
            Real b = (r - q) - aRaw;

            if (quantoHelper_ != nullptr) {
                const Array qa = quantoHelper_->quantoAdjustment(
                    Array(1, std::sqrt(v)), t1, t2);
                QL_REQUIRE(!qa.empty(), "empty quanto adjustment");
                b -= qa[0];
            }

            for (const auto& iter : *mesher_->layout()) {
                const Size i = iter.index();
                convection[i] = b;

                const Real h = interiorSpacing(mesher_, iter, direction_);
                if (h == Null<Real>()) {
                    fittedDiffusion[i] = a; // boundary node: skip fitting
                    continue;
                }

                const Real Pe = b*h/(2.0*a);
                Real rho = fittingFactor(Pe);
                rho = std::max<Real>(1.0, rho);

                fittedDiffusion[i] = a*rho;
            }

#ifdef QL_EXTRA_SAFETY_CHECKS
            // For constant coefficients on a uniform log-mesh, Pe should be constant
            // across interior nodes. We only check this in the truly-uniform case.
            Real hRef = Null<Real>();
            bool uniform = true;

            for (const auto& iter : *mesher_->layout()) {
                if (!isInteriorNode(mesher_, iter, direction_))
                    continue;

                const Real h = interiorSpacing(mesher_, iter, direction_);
                if (h == Null<Real>())
                    continue;

                if (hRef == Null<Real>()) {
                    hRef = h;
                } else {
                    const Real tolH =
                        50.0*std::sqrt(QL_EPSILON)*std::max(1.0, std::fabs(hRef));
                    if (std::fabs(h - hRef) > tolH) {
                        uniform = false;
                        break;
                    }
                }
            }

            if (uniform && hRef != Null<Real>()) {
                bool peInit = false;
                Real peRef = 0.0;

                for (const auto& iter : *mesher_->layout()) {
                    if (!isInteriorNode(mesher_, iter, direction_))
                        continue;

                    const Real h = interiorSpacing(mesher_, iter, direction_);
                    if (h == Null<Real>())
                        continue;

                    const Real Pe = b*h/(2.0*a);

                    if (!peInit) {
                        peInit = true;
                        peRef = Pe;
                    } else {
                        const Real tolPe =
                            100.0*std::sqrt(QL_EPSILON)*std::max(1.0, std::fabs(peRef));
                        QL_REQUIRE(std::fabs(Pe - peRef) <= tolPe,
                                   "internal error: Pe is not constant on a uniform mesh");
                    }
                }
            }
#endif
        }

        mapT_.axpyb(convection, dxMap_,
                    dxxMap_.mult(fittedDiffusion),
                    Array(1, -r));
    }

    Array FdmFittedBlackScholesOp::apply(const Array& u) const {
        return mapT_.apply(u);
    }

    Array FdmFittedBlackScholesOp::apply_direction(Size direction,
                                                   const Array& r) const {
        if (direction == direction_)
            return mapT_.apply(r);
        else
            return Array(r.size(), 0.0);
    }

    Array FdmFittedBlackScholesOp::apply_mixed(const Array& r) const {
        return Array(r.size(), 0.0);
    }

    Array FdmFittedBlackScholesOp::solve_splitting(Size direction,
                                                   const Array& r,
                                                   Real dt) const {
        if (direction == direction_)
            return mapT_.solve_splitting(r, dt, 1.0);
        else
            return r;
    }

    Array FdmFittedBlackScholesOp::preconditioner(const Array& r,
                                                  Real dt) const {
        return solve_splitting(direction_, r, dt);
    }

    std::vector<SparseMatrix> FdmFittedBlackScholesOp::toMatrixDecomp() const {
        return std::vector<SparseMatrix>(1, mapT_.toMatrix());
    }

    bool FdmFittedBlackScholesOp::mMatrixSatisfied() const {
        return mMatrixViolationCount() == 0U;
    }

    Size FdmFittedBlackScholesOp::mMatrixViolationCount() const {
        const ModTripleBandLinearOp mod(mapT_);

        // Tolerance for tiny negative values due to rounding.
        const Real tol = 100.0*QL_EPSILON;

        Size count = 0U;
        for (const auto& iter : *mesher_->layout()) {
            if (!isInteriorNode(mesher_, iter, direction_))
                continue;

            const Size i = iter.index();
            if (mod.lower(i) < -tol || mod.upper(i) < -tol)
                ++count;
        }
        return count;
    }

}
```

---

## 3) `ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp`

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
    Copyright (C) 2026 QuantLib contributors
*/

/*! \file fdmdiagnostics.hpp
    \brief Lightweight diagnostics utilities for FDM solution verification
*/

#ifndef quantlib_fdm_diagnostics_hpp
#define quantlib_fdm_diagnostics_hpp

#include <ql/math/array.hpp>
#include <ql/types.hpp>

namespace QuantLib {

    struct FdmDiagnosticsReport {
        Real minValue = 0.0;
        Size negativeCount = 0U;
        Real oscillationScore = 0.0;
        Size mMatrixViolationCount = 0U;
        Size nanCount = 0U;
    };

    class FdmDiagnostics {
      public:
        enum Level { Off, Light, Full };

        explicit FdmDiagnostics(Level level = Off);

        FdmDiagnosticsReport checkSolution(const Array& u) const;
        // Light: minValue, negativeCount, nanCount only
        // Full: all fields including oscillationScore

        static Real oscillationScore(const Array& u);
        // Count sign changes in Δu_j = u[j+1]−u[j], ignoring |Δu| < 1e-15.
        // Normalize by max(1, size−2). Score of 0 = monotone; >0.1 = oscillating.
        // Algorithm is O(N) with no heap allocation.

        static FdmDiagnosticsReport merge(const FdmDiagnosticsReport& a,
                                         const FdmDiagnosticsReport& b);

        Level level() const;

      private:
        Level level_;
    };

}

#endif
```

---

## 4) `ql/methods/finitedifferences/utilities/fdmdiagnostics.cpp`

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
    Copyright (C) 2026 QuantLib contributors
*/

#include <ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp>

#include <ql/utilities/null.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace QuantLib {

    FdmDiagnostics::FdmDiagnostics(Level level) : level_(level) {}

    FdmDiagnostics::Level FdmDiagnostics::level() const { return level_; }

    FdmDiagnosticsReport FdmDiagnostics::checkSolution(const Array& u) const {
        FdmDiagnosticsReport rep;

        if (level_ == Off) {
            rep.minValue = Null<Real>();
            rep.negativeCount = 0U;
            rep.nanCount = 0U;
            rep.oscillationScore = Null<Real>();
            rep.mMatrixViolationCount = 0U;
            return rep;
        }

        rep.minValue = std::numeric_limits<Real>::infinity();
        rep.negativeCount = 0U;
        rep.nanCount = 0U;
        rep.oscillationScore = 0.0;
        rep.mMatrixViolationCount = 0U;

        for (Size i = 0U; i < u.size(); ++i) {
            const Real x = u[i];

            if (!std::isfinite(x)) {
                ++rep.nanCount;
                continue;
            }

            rep.minValue = std::min(rep.minValue, x);

            if (x < 0.0)
                ++rep.negativeCount;
        }

        if (u.empty() || rep.minValue == std::numeric_limits<Real>::infinity())
            rep.minValue = Null<Real>();

        if (level_ == Full) {
            rep.oscillationScore =
                (rep.nanCount == 0U) ? FdmDiagnostics::oscillationScore(u) : Null<Real>();
        } else {
            rep.oscillationScore = Null<Real>();
        }

        return rep;
    }

    Real FdmDiagnostics::oscillationScore(const Array& u) {
        const Size n = u.size();
        if (n < 3U)
            return 0.0;

        // Ignore tiny differences to avoid counting noise as oscillations.
        const Real eps = 1e-15;

        Integer prevSign = 0; // -1, 0, +1
        Size changes = 0U;

        for (Size j = 0U; j + 1U < n; ++j) {
            const Real a = u[j];
            const Real b = u[j + 1U];

            if (!std::isfinite(a) || !std::isfinite(b))
                return Null<Real>();

            const Real du = b - a;
            if (std::fabs(du) < eps)
                continue;

            const Integer sign = (du > 0.0) ? 1 : -1;

            if (prevSign != 0 && sign != prevSign)
                ++changes;

            prevSign = sign;
        }

        const Real denom = std::max<Real>(1.0, static_cast<Real>(n - 2U));
        return static_cast<Real>(changes) / denom;
    }

    FdmDiagnosticsReport FdmDiagnostics::merge(const FdmDiagnosticsReport& a,
                                              const FdmDiagnosticsReport& b) {
        FdmDiagnosticsReport r;

        // minValue: take the minimum among available values
        if (a.minValue == Null<Real>())
            r.minValue = b.minValue;
        else if (b.minValue == Null<Real>())
            r.minValue = a.minValue;
        else
            r.minValue = std::min(a.minValue, b.minValue);

        r.negativeCount = a.negativeCount + b.negativeCount;
        r.nanCount = a.nanCount + b.nanCount;
        r.mMatrixViolationCount = a.mMatrixViolationCount + b.mMatrixViolationCount;

        // oscillationScore: keep the worst (largest) available score
        if (a.oscillationScore == Null<Real>())
            r.oscillationScore = b.oscillationScore;
        else if (b.oscillationScore == Null<Real>())
            r.oscillationScore = a.oscillationScore;
        else
            r.oscillationScore = std::max(a.oscillationScore, b.oscillationScore);

        return r;
    }

}
```

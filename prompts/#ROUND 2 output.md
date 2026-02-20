## 1) `ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.hpp`

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
    Copyright (C) 2026 QuantLib contributors
*/

/*! \file fdmcnvariantblackscholesop.hpp
    \brief Crank-Nicolson variant Black-Scholes linear operator (log-space)
*/

#ifndef quantlib_fdm_cn_variant_black_scholes_op_hpp
#define quantlib_fdm_cn_variant_black_scholes_op_hpp

#include <ql/processes/blackscholesprocess.hpp>
#include <ql/utilities/null.hpp>
#include <ql/methods/finitedifferences/utilities/fdmquantohelper.hpp>
#include <ql/methods/finitedifferences/operators/firstderivativeop.hpp>
#include <ql/methods/finitedifferences/operators/triplebandlinearop.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp>

#include <vector>

namespace QuantLib {

    //! Crank-Nicolson variant spatial operator for the log-space Black-Scholes PDE
    /*! This operator implements Milev–Tagliani's Crank–Nicolson variant
        via an *effective diffusion* reformulation in QuantLib's log-space
        convention \(x = \ln(S)\).

        In log-space, QuantLib's Black-Scholes PDE is
        \[
            -u_\tau + b u_x + a u_{xx} - r u = 0,
        \]
        where \(a=\sigma^2/2\) and \(b=(r-q)-\sigma^2/2\).

        The Milev–Tagliani CN variant modifies the reaction term stencil. When
        mapped to QuantLib's standard CN time-stepping (which halves the operator),
        this is equivalent to using an enhanced diffusion coefficient at each node:
        \[
            a_{\text{eff}} = \frac{\sigma^2}{2} + \frac{r^2 h^2}{8\sigma^2},
        \]
        where \(h\) is the (local) log-space mesh spacing.

        \note On nonuniform meshes, a per-node effective spacing
              \(h_j = \tfrac{1}{2}(d^-_j + d^+_j)\) is used, skipping boundary nodes
              where \(d^\pm\) is Null<Real>().

        \warning With local volatility enabled, \(\sigma=\sigma(S,t)\) and the
                 CN-variant parameter \(\omega=-r/(16\sigma^2)\) becomes spatially
                 varying; the scheme is still assembled via \(a_{\text{eff}}\) per node,
                 but the theoretical guarantees and timestep constraint become
                 node-dependent.
    */
    class FdmCNVariantBlackScholesOp : public FdmLinearOpComposite {
      public:
        FdmCNVariantBlackScholesOp(
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

        //! Returns the last computed \f$\omega=-r/(16\sigma^2)\f$ (scalar reference value).
        /*! For constant volatility, this is exact for the current time step.
            For local volatility, it is computed at spot = strike_ (diagnostic only).
        */
        Real omega() const;

        //! Returns an estimated max timestep for non-negative explicit-side diagonal (CE-20).
        /*! For constant volatility this reduces to the (single) worst-case interior node.
            On concentrated meshes, very small local spacing can make the returned value
            impractically tight.
        */
        Real timestepConstraint() const;

        //! True if all interior-node off-diagonals of the current mapT_ are non-negative.
        bool mMatrixSatisfied() const;

        //! Number of interior nodes where lower or upper off-diagonal is negative.
        Size mMatrixViolationCount() const;

      private:
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

        Real omega_;
        Real timestepConstraint_;
    };
}

#endif
```

---

## 2) `ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.cpp`

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
    Copyright (C) 2026 QuantLib contributors
*/

#include <ql/methods/finitedifferences/operators/fdmcnvariantblackscholesop.hpp>

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

    FdmCNVariantBlackScholesOp::FdmCNVariantBlackScholesOp(
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
      quantoHelper_(std::move(quantoHelper)),
      omega_(Null<Real>()),
      timestepConstraint_(Null<Real>()) {}

    Size FdmCNVariantBlackScholesOp::size() const { return 1U; }

    void FdmCNVariantBlackScholesOp::setTime(Time t1, Time t2) {

        const Rate r = rTS_->forwardRate(t1, t2, Continuous).rate();
        const Rate q = qTS_->forwardRate(t1, t2, Continuous).rate();

        const Size n = mesher_->layout()->size();
        Array convection(n), effectiveDiffusion(n);

        // Guard to avoid division by (near) zero in the r^2 h^2 / (8 sigma^2) term.
        // The scheme is intended for low-volatility, but not for exactly zero volatility.
        const Real minSigma2 = 1e-20;

        Real hMin = Null<Real>();
        Real aEffAtHMin = Null<Real>();

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

                const Real sigma2Raw = variance[i];
                const Real sigma2Den = std::max(sigma2Raw, minSigma2);

                const Real aRaw = 0.5*sigma2Raw;     // physical a = σ^2/2 in log-space
                Real b = (r - q) - aRaw;             // log-space drift coefficient

                if (quantoHelper_ != nullptr)
                    b -= (quantoAdjIsScalar ? quantoAdjScalar : quantoAdj[i]);

                convection[i] = b;

                const Real h = interiorSpacing(mesher_, iter, direction_);
                if (h == Null<Real>()) {
                    effectiveDiffusion[i] = aRaw; // boundary node: h not defined
                    continue;
                }

                const Real aEff = aRaw + (r*r*h*h)/(8.0*sigma2Den);
                effectiveDiffusion[i] = aEff;

                if (hMin == Null<Real>() || h < hMin) {
                    hMin = h;
                    aEffAtHMin = aEff;
                }
            }

            // Scalar diagnostic omega: compute at spot = strike_.
            Real sigma2Strike = Null<Real>();
            if (illegalLocalVolOverwrite_ < 0.0) {
                sigma2Strike = squared(localVol_->localVol(tMid, strike_, true));
            } else {
                try {
                    sigma2Strike = squared(localVol_->localVol(tMid, strike_, true));
                } catch (Error&) {
                    sigma2Strike = squared(illegalLocalVolOverwrite_);
                }
            }
            omega_ = -r/(16.0*std::max(sigma2Strike, minSigma2));

        } else {

            const Real v = volTS_->blackForwardVariance(t1, t2, strike_)/(t2 - t1);
            const Real sigma2Den = std::max(v, minSigma2);

            const Real aRaw = 0.5*v;
            Real b = (r - q) - aRaw;

            if (quantoHelper_ != nullptr) {
                const Array qa =
                    quantoHelper_->quantoAdjustment(Array(1, std::sqrt(v)), t1, t2);
                QL_REQUIRE(!qa.empty(), "empty quanto adjustment");
                QL_REQUIRE(qa.size() == 1U, "unexpected quanto adjustment size");
                b -= qa[0];
            }

            omega_ = -r/(16.0*sigma2Den);

            for (const auto& iter : *mesher_->layout()) {
                const Size i = iter.index();
                convection[i] = b;

                const Real h = interiorSpacing(mesher_, iter, direction_);
                if (h == Null<Real>()) {
                    effectiveDiffusion[i] = aRaw; // boundary node: h not defined
                    continue;
                }

                const Real aEff = aRaw + (r*r*h*h)/(8.0*sigma2Den);
                effectiveDiffusion[i] = aEff;

                if (hMin == Null<Real>() || h < hMin) {
                    hMin = h;
                    aEffAtHMin = aEff;
                }
            }
        }

        // CE-20 timestep constraint (log-space): require explicit-side diagonal >= 0
        //  1 - 0.5*dt*(2*a_eff/h^2 + r) >= 0  -> dt <= 1 / (a_eff/h^2 + r/2).
        //
        // We use the smallest interior-node h as the practical worst case on typical
        // meshes; on concentrated meshes this can be extremely restrictive.
        if (hMin != Null<Real>() && aEffAtHMin != Null<Real>() && hMin > 0.0) {
            const Real denom = (aEffAtHMin/(hMin*hMin)) + 0.5*r;
            timestepConstraint_ = (denom > 0.0) ? 1.0/denom : Null<Real>();
        } else {
            timestepConstraint_ = Null<Real>();
        }

        mapT_.axpyb(convection, dxMap_,
                    dxxMap_.mult(effectiveDiffusion),
                    Array(1, -r));
    }

    Array FdmCNVariantBlackScholesOp::apply(const Array& r) const {
        return mapT_.apply(r);
    }

    Array FdmCNVariantBlackScholesOp::apply_direction(Size direction,
                                                      const Array& r) const {
        if (direction == direction_)
            return mapT_.apply(r);
        else
            return Array(r.size(), 0.0);
    }

    Array FdmCNVariantBlackScholesOp::apply_mixed(const Array& r) const {
        return Array(r.size(), 0.0);
    }

    Array FdmCNVariantBlackScholesOp::solve_splitting(Size direction,
                                                      const Array& r,
                                                      Real dt) const {
        if (direction == direction_)
            return mapT_.solve_splitting(r, dt, 1.0);
        else
            return r;
    }

    Array FdmCNVariantBlackScholesOp::preconditioner(const Array& r,
                                                     Real dt) const {
        return solve_splitting(direction_, r, dt);
    }

    std::vector<SparseMatrix> FdmCNVariantBlackScholesOp::toMatrixDecomp() const {
        return std::vector<SparseMatrix>(1, mapT_.toMatrix());
    }

    Real FdmCNVariantBlackScholesOp::omega() const {
        return omega_;
    }

    Real FdmCNVariantBlackScholesOp::timestepConstraint() const {
        return timestepConstraint_;
    }

    bool FdmCNVariantBlackScholesOp::mMatrixSatisfied() const {
        return mMatrixViolationCount() == 0U;
    }

    Size FdmCNVariantBlackScholesOp::mMatrixViolationCount() const {
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

## 3) `ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp`

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
    Copyright (C) 2026 QuantLib contributors
*/

/*! \file fdmbarrierprojectioncondition.hpp
    \brief Discrete barrier projection (corridor indicator) step condition
*/

#ifndef quantlib_fdm_barrier_projection_condition_hpp
#define quantlib_fdm_barrier_projection_condition_hpp

#include <ql/methods/finitedifferences/stepcondition.hpp>
#include <ql/shared_ptr.hpp>
#include <ql/types.hpp>
#include <ql/utilities/null.hpp>

#include <vector>

namespace QuantLib {

    class FdmMesher;

    //! Discrete barrier projection condition: U <- U * 1_[L,U](S)
    /*! Implements discrete monitoring of a double (or single) barrier by
        projecting the solution values to zero outside a corridor at monitoring
        times:
        \[
            U_j \leftarrow U_j \cdot \mathbf{1}_{[L,U]}(S_j).
        \]

        \warning QuantLib's Black-Scholes mesher stores x = ln(S) as spatial
                 coordinate (see FdmBlackScholesMesher and FdmBlackScholesOp).
                 This condition therefore converts S-space barriers to log-space
                 via ln(L), ln(U) and compares against mesher->location(iter, direction).

        Grid alignment note:
        For accurate enforcement, ln(L) and ln(U) should ideally coincide with
        mesh points. This condition uses strict comparisons against the mesh
        coordinates, so misalignment can leave small residual values near barriers.
    */
    class FdmBarrierProjectionCondition : public StepCondition<Array> {
      public:
        FdmBarrierProjectionCondition(
            std::vector<Time> monitoringTimes,
            Real lowerBarrier,  // in S-space; 0.0 means "no lower barrier"
            Real upperBarrier,  // in S-space; Null<Real>() means "no upper barrier"
            const ext::shared_ptr<FdmMesher>& mesher,
            Size direction = 0);

        void applyTo(Array& a, Time t) const override;

        const std::vector<Time>& monitoringTimes() const;

      private:
        bool isMonitoringTime(Time t) const;

        std::vector<Time> monitoringTimes_;
        std::vector<Size> outsideIndices_;
        Size size_;
    };

}

#endif
```

---

## 4) `ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.cpp`

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
    Copyright (C) 2026 QuantLib contributors
*/

#include <ql/methods/finitedifferences/stepconditions/fdmbarrierprojectioncondition.hpp>

#include <ql/errors.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/utilities/null.hpp>

#include <algorithm>
#include <cmath>
#include <utility>

namespace QuantLib {

    FdmBarrierProjectionCondition::FdmBarrierProjectionCondition(
        std::vector<Time> monitoringTimes,
        Real lowerBarrier,
        Real upperBarrier,
        const ext::shared_ptr<FdmMesher>& mesher,
        Size direction)
    : monitoringTimes_(std::move(monitoringTimes)),
      outsideIndices_(),
      size_(mesher->layout()->size()) {

        QL_REQUIRE(lowerBarrier >= 0.0, "lower barrier must be non-negative");
        if (upperBarrier != Null<Real>())
            QL_REQUIRE(upperBarrier > 0.0, "upper barrier must be positive");

        const bool hasLower = (lowerBarrier > 0.0);
        const bool hasUpper = (upperBarrier != Null<Real>());

        if (hasLower && hasUpper)
            QL_REQUIRE(upperBarrier > lowerBarrier,
                       "upper barrier (" << upperBarrier
                                         << ") must be greater than lower barrier ("
                                         << lowerBarrier << ")");

        // Normalize monitoring times (sorted, unique).
        std::sort(monitoringTimes_.begin(), monitoringTimes_.end());
        monitoringTimes_.erase(
            std::unique(monitoringTimes_.begin(), monitoringTimes_.end()),
            monitoringTimes_.end());

        for (const Time t : monitoringTimes_)
            QL_REQUIRE(t >= 0.0, "negative monitoring time given");

        // Convert S-space barriers to log-space (x = ln(S)).
        Real lnLower = Null<Real>();
        Real lnUpper = Null<Real>();

        if (hasLower) {
            // lowerBarrier > 0, enforced above
            lnLower = std::log(lowerBarrier);
        }
        if (hasUpper) {
            // upperBarrier > 0, enforced above
            lnUpper = std::log(upperBarrier);
        }

        const Real tol = 1e-12;

        // Precompute indices outside the barrier corridor.
        if (hasLower || hasUpper) {
            outsideIndices_.reserve(size_);

            for (const auto& iter : *mesher->layout()) {
                const Size i = iter.index();
                const Real x = mesher->location(iter, direction);

                bool outside = false;

                if (hasLower && x < lnLower - tol)
                    outside = true;
                if (hasUpper && x > lnUpper + tol)
                    outside = true;

                if (outside)
                    outsideIndices_.push_back(i);
            }
        }
    }

    const std::vector<Time>&
    FdmBarrierProjectionCondition::monitoringTimes() const {
        return monitoringTimes_;
    }

    bool FdmBarrierProjectionCondition::isMonitoringTime(Time t) const {
        if (monitoringTimes_.empty())
            return false;

        // Matching tolerance required by spec.
        const Real tol = 1e-10;

        const auto it = std::lower_bound(monitoringTimes_.begin(),
                                         monitoringTimes_.end(), t);

        if (it != monitoringTimes_.end() && std::fabs(*it - t) <= tol)
            return true;

        if (it != monitoringTimes_.begin()) {
            const auto itPrev = std::prev(it);
            if (std::fabs(*itPrev - t) <= tol)
                return true;
        }

        return false;
    }

    void FdmBarrierProjectionCondition::applyTo(Array& a, Time t) const {

        if (outsideIndices_.empty() || monitoringTimes_.empty())
            return;

        QL_REQUIRE(a.size() == size_,
                   "inconsistent solution size (" << a.size()
                                                  << ", expected " << size_ << ")");

        if (!isMonitoringTime(t))
            return;

        for (Size k = 0U; k < outsideIndices_.size(); ++k)
            a[outsideIndices_[k]] = 0.0;
    }

}
```

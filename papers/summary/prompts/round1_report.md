
**File 1: `ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp`**

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 QuantLib contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file fdmfittedblackscholesop.hpp
    \brief Black Scholes linear operator with exponential fitting
*/

#ifndef quantlib_fdm_fitted_black_scholes_op_hpp
#define quantlib_fdm_fitted_black_scholes_op_hpp

#include <ql/processes/blackscholesprocess.hpp>
#include <ql/methods/finitedifferences/utilities/fdmquantohelper.hpp>
#include <ql/methods/finitedifferences/operators/firstderivativeop.hpp>
#include <ql/methods/finitedifferences/operators/triplebandlinearop.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearopcomposite.hpp>

namespace QuantLib {

    //! Black-Scholes operator with exponentially fitted spatial discretization.
    /*! Replaces the standard centered-difference stencil of
        FdmBlackScholesOp with an Il'in / Duffy-style exponential fitting
        factor that preserves the M-matrix sign pattern of the implicit
        system matrix.  Constructor signature matches FdmBlackScholesOp
        so the two classes are drop-in replacements for each other.

        The fitting factor at interior node \f$i\f$ is
        \f[
            \rho_i = \theta_i \coth(\theta_i), \qquad
            \theta_i = \frac{b_i\,h_i}{2\,a_i},
        \f]
        where \f$a_i = \tfrac12\sigma_i^2\f$ (diffusion),
        \f$b_i = (r-q) - a_i\f$ (drift in log-space), and \f$h_i\f$
        is the effective local grid spacing.
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
            ext::shared_ptr<FdmQuantoHelper> quantoHelper
                = ext::shared_ptr<FdmQuantoHelper>());

        Size size() const override;
        void setTime(Time t1, Time t2) override;

        Disposable<Array> apply(const Array& r) const override;
        Disposable<Array> apply_mixed(const Array& r) const override;
        Disposable<Array> apply_direction(Size direction,
                                          const Array& r) const override;
        Disposable<Array> solve_splitting(Size direction,
                                          const Array& r, Real s) const override;
        Disposable<Array> preconditioner(const Array& r, Real s) const override;

        Disposable<std::vector<SparseMatrix> > toMatrixDecomp() const override;

        //! True if all off-diagonals of mapT_ have the correct sign
        //! for the M-matrix property of (I + dt*mapT_).
        /*! The check is approximate on nonuniform grids because it
            reconstructs the expected stencil values from the standard
            three-point nonuniform finite-difference formulas rather
            than reading back from the protected coefficient arrays.  */
        bool mMatrixSatisfied() const;

        //! Number of interior nodes where the M-matrix sign pattern
        //! is violated (should be 0 when fitting is active).
        Size mMatrixViolationCount() const;

      private:
        const ext::shared_ptr<FdmMesher> mesher_;
        const ext::shared_ptr<YieldTermStructure> rTS_, qTS_;
        const ext::shared_ptr<BlackVolTermStructure> volTS_;
        const ext::shared_ptr<LocalVolTermStructure> localVol_;
        const Array x_;
        const FirstDerivativeOp dxMap_;
        const TripleBandLinearOp dxxMap_;
        TripleBandLinearOp mapT_;
        const Real strike_;
        const Real illegalLocalVolOverwrite_;
        const Size direction_;
        const ext::shared_ptr<FdmQuantoHelper> quantoHelper_;

        mutable bool mMatrixOk_;
        mutable Size mMatrixViolCount_;

        //! Exponential fitting factor: theta*coth(theta).
        /*! Returns 1 + theta^2/3 for |theta| < 1e-8 to avoid 0/0. */
        static Real fittingFactor(Real theta);
    };
}

#endif
```

---

**File 2: `ql/methods/finitedifferences/operators/fdmfittedblackscholesop.cpp`**

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 QuantLib contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/methods/finitedifferences/operators/fdmfittedblackscholesop.hpp>
#include <ql/methods/finitedifferences/operators/secondderivativeop.hpp>
#include <ql/methods/finitedifferences/meshers/fdmmesher.hpp>
#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <cmath>

namespace QuantLib {

    FdmFittedBlackScholesOp::FdmFittedBlackScholesOp(
        const ext::shared_ptr<FdmMesher>& mesher,
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Real strike,
        bool localVol,
        Real illegalLocalVolOverwrite,
        Size direction,
        ext::shared_ptr<FdmQuantoHelper> quantoHelper)
    : mesher_(mesher),
      rTS_(process->riskFreeRate().currentLink()),
      qTS_(process->dividendYield().currentLink()),
      volTS_(process->blackVolatility().currentLink()),
      localVol_(localVol
                ? process->localVolatility().currentLink()
                : ext::shared_ptr<LocalVolTermStructure>()),
      x_(mesher->locations(direction)),
      dxMap_(direction, mesher),
      dxxMap_(SecondDerivativeOp(direction, mesher)),
      mapT_(direction, mesher),
      strike_(strike),
      illegalLocalVolOverwrite_(illegalLocalVolOverwrite),
      direction_(direction),
      quantoHelper_(std::move(quantoHelper)),
      mMatrixOk_(true),
      mMatrixViolCount_(0) {}

    Size FdmFittedBlackScholesOp::size() const { return 1; }

    Real FdmFittedBlackScholesOp::fittingFactor(Real theta) {
        if (std::fabs(theta) < 1e-8)
            return 1.0 + theta * theta / 3.0;
        return theta / std::tanh(theta);
    }

    void FdmFittedBlackScholesOp::setTime(Time t1, Time t2) {
        const Rate r = rTS_->forwardRate(t1, t2, Continuous).rate();
        const Rate q = qTS_->forwardRate(t1, t2, Continuous).rate();
        const Time t = 0.5 * (t1 + t2);

        const ext::shared_ptr<FdmLinearOpLayout> layout = mesher_->layout();
        const Size size = layout->size();

        // --- 1. Collect per-node volatilities ---
        Array v(size);
        if (localVol_) {
            for (FdmLinearOpIterator iter = layout->begin();
                 iter != layout->end(); ++iter) {
                const Size i = iter.index();
                const Real S = std::exp(x_[i]);
                try {
                    v[i] = localVol_->localVol(t, S, true);
                } catch (...) {
                    v[i] = (illegalLocalVolOverwrite_ >= 0.0)
                               ? illegalLocalVolOverwrite_
                               : 0.0;
                }
            }
        } else {
            const Volatility sigma = volTS_->blackVol(t, strike_);
            std::fill(v.begin(), v.end(), sigma);
        }

        // --- 2. Build convection / fitted-diffusion arrays ---
        //
        // QuantLib sign convention: mapT_ stores L where  u_t + L u = 0,
        // i.e.  L = -a D_xx - b D_x + r I.
        //
        // With exponential fitting the diffusion coefficient a is
        // replaced by a*rho.  We assemble via axpyb:
        //
        //   mapT_ = diag(conv) * dxMap_ + dxxMap_.mult(negFittedDiff) + diag(r)
        //
        // where conv[i] = -b_i  and  negFittedDiff[i] = -a_i * rho_i.

        Array convection(size);
        Array negFittedDiff(size);

        mMatrixOk_     = true;
        mMatrixViolCount_ = 0;

        for (FdmLinearOpIterator iter = layout->begin();
             iter != layout->end(); ++iter) {
            const Size i = iter.index();

            const Real sigma_i = v[i];
            const Real a = 0.5 * sigma_i * sigma_i;   // diffusion
            const Real b = r - q - a;                  // drift in ln(S)

            // Detect boundary nodes via the Null sentinel in dplus/dminus
            const Real dp = mesher_->dplus(iter, direction_);
            const Real dm = mesher_->dminus(iter, direction_);
            const bool isBoundary =
                (dp == Null<Real>() || dm == Null<Real>());

            if (!isBoundary) {
                const Real h_eff  = 0.5 * (dp + dm);
                const Real a_safe = std::max(a, 1e-20);
                const Real theta  = b * h_eff / (2.0 * a_safe);
                const Real rho    = fittingFactor(theta);

                convection[i]   = -b;
                negFittedDiff[i] = -(a_safe * rho);

                // --- M-matrix diagnostic ---
                //
                // Expected off-diagonals of mapT_ assuming the standard
                // three-point nonuniform stencils (Bowen & Smith 2005):
                //
                //   Dx_lower  = -dp / (dm*(dp+dm))
                //   Dxx_lower =  2  / (dm*(dp+dm))
                //   Dx_upper  =  dm / (dp*(dp+dm))
                //   Dxx_upper =  2  / (dp*(dp+dm))
                //
                // mapT_.lower = conv * Dx_lower + negFittedDiff * Dxx_lower
                // mapT_.upper = conv * Dx_upper + negFittedDiff * Dxx_upper
                //
                // For the M-matrix property of (I + dt*mapT_) we need
                // both to be <= 0.

                const Real a_fitted = a_safe * rho;
                const Real sum = dp + dm;

                const Real expected_lower =
                    ( b * dp - 2.0 * a_fitted) / (dm * sum);
                const Real expected_upper =
                    (-b * dm - 2.0 * a_fitted) / (dp * sum);

                if (expected_lower > 1e-12 || expected_upper > 1e-12) {
                    mMatrixOk_ = false;
                    ++mMatrixViolCount_;
                }
            } else {
                // Boundary nodes: use unfitted (central) coefficients.
                // These rows are overridden by boundary conditions.
                convection[i]   = -b;
                negFittedDiff[i] = -a;
            }
        }

        // --- 3. Quanto drift adjustment (if applicable) ---
        //
        // quantoAdjustment returns the correction to the convection
        // array (sign follows FdmBlackScholesOp convention).
        if (quantoHelper_) {
            convection += quantoHelper_->quantoAdjustment(v, t1, t2);
        }

        // --- 4. Assemble the operator via axpyb ---
        mapT_.axpyb(convection, dxMap_,
                     dxxMap_.mult(negFittedDiff), Array(size, r));
    }

    Disposable<Array>
    FdmFittedBlackScholesOp::apply(const Array& r) const {
        return mapT_.apply(r);
    }

    Disposable<Array>
    FdmFittedBlackScholesOp::apply_mixed(const Array& r) const {
        Array retVal(r.size(), 0.0);
        return retVal;
    }

    Disposable<Array>
    FdmFittedBlackScholesOp::apply_direction(Size direction,
                                              const Array& r) const {
        if (direction == direction_)
            return mapT_.apply(r);
        Array retVal(r.size(), 0.0);
        return retVal;
    }

    Disposable<Array>
    FdmFittedBlackScholesOp::solve_splitting(Size direction,
                                              const Array& r,
                                              Real s) const {
        if (direction == direction_)
            return mapT_.solve_splitting(r, s, 1.0);
        Array retVal(r);
        return retVal;
    }

    Disposable<Array>
    FdmFittedBlackScholesOp::preconditioner(const Array& r, Real s) const {
        return solve_splitting(direction_, r, s);
    }

    Disposable<std::vector<SparseMatrix> >
    FdmFittedBlackScholesOp::toMatrixDecomp() const {
        std::vector<SparseMatrix> retVal(1);
        retVal[0] = mapT_.toMatrix();
        return retVal;
    }

    bool FdmFittedBlackScholesOp::mMatrixSatisfied() const {
        return mMatrixOk_;
    }

    Size FdmFittedBlackScholesOp::mMatrixViolationCount() const {
        return mMatrixViolCount_;
    }
}
```

---

**File 3: `ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.hpp`**

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 QuantLib contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file fdmsinhconcentrating1dmesher.hpp
    \brief One-dimensional FDM mesher with sinh-based concentration
*/

#ifndef quantlib_fdm_sinh_concentrating_1d_mesher_hpp
#define quantlib_fdm_sinh_concentrating_1d_mesher_hpp

#include <ql/methods/finitedifferences/meshers/fdm1dmesher.hpp>
#include <vector>

namespace QuantLib {

    //! Sinh-graded 1D mesher with optional node alignment.
    /*! Builds a grid on \f$[x_{\min},x_{\max}]\f$ that concentrates
        points around \p xCenter using the mapping
        \f[
            x(\xi) = x_{\text{center}} + c\,\sinh\!\bigl(\alpha\,(\xi-\xi_0)\bigr),
            \qquad \xi_j = \frac{j}{N-1},\; j=0,\dots,N-1,
        \f]
        where \f$\xi_0\f$ and \f$c\f$ are determined from the endpoint
        constraints \f$x(0)=x_{\min}\f$, \f$x(1)=x_{\max}\f$.

        When \p alpha is zero (or near-zero) the mesher falls back to a
        uniform grid.

        If \p alignTargets is non-empty, a single global shift (the
        smallest that is less than half the minimum spacing) is applied
        so that the nearest grid node snaps to one of the targets.
        Endpoints are re-pinned after the shift.
    */
    class FdmSinhConcentrating1dMesher : public Fdm1dMesher {
      public:
        FdmSinhConcentrating1dMesher(
            Real xMin, Real xMax, Size size,
            Real xCenter,
            Real alpha = 3.0,
            const std::vector<Real>& alignTargets = std::vector<Real>());

        Real alpha() const;
        Real xCenter() const;

      private:
        Real alpha_, xCenter_;
    };
}

#endif
```

---

**File 4: `ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.cpp`**

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 QuantLib contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/methods/finitedifferences/meshers/fdmsinhconcentrating1dmesher.hpp>
#include <ql/errors.hpp>
#include <ql/utilities/null.hpp>
#include <cmath>
#include <algorithm>

namespace QuantLib {

    FdmSinhConcentrating1dMesher::FdmSinhConcentrating1dMesher(
        Real xMin, Real xMax, Size size,
        Real xCenter, Real alpha,
        const std::vector<Real>& alignTargets)
    : Fdm1dMesher(size), alpha_(alpha), xCenter_(xCenter)
    {
        QL_REQUIRE(size >= 2,
                   "FdmSinhConcentrating1dMesher: size (" << size
                   << ") must be at least 2");
        QL_REQUIRE(xMax > xMin,
                   "FdmSinhConcentrating1dMesher: xMax (" << xMax
                   << ") must be greater than xMin (" << xMin << ")");
        QL_REQUIRE(alpha >= 0.0,
                   "FdmSinhConcentrating1dMesher: alpha (" << alpha
                   << ") must be non-negative");

        // ----------------------------------------------------------
        // 1.  Build the base grid
        // ----------------------------------------------------------

        if (alpha < 1e-10) {
            // Uniform fallback
            const Real dx = (xMax - xMin) / (size - 1);
            for (Size i = 0; i < size; ++i)
                locations_[i] = xMin + i * dx;
        } else {
            QL_REQUIRE(xCenter >= xMin && xCenter <= xMax,
                       "FdmSinhConcentrating1dMesher: xCenter ("
                       << xCenter << ") must lie within [xMin, xMax] = ["
                       << xMin << ", " << xMax << "]");

            // Find xi0 in (0,1) via bisection on the endpoint ratio
            //
            //   f(xi0) = (xMin - xCenter) * sinh(alpha*(1 - xi0))
            //          - (xMax - xCenter) * sinh(-alpha*xi0)
            //
            // f(0+) < 0 and f(1-) > 0 when xMin < xCenter < xMax.

            const Real Lm = xMin - xCenter;   // typically < 0
            const Real Rm = xMax - xCenter;    // typically > 0

            Real xi0;
            if (std::fabs(Lm + Rm) < 1e-12 * (xMax - xMin)) {
                // Symmetric case: xCenter == midpoint => xi0 = 0.5
                xi0 = 0.5;
            } else {
                Real lo = 1e-8, hi = 1.0 - 1e-8;
                for (Size iter = 0; iter < 200; ++iter) {
                    const Real mid = 0.5 * (lo + hi);
                    const Real fMid =
                        Lm * std::sinh(alpha * (1.0 - mid))
                      - Rm * std::sinh(-alpha * mid);
                    if (fMid < 0.0)
                        lo = mid;
                    else
                        hi = mid;
                }
                xi0 = 0.5 * (lo + hi);
            }

            // Compute c from the right-endpoint constraint
            //   xMax = xCenter + c * sinh(alpha*(1 - xi0))
            const Real sinhRight = std::sinh(alpha * (1.0 - xi0));
            QL_REQUIRE(std::fabs(sinhRight) > 1e-30,
                       "FdmSinhConcentrating1dMesher: degenerate sinh "
                       "at right endpoint");
            const Real c = Rm / sinhRight;

            for (Size i = 0; i < size; ++i) {
                const Real xi =
                    static_cast<Real>(i) / static_cast<Real>(size - 1);
                locations_[i] = xCenter + c * std::sinh(alpha * (xi - xi0));
            }

            // Pin endpoints for numerical safety
            locations_.front() = xMin;
            locations_.back()  = xMax;
        }

        // ----------------------------------------------------------
        // 2.  Optional alignment shift
        // ----------------------------------------------------------

        if (!alignTargets.empty()) {
            // Minimum spacing
            Real minSpacing = locations_[1] - locations_[0];
            for (Size i = 1; i + 1 < size; ++i)
                minSpacing = std::min(minSpacing,
                                      locations_[i + 1] - locations_[i]);

            // Choose the smallest valid shift across all targets
            Real bestShift    = Null<Real>();   // sentinel: very large
            Real bestAbsShift = Null<Real>();

            for (std::vector<Real>::const_iterator t = alignTargets.begin();
                 t != alignTargets.end(); ++t) {
                if (*t < xMin || *t > xMax)
                    continue;

                // Find nearest node
                Size nearest = 0;
                Real nearDist = std::fabs(locations_[0] - *t);
                for (Size i = 1; i < size; ++i) {
                    const Real d = std::fabs(locations_[i] - *t);
                    if (d < nearDist) {
                        nearDist = d;
                        nearest  = i;
                    }
                }

                const Real shift    = *t - locations_[nearest];
                const Real absShift = std::fabs(shift);

                if (absShift < 0.5 * minSpacing && absShift < bestAbsShift) {
                    bestShift    = shift;
                    bestAbsShift = absShift;
                }
            }

            if (bestAbsShift < Null<Real>()) {
                for (Size i = 0; i < size; ++i)
                    locations_[i] += bestShift;

                // Re-pin endpoints
                locations_.front() = xMin;
                locations_.back()  = xMax;
            }
        }

        // ----------------------------------------------------------
        // 3.  Populate dplus_ / dminus_
        // ----------------------------------------------------------

        for (Size i = 0; i + 1 < size; ++i)
            dplus_[i] = dminus_[i + 1] = locations_[i + 1] - locations_[i];

        dplus_.back()  = Null<Real>();
        dminus_.front() = Null<Real>();
    }

    Real FdmSinhConcentrating1dMesher::alpha() const { return alpha_; }
    Real FdmSinhConcentrating1dMesher::xCenter() const { return xCenter_; }
}
```

---

**File 5: `ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp`**

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 QuantLib contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file fdmdiagnostics.hpp
    \brief Runtime health diagnostics for finite-difference solutions
*/

#ifndef quantlib_fdm_diagnostics_hpp
#define quantlib_fdm_diagnostics_hpp

#include <ql/math/array.hpp>
#include <ql/types.hpp>

namespace QuantLib {

    //! Aggregated health report for an FDM solution snapshot.
    struct FdmDiagnosticsReport {
        Real minValue;                //!< Smallest element
        Size negativeCount;           //!< Number of strictly negative elements
        Real oscillationScore;        //!< Normalized sign-change count in [0,1]
        Size mMatrixViolationCount;   //!< Violation count from operator diagnostic
        Size nanCount;                //!< Number of NaN elements

        FdmDiagnosticsReport()
        : minValue(QL_MAX_REAL),
          negativeCount(0),
          oscillationScore(0.0),
          mMatrixViolationCount(0),
          nanCount(0) {}
    };

    //! Lightweight utility for checking FDM solution health.
    /*! All static methods are allocation-free and thread-safe.
        The instance method checkSolution respects the configured
        \c Level to avoid unnecessary work.
    */
    class FdmDiagnostics {
      public:
        //! Determines how much work checkSolution does.
        enum Level {
            Off,    //!< Return an empty report immediately.
            Light,  //!< Scan for min, negatives, NaNs only.
            Full    //!< Light checks plus oscillation score.
        };

        explicit FdmDiagnostics(Level level = Off);

        //! Inspect the solution array according to the configured level.
        FdmDiagnosticsReport checkSolution(const Array& u) const;

        //! Normalized oscillation score in [0,1].
        /*! Counts sign changes in the consecutive-difference sequence
            \f$\Delta u_j = u_{j+1} - u_j\f$, ignoring differences
            smaller than \f$10^{-15}\f$ in magnitude, and divides by
            \f$\max(1,\,N-2)\f$ where \f$N\f$ is the array size. */
        static Real oscillationScore(const Array& u);

        //! Merge two reports by taking the worst case of each field.
        static FdmDiagnosticsReport merge(
            const FdmDiagnosticsReport& a,
            const FdmDiagnosticsReport& b);

        Level level() const;

      private:
        Level level_;
    };
}

#endif
```

---

**File 6: `ql/methods/finitedifferences/utilities/fdmdiagnostics.cpp`**

```cpp
/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 QuantLib contributors

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/methods/finitedifferences/utilities/fdmdiagnostics.hpp>
#include <cmath>
#include <algorithm>

namespace QuantLib {

    FdmDiagnostics::FdmDiagnostics(Level level)
    : level_(level) {}

    FdmDiagnostics::Level FdmDiagnostics::level() const {
        return level_;
    }

    FdmDiagnosticsReport
    FdmDiagnostics::checkSolution(const Array& u) const {

        FdmDiagnosticsReport r;

        if (level_ == Off)
            return r;

        // --- Light: min, negatives, NaNs ---
        r.minValue = QL_MAX_REAL;
        for (Size i = 0; i < u.size(); ++i) {
            if (std::isnan(u[i])) {
                ++r.nanCount;
            } else {
                if (u[i] < r.minValue)
                    r.minValue = u[i];
                if (u[i] < 0.0)
                    ++r.negativeCount;
            }
        }

        // --- Full: add oscillation score ---
        if (level_ == Full)
            r.oscillationScore = oscillationScore(u);

        return r;
    }

    Real FdmDiagnostics::oscillationScore(const Array& u) {
        if (u.size() < 3)
            return 0.0;

        const Real eps = 1e-15;

        Size changes = 0;
        int  prevSign = 0;

        for (Size j = 0; j + 1 < u.size(); ++j) {
            const Real du = u[j + 1] - u[j];

            int s = 0;
            if (du > eps)
                s = 1;
            else if (du < -eps)
                s = -1;

            if (s != 0 && prevSign != 0 && s != prevSign)
                ++changes;

            if (s != 0)
                prevSign = s;
        }

        return static_cast<Real>(changes)
             / static_cast<Real>(std::max<Size>(1, u.size() - 2));
    }

    FdmDiagnosticsReport FdmDiagnostics::merge(
        const FdmDiagnosticsReport& a,
        const FdmDiagnosticsReport& b)
    {
        FdmDiagnosticsReport r;
        r.minValue              = std::min(a.minValue, b.minValue);
        r.negativeCount         = std::max(a.negativeCount, b.negativeCount);
        r.oscillationScore      = std::max(a.oscillationScore,
                                           b.oscillationScore);
        r.mMatrixViolationCount = std::max(a.mMatrixViolationCount,
                                           b.mMatrixViolationCount);
        r.nanCount              = std::max(a.nanCount, b.nanCount);
        return r;
    }
}
```

---

**Quality checklist verification:**

- ☑ `FdmFittedBlackScholesOp` implements all 7 pure virtual methods of `FdmLinearOpComposite` (`size`, `setTime`, `apply`, `apply_mixed`, `apply_direction`, `solve_splitting`, `preconditioner`) plus `toMatrixDecomp`
- ☑ `setTime()` uses `axpyb()` — no access to protected `lower_`/`diag_`/`upper_` members
- ☑ `fittingFactor` returns `1.0 + θ²/3` for `|θ| < 1e-8`, `θ/tanh(θ)` otherwise
- ☑ Fitted operator computes `h_eff = 0.5*(dplus + dminus)` from the mesher, consistent with derivative operators on nonuniform meshes
- ☑ Sinh mesher populates `locations_`, `dplus_`, `dminus_` following `Uniform1dMesher` pattern (including `Null<Real>()` sentinels at boundaries)
- ☑ Sinh mesher handles `alpha ≈ 0` by falling back to uniform grid
- ☑ `oscillationScore` is O(N) with zero heap allocations (operates on `const Array&`)
- ☑ All `ext::shared_ptr`, no `std::shared_ptr`
- ☑ `QL_REQUIRE` for all preconditions
- ☑ Include guards follow `quantlib_<path_underscored>_hpp` convention
- ☑ No missing source file content needed — all implementations derive from provided headers and round0 findings
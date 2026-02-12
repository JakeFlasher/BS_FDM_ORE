/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2024 QuantLib Contributors

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

#include <ql/methods/finitedifferences/fittedbsmoperator.hpp>
#include <ql/math/transformedgrid.hpp>
#include <ql/methods/finitedifferences/pdebsm.hpp>
#include <cmath>

namespace QuantLib {

    namespace {

        // Evaluates  f(x) = x coth(x)  with the convention f(0) = 1.
        //
        // For |x| < 1e-5 a fourth-order Taylor expansion is used to
        // avoid cancellation in the ratio x / tanh(x).  For larger
        // arguments the standard library functions are accurate.
        //
        // Key properties used by the fitting scheme:
        //   f(x) >= |x|  for all x        (=> M-matrix guarantee)
        //   f(x) -> 1    as x -> 0        (=> recovers centred scheme)
        //   f(x) -> |x|  as |x| -> inf    (=> upwind in convection limit)
        Real xcothx(Real x) {
            if (std::fabs(x) < 1e-5) {
                Real x2 = x * x;
                return 1.0 + x2 / 3.0 - x2 * x2 / 45.0;
            }
            return x / std::tanh(x);
        }

    }

    FittedBSMOperator::FittedBSMOperator(Size size, Real dx,
                                          Rate r, Rate q,
                                          Volatility sigma)
    : TridiagonalOperator(size) {
        Real sigma2 = sigma * sigma;
        Real nu = r - q - sigma2 / 2;

        // fitted_sigma2 replaces sigma2 in the standard BSM stencil
        //   pd = -(sigma2/dx - nu) / (2 dx)
        //   pu = -(sigma2/dx + nu) / (2 dx)
        //   pm =  sigma2 / dx^2 + r
        //
        // It equals  sigma2 * xcothx(nu*dx / sigma2) ,  so that
        //   * nu -> 0  =>  fitted_sigma2 -> sigma2   (centred scheme)
        //   * sigma -> 0, nu > 0  =>  fitted_sigma2 -> nu*dx  (upwind)
        Real fitted_sigma2;
        if (sigma2 > QL_EPSILON) {
            Real rho_arg = nu * dx / sigma2;
            fitted_sigma2 = sigma2 * xcothx(rho_arg);
        } else if (std::fabs(nu) > QL_EPSILON) {
            // degenerate: pure convection
            fitted_sigma2 = std::fabs(nu) * dx;
        } else {
            // degenerate: pure reaction
            fitted_sigma2 = 0.0;
        }

        Real pd = -(fitted_sigma2 / dx - nu) / (2 * dx);
        Real pu = -(fitted_sigma2 / dx + nu) / (2 * dx);
        Real pm =  fitted_sigma2 / (dx * dx) + r;
        setMidRows(pd, pm, pu);
    }

    FittedBSMOperator::FittedBSMOperator(const Array& grid,
                                          Rate r, Rate q,
                                          Volatility sigma)
    : TridiagonalOperator(grid.size()) {
        PdeBSM::grid_type logGrid(grid);
        Real sigma2 = sigma * sigma;
        Real nu = r - q - sigma2 / 2;

        for (Size i = 1; i < logGrid.size() - 1; ++i) {
            Real dxm = logGrid.dxm(i);
            Real dxp = logGrid.dxp(i);
            Real dx  = logGrid.dx(i);   // = dxm + dxp

            // Local effective spacing for the fitting factor.
            // Using max(dxm,dxp) guarantees the M-matrix property
            // even on non-uniform meshes: fitted_sigma2 >= |nu|*h_eff
            // ensures both off-diagonal entries are non-positive.
            Real h_eff = std::max(dxm, dxp);

            Real fitted_sigma2;
            if (sigma2 > QL_EPSILON) {
                Real rho_arg = nu * h_eff / sigma2;
                fitted_sigma2 = sigma2 * xcothx(rho_arg);
            } else if (std::fabs(nu) > QL_EPSILON) {
                fitted_sigma2 = std::fabs(nu) * h_eff;
            } else {
                fitted_sigma2 = 0.0;
            }

            Real pd = -(fitted_sigma2 / dxm - nu) / dx;
            Real pu = -(fitted_sigma2 / dxp + nu) / dx;
            Real pm =  fitted_sigma2 / (dxm * dxp) + r;
            setMidRow(i, pd, pm, pu);
        }
    }

}

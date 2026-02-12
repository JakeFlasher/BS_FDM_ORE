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

/*! \file fittedbsmoperator.hpp
    \brief exponentially fitted differential operator for
           Black-Scholes-Merton equation
*/

#ifndef quantlib_fitted_bsm_operator_hpp
#define quantlib_fitted_bsm_operator_hpp

#include <ql/methods/finitedifferences/tridiagonaloperator.hpp>

namespace QuantLib {

    //! Exponentially fitted Black-Scholes-Merton differential operator
    /*! Replaces the centered second-derivative discretization of
        BSMOperator with a fitted diffusion coefficient
        \f[
          \rho = \frac{\nu h}{2}\coth\!\Bigl(\frac{\nu h}{\sigma^2}\Bigr)
        \f]
        where \f$ \nu = r - q - \tfrac12\sigma^2 \f$ is the log-space
        drift and \f$ h \f$ is the local grid spacing.  The resulting
        tridiagonal system is an M-matrix for any mesh Péclet number,
        eliminating the spurious oscillations that standard
        centred-in-space discretizations can produce when the diffusion
        coefficient is small relative to convection.

        In the limit \f$ \sigma\to 0 \f$ the scheme degrades
        gracefully to an implicit upwind discretization; in the limit
        \f$ \nu\to 0 \f$ it reduces identically to BSMOperator.

        Reference: D. J. Duffy, "A Critique of the Crank Nicolson
        Scheme," Wilmott, July 2004.

        \ingroup findiff
    */
    class FittedBSMOperator : public TridiagonalOperator {
      public:
        FittedBSMOperator() = default;
        /*! Uniform-grid constructor.
            \param size   number of grid points
            \param dx     uniform log-space grid step
            \param r      risk-free rate
            \param q      dividend yield
            \param sigma  Black-Scholes volatility
        */
        FittedBSMOperator(Size size, Real dx,
                          Rate r, Rate q, Volatility sigma);
        /*! Non-uniform-grid constructor (grid contains S-values). */
        FittedBSMOperator(const Array& grid,
                          Rate r, Rate q, Volatility sigma);
    };

}


#endif

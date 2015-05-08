/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Ferdinando Ametrano
 Copyright (C) 2015 Paolo Mazzocchi

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

#ifndef quantlib_polynomial_math_function_hpp
#define quantlib_polynomial_math_function_hpp

#include <ql/math/matrix.hpp>

#include <vector>

namespace QuantLib {
    
    //! %Cubic functional form
    /*! \f[ f(t) = c_0 + c_1*t + c_2*t^2 + c_3*t^3 \f]*/
    class PolynomialFunction : public std::unary_function<Time, Real> {

    public:
        PolynomialFunction(const std::vector<Real>& coeff);

        //! function value at time t: \f[ f(t) \f]
        Real operator()(Time t) const;

        /*! first derivative of the function at time t
        \f[ f'(t)dt = b + 2*c*t + 3*d*t^2 \f] */
        Real derivative(Time t) const;

        /*! indefinite integral of the function at time t
        \f[ \int f(t)dt = a*t + b*x^2/2 + c*t^3/3 + d*t^4/4 \f] */
        Real primitive(Time t) const;

        /*! definite integral of the function between t1 and t2
        \f[ \int_{t1}^{t2} f(t)dt \f] */
        Real definiteIntegral(Time t1,
                              Time t2) const;

        /*! Inspectors */
        const std::vector<Real>& coefficients() { return c_; }
        const std::vector<Real>& derivativeCoefficients() { return derC_; }
        const std::vector<Real>& primitiveCoefficients() { return prC_; }
        /*! coefficients of definite integral on a rolling window of tau, with tau = t2-t */
        std::vector<Real> definiteIntegralCoefficients(Time t,
                                                       Time t2) const;
        std::vector<Real> definiteDerivativeCoefficients(Time t,
                                                         Time t2) const;


    private:
        Size order_;
        std::vector<Real> c_;
        mutable Matrix eqs_;
        std::vector<Real> derC_, prC_;
        void initializeEqs_(Time t,Time t2) const;
    };

}

#endif

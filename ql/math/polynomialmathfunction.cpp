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

#include <ql/math/polynomialmathfunction.hpp>
#include <ql/math/comparison.hpp>
#include <ql/math/tartaglia.hpp>
#include <ql/math/matrix.hpp>

namespace QuantLib {

    PolynomialFunction::PolynomialFunction(const std::vector<Real>& coeff)
    : c_(coeff) {
        QL_REQUIRE(!c_.empty(), "empty vector");

        order_ = c_.size();
        derC_ = std::vector<Real>(order_ - 1), prC_ = std::vector<Real>(order_);
        Size i;
        for (i=0; i<order_-1; ++i) {
            prC_[i] = c_[i]/(i+1);
            derC_[i] = c_[i+1]*(i+1);
        }
        prC_[i] = c_[i]/(i + 1);
    }

    Real PolynomialFunction::definiteIntegral(Time t1,
                                              Time t2) const {
        return primitive(t2)-primitive(t1);
    }

    std::vector<Real>
    PolynomialFunction::definiteIntegralCoefficients(Time t,
                                                     Time t2) const {
        Time dt = t2 - t;
        std::vector<Real> result(order_, 0);
        Real coef, tau;
        for (Size i=0; i<order_; ++i) {
            tau = 1.0;
            for (Size j=i; j<order_; ++j) {
                 coef = Tartaglia::get(j+1)[i];
                 tau *= dt;
                 result[i] += prC_[j]*coef*tau; 
            }
        }
        return result;
    }

    std::vector<Real>
    PolynomialFunction::definiteDerivativeCoefficients(Time t,
                                                       Time t2) const {
        Time dt = t2 - t;
        Matrix eqs(order_, order_);
        Array k(order_);
        Real tau;
        for (Size i=0; i<order_; ++i) {
            k[i] = c_[i];
            tau= 1.0;
            for (Size j=0; j<order_; ++j) {
                if (j<i)
                    eqs[i][j]=0;
                else {
                    tau *= dt;
                    eqs[i][j] = (tau * Tartaglia::get(j + 1)[i])/(j + 1);
                }
            }
        }
        Array coeff = inverse(eqs) * k;
        std::vector<Real> result(coeff.begin(), coeff.end());
        return result;
    }

}

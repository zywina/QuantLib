/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007, 2015 Ferdinando Ametrano
 Copyright (C) 2006 Cristina Duminuco
 Copyright (C) 2005, 2006 Klaus Spanderen
 Copyright (C) 2007 Giorgio Facchinetti
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

#include <ql/math/abcdmathfunction.hpp>
#include <ql/math/comparison.hpp>

namespace QuantLib {

    AbcdMathFunction::AbcdMathFunction(Real a, Real b, Real c, Real d)
    : a_(a), b_(b), c_(c), d_(d) {
        validateAbcdParameters(a, b, c, d);
        da_ = b_ - c_*a_;
        db_ = -c_*b_;

        pa_ = -(a_ + b_/c_)/c_;
        pb_ = -b_/c_;
        K_ = 0.0;

        dibc_ = b_/c_;
        diacplusbcc_ = a_/c_ + dibc_/c_;
    }

    Time AbcdMathFunction::maximumLocation() const {
        if (b_<=0.0)
            return 0.0;
        else if (c_==0.0)
            return QL_MAX_REAL;
        else {
            Real temp = (b_-c_*a_)/(c_*b_);
            if (temp>0)
                return temp;
            else
                return 0.0;
        }
    }

    Real AbcdMathFunction::definiteIntegral(Time t1, Time t2) const {
        QL_REQUIRE(t2>=t1, "final time (" << t2 << ") must be greater "
                           "than initial time (" << t1 << ")");

        Time dt = t2 - t1;
        Real expcdt = std::exp(-c_*dt);
        Real dia = - (diacplusbcc_ + dibc_*dt)*expcdt + diacplusbcc_;
        Real dib = dibc_ * (1.0 - expcdt);
        return (dia + dib*t1)*std::exp(-c_*t1) + d_*dt;
    }

    PolynomialFunction::PolynomialFunction(const std::vector<Real>& coeff)
    : c_(coeff) {
        QL_REQUIRE(!c_.empty(), "empty vector");

        order_ = c_.size();
        derC_ = std::vector<Real>(order_ - 1), prC_ = std::vector<Real>(order_);
        Size i;
        for (i = 0; i<order_-1; ++i) {
            prC_[i] = c_[i]/(i+1);
            derC_[i] = c_[i+1]*(i+1);
        }
        prC_[i] = c_[i] / (i + 1);
    }

    Real PolynomialFunction::definiteIntegral(Time t1, Time t2) const {
        QL_REQUIRE(t2 >= t1, "final time (" << t2 << ") must be greater " 
                             "than initial time (" << t1 << ")");

        return primitive(t2)-primitive(t1);
    }

    const std::vector<Real>& PolynomialFunction::
                                definitiveIntegralCoefficients(Time t1,
                                                               Time t2) const {
        QL_REQUIRE(t2 >= t1, "final time (" << t2 << ") must be greater "
                             "than initial time (" << t1 << ")");
        //order_ = c_.size();
        Time dt = t2 - t1;
        std::vector<Real> diCoef(order_, 0);
        Size i,j;
        for (i = 0; i<order_; ++i) {
            for (j = i; j<order_; ++j){
                diCoef[i] += prC_[j]*1.0*std::pow(dt,i-j); 
            }
        }
        return diCoef;
    }

}

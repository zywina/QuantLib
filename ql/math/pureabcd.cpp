/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007, 2015 Ferdinando Ametrano
 Copyright (C) 2006 Cristina Duminuco
 Copyright (C) 2005, 2006 Klaus Spanderen
 Copyright (C) 2007 Giorgio Facchinetti

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

#include <ql/math/pureabcd.hpp>
#include <ql/math/comparison.hpp>

namespace QuantLib {

    PureAbcdFunction::PureAbcdFunction(Real a, Real b, Real c, Real d)
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

    Real PureAbcdFunction::maximumLocation() const {
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

    Real PureAbcdFunction::definiteIntegral(Time t1, Time t2) const {
        QL_REQUIRE(t2>=t1, "final time (" << t2 << ") must be greater "
                           "than intial time (" << t1 << ")");

        Time dt = t2 - t1;
        Real expcdt = std::exp(-c_*dt);
        Real dia = - (diacplusbcc_ + dibc_*dt)*expcdt + diacplusbcc_;
        Real dib = dibc_ * (1.0 - expcdt);
        return (dia + dib*t1)*std::exp(-c_*t1) + d_*dt;
    }
}

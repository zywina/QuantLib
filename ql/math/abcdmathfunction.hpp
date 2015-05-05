/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007, 2015 Ferdinando Ametrano
 Copyright (C) 2006 Cristina Duminuco
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

#ifndef quantlib_abcd_math_function_hpp
#define quantlib_abcd_math_function_hpp

#include <ql/types.hpp>
#include <ql/errors.hpp>

#include <vector>

namespace QuantLib {
    
    inline void validateAbcdParameters(Real a,
                                       Real, // no condition on b
                                       Real c,
                                       Real d) {
        QL_REQUIRE(a+d>=0,
                   "a (" << a << ") + d (" << d << ") must be non negative");
        QL_REQUIRE(c>=0,
                   "c (" << c << ") must be non negative");
        QL_REQUIRE(d>=0,
                   "d (" << d << ") must be non negative");
    }

    //! %Abcd functional form
    /*! \f[ f(t) = [ a + b*t ] e^{-c*t} + d \f]
        following Rebonato's notation. */
    class AbcdMathFunction : public std::unary_function<Time, Real> {

      public:
        AbcdMathFunction(Real a = 0.002,
                         Real b = 0.001, 
                         Real c = 0.16,
                         Real d = 0.0005);

        //! function value at time t: \f[ f(t) \f]
        Real operator()(Time t) const;

        //! time at which the function reaches maximum (if any)
        Time maximumLocation() const;

        //! maximum value of the function
        Real maximumValue() const;

        //! function value at time +inf: \f[ f(\inf) \f]
        Real longTermValue() const { return d_; }

        /*! first derivative of the function at time t
            \f[ f'(t)dt = [ (b-c*a) + (-c*b)*t) ] e^{-c*t} \f] */
        Real derivative(Time t) const;
        
        /*! indefinite integral of the function at time t
            \f[ \int f(t)dt = [ (-a/c-b/c^2) + (-b/c)*t ] e^{-c*t} + d*t \f] */
        Real primitive(Time t) const;
        
        /*! definite integral of the function between t1 and t2
            \f[ \int_{t1}^{t2} f(t)dt \f] */
        Real definiteIntegral(Time t1, Time t2) const;

        /*! Inspectors */
        Real a() const { return a_; }
        Real b() const { return b_; }
        Real c() const { return c_; }
        Real d() const { return d_; }

        Real derivativeA() const { return da_; }
        Real derivativeB() const { return db_; }
        Real derivativeC() const { return c_; }
        Real derivativeD() const { return 0.0; }

        Real definiteIntegralA(Time t1, Time t2) const;
        Real definiteIntegralB(Time t1, Time t2) const;
        Real definiteIntegralC(Time,    Time   ) const { return c_; }
        Real definiteIntegralD(Time t1, Time t2) const;

      protected:
        Real a_, b_, c_, d_;
      private:
        Real da_, db_;
        Real pa_, pb_, K_;

        Real dibc_, diacplusbcc_;
    };

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
        Real definiteIntegral(Time t1, Time t2) const;

        /*! coefficients of definite integral on a rolling window of tau, with tau = t2-t1 */
        std::vector<Real> definiteIntegralCoefficients(Time t1,
                                                       Time t2) const;

        /*! Inspectors */
        const std::vector<Real>& coefficients() { return c_; }
        const std::vector<Real>& derivativeCoefficients() { return derC_; }
        const std::vector<Real>& primitiveCoefficients() { return prC_; }

    protected:
        std::vector<Real> c_;
        Size order_;
    private:
        std::vector<Real> derC_, prC_;
    };

    // inline AbcdMathFunction
    inline Real AbcdMathFunction::operator()(Time t) const {
        //return (a_ + b_*t)*std::exp(-c_*t) + d_;
        return t<0 ? 0.0 : (a_ + b_*t)*std::exp(-c_*t) + d_;
    }

    inline Real AbcdMathFunction::derivative(Time t) const {
        //return (da_ + db_*t)*std::exp(-c_*t);
        return t<0 ? 0.0 : (da_ + db_*t)*std::exp(-c_*t);
    }

    inline Real AbcdMathFunction::primitive(Time t) const {
        //return (pa_ + pb_*t)*std::exp(-c_*t) + d_*t + K_;
        return t<0 ? 0.0 : (pa_ + pb_*t)*std::exp(-c_*t) + d_*t + K_;
    }

    inline Real AbcdMathFunction::maximumValue() const {
        if (c_==0.0)
            return QL_MAX_REAL;

        return this->operator()(maximumLocation());
    }

    inline Real AbcdMathFunction::definiteIntegralA(Time t1, Time t2) const {
        Time dt = t2 - t1;
        Real expcdt = std::exp(-c_*dt);
        return diacplusbcc_ - (diacplusbcc_ + dibc_*dt)*expcdt;
    }

    inline Real AbcdMathFunction::definiteIntegralB(Time t1, Time t2) const {
        Time dt = t2 - t1;
        Real expcdt = std::exp(-c_*dt);
        return dibc_ * (1.0 - expcdt);
    }

    inline Real AbcdMathFunction::definiteIntegralD(Time t1, Time t2) const {
        Time dt = t2 - t1;
        return d_*dt;
    }

    // inline PolynomialFunction
    inline Real PolynomialFunction::primitive(Time t) const {
        Real result = 0.0, tPower = t;
        for (Size i = 0; i<order_; ++i) {
            result += prC_[i]*tPower;
            tPower *= t;
        }
        return result;
    }

    inline Real PolynomialFunction::operator()(Time t) const {
        Real result = 0.0, tPower = 1.0;
        for (Size i = 0; i<order_; ++i) {
            result += c_[i] * tPower;
            tPower *= t;
        }
        return result;
    }

    inline Real PolynomialFunction::derivative(Time t) const {
        Real result = 0.0, tPower = 1.0;
        for (Size i = 0; i<order_-1; ++i) {
            result += derC_[i] * tPower;
            tPower *= t;
        }
        return result;
    }

}

#endif

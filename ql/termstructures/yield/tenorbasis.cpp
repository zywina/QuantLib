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

#include <ql/termstructures/yield/tenorbasis.hpp>
#include <ql/indexes/iborindex.hpp>

using boost::shared_ptr;
using std::unary_function;

namespace QuantLib {

    TenorBasis::TenorBasis(Date settlementDate,
                           shared_ptr<IborIndex> iborIndex,
                           const Handle<YieldTermStructure>& baseCurve)
    : settlementDate_(settlementDate), index_(iborIndex),
      baseCurve_(baseCurve) {
        // TODO: check iborIndex and contBasis pointers
        // TODO: check iborIndex dayCounter
        dc_ = index_->dayCounter();
        bdc_ = index_->businessDayConvention();
        eom_ = index_->endOfMonth();
        cal_ = index_->fixingCalendar();
        tenor_ = index_->tenor();
        Date endDate = cal_.advance(settlementDate, tenor_, bdc_, eom_);
        dt_ = dc_.yearFraction(settlementDate, endDate);
        time2date_ = (endDate - settlementDate)/dt_;
    }

    AbcdTenorBasis::AbcdTenorBasis(Date settlementDate,
                                   shared_ptr<IborIndex> iborIndex,
                                   const Handle<YieldTermStructure>& baseCurve,
                                   shared_ptr<AbcdMathFunction> abcd)
    : TenorBasis(settlementDate, iborIndex, baseCurve), basis_(abcd) {

        Real a = basis_->definiteDerivativeA(0.0, dt_) * dt_;
        Real b = basis_->definiteDerivativeB(0.0, dt_) * dt_;
        Real c = basis_->definiteDerivativeC(0.0, dt_);
        Real d = basis_->definiteDerivativeD(0.0, dt_) * dt_;
        instBasis_ = shared_ptr<AbcdMathFunction>(new
                                                AbcdMathFunction(a, b, c, d));
}

    PolynomialTenorBasis::PolynomialTenorBasis(
                                    Date settlementDate,
                                    shared_ptr<IborIndex> iborIndex,
                                    const Handle<YieldTermStructure>& baseCurve,
                                    shared_ptr<PolynomialFunction> p)
    : TenorBasis(settlementDate, iborIndex, baseCurve), basis_(p) {}

    IntegralTenorBasis::IntegralTenorBasis(
                                Date settlementDate,
                                shared_ptr<IborIndex> iborIndex,
                                const Handle<YieldTermStructure>& baseCurve,
                                shared_ptr<unary_function<Real, Real> > b)
    : TenorBasis(settlementDate, iborIndex, baseCurve), instBasis_(b) {}


    Real IntegralTenorBasis::value(Time t) const {
        Date d = dateFromTime(t);
        Date d2 = cal_.advance(d, tenor_, bdc_, eom_);
        Time t2 = timeFromSettlementDate(d2);
        return value(t, t2);
    }

    Real IntegralTenorBasis::value(Time t1, 
                                   Time t2) const {
        Real bigDelta = integrate(t1, t2);
        Date d1 = dateFromTime(t1);
        Date d2 = dateFromTime(t2);
        Time dt = t2 - t1;

        //DiscountFactor disc1 = baseCurve_->discount(t1);
        //DiscountFactor disc2 = baseCurve_->discount(t2);
        //Real discRatio = disc1 / disc2;
        //Rate fwd = (discRatio*std::exp(bigDelta) - 1.0) / dt;
        //Rate baseCurveFwd = (discRatio - 1.0) / dt;

        Rate baseCurveFwd = baseCurve_->forwardRate(d1, d2, dc_, Simple, Annual, 0);
        Rate fwd = ((1.0 + baseCurveFwd*dt)*std::exp(bigDelta) - 1.0) / dt;
        return fwd - baseCurveFwd;
    }

    Real IntegralTenorBasis::integrate(Date d) const {
        Time t = dc_.yearFraction(settlementDate_, d);
        Date d2 = cal_.advance(d, tenor_, bdc_, eom_);
        Time t2 = dc_.yearFraction(settlementDate_, d2);
        Real result = this->integrate(t, t2);
        return result;
    }

    AbcdIntegralTenorBasis::AbcdIntegralTenorBasis(
                                Date settlementDate,
                                shared_ptr<IborIndex> iborIndex,
                                const Handle<YieldTermStructure>& baseCurve,
                                shared_ptr<AbcdMathFunction> abcd)
    : IntegralTenorBasis(settlementDate, iborIndex, baseCurve, abcd),
      instBasis_(abcd) {

        Real a = instBasis_->definiteIntegralA(0.0, dt_)/dt_;
        Real b = instBasis_->definiteIntegralB(0.0, dt_)/dt_;
        Real c = instBasis_->definiteIntegralC(0.0, dt_);
        Real d = instBasis_->definiteIntegralD(0.0, dt_)/dt_;
        basis_ = shared_ptr<AbcdMathFunction>(new 
                                                AbcdMathFunction(a, b, c, d));
    }

    PolynomialIntegralTenorBasis::PolynomialIntegralTenorBasis(
                                Date settlementDate,
                                shared_ptr<IborIndex> iborIndex,
                                const Handle<YieldTermStructure>& baseCurve,
                                shared_ptr<PolynomialFunction> p)
    : IntegralTenorBasis(settlementDate, iborIndex, baseCurve, p),
      instBasis_(p) {

        std::vector<Real> coef = 
                        instBasis_->definiteIntegralCoefficients(0.0, dt_);///dt_;
        for (Size i = 0; i < coef.size(); ++i){
            coef[i] = coef[i]/dt_;
        }
        basis_ = shared_ptr<PolynomialFunction>(new PolynomialFunction(coef));
    }

}

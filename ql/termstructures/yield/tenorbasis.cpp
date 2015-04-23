/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Ferdinando Ametrano

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
                                   shared_ptr<PureAbcdFunction> abcd)
    : TenorBasis(settlementDate, iborIndex, baseCurve), abcd_(abcd) {}

    PolynomialTenorBasis::PolynomialTenorBasis(
                                    Date settlementDate,
                                    shared_ptr<IborIndex> iborIndex,
                                    const Handle<YieldTermStructure>& baseCurve,
                                    shared_ptr<PolynomialFunction> p)
    : TenorBasis(settlementDate, iborIndex, baseCurve), p_(p) {}

    IntegralTenorBasis::IntegralTenorBasis(
                                Date settlementDate,
                                shared_ptr<IborIndex> iborIndex,
                                const Handle<YieldTermStructure>& baseCurve,
                                shared_ptr<unary_function<Real, Real> > b)
    : TenorBasis(settlementDate, iborIndex, baseCurve), basis_(b) {}

    //Real IntegralTenorBasis::forwardRate(Time t,
    //                                     Time t2) const {
    //    Real bigDelta = integrate(t, t2);
    //    DiscountFactor disc1 = baseCurve_->discount(t);
    //    DiscountFactor disc2 = baseCurve_->discount(t2);
    //    return (disc2/disc1*std::exp(bigDelta) - 1.0)/(t2-t);
    //}

    Real IntegralTenorBasis::value(Time t) const {
        Date d = dateFromTime(t);
        Date d2 = cal_.advance(d, tenor_, bdc_, eom_);
        Time t2 = timeFromSettlementDate(d2);
        Real bigDelta = integrate(t, t2);
        Time dt = t2-t;
        DiscountFactor disc1 = baseCurve_->discount(t);
        DiscountFactor disc2 = baseCurve_->discount(t2);
        Real discRatio = disc2/disc1;
        Rate fwd = (discRatio*std::exp(bigDelta) - 1.0)/dt;
        Rate fwdBase = (discRatio - 1.0)/dt;
        return fwd - fwdBase;
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
                                shared_ptr<PureAbcdFunction> abcd)
    : IntegralTenorBasis(settlementDate, iborIndex, baseCurve, abcd),
      abcd_(abcd) {

        Real a = abcd_->definiteIntegralA(0.0, dt_)/dt_;
        Real b = abcd_->definiteIntegralB(0.0, dt_)/dt_;
        Real c = abcd_->definiteIntegralC(0.0, dt_);
        Real d = abcd_->definiteIntegralD(0.0, dt_)/dt_;
        integratedBasis_ = shared_ptr<PureAbcdFunction>(new 
                                                PureAbcdFunction(a, b, c, d));
    }

    PolynomialIntegralTenorBasis::PolynomialIntegralTenorBasis(
                                Date settlementDate,
                                shared_ptr<IborIndex> iborIndex,
                                const Handle<YieldTermStructure>& baseCurve,
                                shared_ptr<PolynomialFunction> p)
    : IntegralTenorBasis(settlementDate, iborIndex, baseCurve, p),
      p_(p) {

    }

}

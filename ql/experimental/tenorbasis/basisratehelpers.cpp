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

#include <ql/experimental/tenorbasis/basisratehelpers.hpp>
#include <ql/experimental/tenorbasis/tenorbasis.hpp>
#include <ql/indexes/iborindex.hpp>

using boost::shared_ptr;

namespace QuantLib {

    namespace {
        void no_deletion(TenorBasis*) {}
    }

    BasisRateHelper::BasisRateHelper(const Handle<Quote>& basis,
                                     const Date& d)
    : BasisHelper(basis){
        earliestDate_ = d;
    }

    BasisRateHelper::BasisRateHelper(Rate basis,
                                     const Date& d)
    : BasisHelper(basis) {
        earliestDate_ = d;
    }

    void BasisRateHelper::setTermStructure(TenorBasis* t) {

        iborIndex_ = t->iborIndex();
        dc_ = iborIndex_->dayCounter();
        bdc_ = iborIndex_->businessDayConvention();
        eom_ = iborIndex_->endOfMonth();
        cal_ = iborIndex_->fixingCalendar();
        tenor_ = iborIndex_->tenor();
        latestDate_ = cal_.advance(earliestDate_, tenor_, bdc_, eom_);
        tau_ = dc_.yearFraction(earliestDate_, latestDate_);

        // do not set the relinkable handle as an observer -
        // force recalculation when needed
        bool observer = false;

        shared_ptr<TenorBasis> temp(t, no_deletion);
        termStructureHandle_.linkTo(temp, observer);

        baseCurveHandle_ = Handle<YieldTermStructure>(t->baseCurve());
        baseCurveRelinkableHandle_.linkTo(*baseCurveHandle_, observer);

        BasisHelper::setTermStructure(t);
    }

    Real BasisRateHelper::impliedQuote() const {
        QL_REQUIRE(termStructure_ != 0, "tenor basis not set");
        Rate result = termStructure_->tenorForwardRate(earliestDate_);
        return result;
    }

    void BasisRateHelper::accept(AcyclicVisitor& v) {
        Visitor<BasisRateHelper>* v1 =
            dynamic_cast<Visitor<BasisRateHelper>*>(&v);
        if (v1 != 0)
            v1->visit(*this);
        else
            BasisHelper::accept(v);
    }

}

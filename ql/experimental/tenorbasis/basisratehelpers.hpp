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

/*! \file ratehelpers.hpp
    \brief deposit, FRA, futures, and swap rate helpers
*/

#ifndef quantlib_ratehelpers_hpp
#define quantlib_ratehelpers_hpp

#include <ql/termstructures/bootstraphelper.hpp>

namespace QuantLib {

    class Quote;
    class TenorBasis;
    class IborIndex;

    typedef BootstrapHelper<TenorBasis> BasisHelper;

    class BasisRateHelper : public BasisHelper {
      public:
          BasisRateHelper(const Handle<Quote>& basis,
                          const Date& d);
          BasisRateHelper(Rate basis,
                          const Date& d);
        //! \name BasisHelper interface
        //@{
        Real impliedQuote() const;
        void setTermStructure(TenorBasis*);
        //@}
        //! \name Visitability
        //@{
        void accept(AcyclicVisitor&);
        //@}
      private:
        DayCounter dc_;
        BusinessDayConvention bdc_;
        bool eom_;
        Calendar cal_;
        Period tenor_;
        Time tau_;
        boost::shared_ptr<IborIndex> iborIndex_;
        RelinkableHandle<TenorBasis> termStructureHandle_;

        Handle<YieldTermStructure> baseCurveHandle_;
        RelinkableHandle<YieldTermStructure> baseCurveRelinkableHandle_;
    };

}

#endif

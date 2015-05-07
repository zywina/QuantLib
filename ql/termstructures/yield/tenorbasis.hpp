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

#ifndef quantlib_tenor_basis_hpp
#define quantlib_tenor_basis_hpp

#include <ql/math/abcdmathfunction.hpp>
#include <ql/math/polynomialmathfunction.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/calendar.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/handle.hpp>

namespace QuantLib {

    class IborIndex;

    //! Tenor (simple) basis between a given forwarding curve and a base curve
    /*!  */
    class TenorBasis {
      public:
        TenorBasis(Date settlementDate,
                   boost::shared_ptr<IborIndex> iborIndex,
                   const Handle<YieldTermStructure>& baseCurve);

        //! \name Interface
        //@{
        //! tenor basis as function of Date
        Spread value(Date d) const;
        //! tenor basis avoiding approximation as much as possible
        Spread exactValue(Date d) const;
        //! tenor basis as function of Time
        virtual Spread value(Time t) const = 0;

        //! return the value of the forward rate at date d
        Rate forwardRate(Date d) const;
        //! return the value of the forward rate at time t
        Rate forwardRate(Time t) const;

        virtual const std::vector<Real>& coefficients() = 0;
        virtual const std::vector<Real>& instCoefficients() = 0;
        //@}

        //! \name Improper (non-tenor) basis functions
        //@{
        //! tenor basis as integral of instantaneous continuous basis
        Spread value(Date d1,
                     Date d2) const;
        //! tenor basis as integral of instantaneous continuous basis
        Spread value(Time t1,
                     Time t2) const;
        //@}

        //! \name Improper (non-tenor) forward rate functions
        //@{
        // synthetic forward rate obtained as integral of inst. cont. basis
        Rate syntheticRate(Date d1,
                           Date d2) const;
        // synthetic forward rate obtained as integral of inst. cont. basis
        Rate syntheticRate(Time t1,
                           Time t2) const;
        //@}


        //! \name Dates and Time
        //@{
        //! the day counter used for date/time conversion
        DayCounter dayCounter() const { return dc_; }
        //! date-to-time conversion
        Time timeFromSettlementDate(Date d) const;
        //! time-to-date conversion
        Date dateFromTime(Time t) const;
        //@}

        //! Inspectors
        //@{
        //! settlement date for which t=0
        Date settlementDate() const { return settlementDate_; }
        //! IborIndex proving the forwarding curve
        const boost::shared_ptr<IborIndex>& iborIndex() const;
        //! Base curve used as reference for the basis
        const Handle<YieldTermStructure>& baseCurve() const;
        //@}
      protected:
        //! \name Integral functions
        //@{

        /*! \f[ I(d) = \int_{d}^{d+\tau} b(s)ds \f]
            with \f[ b(t) \f] being the instantaneous continuous basis
             and \f[ \tau \f] being the iborIndex tenor */
        Real integrate_(Date d) const;

        /*! \f[ I(d1, d2) = \int_{d1}^{d2} b(s)ds \f]
            with \f[ b(t) \f] being the instantaneous continuous basis */
        Real integrate_(Date d1,
                        Date d2) const;

        /*! \f[ I(t1, t2) = \int_{t1}^{t2} b(s)ds \f]
            with \f[ b(t) \f] being the instantaneous continuous basis
            TODO: possibly implement numerical integration as default */
        virtual Real integrate_(Time t1,
                                Time t2) const = 0;
        //@}
        Date settlementDate_;
        boost::shared_ptr<IborIndex> index_;
        Handle<YieldTermStructure> baseCurve_;

        DayCounter dc_;
        BusinessDayConvention bdc_;
        bool eom_;
        Calendar cal_;
        Period tenor_;
        Time dt_;
        Real time2date_;
    };

    class AbcdTenorBasis : public TenorBasis {
      public:
        AbcdTenorBasis(Date settlementDate,
                       boost::shared_ptr<IborIndex> iborIndex,
                       const Handle<YieldTermStructure>& baseCurve,
                       bool isSimple,
                       boost::shared_ptr<AbcdMathFunction> f);
        //! \name TenorBasis Interface
        //@{
        Spread value(Time t) const { return (*basis_)(t); }
        const std::vector<Real>& coefficients();
        const std::vector<Real>& instCoefficients();
        //@}

        //! date at which the simple tenor basis reaches maximum, if any
        Date maximumLocation() const;
        //! maximum values for the simple tenor basis, if any
        Spread maximumValue() const;
        Spread longTermValue() const { return basis_->d(); }

      protected:
        //! \name TenorBasis Interface
        //@{
        Real integrate_(Time t1,
                        Time t2) const;
        //@}
        boost::shared_ptr<AbcdMathFunction> basis_, instBasis_;
    };

    class PolynomialTenorBasis : public TenorBasis {
      public:
        PolynomialTenorBasis(Date settlementDate,
                             boost::shared_ptr<IborIndex> iborIndex,
                             const Handle<YieldTermStructure>& baseCurve,
                             bool isSimple,
                             boost::shared_ptr<PolynomialFunction> f);
        //! \name TenorBasis Interface
        //@{
        Spread value(Time t) const { return (*basis_)(t); }
        const std::vector<Real>& coefficients();
        const std::vector<Real>& instCoefficients();
        //@}
      protected:
        //! \name TenorBasis Interface
        //@{
        Real integrate_(Time t1,
                        Time t2) const;
        //@}
        boost::shared_ptr<PolynomialFunction> basis_, instBasis_;
    };

    // inline

}

#endif

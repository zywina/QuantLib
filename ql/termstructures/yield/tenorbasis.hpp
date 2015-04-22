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

#ifndef quantlib_tenor_basis_hpp
#define quantlib_tenor_basis_hpp

#include <ql/math/pureabcd.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/calendar.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/handle.hpp>

namespace QuantLib {

    class IborIndex;

    //! Tenor basis as simple basis between a given forwarding curve and a base curve
    class TenorBasis : public std::unary_function<Date, Real> {

    public:
        TenorBasis(Date settlementDate,
                   boost::shared_ptr<IborIndex> iborIndex,
                   const Handle<YieldTermStructure>& baseCurve);

        //! \name Interface
        //@{
        Real operator() (Date d) const;
        virtual Real value(Time t,
                           Time t2) const = 0;
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

        /*! Inspectors */
        Date settlementDate() { return settlementDate_; }
        boost::shared_ptr<IborIndex> iborIndex() { return index_; }
        //! the calendar used for date calculation
        Calendar calendar() const { return cal_; }
      protected:
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

    //! Tenor basis as definite integral of an instantaneous continuous basis
    class IntegralTenorBasis : public TenorBasis {

    public:
        IntegralTenorBasis(
                        Date settlementDate,
                        boost::shared_ptr<IborIndex> iborIndex,
                        const Handle<YieldTermStructure>& baseCurve,
                        boost::shared_ptr<std::unary_function<Real, Real> > b);

        //! \name TenorBasis Interface
        //@{
        Real value(Time t, Time t2) const;
        //@}

        //! \name integral functions
        //@{
        /*! \f[ G(d) = \int_{d}^{d+\tau} f(s)ds \f]
            with \f[ f(t) \f] being the instantaneous continuous basis
             and \f[ \tau \f] being the iborIndex tenor */
        Real integrate(Date d) const;

        /*! \f[ G(d1, d2) = \int_{d1}^{d2} f(s)ds \f]
            with \f[ f(t) \f] being the instantaneous continuous basis */
        Real integrate(Date d1, Date d2) const;

        /*! \f[ G(t1, t2) = \int_{t1}^{t2} f(s)ds \f]
            with \f[ f(t) \f] being the instantaneous continuous basis
            TODO: implement numerical integration as default */
        virtual Real integrate(Time t1, Time t2) const = 0;
        //@}

        boost::shared_ptr<std::unary_function<Real, Real> > basis() {
                                                            return basis_; }
      protected:
        boost::shared_ptr<std::unary_function<Real, Real> > basis_;
    };

    //! Tenor basis as definite integral of an instantaneous PureAbcdFunction basis
    /*! \f[ G(d) = \int_{d}^{d+\tau} f(s)ds \f]
        with \f[ \tau \f] being the iborIndex tenor
         and \f[ f(t) = [ a + b*t ] e^{-c*t} + d \f] */
    class AbcdIntegralTenorBasis : public IntegralTenorBasis {

      public:
        AbcdIntegralTenorBasis(Date settlementDate,
                               boost::shared_ptr<IborIndex> iborIndex,
                               const Handle<YieldTermStructure>& baseCurve,
                               boost::shared_ptr<PureAbcdFunction> abcd);
        //! \name IntegralTenorBasis Interface
        //@{
        Real integrate(Time t1, Time t2) const {
                                    return abcd_->definiteIntegral(t1, t2); }
        //@}

        /*! Inspectors */
        //! instantaneous continuous tenor basis
        boost::shared_ptr<PureAbcdFunction> instantaneousBasis() {
                                                                return abcd_; }

        //! integrated instantaneous continuous tenor basis
        boost::shared_ptr<PureAbcdFunction> integratedBasis() {
                                            return integratedBasis_; }
        Real a() const { return integratedBasis_->a(); }
        Real b() const { return integratedBasis_->b(); }
        Real c() const { return integratedBasis_->c(); }
        Real d() const { return integratedBasis_->d(); }

        /*! approximated date at which the integrated continuous
            tenor basis (i.e. simple basis) reaches maximum (if any) */
        Date maximumLocation() const;

        /*! approximated maximum values for the integrated continuous
            tenor basis (i.e. simple basis) */
        Real maximumValue() const;

        /*! value of the integrated continuous tenor basis (i.e. simple basis)
            at time 0: \f[ f(0) \f] */
        Real shortTermValue() const;

        /*! value of the integrated continuous tenor basis (i.e. simple basis)
            at time +inf: \f[ f(\inf) \f] */
        Real longTermValue() const;

      private:
        boost::shared_ptr<PureAbcdFunction> abcd_;
        boost::shared_ptr<PureAbcdFunction> integratedBasis_;
    };


    // inline

    inline Real TenorBasis::operator() (Date d) const {
        Date d2 = cal_.advance(d, tenor_, bdc_, eom_);
        Time t = timeFromSettlementDate(d);
        Time t2 = timeFromSettlementDate(d2);
        return value(t, t2);
    }

    inline Time TenorBasis::timeFromSettlementDate(Date d) const {
        return dc_.yearFraction(settlementDate_, d);
    }

    inline Date TenorBasis::dateFromTime(Time t) const {
        BigInteger result = settlementDate_.serialNumber() + BigInteger(t*time2date_);
        if (result >= Date::maxDate().serialNumber())
            return Date::maxDate();
        return Date(result);
    }

    inline Real IntegralTenorBasis::integrate(Date d1, Date d2) const {
        Time t1 = timeFromSettlementDate(d1);
        Time t2 = timeFromSettlementDate(d2);
        return integrate(t1, t2);
    }

    inline Date AbcdIntegralTenorBasis::maximumLocation() const {
        Time maximumLocation = integratedBasis_->maximumLocation();
        return dateFromTime(maximumLocation);
    }

    inline Real AbcdIntegralTenorBasis::maximumValue() const {
        Date d = maximumLocation();
        return TenorBasis::operator()(d);
    }

    inline Real AbcdIntegralTenorBasis::shortTermValue() const {
        return integratedBasis_->shortTermValue();
    }

    inline Real AbcdIntegralTenorBasis::longTermValue() const {
        return integratedBasis_->longTermValue();
    }

}

#endif

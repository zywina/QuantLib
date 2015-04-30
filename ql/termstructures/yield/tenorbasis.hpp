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
#include <ql/time/daycounter.hpp>
#include <ql/time/calendar.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/handle.hpp>

namespace QuantLib {

    class IborIndex;

    //! Tenor (simple) basis between a given forwarding curve and a base curve
    /*! \f[ B(d) = iborIndex_{\tau}(d)-F_{\tau}(d) \f] where
        \f[ F_{\tau}(d) \f] is the forward rate calculated on the base curve
        and \f[ \tau \f] is the iborIndex tenor */
    class TenorBasis : public std::unary_function<Date, Real> {
      public:
        TenorBasis(Date settlementDate,
                   boost::shared_ptr<IborIndex> iborIndex,
                   const Handle<YieldTermStructure>& baseCurve);

        //! \name Interface
        //@{
        //! simple tenor basis as function of Date
        Real operator() (Date d) const;

        //! simple tenor basis as function of Time
        virtual Real value(Time t) const = 0;
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
        Date settlementDate() { return settlementDate_; }
        //! IborIndex proving the forwarding curve
        const boost::shared_ptr<IborIndex>& iborIndex() { return index_; }
        //! Base curve used as reference for the basis
        const Handle<YieldTermStructure>& baseCurve() { return baseCurve_; }
        //! the calendar used for date calculation
        Calendar calendar() const { return cal_; }
        //@}
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

    //! Tenor basis as abcd parametrization of simple basis
    /*! \f[ B(t) = [ a + b*t ] e^{-c*t} + d \f] */
    class AbcdTenorBasis : public TenorBasis {
      public:
        AbcdTenorBasis(Date settlementDate,
                       boost::shared_ptr<IborIndex> iborIndex,
                       const Handle<YieldTermStructure>& baseCurve,
                       boost::shared_ptr<AbcdMathFunction> basis);
        //! \name TenorBasis Interface
        //@{
        Real value(Time t) const { return (*basis_)(t); }
        //@}
        //! Inspectors
        //@{
        const boost::shared_ptr<AbcdMathFunction>& basis() { return basis_; }
        //@}
      protected:
        boost::shared_ptr<AbcdMathFunction> basis_;
    };

    //! Tenor basis as polynomial parametrization of simple basis
    /*! \f[ B(t) = \sum_0^n{c_i*t^i} \f] */
    class PolynomialTenorBasis : public TenorBasis {
      public:
        PolynomialTenorBasis(Date settlementDate,
                             boost::shared_ptr<IborIndex> iborIndex,
                             const Handle<YieldTermStructure>& baseCurve,
                             boost::shared_ptr<PolynomialFunction> basis);
        //! \name TenorBasis Interface
        //@{
        Real value(Time t) const { return (*basis_)(t); }
        //@}
        //! Inspectors
        //@{
        const boost::shared_ptr<PolynomialFunction>& basis() { return basis_; }
        //@}
     protected:
        boost::shared_ptr<PolynomialFunction> basis_;
    };


    //! Tenor basis as definite integral of an instantaneous continuous basis
    /*! \f (((1+F_{\tau}(d)\tau) \exp{I(d)} - 1)/\tau - F_{\tau}(d)) \f] where
        \f[ I(d) = \int_{d}^{d+\tau} b(s)ds \f] and 
        \f[ \tau \f] being the iborIndex tenor */
    class IntegralTenorBasis : public TenorBasis {
      public:
        IntegralTenorBasis(
                        Date settlementDate,
                        boost::shared_ptr<IborIndex> iborIndex,
                        const Handle<YieldTermStructure>& baseCurve,
                        boost::shared_ptr<std::unary_function<Real, Real> >);
        //! \name TenorBasis Interface
        //@{
        Real value(Time t) const;
        //@}

        //! \name Simple Basis functions
        //@{
        //! simple tenor basis as integral function between two times
        Real value(Time t1,
                   Time t2) const;

        //! simple tenor basis as integral function between two dates
        Real value(Date d1,
                   Date d2) const;
        //@}

        //! \name Forward Rate functions
        //@{
        // return the value of the forward rate between d1 and d2
        Real fwdRate(Date d1,
                     Date d2) const;

        //! simple tenor basis as integral function between d and d+tau
        Real fwdRate(Date d) const;
        //@}

        //! \name Integral functions
        //@{

        /*! \f[ I(d) = \int_{d}^{d+\tau} b(s)ds \f]
            with \f[ b(t) \f] being the instantaneous continuous basis
             and \f[ \tau \f] being the iborIndex tenor */
        Real integrate(Date d) const;

        /*! \f[ I(d1, d2) = \int_{d1}^{d2} b(s)ds \f]
            with \f[ b(t) \f] being the instantaneous continuous basis */
        Real integrate(Date d1,
                       Date d2) const;

        /*! \f[ I(t1, t2) = \int_{t1}^{t2} b(s)ds \f]
            with \f[ b(t) \f] being the instantaneous continuous basis
            TODO: possibly implement numerical integration as default */
        virtual Real integrate(Time t1,
                               Time t2) const = 0;
        //@}

        const boost::shared_ptr<std::unary_function<Time, Real> >& instBasis();
      protected:
        boost::shared_ptr<std::unary_function<Time, Real> > instBasis_;
    };

    //! Tenor basis as definite integral of an instantaneous abcd continuous basis
    /*! \f (((1+F_{\tau}(t)\tau) \exp{I(t)} - 1)/\tau - F_{\tau}(t)) \f] where
        \f[ I(t) = \int_{t}^{t+\tau} b(s)ds \f], 
        \f[ \tau \f] being the iborIndex tenor,
         and \f[ b(t) = [ a + b*t ] e^{-c*t} + d \f] */
    class AbcdIntegralTenorBasis : public IntegralTenorBasis {
      public:
        AbcdIntegralTenorBasis(Date settlementDate,
                               boost::shared_ptr<IborIndex> iborIndex,
                               const Handle<YieldTermStructure>& baseCurve,
                               boost::shared_ptr<AbcdMathFunction> instBasis);
        //! \name IntegralTenorBasis Interface
        //@{
        Real integrate(Time t1,
                       Time t2) const;
        //@}
        /*! Inspectors */
        //@{
        //! instantaneous continuous tenor basis
        const boost::shared_ptr<AbcdMathFunction>& instBasis();

        //! simple basis, i.e. integrated instantaneous continuous basis
        const boost::shared_ptr<AbcdMathFunction>& basis();

        /*! parameters of the simple basis, i.e. of the integrated
            instantaneous continuous basis.  Not to be confused with
            the {a, b, c, d} parameters of the instantaneous basis */
        Real a() const { return basis_->a(); }
        Real b() const { return basis_->b(); }
        Real c() const { return basis_->c(); }
        Real d() const { return basis_->d(); }
        //@}

        /*! approximated date at which the simple tenor basis (i.e.
            the integrated continuous basis) reaches maximum, if any */
        Date maximumLocation() const;

        /*! approximated maximum values for the simple tenor basis (i.e.
            the integrated continuous basis), if any */
        Real maximumValue() const;

      private:
        boost::shared_ptr<AbcdMathFunction> instBasis_, basis_;
    };

    //! Tenor basis as definite integral of an instantaneous polynomial continuous basis
    /*! \f (((1+F_{\tau}(t)\tau) \exp{I(t)} - 1)/\tau - F_{\tau}(t)) \f] where
        \f[ I(t) = \int_{t}^{t+\tau} b(s)ds \f], 
        \f[ \tau \f] being the iborIndex tenor,
         and \f[ b(t) = \sum_0^n{c_i*t^i} \f] */
    class PolynomialIntegralTenorBasis : public IntegralTenorBasis {
      public:
        PolynomialIntegralTenorBasis(Date settlementDate,
                                     boost::shared_ptr<IborIndex> iborIndex,
                                     const Handle<YieldTermStructure>& b,
                                     boost::shared_ptr<PolynomialFunction>);
        //! \name IntegralTenorBasis Interface
        //@{
        Real integrate(Time t1,
                       Time t2) const;
        //@}

        /*! Inspectors */
        //@{
        //! instantaneous continuous tenor basis
        const boost::shared_ptr<PolynomialFunction>& instBasis();

        //! simple basis, i.e. integrated instantaneous continuous basis
        // const boost::shared_ptr<PolynomialFunction>& basis();
        //@}
      private:
        // boost::shared_ptr<PolynomialFunction> instBasis_, basis_;
        boost::shared_ptr<PolynomialFunction> instBasis_;
    };


    // inline

    inline Real TenorBasis::operator() (Date d) const {
        Time t = timeFromSettlementDate(d);
        return value(t);
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



    inline Real IntegralTenorBasis::integrate(Date d1,
                                              Date d2) const {
        Time t1 = timeFromSettlementDate(d1);
        Time t2 = timeFromSettlementDate(d2);
        return integrate(t1, t2);
    }

    // return the value of the simple basis between d1 and d2
    inline Real IntegralTenorBasis::value(Date d1,
                                          Date d2) const {
        Time t1 = timeFromSettlementDate(d1);
        Time t2 = timeFromSettlementDate(d2);
    //    Rate fwdBase = baseCurve_->forwardRate(d1, d2, dc_, Simple, Annual, 0);
        return value(t1, t2);
    }

    // return the value of the forward rate between d1 and d2
    inline Real IntegralTenorBasis::fwdRate(Date d1,
                                            Date d2) const {
        Time t1 = timeFromSettlementDate(d1);
        Time t2 = timeFromSettlementDate(d2);
        Rate baseCurveFwd =
            baseCurve_->forwardRate(d1, d2, dc_, Simple, Annual, 0);
        return value(t1, t2) + baseCurveFwd;
    }

    // return the value of the forward rate basis between d and d+tau
    inline Real IntegralTenorBasis::fwdRate(Date d) const {
        Date d2 = cal_.advance(d, tenor_, bdc_, eom_);
        return fwdRate(d, d2);
}

    inline const boost::shared_ptr<std::unary_function<Time, Real> >&
    IntegralTenorBasis::instBasis() {
        return instBasis_;
    }



    inline Real AbcdIntegralTenorBasis::integrate(Time t1,
                                                  Time t2) const {
        return instBasis_->definiteIntegral(t1, t2);
    }

    inline const boost::shared_ptr<AbcdMathFunction>&
    AbcdIntegralTenorBasis::instBasis() {
        return instBasis_;
    }

    inline const boost::shared_ptr<AbcdMathFunction>&
    AbcdIntegralTenorBasis::basis() {
        return basis_;
    }

    inline Date AbcdIntegralTenorBasis::maximumLocation() const {
        Time maximumLocation = basis_->maximumLocation();
        return dateFromTime(maximumLocation);
    }

    inline Real AbcdIntegralTenorBasis::maximumValue() const {
        Date d = maximumLocation();
        return TenorBasis::operator()(d);
    }


    inline Real PolynomialIntegralTenorBasis::integrate(Time t1,
                                                        Time t2) const {
        return instBasis_->definiteIntegral(t1, t2);
    }

    inline const boost::shared_ptr<PolynomialFunction>&
    PolynomialIntegralTenorBasis::instBasis() {
        return instBasis_;
    }

    //inline const boost::shared_ptr<PolynomialFunction>&
    //PolynomialIntegralTenorBasis::basis() {
    //    return basis_;
    //}

}

#endif

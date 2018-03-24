/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Roy Zywina

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

/*! \file powell.hpp
    \brief Powell's optimization method
*/

#ifndef quantlib_optimization_powell_hpp
#define quantlib_optimization_powell_hpp

#include <stdio.h>
#include <ql/math/comparison.hpp>
#include <ql/utilities/null.hpp>
#include <ql/patterns/curiouslyrecurring.hpp>
#include <ql/errors.hpp>
#include <ql/math/optimization/linesearchbasedmethod.hpp>
#include <ql/math/optimization/linesearch.hpp>
#include <ql/math/optimization/problem.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace QuantLib{

  template <class Impl>
  class Minimizer1D : public CuriouslyRecurringTemplate<Impl> {
  public:
    Minimizer1D(){}

    template <class F>
    Real minimize(const F& f,
               Real accuracy,
               Real guess,
               Real step)
    {
      xMin_=guess;
      xMax_=guess+step;
      bracket(f,xMin_,xMax_,xMid_,fxMin_,fxMax_,fxMid_);
      return this->impl().minimizeImpl(f, accuracy);
    }

    template <class F>
    void bracket(const F& f,
      Real &ax,Real &bx,Real &cx,
      Real &fa,Real &fb,Real &fc)
    {
      fa = f(ax);
      fb = f(bx);
      if (fb > fa){
        std::swap(fa,fb);
        std::swap(ax,bx);
      }
      const Real growth = 1.6;
      cx = bx + growth*(bx-ax);
      fc = f(cx);
      while (fb > fc){
        // TODO: implement parabolic optimization
        Real u = cx + growth*(cx-bx);
        Real fu = f(u);
        // remove oldest point
        ax=bx; bx=cx; cx=u;
        fa=fb; fb=fc; fc=fu;
      }
      if (ax > bx){
        std::swap(fa,fb);
        std::swap(ax,bx);
      }
    }
  protected:
    Real xMin_,xMax_,xMid_,fxMin_,fxMax_,fxMid_;
  };

  class GoldenSection : public Minimizer1D<GoldenSection>{
  public:
    template <class F>
    Real minimizeImpl(const F& f, Real xAccuracy) {
      const Real gold = (sqrt(5.0)+1)/2;
      Real xc,xd,yc,yd;
      do{
        xc = xMax_ - (xMax_ - xMin_) / gold;
        xd = xMin_ + (xMax_ - xMin_) / gold;
        yc = f(xc);
        yd = f(xd);
        if (yc < yd){
          xMax_ = xd;
          fxMax_ = yd;
          xMid_ = xc;
          fxMid_ = yc;
        }else{
          xMin_ = xc;
          fxMin_ = yc;
          xMid_ = xd;
          fxMid_ = yd;
        }
        //printf("%.8f\n",std::abs(xc-xd));
      }while (std::abs(xc-xd)>xAccuracy);
      return (xMin_+xMax_)/2;
    }
  };

  class BrentMinimizer : public Minimizer1D<BrentMinimizer>{
  public:
    template <class F>
    Real minimizeImpl(const F& f, Real xAccuracy) {
      Real a,b,d,x,w,v,fx,fv,fw,xm,p,q,r,u,fu;
      Real tol1,tol2,tol=0.001,e=0,etemp;
      a = xMin_;
      b = xMax_;
      x=w=v= xMid_;
      fw=fv=fx = fxMid_;
      const Real CGOLD = 0.3819660;
      for (int iter=0; iter<100; iter++){
        xm = (a+b)/2;
        tol1 = tol*std::abs(x)+QL_EPSILON;
        tol2 = 2*tol1;
        if (std::abs(x-xm) <= (tol2-0.5*(b-a))){
          return x;
        }
        if (std::abs(e) > tol1){
          r = (x-w)*(fx-fv);
          q = (x-v)*(fx-fw);
          p = (x-v)*q-(x-w)*r;
          q=2.0*(q-r);
          if (q>0.0) p = -p;
          q = std::abs(q);
          etemp=e;
          e=d;
          if (std::abs(p) >= std::abs(0.5*q*etemp) ||
            p <= q*(a-x) || p >= q*(b-x))
          {
            e = x >= xm ? a-x : b-x;
            d = CGOLD * e;
          }else{
            d = p/q;
            u = x+d;
            if (u-a < tol2 || b-u < tol2)
              d = sign(tol1,xm-x);
          }
        }else{
          e = x >= xm ? a-x : b-x;
          d = CGOLD * e;
        }
        u = std::abs(d) >= tol1 ? x+d : x+sign(tol1,d);
        fu = f(u);
        //printf("Brent f(%f) = %f\n",u,fu);
        if (fu <= fx){
          if (u>=x)
            a=x;
          else
            b=x;
          v=w; w=x; x=u;
          fv=fw; fw=fx; fx=fu;
        }else{
          if (u<x) a=u; else b=u;
          if (fu <= fw || w==x){
            v=w;
            w=u;
            fv=fw;
            fw=fu;
          }else if (fu<=fv || v==x || v==w){
            v=u;
            fv=fu;
          }
        }
      }
      return x;
    }
    Real sign(Real a, Real b) const {
        return b >= 0.0 ? std::fabs(a) : -std::fabs(a);
    }
  };

  //template <class Minimizer = BrentMinimizer>
  class Powell : public OptimizationMethod {
    public:
      Powell(
        Real step=0,
        bool positiveOptimization = true,
        std::vector<Real> steps = std::vector<Real>()
      ) :
        OptimizationMethod(),
        step_(step),
        steps_(steps),
        positiveOptimization_(positiveOptimization)
        {}
      Powell(
        const std::vector<Array> &directions,
        const std::vector<Real> &steps,
        bool positiveOptimization = true
      ) :
        OptimizationMethod(),
        direction_(directions),
        step_(0.0),
        steps_(steps),
        positiveOptimization_(positiveOptimization)
      {
        QL_REQUIRE(steps.size()==direction_.size()+1,
          "Powell requires one more step size than direciton count");
      }

      virtual EndCriteria::Type minimize(Problem& P,
            const EndCriteria& endCriteria)
      {
        EndCriteria::Type ecType = EndCriteria::None;
        P.reset();
        current_ = P.currentValue();
        Size N = current_.size();
        if (direction_.empty()){
          rebuildDirections(N);
        }

        if (step_==0 && steps_.empty()){
          QL_FAIL("No step size given in Powell");
        }else if (step_!=0){
          steps_.resize(N+1, step_);
        }else{
          QL_REQUIRE(steps_.size()==direction_.size()+1,
            "Powel requires N+1 step sizes");
        }
        Real xtol = endCriteria.rootEpsilon();

        boost::function<Real(Real)> f =
          boost::bind(&Powell::calculate,this,boost::ref(P),_1);

        Size statStateIterations = 0;
        dir_ = Array(N,0.0);
        Real fPrev = f(0.0);
        endCriteria.checkStationaryFunctionAccuracy(
          P.functionValue(),positiveOptimization_,ecType);

        for (Size iteration=0; ecType==EndCriteria::None; iteration++){
          Array average(N,0);
          for (size_t j=0; j<direction_.size(); j++){
            current_ = P.currentValue();
            dir_ = direction_[j];
            Real x = BrentMinimizer().minimize(
              f,xtol, 0.0, steps_[j]);
            dir_ *= x;
            average += dir_;
          }
          current_ = P.currentValue();
          dir_ = average;
          // search the average movement vector
          BrentMinimizer().minimize(
            f,xtol, 0.0, steps_[steps_.size()-1]);

          endCriteria.checkMaxIterations(iteration,ecType);
          endCriteria.checkStationaryFunctionAccuracy(
            P.functionValue(),positiveOptimization_,ecType);
          endCriteria.checkStationaryFunctionValue(
            fPrev, P.functionValue(), statStateIterations, ecType);
          fPrev = P.functionValue();
        }
        return ecType;
      }
    private:
      // direction vectors
      std::vector<Array> direction_;
      Array current_, dir_;
      std::vector<Real> steps_;
      Real step_;
      bool positiveOptimization_;

      // default to unit direction set
      void rebuildDirections(Size N){
        direction_.resize(N);
        for (Size i=0; i<direction_.size(); i++){
          direction_[i] = Array(N,0.0);
          direction_[i][i] = 1.0;
        }
      }

      Real calculate(Problem &P, Real x){
        Array xx = dir_;
        xx *= x;
        xx += current_;
        Real err = P.value(xx);
        if (err < P.functionValue()){
          P.setFunctionValue(err);
          P.setCurrentValue(xx);
          std::cout << xx << " => " << err << std::endl;
        }
        return err;
      }
  };

  /*
  Powell::Powell() :
    OptimizationMethod()
  {
  }

  EndCriteria::Type Powell::minimize(Problem& P, const EndCriteria& endCriteria)
  {
    EndCriteria::Type ecType = EndCriteria::None;
    P.reset();
    current_ = P.currentValue();
    Size iterationNumber_ = 0;
    Size N = current_.size();
    if (direction_.size()!=N)
      rebuildDirections(N);

    boost::function<Real(Real)> f =
      boost::bind(&Powell::calculate,this,&P,_1);
    for (int i=0; i<2; i++){
      Array avg(N,0);
      for (size_t j=0; j<N; j++){
        current_ = P.currentValue();
        dir_ = direction_[j];
        Real x = BrentMinimizer().minimize(
          boost::bind(&Powell::calculate,this,&P,_1),
          1e-6, 0.0, 0.01);
        dir_ *= x;
        avg += dir_;
      }
      current_ = P.currentValue();
      dir_ = avg;
      // search the average movement vector
      BrentMinimizer().minimize(
        boost::bind(&Powell::calculate,this,&P,_1),
        1e-6, 0.0, 0.01);
    }

    return ecType;
  }

  Real Powell::calculate(Problem *P, Real x){
    Array xx = dir_;
    xx *= x;
    xx += current_;
    //for (size_t i=0; i<current_.size(); i++){
    //  xx[i] = current_[i] + dir_[i]*x;
    //}
    Real err = P->value(xx);
    if (err < P->functionValue()){
      P->setFunctionValue(err);
      P->setCurrentValue(xx);
      for (Size i=0; i<xx.size(); i++)
        printf("%f ",xx[i]);
      printf(" =  %.8f\n",err);
    }
    return err;
  }

  void Powell::rebuildDirections(Size N){
    // default to unit direction vectors
    direction_.resize(N);
    for (Size i=0; i<direction_.size(); i++){
      direction_[i] = Array(N,0.0);
      direction_[i][i] = 1.0;

      for (Size j=0; j<N; j++){
        direction_[i][j] = cos(i*M_PI*j/(Real)(N-1));
      }
    }
  }
  */

};

#endif

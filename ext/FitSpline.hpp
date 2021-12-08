#ifndef FITSPLINE_HPP_
#define FITSPLINE_HPP_

#include <memory>
#include <valarray>

#include "skynet/Utilities/FunctionVsTime.hpp"

class
FitSpline : public FunctionVsTime<double> {
public:
  FitSpline(
    const std::valarray<double>& times, 
    const std::valarray<double>& values);
  
  double operator()(const double time) const;

  std::unique_ptr<FunctionVsTime<double>> MakeUniquePtr() const;

  void print_coefs() const;

private:
  std::valarray<double> mTimes;
  std::valarray<double> mCfsA;
  std::valarray<double> mCfsB;
  std::valarray<double> mCfsC;
  std::valarray<double> mCfsD;
  
  int interpolate(const std::valarray<double>& values);

};

#endif // FITSPLINE_HPP_

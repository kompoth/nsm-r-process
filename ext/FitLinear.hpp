#ifndef FITLINEAR_HPP_
#define FITLINEAR_HPP_

#include <memory>
#include <valarray>

#include "skynet/Utilities/FunctionVsTime.hpp"

class
FitLinear : public FunctionVsTime<double> {
public:
  FitLinear(
    const std::valarray<double>& times, 
    const std::valarray<double>& values);
   
  double operator()(const double time) const;

  std::unique_ptr<FunctionVsTime<double>> MakeUniquePtr() const;

  void print_coefs() const;

private:
  std::valarray<double> mTimes;
  std::valarray<double> mCfsA;
  std::valarray<double> mCfsB;
  
  FitLinear(
    const std::valarray<double>& times, 
    const std::valarray<double>& coefs_a,
    const std::valarray<double>& coefs_b);
  
  int interpolate(const std::valarray<double>& values);
};

#endif // FITLINEAR_HPP_

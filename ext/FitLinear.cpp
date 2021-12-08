#include "FitLinear.hpp"

#include <cmath>
#include <iterator>
#include <stdexcept>

/* Constructors */
FitLinear::FitLinear(
  const std::valarray<double>& times, 
  const std::valarray<double>& values):
  mTimes(times) {
  if (mTimes.size() != values.size()) 
    throw std::invalid_argument("Data arrays must be of the same size\n");
  if (mTimes.size() < 2) 
    throw std::invalid_argument("At least 2 data points are required\n");
  interpolate(values);
}

FitLinear::FitLinear(
  const std::valarray<double>& times, 
  const std::valarray<double>& coefs_a,
  const std::valarray<double>& coefs_b):
  mTimes(times),
  mCfsA(coefs_a),
  mCfsB(coefs_b) {
  if (mTimes.size() != mCfsA.size() && mTimes.size() != mCfsB.size()) 
    throw std::invalid_argument("Data arrays must be of the same size\n");
  if (mTimes.size() < 2) 
    throw std::invalid_argument("At least 2 data points are required\n");
}

/* Perform linear interpolation */
int
FitLinear::interpolate(const std::valarray<double>& values) {
  size_t npts = mTimes.size();
  double h;
  mCfsA = std::valarray<double>(npts);
  mCfsB = std::valarray<double>(npts);

  for (size_t k = 1; k < npts; k++) {
    h = mTimes[k] - mTimes[k-1];
    if (h < 0.)
      throw std::invalid_argument("Time values decreasing\n");
    mCfsA[k] = values[k-1];
    mCfsB[k] = (values[k] - values[k-1]) / h;  
  }
   
  return 0;
}

/* Base class functions */
double 
FitLinear::operator()(const double time) const {
  if (time > mTimes.max() || time < mTimes.min())
    throw std::invalid_argument("No data for given time " 
      + std::to_string(time) + "\n");
  
  auto left = std::begin(mTimes);
  auto right = std::end(mTimes);
  
  while (left + 1 != right) {
    auto middle = left + std::distance(left, right) / 2;
    if (*middle > time)
      right = middle;  
    else
      left = middle;
  }

  size_t k = right - std::begin(mTimes);
  double dt = time - *left;
  return mCfsA[k] + mCfsB[k] * dt;
}

std::unique_ptr<FunctionVsTime<double>> 
FitLinear::MakeUniquePtr() const {
  return std::unique_ptr<FunctionVsTime<double>>(
    new FitLinear(mTimes, mCfsA, mCfsB));
}

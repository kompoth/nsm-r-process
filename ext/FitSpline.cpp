#include "FitSpline.hpp"

#include <cmath>
#include <iterator>
#include <stdexcept>

/* Constructor */
FitSpline::FitSpline(
  const std::valarray<double>& times, 
  const std::valarray<double>& values):
  mTimes(times) {
  // Check size
  if (mTimes.size() != values.size()) 
    throw std::invalid_argument("Data arrays must be of the same size\n");
  interpolate(values);
}

/* TODO: Need another constructor for already interpolated data */
/* TODO: Sort data by time */

/* Perform spline interpolation */
int
FitSpline::interpolate(const std::valarray<double>& values) {
  size_t k;
  size_t npts = mTimes.size();
  std::valarray<double> h(npts), l(npts);
  std::valarray<double> H(npts), F(npts);
  double mul;
  
  mCfsA = values.shift(-1);
  mCfsB = std::valarray<double>(npts);
  mCfsC = std::valarray<double>(npts);
  mCfsD = std::valarray<double>(npts);
  
  for (k = 1; k < npts; k++) {
    h[k] =  mTimes[k] - mTimes[k-1];
    l[k] =  mCfsA[k] - mCfsA[k-1];
  }
  
  for (k = 1; k < npts - 1; k++) {
    H[k] = 2. * (h[k] + h[k+1]);
    F[k] = 3. * (h[k] * l[k+1] - h[k+1] * l[k]) / (h[k] * h[k+1]);
  }

  for (k = 2; k < npts - 1; k++) {
    mul = h[k] / H[k-1];
    H[k] -= mul * h[k];
    F[k] -= mul * F[k-1];
  }
  
  mCfsC[1] = 0.;
  mCfsC[npts-1] = F[npts-2] / H[npts-2];
  for (k = npts - 2; k > 1; k--) {
    mCfsC[k] = (F[k-1] -  mCfsC[k+1] * h[k]) / H[k-1]; 
  } 
  
  for (k = 1; k < npts; k++) {
    mCfsB[k] = l[k] / h[k] - h[k] * h[k] * (mCfsC[k+1] - mCfsC[k]) / 3.;
    mCfsD[k] = (mCfsC[k+1] - mCfsC[k]) / (3. * h[k]);
  }

  return 0;
}

/* Print interpolation coefficients */
void
FitSpline::print_coefs() const {
  size_t k;
  printf("# S = a + b * dt + c * dt^2 + d * dt^3\n");  
  printf("#%16s%16s%16s%16s%16s%16s\n", "t1", "t2", "a", "b", "c", "d");  
  for (k = 1; k < mTimes.size(); k++)
    printf(" %16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n", 
      mTimes[k-1], mTimes[k], mCfsA[k], mCfsB[k], mCfsC[k], mCfsD[k]);  
}

/* Base class functions */
double 
FitSpline::operator()(const double time) const {
  if (time > mTimes.max() || time < mTimes.min())
    throw std::invalid_argument("No data for given time\n");
  
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
  return mCfsA[k] + mCfsB[k] * dt + mCfsC[k] * pow(dt, 2) \
    + mCfsD[k] * pow(dt, 3);
}

std::unique_ptr<FunctionVsTime<double>> 
FitSpline::MakeUniquePtr() const {
  return std::unique_ptr<FunctionVsTime<double>>(
    new FitSpline(mTimes, mCfsA));
}

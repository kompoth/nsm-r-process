#include "LinearExpansion.hpp"

#include <cmath>
#include <iterator>
#include <stdexcept>

LinearExpansion::LinearExpansion(const double rho0, const double tau):
  mRho0(rho0),
  mTau(tau) {
  if (mRho0 <= 0.0)
    throw std::invalid_argument("Rho0 for an LinearExpansion density profile must "
      "be positive");
  if (mTau <= 0.0)
    throw std::invalid_argument("Tau for an LinearExpansion density profile must "
      "be positive");
}

double LinearExpansion::operator()(const double time) const {
  if (time < 1e-20)
    return mRho0;
  return mRho0 * pow(mTau / time, 3.);
}

std::unique_ptr<FunctionVsTime<double>> LinearExpansion::MakeUniquePtr() const {
  return std::unique_ptr<FunctionVsTime<double>>(new LinearExpansion(mRho0, mTau));
}

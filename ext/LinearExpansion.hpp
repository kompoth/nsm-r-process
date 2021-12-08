#ifndef LINEAREXPANSION_HPP_
#define LINEAREXPANSION_HPP_

#include <memory>

#include "skynet/Utilities/FunctionVsTime.hpp"

class LinearExpansion : public FunctionVsTime<double> {
public:
  LinearExpansion(const double rho0, const double tau);

  double operator()(const double time) const;

  std::unique_ptr<FunctionVsTime<double>> MakeUniquePtr() const;

private:
  double mRho0; // in g / cm^3
  double mTau; // in seconds
};

#endif // LINEAREXPANSION_HPP_

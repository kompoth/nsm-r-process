#include <fstream>
#include <vector>
#include <valarray>
#include <string>

#include "FitLinear.hpp"
#include "FitSpline.hpp"
#include "LinearExpansion.hpp"

int
main(int argc, char *argv[]) {
  std::vector<double> times;
  std::vector<double> values;
  double t, val;
  std::string path;

  if (argc > 1)
    path = argv[1];

  std::ifstream temp(path);
  if (!temp)
    throw std::runtime_error("error loading data file");
  while (temp >> t >> val) {
    times.push_back(t);
    values.push_back(val);
  }
  temp.close();
  
  std::valarray<double> arr_times(times.data(), times.size());
  std::valarray<double> arr_values(values.data(), values.size());
  FitSpline profile (arr_times, arr_values); 
  
  int npts = 1000;
  double dt = 1. / npts; 
  t = 0.;
  for (int i = 0; i < npts; i++) {
    printf("%16.6e%16.6e\n", t, profile(t));
    t += dt;
  }

  return 0;
}

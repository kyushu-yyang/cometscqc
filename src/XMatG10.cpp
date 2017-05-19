#include <cmath>
#include "XQuenchLogger.hpp"
#include "XQuenchExcept.hpp"
#include "XMatG10.hpp"

using Quench::XQuenchLogger;

XMatG10 :: XMatG10()
    : fDirection("Normal")
{}


XMatG10 :: ~XMatG10()
{}


void XMatG10 :: SetDirection(const std::string &dir)
{
  fDirection = dir;

  if (fDirection!="Wrap" || fDirection!="Normal") {
    XQuenchExcept except("this direction did not exist.");
    throw except;
  }
}


double XMatG10 :: GetConductivity()
{
  const double k = GetNistConductivity(fTemp);
  return k;
}


double XMatG10 :: GetNistConductivity(const double T)
{
  // parameters
  double p[9] = {-4.1236, 13.788, -26.068, 26.272, -14.663,
                 4.4954, -0.6905, 0.0397, 0.};

  if (fDirection=="Wrap") {
    p[0] = -2.64827;
    p[1] = 8.80228;
    p[2] = -24.8998;
    p[3] = 41.1625;
    p[4] = -39.8754;
    p[5] = 23.1778;
    p[6] = -7.95635;
    p[7] = 1.48806;
    p[8] = -0.11701;
  }

  double ax = 0.;

  for (int i=0; i<9; i++) {
    ax += p[i]*pow(log10(T),i);
  }

  const double k = pow(10., ax);

  return k;
}


double XMatG10 :: GetCapacity()
{
  const double C = GetNistCapacity(fTemp);
  return C;
}


double XMatG10 :: GetNistCapacity(const double T)
{
  const double p[9] = {-2.4083, 7.6006, -8.2982, 7.3301, -4.2386,
                       1.4294, -0.24396, 0.015236, 0.};

  double ax = 0.;

  for (int i=0; i<9; i++) {
    ax += p[i]*pow(log10(T), i);
  }

  const double C = pow(10., ax);

  return C;
}
 

#include <cmath>
#include "IFdmUnits.hpp"
#include "XQuenchLogger.hpp"
#include "XQuenchExcept.hpp"
#include "XMatAl5083.hpp"

using Quench::XQuenchLogger;

XMatAl5083 :: XMatAl5083()
{}


XMatAl5083 :: ~XMatAl5083()
{}


double XMatAl5083 :: GetConductivity()
{
  const double k = GetNistConductivity(fTemp);
  return k;
}


double XMatAl5083 :: GetNistConductivity(const double T)
{
  // parameters
  const double p[9] = {-0.90933, 5.751, -11.112, 13.612, -9.3977,
                       3.6873, -0.77295, 0.067336, 0.};

  double ax = 0.;

  for (int i=0; i<9; i++) {
    ax += p[i]*pow(log10(T),i);
  }

  const double k = pow(10., ax);

  return k;
}


double XMatAl5083 :: GetCapacity()
{
  const double C = GetNistCapacity(fTemp);
  return C;
}


double XMatAl5083 :: GetNistCapacity(const double T)
{
  // parameters
  const double p[9] = {46.6467, -314.292, 866.662, -1298.3, 1162.27,
                       -637.795, 210.351, -38.3094, 2.96344};
  
  double ax = 0.;

  for (int i=0; i<9; i++) {
    ax += p[i]*pow(log10(T),i);
  }

  const double k = pow(10., ax);

  return k;
}


double XMatAl5083 :: GetResistivity()
{
  const double rho = Lwf * fTemp / GetNistConductivity(fTemp);
  return rho;
}

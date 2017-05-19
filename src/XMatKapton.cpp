#include <cmath>
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XMatKapton.hpp"

using Quench::XQuenchLogger;

XMatKapton :: XMatKapton()
{}

double XMatKapton :: GetCapacity()
{
  double C = 0.;
  calcapacity(C);

  return C;
}

double XMatKapton :: GetConductivity()
{
  double k = 0.;
  calconductivity(k);

  return k;
}

void XMatKapton :: calconductivity(double &k)
{
  const int n = 8;
  const double p[n] = {5.73101, -39.5199, 79.9313, -83.8572, 50.9157,
                       -17.9835, 3.42413, -0.27133};

  double ax = 0.;

  if (fTemp<4.3)
    k = 0.0378 + 0.00161*fTemp;
  else {
    for (int i=0; i<n; i++)
      ax += p[i] * pow(log10(fTemp),i);
    //
    k = pow(10, ax);
  }
}

void XMatKapton :: calcapacity(double &C)
{
  // fitting parameters
  const int n = 8;
  const double p[n] = {-1.3684,  0.65892, 2.8719, 0.42651, -3.0088,
                        1.9558, -0.51998, 0.051574};

  double ax = 0.;
  for (int i=0; i<n; i++)
    ax += p[i] * pow(log10(fTemp),i);

  C = pow(10, ax);
}

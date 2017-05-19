#include <iostream>
#include <cmath>
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XMatCopper.hpp"

using Quench::XQuenchLogger;


XMatCopper :: XMatCopper()
    : fRhoRT(1.553e-8)
{}


XMatCopper :: ~XMatCopper()
{}


double XMatCopper :: GetCapacity()
{
  const double C = calcapacity(fTemp);

  return C;
}

double XMatCopper :: GetConductivity()
{
  const double k = calconduct(fTemp);

  if (fFld > 0.)
    calmagcdt(k);

  return k;
}

double XMatCopper :: GetResistivity()
{
  double res = calresist(fTemp);

  if ( fFld>0. )
    res = calmagresist(res, fFld);

  return res;
}

double XMatCopper :: calresist(const double T) const
{
  // fitting parameter
  const int n = 7;
  const double p[n] = {1.171e-17, 4.49, 3.841e+10, 1.14, 50.0,
                       6.428, 0.4531};

  // resistivity at room temperature [Ohm*m]
  const double rho0  = fRhoRT / fRRR;
  const double rhoi  = p[0] * pow(T,p[1]) / (1 + p[0] * p[2] * pow(T, p[1]-p[3]) * exp(-pow(p[4]/T,p[5])));
  const double rhoi0 = p[6] * rhoi * rho0 / (rho0 + rhoi);

  const double rho   = rho0 + rhoi + rhoi0;

  return rho;
}

double XMatCopper :: calmagresist(double rho, double B) const
{
  const int nb = 5;
  const double pb[nb] = {-2.662, 0.3168, 0.6229, -0.1839, 0.01827};

  const double ax = pb[0]*pow(log10(fRRR*B),0) + pb[1]*pow(log10(fRRR*B),1) 
                  + pb[2]*pow(log10(fRRR*B),2) + pb[3]*pow(log10(fRRR*B),3) 
                  + pb[4]*pow(log10(fRRR*B),4); 
  rho *= (1 + pow(10,ax));

  return rho;
}

double XMatCopper :: calconduct(const double T) const
{
  const double beta  = 0.634 / fRRR;
  const double betar = beta / 0.0003;

  // fitting parameter
  const int n = 7;
  const double p[n] = {1.754e-8, 2.763, 1102.0, -0.165, 70.0,
                       1.756, 0.838/pow(betar,0.1661)};

  //
  const double W0  = beta / T;
  const double Wi  = p[0]*pow(T,p[1]) / (1 + p[0]*p[2]*pow(T,(p[1]+p[3]))*exp(-pow(p[4]/T,p[5])));
  const double Wi0 = p[6] * Wi * W0 / (Wi + W0);
  const double k   = 1. / (W0 + Wi + Wi0);

  return k;
}

double XMatCopper :: calmagcdt(double k) const
{
  const double res = calresist(fTemp);
  if ( fFld>0. )
    calmagresist(res, fFld);

  const double zero = calresist(fTemp);
  k = zero * k / res;

  return k;
}

double XMatCopper :: calcapacity(const double T) const
{
  // fitting parameter
  const int n = 8;
  const double p[n] = {-1.91844, -0.15973, 8.61013, -18.996,
                        21.9661, -12.7328, 3.54322, -0.3797};

  //
  double ax = 0.;
  for (int i=0; i<n; i++)
    ax += p[i] * pow(log10(T), i);

  const double C = pow(10, ax);

  return C;
}

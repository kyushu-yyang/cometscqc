#include <cmath>
#ifndef IFdmUnits_HH
#include "IFdmUnits.hpp"
#endif
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XMatAluminium.hpp"

using Quench::XQuenchLogger;

XMatAluminium :: XMatAluminium()
    : fRhoRT(2.75e-8)
{}

XMatAluminium :: ~XMatAluminium()
{}

double XMatAluminium :: calresist() const
{
  double* par = new double[7];

  par[0] = 1.671e-17;
  par[1] = 4.36;
  par[2] = 2.841e+10;
  par[3] = 1.18;
  par[4] = 64.0;
  par[5] = 4.428;
  par[6] = 1.2031;

  const double res = evalresist(par);

  delete [] par;

  return res;
}


double XMatAluminium :: evalresist(double* par) const
{
  const double rho0  = fRhoRT / fRRR;
  const double rhoi  = par[0] * pow(fTemp,par[1]) / 
                       ( 1 + par[0]*par[2]*pow(fTemp,par[1]-par[3])*exp(-pow(par[4]/fTemp,par[5])) );
  const double rhoi0 = par[6] * rhoi * rho0 / (rho0 + rhoi);

  double rho = rho0 + rhoi + rhoi0;

  return rho;
}


double XMatAluminium :: calmagres(double res) const
{
  // fitting parameter for Al magnetoresistance
  const int nb = 5;
  const double pb[nb] = {3.62857, 2.90419e-5, 3.79649e+6, 10975.9, 0.761609};

  const double h    = fFld * 10. * fRhoRT / res;
  res += h*h*(pb[0] - pb[1]*h)*res / (pb[2] + pb[3]*h + pb[4]*h*h);

  return res;
}


double XMatAluminium :: GetResistivity()
{
  double res = calresist();

  if ( fFld>0. )
    res = calmagres(res);

  return res;
}


double XMatAluminium :: GetConductivity()
{
  const double res = GetResistivity();
  const double k   = Lwf * fTemp / res;

  return k;
}


double XMatAluminium :: GetCapacity()
{
  const double C = calcapacity(fTemp);

  return C;
}


double XMatAluminium :: calcapacity(const double T) const
{
  // switch point
  const double sp1 = 22.67;
  const double sp2 = 46.0;

  // fitting parameters
  const int n = 4;
  const double p1[n] = {-0.207489, 0.165759, -0.0142572, 0.00146459};
  const double p2[n] = {7.88e-13, 6.93201, -0.07139, 46.4363};
  const double p3[n] = {6.273517, -0.5469, 0.000925, -156.932};

  double C;

  if (T>0. && T<sp1)
    C = p1[0] + p1[1]*T + p1[2]*pow(T,2) + p1[3]*pow(T,3);
  else if (T>=sp1 && T<sp2)
    C = p2[0] * pow(T,p2[1]) * exp(p2[2]*T) * exp(p2[3]/T) * 4.186e+3;
  else if (T>=sp2)
    C = p3[0] * pow(T,p3[1]) * exp(p3[2]*T) * exp(p3[3]/T) * 4.186e+3;
  else {
    QuenchError( XQuenchLogger::ERROR, "Temperature: " << T << " is out of range for capacity of Al." );
    XQuenchExcept except("Temperature is out of range!");
    throw except;
  }

  return C;
}

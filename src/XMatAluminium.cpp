#include <cmath>
#ifndef IFdmUnits_HH
#include "IFdmUnits.hpp"
#endif
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XMatAluminium.hpp"

using Quench::XQuenchLogger;

XMatAluminium :: XMatAluminium()
    : fRhoRT(2.7e-8)
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

/*
  par[0] = 0.09052e-16;
  par[1] = 4.551;
  par[2] = 5.173e+10;
  par[3] = 1.26;
  par[4] = 40.;
  par[5] = 13.64;
  par[6] = 0.7416;
*/

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
  const double pb[nb] = {3.62857, -2.90419e-5, 3.79649e+6, 10975.9, 0.761609};
  //const double pb[nb] = {1., 0.00177, 1.8, 1.6, 0.53};

  const double rhoRT = 2.75e-8;
  const double h    = fFld*10. * rhoRT / res;
  res = res + pow(h,2)*(pb[0]+pb[1]*h)*res / (pb[2]+pb[3]*h+pb[4]*pow(h,2));

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

/*
double XMatAluminium :: kohler_plot(const double RRR, const double B) const
{
  // my fitting parameter
  //const double b[5] = {31.1133, 1.834e-4, 4.1053e+5, 6.3665e+3, 14.9814};
  const double b[5] = {8.7744e+2, 5.1726e-3, 1.1578e+7, 1.7954e+5, 4.225e+2};
  //const double b[5] = {1., 0.00177, 1.8, 1.6, 0.53};

  
  const double h     = B * RRR;
  const double Rb_R  = pow(h,2)*(b[0] - b[1]*h) / (b[2] + b[3]*h + b[4]*pow(h,2)) + 1.;

  const double RRR_eq = RRR / Rb_R;
  return RRR_eq;
}
*/

double XMatAluminium :: kohler_plot(const double rhoB0, const double B) const
{
  // my fitting parameter
  //const double b[5] = {31.1133, 1.834e-4, 4.1053e+5, 6.3665e+3, 14.9814};
  //const double b[5] = {8.7744e+2, 5.1726e-3, 1.1578e+7, 1.7954e+5, 4.225e+2};
  const double b[5] = {903.5, -0.0001985, 2.928e+7, 1.451e+5, 417.};
  //const double b[5] = {1., 0.00177, 1.8, 1.6, 0.53};

  const double h    = B * 2.70e-8 / rhoB0;
  const double drho = pow(h,2)*(b[0] - b[1]*h) / (b[2] + b[3]*h + b[4]*pow(h,2)) + 1.;

  const double rhoB = rhoB0 * drho;
  return rhoB;
}

double XMatAluminium :: hust_eq_resist(const double T, double RRR, const double B) const
{
  // hust's fitting parameter
  const double p[7] = {0.9052e-17, 4.551, 5.173e+10, 1.26, 40., 13.64, 0.7416};

  const double rhoRT = 2.70e-8;
  const double rho0  = rhoRT / RRR;
  const double rhoi  = p[0]*pow(T,p[1]) / (1 + p[0]*p[2]*pow(T,p[1]-p[3])*exp(-pow(p[4]/T,p[5])));
  const double rhoi0 = p[6] * rhoi * rho0 / (rho0 + rhoi);
  double rho   = rho0 + rhoi + rhoi0;

  if (B>0.)
    rho = kohler_plot(rho, B);

  return rho;
}

double XMatAluminium :: hust_eq_therm(const double T, double RRR, const double B) const
{
  // hust's fitting parameter
  // ref.: hust et al., national bureau of standards
  const double p[7]  = {4.716e-8, 2.446, 623.6, -0.16, 130.9, 2.5, 0.8168};

  //if (B>0.)
  //  RRR = kohler_plot(RRR, B);
  double scale = 1.;
  if (B>0.)
    scale = hust_eq_resist(T,RRR,0.) / hust_eq_resist(T,RRR,B);

  const double rhoRT = 2.70e-8;
  const double rho0  = rhoRT / RRR;
  const double beta  = rho0 / Lwf;

  const double W0  = beta / T;
  const double Wc  = -0.0005*log(T/330.)*exp(-pow(-log(T/380.)/0.6,2)) - 0.0013*log(T/110.)*exp(-pow(log(T/94.)/0.5,2));
  const double Wi  = p[0] * pow(T,p[1]) * pow(1 + p[0]*p[2]*pow(T,p[1]+p[3])*exp(-pow(p[4]/T,p[5])), -1.) + Wc;
  const double Wi0 = p[6] * Wi * W0 / (Wi+W0);

  const double k = scale * 1. / (W0+Wi+Wi0);
  return k;
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

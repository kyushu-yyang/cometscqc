#include <cmath>
#include "XQuenchLogger.hpp"
#include "XMatNbTi.hpp"

using Quench::XQuenchLogger;

XMatNbTi :: XMatNbTi()
    : fIc5(14.2e+3),
      fPar(3),
      XMaterial()
{}

void XMatNbTi :: SetIcAt5Tesla(const double Ic)
{
  fIc5 = Ic;
//  QuenchInfo( "Ic at 5T and 4.2K: " << fIc5 );
}

double XMatNbTi :: GetCapacity()
{
  double C = calcapacity(fTemp, fFld);
  return C;
}

double XMatNbTi :: GetCriticalI()
{
  // wilson's book page 3
  // Jc ~ 2.4e+9 A/m2 at 4.2K, 5T
  // 2400 A/mm2 * 4.73mm * 15mm -> I0
  //const double I0 = 2400.*4.73*15.;    // normalized factor
  // comet stabilized cable
  double Ic = calcriticalcurrent(fTemp, fFld, fIc5);
  return Ic;
}

double XMatNbTi :: GetCriticalT()
{
  double Tc = calTc(fFld);
  return Tc;
}


double XMatNbTi :: GetSharingT(const double I)
{
  const double Tc = GetCriticalT();
  const double Ic = GetCriticalI();

  if (Ic==0.)
    return 0.;

  const double Tcs = fTemp + (Tc - fTemp) * (1 - I/Ic);
  //const double Tcs = Tc - (Tc - T0) * I / Ic;

  return Tcs;
}


double XMatNbTi :: GetSharingT(const double I, const double T)
{
  const double Tc = GetCriticalT();
  const double Ic = GetCriticalI();

  if (Ic==0.)
    return 0.;

  const double Tcs = T + (Tc - T) * (1 - I/Ic);
  return Tcs; 
}


double XMatNbTi :: calcapacity(const double T, const double B) const
{
  // fitting parameter
  const int n = 5;
  const double p1[n] = {     0.0,   64.*B,    0.0,    49.1,      0.0};
  const double p2[n] = {     0.0,    928.,    0.0,   16.24,      0.0};
  const double p3[n] = {  41383., -7846.1, 553.71, 11.9838,  -0.2177};
  const double p4[n] = {-1.53e+6,  83022., -716.3,   2.976, -0.00482};
  const double p5[n] = { 1.24e+6,  13706., -51.66, 0.09296, -6.29e-5};
  const double p6[n] = { 2.45e+6,   955.5, -0.257,      0.,       0.};

  //const double Tc  = 9.4;
  const double Tc  = calTc(B);
  double C = 0.;

  if (T>0. && T<Tc) {
    for (int i=0; i<n; i++)
      C += p1[i] * pow(T, i); 
  }
  else if (T>=Tc && T<28.358) {
    for (int i=0; i<n; i++)
      C += p2[i] * pow(T, i);
  }
  else if (T>=28.358 && T<50.99) {
    for (int i=0; i<n; i++)
      C += p3[i] * pow(T, i);
  }
  else if (T>=50.99 && T<165.8) {
    for (int i=0; i<n; i++)
      C += p4[i] * pow(T, i);
  }
  else if (T>=165.8 && T<496.54) {
    for (int i=0; i<n; i++)
      C += p5[i] * pow(T, i);
  }
  else if (T>=496.54) {
    for (int i=0; i<n; i++)
      C += p6[i] * pow(T, i);
  }

  return C / GetDensity();
}


void XMatNbTi :: SetIcParameter(const int par)
{
  fPar = par;
  double Tc0, Bc20, C0, alpha, beta, gamma;
  GetIcPar(Tc0, Bc20, C0, alpha, beta, gamma);

  QuenchInfo( "Ic parameter -> Tc: " << Tc0 << " , Bc20: " << Bc20 << ", C0: " << C0 << 
              " , alpha: " << alpha << " , beta: " << beta << " , gamma: " << gamma );
}


void XMatNbTi :: GetIcPar(double &Tc0, double &Bc20, double &C0, double &alpha, double &beta, double &gamma)
{
  // fitting equation from L. Bottura's paper
  const int    n = 5;
  const double T[n] = { 9.2,  8.5,  8.9,  9.2,  9.35};
  const double B[n] = {14.5, 14.2, 14.4, 14.4, 14.25};
  const double C[n] = {23.8, 28.6, 28.5, 37.7,  28.4};
  const double a[n] = {0.57, 0.76, 0.64, 0.89,  0.80};
  const double b[n] = {0.90, 0.85, 0.75, 1.10,  0.89};
  const double g[n] = {1.90, 1.76, 2.30, 2.09,  1.87};

  Tc0   = T[fPar];
  Bc20  = B[fPar];
  C0    = C[fPar];
  alpha = a[fPar];
  beta  = b[fPar];
  gamma = g[fPar];
}


double XMatNbTi :: calcriticalcurrent(const double T, const double B, const double I0) const
{
  // fitting equation from L. Bottura's paper
  const double n = 1.7;
  const double Tc0 [5] = { 9.2,  8.5,  8.9,  9.2,  9.35};
  const double Bc20[5] = {14.5, 14.2, 14.4, 14.4, 14.25};
  const double C0  [5] = {23.8, 28.6, 28.5, 37.7,  28.4};
  const double alp [5] = {0.57, 0.76, 0.64, 0.89,  0.80};
  const double beta[5] = {0.90, 0.85, 0.75, 1.10,  0.89};
  const double gam [5] = {1.90, 1.76, 2.30, 2.09,  1.87};

  const double t   = T / Tc0[fPar];
  const double Bc2 = Bc20[fPar] * (1 - pow(t,n));
  const double b   = B / Bc2;

  if (b>1)
    return 0.;

  // normalized critical current density
  double Ic = C0[fPar] * pow(b,alp[fPar]) * pow((1-b),beta[fPar]) * pow((1-pow(t,n)),gam[fPar]) / B;
  Ic *= I0;

  return Ic;
}


double XMatNbTi :: calTc(const double B) const
{
  const double Tc0 [5] = { 9.2,  8.5,  8.9,  9.2,  9.35};
  const double Bc20[5] = {14.5, 14.2, 14.4, 14.4, 14.25};

  const double n = 1.7;
  double Tc = Tc0[fPar] * pow((1 - B/Bc20[fPar]), 1/n);

  return Tc;
}

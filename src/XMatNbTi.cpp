#include <cmath>
#include "XMatNbTi.hpp"

XMatNbTi :: XMatNbTi()
{}

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
  const double I0 = 2400.*4.73*15.;    // normalized factor
  double Ic = calcriticalcurrent(fTemp, fFld, I0);

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
  const double T0 = 4.5;
  const double Ic = GetCriticalI();
  //const double Tcs = Tc + (Tc - T0) * (1 - I/Ic);
  const double Tcs = Tc - (Tc - T0) * I / Ic;

  return Tcs;
}


double XMatNbTi :: calcapacity(const double T, const double B) const
{
  // fitting parameter
  const int n = 5;
  const double p1[n] = {0.0, 64.*B, 0.0, 49.1, 0.0};
  const double p2[n] = {0.0, 928., 0.0, 16.24, 0.0};
  const double p3[n] = {41383., -7846.1, 553.71, 11.9838, -0.2177};
  const double p4[n] = {-1.53e+6, 83022., -716.3, 2.976, -0.00482};
  const double p5[n] = {1.24e+6, 13706., -51.66, 0.09296, -6.29e-5};
  const double p6[n] = {2.45e+6, 955.5, -0.257, 0., 0.};

  const double Tc  = 9.4;
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


double XMatNbTi :: calcriticalcurrent(const double T, const double B, const double I0) const
{
  // fitting equation from L. Bottura's paper
  const int ni = 5;
  const double n = 1.7;
  const double Tc0 [ni] = { 9.2,  8.5,  8.9,  9.2,  9.35};
  const double Bc20[ni] = {14.5, 14.2, 14.4, 14.4, 14.25};
  const double C0  [ni] = {23.8, 28.6, 28.5, 37.7,  28.4};
  const double alp [ni] = {0.57, 0.76, 0.64, 0.89,  0.80};
  const double beta[ni] = {0.90, 0.85, 0.75, 1.10,  0.89};
  const double gam [ni] = {1.90, 1.76, 2.30, 2.09,  1.87};

  const int    m = 0;
  const double t   = T / Tc0[m];
  const double Bc2 = Bc20[m] * (1 - pow(t,n));
  const double b   = B / Bc2;

  // normalized critical current density
  double Ic = C0[m] * pow(b,alp[m]) * pow((1-b),beta[m]) * pow((1-pow(t,n)),gam[m]) / B;
  Ic *= I0;

  return Ic;
}


double XMatNbTi :: calTc(const double B) const
{
  const int m = 5;
  double Tc0 [m] = { 9.2,  8.5,  8.9,  9.2,  9.35};
  double Bc20[m] = {14.5, 14.2, 14.4, 14.4, 14.25};

  const int data = 0;
  const double n = 1.7;
  double Tc = Tc0[data] * pow((1 - B/Bc20[data]), 1/n);

  return Tc;
}

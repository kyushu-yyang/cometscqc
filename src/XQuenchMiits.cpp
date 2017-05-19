#include <iostream>
#include <TSpline.h>
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XMatAluminium.hpp"
#include "XMatCopper.hpp"
#include "XMatNbTi.hpp"
#include "XQuenchMiits.hpp"

using Quench::XQuenchLogger;

XQuenchMiits :: XQuenchMiits()
    : fResist(0.182),
      fIndct(12.69),
      fCurr(2700.),
      fField(5.5),
      fRRR(100.),
      fTemp0(4.5),
      fTempf(400.),
      fCoil(NULL)
{}


XQuenchMiits :: ~XQuenchMiits()
{
  if (fCoil) delete fCoil;
}

void XQuenchMiits :: SetDumpResistor(const double R)
{
  fResist = R;

  QuenchInfo("resistance of dump resistor: " << fResist);
}


void XQuenchMiits :: SetInductance(const double L)
{
  fIndct = L;

  QuenchInfo("magnet total inductance: " << fIndct);
}


double XQuenchMiits :: GetTimeConstance() const
{
  const double tau = fIndct / fResist;
  QuenchInfo("time constance: " << tau << " sec");

  return tau;
}


void XQuenchMiits :: SetCurrent(const double I0)
{
  fCurr = I0;

  QuenchInfo("intial current: " << fCurr << " A");
}


void XQuenchMiits :: SetField(const double B)
{
  fField = B;

  QuenchInfo("(MIITs) setting field: " << fField << " Tesla");
}


void XQuenchMiits :: SetRRR(const double RRR)
{
  fRRR = RRR;

  QuenchInfo("(MIITs) setting RRR: " << fRRR);
}


void XQuenchMiits :: SetInitialTemp(const double T0)
{
  fTemp0 = T0;

  QuenchInfo("initial temperature: " << fTemp0 << " K");
}


void XQuenchMiits :: SetFinalTemp(const double Tf)
{
  fTempf = Tf;

  QuenchInfo("final temperature: " << fTempf << " K");
}


void XQuenchMiits :: SetTempRange(const double T0, const double Tf)
{
  SetInitialTemp(T0);
  SetFinalTemp(Tf);
}


double XQuenchMiits :: GetMiits() const
{
  //const double Atot = fCoil->GetCoilParts(kConductor)->GetArea();
  //const double J0 = fCurr / Atot;
  const double I0 = fCurr;
  //const double miits = fIndct * pow(J0, 2) / fResist / 2.;
  const double miits = fIndct * pow(I0, 2) / fResist / 2.;

  //QuenchInfo("initial current density: " << J0 << " A/m2");
  QuenchInfo("initial current: " << I0 << " A");
  QuenchInfo("Miits: " << miits);

  return miits;
}


std::map<double,double> XQuenchMiits :: GetMiitsDecay(const double t0, 
                                                      const double tf, const double T)
{
  const int    nt = 800;
  const double dt = (tf - t0) / nt;

  //const double Atot = fCoil->GetCoilParts(kConductor)->GetArea();
  //const double J0 = fCurr / Atot;
  const double I0 = fCurr;

  double t = t0;
  double miits = 0.;
  double res = fResist;
  const double v = GetVelocity(fCurr);

  std::map<double,double> decay;

  while (t<=tf) {
    res = fResist + GetAvgResistance(T, fRRR, fField, v*t);
    //res = fResist;
    miits += pow(I0,2) * exp(-2 * res * t / fIndct) * dt;
    decay.insert( std::map<double,double>::value_type(t,miits) );

    t += dt;
  }

  return decay;
}


std::map<double, double> XQuenchMiits :: Eval(const double T0, const double Tf)
{
  const int    nT = 1500;
  const double dT = (Tf - T0) / nT;

  double T = T0;
  std::map<double, double> miits;

  const double Atot = fCoil->GetCoilParts(kConductor)->GetArea();
  const double rho_avg = 4000.;
  double C_avg = 0.;
  double res = 0.;

  double mii = 0.;

  while (T<=Tf) {
    C_avg = GetAvgCapacity(T, fRRR, fField);
    res = GetAvgResistance(T, fRRR, fField);

    mii += Atot * rho_avg * C_avg * dT / res;

    miits.insert( std::map<double, double>::value_type(T,mii) );

    T += dT;
  }

  return miits;
}


double XQuenchMiits :: interpolate(std::map<double, double> hc, const double x) const
{
  double x0 = 0.;
  double dx = 0.;

  // find t0 and dt
  int cnt = 0;
  for (std::map<double, double>::iterator it=hc.begin(); it!=hc.end(); it++) {
    if (cnt==0)
      x0 = it->first;
    if (cnt==1) {
      dx = it->first - x0;
      break;
    }
    cnt++;
  }

  // find data x position
  const int entry = (x - x0) / dx;

  double xa = 0.;
  double xb = 0.;
  double ya = 0.;
  double yb = 0.;
  cnt = 0;

  for (std::map<double, double>::iterator it=hc.begin(); it!=hc.end(); it++) {
    if (cnt==entry-1) {
      xa = it->first;
      ya = it->second;
    }
    if (cnt==entry+1) {
      xb = it->first;
      yb = it->second;
      break;
    }
    cnt++;
  }

  // at point (x, y)
  const double y = ya + (yb-ya)*(x-xa)/(xb-xa);

  return y;
}


double XQuenchMiits :: CalMiits2(const double Temp)
{
  std::map<double, double> miit = Eval(0.1, fTempf); 
  /*
  int cnt = 0;

  const int n = miit.size();
  double T[n];
  double M[n];

  for (std::map<double, double>::iterator it=miit.begin(); it!=miit.end(); it++) {
    M[cnt] = it->second;
    T[cnt] = it->first;
    cnt++;
  }

  TSpline3 sp("miits2", T, M, n);
  */
  //double miits = sp.Eval(Temp);
  double miits = interpolate(miit, Temp);

  return miits;
}


double XQuenchMiits :: GetMaxTemperature(const double miits)
{
  std::map<double, double> miit = Eval(fTemp0, fTempf); 
  int cnt = 0;

  const int n = miit.size();
  double T[n];
  double M[n];

  for (std::map<double, double>::iterator it=miit.begin(); it!=miit.end(); it++) {
    M[cnt] = it->second;
    T[cnt] = it->first;
    cnt++;
  }

  TSpline3 sp("miits", M, T, n);
  double maxtemp = sp.Eval(miits);

  return maxtemp;
}


double XQuenchMiits :: GetMPZ() const
{
  XMatNbTi nbti;
  XMatAluminium al;
  //XMatCopper cu;

  nbti.SetField(fField);
  const double Tc = nbti.GetCriticalT();
  const double Atot = fCoil->GetCoilParts(kConductor)->GetArea();
  // take the average temperature for setting the temperature property for resistivity etc.
  //const double Tavg = (Tc + fTemp0) / 2.;
  const double Tavg = fTemp0;

  //cu.SetMaterialProperty(Tavg, fRRR, fField);
  al.SetMaterialProperty(Tavg, fRRR, fField);
  nbti.SetTemperature(Tavg);

  const double k = al.GetConductivity();
  const double rho = GetAvgResistance(Tavg, fRRR, fField) * Atot;
  //const double rho = Lwf * Tavg / k;
  //double Jc = nbti.GetCriticalI() / Atot;
  // convert to the effective critical current
  //Jc *= fCoil->GetMaterialRatio(iNbTi);
  const double Jc = fCurr / Atot;

  const double mpz = sqrt(2. * k * (Tc - fTemp0) / pow(Jc,2) / rho);

  //std::cout << Tc << "   " << k/rho << "  " << Jc << std::endl;
  QuenchInfo("Jc: " << Jc << " A/m2, Tc: " << Tc << " K, " <<
             "Tavg: " << Tavg << " K, k:" << k << " W/m/K, " <<
             "rho: " << rho << " Ohm*m, " <<
             "MPZ: " << mpz << " m");

  return mpz;
}


std::map<double,double> XQuenchMiits :: GetTimeTemp(const double T0, const double Tf)
{
  const int    nT = 500;
  const double dT = (Tf - T0) / nT;

  std::map<double, double> temp = Eval(T0, Tf+5.);
  std::map<double, double> curr;
  std::map<double, double> collect;

  double miit = 0.;
  double T = 0.;
  double time = 0.;

  for (unsigned int i=0; i<nT; i++) {
    T = T0 + i*dT;
    curr = GetMiitsDecay(0., 200., T);
    miit = interpolate(temp, T);
    time = findtime(curr, miit);
    collect.insert( std::map<double,double>::value_type(time, T) );
  }

  return collect;
}


std::map<double, double> XQuenchMiits :: GetTimeCurrent(const double T0, const double Tf)
{
  std::map<double,double> timetemp = GetTimeTemp(T0, Tf);
  std::map<double,double> timecurr;

  const double I0 = fCurr;
  double I = I0;
  const double v = GetVelocity(I0);
  double T = 0.;
  double t = 0.;

  for (std::map<double,double>::iterator it=timetemp.begin(); it!=timetemp.end(); it++) {
    t = it->first;
    T = it->second;
    I = I0 * exp(- (fResist + GetAvgResistance(T,fRRR,fField,v*t))*t/fIndct);
    timecurr.insert( std::map<double,double>::value_type(t,I) );
  }

  return timecurr;
}


std::map<double, double> XQuenchMiits :: GetTimeResist(const double T0, const double Tf)
{
  std::map<double, double> timetemp = GetTimeTemp(T0, Tf);
  std::map<double, double> timeres;

  double R = 0.;
  double t = 0.;
  double T = 0.;
  const double v = GetVelocity(fCurr);

  for (std::map<double,double>::iterator it=timetemp.begin(); it!=timetemp.end(); it++) {
    t = it->first;
    T = it->second;
    R = GetAvgResistance(T,fRRR,fField,v*t);
    timeres.insert( std::map<double,double>::value_type(t,R) );
  }

  std::cout << "Estimated quench velocity: " << v << std::endl;

  return timeres;
}


std::map<double, double> XQuenchMiits :: GetTimeVolt(const double T0, const double Tf)
{
  std::map<double, double> timecurr = GetTimeCurrent(T0, Tf);
  std::map<double, double> timeres  = GetTimeResist(T0, Tf);
  std::map<double, double> timevolt;

  double t = 0.;
  double I = 0.;
  double R = 0.;
  double V = 0.;

  for (std::map<double,double>::iterator it=timecurr.begin(); it!=timecurr.end(); it++) {
    t = it->first;
    I = it->second;
    R = timeres[t];
    V = I*R;
    timevolt.insert( std::map<double,double>::value_type(t,V) );
  }

  return timevolt;
}


double XQuenchMiits :: GetAvgCapacity(const double T, const double RRR, const double B)
{
  XMatNbTi      sc;
  XMatCopper    cu;
  XMatAluminium al;

  sc.SetMaterialProperty(T, RRR, B);
  cu.SetMaterialProperty(T, 50., B);
  al.SetMaterialProperty(T, RRR, B);

  const double r_Al = fCoil->GetMaterialRatio(iAluminium);
  const double r_Cu = fCoil->GetMaterialRatio(iCopper);
  const double r_SC = fCoil->GetMaterialRatio(iNbTi);

  const double rho_Cu = cu.GetDensity();
  const double rho_Al = al.GetDensity();
  const double rho_SC = sc.GetDensity();
  const double rho_avg = 4000.;

  const double C_Al = r_Al * rho_Al * al.GetCapacity() / rho_avg;
  const double C_Cu = r_Cu * rho_Cu * cu.GetCapacity() / rho_avg;
  const double C_SC = r_SC * rho_SC * sc.GetCapacity() / rho_avg;
  const double C_avg = C_Al + C_Cu + C_SC;

  return C_avg;
}


double XQuenchMiits :: GetAvgResistance(double T, double RRR, double B, double l) const
{
  XMatNbTi      sc;
  XMatCopper    cu;
  XMatAluminium al;

  sc.SetMaterialProperty(T, RRR, B);
  cu.SetMaterialProperty(T, 50., B);
  al.SetMaterialProperty(T, RRR, B);

  const double r_Al = fCoil->GetMaterialRatio(iAluminium);
  const double r_Cu = fCoil->GetMaterialRatio(iCopper);

  const double Atot = fCoil->GetCoilParts(kConductor)->GetArea();
  const double A_Al = r_Al * Atot;
  const double A_Cu = r_Cu * Atot;

  const double R_Al  = al.GetResistivity() * l / A_Al;
  const double R_Cu  = cu.GetResistivity() * l / A_Cu;
  const double R_avg = pow( (1./R_Al + 1./R_Cu), -1. );
  //const double res   = Atot * R_avg / l;
  //
  //std::cout << RRR << " " << al.GetResistivity();
  //std::cout << " " << cu.GetResistivity() << std::endl;

  //return res;
  return R_avg;
}


double XQuenchMiits :: findtime(std::map<double,double> time, double miits)
{
  std::map<double, double> collect;

  double m = 0.;
  int cnt = 0;
  for (std::map<double,double>::iterator it=time.begin(); it!=time.end(); it++) {
    if (cnt==time.size()-1) {
      m = it->second;
      break;
    }
    cnt++;
  }

  if (miits > m) {
    return 0.;
  }

  for (std::map<double,double>::iterator it=time.begin(); it!=time.end(); it++) {
    collect.insert( std::map<double,double>::value_type(it->second, it->first) );
  }

  return interpolate(collect, miits);
}


double XQuenchMiits :: GetVelocity(const double Iop)
{
  XMatNbTi sc;
  sc.SetMaterialProperty(fTemp0, fRRR, fField);

  const double Tcs = sc.GetSharingT(Iop);
  
  const double Atot = fCoil->GetCoilParts(kConductor)->GetArea();
  const double Jop = Iop / Atot;

  const double rho_avg = 4000.;
  const double C_avg = GetAvgCapacity((Tcs+fTemp0)/2, fRRR, fField);
  //const double C_avg = 0.8;

  const double vqch = Jop * sqrt(Lwf * Tcs / (Tcs-fTemp0)) / rho_avg / C_avg;
  //std::cout << C_avg << " " << Tcs << std::endl;

  return vqch;
}


double XQuenchMiits :: GetQuenchRes(const double t)
{
  const double T0    = fTemp0;
  const double Atot  = fCoil->GetCoilParts(kConductor)->GetArea();
  const double rho0  = GetAvgResistance(T0, fRRR, fField)/Atot;
  const double alpha = 1.;
  const double J0    = fCurr / Atot;
  const double v     = GetVelocity(fCurr);
  //const double U0    = CalMiits2(T0);
  const double U0    = 1e+12;

  double R = 4*M_PI*rho0*pow(alpha,2)*pow(J0,4)*pow(v,3)*pow(t,5)/(30.*pow(Atot,2)*pow(U0,2));

  return R;
}

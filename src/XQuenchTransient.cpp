#include <iostream>
#include "XQuenchLogger.hpp"
#include "XQuenchExcept.hpp"
#include "XQuenchTransient.hpp"

using Quench::XQuenchLogger;

XQuenchTransient :: XQuenchTransient()
    : fDumpRes(0.125),
      fInduct(12.69),
      fPreI(2700.),
      fCurr(2700.*Amp),
      fVth(0.1),
      fDetTime(1.*sec),
      fDiode(0.7),
      fHotZ(1),
      fHotPhi(1),
      fHotR(1),
      fHotSpotHeat(5000.*100.),
      fMagConnect(NULL),
      fShellConnect(NULL),
      XTransientLoop()
{}


XQuenchTransient :: ~XQuenchTransient()
{
  if (fMagConnect) delete [] fMagConnect;
  if (fShellConnect) delete [] fShellConnect;
}


void XQuenchTransient :: SetHotSpot(const int z, const int phi, const int r, const double q)
{
  fHotZ = z;
  fHotPhi = phi;
  fHotR = r;
  fHotSpotHeat = q;

  QuenchInfo(" set hot spot on (" << z << ", " << phi << ", " << r << ")");
}


void XQuenchTransient :: SetDumpResistor(const double R)
{
  fDumpRes = R;
  QuenchInfo(" set dump resistor resistance: " << fDumpRes << " [Ohm]");
}


void XQuenchTransient :: SetInductance(const double L)
{
  fInduct = L;
  QuenchInfo(" set coil inductance: " << fInduct << " [H]");
}


void XQuenchTransient :: SetCurrent(const double I)
{
  fCurr = I;
  fPreI = I;
  QuenchInfo(" set operating current: " << fCurr << " [A]");
}


void XQuenchTransient :: SetVoltage(const double V)
{
  fDumpRes = V / fCurr;
  QuenchInfo(" calculating the dump resistance: " << fDumpRes << " [Ohm]"
             << " (I=" << fCurr << "[A], V=" << V << "[V]");
}


void XQuenchTransient :: SetThreshold(const double Vth)
{
  fVth = Vth;
  QuenchInfo(" set threshold voltage: " << fVth << " [V]");
}


void XQuenchTransient :: SetDetectTime(const double time)
{
  fDetTime = time;
  QuenchInfo(" set detection time: " << fDetTime << " [sec]");
}


void XQuenchTransient :: SetDiode(const double V)
{
  fDiode = V;
  QuenchInfo(" set diode turn on voltage: " << fDiode << " [V]");
}


double XQuenchTransient :: GetDiodeVoltage(const double I) const 
{
  const double Is = 1e-12;
  const double Vt = 0.026;

  const double V = Vt * log( I/Is + 1 );
  return V;
}


double XQuenchTransient :: CalCurrentDecay(const double preI, const double res, const double dt)
{
  fDiode = GetDiodeVoltage(preI);

  double I = (preI*fInduct - fDiode*dt) / (fInduct + dt*(fDumpRes+res));
  return I;
}


void XQuenchTransient :: CalFieldDecay(XThermalSolver* solver)
{
  double preB = 0.;
  double dI = fCurr - fPreI;

  for (unsigned int i=0; i<solver->GetProcess()->GetMaterialEntries(); i++) {
    preB = solver->GetProcess()->GetMaterialEntry(i)->GetField();
    if (preB>0.) {
      solver->GetProcess()->GetMaterialEntry(i)->SetField( preB*(1+dI/fPreI) );
    }
    else {
      solver->GetProcess()->GetMaterialEntry(i)->SetField( 0. );
    }
  }
}


int XQuenchTransient :: GetTotalConductor(XThermalSolver* solver)
{
  int idz = 0;
  int idp = 0;
  int idr = 0;

  const int mshz = solver->GetProcess()->GetMesh(iZ);
  const int mshp = solver->GetProcess()->GetMesh(iPhi);
  const int mshr = solver->GetProcess()->GetMesh(iR);

  int cnt = 0;

  for (unsigned int i=0; i<solver->GetProcess()->GetMaterialEntries(); i++) {
    idz = solver->GetProcess()->GetDimensionEntry(i)->GetId(iZ);
    idp = solver->GetProcess()->GetDimensionEntry(i)->GetId(iPhi);
    idr = solver->GetProcess()->GetDimensionEntry(i)->GetId(iR);
    
    if (solver->GetProcess()->GetDimensionEntry(i)->GetGeometry()==kConductor &&
        idz>0 && idz<mshz+1 &&
        idp>0 && idp<mshp+1 &&
        idr>0 && idr<mshr+1)
      cnt ++;
  }

  return cnt;
}


int XQuenchTransient :: GetQuenchConductor(XThermalSolver* solver, QuenchStatus qch)
{
  int idz = 0;
  int idp = 0;
  int idr = 0;

  const int mshz = solver->GetProcess()->GetMesh(iZ);
  const int mshp = solver->GetProcess()->GetMesh(iPhi);
  const int mshr = solver->GetProcess()->GetMesh(iR);

  int cnt = 0;

  for (unsigned int i=0; i<solver->GetProcess()->GetMaterialEntries(); i++) {
    idz = solver->GetProcess()->GetDimensionEntry(i)->GetId(iZ);
    idp = solver->GetProcess()->GetDimensionEntry(i)->GetId(iPhi);
    idr = solver->GetProcess()->GetDimensionEntry(i)->GetId(iR);
    
    if (solver->GetProcess()->GetDimensionEntry(i)->GetGeometry()==kConductor &&
        solver->GetProcess()->GetMaterialEntry(i)->GetStatus()==qch &&
        idz>0 && idz<mshz+1 &&
        idp>0 && idp<mshp+1 &&
        idr>0 && idr<mshr+1)
      cnt ++;
  }

  return cnt;
}


void XQuenchTransient :: Begin()
{
  std::cout << "running the quench transient loop ... " << std::endl;
}


void XQuenchTransient :: Run()
{

}


void XQuenchTransient :: End()
{

}

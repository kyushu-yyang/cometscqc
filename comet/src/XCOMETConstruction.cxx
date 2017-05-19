#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <TString.h>
#include "XMatCopper.hpp"
#include "XMatAluminium.hpp"
#include "XMatNbTi.hpp"

#include "XCoilConductor.hpp"
#include "XCoilStrip.hpp"
#include "XCoilShell.hpp"

#include "XRadiationHandle.hpp"
#include "XRootOutput.hpp"
#include "XQuenchOutput.hpp"
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XCOMETConstruction.h"


XCOMETConstruction :: XCOMETConstruction()
    : fDay(90*day),
      fFld(NULL),
      fCS0(NULL),
      fCS1(NULL),
      fMS1(NULL),
      fMS2(NULL),
      XQuenchTransient()
{
  if (!fFld)  fFld = new XFieldHandle();
  ConstructField(fFld);
}


XCOMETConstruction :: ~XCOMETConstruction()
{
  if (fFld)  delete fFld;
  if (fCS0)  delete fCS0;
  if (fCS1)  delete fCS1;
  if (fMS1)  delete fMS1;
  if (fMS2)  delete fMS2;
}


void XCOMETConstruction :: ConstructCS0(const std::string& radfile, const char* tempfile)
{
  XRadiationHandle* rad = new XRadiationHandle(radfile);
  rad->SetIrrTime(fDay);

  XInitialTemperature* Tinit = new XInitialTemperature;
  Tinit->Load(tempfile);

  const double r = (672.+823.65)/2.*mm;
  const std::string name = "CS0";
  
  XCoilHandle* coil = new XCoilHandle();
  coil->SetName(name);
  coil->SetCoilSize(0., 2.*M_PI*r, 0.);
  coil->SetMesh(35, 4, 19);
  coil->SetCoilLayers(9);
  coil->SetCoilTurns(35);
  coil->SetMaterialRatio(7.3, 1, 0.9);

  // set coil structure
  coil->AddLayer(1, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(2, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(3, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(4, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(5, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(6, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(7, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(8, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(9, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(10, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(11, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(12, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(13, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(14, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(15, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(16, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(17, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(18, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(19, kShell, GetShell(), kAdiabatic, 0.*cm);

  // set processor
  XProcessManager* pro = new XProcessManager();
  pro->SetCoilHandler(coil);
  pro->Initialize();
  
  // get magnetic field map
  fFld->SetTarget(name);
  fFld->Run();
  pro->SetFieldHandler(fFld);
  pro->SetRadiationHandler(rad);
  pro->SetInitialTemperature(Tinit);
  
  // uniform RRR and magnetic field
  //pro->SetUniformField(5.0);
  //pro->SetUniformRRR(kConductor, 113.);
  //pro->SetUniformRRR(kStrip, 122.);
  //pro->SetUniformRRR(kShell, 1.);

  // fill the materail info into the coil solver
  if (!fCS0) fCS0 = new XThermalSolver();
  fCS0->SetProcessHandle(pro);

  // write the geometry
  XQuenchOutput* geo = new XQuenchOutput("geoCS0.dat", iOfstream);
  geo->WriteGeometry(pro);
  geo->Close();

  std::cout << "constructing " << name << "..." << std::endl;
}


void XCOMETConstruction ::ConstructCS1(const std::string& radfile, const char* tempfile)
{
  XRadiationHandle* rad = new XRadiationHandle(radfile);
  rad->SetIrrTime(fDay);

  XInitialTemperature* Tinit = new XInitialTemperature;
  Tinit->Load(tempfile);

  const double r = (672.+823.65)/2.*mm;
  const std::string name = "CS1";
  
  XCoilHandle* coil = new XCoilHandle();
  coil->SetName(name);
  coil->SetCoilSize(0., 2.*M_PI*r, 0.);
  coil->SetMesh(45, 4, 19);
  coil->SetCoilLayers(9);
  coil->SetCoilTurns(270);
  coil->SetMaterialRatio(7.3, 1, 0.9);

  // set coil structure
  coil->AddLayer(1, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(2, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(3, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(4, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(5, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(6, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(7, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(8, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(9, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(10, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(11, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(12, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(13, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(14, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(15, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(16, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(17, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(18, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(19, kShell, GetShell(), kAdiabatic, 0.*cm);

  // set processor
  XProcessManager* pro = new XProcessManager();
  pro->SetCoilHandler(coil);
  pro->Initialize();
  
  // get magnetic field map
  fFld->SetTarget(name);
  fFld->Run();
  pro->SetFieldHandler(fFld);
  pro->SetRadiationHandler(rad);
  pro->SetInitialTemperature(Tinit);
  
  // uniform RRR and magnetic field
  //pro->SetUniformField(5.0);
  //pro->SetUniformRRR(kConductor, 209.);
  //pro->SetUniformRRR(kStrip, 241.);
  //pro->SetUniformRRR(kShell, 1.);

  // fill the materail info into the coil solver
  if (!fCS1) fCS1 = new XThermalSolver();
  fCS1->SetProcessHandle(pro);

  // write the geometry
  XQuenchOutput* geo = new XQuenchOutput("geoCS1.dat", iOfstream);
  geo->WriteGeometry(pro);
  geo->Close();

  std::cout << "constructing " << name << "..." << std::endl;
}


void XCOMETConstruction ::ConstructMS1(const std::string& radfile, const char* tempfile)
{
  XRadiationHandle* rad = new XRadiationHandle(radfile);
  rad->SetIrrTime(fDay);

  XInitialTemperature* Tinit = new XInitialTemperature;
  Tinit->Load(tempfile);

  const double r = (672.+756.25)/2.*mm;
  const std::string name = "MS1";
  
  XCoilHandle* coil = new XCoilHandle();
  coil->SetName(name);
  coil->SetCoilSize(0., 2.*M_PI*r, 0.);
  coil->SetMesh(45, 4, 11);
  coil->SetCoilLayers(5);
  coil->SetCoilTurns(285);
  coil->SetMaterialRatio(7.3, 1, 0.9);

  // set coil structure
  coil->AddLayer(1, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(2, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(3, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(4, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(5, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(6, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(7, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(8, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(9, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(10, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(11, kShell, GetShell(), kAdiabatic, 0.*cm);

  // set processor
  XProcessManager* pro = new XProcessManager();
  pro->SetCoilHandler(coil);
  pro->Initialize();
  
  // get magnetic field map
  fFld->SetTarget(name);
  fFld->Run();
  pro->SetFieldHandler(fFld);
  pro->SetRadiationHandler(rad);
  pro->SetInitialTemperature(Tinit);
  
  // uniform RRR and magnetic field
  //pro->SetUniformField(4.0);
  //pro->SetUniformRRR(kConductor, 250.);
  //pro->SetUniformRRR(kStrip, 297.);
  //pro->SetUniformRRR(kShell, 1.);

  // fill the materail info into the coil solver
  if (!fMS1) fMS1 = new XThermalSolver();
  fMS1->SetProcessHandle(pro);

  // write the geometry
  XQuenchOutput* geo = new XQuenchOutput("geoMS1.dat", iOfstream);
  geo->WriteGeometry(pro);
  geo->Close();

  std::cout << "constructing " << name << "..." << std::endl;
}


void XCOMETConstruction ::ConstructMS2(const std::string& radfile, const char* tempfile)
{
  XRadiationHandle* rad = new XRadiationHandle(radfile);
  rad->SetIrrTime(fDay);

  XInitialTemperature* Tinit = new XInitialTemperature;
  Tinit->Load(tempfile);

  const double r = (672.+789.95)/2.*mm;
  const std::string name = "MS2";
  
  XCoilHandle* coil = new XCoilHandle();
  coil->SetName(name);
  coil->SetCoilSize(0., 2.*M_PI*r, 0.);
  coil->SetMesh(50, 4, 15);
  coil->SetCoilLayers(7);
  coil->SetCoilTurns(140);
  coil->SetMaterialRatio(7.3, 1, 0.9);

  // set coil structure
  coil->AddLayer(1, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(2, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(3, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(4, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(5, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(6, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(7, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(8, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(9, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(10, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(11, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(12, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(13, kStrip, GetStrip(), kAdiabatic, 20.*cm);
  coil->AddLayer(14, kConductor, GetConductor(), kAdiabatic, 0.*cm);
  coil->AddLayer(15, kShell, GetShell(), kAdiabatic, 0.*cm);

  // set processor
  XProcessManager* pro = new XProcessManager();
  pro->SetCoilHandler(coil);
  pro->Initialize();
  
  // get magnetic field map
  fFld->SetTarget(name);
  fFld->Run();
  pro->SetFieldHandler(fFld);
  pro->SetRadiationHandler(rad);
  pro->SetInitialTemperature(Tinit);
  
  // uniform RRR and magnetic field
  //pro->SetUniformField(3.5);
  //pro->SetUniformRRR(kConductor, 400.);
  //pro->SetUniformRRR(kStrip, 916.);
  //pro->SetUniformRRR(kShell, 1.);

  // fill the materail info into the coil solver
  if (!fMS2) fMS2 = new XThermalSolver();
  fMS2->SetProcessHandle(pro);

  // write the geometry
  XQuenchOutput* geo = new XQuenchOutput("geoMS2.dat", iOfstream);
  geo->WriteGeometry(pro);
  geo->Close();

  std::cout << "constructing " << name << "..." << std::endl;
}


void XCOMETConstruction :: ConstructField(XFieldHandle* fld)
{
  fld->SetCurrent(2700.*Amp);
  fld->AddCoil( "CS0", 857.88*mm, 1038.12*mm, 672.*mm, 823.65*mm );
  fld->SetMesh( "CS0", 35, 18 );

  fld->AddCoil( "CS1", -595.25*mm, 795.25*mm, 672.*mm, 823.65*mm );
  fld->SetMesh( "CS1", 45, 18 );

  fld->AddCoil( "MS1", -2121.375*mm, -653.625*mm, 672.*mm, 756.25*mm );
  fld->SetMesh( "MS1", 45, 10 );

  fld->AddCoil( "MS2", -2910.5*mm, -2189.5*mm, 672.*mm, 789.95*mm );
  fld->SetMesh( "MS2", 50, 14 );
}


XCoilBase* XCOMETConstruction :: GetConductor()
{
  XCoilConductor* cdt = new XCoilConductor();
  cdt->SetDimension( 4.73*mm, 15.*mm );
  cdt->SetInsSize( 0.15*mm, 0.15*mm );

  return dynamic_cast<XCoilBase*>(cdt);
}


XCoilBase* XCOMETConstruction :: GetStrip(const double thick)
{
  XCoilStrip* strip = new XCoilStrip();
  strip->SetDimension( 4.73*mm+0.15*2*mm, thick );
  strip->SetInsSize( 0., 0.25*mm );

  return dynamic_cast<XCoilBase*>(strip);
}


XCoilBase* XCOMETConstruction :: GetShell(const double thick)
{
  XCoilShell * shell = new XCoilShell();
  shell->SetDimension( 4.73*mm+0.15*2*mm, thick );
  shell->SetInsSize( 0., 3.*mm );

  return dynamic_cast<XCoilBase*>(shell);
}


double XCOMETConstruction :: GetCoilResistance(XThermalSolver* solve)
{
  double res = 0.;
  int idz = 0;
  int idp = 0;
  int idr = 0;

  const int mshz = solve->GetProcess()->GetMesh(iZ);
  const int mshp = solve->GetProcess()->GetMesh(iPhi);
  const int mshr = solve->GetProcess()->GetMesh(iR);
  
  for (unsigned int i=0; i<solve->GetProcess()->GetMaterialEntries(); i++) {
    idz = solve->GetProcess()->GetDimensionEntry(i)->GetId(iZ);
    idp = solve->GetProcess()->GetDimensionEntry(i)->GetId(iPhi);
    idr = solve->GetProcess()->GetDimensionEntry(i)->GetId(iR);

    if (solve->GetProcess()->GetDimensionEntry(i)->GetGeometry()==kConductor &&
        idz>0 && idz<mshz+1 &&
        idp>0 && idp<mshp+1 &&
        idr>0 && idr<mshr+1 )
      res += solve->GetProcess()->GetMaterialEntry(i)->GetResistance();
  }

  return res;
}


void XCOMETConstruction :: UpdateQuench(XThermalSolver* solve)
{
  XMatCopper    cu;
  XMatAluminium al;
  XMatNbTi      sc;

  double T     = 0.;
  double RRR   = 0.;
  double B     = 0.;
  double R_Al  = 0.;
  double R_Cu  = 0.;
  double R_avg = 0.;
  double Tcs   = 0.;
  double Tc    = 0.;
  double Rcs   = 0.;

  int idz = 0;
  int idp = 0;
  int idr = 0;

  const int mshz = solve->GetProcess()->GetMesh(iZ);
  const int mshp = solve->GetProcess()->GetMesh(iPhi);
  const int mshr = solve->GetProcess()->GetMesh(iR);

  const double factor = solve->GetProcess()->GetCoilHandler()->GetApproachZ();

  const double l_Phi = solve->GetProcess()->GetCoilHandler()->GetCoilSize(iPhi) / solve->GetProcess()->GetMesh(iPhi);
  const double r_Cu  = solve->GetProcess()->GetCoilHandler()->GetMaterialRatio(iCopper);
  const double r_Al  = solve->GetProcess()->GetCoilHandler()->GetMaterialRatio(iAluminium);
  //const double A_cdt = solve->GetProcess()->GetCoilHandler()->GetCoilParts(kConductor)->GetArea();
  const double A_cdt = solve->GetProcess()->GetCoilHandler()->GetCoilType(2)->GetArea();

  const double A_Cu  = A_cdt * r_Cu;
  const double A_Al  = A_cdt * r_Al;

  //double Volume = solve->GetProcess()->GetCoilHandler()->GetCoilParts(kConductor)->GetTotalArea() * l_Phi;
  double Volume = solve->GetProcess()->GetCoilHandler()->GetCoilType(2)->GetTotalArea() * l_Phi;

  for (unsigned int i=0; i<solve->GetProcess()->GetMaterialEntries(); i++) {
    idz = solve->GetProcess()->GetDimensionEntry(i)->GetId(iZ);
    idp = solve->GetProcess()->GetDimensionEntry(i)->GetId(iPhi);
    idr = solve->GetProcess()->GetDimensionEntry(i)->GetId(iR);

    if ( solve->GetProcess()->GetDimensionEntry(i)->GetGeometry()==kConductor &&
         idz>0 && idz<mshz+1 &&
         idp>0 && idp<mshp+1 &&
         idr>0 && idr<mshr+1 ) {
      T   = solve->GetProcess()->GetMaterialEntry(i)->GetTemperature();
      RRR = solve->GetProcess()->GetMaterialEntry(i)->GetRRR();
      B   = solve->GetProcess()->GetMaterialEntry(i)->GetField();

      // update Temperature, RRR and Field
      cu.SetMaterialProperty(T, 50., B);
      al.SetMaterialProperty(T, RRR, B);
      sc.SetMaterialProperty(T, RRR, B);

      // calculate Tcs and Tc
      Tcs = sc.GetSharingT( fCurr );
      Tc  = sc.GetCriticalT();

      // calculate resistance
      R_Al  = al.GetResistivity() * l_Phi / A_Al;
      R_Cu  = cu.GetResistivity() * l_Phi / A_Cu;
      R_avg = pow(1./R_Al + 1./R_Cu, -1);
      
      // setup quench status
      if ( T<Tcs ) {
        solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kSuperconduct);
        solve->GetProcess()->GetMaterialEntry(i)->SetResistance(0.);
        solve->GetProcess()->GetMaterialEntry(i)->SetVoltage(0.);
        solve->GetProcess()->GetMaterialEntry(i)->SetHeat(0.);
      }
      else if ( T>=Tcs && T<Tc ) {
        Rcs = R_avg * (T-Tcs) / (Tc-Tcs);
        solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kTransition);
        solve->GetProcess()->GetMaterialEntry(i)->SetResistance( Rcs*factor );
        solve->GetProcess()->GetMaterialEntry(i)->SetVoltage( Rcs*fCurr*factor );
        solve->GetProcess()->GetMaterialEntry(i)->SetHeat( pow(fCurr,2)*Rcs/Volume );
      }
      else if ( T>=Tc ) {
        solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kNormal);
        solve->GetProcess()->GetMaterialEntry(i)->SetResistance( R_avg*factor );
        solve->GetProcess()->GetMaterialEntry(i)->SetVoltage( R_avg*fCurr*factor );
        solve->GetProcess()->GetMaterialEntry(i)->SetHeat( pow(fCurr,2)*R_avg/Volume );
      }

    }
    else {
      solve->GetProcess()->GetMaterialEntry(i)->SetResistance(0.);
      solve->GetProcess()->GetMaterialEntry(i)->SetVoltage(0.);
      solve->GetProcess()->GetMaterialEntry(i)->SetHeat(0.);
    }
  }
}


void XCOMETConstruction :: Begin()
{
  std::cout << "start the quench simulation for COMET. " << std::endl;
  std::cout << "....................................................." << std::endl;
  std::cout << " Running" << std::endl;
  std::cout << "....................................................." << std::endl;

  time_t now = time(0);
  tm* local = localtime(&now);

  std::cout << "time: " << asctime(local) << std::endl;
}


void XCOMETConstruction :: Run()
{
  const double Tcool = 4.5*K;

  double dt = fdt;
  fCS0->SetTimeInterval(dt);
  fCS1->SetTimeInterval(dt);
  fMS1->SetTimeInterval(dt);
  fMS2->SetTimeInterval(dt);
  
  /*
  fCS0->UseCylinderConnect();
  fCS1->UseCylinderConnect();
  fMS1->UseCylinderConnect();
  fMS2->UseCylinderConnect();
  */

  double time = fTime0;
  int cnt = 0;
  int ocnt = 0;
  double CoilRes = 0.;
  double qchtime = fTimef;
  bool   quenched = false;
  bool   preqch   = false;
  double timestep[4];

  const int numcdt = GetTotalConductor(fCS0) + GetTotalConductor(fCS1) +
                     GetTotalConductor(fMS1) + GetTotalConductor(fMS2);
  int qchcdt = 0;

  while (time<fTimef) {
    
    // 1. update material thermal conductivity, capacity
    fCS0->GetProcess()->SetMaterial();
    fCS1->GetProcess()->SetMaterial();
    fMS1->GetProcess()->SetMaterial();
    fMS2->GetProcess()->SetMaterial();

    // 2. update material resistivity, voltage
    UpdateQuench(fCS0); 
    UpdateQuench(fCS1);
    UpdateQuench(fMS1);
    UpdateQuench(fMS2);

    // 3. calculate coil resistance
    CoilRes = 0.;
    CoilRes += GetCoilResistance(fCS0);
    CoilRes += GetCoilResistance(fCS1);
    CoilRes += GetCoilResistance(fMS1);
    CoilRes += GetCoilResistance(fMS2);

    preqch = quenched;

    if ( fVth<CoilRes*fCurr )
      quenched = true;
    else
      quenched = false;

    if ( preqch==false && quenched==true )
      qchtime = time + fDetTime;

    // 4. calculate the current decay
    if ( quenched==true && time>qchtime ) {
      fCurr = CalCurrentDecay(fPreI, CoilRes, dt);
      fPreI = fCurr;

    // 5. calculate the field decay
      CalFieldDecay(fCS0);
      CalFieldDecay(fCS1);
      CalFieldDecay(fMS1);
      CalFieldDecay(fMS2);
    }

    // set heat generation before quench
    if ( quenched==false )
      fCS1->GetProcess()->GetMaterialEntry(fCS1->GetProcess()->Id(fHotZ,fHotPhi,fHotR))->SetHeat(fHotSpotHeat);

    // 6. solve the thermal equation
    fCS0->Solve(dt);
    fCS0->SetBoundary();

    fCS1->Solve(dt);
    fCS1->SetBoundary();

    fMS1->Solve(dt);
    fMS1->SetBoundary();
    
    fMS2->Solve(dt);
    fMS2->SetBoundary();

    ConnectShell(fCS0, fCS1, 5.*cm);
    ConnectShell(fCS1, fMS1, 10.*cm);
    ConnectShell(fMS1, fMS2, 10.*cm);

    ConnectMagnet(fCS0, fCS1, 50.*cm);
    ConnectMagnet(fCS1, fMS1, 50.*cm);
    ConnectMagnet(fMS1, fMS2, 50.*cm);

    // count quenched conductors
    qchcdt = GetQuenchConductor(fCS0) + GetQuenchConductor(fCS1) +
             GetQuenchConductor(fMS1) + GetQuenchConductor(fMS2);

    // print out
    if (dt>0.01) fDisplay=1;

    if (cnt%fDisplay==0) {
      std::cout << "time: " << time << " [sec], step: " << dt << " [sec], Rtot: "
                << CoilRes  << " [Ohm], Vtot: " << CoilRes*fCurr << " [V], I: "
                << fCurr << " [A], Ratio: " << static_cast<double>(qchcdt)/numcdt*100. << " %, diode: "
                << fDiode << " [V], ";
      fCS1->Print(fHotZ, fHotPhi, fHotR);
    }

    // fill data into file
    if (cnt%(fDisplay*8)==0) {
      XRootOutput output;
      if (ocnt==0)
        output.SetPath("./output");
      output.SetFilename(Form("./output/qchout%i.root",ocnt));
      output.SetSubDirectory("CS0");
      output.SetSubDirectory("CS1");
      output.SetSubDirectory("MS1");
      output.SetSubDirectory("MS2");
      output.SetHeader(cnt, time, fCurr, CoilRes, CoilRes*fCurr);
      output.Fill("CS0", fCS0->GetProcess());
      output.Fill("CS1", fCS1->GetProcess());
      output.Fill("MS1", fMS1->GetProcess());
      output.Fill("MS2", fMS2->GetProcess());
      output.Close();

      ocnt ++;
    }

    // find minimum time step
    timestep[0] = fCS0->FindTimeStep();
    timestep[1] = fCS1->FindTimeStep();
    timestep[2] = fMS1->FindTimeStep();
    timestep[3] = fMS2->FindTimeStep();

    dt = timestep[0];

    for (int i=1; i<4; i++) {
      if (timestep[i]<dt)
        dt = timestep[i];
    }
    
    time += dt;
    cnt ++;
  }

}


void XCOMETConstruction :: ConnectShell(XThermalSolver* mag1, XThermalSolver* mag2, const double l)
{
  const int mshz2 = mag2->GetProcess()->GetMesh(iZ);
  const int mshr2 = mag2->GetProcess()->GetMesh(iR);
  const int mshr1 = mag1->GetProcess()->GetMesh(iR);

  int id;
  double temp, k, q;
  // connect 1st shell
  for (int i=1; i<5; i++) {
    id = mag2->GetProcess()->Id(mshz2,i,mshr2);
    temp = mag2->GetProcess()->GetMaterialEntry(id)->GetTemperature();

    k = mag2->GetProcess()->GetMaterialEntry(id)->GetConductivity(iPhi);
    q = mag2->GetProcess()->GetMaterialEntry(id)->GetHeatFlux(iPhi);
    temp += l * q / k;

    id = mag1->GetProcess()->Id(0, i, mshr1);
    mag1->GetProcess()->GetMaterialEntry(id)->SetTemperature(temp);
  
    // connect 2nd coil
    id = mag1->GetProcess()->Id(1, i, mshr1);
    temp = mag1->GetProcess()->GetMaterialEntry(id)->GetTemperature();

    q = mag1->GetProcess()->GetMaterialEntry(id)->GetHeatFlux(iPhi);
    temp += l * q / k;

    id = mag2->GetProcess()->Id(mshz2+1, 0, mshr2);
    mag2->GetProcess()->GetMaterialEntry(id)->SetTemperature(temp);
  }
}

void XCOMETConstruction :: ConnectMagnet(XThermalSolver* mag1, XThermalSolver* mag2, const double l)
{
  // connect magnet1 with magnet2
  //const int mshz1 = mag1->GetProcess()->GetMesh(iZ);
  //const int mshr1 = mag1->GetProcess()->GetMesh(iR);

  const int mshz2 = mag2->GetProcess()->GetMesh(iZ);
  const int mshr2 = mag2->GetProcess()->GetMesh(iR);

  // connect 1st coil
  int id = mag2->GetProcess()->Id(mshz2,1,mshr2-1);
  double temp = mag2->GetProcess()->GetMaterialEntry(id)->GetTemperature();

  double k = mag2->GetProcess()->GetMaterialEntry(id)->GetConductivity(iPhi);
  //const double l = 50.*cm;
  double q = mag2->GetProcess()->GetMaterialEntry(id)->GetHeatFlux(iPhi);
  temp += l * q / k;

  id = mag1->GetProcess()->Id(1, 5, 2);
  mag1->GetProcess()->GetMaterialEntry(id)->SetTemperature(temp);
  
  // connect 2nd coil
  id = mag1->GetProcess()->Id(1, 4, 2);
  temp = mag1->GetProcess()->GetMaterialEntry(id)->GetTemperature();

  q = mag1->GetProcess()->GetMaterialEntry(id)->GetHeatFlux(iPhi);
  temp += l * q / k;

  id = mag2->GetProcess()->Id(mshz2, 0, mshr2-1);
  mag2->GetProcess()->GetMaterialEntry(id)->SetTemperature(temp);
}


void XCOMETConstruction :: End()
{
  std::cout << "....................................................." << std::endl;
  std::cout << " Finished" << std::endl;
  std::cout << "....................................................." << std::endl;

  time_t now = time(0);
  tm* local = localtime(&now);

  std::cout << "time: " << asctime(local) << std::endl;
}



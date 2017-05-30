#include <iostream>
#include <cmath>
#include <TString.h>
#include "XMatCopper.hpp"
#include "XMatAluminium.hpp"
#include "XMatNbTi.hpp"

#include "XCoilConductor.hpp"
#include "XCoilStrip.hpp"
#include "XCoilShell.hpp"

#include "XRootOutput.hpp"
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XQuenchOutput.hpp"
#include "XCOMETConstruction.h"


XCOMETConstruction :: XCOMETConstruction()
    : fFld(NULL),
      fCS(NULL),
      XQuenchTransient()
{
  if (!fFld)  fFld = new XFieldHandle();
  ConstructField(fFld);
}


XCOMETConstruction :: ~XCOMETConstruction()
{
  if (fFld)  delete fFld;
  if (fCS)   delete fCS;
}


void XCOMETConstruction ::ConstructAtlas()
{
  const double r = (1229.+45./2)*mm;
  const std::string name = "CS";
  
  XCoilHandle* coil = new XCoilHandle();
  coil->SetName(name);
  coil->SetCoilSize(0., 2.*M_PI*r, 0.);
  coil->SetMesh(800, 1, 3);
  coil->SetCoilLayers(1);
  coil->SetCoilTurns(1173);
  coil->SetMaterialRatio(14, 0.9, 1.);

  // set coil structure
  coil->AddLayer(1, kStrip, GetStrip(), kAdiabatic, 10.*cm);
  coil->AddLayer(2, kConductor, GetConductor(), kAdiabatic, 0.*mm);
  coil->AddLayer(3, kShell, GetShell(), kAdiabatic, 0.*mm);

  // set processor
  XProcessManager* pro = new XProcessManager();
  pro->SetCoilHandler(coil);
  pro->Initialize();
  pro->SetNbTiIc(22.92e+3);
  
  // get magnetic field map
  fFld->SetTarget(name);
  fFld->Run();
  pro->SetFieldHandler(fFld);
  
  // uniform RRR and magnetic field
  //pro->SetUniformField(2.0);
  pro->SetUniformRRR(kConductor, 570.);
  pro->SetUniformRRR(kStrip, 3000.);
  //pro->SetUniformRRR(kShell, 10.);

  // fill the materail info into the coil solver
  if (!fCS) fCS = new XThermalSolver();
  fCS->SetProcessHandle(pro);
}


void XCOMETConstruction :: ConstructField(XFieldHandle* fld)
{
  fld->SetCurrent(7600.*Amp);
  fld->AddCoil( "CS", 0.*mm, 5300.*mm, 1229.*mm, (1229.+33.)*mm );
  fld->SetMesh( "CS", 800, 1 );
}


XCoilBase* XCOMETConstruction :: GetConductor()
{
  XCoilConductor* cdt = new XCoilConductor();
  cdt->SetDimension( 4.25*mm, 30.*mm );
  cdt->SetInsSize( 0.1*mm, 0.1*mm );

  return dynamic_cast<XCoilBase*>(cdt);
}


XCoilBase* XCOMETConstruction :: GetStrip()
{
  XCoilStrip* strip = new XCoilStrip();
  strip->SetDimension( 4.25*mm+0.1*2*mm, 1.*mm );
  strip->SetInsSize( 0., 0.5*mm );

  return dynamic_cast<XCoilBase*>(strip);
}


XCoilBase* XCOMETConstruction :: GetShell()
{
  XCoilShell * shell = new XCoilShell();
  shell->SetDimension( 4.25*mm+0.15*2*mm, 12.*mm );
  shell->SetInsSize( 0., 0.5*mm );

  return dynamic_cast<XCoilBase*>(shell);
}


void XCOMETConstruction :: SetQuenchHeating(XThermalSolver* solve)
{
  solve->GetProcess()->GetMaterialEntry(solve->GetProcess()->Id(fHotZ,fHotPhi,fHotR))->SetHeat(6.4 * 15./0.0015);
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


void XCOMETConstruction :: UpdateQuench(XThermalSolver* solve, const double time)
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
      Tcs = sc.GetSharingT( fCurr, T );
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
        solve->GetProcess()->GetMaterialEntry(i)->SetResistance(Rcs*factor);
        solve->GetProcess()->GetMaterialEntry(i)->SetVoltage( factor*Rcs*fCurr );
        solve->GetProcess()->GetMaterialEntry(i)->SetHeat( pow(fCurr,2)*Rcs/Volume );
        if ( solve->GetProcess()->GetMaterialEntry(i)->GetQuenchTime()<0. )
          solve->GetProcess()->GetMaterialEntry(i)->SetQuenchTime(time);
      }
      else if ( T>=Tc ) {
        solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kNormal);
        solve->GetProcess()->GetMaterialEntry(i)->SetResistance(R_avg*factor);
        solve->GetProcess()->GetMaterialEntry(i)->SetVoltage( R_avg*fCurr*factor );
        solve->GetProcess()->GetMaterialEntry(i)->SetHeat( pow(fCurr,2)*R_avg/Volume );
      }

    }
    else {
      solve->GetProcess()->GetMaterialEntry(i)->SetStatus(kSuperconduct);
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

  XQuenchOutput* geo = new XQuenchOutput("geometry.dat", iOfstream);
  geo->WriteGeometry(fCS->GetProcess());
  geo->Close();
}


void XCOMETConstruction :: Run()
{
  const double Tcool = 4.5*K;

  double dt = fdt;
  fCS->SetTimeInterval(dt);
  fCS->SetAccelerateFactor(1.);

  double time = fTime0;
  int cnt = 0;
  int ocnt = 0;
  double CoilRes = 0.;
  double qchtime = fTimef;
  bool   quenched = false;
  bool   preqch   = false;

  while (time<fTimef) {
    
    // 1. update material thermal conductivity, capacity
    fCS->GetProcess()->SetMaterial();

    // 2. update material resistivity, voltage
    UpdateQuench(fCS, time); 

    // 3. calculate coil resistance
    CoilRes = 0.;
    CoilRes += GetCoilResistance(fCS);

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
      CalFieldDecay(fCS);
    }

    // set heat generation before quench
    //if ( quenched==false )
    if ( time<=1. )
      SetQuenchHeating(fCS);
    //if ( quenched==false )
    //  fCS->GetProcess()->GetMaterialEntry(fCS->GetProcess()->Id(fHotZ,fHotPhi,fHotR))->SetHeat(6400.*4);

    // 6. solve the thermal equation
    fCS->Solve(dt);
    
    fCS->SetBoundary();

    if (dt>0.01) fDisplay=1;

    //if (cnt%fDisplay==0 && (int)(time*10000)%10==0) {
    if ((int)(time*1000000)%5000==0) {
      std::cout << std::setprecision(4) << "time: " << time << " [sec], step: " << dt << " [sec], Rtot: "
                << CoilRes  << " [Ohm], Vtot: " << CoilRes*fCurr << " [V], I: "
                << fCurr << " [A]";
      fCS->Print(fHotZ,fHotPhi,fHotR);
    }

    //if (cnt%(fDisplay*20)==0) {
    if ((int)(time*1000000)%50000==0) {
      XRootOutput output( Form("./output/qchout%i.root",ocnt) );
      output.SetSubDirectory("CS");
      output.SetHeader(cnt, time, fCurr, CoilRes, CoilRes*fCurr);
      output.Fill("CS", fCS->GetProcess());
      output.Close();
      ocnt++;
    }

    dt = fCS->FindTimeStep();
    
    time += dt;
    cnt ++;
  }

}


void XCOMETConstruction :: End()
{
  std::cout << "....................................................." << std::endl;
  std::cout << " Finished" << std::endl;
  std::cout << "....................................................." << std::endl;

  XQuenchOutput* out = new XQuenchOutput("output.dat", iOfstream);
  out->Write(fCS->GetProcess());
  out->Close();
}



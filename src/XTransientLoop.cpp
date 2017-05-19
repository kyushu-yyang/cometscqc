#include <iostream>
#include "XQuenchLogger.hpp"
#include "XQuenchExcept.hpp"
#include "XThermalSolver.hpp"
#include "XTransientLoop.hpp"

using Quench::XQuenchLogger;

XTransientLoop :: XTransientLoop()
    : fSolver(NULL),
      fTime0(0.*sec),
      fTimef(60.*sec),
      fdt(1.*msec),
      fDisplay(500)
{
  if (!fSolver) fSolver = new XThermalSolver();
}

XTransientLoop :: ~XTransientLoop()
{
  if (fSolver) delete fSolver;
}

void XTransientLoop :: SetTime(const double t0, const double tf, const double dt)
{
  fTime0 = t0;
  fTimef = tf;
  fdt    = dt;

  QuenchError( XQuenchLogger::INFO, "setup time -> t0:" << fTime0 << "sec, tf:"
                                    << fTimef << "sec, dt:" << fdt << "sec" );
}

void XTransientLoop :: Begin()
{
  std::cout << "running the thermal transient loop..." << std::endl;
}

void XTransientLoop :: SetProcess(Quench::XProcessManager* hand)
{
  fSolver->SetProcessHandle(hand);
}

void XTransientLoop :: Run()
{
  double dt = fdt;
  fSolver->SetTimeInterval(dt);

  double time = fTime0;
  const double Tcool = 4.5*K;
  int cnt = 0;

  //double print_t = 0.1;

  while (time<fTimef) {

    //fSolver->GetProcess()->SetMaterial();

    fSolver->Solve(dt);

    fSolver->SetBoundary();

    for (unsigned int k=1; k<fSolver->GetProcess()->GetMesh(iR)+1; k++)
      fSolver->SetCoolingPath( k, Tcool, fSolver->GetProcess()->GetCoilHandler()->GetCoolingConfigure(k) );

    // setup cooling point on shell
    for (unsigned int i=1; i<fSolver->GetProcess()->GetMesh(iZ)+1; i++) {
      fSolver->SetLastCoolingPoint(i, Tcool);
    }

    if (cnt%fDisplay==0) {
      std::cout << "time: " << time << " [sec], step: " << dt << " [sec] ";
      fSolver->Print(fSolver->GetProcess()->GetMesh(iZ)*4./7., 3, 2);
    }

    dt = fSolver->FindTimeStep();

    time += dt;
    cnt ++;
  }

}

void XTransientLoop :: End()
{
  for (int i=0; i<fSolver->GetProcess()->GetMaterialEntries(); i++) {
    std::cout << fSolver->GetProcess()->GetDimensionEntry(i)->GetId(iZ) << "  "
              << fSolver->GetProcess()->GetDimensionEntry(i)->GetId(iPhi) << "  "
              << fSolver->GetProcess()->GetDimensionEntry(i)->GetId(iR) << "  "
              << fSolver->GetProcess()->GetMaterialEntry(i)->GetTemperature() << std::endl;
  }
}

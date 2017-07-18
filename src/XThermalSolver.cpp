#include <iostream>
#include <cmath>
#include "XQuenchLogger.hpp"
#include "XQuenchExcept.hpp"
#include "XThermalSolver.hpp"

using Quench::XQuenchLogger;

XThermalSolver :: XThermalSolver()
    : fProcess(NULL),
      fMshZ(0),
      fMshP(0),
      fMshR(0),
      fAcce(1.),
      fdt(0.1*sec),
      fCylinder(false)
{}

XThermalSolver :: ~XThermalSolver()
{
  if (fProcess)  delete fProcess;
}

void XThermalSolver :: SetProcessHandle(Quench::XProcessManager* hand)
{
  if (!hand) {
    QuenchError( XQuenchLogger::ERROR, "process handler is null." );
    XQuenchExcept except("process handler is null.");
    throw except;
  }

  fProcess = hand;

  fMshZ = fProcess->GetMesh(iZ);
  fMshP = fProcess->GetMesh(iPhi);
  fMshR = fProcess->GetMesh(iR);
}


void XThermalSolver :: SetAccelerateFactor(double acc)
{
  fAcce = acc;
  QuenchInfo("set the factor on time step: " << fAcce );
}


void XThermalSolver :: SetTimeInterval(const double dt)
{
  fdt = dt;

  int id = 0;

  for (int k=0; k<fMshR+2; k++) {
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) {
        id = fProcess->Id(i,j,k);
        fProcess->GetMaterialEntry(id)->SetTimeStep(fdt);
      }
    }
  }
}

void XThermalSolver :: Initial()
{
  int id = 0;
  double T = 4.5;

  for (int k=0; k<fMshR+2; k++) {
    for (int j=0; j<fMshP+2; j++) {
      for (int i=0; i<fMshZ+2; i++) {
        id = fProcess->Id(i,j,k);
        T = fProcess->GetMaterialEntry(id)->GetTemperature();
        fProcess->GetMaterialEntry(id)->SetPreTemp(T);
      }
    }
  }
}

double XThermalSolver :: FindTimeStep() const
{
  double minstep = 10000.;
  double step; 
  int    id;

  for (int k=1; k<fMshR+1; k++) {
    for (int j=1; j<fMshP+1; j++) {
      for (int i=1; i<fMshZ+1; i++) {
        id = fProcess->Id(i,j,k);
        step = fProcess->GetMaterialEntry(id)->GetTimeStep();
        minstep = step < minstep ? step : minstep;
      }
    }
  }

  return minstep;
}

void XThermalSolver :: Solve(const double dt)
{
  fdt = dt;
  
  Initial();

  // thermal calculation
  for (int k=1; k<fMshR+1; k++) {
    for (int j=1; j<fMshP+1; j++) {
      for (int i=1; i<fMshZ+1; i++) {
        InTheLoop(i, j, k);
      }
    }
  }

}

void XThermalSolver :: InTheLoop(const int i, const int j, const int k)
{
  // calculate id
  const int zpr   = fProcess->Id(  i,   j,   k);
  const int preZ  = fProcess->Id(i-1,   j,   k);
  const int postZ = fProcess->Id(i+1,   j,   k);
  const int preP  = fProcess->Id(  i, j-1,   k);
  const int postP = fProcess->Id(  i, j+1,   k);
  const int preR  = fProcess->Id(  i,   j, k-1);
  const int postR = fProcess->Id(  i,   j, k+1);

  // get pre-temperature
  const double T      = fProcess->GetMaterialEntry(  zpr)->GetPreTemp();
  const double preTz  = fProcess->GetMaterialEntry( preZ)->GetPreTemp();
  const double postTz = fProcess->GetMaterialEntry(postZ)->GetPreTemp();
  const double preTp  = fProcess->GetMaterialEntry( preP)->GetPreTemp();
  const double postTp = fProcess->GetMaterialEntry(postP)->GetPreTemp();
  const double preTr  = fProcess->GetMaterialEntry( preR)->GetPreTemp();
  const double postTr = fProcess->GetMaterialEntry(postR)->GetPreTemp();

  // get distance
  const double dpreZ  = fProcess->GetDimensionEntry(zpr)->GetPosition(iZ) -
                        fProcess->GetDimensionEntry(zpr)->GetPrePosition(iZ);
  const double dpostZ = fProcess->GetDimensionEntry(zpr)->GetPostPosition(iZ) - 
                        fProcess->GetDimensionEntry(zpr)->GetPosition(iZ);

  const double dpreP  = fProcess->GetDimensionEntry(zpr)->GetPosition(iPhi) -
                        fProcess->GetDimensionEntry(zpr)->GetPrePosition(iPhi);
  const double dpostP = fProcess->GetDimensionEntry(zpr)->GetPostPosition(iPhi) - 
                        fProcess->GetDimensionEntry(zpr)->GetPosition(iPhi);

  const double dpreR  = fProcess->GetDimensionEntry(zpr)->GetPosition(iR) -
                        fProcess->GetDimensionEntry(zpr)->GetPrePosition(iR);
  const double dpostR = fProcess->GetDimensionEntry(zpr)->GetPostPosition(iR) - 
                        fProcess->GetDimensionEntry(zpr)->GetPosition(iR);

  // get cell size
  const double lz = fProcess->GetDimensionEntry(zpr)->GetCellSize(iZ);
  const double lp = fProcess->GetDimensionEntry(zpr)->GetCellSize(iPhi);
  const double lr = fProcess->GetDimensionEntry(zpr)->GetCellSize(iR);

  // get density
  const double rho = fProcess->GetMaterialEntry(zpr)->GetDensity();

  // get heat capacity
  const double C = fProcess->GetMaterialEntry(zpr)->GetCapacity();

  // get thermal conductivity
  const double kz = fProcess->GetMaterialEntry(zpr)->GetConductivity(iZ);
  const double kp = fProcess->GetMaterialEntry(zpr)->GetConductivity(iPhi);
  const double kr = fProcess->GetMaterialEntry(zpr)->GetConductivity(iR);

  // calculate pre and post thermal conductivity
  double kpreZ = fProcess->GetMaterialEntry(preZ)->GetConductivity(iZ);
  if (i==1) kpreZ = kz;
  double kpostZ = fProcess->GetMaterialEntry(postZ)->GetConductivity(iZ);
  if (i==fMshZ) kpostZ = kz;

  double kpreP = fProcess->GetMaterialEntry(preZ)->GetConductivity(iPhi);
  if (j==1) kpreP = fProcess->GetMaterialEntry(fProcess->Id(i,fMshP+1,k))->GetConductivity(iPhi);
  double kpostP = fProcess->GetMaterialEntry(postZ)->GetConductivity(iPhi);
  if (j==fMshP) kpostP = fProcess->GetMaterialEntry(fProcess->Id(i,1,k))->GetConductivity(iPhi);
  
  double kpreR = fProcess->GetMaterialEntry(preR)->GetConductivity(iR);
  if (k==1) kpreR = kr;
  double kpostR = fProcess->GetMaterialEntry(postR)->GetConductivity(iR);
  if (k==fMshR) kpostR = kr;

  if (fProcess->GetDimensionEntry(zpr)->GetGeometry()==kShell)
    kpostR = fProcess->GetMaterialEntry(zpr)->GetConductivity(iZ);

  const double prekZ  = 2*pow(1./kpreZ + 1./kz, -1);
  const double postkZ = 2*pow(1./kpostZ + 1./kz, -1);

  const double prekR  = dpreR / (lr/2./kr + (dpreR-lr/2.)/kpreR);
  const double postkR = dpostR / (lr/2./kr + (dpostR-lr/2.)/kpostR);

  // get diffusion velocity
  double az = kz / rho / C;
  double ap = kp / rho / C;
  double ar = kr / rho / C;

  // calculate time step
  double step = 1. / (az/dpreZ/lz + ap/dpreP/lp + ar/dpreR/lr) / 2.;
  step *= fAcce;

  fProcess->GetMaterialEntry(zpr)->SetTimeStep(step);

  // get heat generation 
  //const double gen = fProcess->GetMaterialEntry(zpr)->GetDeposit() * 4000.;
  const double gen = fProcess->GetMaterialEntry(zpr)->GetHeat();

  // calculate heat flux grad [W/m3]
  //double qz = kz * ( (postTz-T)/dpostZ - (T-preTz)/dpreZ ) / lz;
  double qz = postkZ * (postTz-T)/dpostZ/lz - prekZ * (T-preTz)/dpreZ/lz;
  //double qp = kp * ( (postTp-T)/dpostP - (T-preTp)/dpreP ) / lp;
  double qp = kp * ( (postTp-T)/dpostP - (T-preTp)/dpreP ) / lp;
  //double qr = kr * ( (postTr-T)/dpostR - (T-preTr)/dpreR ) / lr;
  double qr = postkR * (postTr-T)/dpostR/lr - prekR * (T-preTr)/dpreR/lr;

  // calculate heat flux (center differential)
  const double fluz = -kz * (postTz - preTz) / (dpostZ + dpreZ);
  const double flup = -kp * (postTp - preTp) / (dpostP + dpreP);
  const double flur = -kr * (postTr - preTr) / (dpostR + dpreR);

  // set heat flux
  fProcess->GetMaterialEntry(zpr)->SetHeatFlux( fluz, flup, flur );

  double Q = (qz + qp + qr + gen) / rho / C;
  double Temp = T + fdt*Q;

  //std::cout << Temp << " " << Q << " " << qz << " " << qp << " " << qr << " "
            //<< gen << " " << kz << " " << kp << " " << kr << " "
            //<< dpostZ << " " << dpreZ << " " << lz << " "
            //<< i << " " << j << " " << k << " "
            //<< std::endl;

  fProcess->GetMaterialEntry(zpr)->SetTemperature(Temp);
}

void XThermalSolver :: SetBoundary()
{
  ///////////////////////////////////////
  // R DIRECTION
  ///////////////////////////////////////
  int id_bdy  = 0;
  int id_edge = 0;
  double T = 0.;

  // adiabatic boundary
  for (int j=0; j<fMshP+2; j++) {
    for (int i=0; i<fMshZ+2; i++) {
      id_bdy  = fProcess->Id(i,j,0);
      id_edge = fProcess->Id(i,j,1);
      T = fProcess->GetMaterialEntry(id_edge)->GetTemperature();
      fProcess->GetMaterialEntry(id_bdy)->SetTemperature(T);

      id_bdy  = fProcess->Id(i,j,fMshR+1);
      id_edge = fProcess->Id(i,j,fMshR);
      T = fProcess->GetMaterialEntry(id_edge)->GetTemperature();
      fProcess->GetMaterialEntry(id_bdy)->SetTemperature(T);
    }
  }

  ///////////////////////////////////////
  // Z DIRECTION
  ///////////////////////////////////////
  // abiabatic boundary
  for (int k=0; k<fMshR+2; k++) {
    for (int j=0; j<fMshP+2; j++) {
      id_bdy  = fProcess->Id(0,j,k);
      id_edge = fProcess->Id(1,j,k);
      T = fProcess->GetMaterialEntry(id_edge)->GetTemperature();
      fProcess->GetMaterialEntry(id_bdy)->SetTemperature(T);

      id_bdy  = fProcess->Id(fMshZ+1,j,k);
      id_edge = fProcess->Id(fMshZ,j,k);
      T = fProcess->GetMaterialEntry(id_edge)->GetTemperature();
      fProcess->GetMaterialEntry(id_bdy)->SetTemperature(T);
    }
  }

  ///////////////////////////////////////
  // PHI DIRECTION
  ///////////////////////////////////////
  if (fCylinder==true)
    SetCylinderPhi();
  else
    SetConductorPhi();
 
}


void XThermalSolver :: SetCylinderPhi() {
  int id_bdy = 0;
  int id_edge = 0;
  double T = 0.;

  for (int k=0; k<fMshR+2; k++) {
    for (int i=0; i<fMshZ+2; i++) {
      id_bdy  = fProcess->Id(i, 0, k);
      id_edge = fProcess->Id(i, fMshP, k);
      T = fProcess->GetMaterialEntry(id_edge)->GetTemperature();
      fProcess->GetMaterialEntry(id_bdy)->SetTemperature(T);

      id_bdy  = fProcess->Id(i, fMshP+1, k);
      id_edge = fProcess->Id(i, 1, k);
      T = fProcess->GetMaterialEntry(id_edge)->GetTemperature();
      fProcess->GetMaterialEntry(id_bdy)->SetTemperature(T);
    }
  }
}


void XThermalSolver :: SetConductorPhi() {
  int id_bdy = 0;
  int id_edge = 0;
  double T = 0.;

  for (int k=0; k<fMshR+2; k++) {
    for (int i=0; i<fMshZ+2; i++) {
      id_bdy  = fProcess->Id(i,0,k);
      id_edge = fProcess->Id(i,fMshP,k);
      Connect(id_edge, id_bdy);

      id_bdy  = fProcess->Id(i,fMshP+1,k);
      id_edge = fProcess->Id(i,1,k);
      Connect(id_edge, id_bdy);

      // leave the first and last turn
      if ( i>2 && i<fMshZ && fProcess->GetDimensionEntry(id_bdy)->GetGeometry()==kConductor ) {
        id_bdy  = fProcess->Id(i,0,k);
        id_edge = fProcess->Id(i-1,fMshP,k);    // privous turn
        T = fProcess->GetMaterialEntry(id_edge)->GetTemperature();
        fProcess->GetMaterialEntry(id_bdy)->SetTemperature(T);

        id_bdy  = fProcess->Id(i-1,fMshP+1,k);
        id_edge = fProcess->Id(i,1,k);
        Connect(id_edge, id_bdy);
      }

      if ( i==1 && fProcess->GetDimensionEntry(id_bdy)->GetGeometry()==kConductor ) {
        id_bdy  = fProcess->Id(i, 0, k);
        id_edge = fProcess->Id(i, 1, k);
        T = fProcess->GetMaterialEntry(id_edge)->GetTemperature();
        fProcess->GetMaterialEntry(id_bdy)->SetTemperature(T);
      }

      if ( i==fMshZ && fProcess->GetDimensionEntry(id_bdy)->GetGeometry()==kConductor ) {
        id_bdy  = fProcess->Id(i, fMshP+1, k);
        id_edge = fProcess->Id(i, fMshP, k);
        T = fProcess->GetMaterialEntry(id_edge)->GetTemperature();
        fProcess->GetMaterialEntry(id_bdy)->SetTemperature(T);
      }
    }
  }
}


void XThermalSolver :: Connect(const int from_z, const int from_p, const int from_r,
                               const int to_z,   const int to_p,   const int to_r)
{
  int from_id = fProcess->Id(from_z, from_p, from_r);
  int to_id   = fProcess->Id(to_z, to_p, to_r);
  
  double T = fProcess->GetMaterialEntry(from_id)->GetTemperature();
  fProcess->GetMaterialEntry(to_id)->SetTemperature(T);
}


void XThermalSolver :: Connect(const int from_id, const int to_id)
{
  double T = fProcess->GetMaterialEntry(from_id)->GetTemperature();
  fProcess->GetMaterialEntry(to_id)->SetTemperature(T);
}


void XThermalSolver :: SetCoolingPath(const int r, const double T, const Cooling opt)
{
  int bdy = fProcess->Id(1,1,r);

  /*
  if (fProcess->GetDimensionEntry(bdy)->GetGeometry()!=kStrip) 
    QuenchError( XQuenchLogger::WARNING, "layer: " << r << " is not strip.");
  */

  // set cooling path
  for (int j=0; j<fMshP+2; j++) {
    if (opt==kLeft) {
      bdy  = fProcess->Id(0,j,r);
      fProcess->GetMaterialEntry(bdy)->SetTemperature(T);
      //QuenchInfo( "layer: " << r << " is cooling from left." );
    }
    else if (opt==kRight) {
      bdy = fProcess->Id(fMshZ+1,j,r);
      fProcess->GetMaterialEntry(bdy)->SetTemperature(T);
      //QuenchInfo( "layer: " << r << " is cooling from right." );
    } 
    else if (opt==kSide) {
      bdy  = fProcess->Id(0,j,r);
      fProcess->GetMaterialEntry(bdy)->SetTemperature(T);

      bdy = fProcess->Id(fMshZ+1,j,r);
      fProcess->GetMaterialEntry(bdy)->SetTemperature(T);
      //QuenchInfo( "layer: " << r << " is cooling from two side." );
    }
    //else
      //QuenchInfo( "layer: " << r << " is adiabatic condition." );
  }
}


void XThermalSolver :: SetLastCoolingPoint(const int z, const double T)
{
  int id = 0;

  for (int j=1; j<fMshP+1; j++) {
    //id = fProcess->Id(z, j, fMshR);
    //if (fProcess->GetDimensionEntry(id)->GetGeometry()!=kShell)
      //throw
    id = fProcess->Id(z, j, fMshR+1);
    fProcess->GetMaterialEntry(id)->SetTemperature(T);
  }
}


void XThermalSolver :: SetFirstCoolingPoint(const int z, const double T)
{
  int id = 0;

  for (int j=1; j<fMshP+1; j++) {
    id = fProcess->Id(z, j, 0);
    fProcess->GetMaterialEntry(id)->SetTemperature(T);
  }

}


void XThermalSolver :: Print(const int z, const int phi, const int r)
{
  if (fPrint.size()==0)
    fPrint.push_back( fProcess->Id(z, phi, r) );

  double Temp = 0.;
  double preTemp = 0.;
  double step = 0.;

  for (std::vector<int>::size_type i=0; i<fPrint.size(); i++) {
    Temp    = fProcess->GetMaterialEntry(fPrint.at(i))->GetTemperature();
    preTemp = fProcess->GetMaterialEntry(fPrint.at(i))->GetPreTemp();
    step    = fProcess->GetMaterialEntry(fPrint.at(i))->GetTimeStep();

    std::cout << " T: "     << Temp << " [K], " 
              << " dT/dt: " << (Temp-preTemp)/step << " [K/sec] ";

    //QuenchInfo(" T: "     << Temp << " [K], " << " dT/dt: " << (Temp-preTemp)/step << " [K/sec] ");

    if (i==fPrint.size()-1) {
      std::cout << "\n";
    }
  }
}


void XThermalSolver :: AddOutput(const int z, const int phi, const int r)
{
  fPrint.push_back( fProcess->Id(z, phi, r) );
  
  QuenchInfo("set output temperature at mesh (" << z << ", " << phi << ", " << r << ")");
}


void XThermalSolver :: UseCylinderConnect()
{
  fCylinder = true;

  QuenchInfo("using the cylinderial connection at phi direction.");
}


void XThermalSolver :: UseConductorConnect()
{
  fCylinder = false;

  QuenchInfo("using the conductor connection at phi direction");
}

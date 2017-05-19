#include <iostream>
#include <cmath>
#include "XQuenchLogger.hpp"
#include "XMagneticField.hpp"

using Quench::XFieldContainer;
using Quench::XBiotSavart;

XBiotSavart :: XBiotSavart()
    : fMapMesh(NULL), fMapRange(NULL),
      fCurrent(2700.), fCdtA(4.75*mm * 15.05*mm),
      fAMsh(NULL), fSolenoid(NULL)
{
  if (!fAMsh)  fAMsh = new int[3];

  fAMsh[0] = 10;  /// R
  fAMsh[1] = 80;  /// Phi
  fAMsh[2] = 80;  /// Z
}

XBiotSavart :: ~XBiotSavart()
{
  if (fMapMesh)  delete [] fMapMesh;
  if (fMapRange) delete [] fMapRange;
  if (fAMsh)     delete [] fAMsh;
  if (fSolenoid) delete [] fSolenoid;
}

void XBiotSavart :: SetSolenoidMesh(const int mz, const int mp, const int mr)
{
  if (!fAMsh) fAMsh = new int[3];

  fAMsh[0] = mr;
  fAMsh[1] = mp;
  fAMsh[2] = mz;
}

void XBiotSavart :: SetSolenoid(const double z0, const double z1,
                                const double r0, const double r1)
{
  if (!fSolenoid)  fSolenoid = new double[4];

  fSolenoid[0] = r0;
  fSolenoid[1] = r1;
  fSolenoid[2] = z0;
  fSolenoid[3] = z1;
}

void XBiotSavart :: SetMapMesh(const int mz, const int mr)
{
  if (!fMapMesh)  fMapMesh = new int[2];

  fMapMesh[0] = mr;
  fMapMesh[1] = mz;
}

void XBiotSavart :: SetMapRange(const double z0, const double z1,
                                const double r0, const double r1)
{
  if (!fMapRange)  fMapRange = new double[4];

  fMapRange[0] = r0;
  fMapRange[1] = r1;
  fMapRange[2] = z0;
  fMapRange[3] = z1;
}

void XBiotSavart :: calfield()
{
  const double dr = (fMapRange[1] - fMapRange[0]) / fMapMesh[0];
  const double dz = (fMapRange[3] - fMapRange[2]) / fMapMesh[1];

  double R, Z;
  double B[2];
  double A[ fMapMesh[0]+1 ][ fMapMesh[1]+1 ];
  XFieldContainer* fld = NULL;

  if (fHC.size() != 0) {
    QuenchError( XQuenchLogger::DEBUG, "clear the field container." );
    fHC.clear();
  }
  
  for (int i=0; i<fMapMesh[0]+1; i++) {
    R = fMapRange[0] + i*dr + dr/2.;

    for (int j=0; j<fMapMesh[1]+1; j++) {
      Z = fMapRange[2] + j*dz + dz/2.;
      // calculate A_phi
      A[i][j] = calpotential(Z, R);

      if (i>0 && j>0) {
        B[0] = -(A[i][j] - A[i][j-1]) / dz;                    // Br
        B[1] = A[i][j] / R + (A[i][j] - A[i-1][j]) / dr;       // Bz

        fld = new XFieldContainer();
        fld->SetId(j, i);
        fld->SetPosition(Z, R);
        fld->SetField(B[1], B[0]);
        fld->SetPotential(A[i][j]);
        fHC.push_back(fld);
      }
      else {
        //B[0] = -(A[i][j] - 0.) / dz;                   // Br
        //B[1] = A[i][j] / R + (A[i][j] - 0.) / dr;      // Bz
        B[0] = 0.;
        B[1] = 0.;
      }
    }
  }

  QuenchError( XQuenchLogger::DEBUG, "finished the field calculation." );
}

double XBiotSavart :: calpotential(const double z, const double r) const
{
  const double dr = (fSolenoid[1] - fSolenoid[0]) / fAMsh[0];
  const double dp = (2 * M_PI - 0.) / fAMsh[1];
  const double dz = (fSolenoid[3] - fSolenoid[2]) / fAMsh[2];

  double R, Z, PHI, A;
  double rPQ = 0.;
  double Aphi = 0.;

  for (int i=0; i<fAMsh[0]; i++) {    // R
    R = fSolenoid[0] + i*dr + dr/2.;

    for (int j=0; j<fAMsh[1]; j++) {  // PHI
      PHI = 0. + j*dp + dp/2.;

      for (int k=0; k<fAMsh[2]; k++) {  // Z
        Z = fSolenoid[2] + k*dz + dz/2.;  
        
        rPQ = sqrt( pow(Z-z,2) + pow(R,2) + pow(r,2) - 2*r*R*cos(PHI) ); 

        if ( rPQ==0. )
          A = 0.;
        else 
          A = (mu0 * fCurrent / 4 / M_PI) * (R * cos(PHI) * dr * dz * dp / rPQ);

        Aphi += A;
      }
    }
  }

  return Aphi;
}


XFieldContainer* XBiotSavart :: GetFieldEntry(const int iz, const int jr)
{
  if (fHC.size()==0) {
    QuenchError( XQuenchLogger::DEBUG, "started the field calculation." );
    calfield();
  }

  return fHC.at( jr*fMapMesh[1] + iz );
}

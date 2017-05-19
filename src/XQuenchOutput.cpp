#include <iostream>
#include <TFile.h>

#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XQuenchOutput.hpp"

using Quench::XQuenchLogger;
using Quench::XQuenchOutput;


XQuenchOutput :: XQuenchOutput(const std::string &filename, const FileOutput opt)
    : fFilename(filename),
      fRootfile(NULL)
{
  init(opt);
}


XQuenchOutput :: ~XQuenchOutput()
{
  if (fRootfile) fRootfile->Close();
  if (fNormfile.is_open()) fNormfile.close();
}


void XQuenchOutput :: init(const FileOutput opt)
{
  switch (opt) {
    case iROOT: 
      fRootfile = new TFile(fFilename.c_str(), "recreate");
      QuenchError( XQuenchLogger::CONFIG, "write data to the ROOT file." );
      break;
    case iOfstream:
      fNormfile.open( fFilename, std::ios::out );
      QuenchError( XQuenchLogger::CONFIG, "write data to the stream file.");
      break;
    default:
      QuenchError( XQuenchLogger::ERROR, "option " << opt << " is not available." );
      XQuenchExcept except("output option is not available.");
      throw except;
  }
}


void XQuenchOutput :: WriteGeometry(XProcessManager* man)
{
  if (!man) {
    QuenchError( XQuenchLogger::ERROR, "processing manager is null.");
  }

  const int mshz = man->GetMesh(iZ);
  const int mshp = man->GetMesh(iPhi);
  const int mshr = man->GetMesh(iR);

  int idz = 0;
  int idp = 0;
  int idr = 0;

  double z   = 0.;
  double phi = 0.;
  double r   = 0.;

  double lz   = 0.;
  double lphi = 0.;
  double lr   = 0.;

  for (unsigned int i=0; i<man->GetEntries(); i++) {

    idz = man->GetDimensionEntry(i)->GetId(iZ);
    idp = man->GetDimensionEntry(i)->GetId(iPhi);
    idr = man->GetDimensionEntry(i)->GetId(iR);

    if ( idz>0 && idz<mshz+1 &&
         idp>0 && idp<mshp+1 &&
         idr>0 && idr<mshr+1 ) {
    
      lz   = man->GetDimensionEntry(i)->GetCellSize(iZ);
      lphi = man->GetDimensionEntry(i)->GetCellSize(iPhi);
      lr   = man->GetDimensionEntry(i)->GetCellSize(iR);

      z   = man->GetDimensionEntry(i)->GetPosition(  iZ) - lz/2.;
      phi = man->GetDimensionEntry(i)->GetPosition(iPhi) - lphi/2.;
      r   = man->GetDimensionEntry(i)->GetPosition(  iR) - lr/2.;

      if (idz==1 && man->GetDimensionEntry(i)->GetGeometry()==kStrip)
        z = man->GetDimensionEntry(i)->GetPrePosition(iZ);

      if (idz==mshz && man->GetDimensionEntry(i)->GetGeometry()==kStrip)
        z = man->GetDimensionEntry(i)->GetPrePosition(iZ)
           +man->GetDimensionEntry(i-1)->GetCellSize(iZ)/2.;

      fNormfile << i << " "
                << idz << " " << idp << " " << idr << " "
                << z << " " << phi << " " << r << " "
                << lz << " " << lphi << " " << lr << " "
                << man->GetDimensionEntry(i)->GetGeometry() << " "
                <<"\n";
    }
  }
}


void XQuenchOutput :: SetHeader(XProcessManager* man)
{
  // write the header info for data recording

}


void XQuenchOutput :: Write(XProcessManager* man)
{
  if (!man) {
    QuenchError( XQuenchLogger::ERROR, "processing manager is null.");
  }

  const int mshz = man->GetMesh(iZ);
  const int mshp = man->GetMesh(iPhi);
  const int mshr = man->GetMesh(iR);

  int idz = 0;
  int idp = 0;
  int idr = 0;

  for (unsigned int i=0; i<man->GetEntries(); i++) {

    idz = man->GetDimensionEntry(i)->GetId(iZ);
    idp = man->GetDimensionEntry(i)->GetId(iPhi);
    idr = man->GetDimensionEntry(i)->GetId(iR);

    if ( idz>0 && idz<mshz+1 &&
         idp>0 && idp<mshp+1 &&
         idr>0 && idr<mshr+1 ) {
      fNormfile << i << " "
                << man->GetMaterialEntry(i)->GetTemperature()      << " "
                << man->GetMaterialEntry(i)->GetRRR()              << " "
                << man->GetMaterialEntry(i)->GetField()            << " "
                << man->GetMaterialEntry(i)->GetCapacity()         << " "
                << man->GetMaterialEntry(i)->GetConductivity(iZ)   << " "
                << man->GetMaterialEntry(i)->GetConductivity(iPhi) << " "
                << man->GetMaterialEntry(i)->GetConductivity(iR)   << " "
                << man->GetMaterialEntry(i)->GetHeat()             << " "
                <<"\n";

    }
  }
}

void XQuenchOutput :: Close()
{
  if (fRootfile) fRootfile->Close();
  if (fNormfile.is_open()) fNormfile.close();
}

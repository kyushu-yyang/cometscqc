#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XRootOutput.hpp"

using Quench::XQuenchLogger;
using Quench::XProcessManager;


XRootOutput :: XRootOutput()
    : fFile(NULL),
      fInfo(NULL),
      fName(""),
      fMshZ(0),
      fMshP(0),
      fMshR(0)
{}


XRootOutput :: XRootOutput(const char* filename)
    : fFile(NULL),
      fInfo(NULL),
      fName(""),
      fMshZ(0),
      fMshP(0),
      fMshR(0)
{
  SetFilename(filename);
}


XRootOutput :: ~XRootOutput()
{
  //if (fInfo) delete fInfo;
}


void XRootOutput :: SetFilename(const char* filename)
{
  fName = std::string(filename);

  if (!fFile)  fFile = new TFile(filename, "RECREATE");
  if (!fInfo)  fInfo = new TTree("info", "data info");
}


void XRootOutput :: SetHeader(const int i, const double t, const double I, const double R, const double V)
{
  TTree* tree = new TTree("head", "data info");

  int id;
  double time, curr, res, volt;

  tree->Branch("id", &id, "id/I");
  tree->Branch("time", &time, "time/D");
  tree->Branch("I", &curr, "I/D");
  tree->Branch("R", &res, "R/D");
  tree->Branch("V", &volt, "V/D");

  id = i;
  time = t;
  curr = I;
  res = R;
  volt = V;

  tree->Fill();
  tree->Write();

  fInfo->Branch("name", &fName);
  fInfo->Branch("mshz", &fMshZ, "mshz/I");
  fInfo->Branch("mshp", &fMshP, "mshp/I");
  fInfo->Branch("mshr", &fMshR, "mshr/I");
}


void XRootOutput :: SetSubDirectory(const char* dirname)
{
  TDirectory* dir = fFile->mkdir(dirname);
  fDir.insert( std::map<const char*, TDirectory*>::value_type(dirname,dir) );
}


void XRootOutput :: Fill(const char* name, XProcessManager* man)
{
  fName = std::string(name);
  fMshZ = man->GetMesh(iZ);
  fMshP = man->GetMesh(iPhi);
  fMshR = man->GetMesh(iR);
  fInfo->Fill();
  fInfo->Write();

  // access the given directory
  fDir.at(name)->cd();

  // create tree file
  TTree* tree = new TTree("tree", name);

  double Temp, preTemp, RRR, B, C, R, Q;
  int    node;
  double k  [3];
  int    id [3];
  double pos[3];
  double flx[3];
  int    status;

  tree->Branch("node", &node, "node/I");
  tree->Branch("id", id, "id[3]/I");
  tree->Branch("pos", pos, "pos[3]/D");
  tree->Branch("C", &C, "C/D");
  tree->Branch("k", k, "k[3]/D");
  tree->Branch("T", &Temp, "T/D");
  tree->Branch("preT", &preTemp, "preT/D");
  tree->Branch("RRR", &RRR, "RRR/D");
  tree->Branch("B", &B, "B/D");
  tree->Branch("R", &R, "R/D");
  tree->Branch("Q", &Q, "Q/D");
  tree->Branch("q", flx, "q[3]/D");
  tree->Branch("status", status, "status/I");

  for (unsigned int i=0; i<man->GetEntries(); i++) {

    node = man->GetDimensionEntry(i)->GetNodeId();

    id[0] = man->GetDimensionEntry(i)->GetId(iZ);
    id[1] = man->GetDimensionEntry(i)->GetId(iPhi);
    id[2] = man->GetDimensionEntry(i)->GetId(iR);

    pos[0] = man->GetDimensionEntry(i)->GetPosition(iZ);
    pos[1] = man->GetDimensionEntry(i)->GetPosition(iPhi);
    pos[2] = man->GetDimensionEntry(i)->GetPosition(iR);

    k[0] = man->GetMaterialEntry(i)->GetConductivity(iZ);
    k[1] = man->GetMaterialEntry(i)->GetConductivity(iPhi);
    k[2] = man->GetMaterialEntry(i)->GetConductivity(iR);

    flx[0] = man->GetMaterialEntry(i)->GetHeatFlux(iZ);
    flx[1] = man->GetMaterialEntry(i)->GetHeatFlux(iPhi);
    flx[2] = man->GetMaterialEntry(i)->GetHeatFlux(iR);

    C = man->GetMaterialEntry(i)->GetCapacity();
    Temp = man->GetMaterialEntry(i)->GetTemperature();
    preTemp = man->GetMaterialEntry(i)->GetPreTemp();
    RRR = man->GetMaterialEntry(i)->GetRRR();
    B = man->GetMaterialEntry(i)->GetField();
    R = man->GetMaterialEntry(i)->GetResistance();
    Q = man->GetMaterialEntry(i)->GetHeat();
    status = man->GetMaterialEntry(i)->GetStatus();

    if ( id[0]>0 && id[0]<fMshZ+1 &&
         id[1]>0 && id[1]<fMshP+1 &&
         id[2]>0 && id[2]<fMshR+1 ) {
      tree->Fill();
    }
  }

  tree->Write();
}


void XRootOutput :: Close()
{
  if (fFile) {
    fFile->Write();
    fFile->Close();
  }
}


void XRootOutput :: SetPath(const char* path)
{
  struct stat info;
  mode_t mode = 0755;

  if ( stat(path, &info)!=0 ) {
    // cannot access this directory
    mkdir(path, mode);
    QuenchInfo("create new directory: " << path );
  }
  else if ( info.st_mode & S_IFDIR )
    // it is a directory
    QuenchInfo("find directory: " << path );
  else {
    mkdir(path, mode);
    QuenchInfo("create new directory: " << path );
  }
}

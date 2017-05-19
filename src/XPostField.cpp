#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2F.h>
#include "XQuenchLogger.hpp"
#include "XPostField.hpp"

using Quench::XQuenchLogger;
using Quench::XFieldHandle;

XPostField :: XPostField()
    : fMagnet(""), fInfo(NULL)
{}

XPostField :: XPostField(XFieldHandle* fld)
    : fMagnet(""), fInfo(NULL)
{
  SetFieldHandler(fld);
}

XPostField :: ~XPostField()
{
  if (fInfo) delete fInfo;
}

void XPostField :: SetFieldHandler(XFieldHandle* fld)
{
  fMagnet  = fld->GetTarget();
  fInfo    = fld->GetInfoEntry(fMagnet);
  fCollect = fld->GetFieldCollection();
}

void XPostField :: Plot()
{
  TCanvas* c0 = new TCanvas(fMagnet.c_str(), fMagnet.c_str(), 750, 420);
  c0->SetTicks(1,1);
  c0->SetRightMargin(0.13);

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(99);

  TH2F* hist = NULL;
  plot2d(hist);
}

void XPostField :: plot2d(TH2F* hist)
{
  const int xbin = fInfo->GetMesh()[1];     // z direction
  const int ybin = fInfo->GetMesh()[0];     // r direction

  const double xmin = fInfo->GetDimension()[2] / mm;
  const double xmax = fInfo->GetDimension()[3] / mm;
  const double ymin = fInfo->GetDimension()[0] / mm;
  const double ymax = fInfo->GetDimension()[1] / mm;

  hist = new TH2F(fMagnet.c_str(), fMagnet.c_str(), xbin, xmin, xmax,
                                                    ybin, ymin, ymax);

  double Bz, Br;
  double B;
  double Bmax = 0.;
  double z, r;

  for (std::vector<int>::size_type i=0; i<fCollect.size(); i++) {
    z  = fCollect.at(i)->GetPosition()[0] / mm;
    r  = fCollect.at(i)->GetPosition()[1] / mm;
    Bz = fCollect.at(i)->GetField()[0];
    Br = fCollect.at(i)->GetField()[1];
    B  = sqrt( pow(Bz,2) + pow(Br,2) );

    hist->Fill(z, r, B);
    Bmax = B > Bmax ? B : Bmax;
  }

  std::cout << " Maximum field: " << Bmax << " [Tesla]" << std::endl;
  QuenchError( XQuenchLogger::INFO, " Maximum field: " << Bmax << " [Tesla]" );

  hist->SetTitle( Form("%s; Z [mm]; R [mm]; Magnetic Field [Tesla]",fMagnet.c_str()) );
  hist->Draw("colz");
}

#include <iostream>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TStyle.h>

#include "XQuenchInfo.hpp"
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XPostInfoPlot.hpp"

using Quench::XDimensionInfo;
using Quench::XMaterialInfo;
using Quench::XQuenchLogger;

XPostInfoPlot :: ~XPostInfoPlot()
{
  delete fProcess;
}


void XPostInfoPlot :: SetProcessManager(Quench::XProcessManager* pro)
{
  if (!pro) {
    QuenchError( XQuenchLogger::WARNING, "XProcessManager is null." );
  }

  fProcess = pro;
}


void XPostInfoPlot :: PlotNode(const bool save)
{
  TCanvas* c0 = new TCanvas("node", "node", 640, 450);
  c0->SetTicks(1,1);

  int idOfPhi = 1;
  GetNodeGraph(idOfPhi)->Draw("ap");

  if (save==true)
    c0->Print("NodePosInfo.pdf");
}


void XPostInfoPlot :: PlotField(const bool save)
{
  TCanvas* c0 = new TCanvas("field", "field", 640, 450);
  c0->SetTicks(1,1);
  c0->SetRightMargin(0.13);

  gStyle->SetNumberContours(75);
  //gStyle->SetPalette(51);

  //GetFieldGraph()->Draw("cont4z");
  GetFieldGraph()->Draw("colz");

  if (save==true)
    c0->Print("MagField.pdf");
}


TGraph* XPostInfoPlot :: GetNodeGraph(const int phi)
{
  if (!fProcess) {
    QuenchError( XQuenchLogger::ERROR, "XProcessManager is not set." );
    XQuenchExcept except("please set the XProcessManager first." );
    throw except;
  }

  std::vector<XDimensionInfo*> dim = fProcess->GetDimensionContainer();

  TGraph* gr = new TGraph();
  int cnt = 0;

  for (std::vector<int>::size_type i=0; i<dim.size(); i++) {
    if ( dim.at(i)->GetId(iPhi)==phi ) {
      gr->SetPoint( cnt, dim.at(i)->GetPosition(iZ), dim.at(i)->GetPosition(iR) );
      cnt ++;
    }
  }

  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kRed);
  gr->SetTitle( Form("#phi = %i; Z [m]; R [m]", phi) );

  return gr;
}


TGraph2D* XPostInfoPlot :: GetFieldGraph()
{
  if (!fProcess) {
    QuenchError( XQuenchLogger::ERROR, "XProcessManager is not set." );
    XQuenchExcept except("please set the XProcessManager first." );
    throw except;
  }

  std::vector< XMaterialInfo*> mat = fProcess->GetMaterialContainer();
  std::vector<XDimensionInfo*> dim = fProcess->GetDimensionContainer();

  TGraph2D* gr = new TGraph2D();
  int cnt = 0;

  for (std::vector<int>::size_type i=0; i<mat.size(); i++) {
    if ( dim.at(i)->GetId(iPhi)==1 ) {
      gr->SetPoint( cnt, dim.at(i)->GetId(iZ), dim.at(i)->GetId(iR), mat.at(i)->GetField() );
      cnt++;
    }
  }

  gr->SetTitle("; Z [m]; R[m]; Magnetic Field [Tesla]");
  return gr;
}

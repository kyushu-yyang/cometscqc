#include <iostream>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TLegend.h>
#include "XQuenchLogger.hpp"
#include "XPlotMaterial.hpp"

using Quench::XQuenchLogger;

XPlotMaterial :: XPlotMaterial()
    : fOpt("R"), fRange(NULL),
      fMat(NULL), fMg(new TMultiGraph()),
      fPlots(0),
      fLg(NULL)
{
  fLg = new TLegend(0.6, 0.6, 0.89, 0.89);
  fLg->SetBorderSize(0);
  fLg->SetFillStyle(0);
  fLg->SetTextSize(0.043);

  gStyle->SetTitleSize(0.043, "xy");
  gStyle->SetLabelSize(0.043, "xy");
  gStyle->SetTitleOffset(1.3, "y");
}


XPlotMaterial :: ~XPlotMaterial()
{
  if (fRange)  delete [] fRange;
  if (fMat)    delete fMat;
  if (fMg)     delete fMg;
  if (fLg)     delete fLg;
}


void XPlotMaterial :: SetRange(const double min, const double max)
{
  if (!fRange)  fRange = new double[2];

  fRange[0] = min;
  fRange[1] = max;
}


void XPlotMaterial :: Add(const std::string& opt, const double var)
{
  TGraph* gr = new TGraph();

  if (opt=="RRR") {
    fMat->SetRRR(var);
    fLg->AddEntry( gr, Form("RRR: %.1f", var), "l" );
  }
  else if (opt=="B" || opt=="field") {
    fMat->SetField(var);
    fLg->AddEntry( gr, Form("B = %.1f Tesla", var), "l" );
  }
  else { 
    QuenchError( XQuenchLogger::WARNING, "this option " << opt << " did not exist." );
  }

  const double T0 = fRange[0];
  const double Tf = fRange[1];
  const int    nT = 200;
  const double dT = (Tf - T0) / nT;

  double T = T0;
  double par = 0.;

  for (int i=0; i<nT; i++) {
    fMat->SetTemperature(T);
    if (fOpt=="R" || fOpt=="resist")
      par = fMat->GetResistivity();
    else if (fOpt=="k" || fOpt=="conduct")
      par = fMat->GetConductivity();
    gr->SetPoint(i, T, par);
    T += dT;
  }

  gr->SetLineWidth(2.);
  gr->SetLineColor(fPlots+kSpring+1);
  fMg->Add(gr, "l");
  fPlots++;
}


void XPlotMaterial :: Plot()
{
  if (fOpt=="R" || fOpt=="resist")
    fMg->SetTitle("; Temperature [K]; Resistivity [#Omega #upoint m]");
  else if (fOpt=="k" || fOpt=="conduct")
    fMg->SetTitle("; Temperature [K]; Thermal Conductivity [W/m/K]");

  TCanvas* c0 = new TCanvas("mat", "mat", 650, 480);
  c0->SetTicks(1,1);
  c0->SetLogy();

  gStyle->SetTitleSize(0.043, "xy");
  gStyle->SetLabelSize(0.043, "xy");
  gStyle->SetTitleOffset(1.3, "y");

  fMg->Draw("a");
  fLg->Draw();
}


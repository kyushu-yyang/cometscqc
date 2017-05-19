#include <iostream>
#include "XCoilBase.hpp"
#include "XQuenchLogger.hpp"
#include "XQuenchExcept.hpp"
#include "XInitialTemperature.hpp"

using Quench::XQuenchLogger;

XInitialTemperature :: XInitialTemperature()
    : fIntialTemp(4.5),
      fFactor(1.)
{}

XInitialTemperature :: ~XInitialTemperature()
{}

void XInitialTemperature :: SetFactor(const double factor)
{
  fFactor = factor;
  QuenchInfo("set the factor: " << fFactor << " on initial temperature." );
}

void XInitialTemperature :: Load(const char* filename)
{
  if (fTemp.size()!=0)
    fTemp.clear();

  std::ifstream file(filename);

  if (!file) {
    QuenchFatal( "cannot load this file:" << filename << "." );
    XQuenchExcept except("failed to load this file.");
    throw except;
  }

  int id;
  double T;
  double RRR;
  double B;
  double C;
  double k[3];
  double q;

  XInitTempContainer* temp = NULL;

  while (true) {
    file >> id >> T >> RRR >> B >> C >> k[0] >> k[1] >> k[2] >> q;
    if (!file)  break;
    temp = new XInitTempContainer;
    temp->SetNodeId( id );
    temp->SetTemperature( T );
    fTemp.push_back( temp );
  }
}



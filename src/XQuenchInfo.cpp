#include <iostream>
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XQuenchInfo.hpp"

using Quench::XMaterialInfo;
using Quench::XDimensionInfo;
using Quench::XQuenchLogger;

XMaterialInfo :: XMaterialInfo()
    : fField(0.), 
      fPreTemp(4.5),
      fTemp(4.5), 
      fCapacity(0.),
      fRho(0.), 
      fRRR(100.),
      fStatus(kSuperconduct), 
      fGen(0.),
      fDose(0.),
      fStep(0.1),
      fVolt(0.),
      fStack(0.),
      fQcht(-1.),
      fTcs(-1.)
{}


void XMaterialInfo :: SetField(const double &fld)
{
  if (fld<0) {
    QuenchError( XQuenchLogger::WARNING, "field is " << fld );
    XQuenchExcept except("field is negative!");
    throw except;
  }

  fField = fld;
}

void XMaterialInfo :: SetTemperature(const double temp)
{
  fTemp = temp;
}

void XMaterialInfo :: SetCapacity(const double C)
{
  // check input capacity, if the capacity is negative, then throw it to exception
  if (C<0) {
    QuenchError( XQuenchLogger::ERROR, "heat capacity is {" << C );
    XQuenchExcept except("heat capacity is negative!");
    throw except;
  }

  fCapacity = C;
}


void XMaterialInfo :: SetHeatFlux(const double Qx, const double Qy, const double Qz)
{
  fHeat.at(0) = Qx;
  fHeat.at(1) = Qy;
  fHeat.at(2) = Qz;
}


void XMaterialInfo :: SetConductivity(const double kx, const double ky, const double kz)
{
  if (kx<0 || ky<0 || kz<0) {
    QuenchError(XQuenchLogger::ERROR, "thermal conductivity is {" << kx
                << ", " << ky << ", " << kz << "}");
    XQuenchExcept except("thermal conductivity is negative!");
    throw except;
  }

  fk.at(0) = kx;
  fk.at(1) = ky;
  fk.at(2) = kz;
}



/*****************************************************************************************************/
void XDimensionInfo :: SetId(const int i, const int j, const int k)
{
  if (i<0 || j<0 || k<0) {
    QuenchError(XQuenchLogger::ERROR, "cell id is {" << i << ", " << j << ", " << k << "}" );
    XQuenchExcept except("cell id is negative!");
    throw except;
  }

  fId.at(0) = i;
  fId.at(1) = j;
  fId.at(2) = k;
}


void XDimensionInfo :: SetPrePosition(const double x, const double y, const double z)
{
  fPrePos.at(0) = x;
  fPrePos.at(1) = y;
  fPrePos.at(2) = z;
}


void XDimensionInfo :: SetPosition(const double x, const double y, const double z)
{
  fPos.at(0) = x;
  fPos.at(1) = y;
  fPos.at(2) = z;
}


void XDimensionInfo :: SetPostPosition(const double x, const double y, const double z)
{
  fPostPos.at(0) = x;
  fPostPos.at(1) = y;
  fPostPos.at(2) = z;
}


void XDimensionInfo :: SetNodeId(const int node)
{
  if (node<0) {
    QuenchError(XQuenchLogger::ERROR, "node id is " << node);
    XQuenchExcept except("node id is negative!");
    throw except;
  }

  fNode = node;
}


void XDimensionInfo :: SetGeometry(const Geometry geo)
{
  fGeo = geo;
}


void XDimensionInfo :: SetCellSize(const double lx, const double ly, const double lz)
{
  fCell.at(0) = lx;
  fCell.at(1) = ly;
  fCell.at(2) = lz;
}

void XDimensionInfo :: SetDistance(const double dx, const double dy, const double dz)
{
  fDistance.at(0) = dx;
  fDistance.at(1) = dy;
  fDistance.at(2) = dz;
}



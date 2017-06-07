#include <iostream>
#include "XQuenchLogger.hpp"
#include "XQuenchExcept.hpp"
#include "XCoilHandle.hpp"

using Quench::XQuenchLogger;
using Quench::XCoilHandle;

XCoilHandle :: XCoilHandle()
    : fName(""),
      fCoil(NULL),
      fMsh(NULL),
      fLayer(0),
      fTurn(0),
      fEdge(20.*cm),
      fCdt(NULL),
      fStrip(NULL),
      fShell(NULL)
{
  // initialize the material ratio
  fRatio.insert( std::map<const Material, double>::value_type(iAluminium, 7.) );
  fRatio.insert( std::map<const Material, double>::value_type(   iCopper, 1.) );
  fRatio.insert( std::map<const Material, double>::value_type(     iNbTi, 1.) );
}

XCoilHandle :: ~XCoilHandle()
{
  if (fCoil) delete [] fCoil;
  if (fMsh ) delete [] fMsh;
  if (fCdt)  delete fCdt;
  if (fStrip) delete fStrip;
  if (fShell) delete fShell;
}


void XCoilHandle :: SetMaterialRatio(const double Al, const double Cu, const double SC)
{
  double tot = Al + Cu + SC;
  fRatio[iAluminium] = Al / tot;
  fRatio[   iCopper] = Cu / tot;
  fRatio[     iNbTi] = SC / tot;

  QuenchError( XQuenchLogger::INFO, "Magnet: " << fName << " : material ratio :: Al:Cu:SC = " 
                                                        << fRatio[iAluminium] 
                                                        << ":" << fRatio[iCopper]
                                                        << ":" << fRatio[iNbTi] );
}


void XCoilHandle :: SetAlPercent(const double perc)
{
  fRatio[iAluminium] = perc;
  QuenchError( XQuenchLogger::INFO, fName << " : material ratio :: Al = " << fRatio[iAluminium] );
}


void XCoilHandle :: SetCuPercent(const double perc)
{
  fRatio[iCopper] = perc;
  QuenchError( XQuenchLogger::INFO, fName << " : material ratio :: Cu = " << fRatio[iCopper] );
}


void XCoilHandle :: SetScPercent(const double perc)
{
  fRatio[iNbTi] = perc;
  QuenchError( XQuenchLogger::INFO, fName << " : material ratio :: SC = " << fRatio[iNbTi] );
}


double XCoilHandle :: GetMaterialRatio(const Material mat) const
{
  if (fRatio.find(mat)==fRatio.end()) {
    QuenchError( XQuenchLogger::ERROR, "not found material: " << mat );
    XQuenchExcept except("XCoilHandle: not found this material.");
    throw except;
  }
  else {
    return fRatio.find(mat)->second;
  }
}


void XCoilHandle :: SetCoilSize(const double lz, const double lp, const double lr)
{
  if (!fCoil) fCoil = new double[3];
  fCoil[  iZ] = lz;
  fCoil[iPhi] = lp;
  fCoil[  iR] = lr;

  QuenchError( XQuenchLogger::INFO, "coil " << fName << " size :: lz:" << fCoil[iZ] 
                                            << ", lp:" << fCoil[iPhi]
                                            << ", lr:" << fCoil[iR] );
}


double XCoilHandle :: GetCoilSize(const Coil dim) const
{
  if (dim>=3) {
    QuenchError( XQuenchLogger::ERROR, "the dimension" << dim << " is not correct." );
    XQuenchExcept except("XCoilHandle: the coil dimension is not correct.");
    throw except;
  }
  else {
    return fCoil[dim];
  }
}


void XCoilHandle :: SetMesh(const int mz, const int mp, const int mr)
{
  if (!fMsh) fMsh = new int[3];
  fMsh[  iZ] = mz;
  fMsh[iPhi] = mp;
  fMsh[  iR] = mr;

  QuenchError( XQuenchLogger::INFO, "coil mesh :: z:" << fMsh[iZ] << ", phi:" << fMsh[iPhi]
                                            << ", r:" << fMsh[iR] );
}


int XCoilHandle :: GetMesh(const Coil dim) const
{
  if (dim>=3) {
    QuenchError( XQuenchLogger::ERROR, "mesh dimension " << dim << " is wrong." );
    XQuenchExcept except("XCoilHandle: the mesh dimension is not right.");
    throw except;
  }
  else {
    return fMsh[dim];
  }
}


void XCoilHandle :: SetStripEdge(const double length)
{
  fEdge = length;

  QuenchError( XQuenchLogger::INFO, "the length of edge of aluminium: " << fEdge/cm << " cm" );
}


double XCoilHandle :: GetApproachZ() const 
{
  const double app = static_cast<double>(fTurn) / static_cast<double>(fMsh[iZ]);
  return app;
}


void XCoilHandle :: SetCoilParameters(const int turn, const int layer)
{
  SetCoilTurns(turn);
  SetCoilLayers(layer);
}


void XCoilHandle :: SetCoilLayers(const int layer)
{
  fLayer = layer;
  QuenchError( XQuenchLogger::INFO, "coil layer: " << fLayer );
}


void XCoilHandle :: SetCoilTurns(const int turn)
{
  fTurn = turn;
  QuenchError( XQuenchLogger::INFO, "coil turns: " << fTurn );
}


const Quench::XCoilBase* XCoilHandle :: GetCoilParts(const Geometry geo)
{
  switch (geo) {
    case kConductor: return fCdt;
    case kStrip: return fStrip;
    case kShell: return fShell;
    default: return NULL;
  }
}


void XCoilHandle :: AddLayer(const int layer, const Geometry geo)
{
  if ( !is_exist(layer) || !IsOutRange(layer) ) {
    fLayerGeo.insert( std::map<const int, const Geometry>::value_type(layer, geo) );
    QuenchError( XQuenchLogger::INFO, "coil layer = " << layer << ", geometry = " << GetGeometryName(geo) );
  }
  else {
    QuenchError( XQuenchLogger::ERROR, "input layer: " << layer << " already existed." );
  }
}


void XCoilHandle :: AddLayer(const int layer, const Geometry geo, const XCoilBase* coil, const Cooling cool, const double edgesize)
{
  if (!is_exist(layer) || !IsOutRange(layer)) {
    fLayerGeo.insert( std::map<const int, const Geometry>::value_type(layer, geo) );
    fCoilType.insert( std::map<const int, const XCoilBase*>::value_type(layer, coil) );
    fCooling .insert( std::map<const int, const Cooling>::value_type(layer, cool) );
    fEdgesize.insert( std::map<const int, const double>::value_type(layer, edgesize) );

    //QuenchInfo( "layer: " << layer << " => geometry: " << GetGeometryName(geo) );
  }
  else {
    QuenchInfo("this layer: " << layer << " existed already.");
  }
}


const Quench::XCoilBase* XCoilHandle :: GetCoilType(const int layer)
{
  if (is_exist(layer)) {
    return fCoilType.at(layer);
  }
  else {
    XQuenchExcept except("this layer does not exist.");
    throw except;
  }
}


const double XCoilHandle :: GetStripEdgeSize(const int layer)
{
  if (is_exist(layer)) {
    return fEdgesize.at(layer);
  }
  else {
    XQuenchExcept except("this layer does not exist.");
    throw except;
  }
}


const Cooling XCoilHandle :: GetCoolingConfigure(const int layer)
{
  if (is_exist(layer)) {
    return fCooling.at(layer);
  }
  else {
    XQuenchExcept except("this layer does not exist.");
    throw except;
  }
}


const Quench::XCoilBase* XCoilHandle :: GetCoilLayout(const int layer)
{
  if (is_exist(layer)) {
    switch ( fLayerGeo.find(layer)->second ) {
      case kConductor:  return fCdt;
      case kStrip:      return fStrip;
      case kShell:      return fShell;
      default: throw;
    }
  }
  else {
    QuenchError( XQuenchLogger::ERROR, "this layer: " << layer << " does not exist." );
    XQuenchExcept except("this layer does not exist.");
    throw except;
  }
}


std::vector<int> XCoilHandle :: GetLayerId(const Geometry geo)
{
  std::vector<int> id;

  for (std::map<const int, const Geometry>::const_iterator it=fLayerGeo.begin(); it!=fLayerGeo.end(); it++) {
    if ( it->second == geo ) {
      id.push_back( it->first );
    }
  }

  return id;
}


bool XCoilHandle :: IsOutRange(const int layer) const
{
   bool outrange = false;
   if (layer>fMsh[iR])  outrange = true;

   return outrange;
}


bool XCoilHandle :: is_out_range() const
{
  bool outrange = false;

  for (std::map<const int, const Geometry>::const_iterator it=fLayerGeo.begin(); it!=fLayerGeo.end(); it++) {
    if ( IsOutRange(it->first) )  {
      outrange = true;
      break;
    }
  }

  return outrange;
}


bool XCoilHandle :: is_exist(const int layer) const
{
  bool exist = false;

  for (std::map<const int, const Geometry>::const_iterator it=fLayerGeo.begin(); it!=fLayerGeo.end(); it++) {
    if (it->first==layer) {
      exist = true;
      break;
    }
  }

  return exist;
}


int XCoilHandle :: min_layer() const
{
  int min = 99999999;

  for (std::map<const int, const Geometry>::const_iterator it=fLayerGeo.begin(); it!=fLayerGeo.end(); it++) {
    if (it->first<min)
      min = it->first;
  }

  return min;
}


int XCoilHandle :: max_layer() const
{
  int max = -99999999;

  for (std::map<const int, const Geometry>::const_iterator it=fLayerGeo.begin(); it!=fLayerGeo.end(); it++) {
    if (it->first>max)
      max = it->first;
  }

  return max;
}


std::string XCoilHandle :: GetGeometryName(const Geometry geo)
{
  std::string name = "";
  switch (geo) {
    case kConductor: name = "Conductor"; break;
    case     kStrip: name = "Strip";     break;
    case     kShell: name = "Shell";     break;
    case       kG10: name = "G10";       break;
    case     kA5083: name = "A5083";     break;
    default: break;
  }

  return name;
}


std::string XCoilHandle :: GetMaterialName(const Material mat)
{
  std::string name = "";
  switch (mat) {
    case iAluminium: name = "Aluminium"; break;
    case    iCopper: name = "Copper";    break;
    case      iNbTi: name = "NbTi";      break;
    default: break;
  }

  return name;
}

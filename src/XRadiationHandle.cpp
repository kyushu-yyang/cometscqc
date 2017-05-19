#include <iostream>
#include <fstream>
#include "XQuenchLogger.hpp"
#include "XQuenchExcept.hpp"
#include "XRadiationHandle.hpp"

using Quench::XQuenchLogger;

void XRadiationContainer :: SetId(const int z, const int p, const int r)
{
  fId.at(iZ)   = z;
  fId.at(iPhi) = p;
  fId.at(iR)   = r;
}


void XRadiationContainer :: SetPosition(const double z, const double p, const double r)
{
  fPos.at(iZ)   = z;
  fPos.at(iPhi) = p;
  fPos.at(iR)   = r;
}


XRadiationHandle :: XRadiationHandle()
    : fName(""), 
      fIrrad(1.*year),
      fMshz(0),
      fMshp(0),
      fMshr(0)
{}

XRadiationHandle :: XRadiationHandle(const std::string& filename)
    : fName(""),
      fIrrad(1.*year),
      fMshz(0),
      fMshp(0),
      fMshr(0)
{
  Load(filename);
}

bool XRadiationHandle :: IsOverRange(const int id) const
{
  bool over_range = false;
  const int size = fMshz * fMshp * fMshr;
  if (id>=size)  over_range = true;

  return over_range;
}

const int XRadiationHandle :: Id(const int i, const int j, const int k)
{
  //const int id = i*fMshr*fMshp + j*fMshr + k;
  const int id = k*fMshp*fMshz + i*fMshp + j;
  if (IsOverRange(id)) {
    QuenchError( XQuenchLogger::ERROR, "local id is over range, z:" << i
                                       << ", phi: " << j << ", r:" << k
                                       << ", id: " << id );
    XQuenchExcept except("local id is over range.");
    throw except;
  }

  return id;
}

void XRadiationHandle :: SetIrrTime(const double time)
{
  fIrrad = time;

  if (fRC.size()>3) {
    QuenchError( XQuenchLogger::INFO, "set irradiation time: " << fIrrad << " sec <-> "
                                      << fIrrad/day << " day.");
    const double factor = fIrrad / (1*year);
    for (std::vector<int>::size_type i=0; i<fRC.size(); i++) 
      fRC.at(i)->SetNeutron( fRC.at(i)->GetNeutron() * factor );
  }
}

void XRadiationHandle :: Load(const std::string& filename)
{
  std::ifstream file(filename);

  if (!file) {
    QuenchError( XQuenchLogger::ERROR, "file: " << filename << " did not exist." );
    XQuenchExcept except("file did not exist.");
    throw except;
  }

  // factor to convert n/year to n/?
  const double factor = fIrrad / (1*year);

  int    z, phi, r;
  double dbuff[3];
    
  XRadiationContainer* rad = NULL;
  QuenchError( XQuenchLogger::INFO, "loading file: " << filename );

  while (true) {
    file >> z >> phi >> r >> dbuff[0] >> dbuff[1] >> dbuff[2];
    if (!file) break;
    rad = new XRadiationContainer;
    rad->SetId( z, phi, r );
    rad->SetNeutron( dbuff[0] * factor );
    rad->SetDose( dbuff[1] );
    fRC.push_back( rad );
  }

  fMshz = GetMesh(iZ);
  fMshp = GetMesh(iPhi);
  fMshr = GetMesh(iR);
}


const int XRadiationHandle :: GetMesh(const Coil dim)
{
  const int msh = findmax(dim) + 1;
  return msh;
}


double XRadiationHandle :: GetRRR(const Geometry geo, const double neu) const 
{
  const double irr = 2.7e-22;
  const double rho0_Cdt = 0.0675;
  const double rho0_Al  = 0.0135;

  double RRR = 1.;

  switch (geo) {
    case kConductor:
      RRR = 27. / (rho0_Cdt + neu*irr);
      break;
    case kStrip:
      RRR = 27. / (rho0_Al + neu*irr);
      break;
    default:
      break;
  }

  return RRR;
}


int XRadiationHandle :: findmax(const Coil dim) const
{
  int max = -999;
  int id;

  for (std::vector<int>::size_type i=0; i<fRC.size(); i++) {
    id = fRC.at(i)->GetId(dim);
    max = id > max ? id : max;
  }

  return max;
}

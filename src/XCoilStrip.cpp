#include <iostream>
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XCoilStrip.hpp"

using Quench::XCoilStrip;
using Quench::XQuenchLogger;


XCoilStrip :: XCoilStrip()
    : fSize(new double[2]),
      fIns(new double[2]),
      fType(kStrip)
{}


XCoilStrip :: ~XCoilStrip()
{
  if (fSize)  delete [] fSize;
  if (fIns)   delete [] fIns;
}


void XCoilStrip :: SetDimension(const double lz, const double lr)
{
  if (!fSize)  fSize = new double[2];

  fSize[0] = lz;
  fSize[1] = lr;
  QuenchError( XQuenchLogger::INFO, "strip::lz: " << lz << ", lr: " << lr );
}


double XCoilStrip :: GetDimension(const Coil dim) const
{
  switch (dim) {
    case iZ: return fSize[0];
    case iR: return fSize[1];
    default:
      QuenchError( XQuenchLogger::WARNING, "no this dimension: " << dim );
      return 0.;
  }
}


double XCoilStrip :: GetTotalLength(const Coil dim) const
{
  switch (dim) {
    case iZ: return fSize[0];
    case iR: return fSize[1] + 2*fIns[1];
    default:
      QuenchError( XQuenchLogger::WARNING, "no this dimension: " << dim );
      return 0.;
  }
}


double XCoilStrip :: GetInsSize(const Coil dim) const
{
  switch (dim) {
    case iZ: return fIns[0];
    case iR: return fIns[1];
    default:
      QuenchError( XQuenchLogger::WARNING, "no this dimension: " << dim );
      return 0.;
  }
}


void XCoilStrip :: SetInsSize(const double lz, const double lr)
{
  if (!fIns)  fIns = new double[2];

  fIns[0] = lz;
  fIns[1] = lr;
  QuenchError( XQuenchLogger::INFO, "insulation::strip::lz: " << lz << ", lr: " << lr );
}


double XCoilStrip :: GetTotalArea() const
{
  if (!fIns || !fSize) {
    XQuenchExcept except("insulation and strip size is not set.");
    throw except;
  }

  double lz = fSize[0];
  double lr = fSize[1] + 2*fIns[1];
  return lz * lr;
}


double XCoilStrip :: GetArea() const
{
  if (!fSize) {
    XQuenchExcept except("strip size is not set.");
    throw except;
  }

  return fSize[0] * fSize[1];
}


double XCoilStrip :: GetInsArea() const
{
  double total = GetTotalArea();
  double strip = GetArea();

  return total - strip;
}

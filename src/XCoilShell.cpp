#include <iostream>
#include "XQuenchExcept.hpp"
#include "XQuenchLogger.hpp"
#include "XCoilShell.hpp"

using Quench::XQuenchLogger;
using Quench::XCoilShell;

XCoilShell :: XCoilShell()
    : fSize(new double[2]),
      fIns(new double[2]),
      fType(kShell)
{}


XCoilShell :: ~XCoilShell()
{
  if (fSize)  delete [] fSize;
  if (fIns)   delete [] fIns;
}


void XCoilShell :: SetDimension(const double lz, const double lr)
{
  if (!fSize) fSize = new double[2];
  fSize[0] = lz;
  fSize[1] = lr;
  QuenchError( XQuenchLogger::INFO, "shell size -> lz: " << lz << ", lr: " << lr );
}


double XCoilShell :: GetDimension(const Coil dim) const
{
  switch (dim) {
    case iZ: return fSize[0];
    case iR: return fSize[1];
    default:
      QuenchError( XQuenchLogger::WARNING, "no this dimension: " << dim );
      return 0.;
  }
}


double XCoilShell :: GetTotalLength(const Coil dim) const
{
  switch (dim) {
    case iZ: return fSize[0];
    case iR: return fSize[1] + fIns[1];
    default:
      QuenchError( XQuenchLogger::WARNING, "no this dimension: " << dim );
      return 0.;
  }
}


void XCoilShell :: SetInsSize(const double lz, const double lr)
{
  if (!fIns)  fIns = new double[2];;
  fIns[0] = lz;
  fIns[1] = lr;
  QuenchError( XQuenchLogger::INFO, "insulation size -> lz: " << lz << ", lr: " << lr );
}


double XCoilShell :: GetInsSize(const Coil dim) const
{
  switch (dim) {
    case iZ: return fIns[0];
    case iR: return fIns[1];
    default:
      QuenchError( XQuenchLogger::WARNING, "no this dimension: " << dim );
      return 0.;
  }
}


double XCoilShell :: GetTotalArea() const
{
  double lz = fSize[0];        //  no insulation along z direction
  double lr = fIns[1] + fSize[1];   // 1 layer insulation
  return lz * lr;
}


double XCoilShell :: GetArea() const
{
  return fSize[0] * fSize[1];
}


double XCoilShell :: GetInsArea() const
{
  double total = GetTotalArea();
  double shell = GetArea();
  return total - shell;
}
